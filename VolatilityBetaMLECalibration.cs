/// <author>
/// Justin Chan
/// </author>
/// <owner>
/// Phil Koop
/// </owner>
/// <summary>
/// Calibration for the Beta Asset Price Volatility model using maximum likelihood (MLE) techniques. 
/// There are actually two calibrations, the first for a garden variety name volatility that 
/// will be driven by some systemic factor via a beta. This calibration is conceptually a 
/// combination of the Simple Volatility model calibration and the Multi-GBM Asset Price 
/// beta calibration in MLE framework. The second calibration is for the systemic factors 
/// themselves; it is implemented using the Simple Volatility model MLE calibration adapted to 
/// the Beta Asset Price Volatility model by setting the beta of the factor to itself equal to one.
/// </summary>
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibrate a Beta Asset Price Vol model using a beta to a systemic factor.
    /// </summary>
    [DisplayName("Beta Asset Price Volatility Beta Calibration")]
    public class BetaAssetPriceVolBetaCalibration : BetaAssetPriceVolSystemicCalibration
    {
        // Tag name used for IPriceFactorGrouping; the definition of systemic driver may be given in a grouping file, under column "Benchmark".
        private const string DriverMappingTagName = "Benchmark";

        private double fSignificanceThreshold;

        private List<TDate> fDriverHistoricalDates;
        private double[] fDriverProcessValues;
        private double[] fDriverProcessValuesAtHorizon;

        /// <summary>
        /// Initializes a new instance of the BetaAssetPriceVolBetaCalibration class.
        /// </summary>
        public BetaAssetPriceVolBetaCalibration()
        {
            Calibration_Method      = CalibrationMethod.MLE;
            MLE_Parameters          = new BetaAssetPriceVolBetaMLEParameters<AssetPriceVolSurface>(fVolatilityWindow);
            fSignificanceThreshold  = CalcUtils.DefaultMinProcessWeight;
        }

        /// <summary>
        /// Calibration method. The base property is overridden to ensure this property appears first as in similar calibrations.
        /// </summary>
        public new CalibrationMethod Calibration_Method { get; set; }

        /// <summary>
        /// The Driver used to lookup beta values.
        /// </summary>
        /// <remarks>
        /// If this is not set then the driver specified in the statistics set will be used.
        /// </remarks>
        public string Systemic_Factor                   { get; set; }

        /// <summary>
        /// Specify whether the calibration will allow an idiosyncratic factor.
        /// </summary>
        public YesNo Idiosyncratic_Factor               { get; set; }

        /// <summary>
        /// Specify the significance threshold of the model.
        /// </summary>
        public double Significance_Threshold
        {
            get
            {
                return fSignificanceThreshold;
            }
            set
            {
                if (value > 0.0)
                    fSignificanceThreshold = value;
                else
                    throw new AnalyticsException(string.Format("Significance threshold {0} cannot be non-positive", value));
            }
        }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new BetaAssetPriceVolBetaStatisticsCalibration
                {
                    Blend_Method            = Blend_Method,
                    Systemic_Factor         = Systemic_Factor,
                    Idiosyncratic_Factor    = Idiosyncratic_Factor,
                    Significance_Threshold  = Significance_Threshold
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            // Step(1): Validation and Initialisation
            BetaAssetPriceVolModel model        = CalibrationHelper.Validate<BetaAssetPriceVolModel>(calibrationData, priceModel, output);
            IPriceFactor priceFactor            = CalibrationHelper.GetPriceFactor(model);
            FactorID driverID                   = GetDriverID(calibrationData, priceFactor);
            ISpotProcessVol driverPriceFactor   = GetDriverPriceFactor(calibrationData, priceFactor, driverID);

            // Step(2): Base Case: calibrate as a systemic risk; future steps are mainly to figure out correlation
            // Most calibration work is done here. The rest of the calibration code is just to figure out correlation with the systemic driver.
            base.Calibrate(calibrationData, priceModel, output, errorList); 

            model.Systemic_Factor = driverID.ToCode();

            if (Idiosyncratic_Factor == YesNo.No || Math.Abs(model.Vol_of_Vol) < CalcUtils.TINY)
            {
                // If vol of vol == 0, then don't even bother computing correlation; this can happen when practically data are not updated and are kept constant.
                base.GetUpdatedCurves(priceModel, output);
                return;
            }

            TimeSeriesStructure<VolSurface.Point> driverVolSurfaceStructure;
            if (!GenericDataRetrieval.TryGetHistoricalVolSurface<AssetPriceVolSurface>(calibrationData, driverPriceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out driverVolSurfaceStructure))
            {
                return;
            }

            if (MLE_Parameters.Detrended_Time_Series == YesNo.Yes)
            {
                double[] drifts;
                IHolidayCalendar calendar = CalibrationHelper.GetHolidayCalendar(calibrationData, errorList, MLE_Parameters.Data_Retrieval_Parameters.Calendar);
                if (!TimeSeriesAnalysis.TryDetrendExponential(driverVolSurfaceStructure, errorList, calendar, out drifts))
                {
                    return;
                }
            }

            fDriverHistoricalDates = new List<TDate>(driverVolSurfaceStructure.Dates);

            var validDiverCoordinates = fVolatilityWindow.GetPointsInWindow(driverVolSurfaceStructure.TimeSeriesIDs.ToList());
            if (!validDiverCoordinates.Any())
                throw new CalibrationException(string.Format("No driver historical data of {1} for driver of {0} in selected window.", priceFactor.GetKey(), driverPriceFactor.GetKey()));

            var isValidDriverCoordinate = driverVolSurfaceStructure.TimeSeriesIDs.Select(point => validDiverCoordinates.Contains(point)).ToArray();

            // Step(3): calculate observations of average driver process
            fDriverProcessValues          = new double[driverVolSurfaceStructure.Length];
            fDriverProcessValuesAtHorizon = new double[driverVolSurfaceStructure.Length];
            for (int i = 0; i < driverVolSurfaceStructure.Length; i++)
            {
                double sum          = 0.0;
                double sumAtHorizon = 0.0;
                for (int j = 0; j < driverVolSurfaceStructure.NumTimeSeries; j++)
                {
                    if (!isValidDriverCoordinate[j])
                        continue;

                    sum          += Math.Log(driverVolSurfaceStructure.Values[j][i]);
                    sumAtHorizon += Math.Log(driverVolSurfaceStructure.ValuesAtHorizon[j][i]);
                }

                fDriverProcessValues[i]          = sum / validDiverCoordinates.Count;
                fDriverProcessValuesAtHorizon[i] = sumAtHorizon / validDiverCoordinates.Count;
            }

            var parameters = (BetaAssetPriceVolBetaMLEParameters<AssetPriceVolSurface>)MLE_Parameters;
            if (parameters == null)
                throw new CalibrationException("Unexpected MLE parameters type.");

            // Step(4): calculate alpha of average driver process
            var driverAlphaFixed = parameters.Condition_On_Systemic_Factor_Reversion_Speed == YesNo.Yes ? (double?)parameters.Systemic_Factor_Reversion_Speed : null;
            MeanReversionStatistics driverMrs;
            if (!ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(fDriverProcessValues, fDriverProcessValuesAtHorizon, ReturnMethod.Diff,
                                                                                    driverVolSurfaceStructure.DeltaT, driverAlphaFixed, null, false, errorList,
                                                                                    out driverMrs))
            {
                errorList.Add(ErrorLevel.Error, "Calibration Error");
                return;
            }

            var driverAlpha = driverMrs.Alpha;

            // Step(5): get observations of average process and average driver process on common set of dates
            IEnumerable<TDate> commonDates = fHistoricalDates.Where(date => fDriverHistoricalDates.Contains(date));

            var processValues                = new List<double>(commonDates.Count());
            var processValuesAtHorizon       = new List<double>(commonDates.Count());
            var driverProcessValues          = new List<double>(commonDates.Count());
            var driverProcessValuesAtHorizon = new List<double>(commonDates.Count());

            foreach (var date in commonDates)
            {
                var index = fHistoricalDates.IndexOf(date); 
                var driverIndex = fDriverHistoricalDates.IndexOf(date);
                processValues.Add(fProcessValues[index]);
                processValuesAtHorizon.Add(fProcessValuesAtHorizon[index]);
                driverProcessValues.Add(fDriverProcessValues[driverIndex]);
                driverProcessValuesAtHorizon.Add(fDriverProcessValuesAtHorizon[driverIndex]);
            }

            // Step(6): calculate correlation between driver process and average driver process
            var values          = new double[][] { processValues.ToArray(), driverProcessValues.ToArray() };
            var valuesAtHorizon = new double[][] { processValuesAtHorizon.ToArray(), driverProcessValuesAtHorizon.ToArray() };
            var alphas          = new double[] { model.Reversion_Speed, driverAlpha };

            var covariance = ReversionStatisticsCalculationHelper.ComputeCovarianceMatrix(values, valuesAtHorizon, alphas, null, ReturnMethod.Diff, driverVolSurfaceStructure.DeltaT, false);
            double correlation = CalcUtils.IsTiny(covariance[1, 1]) ? 0.0 : covariance[1, 0] / Math.Sqrt(covariance[0, 0] * covariance[1, 1]);

            model.Systemic_Weight = GetSystemicWeight(correlation);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        public override void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            if (fHistoricalDates == null || fHistoricalDates.Count < StatisticsCalculationHelper.MinTimeSeriesLength)
                return;

            var graphData = new ResultSeriesList();
            var dependentOUProcess  = new Curve();

            graphData.fXAxis.fLabel = "Date";
            graphData.fYAxis.fLabel = "Underlying Ornstein-Uhlenbeck Processes";

            for (int i = 0; i < fHistoricalDates.Count; i++)
                dependentOUProcess[(float)fHistoricalDates[i]] = fProcessValues[i];

            graphData.Add(new ResultSeries(string.Format("Dependent OU Process"), LineStyle.Solid, dependentOUProcess));

            if (fDriverHistoricalDates != null && fDriverHistoricalDates.Count >= StatisticsCalculationHelper.MinTimeSeriesLength)
            {
                var driverOUProcess = new Curve();
                for (int i = 0; i < fDriverHistoricalDates.Count; i++)
                    driverOUProcess[(float)fDriverHistoricalDates[i]] = fDriverProcessValues[i];

                graphData.Add(new ResultSeries(string.Format("Systemic OU Process"), LineStyle.Solid, driverOUProcess));
            }

            graphData.ToXml(output);
        }

        /// <summary>
        /// Returns the driver factor ID.
        /// </summary>
        /// <remarks>
        /// If the input field "Driver" is populated, the given Driver is used. Otherwise, this method goes into historical data archive (e.g., ADA file),
        /// finds all price factor's corresponding points, and finds their matching pairs in the factor grouping (e.g., Grouping file). Then, the method
        /// determines the systemic driver, and returns its FactorID. This method should never return null; if it cannot find a logical FactorID, it will throw errors.
        /// </remarks>
        private FactorID GetDriverID(ICalibrationData calibrationData, IPriceFactor priceFactor)
        {
            if (!string.IsNullOrEmpty(Systemic_Factor))
            {
                return new FactorID(Systemic_Factor);
            }

            // Reading from grouping file
            IPriceFactorGrouping priceFactorGrouping = calibrationData.Grouping;

            string priceFactorType  = priceFactor.TypeDotSubType();
            FactorID priceFactorID  = priceFactor.GetID();
            FactorID driverFactorID = null;

            var factors = new FactorTypeAndIDList { priceFactor };

            IHistoricalArchive archive = calibrationData.GetHistoricalMarketRates(factors, calibrationData.ArchiveStartDate, calibrationData.ArchiveEndDate);

            foreach (ArchiveItemID item in archive.ItemIDs)
            {
                if (StatisticsData.GetTypeForFactor(item.fPriceFactorSubType) == priceFactorType &&
                    item.fPriceFactorID.CompareTo(priceFactorID) == 0)
                {
                    string riskFactorString = StatisticsData.BuildScalarID(priceFactorType, priceFactorID, item.fPointID);
                    string candidateDriver  = priceFactorGrouping.GetTagValue(riskFactorString, DriverMappingTagName);

                    string candidateDriverFactorType;
                    string candidateDriverPointID;
                    FactorID candidateDriverFactorID;
                    StatisticsData.ParseScalarID(candidateDriver, out candidateDriverFactorType,
                                                 out candidateDriverFactorID, out candidateDriverPointID);

                    if (string.IsNullOrEmpty(candidateDriverFactorID.ToCode()) || priceFactorType != candidateDriverFactorType)
                        continue;

                    if (driverFactorID == null)
                        driverFactorID = candidateDriverFactorID;
                    else if (driverFactorID.ToCode() != candidateDriverFactorID.ToCode())
                        throw new CalibrationException(
                            string.Format(
                                "Ambiguous {0} definition of systemic driver for {1}: mapped at least to both {2} and {3} in grouping file.",
                                DriverMappingTagName, priceFactor.GetKey(), driverFactorID.ToCode(), candidateDriverFactorID.ToCode()));
                }
            }

            if (driverFactorID == null)
                throw new CalibrationException(string.Format(
                        "Info about either systemic driver for {0} in grouping file or historical data for {0} is not available.",
                        priceFactor.GetKey()));

            return driverFactorID;
        }

        /// <summary>
        /// Gets the driver price factor from the market data (e.g., DAT file).
        /// </summary>
        /// <remarks>
        /// In addition to finding a matching price factor, we especially need the current vol surface of the driver to calibrate the driver process and calculate beta.
        /// This method will never return null; if it cannot find a logical driverPriceFactor, it will throw errors.
        /// </remarks>
        private ISpotProcessVol GetDriverPriceFactor(ICalibrationData calibrationData, IPriceFactor priceFactor, FactorID driverID)
        {
            if (driverID == null)
                return null;

            ISpotProcessVol driverPriceFactor = null;
            foreach (IPriceFactor pf in calibrationData.MarketData.PriceFactors)
            {
                if (driverPriceFactor != null)
                    break;

                if (pf.GetID().CompareTo(driverID) == 0)
                    driverPriceFactor = pf as ISpotProcessVol;
            }

            if (driverPriceFactor == null)
                throw new CalibrationException(
                    string.Format("The driver {0} for {1} of type {2} doesn't exist in market data.", Systemic_Factor,
                                  priceFactor.GetKey(), typeof(ISpotProcessVol)));

            return driverPriceFactor;
        }

        /// <summary>
        /// Coerce the given correlation above the Significance Threshold.
        /// </summary>
        private double GetSystemicWeight(double correlation)
        {
            if (Math.Abs(correlation) < Significance_Threshold)
                return 0.0;

            if (Math.Abs(correlation - 1) < Significance_Threshold)
                return 1.0;

            return correlation;
        }
    }

    /// <summary>
    /// Calibrate a Beta Asset Price Vol model as a systemic factor, i.e. without using an input beta.
    /// </summary>
    /// <remarks>
    /// This calibration is equivalent to a simple volatility MLE calibration.
    /// </remarks>
    [DisplayName("Systemic Beta Asset Price Volatility Calibration")]
    public class BetaAssetPriceVolSystemicCalibration : SimpleAssetPriceVolCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new BetaAssetPriceVolSystemicStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            BetaAssetPriceVolModel model = (BetaAssetPriceVolModel)priceModel;

            model.Systemic_Factor   = model.fID.ToCode();
            model.Systemic_Weight   = 1.0;

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(BetaAssetPriceVolSystemicCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public override Type ModelType()
        {
            return typeof(BetaAssetPriceVolModel);
        }
    }

    /// <summary>
    /// Parameters for the MLE method.
    /// </summary>
    /// <typeparam name="TVolSurface">A VolSurface object.</typeparam>
    public class BetaAssetPriceVolBetaMLEParameters<TVolSurface> : SimpleVolMLEParameters<TVolSurface> where TVolSurface : VolSurface, new()
    {
        private double fDriverAlpha;

        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public BetaAssetPriceVolBetaMLEParameters(SurfaceWindow<TVolSurface> surfaceWindow) : base(surfaceWindow)
        {
            Condition_On_Systemic_Factor_Reversion_Speed = YesNo.No;
        }

        /// <summary>
        /// Optionally using conditional Reversion Speed (i.e., Alpha) for the driver (systemic risk) process; output parameters will be conditional on the given Alpha.
        /// </summary>
        public YesNo Condition_On_Systemic_Factor_Reversion_Speed   { get; set; }

        /// <summary>
        /// Optional user-defined Reversion_Speed (i.e., Alpha) for conditional MLE estimation for the driver(systemic risk) process; used in conjunction with Condition_On_Driver_Reversion_Speed.
        /// </summary>
        public double Systemic_Factor_Reversion_Speed
        {
            get
            {
                return fDriverAlpha;
            }
            set
            {
                fDriverAlpha = Math.Abs(value) >= CalcUtils.MinReversionRate ? value : 0.0;
            }
        }
    }
}