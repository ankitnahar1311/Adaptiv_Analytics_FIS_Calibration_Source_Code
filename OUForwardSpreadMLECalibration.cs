/// <author>
/// Justin Chan
/// </author>
/// <owner>
/// Philip Koop
/// </owner>
/// <summary>
/// Maximum likelihood estimation calibrations for OU forward spread models.
/// </summary>
using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Base class for calibratiohn of OU models with maximum likelihood estimation.
    /// </summary>
    public abstract class OUForwardSpreadCalibrationBase : IModelCalibration
    {
        protected TimeSeriesStructure<Period> fTimeSeries;
        protected double[] fPointDates; // List of dates corresponding to points (e.g., 1/1/2011, 6/1/2011, 1/1/2012, etc)

        /// <summary>
        /// Initializes a new instance of the OUForwardSpreadMLECalibrationBase class.
        /// </summary>
        protected OUForwardSpreadCalibrationBase()
        {
            // Use MLE by default.
            Calibration_Method  = CalibrationMethod.MLE;
        }

        /// <summary>
        /// Calibration method can be statistics or MLE.
        /// </summary>
        public CalibrationMethod Calibration_Method             { get; set; }
        
        /// <summary>
        /// Base date. If not set, we use the base date provided in the market data.
        /// </summary>
        public TDate Base_Date              { get; set; }

        /// <summary>
        /// Last rollover date.
        /// </summary>
        public TDate Last_Rollover_Date     { get; set; }

        /// <summary>
        /// Parameters for the MLE method.
        /// </summary>
        public OUForwardSpreadBaseMLEParameters MLE_Parameters  { get; set; }

        /// <summary>
        /// Calibrates the price model from the statistics files.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Step(1): Validation and initialisation
            ForwardPriceModel model     = CalibrationHelper.Validate<ForwardPriceModel>(calibrationData, priceModel, output);
            ForwardPrice currentSpread  = (ForwardPrice)CalibrationHelper.GetPriceFactor(model);
            IPriceFactor priceFactor    = CalibrationHelper.GetPriceFactor(model);

            if (Last_Rollover_Date > Base_Date)
                errorList.Add(ErrorLevel.Warning, string.Format("Last Rollover Date {0} should not be after Base Date {1} for price factor {2}", Last_Rollover_Date, Base_Date, priceFactor.GetKey()));

            if (currentSpread.Curve.Count < 1)
                throw new CalibrationException(string.Format("Current spread curve for {0} must be populated.", priceFactor.GetKey()));

            // Extract a base date.
            TDate baseDate = Base_Date == 0.0 ? (TDate)currentSpread.Curve.X[0] : Base_Date;

            // If these dates are empty, then assume they are the first date on the spread curve.
            TDate lastRolloverDate = Last_Rollover_Date == 0.0 ? (TDate)currentSpread.Curve.X[0] : Last_Rollover_Date;

            // Step(2): Obtain raw data and basic data processing
            if (!GenericDataRetrieval.TryGetHistoricalRateTermStructure(calibrationData, priceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out fTimeSeries))
            {
                return;
            }

            // Currently using simple detrending (i.e., no seasonality adjustment).
            if (MLE_Parameters.Detrended_Time_Series == YesNo.Yes)
            {
                double[] drifts;

                IHolidayCalendar calendar = MLE_Parameters.Data_Retrieval_Parameters.GetHolidayCalendar(calibrationData, errorList);

                if (!TimeSeriesAnalysis.TryDetrendLinear(fTimeSeries, errorList, calendar, out drifts))
                {
                    errorList.Add(ErrorLevel.Error, "Calibration Error");
                    return;
                }
            }

            fPointDates = fTimeSeries.TimeSeriesIDs.Select(point => baseDate + CalcUtils.YearsToDays(point)).ToArray();

            // Step(3): Core model calibration; to be specified at each individual model
            CalibrateModelParameters(lastRolloverDate, model, errorList);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            if (fTimeSeries == null || fTimeSeries.Length < StatisticsCalculationHelper.MinTimeSeriesLength)
                return;

            var  graphData = new ResultSeriesList
                                 {
                                     fXAxis = { fLabel = "Date" },
                                     fYAxis =
                                         {
                                             fLabel = MLE_Parameters.Detrended_Time_Series == YesNo.Yes
                                                     ? "Detrended Price Spread" : "Price Spread"
                                         }
                                 };

            for (int i = 0; i < fTimeSeries.NumTimeSeries; i++)
            {
                // For each time series of price spread at a particular point, plot a curve
                var historicalData = new Curve();
                for (int j = 0; j < fTimeSeries.Length; j++)
                    historicalData[fTimeSeries.Dates[j]] = fTimeSeries.Values[i][j];

                graphData.Add(new ResultSeries(string.Format(CultureInfo.InvariantCulture, "Term = {0}", (float)fTimeSeries.TimeSeriesIDs[i]), LineStyle.Solid, historicalData));
            }

            graphData.ToXml(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        public abstract Type ModelType();

        /// <summary>
        /// Core calibration method to be defined at each concrete child.
        /// </summary>
        protected abstract void CalibrateModelParameters(TDate lastRolloverDate, ForwardPriceModel model, ErrorList errorList);
    }

    /// <summary>
    /// Calibrate the OU Forward Spread TS model with maximum likelihood estimation.
    /// </summary>
    [DisplayName("Ornstein-Uhlenbeck Forward Spread Term Structure Calibration")]
    public class OUForwardSpreadTSCalibration : OUForwardSpreadCalibrationBase
    {
        /// <summary>
        /// Initializes a new instance of the OUForwardSpreadTSMLECalibration class.
        /// </summary>
        public OUForwardSpreadTSCalibration()
        {
            MLE_Parameters = new OUForwardSpreadTSMLEParameters();
        }

        /// <summary>
        /// Calibrates the price model from the statistics files.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new OUForwardSpreadTSStatisticsCalibration
                                                {
                                                    Base_Date           = Base_Date,
                                                    Last_Rollover_Date  = Last_Rollover_Date
                                                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
            }
            else
            {
                base.Calibrate(calibrationData, priceModel, output, errorList);
            }
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public override Type ModelType()
        {
            return typeof(OUForwardSpreadTSModel);
        }

        /// <summary>
        /// Core calibration method to be defined at each concrete child.
        /// </summary>
        protected override void CalibrateModelParameters(TDate lastRolloverDate, ForwardPriceModel model, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            var parameters = (OUForwardSpreadTSMLEParameters)MLE_Parameters;

            // Step(1): Validation and Initialisation
            if (parameters.Max_Reversion_Speed > 0.0 && parameters.Max_Reversion_Speed < parameters.Min_Reversion_Speed)
                errorList.Add(ErrorLevel.Error, string.Format("Maximum Reversion Speed {0} cannot be less than Minimum Reversion Speed {1} for price factor {2}.",
                        parameters.Max_Reversion_Speed, parameters.Min_Reversion_Speed, priceFactor));

            var price = (ForwardPrice)CalibrationHelper.GetPriceFactor(model);

            var reversionSpeed = new PriceCurve();
            var reversionLevel = new PriceCurve();
            var volatility     = new PriceCurve();

            double? alphaMax   = null;
            double? alphaMin   = null;
            double? alphaFixed = null;
            if (parameters.Min_Reversion_Speed > 0.0)
            {
                // Min_Reversion_Speed specified
                if (parameters.Max_Reversion_Speed == parameters.Min_Reversion_Speed)
                {
                    // Calibrate with known alpha
                    alphaFixed = parameters.Max_Reversion_Speed;
                }
                else
                {
                    alphaMin = parameters.Min_Reversion_Speed;
                    if (parameters.Max_Reversion_Speed > parameters.Min_Reversion_Speed)
                        alphaMax = parameters.Max_Reversion_Speed;
                }
            }

            // Mean Reversion Level assumed to be the current term structure
            reversionLevel.Assign(price.Curve);  
            var thetaFixed = new double?[fTimeSeries.NumTimeSeries];
            for (int i = 0; i < fTimeSeries.NumTimeSeries; ++i)
                thetaFixed[i] = reversionLevel[fPointDates[i]];

            // Step(2): Calibrate the OU processes
            var mrsParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Diff, alphaMin, alphaMax, null, alphaFixed, thetaFixed, false);
            CorrelatedMeanReversionStatisticsVectorAlpha results;
            ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(fTimeSeries.Values, fTimeSeries.ValuesAtHorizon, fTimeSeries.DeltaT, mrsParameters, errorList, out results);

            for (int i = 0; i < fTimeSeries.NumTimeSeries; ++i)
            {
                reversionSpeed[fPointDates[i]] = results.MeanReversionSpeeds[i];
                volatility[fPointDates[i]]     = results.MeanReversionVols[i];
            }

            // Step(3): Model parameters assignment
            var spreadModel = model as OUForwardSpreadTSModel;
            if (spreadModel == null)
                throw new CalibrationException(string.Format("This Calibration can only calibrate model of type {0} but not {1} for price factor {2}.", typeof(OUForwardSpreadTSModel), model.GetType(), priceFactor.GetKey()));

            spreadModel.Last_Rollover_Date = lastRolloverDate;
            spreadModel.Reversion_Speed    = reversionSpeed;
            spreadModel.Reversion_Level    = reversionLevel;
            spreadModel.Volatility         = volatility;
        }
    }

    /// <summary>
    /// Calibrate the OU Forward Spread model with maximum likelihood estimation.
    /// </summary>
    [DisplayName("Ornstein-Uhlenbeck Forward Spread Calibration")]
    public class OUForwardSpreadCalibration : OUForwardSpreadCalibrationBase
    {
        /// <summary>
        /// Initializes a new instance of the OUForwardSpreadMLECalibration class.
        /// </summary>
        public OUForwardSpreadCalibration()
        {
            MLE_Parameters = new OUForwardSpreadMLEParameters();
        }

        /// <summary>
        /// Calibrates the price model from the statistics files.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new OUForwardSpreadStatisticsCalibration
                                                {
                                                    Base_Date           = Base_Date,
                                                    Last_Rollover_Date  = Last_Rollover_Date
                                                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
            }
            else
            {
                base.Calibrate(calibrationData, priceModel, output, errorList);
            }
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public override Type ModelType()
        {
            return typeof(OUForwardSpreadModel);
        }

        /// <summary>
        /// Core calibration method to be defined at each concrete child.
        /// </summary>
        protected override void CalibrateModelParameters(TDate lastRolloverDate, ForwardPriceModel model, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);
            ForwardPrice price      = (ForwardPrice)CalibrationHelper.GetPriceFactor(model);

            var parameters = (OUForwardSpreadMLEParameters)MLE_Parameters;

            var reversionLevel = new PriceCurve();
            var volatility     = new PriceCurve();

            // Mean Reversion Level assumed to be the current term structure.
            reversionLevel.Assign(price.Curve); 
            var thetaFixed = new double?[fTimeSeries.NumTimeSeries];
            for (int i = 0; i < fTimeSeries.NumTimeSeries; ++i)
                thetaFixed[i] = reversionLevel[fPointDates[i]];

            var alphaFixed = parameters.Condition_On_Reversion_Speed == YesNo.Yes ? (double?)parameters.Reversion_Speed : null;

            var mrsParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Diff, null, null, null, alphaFixed, thetaFixed, false);
            CorrelatedMeanReversionStatisticsScalarAlpha results;
            ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(fTimeSeries.Values, fTimeSeries.ValuesAtHorizon, fTimeSeries.DeltaT, mrsParameters, errorList, out results);

            for (int i = 0; i < fTimeSeries.NumTimeSeries; ++i)
            {
                volatility[fPointDates[i]] = results.MeanReversionVols[i];
            }

            var spreadModel = model as OUForwardSpreadModel;
            if (spreadModel == null)
                throw new CalibrationException(string.Format("This Calibration can only calibrate model of type {0} but not {1} for price factor {2}.", typeof(OUForwardSpreadModel), model.GetType(), priceFactor.GetKey()));

            spreadModel.Last_Rollover_Date = lastRolloverDate;
            spreadModel.Reversion_Speed    = results.MeanReversionSpeed;
            spreadModel.Reversion_Level    = reversionLevel;
            spreadModel.Volatility         = volatility;
        }
    }

    /// <summary>
    /// Parameters for the MLE method for OUForwardSpreadTSCalibration.
    /// </summary>
    public class OUForwardSpreadTSMLEParameters : OUForwardSpreadBaseMLEParameters
    {
        private double fMinReversionSpeed;
        private double fMaxReversionSpeed;

        /// <summary>
        /// If calibrated reversion speed is lower than specified min reversion speed, then process will be re-calibrated conditional on given min reversion speed.
        /// </summary>
        /// <remark>
        /// Set both Min_Reversion_Speed and Max_Reversion_Speed to 0.0 to ignore reversion speed bounds.
        /// </remark>
        public double Min_Reversion_Speed
        {
            get
            {
                return fMinReversionSpeed;
            }
            set
            {
                if (value >= 0.0)
                    fMinReversionSpeed = value;
                else
                    throw new AnalyticsException(string.Format("Minimum Reversion Speed {0} cannot be non-positive", value));
            }
        }

        /// <summary>
        /// If calibrated reversion speed is lower than specified min reversion speed, then process will be re-calibrated conditional on given min reversion speed.
        /// </summary>
        /// <remark>
        /// Set both Min_Reversion_Speed and Max_Reversion_Speed to 0.0 to ignore reversion speed bounds.
        /// </remark>
        public double Max_Reversion_Speed
        {
            get
            {
                return fMaxReversionSpeed;
            }
            set
            {
                if (value >= 0.0)
                    fMaxReversionSpeed = value;
                else
                    throw new AnalyticsException(string.Format("Maximum Reversion Speed {0} cannot be non-positive", value));
            }
        }
    }

    /// <summary>
    /// Parameters for the MLE method for OUForwardSpreadTSCalibration.
    /// </summary>
    public class OUForwardSpreadMLEParameters : OUForwardSpreadBaseMLEParameters
    {
        private double fAlpha;

        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public OUForwardSpreadMLEParameters()
        {
            Condition_On_Reversion_Speed = YesNo.No;
        }

        /// <summary>
        /// Optionally using conditional Reversion Speed (i.e., Alpha); 
        /// output parameters will be conditional on the given Alpha.
        /// </summary>
        public YesNo Condition_On_Reversion_Speed { get; set; }

        /// <summary>
        /// Optional user-defined Reversion_Speed (i.e., Alpha) for conditional MLE estimation; 
        /// used in conjunction with Condition_On_Reversion_Speed.
        /// </summary>
        public double Reversion_Speed
        {
            get
            {
                return fAlpha;
            }
            set
            {
                fAlpha = Math.Abs(value) >= CalcUtils.MinReversionRate ? value : 0.0;
            }
        }
    }

    /// <summary>
    /// Parameters for the MLE method for base OUForwardSpreadCalibration.
    /// </summary>
    public abstract class OUForwardSpreadBaseMLEParameters : NestedPresentableObject
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        protected OUForwardSpreadBaseMLEParameters()
        {
            Data_Retrieval_Parameters = new CalibratorDataRetrievalParameters<Period>();
            Detrended_Time_Series     = YesNo.Yes;
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<Period> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// Optionally remove the timeseries' drift before MLE calibration; this is a calibration setting.
        /// </summary>
        public YesNo Detrended_Time_Series                               { get; set; }
    }
}