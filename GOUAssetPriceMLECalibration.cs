/// <author>
/// Justin Chan
/// </author>
/// <owner>
/// Phil Koop
/// </owner>
/// <summary>
/// GBM Asset Price Mean Reversion Model Calibration, using maximum likelihood estimation (MLE). 
/// For this particular model, there happens to exist a closed-form solution for the model parameters.
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
    /// Calibration for GOUAssetPriceModel using (conditional) MLE method.
    /// </summary>
    /// <remarks>
    /// Set Calibration_Method = MLE for the MLE method and Calibration_Method = Pre_Computed_Statistics 
    /// to use pre-computed statistics from the statistics file.
    /// </remarks>
    [DisplayName("GOU Asset Price Calibration")]
    public class GOUAssetPriceCalibration : IModelCalibration
    {
        // Data needed for graph
        private TDate[] fHistoricalDates;
        private double[] fHistoricalPrices;
        private double[] fDetrendedHistoricalPrices;
        private double fLongRunMean;

        /// <summary>
        /// Initializes a new instance of the GBMAssetPriceMRMLECalibration class.
        /// </summary>
        /// <remarks>Use MLE by default.</remarks>
        public GOUAssetPriceCalibration()
        {
            // Initialise default values
            Calibration_Method  = CalibrationMethod.MLE;               
            MLE_Parameters      = new GOUAssetPriceMLEParameters();

            Use_Absolute_Theta  = YesNo.Yes;
            Use_Mu              = YesNo.No;
        }

        /// <summary>
        /// Calibration method can be statistics or MLE.
        /// </summary>
        public CalibrationMethod Calibration_Method         { get; set; }

        /// <summary>
        /// Using either absolute(Yes) or relative(No) mean reversion levels; directly passing this attribute to the model.
        /// </summary>
        /// <value>
        /// Yes for using absolute mean reversion level; no for using relative mean reversion level; default is Yes.
        /// </value>
        public YesNo Use_Absolute_Theta                     { get; set; }

        /// <summary>
        /// Using(Yes) or ignoring(No) the drift term; directly passing this attribute to the model.
        /// </summary>
        /// <value>
        /// Yes for using a drift term; no for ignoring the drift term; default is No.
        /// </value>
        public YesNo Use_Mu                                 { get; set; }

        /// <summary>
        /// MLE parameters.
        /// </summary>
        public GOUAssetPriceMLEParameters MLE_Parameters    { get; set; }

        /// <summary>
        /// Calibrates the GBM Asset Price model with mean reversion from the statistics.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new GOUAssetPriceStatisticsCalibration
                {
                    Use_Absolute_Theta  = Use_Absolute_Theta,
                    Use_Mu              = Use_Mu
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            GOUAssetPriceModel model = CalibrationHelper.Validate<GOUAssetPriceModel>(calibrationData, priceModel, output);
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            // TIME SERIES.
            TimeSeriesStructure<string> timeSeries;
            if (!GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData, priceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out timeSeries))
            {
                return;
            }

            fHistoricalDates = timeSeries.Dates;
            model.Mu = 0.0;
            if (MLE_Parameters.Detrended_Time_Series == YesNo.Yes)
            {
                // Copy the values before detrending (for graph)
                fHistoricalPrices = (double[])timeSeries.Values[0].Clone();

                double[] drifts;
                IHolidayCalendar calendar = MLE_Parameters.Data_Retrieval_Parameters.GetHolidayCalendar(calibrationData, errorList);
                if (!TimeSeriesAnalysis.TryDetrendExponential(timeSeries, errorList, calendar, out drifts))
                {
                    return;
                }

                if (Use_Mu == YesNo.Yes)
                    model.Mu = drifts[0];

                fDetrendedHistoricalPrices = timeSeries.Values[0];
            }
            else
            {
                fHistoricalPrices = timeSeries.Values[0];
            }

            double logInitialPrice = Use_Absolute_Theta == YesNo.Yes ? 0.0 : Math.Log(((IAssetPrice)(model.GetPriceFactor())).GetQuotedSpotPrice()[0]);
            var alphaFixed = MLE_Parameters.Condition_On_Alpha == YesNo.Yes ? (double?)MLE_Parameters.Alpha : null;
            var thetaFixed = MLE_Parameters.Condition_On_Theta == YesNo.Yes ? (double?)MLE_Parameters.Theta + logInitialPrice : null;

            MeanReversionStatistics mrs;
            if (!ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(timeSeries.Values[0], timeSeries.ValuesAtHorizon[0], ReturnMethod.Log, timeSeries.DeltaT, alphaFixed, thetaFixed, true, errorList, out mrs))
                return;

            model.Theta               = mrs.Theta - logInitialPrice;
            model.Alpha               = mrs.Alpha;
            model.Sigma               = mrs.Sigma;
            model.Use_Absolute_Theta  = Use_Absolute_Theta;

            // Store long run mean of log(S) for graph
            fLongRunMean = mrs.Theta;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            if (fHistoricalDates == null || fHistoricalDates.Length < StatisticsCalculationHelper.MinTimeSeriesLength)
                return;

            GOUAssetPriceModel model        = (GOUAssetPriceModel)priceModel;

            var historicalData              = new Curve();
            var detrendedHistoricalData     = new Curve();
            var longRunMeanPlot             = new Curve();
            var longRunOneStdevUpperRange   = new Curve();
            var longRunOneStdevLowerRange   = new Curve();

            var times = fHistoricalDates.Select(date => CalcUtils.DaysToYears(date - fHistoricalDates[0])).ToList();
            fLongRunMean = MeanReversionStatistics.GetThetaOfGouProcess(model.Alpha, model.Sigma, fLongRunMean);

            var graphData = new ResultSeriesList
                                {
                                    fXAxis = { fLabel = "Date" },
                                    fYAxis = { fLabel = "Historical Price" }
                                };

            for (int i = 0; i < fHistoricalDates.Length; i++)
            {
                historicalData[fHistoricalDates[i]] = fHistoricalPrices[i];

                if (model.Alpha > 0.0)
                {
                    longRunMeanPlot[fHistoricalDates[i]]           = fLongRunMean * Math.Exp(model.Mu * times[i]);
                    longRunOneStdevUpperRange[fHistoricalDates[i]] = longRunMeanPlot[fHistoricalDates[i]] * Math.Exp(model.Sigma / Math.Sqrt(2.0 * model.Alpha));
                    longRunOneStdevLowerRange[fHistoricalDates[i]] = longRunMeanPlot[fHistoricalDates[i]] * Math.Exp(-model.Sigma / Math.Sqrt(2.0 * model.Alpha));
                }

                if (MLE_Parameters.Detrended_Time_Series == YesNo.Yes && fDetrendedHistoricalPrices != null)
                    detrendedHistoricalData[fHistoricalDates[i]] = fDetrendedHistoricalPrices[i];
            }

            graphData.Add(new ResultSeries(string.Format("Historical Asset Prices"), LineStyle.Solid, historicalData));
            if (detrendedHistoricalData.Count != 0)
                graphData.Add(new ResultSeries(string.Format("Detrended Historical Asset Prices"), LineStyle.Solid, detrendedHistoricalData));

            // Long run mean, and upper/lower range of variation is only defined when mean reversion speed is positive
            if (model.Alpha > 0.0)
            {
                graphData.Add(new ResultSeries(string.Format("Long Run Mean"), LineStyle.Solid, longRunMeanPlot));
                graphData.Add(new ResultSeries(string.Format("Long Run One Standard Deviation Upper Range"), LineStyle.Solid, longRunOneStdevUpperRange));
                graphData.Add(new ResultSeries(string.Format("Long Run One Standard Deviation Lower Range"), LineStyle.Solid, longRunOneStdevLowerRange));
            }

            graphData.ToXml(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GOUAssetPriceModel);
        }
    }

    /// <summary>
    /// Parameters for the MLE method.
    /// </summary>
    public class GOUAssetPriceMLEParameters : NestedPresentableObject
    {
        private double fAlpha;

        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public GOUAssetPriceMLEParameters()
        {
            Data_Retrieval_Parameters = new CalibratorDataRetrievalParameters<string>();

            Detrended_Time_Series       = YesNo.Yes;
            Condition_On_Alpha          = YesNo.No;
            Condition_On_Theta          = YesNo.No;
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<string> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// Optionally remove the timeseries' drift before MLE calibration; this is a calibration setting.
        /// </summary>
        /// <value>
        /// Yes for detrending before calibration; default is No.
        /// </value>
        public YesNo Detrended_Time_Series { get; set; }

        /// <summary>
        /// Optionally using conditional alpha; i.e., output parameters will be conditional on the given Alpha.
        /// </summary>
        /// <value>
        /// Yes for using Alpha for conditional MLE calibration; default is Yes.
        /// </value>
        public YesNo Condition_On_Alpha { get; set; }

        /// <summary>
        /// Optional user-defined Alpha for conditional MLE estimation; used in conjunction with Condition_On_Alpha.
        /// </summary>
        /// <value>
        /// Value of Alpha; since this is a customisation setting; any value, including negative Alpha, is allowed.
        /// </value>
        public double Alpha
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

        /// <summary>
        /// Optionally using conditional theta; i.e., output parameters will be conditional on the given Theta.
        /// </summary>
        /// <value>
        /// Yes for using Theta for conditional MLE calibration; default is No.
        /// </value>
        public YesNo Condition_On_Theta { get; set; }

        /// <summary>
        /// Optional user-defined Theta for conditional MLE estimation; used in conjunction with Conditional_On_Theta.
        /// </summary>
        /// <value>
        /// Value of Theta; since this is a customisation setting; any value is allowed.
        /// </value>
        public double Theta { get; set; }
    }
}