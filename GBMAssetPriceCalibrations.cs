/// <author>
/// Alastair Wilkins, Nathalie Antunes, Perukrishnen Vytelingum 
/// </author>
/// <owner>
/// Nathalie Antunes
/// </owner>
/// <summary>
///  GBMAssetPriceCalibration :         
///          Standard calibration for a single factor GBM model with constant 
///          drift and vol
///      
///   GBMPriceIndexCalibration :
///          Standard calibration for a GBM index price model
///
///   GBMPriceIndexDriftCalibration :
///          Standard calibration for a GBM index price drift model
///
///   GBMAssetPriceMRCalibration :
///          Standard calibration for a GBM mean reverting model
/// </summary>
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration of the GBM Asset Price Model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("GBM Asset Price Calibration")]
    public class GBMAssetPriceCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the GBM Asset Price model using the statistics.
        /// </summary>
        /// <remarks>
        /// Computed from archive data or obtained from the statistics set
        /// </remarks>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            GBMAssetPriceModel model    = CalibrationHelper.Validate<GBMAssetPriceModel>(calibrationData, priceModel, output);
            var priceFactor             = CalibrationHelper.GetPriceFactor(priceModel);

            // Moment statistics.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, model, ReturnMethod.Log, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            // Set model parameters using earliest point. Validating earliest point.
            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for {0}.", priceFactor.GetKey()));
                return;
            }

            model.Drift = momentStatistics[0].GetGbmDrift();
            model.Vol   = momentStatistics[0].Vol;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GBMAssetPriceModel);
        }
    }

    /// <summary>
    /// Calibration of the GBM Price Index Model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("GBM Price Index Calibration")]
    public class GBMPriceIndexCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the GBM Asset Price model from the statistics.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            GBMPriceIndexModel model    = CalibrationHelper.Validate<GBMPriceIndexModel>(calibrationData, priceModel, output);
            var priceFactor             = CalibrationHelper.GetPriceFactor(priceModel);

            // Moment statistics.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, model, ReturnMethod.Log, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            // Set model parameters using earliest point.
            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for {0}.", priceFactor.GetKey()));
                return;
            }

            model.Drift = momentStatistics[0].GetGbmDrift();
            model.Vol   = momentStatistics[0].Vol;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GBMPriceIndexModel);
        }
    }

    /// <summary>
    /// Calibration of the GBM Price Index Drift Model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("GBM Price Index Drift Calibration")]
    public class GBMPriceIndexDriftCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the GBM Price Index Drift model from the statistics
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            GBMPriceIndexDriftModel model   = CalibrationHelper.Validate<GBMPriceIndexDriftModel>(calibrationData, priceModel, output);
            var priceFactor                 = CalibrationHelper.GetPriceFactor(priceModel);

            // Moment statistics.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, model, ReturnMethod.Log, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            // Set model parameter using earliest available point.
            if(!isValidStatistics[0]) 
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for {0}.", priceFactor.GetKey()));
                return;
            }

            model.Vol = momentStatistics[0].Vol;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GBMPriceIndexDriftModel);
        }
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class GOUAssetPriceStatisticsCalibration
    {
        /// <summary>
        /// Initializes a new instance of the GOUAssetPriceCalibration class; mainly for setting default values.
        /// </summary>
        public GOUAssetPriceStatisticsCalibration()
        {
            // Initialise default values
            Use_Absolute_Theta = YesNo.Yes;
            Use_Mu             = YesNo.No;
        }

        public YesNo Use_Absolute_Theta
        {
            get; set;
        }

        public YesNo Use_Mu
        {
            get; set;
        }

        /// <summary>
        /// Calibrates the GBM Asset Price model with mean reversion from the statistics.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            GOUAssetPriceModel model    = CalibrationHelper.Validate<GOUAssetPriceModel>(calibrationData, priceModel, output);
            IStatistics statistics      = calibrationData.Statistics;
            IPriceFactor priceFactor    = CalibrationHelper.GetPriceFactor(priceModel);

            // Commodities may have multiple prices so we need to pick the earliest
            List<int> indexes;
            List<string> points;

            if (statistics.FindEntries(priceFactor.GetType().Name, model.fID, out indexes, out points))
            {
                int index = indexes[0];
                double drift = StatisticsHelper.GetStatistic(statistics, index,
                                                             BasePlugInStatisticsCalculations.DriftString, ReturnMethod.Log);
                double longRunMean = StatisticsHelper.GetStatistic(statistics, index,
                                                                   BasePlugInStatisticsCalculations.LongRunMeanString, ReturnMethod.Log);
                double alpha = StatisticsHelper.GetStatistic(statistics, index,
                                                             BasePlugInStatisticsCalculations.MeanRevSpeedString, ReturnMethod.Log);
                double sigma = StatisticsHelper.GetStatistic(statistics, index,
                                                             BasePlugInStatisticsCalculations.ReversionVolString, ReturnMethod.Log);

                model.Mu    = Use_Mu == YesNo.Yes ? drift : 0.0;
                model.Alpha = alpha;
                model.Sigma = sigma;

                // Convert long run mean
                if (longRunMean > 0.0)
                {
                    // Let M = longRunMean
                    // Let L = log(M) - 0.25 * sigma * sigma / alpha
                    // When Use_Absolute_Theta = No,  Theta = L - log(S(0))
                    // When Use_Absolute_Theta = Yes, Theta = L
                    double spotPrice    = Use_Absolute_Theta == YesNo.No
                                           ? ((IAssetPrice)(model.GetPriceFactor())).GetQuotedSpotPrice()[0] : 1.0;
                    model.Theta         = Math.Log(longRunMean/spotPrice) - 0.25*sigma*sigma/alpha;
                }
                else
                {
                    // Long run price level cannot be nonpositive.
                    errorList.Add(ErrorLevel.Error, string.Format("Long run mean {0} in statistics for {1} cannot be nonnegative.", longRunMean,
                                      priceFactor.GetKey()));
                    return;
                }

                model.Use_Absolute_Theta = Use_Absolute_Theta;

                GetUpdatedCurves(output);
            }
            else
            {
                errorList.Add(ErrorLevel.Error, string.Format("No statistics for {0}.", priceFactor.GetKey()));
            }
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="output">Revised output.</param>
        private void GetUpdatedCurves(XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }
    }
}