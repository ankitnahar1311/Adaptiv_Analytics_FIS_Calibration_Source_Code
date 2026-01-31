/// <author>
/// Andy Hudson, Alastair Wilkins
/// </author>
/// <owner>
/// Nathalie Antunes
/// </owner>
/// <summary>
/// Calibration of Clewlow-Strickland forward price model.
/// </summary>
using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration for Clewlow-Strickland Forward Price model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Clewlow-Strickland Forward Price Calibration")]
    public class CSForwardPriceCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the CS Forward Price model from the statistics.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output,
                                      ErrorList errorList)
        {
            CSForwardPriceModel model   = CalibrationHelper.Validate<CSForwardPriceModel>(calibrationData, priceModel, output);
            IPriceFactor priceFactor    = CalibrationHelper.GetPriceFactor(priceModel);

            // Need to get statistics for the underlying spot price.
            // Get mean-reversion statistics using the Log return.
            MeanReversionStatistics[] meanReversionStatistics;
            TimeSeriesStructure<string> timeSeries = null;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Log, null, null, null, null, false);
            if (!GetMeanReversionStatistics(calibrationData, priceModel, calculationParameters, errorList, ref timeSeries, out meanReversionStatistics, out isValidStatistics))
            {
                return;
            }

            MomentStatistics[] momentStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, priceModel, ReturnMethod.Log, errorList, ref timeSeries, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            if(!isValidStatistics[0]) 
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics at first point for {0}.", priceFactor.GetKey()));
                return;
            }

            // Use first point of forward curve
            model.Drift = momentStatistics[0].GetGbmDrift();
            model.Sigma = meanReversionStatistics[0].Sigma;
            model.Alpha = meanReversionStatistics[0].Alpha;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
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
            return typeof(CSForwardPriceModel);
        }
    }
}