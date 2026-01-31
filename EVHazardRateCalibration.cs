/// <author>
/// Nathalie Antunes, Alastair Wilkins
/// </author>
/// <owner>
/// Nathalie Antunes
/// </owner>
/// <summary>
/// Calibration of exponential Vasicek hazard rate model.
/// </summary>
using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration for Exponential Vasicek Hazard Rate model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Exponential Vasicek Hazard Rate Calibration")]
    public class ExpVasicekHazardRateCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the Exponential Vasicek Hazard Rate Price model from the statistics.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            ExpVasicekHazardRateModel model = CalibrationHelper.Validate<ExpVasicekHazardRateModel>(calibrationData, priceModel, output);
            IPriceFactor priceFactor        = CalibrationHelper.GetPriceFactor(priceModel);

            MeanReversionStatistics[] meanReversionStatistics;
            TimeSeriesStructure<string> timeSeries = null;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Log, null, null, null, null, false);
            if (!GetMeanReversionStatistics(calibrationData, priceModel, calculationParameters, errorList,
                                            ref timeSeries, out meanReversionStatistics, out isValidStatistics))
            {
                return;
            }

            int count       = 0;
            double sumAlpha = 0.0;

            for (int index = 0; index < meanReversionStatistics.Length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                double tenor = meanReversionStatistics[index].ParsePointIDToDouble();

                if (tenor <= 0.0)
                    continue;

                double alpha = meanReversionStatistics[index].Alpha;
                sumAlpha += alpha;
                ++count;

                // First point
                if (count == 1)
                {
                    double sigma = meanReversionStatistics[index].Sigma;
                    double theta = meanReversionStatistics[index].Theta;

                    // First hazard rate is h_1 = I_1 / t_1.
                    // volatility of h_1 = volatility of I_1
                    model.Sigma = sigma;

                    // Convert long run mean of I_1 to long run mean of log(I_1).
                    // Subtract log(t_1) to get long run mean of log(h_1).
                    model.Theta = Math.Log(theta / tenor) - 0.25 * sigma * sigma / alpha;
                }
            }

            // Alpha set to average of mean reversion speed of I_i
            if (count > 0)
            {
                model.Alpha = sumAlpha / count;
            }
            else
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for {0}.", priceFactor.GetKey()));
                return;
            }

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
            return typeof(ExpVasicekHazardRateModel);
        }
    }
}