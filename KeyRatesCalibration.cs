using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration of Gaussian Key Rates model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Gaussian Key Rates Interest Rate Calibration")]
    public class GaussianKeyRatesInterestRateCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Distribution type can be normal or lognormal.
        /// </summary>
        public ProbabilityDistribution Distribution_Type { get; set; }

        /// <summary>
        /// Drift method.
        /// </summary>
        public DriftMethod2 Rate_Drift_Model { get; set; }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            GaussianKeyRatesInterestRateModel gaussianKeyRatesModel = CalibrationHelper.Validate<GaussianKeyRatesInterestRateModel>(calibrationData, priceModel, output);

            gaussianKeyRatesModel.Distribution_Type = Distribution_Type;
            gaussianKeyRatesModel.Rate_Drift_Model = Rate_Drift_Model;

            CalibrateGaussianKeyRatesInterestRateModel(calibrationData, gaussianKeyRatesModel, errorList);

            WriteDriftVolCalibrationDataAsXml(gaussianKeyRatesModel, output, Rate_Drift_Model == DriftMethod2.Historical);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GaussianKeyRatesInterestRateModel);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            WriteDriftVolCalibrationDataAsXml((GaussianKeyRatesInterestRateModel)priceModel, output, Rate_Drift_Model == DriftMethod2.Historical);
        }

        /// <summary>
        /// Writes calibration data to an XML stream.
        /// </summary>
        private static void WriteDriftVolCalibrationDataAsXml(GaussianKeyRatesInterestRateModel gaussianKeyRatesModel, XmlWriter writer, bool isHistoricalDrift)
        {
            var graphData = new ResultSeriesList
            {
                fXAxis = { fLabel = "Time (Years)" },
                fYAxis = { fLabel = "Value" }
            };

            graphData.Add(new ResultSeries("Vol", LineStyle.Solid, gaussianKeyRatesModel.Vol.Clone()));

            if (isHistoricalDrift)
                graphData.Add(new ResultSeries("Drift", LineStyle.Solid, gaussianKeyRatesModel.Drift.Clone()));

            graphData.ToXml(writer);
        }

        /// <summary>
        /// Calculate the model parameters from the statistics.
        /// </summary>
        private void CalibrateGaussianKeyRatesInterestRateModel(ICalibrationData calibrationData, GaussianKeyRatesInterestRateModel gaussianKeyRatesModel, ErrorList errorList)
        {
            ReturnMethod returnMethod = gaussianKeyRatesModel.Distribution_Type == ProbabilityDistribution.Lognormal ? ReturnMethod.Log : ReturnMethod.Diff;

            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, gaussianKeyRatesModel, returnMethod, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            var validPoints = new List<double>();
            var validVols = new List<double>();
            var validDrifts = new List<double>();

            for (int index = 0; index < momentStatistics.Length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                // Ignore tiny volatilities.
                double vol = momentStatistics[index].Vol;
                if (vol < CalcUtils.MinVolatility)
                    continue;

                validPoints.Add(momentStatistics[index].ParsePointIDToDouble());
                validVols.Add(vol);

                if (Rate_Drift_Model == DriftMethod2.Historical)
                    validDrifts.Add(momentStatistics[index].GetGbmDrift());
            }

            // Number of valid points.
            int nPoints = validPoints.Count;

            // Convert point strings to key rate tenors; set volatility and drift.
            gaussianKeyRatesModel.Risk_Factor_Tenors.Clear();
            gaussianKeyRatesModel.Vol.Clear();
            gaussianKeyRatesModel.Drift.Clear();
            for (int indexValid = 0; indexValid < nPoints; ++indexValid)
            {
                double tenor = validPoints[indexValid];
                gaussianKeyRatesModel.Risk_Factor_Tenors.Add(tenor);
                gaussianKeyRatesModel.Vol[tenor] = validVols[indexValid];

                if (Rate_Drift_Model == DriftMethod2.Historical)
                    gaussianKeyRatesModel.Drift[tenor] = validDrifts[indexValid];
            }
        }
    }
}