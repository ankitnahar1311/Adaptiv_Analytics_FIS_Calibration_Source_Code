using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Hull-White Historical Volatility Calibration.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Hull-White Historical Volatility Calibration")]
    public class HullWhiteZeroRateCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        private Curve fHistoricalCurve;

        /// <summary>
        /// Initialises calibration parameters.
        /// </summary>
        public HullWhiteZeroRateCalibration()
        {
            Number_Of_Iterations = 100;
        }

        /// <summary>
        /// Maximum of iterations for the Downhill Simplex optimiser.
        /// </summary>
        public int Number_Of_Iterations { get; set; }

        /// <summary>
        /// Calibrates the Hull White model using historical zero rate volatilities.
        /// Using BasePlugInStatisticsCalculations.VolatilityString statistics to build the historical curve.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            HullWhite1FactorInterestRateModel model = CalibrationHelper.Validate<HullWhite1FactorInterestRateModel>(calibrationData, priceModel, output);

            // Reads historical vols from RiskMetrics into a Curve.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, priceModel, ReturnMethod.Diff, errorList,
                                               ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            fHistoricalCurve = new Curve();
            for (int index = 0; index < momentStatistics.Length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                double term = momentStatistics[index].ParsePointIDToDouble();
                fHistoricalCurve[term] = momentStatistics[index].Vol;
            }

            if (fHistoricalCurve.Count == 0)
                throw new CalibrationException("No historical data in historical curve.");

            double sigma;
            double alpha;
            if (!FitHWParamsToZeroRateVols(fHistoricalCurve, out sigma, out alpha))
                errorList.Add(ErrorLevel.Error, string.Format("Calibration of {0} failed.", model));

            model.Alpha = alpha;
            model.SetConstantSigma(sigma);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// The ModelType that this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(HullWhite1FactorInterestRateModel);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for</param>
        /// <param name="output">Revised output</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            HullWhite1FactorInterestRateModel hullWhiteModel = priceModel as HullWhite1FactorInterestRateModel;

            if (hullWhiteModel == null)
                return;

            Curve modelCurve = GetModelCurve(hullWhiteModel.GetConstantSigma(), hullWhiteModel.Alpha, fHistoricalCurve);
            WriteCalibrationDataToXml(priceModel, fHistoricalCurve, modelCurve, output);
        }

        /// <summary>
        /// Gets the volatility for a given period based on am alpha and sigma value.
        /// </summary>
        /// <param name="sigma">Sigma for which volatility should be solved.</param>
        /// <param name="alpha">Alpha for which volatility should be solved.</param>
        /// <param name="time">Time for which volatility should be solved.</param>
        /// <returns>THe volatility point corresponding to the time, alpha and sigma.</returns>
        private static double HWZeroRateVol(double sigma, double alpha, double time)
        {
            if (Math.Abs(alpha) < CalcUtils.MinReversionRate || time < CalcUtils.MinTime)
            {
                return sigma;
            }

            double alphaTime = alpha * time;
            return sigma * (1 - Math.Exp(-alphaTime)) / alphaTime;
        }

        /// <summary>
        /// Gets the calibrated alpha and sigma values.
        /// </summary>
        /// <param name="zeroRateVols">Historical curve which contains periods.</param>
        /// <param name="sigma">A double to which sigma will be asigned.</param>
        /// <param name="alpha">A double to which alpha will be asigned.</param>
        /// <returns>A boolean indicating success.</returns>
        private bool FitHWParamsToZeroRateVols(Curve zeroRateVols, out double sigma, out double alpha)
        {
            // simplex is just a triangle in (sigma, alpha) space
            var simplex = new double[3, 2];
            const double scale = 0.01;

            // first vertex at (scale, scale)
            simplex[0, 0] = scale;
            simplex[0, 1] = scale;

            // second vertex at (2 * scale, scale)
            simplex[1, 0] = 2 * scale;
            simplex[1, 1] = scale;

            // third vertex at (scale, 2 * scale)
            simplex[2, 0] = scale;
            simplex[2, 1] = 2 * scale;

            int valuationCount;
            double minDistance;
            return Minimizer<Curve>.MinimizeDownhillSimplex(HWZeroRateVolCurveDistance, zeroRateVols, simplex, Number_Of_Iterations, out sigma, out alpha, out minDistance, out valuationCount);
        }

        /// <summary>
        /// Curve distance function.
        /// </summary>
        /// <param name="sigma">Sigma for which to calculate the difference.</param>
        /// <param name="alpha">Alpha for which to calculate the difference.</param>
        /// <param name="zeroRateVolCurve">The initial zero rate vol curve.</param>
        /// <param name="curveDistance">The calculated curve distance.</param>
        /// <returns>bool indicating success.</returns>
        /// <exception cref="ArgumentException"></exception>
        private static bool HWZeroRateVolCurveDistance(double sigma, double alpha, Curve zeroRateVolCurve, out double curveDistance)
        {
            Debug.Assert(zeroRateVolCurve != null, "zeroRateVolCurve != null, required parameter 'zeroRateVolCurve' should not be null");

            curveDistance = 0;

            if (sigma < 0 || alpha <= 0)
                return false;

            int points = zeroRateVolCurve.Count;
            for (int i = 0; i < points; i++)
            {
                double time = zeroRateVolCurve.X[i];
                double zeroRateVol = zeroRateVolCurve.Y[i];

                double diff = zeroRateVol - HWZeroRateVol(sigma, alpha, time);
                curveDistance += diff * diff;
            }

            return true;
        }

        /// <summary>
        /// Constructs the model curve using the calibrated alpha and sigma values.
        /// </summary>
        /// <param name="sigma">Sigma obtained from the calibration.</param>
        /// <param name="alpha">Alpha obtained from the calibration.</param>
        /// <param name="historicalCurve">Historical curve which contains periods needed to construct the model curve.</param>
        /// <returns>The model calibrated curve.</returns>
        private static Curve GetModelCurve(double sigma, double alpha, Curve historicalCurve)
        {
            var modelCurve = new Curve();
            for (int i = 0; i < historicalCurve.X.Count; ++i)
                modelCurve[historicalCurve.X[i]] = (float)HWZeroRateVol(sigma, alpha, historicalCurve.X[i]);

            if (modelCurve.Count <= 0)
                throw new CalibrationException("Calibrated model curve has no data.");

            return modelCurve;
        }

        /// <summary>
        /// Writes calibration data to XML so it can be utilised.
        /// </summary>
        private static void WriteCalibrationDataToXml(PriceModel priceModel, Curve historicalCurve, Curve modelCurve, XmlWriter output)
        {
            try
            {
                const float scale = (float)Percentage.OverPercentagePoint;

                var scaledHistoricalCurve = new Curve(historicalCurve);
                scaledHistoricalCurve.MultiplyBy(scale);

                var scaledModelCurve = new Curve(modelCurve);
                scaledModelCurve.MultiplyBy(scale);

                var graphData = new ResultSeriesList
                {
                    fXAxis = { fLabel = "Time (Years)" },
                    fYAxis = { fLabel = "Zero Rate (%) Annual Standard Deviation" }
                };
                graphData.Add(new ResultSeries("Historical", LineStyle.Solid, scaledHistoricalCurve));
                graphData.Add(new ResultSeries("Model Implied", LineStyle.Solid, scaledModelCurve));
                graphData.ToXml(output);

                CalibrationHelper.WriteModelParametersToXml(priceModel, output);
            }
            catch (XmlException e)
            {
                throw new CalibrationException("Failed to write calibration data to XML.", e);
            }
        }
    }

    /// <summary>
    /// Hull-White Interest rate statistics-based calibration.
    /// </summary>
    public class HWInterestRateStatisticsCalibration
    {
        /// <summary>
        /// Calibrates the HW Interest Rate Price model from the statistics.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            HullWhite1FactorInterestRateModel model = CalibrationHelper.Validate<HullWhite1FactorInterestRateModel>(calibrationData, priceModel, output);

            IStatistics statistics = calibrationData.Statistics;
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(priceModel);

            List<int> indexes;
            List<string> points;
            if (statistics.FindEntries(priceFactor.GetType().Name, model.fID, out indexes, out points))
            {
                // AA Price factor Model calibration.doc states that sigma is the reversion vol of the first point
                // (commented out code below) however the Panorama code seems to be taking the average over all points
                //model.Sigma = sd.Statistics.GetStatistic(indexes[0], Statistics.ReversionVol);
                double sumSigma = 0.0;
                double sumAlpha = 0.0;
                foreach (int index in indexes)
                {
                    sumAlpha += StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.MeanRevSpeedString, ReturnMethod.Diff);
                    sumSigma += StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.ReversionVolString, ReturnMethod.Diff);
                }

                model.Alpha = sumAlpha / indexes.Count;
                model.SetConstantSigma(sumSigma / indexes.Count);

                GetUpdatedCurves(output);
            }
            else
            {
                errorList.Add(ErrorLevel.Error, string.Format("No statistics for {0}.", priceFactor.GetKey()));
            }
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model
        /// </summary>
        /// <param name="output">Revised output</param>
        private void GetUpdatedCurves(XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }
    }

    /// <summary>
    /// Hull-White Hazard Rate Calibration.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Hull-White Hazard Rate Calibration")]
    public class HWHazardRateCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the HW Hazard Rate Price model from the statistics.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<HWHazardRateModel>(calibrationData, priceModel, output);

            MeanReversionStatistics[] meanReversionStatistics;
            TimeSeriesStructure<string> timeSeries = null;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Diff, null, null, null, null, false);
            bool isSuccess = GetMeanReversionStatistics(calibrationData, priceModel, calculationParameters, errorList,
                                                        ref timeSeries, out meanReversionStatistics, out isValidStatistics);

            if (!isSuccess)
                return;

            int length = meanReversionStatistics.Length;
            double alphaSum = 0.0;
            int numValidStatistics = 0;

            for (int index = 0; index < length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                double tenor = meanReversionStatistics[index].ParsePointIDToDouble();

                if (tenor <= 0.0)
                    continue;

                if (numValidStatistics == 0)
                {
                    // First hazard rate is h_1 = I_1 / t_1.
                    // Volatility of h_1 = volatility of I_1 / t_1
                    double logSurvVol = meanReversionStatistics[index].Sigma;
                    model.Sigma = logSurvVol / tenor;
                }

                alphaSum += meanReversionStatistics[index].Alpha;
                ++numValidStatistics;
            }

            if (numValidStatistics > 0)
            {
                // Alpha set to average of mean reversion speed of I_i
                model.Alpha = alphaSum / numValidStatistics;
            }
            else
            {
                errorList.Add(ErrorLevel.Error, string.Format("{0} Calibration - no valid statistics available for calibration.", priceModel));
            }

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for</param>
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
            return typeof(HWHazardRateModel);
        }
    }
}