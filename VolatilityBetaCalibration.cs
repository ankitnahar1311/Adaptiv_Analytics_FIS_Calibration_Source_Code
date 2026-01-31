/// <author>
/// Philip Koop
/// </author>
/// <owner>
/// Philip Koop, Mark Sinclair-McGarvie
/// </owner>
/// <summary>
/// Calibration for the Beta Asset Price Volatility model. There are actually two calibrations,
/// the first for a garden variety name volatility that will be driven by some systemic factor
/// via a beta. This calibration is conceptually a combination of the Simple Volatility model
/// calibration and the Multi-GBM Asset Price beta calibration. The second calibration is
/// for the systemic factors themselves; it is implemented using the Simple Volatility model
/// calibration adapted to the Beta Asset Price Volatility model by setting the beta of
/// the factor to itself equal to one.
/// </summary>
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibrate a Beta Asset Price Vol model using a beta to a systemic factor.
    /// </summary>
    public class BetaAssetPriceVolBetaStatisticsCalibration
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public BetaAssetPriceVolBetaStatisticsCalibration()
        {
            Systemic_Factor         = string.Empty;
            Significance_Threshold  = CalcUtils.DefaultMinProcessWeight;
        }

        /// <summary>
        /// Specify the blending method to be set on the model.
        /// </summary>
        public BlendMethod Blend_Method         { get; set; }

        /// <summary>
        /// The Driver used to lookup beta values.
        /// </summary>
        /// <remarks>
        /// If this is not set then the driver specified in the statistics set will be used.
        /// </remarks>
        public string Systemic_Factor           { get; set; }

        /// <summary>
        /// Specify whether the calibration will allow an idiosyncratic factor.
        /// </summary>
        public YesNo Idiosyncratic_Factor       { get; set; }

        /// <summary>
        /// Significant threshold.
        /// </summary>
        public double Significance_Threshold    { get; set; }

        /// <summary>
        /// Calibrates the price model from the statistics files 
        /// using linear regression.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            BetaAssetPriceVolModel model = CalibrationHelper.Validate<BetaAssetPriceVolModel>(calibrationData, priceModel, output);

            Calculate(calibrationData.Statistics, model, errorList);
            GetUpdatedCurves(output);
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

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(BetaAssetPriceVolModel);
        }

        /// <summary>
        /// Calculate the model parameters.
        /// </summary>
        private void Calculate(IStatistics statistics, BetaAssetPriceVolModel model, ErrorList errorList)
        {
            ISpotProcessVol price = (ISpotProcessVol)CalibrationHelper.GetPriceFactor(model);

            string nameType = price.GetType().Name;

            // Find the indexes in the statistics set of the name-specific price points
            List<int> indexes;
            List<string> points;
            statistics.FindEntries(nameType, model.fID, out indexes, out points);

            // Find the points to use for calibration - for now this is just the points closest to ATM
            SelectCandidatePoints(ref indexes, ref points);
            if (indexes.Count == 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No statistics for {0}.", price.GetKey()));
                return;
            }

            // Find the driver points to use for calibration
            List<int> driverIndexes;
            List<string> driverPoints;
            FactorID driver = GetDriver(statistics, model.fID, nameType, indexes, points, out driverIndexes, out driverPoints);

            // Fetch the name-specific point statistics. Missing name-specific statistics are allowed; in that case, we assume that the name asset is not
            // traded and the statistic is not observable; we will then substitute the driver statistic. However, we insist on all or none: if a particular
            // statistic is missing for any point, we ignore it for all points. The scalar statistics required for the model are derived from the point
            // statistics by simple averaging over the closest-to-ATM points. This logic is by analogy with the MultiGBMAssetPrice Beta calibration.
            double beta = GetMeanBeta(statistics, indexes, driverIndexes);
            double driverVol = GetMeanStatistic(statistics, BasePlugInStatisticsCalculations.ReversionVolString, driverIndexes);
            double nameVol;
            if (!TryGetMeanStatistic(statistics, BasePlugInStatisticsCalculations.ReversionVolString, indexes, out nameVol))
                nameVol = beta*driverVol;

            double reversion;
            if (!TryGetMeanStatistic(statistics, BasePlugInStatisticsCalculations.MeanRevSpeedString, indexes, out reversion))
            {
                if (!TryGetMeanStatistic(statistics, BasePlugInStatisticsCalculations.MeanRevSpeedString, driverIndexes, out reversion))
                {
                    errorList.Add(ErrorLevel.Error, string.Format("Missing Reversion Speeds for both name {0} {1} and index {0} {2}",
                                                nameType, model.fID, driver));
                    return;
                }
            }

            Surface longRunMeans;
            if (!TryGetLongRunMeans(statistics, indexes, points, out longRunMeans))
            {
                if (!TryGetLongRunMeans(statistics, driverIndexes, driverPoints, out longRunMeans))
                {
                    errorList.Add(ErrorLevel.Error, string.Format("Missing Long Run Means for both name {0} {1} and index {0} {2}",
                                                nameType, model.fID, driver));
                    return;
                }
            }

            // compute the correlation
            double correlation = ComputeCorrelation(beta, nameVol, driverVol);
            if (correlation > 1.0 || correlation < -1.0)
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                        "Correlation between {0} {1} and {0} {2} retrieved from beta and volatilities is not valid.",
                        nameType, model.fID, driver));
                return;
            }

            // Set model parameters. The following logic DEVIATES from that of the MultGBMAssetPriceModel Beta calibration: if there is a
            // a name-specific volatility, then the model uses that name-specific volatility even if  there is no idiosyncratic factor. The
            // Beta Volatility Model does not allow the idiosyncratic weight to be set independently.
            model.Systemic_Factor = driver.ToCode();
            model.Vol_of_Vol = nameVol;
            model.Reversion_Speed = reversion;
            model.Blend_Method = Blend_Method;
                    model.Systemic_Weight = correlation;
            PopulateHistoricalSurface(model.Historical_Surface, longRunMeans, price);
        }

        /// <summary>
        /// Filters the lists of indices and points to select those that will be used for calibration.
        /// </summary>
        /// <remarks>
        /// Currently the points used for calibration are those points which are closest to ATM.
        /// </remarks>
        /// <param name="indexes">
        /// Inidices into the stats set of the points on the surface.  When this function returns, this contains a list
        /// of just those points which should be used for calibration.
        /// </param>
        /// <param name="points">
        /// The point strings of the points on the surface with statistics in <paramref name="indexes"/>.  When this
        /// function returns this contians the point strings of just those points which should be calibrated.
        /// </param>
        private static void SelectCandidatePoints(ref List<int> indexes, ref List<string> points)
        {
            double closestToATM = double.MaxValue;

            var closestIndexes  = new List<int>();
            var closestPoints   = new List<string>();

            for (int index = 0; index < points.Count; index++)
            {
                string point = points[index];
                string[] parts = point.Split(',');
                Debug.Assert(parts.Length == 2);
                double moneyness = CalcUtils.ParsePoint(parts[0]);
                double closeness = Math.Abs(moneyness - 1.0);

                if (closeness <= closestToATM)
                {
                    if (closeness < closestToATM)
                    {
                        closestToATM = closeness;
                        closestIndexes.Clear();
                        closestPoints.Clear();
                    }

                    closestIndexes.Add(indexes[index]);
                    closestPoints.Add(points[index]);
                }
            }

            indexes = closestIndexes;
            points = closestPoints;
        }

        /// <summary>
        /// Tries to get the mean value of a statistic over the specified points.
        /// </summary>
        /// <param name="statistics">The statistics set.</param>
        /// <param name="statistic">The statistic.</param>
        /// <param name="indexes">The indexes of the points to average over.</param>
        /// <param name="mean">
        /// When the function returns, if all the required statistics were present then this contains the mean value
        /// of the statistics, otherwise <see cref="double.NaN"/>.
        /// </param>
        /// <returns>
        /// <see langword="true" /> if all the required statistics were present; otherwise <see langword="false" />.
        /// </returns>
        private static bool TryGetMeanStatistic(IStatistics statistics, string statistic, List<int> indexes, out double mean)
        {
            double sum = 0;

            foreach (int index in indexes)
            {
                double current;
                if (!StatisticsHelper.TryGetStatistic(statistics, index, statistic, ReturnMethod.Log, out current))
                {
                    mean = double.NaN;
                    return false;
                }

                sum += current;
            }

            mean = sum/indexes.Count;
            return true;
        }

        /// <summary>
        /// Gets the mean value of a statistic over the specified points, and throws an exception if any statistics 
        /// entries are missing.
        /// </summary>
        /// <param name="statistics">The statistics set.</param>
        /// <param name="statistic">The statistic.</param>
        /// <param name="indexes">The indexes of the points to average over.</param>
        /// <returns>
        /// The mean value for the statistic.
        /// </returns>
        /// <exception cref="CalibrationException">One or more statistics entries were missing.</exception>
        private static double GetMeanStatistic(IStatistics statistics, string statistic, List<int> indexes)
        {
            return indexes.Average(index => StatisticsHelper.GetStatistic(statistics, index, statistic, ReturnMethod.Log));
        }

        /// <summary>
        /// Gets the mean beta value over the specified points.
        /// </summary>
        /// <param name="statistics">The statistics set.</param>
        /// <param name="indexes">The indexes of the points to average over.</param>
        /// <param name="driverIndexes">The indexes of the points on the driver surface</param>
        /// <returns>
        /// The mean beta value.
        /// </returns>
        /// <exception cref="CalibrationException">
        /// One of the required statistics was missing.
        /// </exception>
        private static double GetMeanBeta(IStatistics statistics, List<int> indexes, List<int> driverIndexes)
        {
            // Beta is mandatory - GetBeta will throw an exception if it is missing
            return indexes.Select((t, i) => StatisticsHelper.GetBeta(statistics, t, driverIndexes[i], ReturnMethod.Log)).Average();
        }

        /// <summary>
        /// Tries to get the long run means for the specified points.
        /// </summary>
        /// <param name="statistics">The statistics set.</param>
        /// <param name="indexes">The indexes of the candidate points.</param>
        /// <param name="points">The point strings corresponding to <see cref="indexes"/>.</param>
        /// <param name="longRunMeans">
        /// When the function returns, contains a surface containing the long run means, or <see langword="null" /> if any long run mean statistics are missing.
        /// </param>
        /// <returns>
        /// <see langword="true"/> if the surface was computed successfully, otherwise <see langword="false" />
        /// </returns>
        private static bool TryGetLongRunMeans(IStatistics statistics, List<int> indexes, List<string> points, out Surface longRunMeans)
        {
            longRunMeans = new Surface(1);
            for (int idx = 0; idx < indexes.Count; ++idx)
            {
                int index = indexes[idx];
                string point = points[idx];

                string[] parts = point.Split(',');
                Debug.Assert(parts.Length == 2);
                double expiry = CalcUtils.ParsePoint(parts[1]);

                double longRunMean;
                if (!StatisticsHelper.TryGetStatistic(statistics, index, BasePlugInStatisticsCalculations.LongRunMeanString, ReturnMethod.Log, out longRunMean))
                {
                    longRunMeans = null;
                    return false;
                }

                longRunMeans.SetValue(new[] { expiry }, longRunMean);
            }
            return true;
        }

        /// <summary>
        /// Populates a historical surface from the specified long run means and price factor.
        /// </summary>
        /// <param name="historicalSurface">The historical surface.</param>
        /// <param name="longRunMeans">The long run means.</param>
        /// <param name="price">The price.</param>
        private static void PopulateHistoricalSurface(Surface historicalSurface, Surface longRunMeans, ISpotProcessVol price)
        {
            historicalSurface.Clear();
            Surface surface = price.GetSurface();
            for (int expiryIdx = 0; expiryIdx < surface.X.Count; ++expiryIdx)
            {
                double expiry           = surface.X[expiryIdx];
                Surface moneynessSlice  = surface.Y[expiryIdx];
                double atmVol           = moneynessSlice[1.0].Get();
                double longRunMean      = longRunMeans[expiry].Get();
                double scale            = longRunMean / atmVol;

                for (int i = 0; i < moneynessSlice.Count; ++i)
                {
                    historicalSurface.SetValue(new [] { moneynessSlice.Curve.X[i], expiry },
                                               moneynessSlice.Curve.Y[i]*scale);
                }
            }
        }

        /// <summary>
        /// Identifies the driver factor, and the points from that factor that should be used for calibration.
        /// </summary>
        /// <param name="statistics">The statistics set.</param>
        /// <param name="id">The id of the named factor.</param>
        /// <param name="nameType">The factor type of the named factor.</param>
        /// <param name="indexes">
        /// The indexes of the points from the named factor that should be used for calibration.
        /// </param>
        /// <param name="points">
        /// The point strings of the points from the named factor that should be used for calibration.
        /// </param>
        /// <param name="driverIndexes">
        /// When the function returns, this contains the indexes of the points on the driver factor that should be used
        /// for calibration.
        /// </param>
        /// <param name="driverPoints">
        /// When the function returns, this contains the point strings of the points on the driver factor that should
        /// be used for calibration.
        /// </param>
        /// <remarks>
        /// <para>
        /// If the <see cref="Systemic_Factor"/> property is set then this returns the driving factor, and the points on the
        /// driving factor which correspond to the points specified by <paramref name="points"/>.
        /// </para>
        /// <para>
        /// If the <see cref="Systemic_Factor"/> property is not set then this returns the points specified by the INDEX
        /// statistic for the points specified by <paramref name="points"/> and <paramref name="indexes"/>
        /// </para>
        /// </remarks>
        /// <returns>The ID of the driver factor.</returns>
        /// <exception cref="CalibrationException">
        /// Thrown if neither the <see cref="Systemic_Factor"/> property is set, nor the INDEX statistic, or if the INDEX
        /// statistics for the named factor refer to multiple factors, or there is no statistics available for one of
        /// the points on the driver factor.
        /// </exception>
        private FactorID GetDriver(IStatistics statistics, FactorID id, string nameType, List<int> indexes, List<string> points, 
                                   out List<int> driverIndexes, out List<string> driverPoints)
        {
            Debug.Assert(indexes.Count == points.Count);

            string driverType;
            FactorID driver;

            if (!string.IsNullOrEmpty(Systemic_Factor))
            {
                // if the driver property is specified, then we want to use all the corresponding points on the driver surface.
                driver = new FactorID(Systemic_Factor);
                driverType = nameType;
                driverPoints = points;
            }
            else
            {
                driver = null;
                driverType = null;
                driverPoints = new List<string>(points.Count);

                Debug.Assert(indexes.Count > 0);

                // otherwise read all the INDEX values from the statistics
                for (int i = 0; i < indexes.Count; i++)
                {
                    int index = indexes[i];
                    string point = points[i];

                    string equityIndex = statistics.GetEquityIndex(index);
                    if (string.IsNullOrEmpty(equityIndex))
                        throw new CalibrationException(
                            string.Format("INDEX statistic not defined for point {0} on {1}", point, id));

                    string driverPoint;
                    FactorID driverID;
                    StatisticsData.ParseScalarID(equityIndex, out driverType, out driverID, out driverPoint);

                    if (string.IsNullOrEmpty(driverPoint))
                        driverPoint = point;

                    // make sure all the indexes refer to the same surface
                    if (driver == null)
                        driver = driverID;
                    else if (!driver.Equals(driverID))
                        throw new CalibrationException(
                            string.Format("The INDEX statistic for all points on {0} must refer to the same surface", id));

                    driverPoints.Add(driverPoint);
                }
            }

            // lookup the indexes from the points.
            driverIndexes = new List<int>(driverPoints.Count);
            foreach (string point in driverPoints)
            {
                int driverIndex = statistics.GetIndex(driverType, driver, point);

                if (driverIndex == -1)
                    throw new CalibrationException(
                        string.Format("Failed to find statistics for point {0} on driver {1}", point, driver));

                driverIndexes.Add(driverIndex);
            }

            return driver;
        }

        /// <summary>
        /// Computes the systemic weight for the model.
        /// </summary>
        /// <param name="beta">The beta value.</param>
        /// <param name="nameVol">The reversion volatility of the named surface.</param>
        /// <param name="driverVol">The reversion volatility of the driver surface.</param>
        /// <returns>
        /// The systemic weight for the model.
        /// </returns>
        /// <remarks>
        /// If <see cref="Idiosyncratic_Factor"/> is set to <see cref="YesNo.Yes"/> then the systemic weight is 1.0;
        /// otherwise it is computed as <c><paramref name="beta"/> * <paramref name="driverVol"/> / <paramref name="nameVol"/></c>.
        /// </remarks>
        private double ComputeCorrelation(double beta, double nameVol, double driverVol)
        {
            if (Idiosyncratic_Factor != YesNo.Yes)
                return 1.0;

            double correlation = beta * driverVol / nameVol;

            if (Math.Abs(correlation) < Significance_Threshold)
                return 0.0;

            if (Math.Abs(correlation - 1) < Significance_Threshold)
                return 1.0;

            if (Math.Abs(correlation + 1) < Significance_Threshold)
                return -1.0;

            return correlation;
        }
    }

    /// <summary>
    /// Calibrate a Beta Asset Price Vol model as a systemic factor, i.e. without using an input beta.
    /// </summary>
    /// <remarks>
    /// This calibration is equivalent to a simple volatility calibration.
    /// </remarks>
    public class BetaAssetPriceVolSystemicStatisticsCalibration : SimpleAssetPriceVolStatisticsStatisticsCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            BetaAssetPriceVolModel model = (BetaAssetPriceVolModel)priceModel;
            model.Systemic_Factor = model.fID.ToCode();
            model.Systemic_Weight = 1.0;

            base.Calibrate(calibrationData, priceModel, output, errorList);
        }
    }
}