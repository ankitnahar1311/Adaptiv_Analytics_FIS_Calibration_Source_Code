using System;
using System.Collections.Generic;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Static helper functions for pre-computed statistics.
    /// </summary>
    public static class StatisticsHelper
    {
        /// <summary>
        /// Gets the beta value for the specified index as computed against the specified benchmark
        /// </summary>
        /// <remarks>
        /// This value is computed against the volatilities of, and the correlations between, the index and the benchmark index.
        /// It is not read from the BETA statistic.
        /// </remarks>
        public static double GetBeta(IStatistics statistics, int index, int benchmarkIndex, ReturnMethod method)
        {
            if (statistics == null)
                throw new ArgumentNullException("statistics");

            double vol = GetStatistic(statistics, index, BasePlugInStatisticsCalculations.VolatilityString, method);
            double indexVol = GetStatistic(statistics, benchmarkIndex, BasePlugInStatisticsCalculations.VolatilityString, method);
            double correlation = statistics.Correlation(index, benchmarkIndex);

            return indexVol > 0 ? correlation * vol / indexVol : 0;
        }

        /// <summary>
        /// Get statistic value for given factor index, statistic title and return method.
        /// Throws an exception if result is NaN.
        /// </summary>
        public static double GetStatistic(IStatistics statistics, int index, string statistic, ReturnMethod method)
        {
            double result = statistics.GetStatistic(index, statistic, method);
            if (!double.IsNaN(result))
                return result;

            // Check if value available for other return method
            ReturnMethod otherMethod = method == ReturnMethod.Diff ? ReturnMethod.Log : ReturnMethod.Diff;
            if (statistic != StatisticsData.StatisticBeta && !double.IsNaN(statistics.GetStatistic(index, statistic, otherMethod)))
                throw new CalibrationException(string.Format("{0} requested for factor index {1} and return type {2} but available only for return type {3}.", statistic, index, method, otherMethod));

            throw new CalibrationException(string.Format("No {0} for factor index {1}.", statistic, index));
        }

        /// <summary>
        /// Get statistic value for given factor index, statistic title and return method.
        /// Return false if result is NaN.
        /// </summary>
        public static bool TryGetStatistic(IStatistics statistics, int index, string statistic, ReturnMethod method, out double result)
        {
            result = statistics.GetStatistic(index, statistic, method);
            return !double.IsNaN(result);
        }

        /// <summary>
        /// Get the adjusted drift (with volatility correction for a log-normal distribution) assuming a GBM Model.
        /// </summary>
        public static double GetGbmDrift(IStatistics statistics, int index, ReturnMethod method)
        {
            double drift = statistics.GetStatistic(index, BasePlugInStatisticsCalculations.DriftString, method);

            if (method != ReturnMethod.Log)
            {
                return drift;
            }

            double vol;
            if (TryGetGbmVol(statistics, index, method, out vol))
            {
                // Volatility correction for the log-normal case
                drift += 0.5 * vol * vol;
            }

            return drift;
        }

        /// <summary>
        /// Return the volatility assuming a GBM Model.
        /// </summary>
        public static double GetGbmVol(IStatistics statistics, int index, ReturnMethod method)
        {
            return statistics.GetStatistic(index, BasePlugInStatisticsCalculations.VolatilityString, method);
        }

        /// <summary>
        /// Get statistic value for the adjusted drift.
        /// Return false if result is NaN.
        /// </summary>
        public static bool TryGetGbmDrift(IStatistics statistics, int index, ReturnMethod method, out double drift)
        {
            drift = GetGbmDrift(statistics, index, method);
            return !double.IsNaN(drift);
        }

        /// <summary>
        /// Get statistic value for the volatility.
        /// Return false if result is NaN.
        /// </summary>
        public static bool TryGetGbmVol(IStatistics statistics, int index, ReturnMethod method, out double vol)
        {
            vol = GetGbmVol(statistics, index, method);
            return !double.IsNaN(vol);
        }

        /// <summary>
        /// Get mean reversion statistics, indexes and points strings for this factor type and factor ID after filtering out 
        /// points with zero volatility (less than MinVolatility).
        /// </summary>
        public static void GetMeanReversionStatistics(IStatistics statistics, Type priceFactorType, FactorID factorId, ReturnMethod returnMethod, 
            out List<int> validIndexes, out List<string> validPointIDs, out List<double> reversionVols, out List<double> reversionSpeeds, out List<double> longRunMeans)
        {
            List<int> indexes;
            List<string> pointIDs;
            statistics.FindEntries(priceFactorType.Name, factorId, out indexes, out pointIDs);

            validIndexes    = new List<int>(indexes.Count);
            validPointIDs   = new List<string>(indexes.Count);
            reversionVols   = new List<double>(indexes.Count);
            longRunMeans    = new List<double>(indexes.Count);
            reversionSpeeds = new List<double>(indexes.Count);

            for (int i = 0; i < indexes.Count; ++i)
            {
                int currentIndex = indexes[i];

                double vol = GetStatistic(statistics, currentIndex, BasePlugInStatisticsCalculations.ReversionVolString, returnMethod);
                if (vol < CalcUtils.MinVolatility)
                {
                    continue;
                }

                validIndexes.Add(currentIndex);
                validPointIDs.Add(pointIDs[i]);
                reversionVols.Add(vol);
                reversionSpeeds.Add(GetStatistic(statistics, currentIndex, BasePlugInStatisticsCalculations.MeanRevSpeedString, returnMethod));
                longRunMeans.Add(GetStatistic(statistics, currentIndex, BasePlugInStatisticsCalculations.LongRunMeanString, returnMethod));
            }
        }

        /// <summary>
        /// Returns the point indexes and the points converted to doubles. The points must be single numeric coordinates.
        /// </summary>
        public static void FindEntries(IStatistics statistics, string priceFactorType, FactorID id, out List<int> indexes, out List<double> points)
        {
            List<string> pointStrings;
            statistics.FindEntries(priceFactorType, id, out indexes, out pointStrings);

            points = new List<double>(pointStrings.Count);
            foreach (string pointString in pointStrings)
                points.Add(CalcUtils.ParsePoint(pointString));
        }

        /// <summary>
        /// Calls statistics.GetFactorDetails then convert the points to doubles. The points must be single numeric coordinates.
        /// </summary>
        public static void GetFactorDetails(IStatistics statistics, int index, out string priceFactorType, out FactorID id, out double point)
        {
            string pointString;
            statistics.GetFactorDetails(index, out priceFactorType, out id, out pointString);
            point = CalcUtils.ParsePoint(pointString);
        }
    }
}