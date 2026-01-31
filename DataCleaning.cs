// <author>
//   Perukrishnen Vytelingum
// </author>
// <summary>
//   Data cleaning.
// </summary>

using System;
using System.ComponentModel;
using System.Globalization;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Method to handle an outlier or an invalid value.
    /// </summary>
    public enum ReplacementMethod
    {
        Remove,
        Replace
    }

    /// <summary>
    /// Direction of the window with respect to data point.
    /// </summary>
    public enum WindowDirection
    {
        Center,
        Left,
        Right
    }

    /// <summary>
    /// Helper class.
    /// </summary>
    public static class DataCleaningHelper
    {
        /// <summary>
        /// Output data cleaning diagnostics.
        /// </summary>
        /// <param name="methodName">Name of the data cleaning method.</param>
        /// <param name="errorList">List of infos, warnings and errors. Assume this can be null.</param>
        /// <param name="isCorrected">Two-dimensional flag identifying any corrections.</param>
        public static void OutputDiagnostics(string methodName, ErrorList errorList, bool[][] isCorrected)
        {
            if (errorList == null)
                return;

            int numTimeSeries   = isCorrected.Length;
            var numOutliers     = GetNumOfCleanedData(isCorrected);

            // Check if any data was cleaned.
            bool hasNothingCleaned = true;
            foreach (int numOutlierInTimeSeries in numOutliers)
            {
                if (numOutlierInTimeSeries > 0)
                {
                    hasNothingCleaned = false;
                    break;
                }
            }

            // Do not output any message if nothing required cleaning.
            if (hasNothingCleaned)
                return;

            string output = string.Format(
                "Historical data required cleaning using the {0} - number of corrections {1}", methodName,
                numTimeSeries == 1 ? "= " : "at each point = (");

            for (int tenorIndex = 0; tenorIndex < numTimeSeries - 1; ++tenorIndex)
            {
                output += numOutliers[tenorIndex] + ",";
            }

            output += numOutliers[numTimeSeries - 1] + (numTimeSeries == 1 ? "." : "). ");

            errorList.Add(ErrorLevel.Info, output);
        }

        /// <summary>
        /// Calculate the mean and sample variance of a window of the timeseries, ignoring the value at [tenorIndex, timeIndex] which is
        /// potentially corrupted. Assume value at time index is not null.
        /// </summary>
        /// <param name="windowSize">Window size.</param>
        /// <param name="windowDirection">Whether window is in the middle, left of right of time index.</param>
        /// <param name="timeIndex">Index of time where statistics required.</param>
        /// <param name="historicalZeroRates">Array of values.</param>
        /// <param name="mean">Return the local mean.</param>
        /// <param name="stdev">Return the local standard deviation.</param>
        public static void CalculateLocalMomentsOfTimeSeries(int windowSize, WindowDirection windowDirection, int timeIndex, double[] historicalZeroRates, out double mean, out double stdev)
        {
            int sampleCount     = 0;
            double sumSample    = 0.0;
            double sumSampleSqr = 0.0;

            if (windowDirection == WindowDirection.Left)
            {
                SumOnTheLeft(timeIndex, historicalZeroRates, windowSize, out sampleCount, out sumSample, out sumSampleSqr);
            }
            else if (windowDirection == WindowDirection.Center)
            {
                int halfWindowSize = windowSize / 2;
                int sampleCountLeft;
                double sumSampleLeft;
                double sumSampleSqrLeft;
                SumOnTheLeft(timeIndex, historicalZeroRates, halfWindowSize, out sampleCountLeft, out sumSampleLeft, out sumSampleSqrLeft);

                int sampleCountRight;
                double sumSampleRight;
                double sumSampleSqrRight;
                SumOnTheRight(timeIndex, historicalZeroRates, halfWindowSize, out sampleCountRight, out sumSampleRight, out sumSampleSqrRight);

                sampleCount     = sampleCountLeft + sampleCountRight;
                sumSample       = sumSampleLeft + sumSampleRight;
                sumSampleSqr    = sumSampleSqrLeft + sumSampleSqrRight;
            }

            // Mean
            mean = sumSample / sampleCount;

            // Unbiased standard deviation
            stdev = Math.Sqrt(sumSampleSqr / sampleCount - mean * mean);
        }

        /// <summary>
        /// Logs diagnostics for the cleansed point.
        /// </summary>
        public static void LogDiagnostics(ErrorList errorList, string pointID, DateTime date, double originalValue, double cleansedValue)
        {
            CalibrationHelper.LogError(errorList, ErrorLevel.Info, ErrorLevel.Info, 
                string.Format(CultureInfo.InvariantCulture, "Point '{0}' on '{1}' changed from '{2}' to '{3}'", pointID, date, originalValue, cleansedValue));
        }

        /// <summary>
        /// Compute statistics on the left of the data point. Ignore any NaN data.
        /// </summary>
        /// <param name="windowSize">Window size.</param>
        /// <param name="timeIndex">Index of time where statistics required.</param>
        /// <param name="historicalZeroRates">Array of values.</param>
        /// <param name="sampleCount">Number of samples.</param>
        /// <param name="sumSample">Sum of samples.</param>
        /// <param name="sumSampleSqr">Sum of square of samples.</param>
        private static void SumOnTheRight(int timeIndex, double[] historicalZeroRates, int windowSize, out int sampleCount, out double sumSample, out double sumSampleSqr)
        {
            sampleCount = 0;
            sumSample = 0.0;
            sumSampleSqr = 0.0;

            for (int index = timeIndex + 1; index < historicalZeroRates.Length; ++index)
            {
                var d = historicalZeroRates[index];

                if (!double.IsNaN(d))
                {
                    sumSample       += d;
                    sumSampleSqr    += d * d;

                    sampleCount++;
                }

                if (sampleCount >= windowSize)
                {
                    return;
                }
            }
        }

        /// <summary>
        /// Compute statistics on the right of the data point. Ignore any NaN data.
        /// </summary>
        /// <param name="windowSize">Window size.</param>
        /// <param name="timeIndex">Index of time where statistics required.</param>
        /// <param name="historicalZeroRates">Array of values.</param>
        /// <param name="sampleCount">Number of samples.</param>
        /// <param name="sumSample">Sum of samples.</param>
        /// <param name="sumSampleSqr">Sum of square of samples.</param>
        private static void SumOnTheLeft(int timeIndex, double[] historicalZeroRates, int windowSize, out int sampleCount, out double sumSample, out double sumSampleSqr)
        {
            sampleCount     = 0;
            sumSample       = 0.0;
            sumSampleSqr    = 0.0;

            for (int index = timeIndex - 1; index >= 0; --index)
            {
                var d = historicalZeroRates[index];

                if (!double.IsNaN(d))
                {
                    sumSample       += d;
                    sumSampleSqr    += d * d;

                    sampleCount++;
                }

                if (sampleCount >= windowSize)
                {
                    return;
                }
            }
        }

        /// <summary>
        /// Return the number of a particular point type in each time series.
        /// </summary>
        /// <param name="isCleaned">Two-dimensional flag of corrections.</param>
        /// <returns>Return an array of the number of corrected data in each array.</returns>
        private static int[] GetNumOfCleanedData(bool[][] isCleaned)
        {
            int numTimeSeries = isCleaned.Length;
            var sum = new int[numTimeSeries];

            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                sum[timeSeriesIndex] = 0;
                for (int dateIndex = 0; dateIndex < isCleaned[timeSeriesIndex].Length; ++dateIndex)
                {
                    if (isCleaned[timeSeriesIndex][dateIndex])
                        sum[timeSeriesIndex]++;
                }
            }

            return sum;
        }
    }

    /// <summary>
    /// Data cleaning method used to cap/floor values, e.g., flooring negative values before a log-normal model calibration.
    /// Extends NestedPresentableObject to appear as nested properties of the method in the Calibration Parameters.
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    [Serializable]
    [DisplayName("Data Cleaning Capping / Flooring")]
    public class DataCleaningCapFloorMethod<T> : NestedPresentableObject, IDataCleaningMethod<T>
    {
        private DoubleNull fCap;
        private DoubleNull fFloor;

        /// <summary>
        /// Initialised with default parameters.
        /// </summary>
        public DataCleaningCapFloorMethod()
        {
            Cap = null;
            Floor = null;
            Replacement_Method  = Replacement_Method;
        }

        /// <summary>
        /// User-defined property for a lower bound on zero-rate across all allTenors.
        /// </summary>
        public DoubleNull Floor
        {
            get
            {
                return fFloor;
            }
            set
            {
                if (value.Exists && double.IsNaN(value.Value))
                    throw new ArgumentException("The Floor must be a valid number.");

                fFloor = value;
            }
        }

        /// <summary>
        /// User-defined property for a upper bound on zero-rate across all Tenors.
        /// </summary>
        public DoubleNull Cap
        {
            get
            {
                return fCap;
            }
            set
            {
                if (value.Exists && double.IsNaN(value.Value))
                    throw new ArgumentException("The Cap must be a valid number.");

                fCap = value;
            }
        }

        /// <summary>
        /// Method used to replace a number based on the specified cap/floor.
        /// </summary>
        public ReplacementMethod Replacement_Method { get; set; }

        /// <inheritdoc />
        public void CleanData(TDate[] dates, T[] timeSeriesIDs, FactorTypeAndID factorTypeAndId, ErrorList errorList, ref double[][] timeSeriesStructure)
        {
            int timeSeriesLength    = dates.Length;
            int numTimeSeries       = timeSeriesStructure.Length;

            var isCleaned = new bool[numTimeSeries][];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                isCleaned[index] = new bool[timeSeriesLength];
            }

            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                for (int timeIndex = 0; timeIndex < timeSeriesLength; ++timeIndex)
                {
                    double rate = timeSeriesStructure[timeSeriesIndex][timeIndex];

                    if (double.IsNaN(rate))
                        continue;

                    if (fFloor.Exists && rate < fFloor.Value)
                    {
                        // Flooring required - remove or replace.
                        SetCleanValue(timeSeriesStructure, isCleaned, timeSeriesIndex, timeIndex, fFloor.Value);
                    }
                    else if (fCap.Exists && rate > fCap.Value)
                    {
                        SetCleanValue(timeSeriesStructure, isCleaned, timeSeriesIndex, timeIndex, fCap.Value);
                    }

                    if (isCleaned[timeSeriesIndex][timeIndex])
                    {
                        DataCleaningHelper.LogDiagnostics(errorList, timeSeriesIDs[timeSeriesIndex].ToString(), dates[timeIndex], rate, timeSeriesStructure[timeSeriesIndex][timeIndex]);
                    }
                }
            }

            // Diagnostics.
            DataCleaningHelper.OutputDiagnostics("Capping / Flooring method", errorList, isCleaned);
        }

        /// <inheritdoc />
        public bool ValidateParameters(ErrorList errorList)
        {
            bool validProperties = true;

            if (Cap.Exists && Floor.Exists && Floor.Value > Cap.Value)
            {
                errorList.Add(ErrorLevel.Error, "Floor should be less than or equal to Cap.");
                validProperties = false;
            }

            if (!Cap.Exists && !Floor.Exists)
            {
                errorList.Add(ErrorLevel.Warning, "Cap and Floor are both undefined. Please set at least one of them with the required value to apply data cleaning.");
                validProperties = false;
            }

            return validProperties;
        }

        private void SetCleanValue(double[][] timeSeriesStructure, bool[][] isCleaned, int timeSeriesIndex, int timeIndex, double value)
        {
            timeSeriesStructure[timeSeriesIndex][timeIndex] = Replacement_Method == ReplacementMethod.Replace ? value : double.NaN;
            isCleaned[timeSeriesIndex][timeIndex] = true;
        }
    }

    /// <summary>
    /// Z-Score data cleaning method - removes unlikely values based on a z-score computed from a local mean and variance.
    /// Extends NestedPresentableObject to appear as nested properties of the method in the Calibration Parameters.
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    [Serializable]
    [DisplayName("Data Cleaning ZScore")]
    public class DataCleaningZScoreMethod<T> : NestedPresentableObject, IDataCleaningMethod<T>
    {
        /// <summary>
        /// Initialised with default parameters.
        /// </summary>
        public DataCleaningZScoreMethod()
        {
            Threshold           = 3.09;                         // Default probability of less that 0.1% that an outlier interest rate occurs.
            Window_Size         = 36;                           // Default window size.
            Perturbation        = 0.0001;                       // 0.0001 roughly corresponds to 1 bp.
            Replacement_Method  = ReplacementMethod.Remove;  // By default, remove outlier.
        }

        /// <summary>
        /// Number of standard deviations used as threshold to identify outliers.
        /// </summary>
        public double Threshold                     { get; set; }

        /// <summary>
        /// Window size over which to compute the local statistics for identifying outliers.
        /// </summary>
        public int Window_Size                      { get; set; }

        /// <summary>
        /// Standard deviation perturbation when identifying outliers.
        /// </summary>
        public double Perturbation                  { get; set; }

        /// <summary>
        /// Method used to replace an outlier.
        /// </summary>
        public ReplacementMethod Replacement_Method { get; set; }

        /// <inheritdoc />
        public void CleanData(TDate[] dates, T[] timeSeriesIDs, FactorTypeAndID factorTypeAndId, ErrorList errorList, ref double[][] timeSeriesStructure)
        {
            int timeSeriesLength    = dates.Length;
            int numTimeSeries       = timeSeriesStructure.Length;

            var isCleaned = new bool[numTimeSeries][];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                isCleaned[index] = new bool[timeSeriesLength];
            }

            if (dates.Length < Window_Size)
            {
                // Do not change market data and no outliers identified.
                if (errorList != null)
                    errorList.Add(ErrorLevel.Warning, string.Format("Z-Score method failed. Window size cannot be larger than length of data set."));
                
                return;
            }

            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                for (int timeIndex = 0; timeIndex < timeSeriesLength; ++timeIndex)
                {
                    double rate = timeSeriesStructure[timeSeriesIndex][timeIndex];

                    if (!double.IsNaN(rate))
                    {
                        // Calculate local mean and stdev
                        double mean, stdev;

                        DataCleaningHelper.CalculateLocalMomentsOfTimeSeries(Window_Size, WindowDirection.Center, timeIndex, timeSeriesStructure[timeSeriesIndex], out mean, out stdev);

                        // Calculate z-score
                        double zscore = stdev + Perturbation <= CalcUtils.TINY ? 0.0 : (rate - mean) / (stdev + Perturbation);

                        if (zscore > Threshold || zscore < -Threshold)
                        {
                            // Outlier identified - replacing with expected value or removing data.
                            timeSeriesStructure[timeSeriesIndex][timeIndex] = Replacement_Method == ReplacementMethod.Replace ? mean : double.NaN;
                            isCleaned[timeSeriesIndex][timeIndex]           = true;

                            DataCleaningHelper.LogDiagnostics(
                                errorList,
                                timeSeriesIDs[timeSeriesIndex].ToString(),
                                dates[timeIndex],
                                rate,
                                timeSeriesStructure[timeSeriesIndex][timeIndex]);
                        }
                    }
                }
            }
            
            // Diagnostics.
            DataCleaningHelper.OutputDiagnostics("Z-Score method", errorList, isCleaned);
        }

        /// <inheritdoc />
        public bool ValidateParameters(ErrorList errorList)
        {
            bool validProperties = true;

            if (Perturbation < 0.0)
            {
                errorList.Add(ErrorLevel.Error, "Z-Score Method Data cleaning - Zscore_Perturbation must be greater or equal to zero.");
                validProperties = false;
            }

            if (Window_Size < 2)
            {
                errorList.Add(ErrorLevel.Error, "Z-Score Method Data cleaning - Window_Size must be greater than 1.");
                validProperties = false;
            }

            if (Threshold <= 0.0)
            {
                errorList.Add(ErrorLevel.Error, "Z-Score Method Data cleaning - Outlier_Threshold must be positive.");
                validProperties = false;
            }
            else if (Threshold <= 1.96)
            {
                errorList.Add(ErrorLevel.Info, "Z-Score Method Data cleaning - A reasonable value for Outlier_Threshold would be greater than 1.96 for a 95% confidence interval.");
            }

            return validProperties;
        }
    }
}
