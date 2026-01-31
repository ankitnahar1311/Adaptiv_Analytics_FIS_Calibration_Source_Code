// <summary>
//   Missing data filling.
// </summary>
using System;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Helper class for data filling.
    /// </summary>
    public static class DataFillingHelper
    {
        /// <summary>
        /// Return the first non-null values in a array of doubles.
        /// </summary>
        /// <param name="values">Array of nullable values.</param>
        /// <returns>Return the first non-null value.</returns>
        public static double GetFirstNonNullValue(double[] values)
        {
            for (int dateIndex = 0; dateIndex < values.Length; ++dateIndex)
            {
                var d = values[dateIndex];
                if (!double.IsNaN(d))
                {
                    return d;
                }
            }

            // Return 0.0 if no valid value available. This code is never reached.
            return 0.0;
        }

        /// <summary>
        /// Return the last non-null values in a array of doubles.
        /// </summary>
        /// <param name="values">Array of nullable values.</param>
        /// <returns>Return the last non-null value.</returns>
        public static double GetLastNonNullValue(double[] values)
        {
            for (int dateIndex = values.Length - 1; dateIndex >= 0; --dateIndex)
            {
                var d = values[dateIndex];
                if (!double.IsNaN(d))
                {
                    return d;
                }
            }

            // Return 0.0 if no valid value available. This code is never reached.
            return 0.0;
        }
    }

    /// <summary>
    /// Brownian Bridge method to fill in completely missing rows of data.
    /// This method assumes a standard holiday calendar. 
    /// Extends NestedPresentableObject to appear as nested properties of the method in the Calibration Parameters.
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    [Serializable]
    [DisplayName("Missing Data Filling Correlated Brownian Bridge")]
    public class BrownianBridgeMissingDataFilling<T> : NestedPresentableObject, IDataCleaningMethod<T>
    {
        // Number of business days in a year - used to normalise dates to years.
        private int? fBusinessDaysInYear = CalcUtils.BUSINESS_DAYS_IN_YEAR;
        private double fNumBusinessDaysInYear = CalcUtils.BUSINESS_DAYS_IN_YEAR;
        
        // Annualised length of returns.
        private double fDeltaT = 1.0 / (double)CalcUtils.BUSINESS_DAYS_IN_YEAR;

        /// <summary>
        /// Set the default parameters.
        /// </summary>
        public BrownianBridgeMissingDataFilling()
        {
            Return_Method           = ReturnMethod.Log;
            Random_Seed             = 12345;
            Window_Size             = 36;
            Horizon                 = new Term(1);
            Business_Days_In_Year   = "260";
        }

        /// <summary>
        /// Return method for brownian bridge.
        /// </summary>
        public ReturnMethod Return_Method   { get; set; }

        /// <summary>
        /// Random seed. Leave field blank to not set the seed.
        /// </summary>
        public int Random_Seed              { get; set; }

        /// <summary>
        /// Window size over which to compute the local statistics.
        /// </summary>
        public int Window_Size              { get; set; }

        /// <summary>
        /// Horizon used to compute the local statistics.
        /// </summary>
        public Term Horizon                 { get; set; }

        /// <summary>
        /// Number of businsess days per year.
        /// </summary>
        /// <remarks>
        /// Will be calculated from historical data if empty.
        /// </remarks>
        public string Business_Days_In_Year    
        { 
            get
            {
                return fBusinessDaysInYear.HasValue ? fBusinessDaysInYear.Value.ToString(CultureInfo.InvariantCulture) : string.Empty;
            }
            set
            {
                int d;
                if (!string.IsNullOrEmpty(value) && int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out d) && d > 0)
                    fBusinessDaysInYear = d;
                else
                    fBusinessDaysInYear = null;
            }
        }

        /// <summary>
        ///  Implementation of the correlated brownian bridge method.
        /// </summary>
        /// <param name="dates">Array of dates.</param>
        /// <param name="ids">IDs of time series, e.g., Period, string ids.</param>
        /// <param name="factorTypeAndId">The factor type and id.</param>
        /// <param name="errorList">List of infos, warnings and errors. Assume this can be null.</param>
        /// <param name="timeSeriesStructure">Two-dimensional array of value with potential outliers.</param>
        public void CleanData(TDate[] dates, T[] ids, FactorTypeAndID factorTypeAndId, ErrorList errorList, ref double[][] timeSeriesStructure)
        {
            if (dates.Length < Window_Size)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Warning, "Time series is too short to compute statistics for data filling - no data filling carried out.");

                return;
            }

            int numTimeSeries       = timeSeriesStructure.Length;
            int lengthTimeSeries    = dates.Length;

            InitializeParameters(dates);

            // Find rows with all missing value - correlated sample required then.
            var requiredCorrelatedSample    = new bool[lengthTimeSeries];
            for (int dateIndex = 0; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                // Backfill if all values in row missing.
                bool allMissing = true;
                for (int index = 0; index < numTimeSeries; ++index)
                {
                    if (!double.IsNaN(timeSeriesStructure[index][dateIndex]))
                    {
                        allMissing = false;
                        break;
                    }
                }

                requiredCorrelatedSample[dateIndex] = allMissing;
            }

            // Compute the correlated Wiener process samples.
            SymmetricMatrix correl;
            double[][] correlatedRandomNumbers;
            bool isSuccess = CorrelationCalculationHelper.GetCorrelationMatrix(timeSeriesStructure, out correl);
            if (!isSuccess)
            {
                // Code is unlikely to be reached. Not doing any data filling if it does.
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, "Unable to compute the correlation matrix. Data filling failed.");

                return;
            }

            isSuccess = ModelSimulationHelper.GetCorrelatedRandomNumbers(requiredCorrelatedSample, correl, numTimeSeries, lengthTimeSeries, 
                Random_Seed, out correlatedRandomNumbers);

            if (!isSuccess)
            {
                // Code is unlikely to be reached. Not doing any data filling if it does.
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, "Unable to compute correlated Wiener Process samples. Data filling failed.");

                return;
            }

            var isFilled = FillMissingRows(dates, timeSeriesStructure, correlatedRandomNumbers, errorList);

            // Diagnostics.
            DataCleaningHelper.OutputDiagnostics("Correlated Brownian Bridge Filling method", errorList, isFilled);
        }

        /// <summary>
        /// Initialize data filling parameters.
        /// </summary>
        public void InitializeParameters(TDate[] dates)
        {
            // Get number of business days per year from data filling parameters.
            // If not specifed then calculate number of business days per year from historical data.
            if (fBusinessDaysInYear.HasValue)
                fNumBusinessDaysInYear = fBusinessDaysInYear.Value;
            else
                fNumBusinessDaysInYear = StatisticsCalculator.GetNumBusinessDaysInYear(dates.First(), dates.Last(),
                    new HolidayCalendar());

            fDeltaT = Horizon.ToYears(fNumBusinessDaysInYear);
        }

        /// <inheritdoc />
        public bool ValidateParameters(ErrorList errorList)
        {
            if (Window_Size < 2)
            {
                errorList.Add(ErrorLevel.Error, "Invalid window size - must be at least 2.");
                return false;
            }

            return true;
        }

        /// <summary>
        /// Fill data across individual time series given the correlated Wiener Process samples.
        /// </summary>
        public bool[][] FillMissingRows(TDate[] dates, double[][] timeSeriesStructure, double[][] correlatedRandomNumbers, ErrorList errorList)
        {
            int numTimeSeries = timeSeriesStructure.Length;
            int lengthTimeSeries = dates.Length;

            // Calendar for purposes of computing time duration.
            var timeAnnualised = new double[lengthTimeSeries];
            timeAnnualised[0] = 0.0;
            for (int dateIndex = 1; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                // Assuming there are no missing business days, the number of business days between each consecutive dates is 1.
                timeAnnualised[dateIndex] = timeAnnualised[dateIndex - 1] + 1.0 / fNumBusinessDaysInYear;
            }

            // Flag all filled data.
            var isFilled = new bool[numTimeSeries][];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                isFilled[timeSeriesIndex] = FillData(errorList, correlatedRandomNumbers[timeSeriesIndex], dates, timeAnnualised, ref timeSeriesStructure[timeSeriesIndex]);
            }

            return isFilled;
        }

        /// <summary>
        /// Fill data across individual time series given the correlated Wiener Process samples.
        /// </summary>
        public bool[] FillData(ErrorList errorList, double[] correlatedRandomNumbers, TDate[] dates, double[] timeAnnualised, 
                                ref double[] timeSeries)
        {
            int timeSeriesLength = dates.Length;
            var isFilled = new bool[timeSeriesLength];

            // Assuming a full time grid.
            // Ensure first value is not null to ensure a left clip at the beginning.
            int index = 0;
            FillMissingValueIfNecessary(timeSeries, index, isFilled, DataFillingHelper.GetFirstNonNullValue);

            // Ensure last value is not null to ensure a right clip at the end.
            index = timeSeries.Length - 1;
            FillMissingValueIfNecessary(timeSeries, index, isFilled, DataFillingHelper.GetLastNonNullValue);

            double previous = timeSeries[0];
            int leftClipIndex = 0;
            for (int dateIndex = 1; dateIndex < timeSeriesLength; ++dateIndex)
            {
                double current = timeSeries[dateIndex];

                if (!double.IsNaN(previous) && double.IsNaN(current))
                {
                    // Beginning of bridge.
                    leftClipIndex = dateIndex - 1;
                }
                else if (double.IsNaN(previous) && !double.IsNaN(current))
                {
                    // End of bridge - create bridge here.
                    int rightClipIndex = dateIndex;

                    // Compute local statistics here.
                    MomentStatistics momentStatistics;
                    bool isSuccess = ComputeLocalGBMStatistics(leftClipIndex, rightClipIndex, dates, timeSeries, out momentStatistics);

                    if (!isSuccess)
                    {
                        if (errorList != null)
                            errorList.Add(ErrorLevel.Error, "Unable to compute drift and volatilty of local GBM.");

                        return isFilled;
                    }

                    // Build the brownian bridge between leftClipIndex and rightClipIndex.
                    BuildBrownianBridge(momentStatistics.Drift, momentStatistics.Vol, correlatedRandomNumbers, leftClipIndex, rightClipIndex, 
                                        timeAnnualised, ref timeSeries, ref isFilled);
                }

                previous = current;
            }

            return isFilled;
        }

        /// <summary>
        /// If the value of the time series at the given index is missing attempt to fill it and mark if
        /// it has been filled.
        /// </summary>
        private static void FillMissingValueIfNecessary(double[] timeSeries, int index, bool[] isFilled,
            Func<double[], double> valueRetrieval)
        {
            if (!double.IsNaN(timeSeries[index]))
                return;

            timeSeries[index] = valueRetrieval(timeSeries);
            if (!double.IsNaN(timeSeries[index]))
                isFilled[index] = true;
        }

        /// <summary>
        /// Build a brownian bridge. Note that for data on partially missing rows, we do not fill and the correlated random numbers are NaN.
        /// </summary>
        private void BuildBrownianBridge(double mu, double sigma, double[] correlatedRandomNumbers, int leftClipIndex, int rightClipIndex, 
            double[] timeAnnualised, ref double[] timeSeries, ref bool[] isFilled)
        {
            var t1   = timeAnnualised[leftClipIndex];
            var t2   = timeAnnualised[rightClipIndex];

            var dLeft   = timeSeries[leftClipIndex];
            var dRight  = timeSeries[rightClipIndex];

            if (!double.IsNaN(dLeft) && !double.IsNaN(dRight))
            {
                double xT1 = dLeft;
                double xT2 = dRight;
                double bT = 0.0;
                double w = Return_Method == ReturnMethod.Log ? Math.Log(xT2 / xT1) / sigma - (mu / sigma - 0.5 * sigma) * (t2 - t1) : (xT2 - xT1 - mu * (t2 - t1)) / sigma;
                for (int bridgeIndex = leftClipIndex + 1; bridgeIndex < rightClipIndex; ++bridgeIndex)
                {
                    // If correlated random number is NaN, we have a partially missing row that we are not going to fill.
                    if (double.IsNaN(correlatedRandomNumbers[bridgeIndex]))
                    {
                        timeSeries[bridgeIndex] = double.NaN;
                        continue;
                    }

                    double t                = timeAnnualised[bridgeIndex];
                    double previousT        = timeAnnualised[bridgeIndex - 1];
                    double stochasticValue  = correlatedRandomNumbers[bridgeIndex] * Math.Sqrt((t2 - t) * (t - previousT) / (t2 - previousT));

                    bT = bT * (t2 - t) / (t2 - previousT) + w * (t - previousT) / (t2 - previousT) + stochasticValue;

                    if (Return_Method == ReturnMethod.Log)
                    {
                        if (double.IsNaN(timeSeries[bridgeIndex]))
                        {
                            timeSeries[bridgeIndex] = xT1 * Math.Exp((mu - 0.5 * sigma * sigma) * (t - t1) + sigma * bT);
                            isFilled[bridgeIndex]   = true;
                        }
                    }
                    else
                    {
                        if (double.IsNaN(timeSeries[bridgeIndex]))
                        {
                            timeSeries[bridgeIndex] = xT1 + mu * (t - t1) + sigma * bT;
                            isFilled[bridgeIndex]   = true;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Compute the local drift and volatility.
        /// </summary>
        private bool ComputeLocalGBMStatistics(int leftClipIndex, int rightClipIndex, TDate[] dates, double[] timeSeries, out MomentStatistics momentStatistics)
        {
            int timeSeriesLength = dates.Length;

            int startIndex  = Math.Max(leftClipIndex - Window_Size, 0);
            int endIndex    = Math.Min(rightClipIndex + Window_Size, timeSeriesLength - 1);
            int length      = endIndex - startIndex;

            var localValues             = new double[length];
            var localValuesAtHorizon    = new double[length];
            var localDates              = new TDate[length];
            for (int index = startIndex; index < endIndex; ++index)
            {
                localDates[index - startIndex]              = dates[index];
                localValues[index - startIndex]             = timeSeries[index];
                localValuesAtHorizon[index - startIndex]    = timeSeries[index + 1];
            }

            return MomentStatisticsCalculationHelper.ComputeMomentStatistics(localValues, localValuesAtHorizon, fDeltaT, Return_Method, null, out momentStatistics);
        }
    }
}
