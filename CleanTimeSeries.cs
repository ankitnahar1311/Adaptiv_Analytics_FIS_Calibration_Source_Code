/// <author>
///   Perukrishnen Vytelingum
/// </author>
/// <summary>
///   Creating time series with data cleaning and missing data filling when available.
/// </summary>
using System;
using System.Collections.Generic;
using System.Linq;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Helper class for time series manipulation.
    /// </summary>
    public static class TimeSeriesHelper
    {
        /// <summary>
        /// Clean the historical data by removing any dates at which there is at least one missing or invalid value.
        /// </summary>
        /// <typeparam name="TPoint">Point type.</typeparam>
        public static void Clean<TPoint>(TimeSeriesStructure<TPoint> timeSeriesStructure, bool requirePositive)
        {
            bool hasHorizon      = timeSeriesStructure.Horizon.HasValue;
            int numTimeSeries    = timeSeriesStructure.NumTimeSeries;
            int timeSeriesLength = timeSeriesStructure.Dates.Length;

            var dates          = new List<TDate>(timeSeriesStructure.Dates.Length);
            var datesAtHorizon = hasHorizon ? new List<TDate>(timeSeriesStructure.Dates.Length) : null;

            var values          = new List<List<double>>(numTimeSeries);
            var valuesAtHorizon = hasHorizon ? new List<List<double>>(timeSeriesStructure.NumTimeSeries) : null;
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                values.Add(new List<double>(timeSeriesLength));
                if (hasHorizon)
                    valuesAtHorizon.Add(new List<double>(timeSeriesLength));
            }

            for (int dateIndex = 0; dateIndex < timeSeriesLength; ++dateIndex)
            {
                if (timeSeriesStructure.Values.All(timeSeries => IsValid(timeSeries[dateIndex], requirePositive)) &&
                    (!hasHorizon || timeSeriesStructure.ValuesAtHorizon.All(timeSeries => IsValid(timeSeries[dateIndex], requirePositive))))
                {
                    // All values and horizon values are all valid at this dateIndex                    
                    dates.Add(timeSeriesStructure.Dates[dateIndex]);
                    if (hasHorizon)
                        datesAtHorizon.Add(timeSeriesStructure.DatesAtHorizon[dateIndex]);

                    for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
                    {
                        values[timeSeriesIndex].Add(timeSeriesStructure.Values[timeSeriesIndex][dateIndex]);
                        if (hasHorizon)
                            valuesAtHorizon[timeSeriesIndex].Add(timeSeriesStructure.ValuesAtHorizon[timeSeriesIndex][dateIndex]);
                    }
                }
            }

            timeSeriesStructure.Dates  = dates.ToArray();
            timeSeriesStructure.Values = TimeSeriesHelper.GetDoubleArray(values);
            if (hasHorizon)
            {
                timeSeriesStructure.DatesAtHorizon  = datesAtHorizon.ToArray();
                timeSeriesStructure.ValuesAtHorizon = TimeSeriesHelper.GetDoubleArray(valuesAtHorizon);
            }
        }

        /// <summary>
        /// Converts the list of array of zero-rates across tenors to a multi-dimensional array.
        /// </summary>
        /// <param name="listDouble">List of list of double.</param>
        /// <typeparam name="T">Type of the list considered, typically a double.</typeparam>
        public static T[][] GetDoubleArray<T>(IList<List<T>> listDouble)
        {
            int length = listDouble.Count;
            var values = new T[length][];
            for (int index = 0; index < length; ++index)
            {
                values[index] = listDouble[index].ToArray();
            }

            return values;
        }

        private static bool IsValid(double value, bool requirePositive)
        {
            return !double.IsNaN(value) && (!requirePositive || value > 0.0);
        }
    }

    /// <summary>
    /// Term structure with a number of time series.
    /// </summary>
    /// <typeparam name="T">Type of point ids of market data, e.g., string for asset prices or Period for interest rates.</typeparam>
    public class TimeSeriesStructure<T>
    {
        // Holiday calendar used to identify business days.
        private readonly IHolidayCalendar fHolidayCalendar;

        /// <summary>
        /// Create a clean time series. isOutlier can be set to null if no data cleaning was done prior.
        /// </summary>
        /// <param name="extractedDates">Array of all extracted dates.</param>
        /// <param name="timeSeriesIDs">Array of timeseries IDs.</param>
        /// <param name="extractedValues">Two-dimensional array of all extracted values (NaN where value is missing).</param>
        /// <param name="holidayCalendar">Holiday calendar used to generate clean timeseries.</param>
        /// <param name="dataRetrievalParameters">Data retrieval parameters.</param>
        /// <param name="errorList">Errorlist to output infos, warnings and errors. This can be set to null for no diagnostics.</param>
        public TimeSeriesStructure(TDate[] extractedDates, T[] timeSeriesIDs, double[][] extractedValues, IHolidayCalendar holidayCalendar, 
                                   DataRetrievalParametersBase<T> dataRetrievalParameters, ErrorList errorList)
        {
            Frequency = dataRetrievalParameters.Frequency;
            Horizon = dataRetrievalParameters.GetHorizon();

            fHolidayCalendar    = holidayCalendar;
            TimeSeriesIDs       = timeSeriesIDs;
            ObservedDates       = extractedDates;
            ObservedValues      = extractedValues;
            RatioMissingData    = 0.0;

            // Get number of business days per year from data retrieval parameters.
            // If not specifed then calculate number of business days per year from historical data.
            int? businessDaysInYear = dataRetrievalParameters.GetBusinessDaysInYear();
            if (businessDaysInYear.HasValue)
                NumberOfBusinessDaysInYear = businessDaysInYear.Value;
            else
                NumberOfBusinessDaysInYear = StatisticsCalculator.GetNumBusinessDaysInYear(ObservedDates.First(), ObservedDates.Last(), fHolidayCalendar);

            DeltaT = Horizon.HasValue ? Horizon.Value.ToYears(NumberOfBusinessDaysInYear) : Frequency.ToYears(NumberOfBusinessDaysInYear);
            
            double archiveDays = ObservedDates[ObservedDates.Length - 1] - ObservedDates[0];
            if (Frequency.ToYears(NumberOfBusinessDaysInYear) < 0.75 * CalcUtils.DaysToYears(archiveDays) / (ObservedDates.Length - 1))
                CalibrationHelper.LogError(errorList, ErrorLevel.Warning, dataRetrievalParameters.Diagnostics_Error_Level, string.Format("Frequency of {0} is too low archive with {1} archive dates in {2} days.", Frequency, ObservedDates.Length, archiveDays));

            double[][] fullTermStructure = GetFullStructure();
            RunDataCleaningMethods(timeSeriesIDs, dataRetrievalParameters, ref fullTermStructure,
                "Pre-Calibration Methods for sampled data:", errorList);
            Values = fullTermStructure;

            if (Horizon.HasValue)
            {
                double[][] fullTermStructureAtHorizon = GetFullStructureAtHorizon();
                RunDataCleaningMethods(timeSeriesIDs, dataRetrievalParameters, ref fullTermStructureAtHorizon,
                    "Pre-Calibration Methods for sampled data at horizon:", errorList);
                ValuesAtHorizon = fullTermStructureAtHorizon;
            }
        }

        /// <summary>
        /// User-defined frequency.
        /// </summary>
        public Term Frequency                           { get; set; }

        /// <summary>
        /// Horizon (equals Frequency if not user-defined).
        /// </summary>
        public Term? Horizon                            { get; set; }

        /// <summary>
        /// Horizon year fraction calculated using days / BusinessDaysInYear.
        /// </summary>
        public double DeltaT                            { get; set; }

        /// <summary>
        /// Number of business days in a year given the holiday calendar.
        /// </summary>
        public double NumberOfBusinessDaysInYear        { get; set; }

        /// <summary>
        /// Ratio of missing data.
        /// </summary>
        public double RatioMissingData                  { get; set; }   

        /// <summary>
        /// Returns the number of tenors.
        /// </summary>
        public int NumTimeSeries
        {
            get { return TimeSeriesIDs != null ? TimeSeriesIDs.Length : 0; }
        }

        /// <summary>
        /// Returns time series length.
        /// </summary>
        public int Length
        {
            get { return Dates != null ? Dates.Length : 0; }
        }            

        /// <summary>
        /// List of IDs structure, e.g., tenor.
        /// </summary>
        public T[] TimeSeriesIDs                        { get; set; }

        /// <summary>
        /// Observed time stamps.
        /// </summary>
        public TDate[] ObservedDates                    { get; set; }

        /// <summary>
        /// Observed market values.
        /// </summary>
        public double[][] ObservedValues                { get; set; }

        /// <summary>
        /// All time stamps.
        /// </summary>
        public TDate[] Dates                            { get; set; }

        /// <summary>
        /// All Market Data, with missing data identified.
        /// </summary>
        public double[][] Values                        { get; set; }

        /// <summary>
        /// All return period end dates.
        /// </summary>
        public TDate[] DatesAtHorizon                   { get; set; }

        /// <summary>
        /// All Market Data, with missing data identified.
        /// </summary>
        public double[][] ValuesAtHorizon               { get; set; }

        private void RunDataCleaningMethods(T[] timeseriesIDs, DataRetrievalParametersBase<T> dataRetrievalParameters,
            ref double[][] termStructure, string infoMessage, ErrorList errorList)
        {
            var dataCleaningMethods = dataRetrievalParameters.GetDataCleaningMethods();
            var localErrors = new ErrorList();
            foreach (var dataCleaningMethod in dataCleaningMethods)
                dataCleaningMethod.CleanData(Dates, timeseriesIDs, null, localErrors, ref termStructure);

            if (localErrors.Count == 0)
                return;

            CalibrationHelper.LogError(errorList, ErrorLevel.Info, dataRetrievalParameters.Diagnostics_Error_Level, infoMessage);
            foreach (var error in localErrors)
                CalibrationHelper.LogError(errorList, dataRetrievalParameters.Diagnostics_Error_Level, error);
        }

        /// <summary>
        /// Get a full structure without missing dates, i.e., through all business dates with missing data set as NaN.
        /// </summary>
        private double[][] GetFullStructure()
        {
            var numTimeSeries = ObservedValues.Length;
            var allDates = new List<TDate>();
            var allValues = new List<double>[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                allValues[index] = new List<double>();
            }

            DateTime firstDate = ObservedDates[0];
            DateTime lastDate = ObservedDates[ObservedDates.Length - 1];
            bool frequencyHasDays = Frequency.Days > 0;

            // Calculate first start date
            var startDate = DateAdjuster.AdjustDate(firstDate, fHolidayCalendar, DateAdjustmentMethod.Following);

            // Calculate first horizon date from start date
            TDate horizonDate = Horizon.HasValue ? DateAdjuster.Add(startDate, Horizon.Value, 1, fHolidayCalendar, DateAdjustmentMethod.Modified_Following) : startDate;

            int startDateIndex = 0;
            int numStartDates = 0;
            int numMissingData = 0;
            while (horizonDate <= lastDate)
            {
                allDates.Add(startDate);
                ++numStartDates;

                // Look for startDate in the ObservedDates
                bool foundStart = FindInObservedDates(startDate, ref startDateIndex);

                // Add values or NaN for this start date
                for (int index = 0; index < numTimeSeries; ++index)
                {
                    var startValueToAdd = foundStart ? ObservedValues[index][startDateIndex] : double.NaN;
                    allValues[index].Add(startValueToAdd);
                    if (double.IsNaN(startValueToAdd))
                        ++numMissingData;
                }

                if (Frequency.IsZero())
                    break; // quit having generated a single start date

                // Generate new start date
                if (frequencyHasDays)
                {
                    // Add years, months and weeks and then business days to last start date (adjustment method is not used)
                    startDate = DateAdjuster.Add(startDate, Frequency, 1, fHolidayCalendar, DateAdjustmentMethod.Modified_Following);
                }
                else
                {
                    // Add multiple of years, months and weeks to first date, and then adjust
                    startDate = DateAdjuster.Add(firstDate, Frequency, numStartDates, fHolidayCalendar, DateAdjustmentMethod.Modified_Following);
                }

                // Generate new horizon date from start date
                horizonDate = Horizon.HasValue ? DateAdjuster.Add(startDate, Horizon.Value, 1, fHolidayCalendar, DateAdjustmentMethod.Modified_Following) : startDate;
            }

            // Calculate the ratio of missing data.
            RatioMissingData = (double)numMissingData / (allDates.Count * numTimeSeries);

            Dates = allDates.ToArray();
            return TimeSeriesHelper.GetDoubleArray(allValues);
        }

        /// <summary>
        /// Construct horizon dates and values.
        /// </summary>
        private double[][] GetFullStructureAtHorizon()
        {
            var numTimeSeries = ObservedValues.Length;
            List<TDate> allHorizonDates = null;
            List<double>[] allHorizonValues = null;
            allHorizonDates = new List<TDate>();
            allHorizonValues = new List<double>[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                allHorizonValues[index] = new List<double>();
            }

            int horizonDateIndex = 0;
            foreach (var startDate in Dates)
            {
                var horizonDate = DateAdjuster.Add(startDate, Horizon.Value, 1, fHolidayCalendar, DateAdjustmentMethod.Modified_Following);
                allHorizonDates.Add(horizonDate);

                // Look for horizonDate in the ObservedDates
                var foundHorizon = FindInObservedDates(horizonDate, ref horizonDateIndex);

                // Add values or NaN for this horizon date
                for (int index = 0; index < numTimeSeries; ++index)
                {
                    var horizonValueToAdd = foundHorizon ? ObservedValues[index][horizonDateIndex] : double.NaN;
                    allHorizonValues[index].Add(horizonValueToAdd);
                }
            }

            DatesAtHorizon = allHorizonDates.ToArray();
            return TimeSeriesHelper.GetDoubleArray(allHorizonValues);
        }

        private bool FindInObservedDates(DateTime date, ref int dateIndex)
        {
            for (int i = dateIndex; i < ObservedDates.Length; ++i)
            {
                if (ObservedDates[i] > date)
                    break;

                // advance date index while ObservedDates[i] <= date
                dateIndex = i;

                if (ObservedDates[i] == date)
                    return true;
            }

            return false;
        }
    }
}
