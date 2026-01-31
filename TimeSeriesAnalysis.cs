//-----------------------------------------------------------------------
// <owner>
// Justin Chan, Perukrishnen Vytelingum, Philip Koop
// </owner>
//-----------------------------------------------------------------------

using System;
using System.Linq;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Static functions to analyse and transform time series.
    /// </summary>
    public static class TimeSeriesAnalysis
    {
        /// <summary>
        /// Detrend the time series using the linear method.
        /// </summary>
        /// <typeparam name="TPoint">Point type.</typeparam>
        public static bool TryDetrendLinear<TPoint>(TimeSeriesStructure<TPoint> timeSeriesStructure, ErrorList errors, IHolidayCalendar calendar, out double[] drifts)
        {
            Requires.Assert(timeSeriesStructure != null);

            bool hasHorizon      = timeSeriesStructure.Horizon.HasValue;
            int numTimeSeries    = timeSeriesStructure.NumTimeSeries;
            int timeSeriesLength = timeSeriesStructure.Dates.Length;
            var firstDate        = timeSeriesStructure.Dates.FirstOrDefault();
            double[] times = TryGetTimeDifferencesInBusinessDays<TPoint>(timeSeriesStructure.Dates, firstDate, errors, calendar);
            double[] horizonTimes = hasHorizon ? TryGetTimeDifferencesInBusinessDays<TPoint>(timeSeriesStructure.DatesAtHorizon, firstDate, errors, calendar) : null;

            drifts = new double[numTimeSeries];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                var values = timeSeriesStructure.Values[timeSeriesIndex];

                // Linear regression of values against times from first observation date
                double b, drift;
                if (!TryDoLinearRegression(values, times, out drift, out b))
                {
                    errors.Add(ErrorLevel.Error, "Linear regression failed");
                    return false;
                }
                
                // Set constant so that mean of detrended values equals mean of original values
                b = -drift * times.Average();

                for (int dateIndex = 0; dateIndex < timeSeriesLength; ++dateIndex)
                {
                    timeSeriesStructure.Values[timeSeriesIndex][dateIndex] -= drift * times[dateIndex] + b;
                    if (hasHorizon)
                        timeSeriesStructure.ValuesAtHorizon[timeSeriesIndex][dateIndex] -= drift * horizonTimes[dateIndex] + b;
                }

                drifts[timeSeriesIndex] = drift;
            }

            return true;
        }

        /// <summary>
        /// Detrend the time series using the exponential method.
        /// </summary>
        /// <typeparam name="TPoint">Point type.</typeparam>
        public static bool TryDetrendExponential<TPoint>(TimeSeriesStructure<TPoint> timeSeriesStructure, ErrorList errors,
                                                         IHolidayCalendar calendar, out double[] drifts)
        {
            Requires.Assert(timeSeriesStructure != null);

            bool hasHorizon      = timeSeriesStructure.Horizon.HasValue;
            int numTimeSeries    = timeSeriesStructure.NumTimeSeries;
            int timeSeriesLength = timeSeriesStructure.Dates.Length;
            var firstDate        = timeSeriesStructure.Dates.FirstOrDefault();
            var times            = TryGetTimeDifferencesInBusinessDays<TPoint>(timeSeriesStructure.Dates, firstDate, errors, calendar);
            double[] horizonTimes = hasHorizon ? TryGetTimeDifferencesInBusinessDays<TPoint>(timeSeriesStructure.DatesAtHorizon, firstDate, errors, calendar) : null;

            drifts = new double[numTimeSeries];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                var values = timeSeriesStructure.Values[timeSeriesIndex];
                var valuesWithoutNaNs = values.Where(value => !double.IsNaN(value)).ToArray();
                if (valuesWithoutNaNs.Length == 0)
                {
                    errors.Add(ErrorLevel.Error, "Detrending Error - All values are NaN");
                    return false;
                }

                var mean = valuesWithoutNaNs.Average();
                var logValues = values.Select(x => Math.Log(x)).ToArray();

                // Linear regression of log values against times from first observation date
                double drift, b;
                if (!TryDoLinearRegression(logValues, times, out drift, out b))
                {
                    errors.Add(ErrorLevel.Error, "Linear regression failed");
                    return false;
                }

                for (int dateIndex = 0; dateIndex < timeSeriesLength; ++dateIndex)
                {
                    values[dateIndex] *= Math.Exp(-drift * times[dateIndex]);

                    if (hasHorizon)
                        timeSeriesStructure.ValuesAtHorizon[timeSeriesIndex][dateIndex] *= Math.Exp(-drift * horizonTimes[dateIndex]);
                }

                // Calculate factor that makes mean of detrended values equals mean of original values
                var detrendedValuesWithoutNaNs = values.Where(value => !double.IsNaN(value)).ToArray();
                if (detrendedValuesWithoutNaNs.Length == 0)
                {
                    errors.Add(ErrorLevel.Error, "Detrending Error - All detrended values are NaN");
                    return false;
                }

                var newMean = detrendedValuesWithoutNaNs.Average();
                var factor = CalcUtils.IsTiny(newMean) ? 1.0 : mean / newMean;

                for (int dateIndex = 0; dateIndex < timeSeriesLength; ++dateIndex)
                {
                    values[dateIndex] *= factor;
                    if (hasHorizon)
                        timeSeriesStructure.ValuesAtHorizon[timeSeriesIndex][dateIndex] *= factor;
                }

                drifts[timeSeriesIndex] = drift;
            }

            return true;
        }

        /// <summary>
        /// Linear Regression of the form Y = a * X + b with ordinary least square.
        /// </summary>
        /// <param name="y">List of values of y.</param>
        /// <param name="x">List of values of x.</param>
        /// <param name="a">Linear coefficient a.</param>
        /// <param name="b">Constant coefficient b.</param>
        /// <remarks>
        /// The ordinary least square solution to linear regression is such that a = cov(X,Y)/var(X) and b = E(Y) - a * E(X).
        /// </remarks>
        public static bool TryDoLinearRegression(double[] y, double[] x, out double a, out double b)
        {
            Requires.Assert(y != null);
            Requires.Assert(x != null);
            Requires.Assert(x.Length == y.Length);

            // Special case: technically if only one point is given, linear regression cannot be performed. In this case, it's 
            // reasonable to default the "slope" to zero, and the "y-intercept" to the only data point we have (i.e., y[0]).
            if (x.Length == 1)
            {
                a = 0.0;
                b = y[0];
                return !double.IsNaN(b);
            }

            int numData  = 0;
            double meanX = 0.0;
            double meanY = 0.0;
            double cov   = 0.0;
            double varX  = 0.0;

            for (int i = 0; i < x.Length; i++)
            {
                if (double.IsNaN(x[i]) || double.IsNaN(y[i]))
                    continue;

                meanX += x[i];
                meanY += y[i];
                cov   += x[i] * y[i];
                varX  += x[i] * x[i];
                ++numData;
            }

            if (numData > 0)
            {
                meanX /= numData;
                meanY /= numData;
                cov    = cov / numData - meanX * meanY;
                varX   = varX / numData - meanX * meanX;
                a = cov / varX;
                b = meanY - a * meanX;
                return true;
            }

            a = double.NaN;
            b = double.NaN;
            return false;
        }

        /// <summary>
        /// Linear Regression of the form Y = a * X + b with ordinary least square.
        /// </summary>
        /// <param name="linearRegressionData">Contains data (i.e., X and Y) on which linear regression is perform.</param>
        /// <param name="a">Linear coefficient a(must have already been assigned to a vector cache).</param>
        /// <param name="b">Constant coefficient b(must have already been assigned to a vector cache).</param>
        /// <remarks>
        /// The ordinary least square solution to linear regression is such that a = cov(X,Y)/var(X) and b = E(Y) - a * E(X).
        /// </remarks>
        public static void DoLinearRegression(VectorCurve linearRegressionData, Vector a, Vector b)
        {
            if (linearRegressionData.Count <= 0)
                throw new CalibrationException(string.Format("X and Y must not be empty."));

            var x = linearRegressionData.X;
            var y = linearRegressionData.Y;

            using (var cache = Vector.CacheLike(y[0]))
            {
                // Special case: technically if only one point is given, linear regression cannot be performed. In this case, it's 
                // reasonable to default the "slope" to zero, and the "y-intercept" to the only data point we have (i.e., y[0]).
                if (x.Count == 1)
                {
                    a.Assign(0.0);
                    b.Assign(y[0]);
                    return;
                }

                int    numData = linearRegressionData.Count;
                double meanX   = 0.0;
                double varX    = 0.0;
                Vector meanY   = cache.GetClear();
                Vector cov     = cache.GetClear();

                for (int i = 0; i < numData; i++)
                {
                    meanX += x[i];
                    meanY.Add(y[i]);
                    cov.Add(x[i] * y[i]);
                    varX  += x[i] * x[i];
                }

                meanX /= numData;
                meanY.MultiplyBy(1.0 / numData);
                cov.Assign(cov / numData - meanX * meanY);
                varX   = varX / numData - meanX * meanX;

                a.Assign(cov / varX);
                b.Assign(meanY - a * meanX);
            }
        }

        /// <summary>
        /// Try to get the time differences (in business days) from a time series.
        /// </summary>
        /// <remarks>If the differences in number of business days cannot be retrieved then absolute date differences are returned.</remarks>
        /// <typeparam name="TPoint">Point type.</typeparam>
        private static double[] TryGetTimeDifferencesInBusinessDays<TPoint>(TDate[] dates, TDate firstDate,
                                                                           ErrorList errors, IHolidayCalendar calendar)
        {
            // If these is a calendar available, then use it.
            if (calendar != null)
            {
                var businessDays = dates.Select(date => IHolidayCalendarExtensions.CountBusinessDays(calendar, DateTime.FromOADate(firstDate),
                                                                                                                  DateTime.FromOADate(date), false, true)).ToArray();
                return businessDays.Select(days => (double)days).ToArray();
            }
            else
            {
                errors.Add(ErrorLevel.Warning, "Failed to obtain a valid calendar during time series analysis. Using date difference instead of business day difference within calibration.");
                return dates.Select(date => CalcUtils.DaysToYears(date - firstDate)).ToArray();
            }
        }
    }
}
