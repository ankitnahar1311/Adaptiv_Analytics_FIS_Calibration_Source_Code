using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Runtime.Serialization;
using System.Xml;
using SunGard.Adaptiv.Analytics.Framework;
using SunGard.Adaptiv.Core;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration can be based on statistics or the MLE method.
    /// </summary>
    public enum CalibrationMethod
    {
        Pre_Computed_Statistics,
        MLE
    }

    /// <summary>
    /// Static helper functions for calibration.
    /// </summary>
    public static class CalibrationHelper
    {
        /// <summary>
        /// Calculate long run mean surface from smile of initial surface and 
        /// long run means of ATM vols in statistics.
        /// </summary>
        /// <typeparam name="TVolSurface">A vol surface class of type VolSurface.</typeparam>
        /// <typeparam name="TPriceFactor">A Price factor class of type PriceFactor.</typeparam>
        /// <typeparam name="TRiskFactor">A risk factor class of type OrnsteinUhlenbeckProcess.</typeparam>
        /// <param name="model">Model's long run mean surface smile is to be fit with model's current vol surface smile.</param>
        /// <param name="useSinglePoint">True for there is only one single point on vol surface.</param>
        /// <param name="singlePoint">The single point on vol surface; only used when useSinglePoint is true and ignored when false.</param>
        /// <remarks>
        /// This method is provided because .ada files derived from
        /// A360 market data archives contain only ATM vols.
        /// The flag useATMShape is useful when there is only one particular timeseries involved in calibration, and the long run mean surface contains
        /// only one single point.
        /// </remarks>
        public static void CalculateLongRunMeanSurface<TVolSurface, TPriceFactor, TRiskFactor>(SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor> model, bool useSinglePoint, double[] singlePoint)
            where TVolSurface : VolSurface, new()
            where TPriceFactor : class, IPriceFactor
            where TRiskFactor : OrnsteinUhlenbeckProcess, new()
        {
            TVolSurface theta    = model.Theta();      // current surface derived from statistics
            TVolSurface newTheta = new TVolSurface(); // replacement surface

            double atmMoneyness = theta.AtmMoneyness();

            TVolSurface initial = new TVolSurface();
            initial.Surface     = model.InitialSurface();

            // Iterate over points on initial surface and construct newTheta on these points
            var iter = new SurfaceIter(initial.Surface);
            var args = new double[initial.Surface.Dimension];

            while (!iter.AtEnd)
            {
                double vol0 = 0.0;
                iter.GetPoint(ref args, ref vol0);
                VolSurface.Point point = initial.PointFromArgs(args);

                if (!useSinglePoint)
                {
                    if (point.Moneyness != atmMoneyness)
                    {
                        // Add interpolated long run mean at corresponding ATM point multiplied by smile factor vol0 / vol0Atm
                        VolSurface.Point atmPoint = new VolSurface.Point() { Expiry = point.Expiry, Tenor = point.Tenor, Moneyness = atmMoneyness };
                        double           volAtm   = theta.GetValue(atmPoint);
                        double           vol0Atm  = initial.GetValue(atmPoint);
                        newTheta.SetValue(point, volAtm * vol0 / vol0Atm);
                    }
                    else
                    {
                        // For ATM points on initial surface, just add interpolated long run mean
                        newTheta.SetValue(point, theta.GetValue(point));
                    }
                }
                else
                {
                    VolSurface.Point singlePointOnVolSurface = theta.PointFromArgs(singlePoint);
                    double           volAtSinglePoint        = theta.GetValue(singlePointOnVolSurface);
                    double           vol0AtSinglePoint       = initial.GetValue(singlePointOnVolSurface);
                    newTheta.SetValue(point, volAtSinglePoint * vol0 / vol0AtSinglePoint);
                }

                iter.Next();
            }

            // Set replacement surface on model
            theta.Surface = newTheta.Surface;
        }

        /// <summary>
        /// Validate arguments of Calibrate: calibrationData, priceModel and output (and not statistics).
        /// Return priceModel cast to ModelClass.
        /// </summary>
        /// <typeparam name="TModelClass">Model being calibrated.</typeparam>
        public static TModelClass ValidateWithoutStatistics<TModelClass>(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output)
            where TModelClass : PriceModel
        {
            if (calibrationData == null)
                throw new ArgumentNullException("calibrationData", "calibrationData is null.");

            if (calibrationData.MarketData == null)
                throw new ArgumentNullException("calibrationData.MarketData", "calibrationData.MarketData is null.");

            if (output == null)
                throw new ArgumentNullException("output", "output is null.");

            if (priceModel == null)
                throw new ArgumentNullException("priceModel", "priceModel is null.");

            TModelClass model = priceModel as TModelClass;
            if (model == null)
                throw new ArgumentException("priceModel", String.Format("priceModel is not of type {0}.", typeof(TModelClass).Name));

            return model;
        }

        /// <summary>
        /// Convert list of T items to list of doubles using the converter.
        /// Remove duplicates and sort.
        /// </summary>
        /// <typeparam name="T">Type parameter.</typeparam>
        public static List<double> ToDoubleList<T>(List<T> list, Converter<T, double> converter)
        {
            if (list == null || list.Count == 0)
                return null;

            List<double> convertedList = list.ConvertAll<double>(converter).Distinct().ToList();
            convertedList.Sort();

            return convertedList;
        }

        /// <summary>
        /// Return true if value is within tolerance of an item in the list.
        /// </summary>
        public static bool IsInList(List<double> list, double value, double tolerance)
        {
            return list.Any(x => Math.Abs(x - value) < tolerance);
        }

        /// <summary>
        /// Get the Holiday Calendar.
        /// </summary>
        public static IHolidayCalendar GetHolidayCalendar(ICalibrationData calibrationData, ErrorList errorList, string calendar)
        {
            if (String.IsNullOrEmpty(calendar))
                return calibrationData.GetArchiveCalendar();

            ICalendarData calendarData = calibrationData.CalendarData;
            if (calendarData == null)
            {
                errorList.Add(ErrorLevel.Warning, String.Format("Failed to load calendar {0} - no calendar data.", calendar));
                return null;
            }

            calendarData.Validate(calendar, errorList);

            return calendarData.GetCalendarUnion(calendar);
        }

        /// <summary>
        /// Write model parameters to xml output stream.
        /// </summary>
        public static void WriteModelParametersToXml(Object priceModel, XmlWriter output)
        {
            output.WriteStartElement("ModelParameters");
            output.WriteValue(Property.GetPropPacked(priceModel));
            output.WriteEndElement();
        }

        /// <summary>
        /// Returns the model's price factor, throws if it is null.
        /// </summary>
        public static IPriceFactor GetPriceFactor(PriceModel priceModel)
        {
            IPriceFactor priceFactor = priceModel.GetPriceFactor();
            if (priceFactor == null)
                throw new AnalyticsException("No price factor attached to model {0}".FormatInvariant(priceModel));
            return priceFactor;
        }

        /// <summary>
        /// Returns the model's price factor's subtype. If the price factor is null it returns "Unknown".
        /// </summary>
        public static string SafeSubType(PriceModel model)
        {
            var pf = model.GetPriceFactor();
            if (pf == null)
                return "Unknown";
            return pf.TypeDotSubType();
        }

        /// <summary>
        /// Validate arguments of Calibrate: calibrationData and output.
        /// </summary>
        public static void Validate(ICalibrationData calibrationData, XmlWriter output)
        {
            if (calibrationData == null)
                throw new ArgumentNullException("calibrationData", "calibrationData is null.");

            if (calibrationData.MarketData == null)
                throw new ArgumentNullException("calibrationData.MarketData", "calibrationData.MarketData is null.");

            if (calibrationData.Statistics == null)
                throw new ArgumentNullException("calibrationData.Statistics", "calibrationData.Statistics is null.");

            if (output == null)
                throw new ArgumentNullException("output", "output is null.");
        }

        /// <summary>
        /// Validate arguments of Calibrate: calibrationData, priceModel and output.
        /// </summary>
        public static void Validate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output)
        {
            Validate(calibrationData, output);

            if (priceModel == null)
                throw new ArgumentNullException("priceModel", "priceModel is null.");
        }

        /// <summary>
        /// Validate arguments of Calibrate: calibrationData, priceModel and output.
        /// </summary>
        /// <typeparam name="TModel">Type of price factor model.</typeparam>
        /// <returns>Returns priceModel cast to TModel.</returns>
        public static TModel Validate<TModel>(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output)
            where TModel : PriceModel
        {
            Validate(calibrationData, priceModel, output);

            TModel model = priceModel as TModel;
            if (model == null)
                throw new ArgumentException("priceModel", String.Format("priceModel is not of type {0}.", typeof(TModel).Name));

            return model;
        }

        /// <summary>
        /// Return true is statistics for all time series are invalid.
        /// </summary>
        public static bool AreAllStatisticsInvalid(IEnumerable<bool> validStatistics)
        {
            if (validStatistics == null)
                return true;

            return validStatistics.All(validStatistic => !validStatistic);
        }

        /// <summary>
        /// Logs the given message if error list is not null and the errorLevel is not <see cref="ErrorLevel.None"/> and 
        /// is greater than minimum error level defined by user.
        /// </summary>
        public static void LogError(ErrorList errorList, ErrorLevel errorLevel, ErrorLevel userDefinedMinErrorLevel, string message)
        {
            LogError(errorList, userDefinedMinErrorLevel, new LogEntryDetails() { Severity = errorLevel.ToSeverity(), Message = message });
        }

        /// <summary>
        /// Logs the given error info if error list is not null and the errorLevel is not <see cref="ErrorLevel.None"/> and 
        /// is greater than minimum error level defined by user.
        /// </summary>
        public static void LogError(ErrorList errorList, ErrorLevel userDefinedMinErrorLevel, LogEntryDetails errorInfo)
        {
            if (userDefinedMinErrorLevel != ErrorLevel.None && errorInfo.Level() >= userDefinedMinErrorLevel)
                errorList.Add(errorInfo);
        }

        /// <summary>
        /// Select a subset of the term structure after removing any tenor where too much data is missing.
        /// </summary>
        /// <param name="allHistoricalTenors">An array of all tenors.</param>
        /// <param name="allHistoricalZeroRates">A two-dimension array of all historical zero-rates.</param>
        /// <param name="minTenorToCalibrate">Minimum tenor to calibrate to. Use minimum tenor available if null.</param>
        /// <param name="thresholdToIgnoreTenor">Relative threshold to ignore a tenor, i.e., ignore any tenor with not enough data. Accept all if null.</param>
        /// <param name="errorList">Warning about tenor ignored.</param>
        /// <param name="unselectedTenors">Unselected tenors.</param>
        /// <param name="invalidSelectedTenors">Tenors that would have been selected, but are invalid because of lack of data.</param>
        /// <param name="validSelectedIndexList">Indices of selected tenors.</param>
        /// <param name="selectedTenors">Output the selected tenors</param>
        /// <param name="selectedHistoricalZeroRates">Output the selected zero-rates.</param>
        public static bool TryGetSelectedTermStructureWithoutInvalidTenors(Period[] allHistoricalTenors, double[][] allHistoricalZeroRates,
                                                                        Period minTenorToCalibrate, double? thresholdToIgnoreTenor, ErrorList errorList,
                                                                        out Period[] unselectedTenors, out Period[] invalidSelectedTenors,
                                                                        out List<int> validSelectedIndexList, out Period[] selectedTenors,
                                                                        out double[][] selectedHistoricalZeroRates)
        {
            List<int> selectedIndices;
            List<int> unselectedIndices;
            GetTenorIndicesToCalibrate(allHistoricalTenors, minTenorToCalibrate, out unselectedIndices, out selectedIndices);

            int selectedLength = selectedIndices.Count;

            if (selectedLength == 0)
            {
               // No valid tenor to calibrate to.
                errorList.Add(ErrorLevel.Error, String.Format("No valid tenor to calibrate to - all tenor(s) below Min Tenor {0}.", minTenorToCalibrate));

                unselectedTenors            = null;
                invalidSelectedTenors       = null;
                validSelectedIndexList      = null;
                selectedTenors              = null;
                selectedHistoricalZeroRates = null;

                return false;
            }

            // Calculate the number of valid data for each tenor within the range.
            var validDataCount = new int[selectedLength];
            for (int index = 0; index < selectedLength; ++index)
            {
                double[] timeseries = allHistoricalZeroRates[selectedIndices[index]];

                // Check if timeseries contains enough data.
                validDataCount[index] = 0;
                foreach (double data in timeseries)
                {
                    validDataCount[index] += Double.IsNaN(data) ? 0 : 1;
                }
            }

            double medianDataCountAcrossTenors = GetMedian(validDataCount);

            // Find valid tenors within selected range that contains enough data.
            validSelectedIndexList = new List<int>();
            var invalidSelectedIndexList = new List<int>();
            for (int index = 0; index < selectedLength; ++index)
            {
                int selectedIndex = selectedIndices[index];

                if (thresholdToIgnoreTenor == null ||
                    validDataCount[index] > thresholdToIgnoreTenor * medianDataCountAcrossTenors)
                {
                    validSelectedIndexList.Add(selectedIndex);
                }
                else
                {
                    invalidSelectedIndexList.Add(selectedIndex);
                }
            }

            if (errorList != null && invalidSelectedIndexList.Count > 0)
            {
                string tenorsNotToCalibrateStr = String.Empty;
                for (int i = 0; i < invalidSelectedIndexList.Count; ++i)
                {
                    int index = invalidSelectedIndexList[i];
                    tenorsNotToCalibrateStr += allHistoricalTenors[index] + (i == invalidSelectedIndexList.Count - 1 ? String.Empty : ", ");
                }

                errorList.Add(ErrorLevel.Warning,
                              String.Format("Ignoring Tenor(s) ({0}) due to insufficient data.", tenorsNotToCalibrateStr));
            }

            selectedTenors = new Period[validSelectedIndexList.Count];
            selectedHistoricalZeroRates = new double[validSelectedIndexList.Count][];
            for (int tenorIndex = 0; tenorIndex < validSelectedIndexList.Count; ++tenorIndex)
            {
                int selectedIndex = validSelectedIndexList[tenorIndex];
                selectedTenors[tenorIndex] = allHistoricalTenors[selectedIndex];
                selectedHistoricalZeroRates[tenorIndex] = allHistoricalZeroRates[selectedIndex];
            }

            unselectedTenors = new Period[unselectedIndices.Count];
            for (int tenorIndex = 0; tenorIndex < unselectedIndices.Count; ++tenorIndex)
            {
                int unselectedIndex = unselectedIndices[tenorIndex];
                unselectedTenors[tenorIndex] = allHistoricalTenors[unselectedIndex];
            }

            invalidSelectedTenors = new Period[invalidSelectedIndexList.Count];
            for (int tenorIndex = 0; tenorIndex < invalidSelectedIndexList.Count; ++tenorIndex)
            {
                int invalidSelectedIndex = invalidSelectedIndexList[tenorIndex];
                invalidSelectedTenors[tenorIndex] = allHistoricalTenors[invalidSelectedIndex];
            }

            return true;
        }

        /// <summary>
        /// Find the index of the single tenor to calibrate to.
        /// </summary>
        public static int GetTenorIndexToCalibrate(Period[] allHistoricalTenors, Period tenorToCalibrate)
        {
            int length = allHistoricalTenors.Length;

            for (int tenorIndex = 0; tenorIndex < length - 1; ++tenorIndex)
            {
                if (Math.Abs(allHistoricalTenors[tenorIndex] - tenorToCalibrate) <= CalcUtils.SMALL)
                {
                    return tenorIndex;
                }

                Period diffRight = (allHistoricalTenors[tenorIndex + 1] - CalcUtils.SMALL) - tenorToCalibrate;
                if (diffRight > 0.0)
                {
                    Period diffLeft = tenorToCalibrate - (allHistoricalTenors[tenorIndex] + CalcUtils.SMALL);
                    return diffLeft <= diffRight ? tenorIndex : tenorIndex + 1;
                }
            }

            return length - 1;
        }

        /// <summary>
        /// Check if Price Index Volatility is required.
        /// </summary>
        public static bool IsPriceIndexVolatilityRequired(IEnumerable<CalibrationInflationInstrumentDefinition> instrumentDefinitions)
        {
            return instrumentDefinitions.Any(definition => definition.Market_Price < CalcUtils.SQRT_DOUBLE_PREC);
        }

        /// <summary>
        /// Returns the median.
        /// </summary>
        private static double GetMedian(int[] values)
        {
            int length = values.Length;
            if (length == 1)
                return 0.0;

            var valueList = new List<double>();
            for (int index = 0; index < length; ++index)
            {
                valueList.Add(values[index]);
            }

            // Sort the list.
            valueList.Sort();

            // Median is calculated as the middle value of the sorted list.
            if (length % 2 == 0)
                return 0.5 * valueList[(int)Math.Floor(length / 2.0)] + 0.5 * valueList[(int)Math.Ceiling(length / 2.0)];

            return valueList[length / 2];
        }

        /// <summary>
        /// Find the indices of the tenors to calibrate and the tenors to extrapolate.
        /// </summary>
        private static void GetTenorIndicesToCalibrate(Period[] allHistoricalTenors, Period minTenorToCalibrate,
                                                       out List<int> unselectedIndices, out List<int> selectedIndices)
        {
            int length = allHistoricalTenors.Length;

            selectedIndices = new List<int>();
            unselectedIndices = new List<int>();
            for (int tenorIndex = 0; tenorIndex < length; ++tenorIndex)
            {
                if (allHistoricalTenors[tenorIndex] < minTenorToCalibrate - CalcUtils.SMALL)
                {
                    unselectedIndices.Add(tenorIndex);
                }
                else
                {
                    selectedIndices.Add(tenorIndex);
                }
            }
        }
    }

    /// <summary>
    /// Helper function for calculation correlation matrices.
    /// </summary>
    public static class CorrelationExtrapolationHelper
    {
        // Maximum number of iterations for Downhill simplex optimiser.
        private const int MaxIterations = 1000;

        // Scale used to generate simplex in downhill algorithm.
        private const double DownhillSimplexScale = 0.01;

        /// <summary>
        /// Flat left extrapolation of correlation matrix.
        /// </summary>
        /// <param name="unselectedLowerTenors">Unselected lower tenors to extrapolate.</param>
        /// <param name="selectedTenors">Selected tenors used to extrapolate.</param>
        /// <param name="matrixToExtrapolate">Covariance/correlation matrix to extrapolate.</param>
        /// <param name="extrapolatedSigmas">Extrapolated sigmas used when dealing with a covariance/correlation matrix transform.</param>
        /// <param name="errorList">List of infos, warnings and errors.</param>
        /// <param name="isCorrelationMatrix">Flag determining whether using a correlation or a covariance matrix.</param>
        /// <param name="extrapolatedMatrix">Extrapolated covariance/correlation matrix.</param>
        /// <returns>Return true if extrapolation is successful.</returns>
        /// <remarks>
        /// The last correlation curve is use to extrapolate lower tenors.
        /// </remarks>
        public static bool ExtrapolateFlatLeftCovarianceMatrix(Period[] unselectedLowerTenors, Period[] selectedTenors, SymmetricMatrix matrixToExtrapolate, double[] extrapolatedSigmas, ErrorList errorList, bool isCorrelationMatrix, out SymmetricMatrix extrapolatedMatrix)
        {
            int numRows = selectedTenors.Length;
            int numRowsToExtrapolate = unselectedLowerTenors.Length;

            if (numRowsToExtrapolate == 0)
            {
                extrapolatedMatrix = matrixToExtrapolate;
                return true;
            }

            if (numRowsToExtrapolate < 0)
            {
                extrapolatedMatrix = matrixToExtrapolate;
                errorList.Add(ErrorLevel.Error, "Extrapolation of correlation matrix was unsuccessful.");
                return false;
            }

            int numCombinedRows = numRows + numRowsToExtrapolate;
            extrapolatedMatrix = SymmetricMatrix.Create(numCombinedRows);
            var allTenors = new double[numCombinedRows];

            // Array of all tenors.
            for (int selectedIndexCol = 0; selectedIndexCol < numRows; ++selectedIndexCol)
                allTenors[numRowsToExtrapolate + selectedIndexCol] = selectedTenors[selectedIndexCol];

            for (int unselectedIndexCol = 0; unselectedIndexCol < numRowsToExtrapolate; ++unselectedIndexCol)
                allTenors[unselectedIndexCol] = unselectedLowerTenors[unselectedIndexCol];

            // Extrapolate from the first selected tenor.
            double alpha;
            double ell;
            var rowOfCorrelations = new double[numRows];
            for (int indexCol = 0; indexCol < numRows; ++indexCol)
            {
                rowOfCorrelations[indexCol] = matrixToExtrapolate[0, indexCol] /
                                              (isCorrelationMatrix
                                                   ? 1.0
                                                   : extrapolatedSigmas[indexCol + numRowsToExtrapolate] * extrapolatedSigmas[indexCol + numRowsToExtrapolate]);
            }

            // Calibrated on correlation matrix (and not covariance matrix)
            bool isSuccess = CalibrateCorrelationFunction(selectedTenors, rowOfCorrelations, 0, out alpha, out ell);

            if (!isSuccess)
            {
                extrapolatedMatrix = null;
                errorList.Add(ErrorLevel.Error, "Extrapolation of correlation matrix was unsuccessful.");
                return false;
            }

            errorList.Add(ErrorLevel.Info, string.Format(CultureInfo.InvariantCulture, "Extrapolation of covariance matrix was successful: a = {0} and b = {1}", alpha, ell));

            // Fill in existing correlations/covariances.
            for (int indexRow = 0; indexRow < numRows; ++indexRow)
            {
                for (int indexCol = 0; indexCol < numRows; ++indexCol)
                {
                    extrapolatedMatrix[indexRow + numRowsToExtrapolate, indexCol + numRowsToExtrapolate] = matrixToExtrapolate[indexRow, indexCol];
                }
            }

            // Fill in extrapolated correlations/covariances (using flat-left interpolation of sigmas for covariances).
            for (int indexRow = 0; indexRow < numRowsToExtrapolate; ++indexRow)
            {
                for (int indexCol = 0; indexCol < numCombinedRows; ++indexCol)
                {
                    double d = allTenors[indexRow] - allTenors[indexCol];
                    extrapolatedMatrix[indexRow, indexCol] = RqCovarianceFunction(d, alpha, ell) * 
                                                                        (isCorrelationMatrix
                                                                             ? 1.0
                                                                             : extrapolatedSigmas[indexRow] * extrapolatedSigmas[indexCol]);
                }
            }

            // Fill in symmetric correlations/covariances.
            for (int indexCol = 0; indexCol < numRowsToExtrapolate; ++indexCol)
            {
                for (int indexRow = numRowsToExtrapolate; indexRow < numCombinedRows; ++indexRow)
                {
                    extrapolatedMatrix[indexRow, indexCol] = extrapolatedMatrix[indexCol, indexRow];
                }
            }

            return true;
        }

        /// <summary>
        /// Calibrate a Rational Quadratic correlation function.
        /// </summary>
        /// <param name="tenors">Array of tenors.</param>
        /// <param name="correlations">Covariance values of function.</param>
        /// <param name="currentTenorIndex">Index of tenor used to calibrate correlation function.</param>
        /// <param name="alpha">Return parameter alpha of correlation function.</param>
        /// <param name="ell">Return parameter ell of correlation function.</param>
        /// <returns>Return true of optimisation is successful.</returns>
        private static bool CalibrateCorrelationFunction(Period[] tenors, double[] correlations, int currentTenorIndex, out double alpha, out double ell)
        {
            // simplex is just a triangle in (alpha, ell) space
            var simplex = new double[3, 2];

            simplex[0, 0] = DownhillSimplexScale;
            simplex[0, 1] = DownhillSimplexScale;
            simplex[1, 0] = 2 * DownhillSimplexScale;
            simplex[1, 1] = DownhillSimplexScale;
            simplex[2, 0] = DownhillSimplexScale;
            simplex[2, 1] = 2 * DownhillSimplexScale;

            int valuationCount;
            double minDistance;
            var info = new CovarianceFunctionInfo
            {
                Xvalues = tenors,
                Curve = correlations,
                CurrentXIndex = currentTenorIndex
            };
            return Minimizer<CovarianceFunctionInfo>.MinimizeDownhillSimplex(CorrelationDistanceError, info, simplex, MaxIterations, out alpha, out ell, out minDistance, out valuationCount);
        }

        /// <summary>
        /// Objective function to be minimized to calibrate the correlation function.
        /// </summary>
        /// <param name="alpha">Parameter alpha of correlation function.</param>
        /// <param name="ell">Parameter ell of correlation function.</param>
        /// <param name="info">Parameters required to compute correlation function.</param>
        /// <param name="error">Return the difference between observed and correlation function value.</param>
        /// <returns>Return true if the difference is valid.</returns>
        private static bool CorrelationDistanceError(double alpha, double ell, CovarianceFunctionInfo info, out double error)
        {
            if (alpha <= CalcUtils.SMALL || ell <= CalcUtils.SMALL)
            {
                error = Double.MaxValue;
                return false;
            }

            int length = info.Length;
            error = 0.0;
            for (int index = 0; index < length; ++index)
            {
                double d = info.Xvalues[index] - info.CurrentXValue;
                double correlPredicted = RqCovarianceFunction(d, alpha, ell);
                double difference = correlPredicted - info.Curve[index];
                error += difference * difference;
            }

            return true;
        }

        /// <summary>
        /// Correlation function.
        /// </summary>
        /// <param name="d">Parameter d.</param>
        /// <param name="alpha">Parameter alpha.</param>
        /// <param name="ell">Parameter ell.</param>
        /// <returns>The function value given those parameters.</returns>
        private static double RqCovarianceFunction(double d, double alpha, double ell)
        {
            return Math.Pow(1 + 0.5 * d * d / (alpha * ell * ell), -alpha);
        }

        /// <summary>
        /// Information used to calibrate the correlation function.
        /// </summary>
        private class CovarianceFunctionInfo
        {
            /// <summary>
            /// Correlation function.
            /// </summary>
            public double[] Curve { get; set; }

            /// <summary>
            /// List of all tenors.
            /// </summary>
            public Period[] Xvalues { get; set; }

            /// <summary>
            /// Index of tenor to calibrate to.
            /// </summary>
            public int CurrentXIndex { get; set; }

            /// <summary>
            /// Return the length of the correlation function. 
            /// </summary>
            public int Length
            {
                get { return Curve.Length; }
            }

            /// <summary>
            /// Return the value of tenor used to calibrate.
            /// </summary>
            public double CurrentXValue
            {
                get { return Xvalues[CurrentXIndex]; }
            }
        }
    }

    /// <summary>
    /// Calibration exceptions.
    /// </summary>
    public class CalibrationException : KnownException
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="CalibrationException"/> class.
        /// </summary>
        public CalibrationException()
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CalibrationException"/> class.
        /// </summary>
        public CalibrationException(string message)
            : base(message)
        {
            Debug.Assert(!string.IsNullOrEmpty(message), "!string.IsNullOrEmpty(message), null or empty string passed to exception.");
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CalibrationException"/> class.
        /// </summary>
        public CalibrationException(string message, Exception innerException)
            : base(message, innerException)
        {
            Debug.Assert(!string.IsNullOrEmpty(message), "!string.IsNullOrEmpty(message), null or empty string passed to exception.");
            Debug.Assert(innerException != null, "innerException != null, null argument passed to exception.");
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CalibrationException"/> class.
        /// </summary>
        protected CalibrationException(SerializationInfo info, StreamingContext context)
            : base(info, context)
        {
        }
    }    
    
    /// <summary>
    /// The purpose of this helper class is to find the most appropriate of the available
    /// forward rates to use as the proxy for a specified forward rate. Forward rates are
    /// assumed to be represented by term (time to start) plus tenor measured in years.
    /// </summary>
    public class ForwardRateDictionary
    {
        public const char PointSeparator = ',';

        protected SortedList<double, TermIndexList> fTenorList = new SortedList<double, TermIndexList>();

        /// <summary>
        /// Size of tenor list.
        /// </summary>
        public int Count
        {
            get { return fTenorList.Count; }
        }

        /// <summary>
        /// Clear the dictionary.
        /// </summary>
        public void Clear()
        {
            fTenorList.Clear();
        }

        /// <summary>
        /// Populate the dictionary with available forward rates.
        /// </summary>
        /// <param name="indexes">List of risk factor indexes</param>
        /// <param name="points">List of datapoint strings corresonding to the indexes</param>
        public void SetEntries(List<int> indexes, List<string> points)
        {
            Debug.Assert(indexes.Count == points.Count, "# of indexes != # of points");

            Clear();
            double term, tenor;
            for (int i = 0; i < indexes.Count; i++)
            {
                ParsePoint(points[i], out term, out tenor);
                TermIndexList indexList;
                if (!fTenorList.TryGetValue(tenor, out indexList))
                {
                    indexList = new TermIndexList();
                    fTenorList.Add(tenor, indexList);
                }

                TermIndex index = new TermIndex();
                index.Term = term;
                index.Index = indexes[i];
                indexList.Add(index);
            }

            foreach (TermIndexList indexList in fTenorList.Values)
                indexList.Sort();
        }

        /// <summary>
        /// Find the risk factor index of the available forward rate that is closest
        /// to the requested term and tenor. Rates are implicitly ordered primarily by
        /// tenor and secondarily by term, so that a rate with a closer tenor is considered
        /// a closer match than a rate with a closer term.
        /// </summary>
        public int GetIndex(double term, double tenor)
        {
            TermIndexList list = ClosestTermList(tenor);
            return list.ClosestIndex(term);
        }

        /// <summary>
        /// Extract the term and tenor of the forward rate from the datapoint string.
        /// </summary>
        public void ParsePoint(string point, out double term, out double tenor)
        {
            StringFragment termString;
            StringFragment tenorString;
            point.Partition(PointSeparator, out termString, out tenorString);
            term = CalcUtils.ParsePoint(termString.ToString());
            tenor = CalcUtils.ParsePoint(tenorString.ToString());
        }

        /// <summary>
        /// Find the list of available terms for the available tenor that is closest to
        /// the given tenor
        /// </summary>
        /// <param name="tenor">The tenor requested</param>
        /// <returns>List of available terms with risk indexes</returns>
        protected TermIndexList ClosestTermList(double tenor)
        {
            // No points
            if (Count == 0)
                return null;

            // At or after the last point
            else if (tenor >= fTenorList.Keys[Count - 1])
                return fTenorList.Values[Count - 1];

            // Before the first point
            else if (tenor < fTenorList.Keys[0])
                return fTenorList.Values[0];

            // Find the list for highest tenor <= that requested
            // using binary search on Keys[0] <= tenor < Keys[Count-1]
            else
            {
                int ll = 0;
                int hh = Count - 1;
                int mm = 0;
                do
                {
                    mm = (ll + hh) >> 1;
                    if (fTenorList.Keys[mm] <= tenor) 
                        ll = mm;
                    else 
                        hh = mm;
                } 
                while (ll < hh - 1);

                if (ll < Count - 1)
                    return (tenor - fTenorList.Keys[ll] < fTenorList.Keys[ll + 1] - tenor) ? fTenorList.Values[ll] : fTenorList.Values[ll + 1];
                else
                    return fTenorList.Values[ll];
            }
        }

        /// <summary>
        /// This struct associates a term with a risk factor index.
        /// </summary>
        public struct TermIndex : IComparable<TermIndex>
        {
            /// <summary>
            /// Term string.
            /// </summary>
            public double Term
            {
                get;
                set;
            }

            /// <summary>
            /// Risk factor index.
            /// </summary>
            public int Index
            {
                get;
                set;
            }

            /// <summary>
            /// Compare term strings.
            /// </summary>
            public int CompareTo(TermIndex other)
            {
                return Term.CompareTo(other.Term);
            }
        }

        /// <summary>
        /// This class is a container for all terms available for rates of a particular tenor.
        /// </summary>
        public class TermIndexList : List<TermIndex>
        {
            /// <summary>
            /// Finds the risk index of the available term that is closest to the given term.
            /// </summary>
            /// <param name="term">The term to look for</param>
            /// <returns>A risk factor index, or -1 on failure</returns>
            public int ClosestIndex(double term)
            {
                // No points
                if (Count == 0)
                    return -1;

                // Look for matching term with binary search
                TermIndex termIndex = new TermIndex();
                termIndex.Term = term;
                int itemIdx = BinarySearch(termIndex);
                if (itemIdx < 0)
                {
                    // The bitwise complement of the index is the next-highest term position
                    itemIdx = ~itemIdx;
                    if (itemIdx >= Count)
                        itemIdx = Count - 1;
                    else if (itemIdx == 0)
                        itemIdx = 0;
                    else if (this[itemIdx].Term - term > term - this[itemIdx - 1].Term)
                        itemIdx = itemIdx - 1;
                }

                return this[itemIdx].Index;
            }
        }
    }
}