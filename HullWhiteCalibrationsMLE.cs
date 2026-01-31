// <summary>
// Maximum Likelihood Estimation (MLE) calibration for a Hull-White Model Zero-Rate Model:
// dr(t) = (theta(t) - alpha * r(t))dt + sigma * dWt.
// </summary>

using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Maximum Likelihood Estimation (MLE) calibration for a Hull-White Zero-Rate Model:
    ///     dr(t) = (theta(t) - alpha * r(t))dt + sigma * dWt
    /// where r(t) is the zero-rate process, alpha is the zero-rate mean reversion speed, sigma is the 
    /// zero-rate volatility, theta(t) is a deterministic curve based on the instantaneous forward rate.
    /// </summary>
    [DisplayName("Hull-White Interest Rate Calibration")]
    public class HWInterestRateCalibration : IModelCalibration
    {
        // Scale used to generate simplex in downhill algorithm.
        private const double DownhillSimplexScale = 0.01;

        // Number of business days in years based on calendars.
        private double fNumBusinessDaysInYear;

        // The drift curve based on the calibrated alpha and sigma.
        private double[] fFirstCalibratedDrift;

        // Annualized times.
        private double[] fTimeAnnualised;

        // Initial yield curve.
        private DoubleCurve fInitialYieldCurve;

        // Array of all tenors.
        private Period[] fHistoricalTenors;

        // Array of selected tenors.
        private Period fSelectedTenorToCalibrateTo;

        // Array of clean zero-rates calibrating to.
        private double[] fSelectedZeroRates;

        // Term structure.
        private TimeSeriesStructure<Period> fTermStructure;

        /// <summary>
        /// Initializes a new instance and sets the defaults.
        /// </summary>
        public HWInterestRateCalibration()
        {
            // Initialise default values
            Calibration_Method = CalibrationMethod.MLE;
            MLE_Parameters     = new HullWhiteInterestRateMLEParameters();
        }

        /// <summary>
        /// Calibration method can be statistics or MLE.
        /// </summary>
        public CalibrationMethod Calibration_Method { get; set; }

        /// <summary>
        /// MLE parameters.
        /// </summary>
        public HullWhiteInterestRateMLEParameters MLE_Parameters { get; set; }

        /// <summary>
        /// Calibrates the Hull-White Short-Rate model using the Maximum Likelihood Estimation method.
        /// </summary>
        /// <param name="calibrationData">Market data and statistics.</param>
        /// <param name="priceModel">The price model to calibrate.</param>
        /// <param name="output">The XmlWriter to which any intermediate data should be written.</param>
        /// <param name="errorList">List of calibration errors and warnings.</param>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            Debug.Assert(errorList != null, "Errorlist cannot be null.");

            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new HWInterestRateStatisticsCalibration();
                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            // Validate properties. Fail if any property not properly set.
            if (!ValidateProperties(errorList, priceModel.ToString()))
            {
                return;
            }

            // Validate the data against the model.
            HullWhite1FactorInterestRateModel hwModel = CalibrationHelper.ValidateWithoutStatistics<HullWhite1FactorInterestRateModel>(calibrationData, priceModel, output);

            // Get the price factor for the model: we need the historical values corresponding to this factor.
            IPriceFactor priceFactor = hwModel.GetPriceFactor();
            if (priceFactor == null)
            {
                errorList.Add(ErrorLevel.Error, "The Hull-White Zero-Rate MLE calibration requires the model's price factor to be present in the market data.");
                return;
            }

            // Term structure.
            if (!GenericDataRetrieval.TryGetHistoricalRateTermStructure(calibrationData, priceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out fTermStructure))
            {
                return;
            }

            // Set the number of business days calculated from holiday calendar.
            fNumBusinessDaysInYear = fTermStructure.NumberOfBusinessDaysInYear;

            // TENORS.
            fHistoricalTenors = fTermStructure.TimeSeriesIDs;

            if (fHistoricalTenors == null || fHistoricalTenors.Length == 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Invalid term points - the {0} could not be calibrated.", priceModel));
                return;
            }

            // Get clean data without any NaN.
            TimeSeriesHelper.Clean(fTermStructure, false);

            if (fTermStructure.Length < StatisticsCalculationHelper.MinTimeSeriesLength)
            {
                errorList.Add(ErrorLevel.Error, "No valid data available. Calibration failed.");
                return;
            }

            // CLEAN TIMESERIES.
            int tenorIndexToCalibrate = CalibrationHelper.GetTenorIndexToCalibrate(fHistoricalTenors, MLE_Parameters.Tenor);
            fSelectedTenorToCalibrateTo = fHistoricalTenors[tenorIndexToCalibrate];
            fSelectedZeroRates = fTermStructure.Values[tenorIndexToCalibrate];

            errorList.Add(ErrorLevel.Info, string.Format("Calibrating {0} to Tenor {1}.", priceModel, fSelectedTenorToCalibrateTo));

            // INITIAL YIELD.
            fInitialYieldCurve = GetInitialYieldCurve(fTermStructure.Length / fNumBusinessDaysInYear, fHistoricalTenors, fTermStructure.Values);

            // MLE OPTIMIZATION.
            double initialAlphaEstimate;
            double initialSigmaEstimate;
            IHolidayCalendar holidayCalendar = MLE_Parameters.Data_Retrieval_Parameters.GetHolidayCalendar(calibrationData, errorList);
            GetInitialEstimate(fSelectedZeroRates, errorList, out initialAlphaEstimate, out initialSigmaEstimate);

            // Pre-compute delta times and annualised times to speed up computation.
            var deltaTimeAnnualised = new double[fTermStructure.Length];
            fTimeAnnualised = new double[fTermStructure.Length];
            for (int dateIndex = 0; dateIndex < fTermStructure.Length; ++dateIndex)
            {
                if (dateIndex > 0)
                {
                    deltaTimeAnnualised[dateIndex] = (holidayCalendar == null ? 1.0 : holidayCalendar.CountBusinessDays(DateTime.FromOADate(fTermStructure.Dates[dateIndex - 1]),
                       DateTime.FromOADate(fTermStructure.Dates[dateIndex]), false, true)) / fNumBusinessDaysInYear;
                }

                fTimeAnnualised[dateIndex] = (holidayCalendar == null ? dateIndex : holidayCalendar.CountBusinessDays(DateTime.FromOADate(fTermStructure.Dates[0]),
                    DateTime.FromOADate(fTermStructure.Dates[dateIndex]), false, true)) / fNumBusinessDaysInYear;
            }

            // Pre-compute the unoptimised part of the drift.
            int lengthTimeSeries = fTimeAnnualised.Length;
            double tau = GetAnnualisedTau(fSelectedTenorToCalibrateTo);

            var unoptimisedDrift = new double[lengthTimeSeries];
            for (int dateIndex = 0; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                double tTime = fTimeAnnualised[dateIndex];
                unoptimisedDrift[dateIndex] = ((tau + tTime) * fInitialYieldCurve[tau + tTime] - tTime * fInitialYieldCurve[tTime]) / tau;
            }

            var likelihoodParameters = new ObjectiveFunctionParameters
            {
                Tau = fSelectedTenorToCalibrateTo,
                HistoricalDates = fTermStructure.Dates,
                HistoricalZeroRates = fSelectedZeroRates,
                InitialYieldCurve = fInitialYieldCurve,
                DeltaTimeAnnualised = deltaTimeAnnualised,
                TimeAnnualised = fTimeAnnualised,
                InitialAlphaEstimate = initialAlphaEstimate,
                InitialSigmaEstimate = initialSigmaEstimate,
                UnoptimisedDrift = unoptimisedDrift,
                FixedAlpha = MLE_Parameters.Alpha_Fixed
            };

            double alphaEstimate;
            double sigmaEstimate;
            GetMaximumLikelihoodEstimate(errorList, likelihoodParameters, out alphaEstimate, out sigmaEstimate);

            if (MLE_Parameters.Alpha_Fixed.Exists)
            {
                hwModel.Alpha = MLE_Parameters.Alpha_Fixed.Value;
            }
            else if (MLE_Parameters.Optimization_Parameters.Alpha_Lower_Bound.HasValue && alphaEstimate < MLE_Parameters.Optimization_Parameters.Alpha_Lower_Bound + DownhillSimplexScale)
            {
                hwModel.Alpha = MLE_Parameters.Optimization_Parameters.Alpha_Lower_Bound.Value;
            }
            else
            {
                hwModel.Alpha = alphaEstimate;
            }

            hwModel.SetConstantSigma(sigmaEstimate);

            // Get the calibrated drift for display purposes.
            fFirstCalibratedDrift = GetDrift(hwModel.Alpha, hwModel.GetConstantSigma(), likelihoodParameters);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model and writes parameters to output xml file.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var model = priceModel as HullWhite1FactorInterestRateModel;
            if (model != null)
            {
                var graphData = new ResultSeriesList { fXAxis = { fLabel = "Time (Years)" }, fYAxis = { fLabel = "Rate (%)" } };

                var driftCurve = new Curve();
                var zeroRateCurve = new Curve();
                var initialYieldCurve = new Curve();

                for (int dateIndex = 1; dateIndex < fTermStructure.Length; ++dateIndex)
                {
                    double time = fTimeAnnualised[dateIndex];
                    initialYieldCurve[time] = fInitialYieldCurve[time] * Percentage.OverPercentagePoint;
                }

                graphData.Add(new ResultSeries("Initial Zero Rates", LineStyle.Solid, initialYieldCurve));

                // Plot the drift and zero-rate curve.
                for (int dateIndex = 1; dateIndex < fTermStructure.Length; ++dateIndex)
                {
                    double time = fTimeAnnualised[dateIndex];
                    driftCurve[time] = fFirstCalibratedDrift[dateIndex] * Percentage.OverPercentagePoint;
                    zeroRateCurve[time] = fSelectedZeroRates[dateIndex] * Percentage.OverPercentagePoint;
                }

                graphData.Add(new ResultSeries(string.Format("Calibrated Drift at Tenor {0}", fSelectedTenorToCalibrateTo), LineStyle.Solid, driftCurve));
                graphData.Add(new ResultSeries(string.Format("Zero Rates at Tenor {0}", fSelectedTenorToCalibrateTo), LineStyle.Solid, zeroRateCurve));

                graphData.ToXml(output);
                CalibrationHelper.WriteModelParametersToXml(priceModel, output);
            }
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(HullWhite1FactorInterestRateModel);
        }

        /// <summary>
        /// Properties validation. Provides user with a list of errors and warnings about all properties.
        /// </summary>
        private bool ValidateProperties(ErrorList errorList, string priceModelID)
        {
            bool validProperties = MLE_Parameters.Data_Retrieval_Parameters.ValidateProperties(errorList);

            if (MLE_Parameters.Optimization_Parameters.Alpha_Lower_Bound > MLE_Parameters.Optimization_Parameters.Alpha_Upper_Bound)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Calibration of {0} - alpha bounds incorrectly set. Upperbound must be greater than LowerBound.", priceModelID));
                validProperties = false;
            }

            if (MLE_Parameters.Optimization_Parameters.Sigma_Lower_Bound >
                MLE_Parameters.Optimization_Parameters.Sigma_Upper_Bound || MLE_Parameters.Optimization_Parameters.Sigma_Lower_Bound <= 0.0)
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                                  "Calibration of {0} - sigma bounds incorrectly set.  Upperbound must be greater than LowerBound and bound must be greater than 0.",
                                  priceModelID));
                validProperties = false;
            }

            if (MLE_Parameters.Tenor <= 0.0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Calibration of {0} - the tenor for calibration is not correctly set. It must be greater than 0D.", priceModelID));
                validProperties = false;
            }

            return validProperties;
        }

        /// <summary>
        /// Returns the initial yield curve.
        /// </summary>
        private DoubleCurve GetInitialYieldCurve(double archiveLengthInYears, Period[] allHistoricalTenors, double[][] allHistoricalZeroRates)
        {
            int numTenors = allHistoricalZeroRates.Length;
            var initialYieldCurve = new RateCurve();

            double rate = allHistoricalZeroRates[0][0];
            initialYieldCurve[0] = double.IsNaN(rate) ? 0.0 : rate;

            for (int tenorIndex = 0; tenorIndex < numTenors; ++tenorIndex)
            {
                double tenor = allHistoricalTenors[tenorIndex];
                rate = allHistoricalZeroRates[tenorIndex][0];
                if (!double.IsNaN(rate))
                    initialYieldCurve[tenor] = rate;
            }

            // Set last point on yield curve to extrapolate if projected time is greater than highest tenor
            double lastTenor = allHistoricalTenors[numTenors - 1];
            if (fSelectedTenorToCalibrateTo + archiveLengthInYears > lastTenor)
            {
                rate = allHistoricalZeroRates[numTenors - 1][0];
                if (!double.IsNaN(rate))
                    initialYieldCurve[lastTenor + archiveLengthInYears] = rate;
            }

            return initialYieldCurve;
        }

        /// <summary>
        /// Creates a simplex for the Downhill Simplex Optimisation algorithm given a starting point.
        /// </summary>
        private double[,] GetSimplex(double[] x)
        {
            var simplex = new double[x.Length + 1, x.Length];

            // Populate each vertex
            for (int vertexIndex = 0; vertexIndex < simplex.GetLength(0); ++vertexIndex)
            {
                // Populate each dimension of the current vertex
                for (int dimensionIndex = 0; dimensionIndex < simplex.GetLength(1); ++dimensionIndex)
                {
                    simplex[vertexIndex, dimensionIndex] = x[dimensionIndex];
                }

                // All vertices except one are offset from x by scale in one dimension
                if (vertexIndex < simplex.GetLength(1))
                    simplex[vertexIndex, vertexIndex] += DownhillSimplexScale;
            }

            return simplex;
        }

        /// <summary>
        /// Optimise over alpha and sigma.
        /// </summary>
        private bool MLEObjective(double[] normalisedAlphaSigma, ObjectiveFunctionParameters likelihoodParameter, out double mleValue)
        {
            double alpha = normalisedAlphaSigma[0] * likelihoodParameter.InitialAlphaEstimate;
            double sigma = normalisedAlphaSigma[1] * likelihoodParameter.InitialSigmaEstimate;

            // Handle the case when alpha or sigma is very small, i.e. close to 0
            if (alpha < MLE_Parameters.Optimization_Parameters.Alpha_Lower_Bound ||
                alpha > MLE_Parameters.Optimization_Parameters.Alpha_Upper_Bound || sigma < MLE_Parameters.Optimization_Parameters.Sigma_Lower_Bound)
            {
                mleValue = Double.MaxValue;
                return false;
            }

            double[] rates = likelihoodParameter.HistoricalZeroRates;

            int lengthTimeSeries = likelihoodParameter.HistoricalDates.Length;

            var var = sigma * sigma;

            double loglikelihood = 0.0;
            double[] deltaAnnualised = likelihoodParameter.DeltaTimeAnnualised;

            double tau = GetAnnualisedTau(likelihoodParameter.Tau);
            double termBtau = (1.0 - Math.Exp(-alpha * tau)) / alpha;
            double[] drift = GetDrift(alpha, sigma, likelihoodParameter);

            for (int dateIndex = 1; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                double deltaTimeAnnualised = deltaAnnualised[dateIndex];

                double timeCorrection = Math.Exp(-alpha * deltaTimeAnnualised);

                double sValue = rates[dateIndex - 1];
                double tValue = rates[dateIndex];

                // Ignore missing values.
                if (double.IsNaN(sValue) || double.IsNaN(tValue))
                    continue;

                // Variance of innovations
                double varTransformed = var * termBtau * termBtau * (1.0 - timeCorrection * timeCorrection) /
                    (2.0 * alpha * tau * tau);

                // IID innovations with mean zero and variance varTransformed
                double iidValue = tValue - sValue * timeCorrection - drift[dateIndex] * (1.0 - timeCorrection);

                loglikelihood -= Math.Log(varTransformed) + iidValue * iidValue / (varTransformed);
            }

            // Negative loglikelihood to minimise
            mleValue = -loglikelihood;
            return true;
        }

        /// <summary>
        /// Optimise over alpha and sigma.
        /// </summary>
        private bool MLEObjectiveFixedAlpha(double[] normalisedSigma, ObjectiveFunctionParameters likelihoodParameter, out double mleValue)
        {
            bool fixedAlpha = likelihoodParameter.FixedAlpha.Exists;
            double normalisedFixedAlpha = fixedAlpha ? likelihoodParameter.FixedAlpha.Value / likelihoodParameter.InitialAlphaEstimate : 1.0;

            return MLEObjective(new[] { normalisedFixedAlpha, normalisedSigma[0] }, likelihoodParameter, out mleValue);
        }

        /// <summary>
        /// Returns the drift given an alpha and a sigma.
        /// </summary>
        private double[] GetDrift(double alpha, double sigma, ObjectiveFunctionParameters parameters)
        {
            int lengthTimeSeries = parameters.HistoricalDates.Length;

            Period tau = parameters.Tau;
            double var = sigma * sigma;
            double termBAlphatau = (1.0 - Math.Exp(-alpha * tau)) / alpha;

            var unoptimisedDrift = parameters.UnoptimisedDrift;

            var drift = new double[lengthTimeSeries];
            double[] timeAnnualised = parameters.TimeAnnualised;
            for (int dateIndex = 0; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                double tTime = timeAnnualised[dateIndex];
                double expAlphaTime = Math.Exp(-alpha * tTime);
                double termBAlphat = (1.0 - expAlphaTime) / alpha;
                double termB2AlphaT = (1.0 - expAlphaTime * expAlphaTime) / (2.0 * alpha);

                drift[dateIndex] = unoptimisedDrift[dateIndex] + 0.5 * var * termBAlphatau * (termBAlphatau * termB2AlphaT + termBAlphat * termBAlphat) / tau;
            }

            return drift;
        }

        /// <summary>
        /// Return the optimal parameters of the Vasicek model using MLE.
        /// </summary>
        private void GetMaximumLikelihoodEstimate(ErrorList errorList, ObjectiveFunctionParameters likelihoodParameters,
            out double alphaEstimate, out double sigmaEstimate)
        {
            // Given that variables are normalised, starting position is (1.0, 1.0).
            int valuationCount;
            double optimalObjective;
            bool fixedAlpha = likelihoodParameters.FixedAlpha.Exists;
            bool success;

            if (fixedAlpha)
            {
                var simplex = GetSimplex(new[] { 1.0 });
                double[] optimalNormalisedEstimates;
                success = Minimizer<ObjectiveFunctionParameters>.MinimizeDownhillSimplex(MLEObjectiveFixedAlpha, likelihoodParameters, simplex,
                    MLE_Parameters.Optimization_Parameters.Max_Iterations, MLE_Parameters.Optimization_Parameters.Fractional_Tolerance, out optimalNormalisedEstimates, out optimalObjective, out valuationCount);

                alphaEstimate = likelihoodParameters.FixedAlpha.Value;
                sigmaEstimate = optimalNormalisedEstimates[0] * likelihoodParameters.InitialSigmaEstimate;
            }
            else
            {
                var simplex = GetSimplex(new[] { 1.0, 1.0 });
                double[] optimalNormalisedEstimates;
                success = Minimizer<ObjectiveFunctionParameters>.MinimizeDownhillSimplex(MLEObjective, likelihoodParameters, simplex,
                    MLE_Parameters.Optimization_Parameters.Max_Iterations, MLE_Parameters.Optimization_Parameters.Fractional_Tolerance, out optimalNormalisedEstimates, out optimalObjective, out valuationCount);

                alphaEstimate = optimalNormalisedEstimates[0] * likelihoodParameters.InitialAlphaEstimate;
                sigmaEstimate = optimalNormalisedEstimates[1] * likelihoodParameters.InitialSigmaEstimate;
            }

            if (!success)
            {
                errorList.Add(ErrorLevel.Warning, "The number of iterations in the Downhill Simplex Optimisation was exceeded - using best estimate.");
            }
        }

        /// <summary>
        /// Calculate initial estimates for the mean reversion speed, reversion volatility and the long-run mean.
        /// </summary>
        private void GetInitialEstimate(double[] historicalZeroRates, ErrorList errorList, out double alpha, out double sigma)
        {
            alpha = 0.0;
            sigma = 0.0;

            var historicalZeroRatesAtHorizon = new double[historicalZeroRates.Length];
            Array.Copy(historicalZeroRates, 1, historicalZeroRatesAtHorizon, 0, historicalZeroRates.Length - 1);
            historicalZeroRatesAtHorizon[historicalZeroRates.Length - 1] = double.NaN;

            var calculationParameters = new MeanReversionStatisticsCalculationParameters(ReturnMethod.Diff, null, null, null, null, null);

            MeanReversionStatistics meanReversionStatistics;
            if (ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(historicalZeroRates,
                                                                               historicalZeroRatesAtHorizon,
                                                                               fTermStructure.DeltaT, calculationParameters,
                                                                               errorList, out meanReversionStatistics))
            {
                alpha = meanReversionStatistics.Alpha;
                sigma = meanReversionStatistics.Sigma;
            }
        }

        /// <summary>
        /// Return Tau in years based on the specified number of business days.
        /// </summary>
        private double GetAnnualisedTau(Period period)
        {
            Term term;
            if (Term.TryParse(period.ToString(), out term))
            {
                return term.Years + term.Months / 12.0 + term.Weeks * 5.0 / fNumBusinessDaysInYear + term.Days / fNumBusinessDaysInYear;
            }

            return period;
        }

        /// <summary>
        /// Parameters used to calculate the likelihood funtion of the MLE calibration.
        /// </summary>
        private class ObjectiveFunctionParameters
        {
            /// <summary>
            /// Annualised deltas.
            /// </summary>
            public double[] DeltaTimeAnnualised { get; set; }

            /// <summary>
            /// Annualised times.
            /// </summary>
            public double[] TimeAnnualised { get; set; }

            /// <summary>
            /// User-defined tenor for zero-rate historical data.
            /// </summary>
            public Period Tau { get; set; }

            /// <summary>
            /// Array of historical zero rates of selected tenor. This is an array with a single array if a tenor for calibration is provided.
            /// </summary>
            public double[] HistoricalZeroRates { get; set; }

            /// <summary>
            /// Array of historical dates.
            /// </summary>
            public TDate[] HistoricalDates { get; set; }

            /// <summary>
            /// Initial Yield Curve defined at all times.
            /// </summary>
            public DoubleCurve InitialYieldCurve { get; set; }

            /// <summary>
            /// Scalar to normalise the scalar alpha parameter.
            /// </summary>
            public double InitialAlphaEstimate { get; set; }

            /// <summary>
            /// Scalar to normalise the sigma parameter.
            /// </summary>
            public double InitialSigmaEstimate { get; set; }

            /// <summary>
            /// Pre-computed unoptimised part of the drift.
            /// </summary>
            public double[] UnoptimisedDrift { get; set; }

            /// <summary>
            /// User-defined property for an explicit alpha.
            /// </summary>
            public DoubleNull FixedAlpha { get; set; }
        }
    }

    /// <summary>
    /// Parameters for the MLE method.
    /// </summary>
    public class HullWhiteInterestRateMLEParameters : NestedPresentableObject
    {
        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public HullWhiteInterestRateMLEParameters()
        {
            // Data cleaning (to remove outliers) used by default.
            Data_Retrieval_Parameters = new DataRetrievalParametersWithCalendar<Period>
            {
                Data_Cleaning_Methods = 
                { 
                    new DataCleaningMethodWrapper<Period> 
                    { 
                        Method_Name = Property.DisplayName(typeof(DataCleaningZScoreMethod<Period>)),
                        Parameters = new DataCleaningZScoreMethod<Period>()
                    } 
                }
            };

            Tenor = new Period("6M");
            Optimization_Parameters = new HullWhiteOptimizationParameters();
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public DataRetrievalParametersWithCalendar<Period> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// User-defined floor on zero-rate historical data.
        /// </summary>
        public Period Tenor { get; set; }

        /// <summary>
        /// User-defined property for an explicit alpha.
        /// </summary>
        public DoubleNull Alpha_Fixed { get; set; }

        /// <summary>
        /// Nested optimization parameters.
        /// </summary>
        public HullWhiteOptimizationParameters Optimization_Parameters { get; set; }
    }

    /// <summary>
    /// Optimisation parameters.
    /// </summary>
    public class HullWhiteOptimizationParameters : NestedPresentableObject
    {
        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        public HullWhiteOptimizationParameters()
        {
            Alpha_Lower_Bound = CalcUtils.MinReversionRate;
            Sigma_Lower_Bound = CalcUtils.MinInterestRate;
            Max_Iterations = 1000;
            Fractional_Tolerance = CalcUtils.SQRT_DOUBLE_PREC;
        }

        /// <summary>
        /// User-defined property for a lower bound on alpha.
        /// </summary>
        public double? Alpha_Lower_Bound { get; set; }

        /// <summary>
        /// User-defined property for an upper bound on alpha.
        /// </summary>
        public double? Alpha_Upper_Bound { get; set; }

        /// <summary>
        /// User-defined property for a lower bound on sigma across all tenors.
        /// </summary>
        public double? Sigma_Lower_Bound { get; set; }

        /// <summary>
        /// User-defined property for a upper bound on sigma across all tenors.
        /// </summary>
        public double? Sigma_Upper_Bound { get; set; }

        /// <summary>
        /// User-defined property for the maximum number of iterations in the MLE optimisation.
        /// </summary>
        public int Max_Iterations { get; set; }

        /// <summary>
        /// Fractional tolerance used in the optimiser.
        /// </summary>
        public double Fractional_Tolerance { get; set; }
    }
}
