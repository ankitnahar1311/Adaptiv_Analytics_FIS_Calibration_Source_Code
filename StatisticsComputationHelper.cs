// <summary>
// Helper class to calculate mean reversion statistics in the presence of missing data.
// </summary>
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.Contracts;
using System.Globalization;
using System.IO;
using System.Linq;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Compute the mean-reversion statistics.
    /// </summary>
    public static class ReversionStatisticsCalculationHelper
    {
        /// <summary>
        /// Compute mean-reversion statistics for a list of price factors with data cleaning, filling and detrending with a horizon.
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="factorTypeAndID">Factor requiring statistics.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="calculationParameters">Calculation parameters.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="timeSeries">Extracted market data.</param>
        /// <param name="meanReversionStatistics">Mean-reversion statistics.</param>
        /// <returns>Return an array of bool to flag successful statistics computation.</returns>
        public static bool[] ComputeMleReversionStatistics(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID,
                                                           DataRetrievalParametersWithCalendar<string> dataRetrievalParameters,
                                                           MeanReversionStatisticsMultiPointCalculationParameters calculationParameters,
                                                           ErrorList errorList, ref TimeSeriesStructure<string> timeSeries, 
                                                           out CorrelatedMeanReversionStatisticsVectorAlpha meanReversionStatistics)
        {
            // Create Time Series using default values if not extracted yet.
            if (timeSeries == null)
            {
                if (!GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID,
                    dataRetrievalParameters, errorList, out timeSeries))
                {
                    meanReversionStatistics = null;
                    return null;
                }
            }

            bool[] isAllSuccess = ComputeMleReversionStatistics(timeSeries.Values, timeSeries.ValuesAtHorizon,
                                                                timeSeries.DeltaT, calculationParameters, errorList, out meanReversionStatistics);

            // Set the point IDs.
            for (int timeSeriesIndex = 0; timeSeriesIndex < meanReversionStatistics.GetLength(); ++timeSeriesIndex)
                meanReversionStatistics.PointsID[timeSeriesIndex] = timeSeries.TimeSeriesIDs[timeSeriesIndex];

            return isAllSuccess;
        }

        /// <summary>
        /// Compute exact mean-reversion statistics for a vector of mean-reversion speeds.
        /// </summary>
        /// <param name="timeSeries">Two-dimensional timeseries.</param>
        /// <param name="timeSeriesAtHorizon">Two-dimensional timeseries at horizon.</param>
        /// <param name="horizon">Horizon of returns.</param>
        /// <param name="holidayCalendar">Holiday calendar used for analysing historical market data.</param>
        /// <param name="calculationParameters">Calculation parameters.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="meanReversionStatistics">Mean-reversion statistics.</param>
        /// <returns>Return an array of bool to flag successful statistics computation.</returns>
        public static bool[] ComputeMleReversionStatistics(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT, 
                                                           MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList, 
                                                           out CorrelatedMeanReversionStatisticsVectorAlpha meanReversionStatistics)
        {
            Requires.Assert(timeSeries != null);
            Requires.Assert(timeSeriesAtHorizon != null);
            Requires.Assert(timeSeriesAtHorizon.Length == timeSeries.Length);
            Requires.Assert(timeSeries.All(x => x.Length == timeSeries.First().Length));
            Requires.Assert(calculationParameters != null);

            int numTimeSeries = timeSeries.Length;
            int dataLength = timeSeries.First().Length;

            // Compute exact or approximate alpha.
            double[] alphaEstimates;
            CorrelatedOUProcessesHelper.GetAlphaEstimate(timeSeries, timeSeriesAtHorizon, deltaT, calculationParameters, errorList, out alphaEstimates);

            // Compute thetas using known alphas
            var thetaEstimates = new double[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                var thetaFixed = calculationParameters.ThetaFixed != null ? calculationParameters.ThetaFixed[index] : null;
                MeanReversionStatistics mrs;
                if (ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(timeSeries[index], timeSeriesAtHorizon[index], calculationParameters.ReturnMethod, deltaT, alphaEstimates[index], thetaFixed, false, errorList, out mrs))
                    thetaEstimates[index] = mrs.Theta;
            }

            // Bound alpha if required.
            for (int index = 0; index < numTimeSeries; ++index)
            {
                if (calculationParameters.MinAlpha != null && alphaEstimates[index] < calculationParameters.MinAlpha.Value)
                    alphaEstimates[index] = calculationParameters.MinAlpha.Value;

                if (calculationParameters.MaxAlpha != null && alphaEstimates[index] > calculationParameters.MaxAlpha.Value)
                    alphaEstimates[index] = calculationParameters.MaxAlpha.Value;
            }
            
            // Compute covariance matrix.
            var covariance = ComputeCovarianceMatrix(timeSeries, timeSeriesAtHorizon, alphaEstimates, calculationParameters.ThetaFixed, calculationParameters.ReturnMethod, deltaT, true);

            // Compute sigma.
            var sigmaEstimates = new double[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
                sigmaEstimates[index] = Math.Sqrt(covariance[index, index]);

            // Transform theta if model is log-normal.
            if (calculationParameters.ReturnMethod == ReturnMethod.Log)
            {
                for (int index = 0; index < numTimeSeries; ++index)
                    thetaEstimates[index] = MeanReversionStatistics.GetThetaOfGouProcess(alphaEstimates[index], sigmaEstimates[index], thetaEstimates[index]);
            }

            meanReversionStatistics = new CorrelatedMeanReversionStatisticsVectorAlpha(alphaEstimates, sigmaEstimates, thetaEstimates, covariance);
            var isStatisticsValid = new bool[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                isStatisticsValid[index] = !double.IsNaN(alphaEstimates[index]) && !double.IsNaN(sigmaEstimates[index]) &&
                                           !double.IsNaN(thetaEstimates[index]);
            }

            return isStatisticsValid;
        }

        /// <summary>
        /// Compute exact mean-reversion statistics for a scalar mean-reversion speed.
        /// </summary>
        /// <param name="timeSeries">Two-dimensional timeseries.</param>
        /// <param name="timeSeriesAtHorizon">Two-dimensional timeseries at horizon.</param>
        /// <param name="horizon">Horizon of returns.</param>
        /// <param name="holidayCalendar">Holiday calendar used for analysing historical market data.</param>
        /// <param name="calculationParameters">Calculation parameters.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="meanReversionStatistics">Mean-reversion statistics.</param>
        /// <returns>Return an array of bool to flag successful statistics computation.</returns>
        public static bool[] ComputeMleReversionStatistics(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT, 
                                                           MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList,
                                                           out CorrelatedMeanReversionStatisticsScalarAlpha meanReversionStatistics)
        {
            Requires.Assert(timeSeries != null);
            Requires.Assert(timeSeriesAtHorizon != null);
            Requires.Assert(timeSeriesAtHorizon.Length == timeSeries.Length);
            Requires.Assert(timeSeries.All(x => x.Length == timeSeries.First().Length));
            Requires.Assert(calculationParameters != null);

            int numTimeSeries = timeSeries.Length;

            // Compute exact or approximate alpha.
            double alphaEstimate;
            CorrelatedOUProcessesHelper.GetAlphaEstimate(timeSeries, timeSeriesAtHorizon, deltaT, calculationParameters, errorList, out alphaEstimate);
            
            // Bound alpha if required.
            if (calculationParameters.MinAlpha.HasValue && alphaEstimate < calculationParameters.MinAlpha.Value)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Warning, string.Format(CultureInfo.InvariantCulture, "Flooring mean reversion speed from {0:G} to {1:G}.", alphaEstimate, calculationParameters.MinAlpha.Value));

                alphaEstimate = calculationParameters.MinAlpha.Value;
            }

            // Bound alpha if required.
            if (calculationParameters.MaxAlpha.HasValue && alphaEstimate > calculationParameters.MaxAlpha.Value)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Warning, string.Format(CultureInfo.InvariantCulture, "Ceiling mean reversion speed from {0:G} to {1:G}.", alphaEstimate, calculationParameters.MaxAlpha.Value));

                alphaEstimate = calculationParameters.MaxAlpha.Value;
            }

            // Compute theta using alpha
            var thetaEstimates = new double[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                var thetaFixed = calculationParameters.ThetaFixed != null ? calculationParameters.ThetaFixed[index] : null;
                MeanReversionStatistics mrs;
                if (ComputeMleReversionStatistics(timeSeries[index], timeSeriesAtHorizon[index], calculationParameters.ReturnMethod, deltaT, alphaEstimate, thetaFixed, false, errorList, out mrs))
                    thetaEstimates[index] = mrs.Theta;
            }

            // Compute covariance matrix.
            var covariance = ComputeCovarianceMatrix(timeSeries, timeSeriesAtHorizon, alphaEstimate, calculationParameters.ThetaFixed, calculationParameters.ReturnMethod, deltaT, true);

            // Compute sigma.
            bool isSigmaCeiled = false;
            var sigmaEstimates = new double[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                sigmaEstimates[index] = Math.Sqrt(covariance[index, index]);

                // Bound sigma if required. No warning since it is applied on the whole curve.
                if (calculationParameters.MaxSigma != null && sigmaEstimates[index] > calculationParameters.MaxSigma.Value)
                {
                    isSigmaCeiled = true;
                    sigmaEstimates[index] = calculationParameters.MaxSigma.Value;
                }
            }

            if (errorList != null && calculationParameters.MaxSigma != null && isSigmaCeiled)
                errorList.Add(ErrorLevel.Warning, string.Format(CultureInfo.InvariantCulture, "Ceiling volatility curve to {0:G}.", calculationParameters.MaxSigma.Value));

            // Transform theta if model is log-normal.
            if (calculationParameters.ReturnMethod == ReturnMethod.Log)
            {
                for (int index = 0; index < numTimeSeries; ++index)
                    thetaEstimates[index] = MeanReversionStatistics.GetThetaOfGouProcess(alphaEstimate, sigmaEstimates[index], thetaEstimates[index]);
            }

            meanReversionStatistics = new CorrelatedMeanReversionStatisticsScalarAlpha(alphaEstimate, sigmaEstimates, thetaEstimates, covariance);
            var isStatisticsValid = new bool[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                isStatisticsValid[index] = !double.IsNaN(alphaEstimate) && !double.IsNaN(sigmaEstimates[index]) &&
                                           !double.IsNaN(thetaEstimates[index]);
            }

            return isStatisticsValid;
        }

        /// <summary>
        /// Compute mean-reversion statistics for a time series.
        /// </summary>
        /// <param name="timeSeries">Single time series.</param>
        /// <param name="timeSeriesAtHorizon">Single timeseries at horizon.</param>
        /// <param name="horizon">Horizon of returns.</param>
        /// <param name="holidayCalendar">Holiday calendar used for analysing historical market data.</param>
        /// <param name="calculationParameters">Calculation parameters.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="meanReversionStatistics">Mean-reversion statistics.</param>
        /// <returns>Return true if statistics computation successful.</returns>
        public static bool ComputeMleReversionStatistics(double[] timeSeries, double[] timeSeriesAtHorizon, double deltaT, 
                                                         MeanReversionStatisticsCalculationParameters calculationParameters, ErrorList errorList, 
                                                         out MeanReversionStatistics meanReversionStatistics)
        {
            if (!ComputeMleReversionStatistics(timeSeries, timeSeriesAtHorizon, calculationParameters.ReturnMethod, deltaT, calculationParameters.AlphaFixed, calculationParameters.ThetaFixed, true, errorList, out meanReversionStatistics))
                return false;

            // Transform theta if model is log-normal.
            if (calculationParameters.ReturnMethod == ReturnMethod.Log)
                meanReversionStatistics.TransformToThetaOfGouProcess();

            return true;
        }

        /// <summary>
        /// Compute MLE of mean reversion speed (alpha), long run mean (theta) and volatility (sigma).
        /// </summary>
        /// <returns>Returns false when there is insufficient historical data.</returns>
        /// <remarks>Does not tranform theta when return method is Log.</remarks>
        public static bool ComputeMleReversionStatistics(double[] values, double[] valuesAtHorizon, ReturnMethod returnMethod, double deltaT, double? alphaFixed, double? thetaFixed, bool calculateSigma, ErrorList errorList, out MeanReversionStatistics meanReversionStatistics)
        {
            double covar, var;
            return ComputeMleReversionStatistics(values, valuesAtHorizon, returnMethod, deltaT, alphaFixed, thetaFixed, calculateSigma, errorList, out meanReversionStatistics, out covar, out var);
        }

        /// <summary>
        /// Compute MLE of mean reversion speed (alpha), long run mean (theta) and volatility (sigma).
        /// </summary>
        /// <remarks>Also output covariance of values and valuesAtHorizon, and variance of values (but without division by sampleSize - 1).</remarks>
        /// <returns>Returns false when there is insufficient historical data.</returns>
        /// <remarks>Does not tranform theta when return method is Log.</remarks>
        public static bool ComputeMleReversionStatistics(double[] values, double[] valuesAtHorizon, ReturnMethod returnMethod, double deltaT, double? alphaFixed, double? thetaFixed, bool calculateSigma, ErrorList errorList, out MeanReversionStatistics meanReversionStatistics, out double covar, out double var)
        {
            Requires.Assert(values != null);
            Requires.Assert(valuesAtHorizon != null);
            Requires.Assert(values.Length == valuesAtHorizon.Length);
            Requires.Assert(deltaT > 0.0);

            covar = 0.0;
            var   = 0.0;

            List<double> logValues, logValuesAtHorizon;
            if (!StatisticsCalculationHelper.GetLogValues(values, valuesAtHorizon, returnMethod, false, out logValues, out logValuesAtHorizon))
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("Insufficient historical data - number of valid historical returns is less than {0}.", StatisticsCalculationHelper.MinReturnSeriesLength));

                meanReversionStatistics = null;
                return false;
            }

            int count = logValues.Count;
            double mean, meanAtHorizon;
            if (thetaFixed.HasValue)
            {
                mean          = thetaFixed.Value;
                meanAtHorizon = thetaFixed.Value;
            }
            else
            {
                mean          = logValues.Average();
                meanAtHorizon = logValuesAtHorizon.Average();
            }

            // Calculate alpha
            if (!alphaFixed.HasValue)
            {
                for (int dateIndex = 0; dateIndex < count; ++dateIndex)
                {
                    covar += (logValuesAtHorizon[dateIndex] - meanAtHorizon) * (logValues[dateIndex] - mean);
                    var   += CalcUtils.Sqr(logValues[dateIndex] - mean);
                }
            }

            double alpha, expAlphaDelta;
            ComputeAlpha(alphaFixed, deltaT, covar, var, errorList, out alpha, out expAlphaDelta);

            // Calculate theta
            double theta;
            if (thetaFixed.HasValue)
            {
                theta = thetaFixed.Value;
            }
            else
            {
                theta = alpha == 0.0 ? meanAtHorizon : (meanAtHorizon - expAlphaDelta * mean) / (1.0 - expAlphaDelta);
            }

            // Calculate sigma, if required
            double sigma = 0.0;
            if (calculateSigma)
            {
                double sum = 0.0;
                if (alpha == 0.0 || expAlphaDelta == 0.0)
                {
                    // Large or small alpha: calculate volatility of returns
                    for (int dateIndex = 0; dateIndex < count; ++dateIndex)
                        sum += CalcUtils.Sqr(logValuesAtHorizon[dateIndex] - meanAtHorizon - logValues[dateIndex] + mean);

                    sigma = Math.Sqrt(sum / (count * deltaT));
                }
                else
                {
                    for (int dateIndex = 0; dateIndex < count; ++dateIndex)
                        sum += CalcUtils.Sqr(logValuesAtHorizon[dateIndex] - meanAtHorizon - expAlphaDelta * (logValues[dateIndex] - mean));

                    sigma = Math.Sqrt(2.0 * alpha * sum / (count * (1.0 - expAlphaDelta * expAlphaDelta)));
                }
            }

            meanReversionStatistics = new MeanReversionStatistics(alpha, sigma, theta);
            return true;
        }

        /// <summary>
        /// Compute alpha and expAlphaDelta from either alphaFixed or covar and var.
        /// </summary>
        public static void ComputeAlpha(double? alphaFixed, double deltaT, double covar, double var, ErrorList errorList, out double alpha, out double expAlphaDelta)
        {
            Requires.Assert(deltaT > 0.0);

            if (alphaFixed.HasValue)
            {
                alpha = alphaFixed.Value;
                expAlphaDelta = Math.Exp(-alpha * deltaT);
            }
            else
            {
                if (var < CalcUtils.TINY)
                {
                    // Warn and set alpha to 0, expAlphaDelta to 1
                    alpha = 0.0;
                    expAlphaDelta = 1.0;
                    if (errorList != null)
                        errorList.Add(ErrorLevel.Warning, "Failed to calculate mean reversion speed because variance of data is zero. Setting mean reversion speed to zero.");
                }
                else
                {
                    expAlphaDelta = covar / var;

                    if (expAlphaDelta < CalcUtils.TINY)
                    {
                        // Warn and set alpha to a LARGE and expAlphaDelta to zero
                        alpha = CalcUtils.LARGE;
                        expAlphaDelta = 0.0;
                        if (errorList != null)
                            errorList.Add(ErrorLevel.Warning, string.Format("Failed to calculate mean reversion speed because linear regression coefficient not positive. Setting mean reversion speed to {0}.", CalcUtils.LARGE));
                    }
                    else
                    {
                        alpha = -Math.Log(expAlphaDelta) / deltaT;
                    }
                }
            }

            // Small alpha: set alpha to 0, expAlphaDelta to 1
            if (Math.Abs(alpha) < CalcUtils.MinReversionRate)
            {
                alpha = 0.0;
                expAlphaDelta = 1.0;
            }
        }

        /// <summary>
        /// Compute the covariance matrix.
        /// </summary>
        public static SymmetricMatrix ComputeCovarianceMatrix(double[][] timeSeries, double[][] timeSeriesAtHorizon, double alpha, double?[] thetas, 
                                                      ReturnMethod returnMethod, double deltaT, bool alphaAdjustment)
        {
            Requires.Assert(timeSeries != null);

            int numTimeSeries = timeSeries.Length;
            var alphas = new double[numTimeSeries];

            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
                alphas[timeSeriesIndex] = alpha;

            return ComputeCovarianceMatrix(timeSeries, timeSeriesAtHorizon, alphas, thetas, returnMethod, deltaT, alphaAdjustment);
        }

        /// <summary>
        /// Compute the covariance matrix.
        /// </summary>
        public static SymmetricMatrix ComputeCovarianceMatrix(double[][] timeSeries, double[][] timeSeriesAtHorizon, double[] alphas, double?[] thetas, 
                                                      ReturnMethod returnMethod, double deltaT, bool alphaAdjustment)
        {
            Requires.Assert(timeSeries != null);
            Contract.Requires(timeSeries.All(x => x.Length == timeSeries.First().Length));
            Requires.Assert(timeSeriesAtHorizon != null);
            Contract.Requires(timeSeriesAtHorizon.All(x => x.Length == timeSeriesAtHorizon.First().Length));
            Requires.Assert(timeSeriesAtHorizon.Length == timeSeries.Length);
            Requires.Assert(alphas != null);
            Requires.Assert(alphas.Length == timeSeries.Length);
            Contract.Requires(thetas == null || thetas.Length == timeSeries.Length);

            int numTimeSeries = timeSeries.Length;
            int dataLength = timeSeries.First().Length;

            var covariance = SymmetricMatrix.Create(numTimeSeries);

            var expAlphaDelta      = new double[numTimeSeries];
            var logValues          = new List<double>[numTimeSeries];
            var logValuesAtHorizon = new List<double>[numTimeSeries];
            var mean               = new double[numTimeSeries];
            var meanAtHorizon      = new double[numTimeSeries];

            for (int index = 0; index < numTimeSeries; ++index)
            {
                expAlphaDelta[index] = Math.Exp(-alphas[index] * deltaT);

                bool validSeries = StatisticsCalculationHelper.GetLogValues(timeSeries[index], timeSeriesAtHorizon[index], returnMethod, true, out logValues[index], out logValuesAtHorizon[index]);

                if (thetas != null && thetas[index].HasValue)
                {
                    mean[index]          = thetas[index].Value;
                    meanAtHorizon[index] = thetas[index].Value;
                }
                else if (validSeries)
                {
                    mean[index]          = logValues[index].Where(x => !double.IsNaN(x)).Average();
                    meanAtHorizon[index] = logValuesAtHorizon[index].Where(x => !double.IsNaN(x)).Average();
                }
            }

            for (int timeSeriesIndexA = 0; timeSeriesIndexA < numTimeSeries; ++timeSeriesIndexA)
            {
                for (int timeSeriesIndexB = 0; timeSeriesIndexB < numTimeSeries; ++timeSeriesIndexB)
                {
                    double sumCovar = 0.0;
                    int validDataLength = 0;
                    for (int dateIndex = 0; dateIndex < dataLength; ++dateIndex)
                    {
                        if (!double.IsNaN(logValues[timeSeriesIndexA][dateIndex]) && 
                            !double.IsNaN(logValues[timeSeriesIndexB][dateIndex]))
                        {
                            sumCovar += (logValuesAtHorizon[timeSeriesIndexA][dateIndex] - meanAtHorizon[timeSeriesIndexA] - expAlphaDelta[timeSeriesIndexA] * (logValues[timeSeriesIndexA][dateIndex] - mean[timeSeriesIndexA])) *
                                        (logValuesAtHorizon[timeSeriesIndexB][dateIndex] - meanAtHorizon[timeSeriesIndexB] - expAlphaDelta[timeSeriesIndexB] * (logValues[timeSeriesIndexB][dateIndex] - mean[timeSeriesIndexB]));

                            validDataLength++;
                        }
                    }

                    if (validDataLength > 0)
                        covariance[timeSeriesIndexA, timeSeriesIndexB] = sumCovar / validDataLength;

                    if (alphaAdjustment)
                        covariance[timeSeriesIndexA, timeSeriesIndexB] *= (alphas[timeSeriesIndexA] + alphas[timeSeriesIndexB]) / (1.0 - expAlphaDelta[timeSeriesIndexA] * expAlphaDelta[timeSeriesIndexB]);
                }
            }

            return covariance;
        }
    }

    /// <summary>
    /// Compute moment statistics.
    /// </summary>
    public static class MomentStatisticsCalculationHelper
    {
        /// <summary>
        /// Compute moment statistics for a time series with data cleaning, filling and detrending.
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="factorTypeAndID">Price factor to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="returnMethod">Return method can be diff or log.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="momentStatistics">Moment statistics.</param>
        /// <param name="timeSeriesStructure">Extracted market data.</param>
        /// <returns>Return an array of bool to flag successful statistics computation.</returns>
        public static bool[] ComputeMomentStatistics(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID,
            DataRetrievalParametersWithCalendar<string> dataRetrievalParameters, ReturnMethod returnMethod, ErrorList errorList,
            out MomentStatistics[] momentStatistics, ref TimeSeriesStructure<string> timeSeriesStructure)
        {
            // Create Time Series using default values.
            if (timeSeriesStructure == null)
            {
                if (!GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID,
                    dataRetrievalParameters, errorList, out timeSeriesStructure))
                {
                    momentStatistics = null;
                    return null;
                }
            }

            bool[] isValidStatistics = ComputeMomentStatistics(timeSeriesStructure, 
                                                               returnMethod, errorList, out momentStatistics);

            // Exit if no moment statistics calculated.
            if (isValidStatistics == null || CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                return isValidStatistics;

            for (int timeSeriesIndex = 0; timeSeriesIndex < momentStatistics.Length; ++timeSeriesIndex)
                momentStatistics[timeSeriesIndex].PointID = timeSeriesStructure.TimeSeriesIDs[timeSeriesIndex];

            return isValidStatistics;
        }

        /// <summary>
        /// Get moment statistics given a time series.
        /// </summary>
        public static bool[] ComputeMomentStatistics(TimeSeriesStructure<string> timeSeriesStructure, ReturnMethod returnMethod, 
                                                     ErrorList errorList, out MomentStatistics[] momentStatistics)
        {
            momentStatistics = null;

            if (timeSeriesStructure == null)
                return null;

            return ComputeMomentStatistics(timeSeriesStructure.Values, timeSeriesStructure.ValuesAtHorizon, timeSeriesStructure.DeltaT, 
                                           returnMethod, errorList, out momentStatistics);
        }

        /// <summary>
        /// Compute moment statistics for a time series given a horizon.
        /// </summary>
        /// <param name="timeSeries">Single time series.</param>
        /// <param name="timeSeriesAtHorizon">Single time series at horizon.</param>
        /// <param name="deltaT">Annualised length of return.</param>
        /// <param name="returnMethod">Return method can be diff or log.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="momentStatistics">Moment statistics.</param>
        /// <returns>Return true if statistics computation successful.</returns>
        public static bool ComputeMomentStatistics(double[] timeSeries, double[] timeSeriesAtHorizon, double deltaT, 
                                                   ReturnMethod returnMethod, ErrorList errorList, out MomentStatistics momentStatistics)
        {
            Requires.Assert(timeSeries != null);
            Requires.Assert(timeSeriesAtHorizon != null);

            var returns = new List<double>();
            if (!StatisticsCalculationHelper.GetReturns(timeSeries, timeSeriesAtHorizon, returnMethod, false, out returns))
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("Insufficient historical data - number of valid historical returns is less than {0}.", StatisticsCalculationHelper.MinReturnSeriesLength));

                momentStatistics = null;
                return false;
            }

            double vol         = Math.Sqrt(CalcUtils.SampleVariance(returns) / deltaT);
            double meanReturn  = CalcUtils.Mean(returns) / deltaT;

            momentStatistics = new MomentStatistics(meanReturn, vol);

            return true;
        }

        /// <summary>
        /// Compute moment statistics for a time series structure where the horizon = frequency.
        /// </summary>
        /// <param name="timeSeries">Multiple time series.</param>
        /// <param name="timeSeriesAtHorizon">Multiple time series at horizon.</param>
        /// <param name="horizon">Horizon of returns.</param>
        /// <param name="holidayCalendar">Holiday calendar used for analysing historical market data.</param>
        /// <param name="returnMethod">Return method can be diff or log.</param>
        /// <param name="errorList">List of info, warnings and errors. Can be null to disable diagnostics.</param>
        /// <param name="momentStatistics">Moment statistics.</param>
        /// <returns>Return true if statistics computation successful.</returns>
        public static bool[] ComputeMomentStatistics(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT, 
                                                     ReturnMethod returnMethod, ErrorList errorList, out MomentStatistics[] momentStatistics)
        {
            Requires.Assert(timeSeries != null);
            Requires.Assert(timeSeriesAtHorizon != null);

            int numTimeSeries = timeSeries.Length;

            momentStatistics = new MomentStatistics[numTimeSeries];
            var isSuccess = new bool[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                isSuccess[index] = ComputeMomentStatistics(timeSeries[index], timeSeriesAtHorizon[index], deltaT, returnMethod, errorList, out momentStatistics[index]);
            }

            return isSuccess;
        }
    }

    /// <summary>
    /// Common helper methods for mean-reversion and moments statistics calculation.
    /// </summary>
    public static class StatisticsCalculationHelper
    {
        /// <summary>
        /// Minimum length for time series of historical data.
        /// </summary>
        public const int MinTimeSeriesLength = 3;

        /// <summary>
        /// Minimum length for time series of historical returns.
        /// </summary>
        public const int MinReturnSeriesLength = MinTimeSeriesLength - 1;

        /// <summary>
        /// Calculate return series from value and value-at-horizon series, including missing values.
        /// </summary>
        /// <returns>True if there are at least MinReturnSeriesLength valid returns.</returns>
        public static bool GetReturns(double[] values, double[] valuesAtHorizon, ReturnMethod returnMethod, bool includeMissing, out List<double> returns)
        {
            Requires.Assert(values.Length == valuesAtHorizon.Length, "values and valuesAtHorizon must be of same length");

            int lengthTimeSeries = values.Length;
            returns = new List<double>(lengthTimeSeries);
            int numValidReturns = 0;

            for (int dateIndex = 0; dateIndex < lengthTimeSeries; ++dateIndex)
            {
                double data, dataAtHorizon;
                if (GetLogValue(values[dateIndex], valuesAtHorizon[dateIndex], returnMethod, out data, out dataAtHorizon))
                {
                    returns.Add(dataAtHorizon - data);
                    ++numValidReturns;
                }
                else if (includeMissing)
                {
                    returns.Add(double.NaN);
                }
            }

            return numValidReturns >= MinReturnSeriesLength;
        }

        /// <summary>
        /// Constructs lists of log values (Log) or values (Diff).
        /// </summary>
        /// <returns>Returns true if there are at least MinReturnSeriesLength valid pairs of log values.</returns>
        public static bool GetLogValues(double[] values, double[] valuesAtHorizon, ReturnMethod returnMethod, bool includeMissing, out List<double> logValues, out List<double> logValuesAtHorizon)
        {
            Requires.Assert(values != null);
            Requires.Assert(valuesAtHorizon != null);
            Requires.Assert(values.Length == valuesAtHorizon.Length);

            logValues = new List<double>(values.Length);
            logValuesAtHorizon = new List<double>(values.Length);
            int numValidValues = 0;

            for (int dateIndex = 0; dateIndex < values.Length; ++dateIndex)
            {
                double logValue, logValueAtHorizon;
                if (GetLogValue(values[dateIndex], valuesAtHorizon[dateIndex], returnMethod, out logValue, out logValueAtHorizon))
                {
                    logValues.Add(logValue);
                    logValuesAtHorizon.Add(logValueAtHorizon);
                    ++numValidValues;
                }
                else if (includeMissing)
                {
                    logValues.Add(double.NaN);
                    logValuesAtHorizon.Add(double.NaN);
                }
            }

            return numValidValues >= MinReturnSeriesLength;
        }

        /// <summary>
        /// Get log values for returnMethod == Log or values for returnMethod == Diff.
        /// </summary>
        /// <returns>Returns true if logValue and logValueAtHorizon are valid.</returns>
        private static bool GetLogValue(double value, double valueAtHorizon, ReturnMethod returnMethod, out double logValue, out double logValueAtHorizon)
        {
            if (double.IsNaN(value) || double.IsNaN(valueAtHorizon) ||
                returnMethod == ReturnMethod.Log && (value <= 0.0 || valueAtHorizon <= 0.0))
            {
                logValue = double.NaN;
                logValueAtHorizon = double.NaN;
                return false;
            }

            logValue = returnMethod == ReturnMethod.Log ? Math.Log(value) : value;
            logValueAtHorizon = returnMethod == ReturnMethod.Log ? Math.Log(valueAtHorizon) : valueAtHorizon;
            return true;
        }
    }

    /// <summary>
    /// Compute the correlation of a term structure.
    /// </summary>
    public static class CorrelationCalculationHelper
    {
        /// <summary>
        /// Matrix type (Covariance or Correlation).
        /// </summary>
        private enum MatrixType
        {
            Covariance, Correlation
        }

        /// <summary>
        /// Calculate the covariance matrix of the returns of the time series structure.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one point in the term structure and every pair of points 
        /// has at least two dates where both returns are valid.
        /// </returns>
        /// <typeparam name="T">Type of point ids of market data, e.g., string for asset prices or Period for interest rates.</typeparam>
        public static bool GetCovarianceMatrix<T>(TimeSeriesStructure<T> timeSeries, ReturnMethod returnMethod, out SymmetricMatrix covarianceMatrix)
        {
            return GetCovarianceOrCorrelationMatrix(timeSeries, returnMethod, MatrixType.Covariance, out covarianceMatrix);
        }

        /// <summary>
        /// Calculate the correlation matrix of the returns of the time series structure.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one point in the term structure and every pair of points 
        /// has at least two dates where both returns are valid.
        /// </returns>
        /// <typeparam name="T">Type of point ids of market data, e.g., string for asset prices or Period for interest rates.</typeparam>
        public static bool GetCorrelationMatrix<T>(TimeSeriesStructure<T> timeSeries, ReturnMethod returnMethod, out SymmetricMatrix correlationMatrix)
        {
            return GetCovarianceOrCorrelationMatrix(timeSeries, returnMethod, MatrixType.Correlation, out correlationMatrix);
        }

        /// <summary>
        /// Calculate the covariance or correlation matrix of the returns of the two time series structures.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one point in the term structure and every pair of points 
        /// has at least two dates where both returns are valid.
        /// </returns>
        /// <typeparam name="T">Type of point ids of market data, e.g., string for asset prices or Period for interest rates.</typeparam>
        public static bool GetCrossCorrelationMatrix<T>(TimeSeriesStructure<T> timeSeries1, TimeSeriesStructure<T> timeSeries2, ReturnMethod returnMethod, out SymmetricMatrix matrix)
        {
            int numTimeSeries1 = timeSeries1.NumTimeSeries;
            int numTimeSeries2 = timeSeries1.NumTimeSeries;
            if (numTimeSeries1 == 0 || numTimeSeries2 == 0)
            {
                matrix = null;
                return false;
            }

            var returnsSeries = new double[numTimeSeries1 + numTimeSeries2][];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries1; ++timeSeriesIndex)
            {
                List<double> returns;
                StatisticsCalculationHelper.GetReturns(timeSeries1.Values[timeSeriesIndex], timeSeries1.ValuesAtHorizon[timeSeriesIndex], returnMethod, true, out returns);
                returnsSeries[timeSeriesIndex] = returns.ToArray();
            }

            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries2; ++timeSeriesIndex)
            {
                List<double> returns;
                StatisticsCalculationHelper.GetReturns(timeSeries2.Values[timeSeriesIndex], timeSeries2.ValuesAtHorizon[timeSeriesIndex], returnMethod, true, out returns);
                returnsSeries[numTimeSeries1 + timeSeriesIndex] = returns.ToArray();
            }

            return GetCovarianceOrCorrelationMatrix(returnsSeries, MatrixType.Correlation, timeSeries1.DeltaT, out matrix);
        }

        /// <summary>
        /// Calculate the covariance matrix of an array of return series.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one return series and every pair of return series 
        /// has at least two dates where both returns are valid.
        /// </returns>
        public static bool GetCovarianceMatrix(double[][] returns, out SymmetricMatrix covarianceMatrix)
        {
            return GetCovarianceOrCorrelationMatrix(returns, MatrixType.Covariance, 1.0, out covarianceMatrix);
        }

        /// <summary>
        /// Calculate the correlation matrix of an array of return series.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one return series and every pair of return series 
        /// has at least two dates where both returns are valid.
        /// </returns>
        public static bool GetCorrelationMatrix(double[][] returns, out SymmetricMatrix correlationMatrix)
        {
            return GetCovarianceOrCorrelationMatrix(returns, MatrixType.Correlation, 1.0, out correlationMatrix);
        }

        /// <summary>
        /// Calculate the covariance from a pair of value and value-at-horizon series.
        /// </summary>
        public static bool GetCovariance(double[] valuesA, double[] valuesAtHorizonA, double[] valuesB, double[] valuesAtHorizonB, ReturnMethod returnMethod, double deltaT, out double covariance)
        {
            double correlation;
            return GetCovarianceAndCorrelation(valuesA, valuesAtHorizonA, valuesB, valuesAtHorizonB, returnMethod, deltaT, out covariance, out correlation);
        }

        /// <summary>
        /// Calculate the correlation from a pair of value and value-at-horizon series.
        /// </summary>
        public static bool GetCorrelation(double[] valuesA, double[] valuesAtHorizonA, double[] valuesB, double[] valuesAtHorizonB, ReturnMethod returnMethod, out double correlation)
        {
            double covariance;
            return GetCovarianceAndCorrelation(valuesA, valuesAtHorizonA, valuesB, valuesAtHorizonB, returnMethod, 1.0, out covariance, out correlation);
        }

        /// <summary>
        /// Returns a covariance matrix given a correlation matrix.
        /// </summary>
        public static SymmetricMatrix GetCovarianceMatrix(SymmetricMatrix correlationMatrix, double[] stdev)
        {
            int rows = correlationMatrix.Rows;
            if (rows != correlationMatrix.Cols)
                return null;

            var covarianceMatrix = SymmetricMatrix.Create(rows);
            for (int rowIndex = 0; rowIndex < rows; ++rowIndex)
            {
                covarianceMatrix[rowIndex, rowIndex] = stdev[rowIndex] * stdev[rowIndex];
                for (int colIndex = 0; colIndex < rowIndex; ++colIndex)
                {
                    covarianceMatrix[rowIndex, colIndex] = correlationMatrix[rowIndex, colIndex] * stdev[rowIndex] * stdev[colIndex];
                }
            }

            return covarianceMatrix;
        }

        /// <summary>
        /// Return a matrix with only the selected indices.
        /// </summary>
        public static SymmetricMatrix GetSelectedMatrix(List<int> selectedIndices, SymmetricMatrix matrix)
        {
            int length          = matrix.Rows;
            int selectedLength  = selectedIndices.Count;

            // Validation.
            if (matrix.Cols != matrix.Rows)
            {
                return null;
            }

            foreach (int index in selectedIndices)
            {
                if (index >= length || index < 0)
                    return null;
            }

            var selectedMatrix = SymmetricMatrix.Create(selectedLength);
            for (int i = 0; i < selectedLength; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    selectedMatrix[i, j] = matrix[selectedIndices[i], selectedIndices[j]];
                }
            }

            return selectedMatrix;
        }

        /// <summary>
        /// Calculate the covariance or correlation matrix of the returns of the time series structure.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one point in the term structure and every pair of points 
        /// has at least two dates where both returns are valid.
        /// </returns>
        /// <typeparam name="T">Type of point ids of market data, e.g., string for asset prices or Period for interest rates.</typeparam>
        private static bool GetCovarianceOrCorrelationMatrix<T>(TimeSeriesStructure<T> timeSeries, ReturnMethod returnMethod, MatrixType matrixType, out SymmetricMatrix matrix)
        {
            int numTimeSeries = timeSeries.NumTimeSeries;
            if (numTimeSeries == 0)
            {
                matrix = null;
                return false;
            }

            var returnsSeries = new double[numTimeSeries][];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                List<double> returns;
                StatisticsCalculationHelper.GetReturns(timeSeries.Values[timeSeriesIndex], timeSeries.ValuesAtHorizon[timeSeriesIndex], returnMethod, true, out returns);
                returnsSeries[timeSeriesIndex] = returns.ToArray();
            }

            return GetCovarianceOrCorrelationMatrix(returnsSeries, matrixType, timeSeries.DeltaT, out matrix);
        }

        /// <summary>
        /// Calculate the covariance or correlation matrix of an array of return series.
        /// </summary>
        /// <returns>
        /// Return true if there is at least one return series and every pair of return series 
        /// has at least two dates where both returns are valid.
        /// </returns>
        private static bool GetCovarianceOrCorrelationMatrix(double[][] returns, MatrixType matrixType, double deltaT, out SymmetricMatrix matrix)
        {
            int numTimeSeries = returns.Length;
            if (numTimeSeries == 0)
            {
                matrix = null;
                return false;
            }

            matrix = SymmetricMatrix.Create(numTimeSeries);

            if (matrixType == MatrixType.Correlation && numTimeSeries == 1)
            {
                matrix[0, 0] = 1.0;
                return true;
            }

            // Get correlation matrix.
            for (int indexA = 0; indexA < numTimeSeries; ++indexA)
            {
                for (int indexB = 0; indexB <= indexA; ++indexB)
                {
                    if (matrixType == MatrixType.Correlation && indexA == indexB)
                    {
                        matrix[indexA, indexA] = 1.0;
                        break;
                    }

                    double covariance, correlation;
                    if (!GetCovarianceAndCorrelation(returns[indexA], returns[indexB], deltaT, out covariance, out correlation))
                    {
                        matrix = null;
                        return false;
                    }

                    double value = matrixType == MatrixType.Correlation ? correlation : covariance;
                    matrix[indexA, indexB] = value;
                }
            }

            return true;
        }

        /// <summary>
        /// Calculate the covariance and correlation from a pair of value and value-at-horizon series.
        /// </summary>
        private static bool GetCovarianceAndCorrelation(double[] valuesA, double[] valuesAtHorizonA, double[] valuesB, double[] valuesAtHorizonB, ReturnMethod returnMethod, double deltaT, out double covariance, out double correlation)
        {
            List<double> returnsA, returnsB;
            StatisticsCalculationHelper.GetReturns(valuesA, valuesAtHorizonA, returnMethod, true, out returnsA);
            StatisticsCalculationHelper.GetReturns(valuesB, valuesAtHorizonB, returnMethod, true, out returnsB);

            return GetCovarianceAndCorrelation(returnsA.ToArray(), returnsB.ToArray(), deltaT, out covariance, out correlation);
        }

        /// <summary>
        /// Calculate the covariance and correlation from a pair of return series.
        /// </summary>
        /// <returns>True if there are least two dates for which returnsA and returnsB are both valid.</returns>
        private static bool GetCovarianceAndCorrelation(double[] returnsA, double[] returnsB, double deltaT, out double covariance, out double correlation)
        {
            Requires.Assert(returnsA.Length == returnsB.Length, "Return series must be of the same length");

            covariance  = 0.0;
            correlation = 0.0;

            var validReturnsA = returnsA.Where(x => !double.IsNaN(x)).ToArray();
            var validReturnsB = returnsB.Where(x => !double.IsNaN(x)).ToArray();

            // Require at least 1 valid return to calculate mean.
            if (validReturnsA.Length == 0 || validReturnsB.Length == 0)
                return false;

            double meanA = validReturnsA.Average();
            double meanB = validReturnsB.Average();
            double varA = 0.0;
            double varB = 0.0;

            int numValidA = 0;
            int numValidB = 0;
            int numValidAB = 0;
            for (int dateIndex = 0; dateIndex < returnsA.Length; ++dateIndex)
            {
                double a = returnsA[dateIndex];
                double b = returnsB[dateIndex];

                if (!double.IsNaN(a))
                {
                    varA += CalcUtils.Sqr(a - meanA);
                    ++numValidA;
                }

                if (!double.IsNaN(b))
                {
                    varB += CalcUtils.Sqr(b - meanB);
                    ++numValidB;
                }

                if (!double.IsNaN(a) && !double.IsNaN(b))
                {
                    covariance += (a - meanA) * (b - meanB);
                    ++numValidAB;
                }
            }

            if (varA < CalcUtils.TINY || varB < CalcUtils.TINY || numValidAB < StatisticsCalculationHelper.MinReturnSeriesLength)
                return false;

            correlation = covariance / Math.Sqrt(varA * varB);
            covariance /= deltaT * Math.Sqrt((numValidA - 1) * (numValidB - 1));

            return true;
        }
    }

    /// <summary>
    /// Find alpha vector or scalar that maximises the log likelihood.
    /// </summary>
    public static class CorrelatedOUProcessesHelper
    {
        /// <summary>
        ///  Get optimal vector alpha that maximises the log likelihood.
        /// </summary>
        public static void GetAlphaEstimate(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT,
                                            MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList, 
                                            out double[] alphaEstimates)
        {
            if (calculationParameters.AlphaFixed.HasValue)
            {
                alphaEstimates = new double[timeSeries.Length];
                for (int index = 0; index < timeSeries.Length; ++index)
                    alphaEstimates[index] = calculationParameters.AlphaFixed.Value;

                return;
            }
            
            alphaEstimates = GetInitialAlphaVectorEstimate(timeSeries, timeSeriesAtHorizon, deltaT, calculationParameters, errorList);
            if (!calculationParameters.UseExactSolution)
                return;

            var likelihoodParameters = new ObjectiveParameters
            {
                Delta                  = deltaT,
                TermStructure          = timeSeries,
                TermStructureAtHorizon = timeSeriesAtHorizon,
                CalculationParameters  = calculationParameters
            };

            MeanReversionStatisticsOptimisationParameters optimisationParameters = calculationParameters.OptimisationParameters;
            var simplex = GetSimplex(optimisationParameters.Downhill_Simplex_Scale, alphaEstimates);
            double optimalObjective;
            int valuationCount = 0;
            Minimizer<ObjectiveParameters>.MinimizeDownhillSimplex(ObjectiveFunctionVectorAlpha, likelihoodParameters, simplex, 
                                                    optimisationParameters.Max_Iterations,
                                                    optimisationParameters.Fractional_Tolerance, out alphaEstimates, out optimalObjective, out valuationCount);
        }

        /// <summary>
        ///  Get optimal scalar alpha that maximises the log likelihood.
        /// </summary>
        public static void GetAlphaEstimate(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT,
                                            MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList, 
                                            out double alphaEstimate)
        {
            if (calculationParameters.AlphaFixed.HasValue)
            {
                alphaEstimate = calculationParameters.AlphaFixed.Value;
                return;
            }

            alphaEstimate = GetInitialAlphaScalarEstimate(timeSeries, timeSeriesAtHorizon, calculationParameters, deltaT, errorList);
            if (!calculationParameters.UseExactSolution)
                return;

            var likelihoodParameters = new ObjectiveParameters
            {
                Delta                  = deltaT,
                TermStructure          = timeSeries,
                TermStructureAtHorizon = timeSeriesAtHorizon,
                CalculationParameters  = calculationParameters
            };

            double[] alphaEstimates;
            MeanReversionStatisticsOptimisationParameters optimisationParameters = calculationParameters.OptimisationParameters;
            var simplex = GetSimplex(optimisationParameters.Downhill_Simplex_Scale, new[] { alphaEstimate });
            double optimalObjective;
            int valuationCount = 0;
            Minimizer<ObjectiveParameters>.MinimizeDownhillSimplex(ObjectiveFunctionScalarAlpha, likelihoodParameters, simplex,
                                                    optimisationParameters.Max_Iterations,
                                                    optimisationParameters.Fractional_Tolerance, out alphaEstimates, out optimalObjective, out valuationCount);

            alphaEstimate = alphaEstimates[0];
        }

        /// <summary>
        /// Compute the an approximate value for alpha.
        /// </summary>
        private static double GetInitialAlphaScalarEstimate(double[][] timeSeries, double[][] timeSeriesAtHorizon, MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, double deltaT, ErrorList errorList)
        {
            Requires.Assert(timeSeries != null);
            Requires.Assert(timeSeriesAtHorizon != null);
            Requires.Assert(timeSeries.Length != 0);
            Requires.Assert(timeSeriesAtHorizon.Length != 0);

            double numerator   = 0.0;
            double denominator = 0.0;
            for (int timeSeriesIndex = 0; timeSeriesIndex < timeSeries.Length; ++timeSeriesIndex)
            {
                var thetaFixed = calculationParameters.ThetaFixed != null ? calculationParameters.ThetaFixed[timeSeriesIndex] : null;
                MeanReversionStatistics mrs;
                double covar, var;
                ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(timeSeries[timeSeriesIndex], timeSeriesAtHorizon[timeSeriesIndex], calculationParameters.ReturnMethod, deltaT, calculationParameters.AlphaFixed, thetaFixed, false, null, out mrs, out covar, out var);
                numerator   += covar;
                denominator += var;
            }

            double alpha, expAlphaDelta;
            ReversionStatisticsCalculationHelper.ComputeAlpha(calculationParameters.AlphaFixed, deltaT, numerator, denominator, errorList, out alpha, out expAlphaDelta);
            return alpha;  
        }

        /// <summary>
        /// Compute initial estimate.
        /// </summary>
        private static double[] GetInitialAlphaVectorEstimate(double[][] timeSeries, double[][] timeSeriesAtHorizon, double deltaT, MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList)
        {
            int numTimeSeries = timeSeries.Length;
            var alphaEstimates = new double[numTimeSeries];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                var thetaFixed = calculationParameters.ThetaFixed != null ? calculationParameters.ThetaFixed[index] : null;
                MeanReversionStatistics mrs;
                if (ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(timeSeries[index], timeSeriesAtHorizon[index], calculationParameters.ReturnMethod, deltaT, calculationParameters.AlphaFixed, thetaFixed, false, errorList, out mrs))
                    alphaEstimates[index] = mrs.Alpha;
            }

            return alphaEstimates;
        }

        /// <summary>
        /// Objective function. Return the determinant of the correlation matrix of the IID Samples.
        /// Minimising the log of the determinant to maximise the log likelihood.
        /// </summary>
        private static bool ObjectiveFunctionVectorAlpha(double[] alphas, ObjectiveParameters objParameters, out double objValue)
        {
            var covarMatrix = ReversionStatisticsCalculationHelper.ComputeCovarianceMatrix(objParameters.TermStructure, objParameters.TermStructureAtHorizon, alphas, objParameters.CalculationParameters.ThetaFixed, objParameters.CalculationParameters.ReturnMethod, objParameters.Delta, false);
            double det = covarMatrix.CalculateDeterminant();
            objValue = det <= 0.0 ? double.MinValue : Math.Log(det);

            return true;
        }

        /// <summary>
        /// Objective function. Return the determinant of the correlation matrix of the IID Samples.
        /// Minimising the log of the determinant to maximise the log likelihood.
        /// </summary>
        private static bool ObjectiveFunctionScalarAlpha(double[] alphas, ObjectiveParameters objParameters, out double objValue)
        {
            var covarMatrix = ReversionStatisticsCalculationHelper.ComputeCovarianceMatrix(objParameters.TermStructure, objParameters.TermStructureAtHorizon, alphas[0], objParameters.CalculationParameters.ThetaFixed, objParameters.CalculationParameters.ReturnMethod, objParameters.Delta, false);
            double det = covarMatrix.CalculateDeterminant();
            objValue = det <= 0.0 ? double.MinValue : Math.Log(det);

            return true;
        }

        /// <summary>
        /// Creates a simplex for the Downhill Simplex Optimisation algorithm given a starting point.
        /// </summary>
        private static double[,] GetSimplex(double downhillSimplexScale, double[] x)
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
                    simplex[vertexIndex, vertexIndex] += downhillSimplexScale;
            }

            return simplex;
        }

        /// <summary>
        /// Paramters used to compute objective function.
        /// </summary>
        private class ObjectiveParameters
        {
            /// <summary>
            /// Historical values.
            /// </summary>
            public double[][] TermStructure { get; set; }

            /// <summary>
            /// Historical values at horizon.
            /// </summary>
            public double[][] TermStructureAtHorizon { get; set; }

            /// <summary>
            /// Delta = annualised time between samples.
            /// </summary>
            public double Delta { get; set; }

            /// <summary>
            /// Calculation parameters.
            /// </summary>
            public MeanReversionStatisticsMultiPointCalculationParameters CalculationParameters { get; set; }
        }
    }

    /// <summary>
    /// Mean-reversion statistics container.
    /// </summary>
    public class MeanReversionStatistics
    {
        /// <summary>
        /// A minimum value used for the transform. A very small or negative value imply convergence at infinity or divergence and is incorrect.
        /// </summary>
        private const double MinReversionSpeedForTransform = 0.1;

        /// <summary>
        /// Set statistics.
        /// </summary>
        public MeanReversionStatistics(string pointID, double alpha, double sigma, double theta)
        {
            PointID = pointID;
            Alpha   = alpha;
            Sigma   = sigma;
            Theta   = theta;
        }

        /// <summary>
        /// Set statistics.
        /// </summary>
        public MeanReversionStatistics(double alpha, double sigma, double theta)
            : this(string.Empty, alpha, sigma, theta)
        {
        }

        /// <summary>
        /// Mean-reversion speed.
        /// </summary>
        public double Alpha { get; set; }

        /// <summary>
        /// Mean-reversion vol.
        /// </summary>
        public double Sigma { get; set; }

        /// <summary>
        /// Mean-reversion level.
        /// </summary>
        public double Theta { get; set; }

        /// <summary>
        /// Point ID.
        /// </summary>
        public string PointID { get; set; }

        /// <summary>
        /// Calculate long run mean of geometric OU process from long mean of OU process.
        /// </summary>
        public static double GetThetaOfGouProcess(double alpha, double sigma, double theta)
        {
            return Math.Exp(theta + 0.25 * sigma * sigma / Math.Max(MinReversionSpeedForTransform, alpha));
        }

        /// <summary>
        /// Parse the point IDs to double.
        /// </summary>
        public double ParsePointIDToDouble()
        {
            double pointD;
            try
            {
                pointD = CalcUtils.ParsePoint(PointID);
            }
            catch (ArgumentException)
            {
                throw new CalibrationException("Points in historical archive could not be parsed to a double.");
            }

            return pointD;
        }

        /// <summary>
        /// Transform long run mean of geometric OU process from long mean of OU process.
        /// </summary>
        public void TransformToThetaOfGouProcess()
        {
            Theta = GetThetaOfGouProcess(Alpha, Sigma, Theta);
        }
    }

    /// <summary>
    /// Moment statistics container.
    /// </summary>
    public class MomentStatistics
    {
        /// <summary>
        /// Set statistics.
        /// </summary>
        public MomentStatistics(string pointID, double drift, double vol)
        {
            PointID = pointID;
            Drift   = drift;
            Vol     = vol;
        }

        /// <summary>
        /// Set statistics.
        /// </summary>
        public MomentStatistics(double drift, double vol)
            : this(string.Empty, drift, vol)
        {
        }

        /// <summary>
        /// Drift without a volatility correction.
        /// </summary>
        public double Drift     { get; set; }

        /// <summary>
        /// Volatility statistics.
        /// </summary>
        public double Vol       { get; set; }

        /// <summary>
        /// Point ID.
        /// </summary>
        public string PointID   { get; set; }

        /// <summary>
        /// Return the volatility-corrected drift.
        /// </summary>
        public double GetGbmDrift()
        {
            return Drift + 0.5 * Vol * Vol;
        }

        /// <summary>
        /// Parse the point IDs to double.
        /// </summary>
        public double ParsePointIDToDouble()
        {
            double pointD;
            try
            {
                pointD = CalcUtils.ParsePoint(PointID);
            }
            catch (ArgumentException)
            {
                throw new CalibrationException("Points in historical archive could not be parsed to a double.");
            }

            return pointD;
        }
    }

    /// <summary>
    /// Mean-reversion statistics container.
    /// </summary>
    public class CorrelatedMeanReversionStatisticsVectorAlpha
    {
        /// <summary>
        /// Set statistics.
        /// </summary>
        public CorrelatedMeanReversionStatisticsVectorAlpha(string[] points, double[] alpha, double[] sigma, double[] theta, SymmetricMatrix covarianceMatrix)
        {
            PointsID            = points;
            MeanReversionSpeeds = alpha;
            MeanReversionVols   = sigma;
            LongRunMeans        = theta;
            CovarianceMatrix    = covarianceMatrix;
        }

        /// <summary>
        /// Set statistics.
        /// </summary>
        public CorrelatedMeanReversionStatisticsVectorAlpha(double[] alpha, double[] sigma, double[] theta, SymmetricMatrix covarianceMatrix)
        {
            PointsID            = new string[sigma.Length];
            MeanReversionSpeeds = alpha;
            MeanReversionVols   = sigma;
            LongRunMeans        = theta;
            CovarianceMatrix    = covarianceMatrix;
        }

        /// <summary>
        /// Mean-reversion speed.
        /// </summary>
        public double[] MeanReversionSpeeds { get; set; }

        /// <summary>
        /// Mean-reversion vol.
        /// </summary>
        public double[] MeanReversionVols   { get; set; }

        /// <summary>
        /// Mean-reversion level.
        /// </summary>
        public double[] LongRunMeans        { get; set; }

        /// <summary>
        /// Point ID.
        /// </summary>
        public string[] PointsID            { get; set; }

        /// <summary>
        /// Covariance matrix.
        /// </summary>
        public SymmetricMatrix CovarianceMatrix { get; set; }

        /// <summary>
        /// Return the numner of time series.
        /// </summary>
        public int GetLength()
        {
            return MeanReversionVols.Length;
        }
    }

    /// <summary>
    /// Mean-reversion statistics container.
    /// </summary>
    public class CorrelatedMeanReversionStatisticsScalarAlpha
    {
        /// <summary>
        /// Set statistics.
        /// </summary>
        public CorrelatedMeanReversionStatisticsScalarAlpha(double alpha, double[] sigma, double[] theta, SymmetricMatrix covarianceMatrix)
        {
            PointsID            = new string[sigma.Length];
            MeanReversionSpeed  = alpha;
            MeanReversionVols   = sigma;
            LongRunMeans        = theta;
            CovarianceMatrix    = covarianceMatrix;
        }

        /// <summary>
        /// Mean-reversion speed.
        /// </summary>
        public double MeanReversionSpeed    { get; set; }

        /// <summary>
        /// Mean-reversion vol.
        /// </summary>
        public double[] MeanReversionVols   { get; set; }

        /// <summary>
        /// Mean-reversion level.
        /// </summary>
        public double[] LongRunMeans        { get; set; }

        /// <summary>
        /// Point ID.
        /// </summary>
        public string[] PointsID            { get; set; }

        /// <summary>
        /// Covariance matrix.
        /// </summary>
        public SymmetricMatrix CovarianceMatrix { get; set; }

        /// <summary>
        /// Return the numner of time series.
        /// </summary>
        public int GetLength()
        {
            return MeanReversionVols.Length;
        }

        /// <summary>
        /// Return the correlation matrix.
        /// </summary>
        public SymmetricMatrix GetCorrelationMatrix()
        {
            int dim = CovarianceMatrix.Rows;
            var correlationMatrix = SymmetricMatrix.Create(dim);

            // Compute correlation matrix if required.
            for (int indexA = 0; indexA < dim; ++indexA)
            {
                correlationMatrix[indexA, indexA] = 1.0;
                for (int indexB = indexA + 1; indexB < dim; ++indexB)
                {
                    correlationMatrix[indexA, indexB] = CovarianceMatrix[indexA, indexB] / (MeanReversionVols[indexA] * MeanReversionVols[indexB]);
                }
            }

            return correlationMatrix;
        }
    }

    /// <summary>
    /// Optimisation parameters.
    /// </summary>
    public class MeanReversionStatisticsOptimisationParameters : NestedPresentableObject
    {
        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        public MeanReversionStatisticsOptimisationParameters()
        {
            Max_Iterations          = 1000;
            Fractional_Tolerance    = CalcUtils.SQRT_DOUBLE_PREC;
            Downhill_Simplex_Scale  = 0.005;
        }

        /// <summary>
        /// User-defined property for the maximum number of iterations in the MLE optimisation.
        /// </summary>
        public int Max_Iterations { get; set; }

        /// <summary>
        /// Fractional tolerance used in the optimiser.
        /// </summary>
        public double Fractional_Tolerance { get; set; }

        /// <summary>
        /// Scale used to generate simplex in downhill algorithm.
        /// </summary>
        public double Downhill_Simplex_Scale { get; set; }
    }

    /// <summary>
    /// Single-point Reversion Statistics calculation parameters.
    /// </summary>
    public class MeanReversionStatisticsCalculationParameters : MeanReversionStatisticsCalculationParametersBase
    {
        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        public MeanReversionStatisticsCalculationParameters(ReturnMethod returnMethod, double? minAlpha, double? maxAlpha, double? maxSigma, double? alphaFixed, double? thetaFixed)
            : base(returnMethod, minAlpha, maxAlpha, maxSigma, alphaFixed)
        {
            ThetaFixed = thetaFixed;
        }

        /// <summary>
        /// Flag to use a fixed theta or not.
        /// </summary>
        public double? ThetaFixed { get; set; }
    }

    /// <summary>
    /// Multi-point Reversion Statistics calculation parameters.
    /// </summary>
    public class MeanReversionStatisticsMultiPointCalculationParameters : MeanReversionStatisticsCalculationParametersBase
    {
        private bool fUseExactSolution;

        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        public MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod returnMethod, double? minAlpha,
                                                                      double? maxAlpha, double? maxSigma,
                                                                      double? alphaFixed, double?[] thetaFixed, bool useExactSolution)
            : base(returnMethod, minAlpha, maxAlpha, maxSigma, alphaFixed)
        {
            ThetaFixed       = thetaFixed;
            UseExactSolution = useExactSolution;
        }

        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        public MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod returnMethod, double? minAlpha,
                                                                      double? maxAlpha, double? maxSigma,
                                                                      double? alphaFixed, bool useExactSolution)
            : this(returnMethod, minAlpha, maxAlpha, maxSigma, alphaFixed, null, useExactSolution)
        {
        }

        /// <summary>
        /// Flag to use a fixed theta or not.
        /// </summary>
        public double?[] ThetaFixed 
        { 
            get; set; 
        }

        /// <summary>
        /// Whether to use the approximation or not.
        /// </summary>
        public bool UseExactSolution
        {
            get
            {
                return fUseExactSolution;
            }
            set
            {
                OptimisationParameters  = value ? new MeanReversionStatisticsOptimisationParameters() : null;
                fUseExactSolution       = value;
            }
        }

        /// <summary>
        /// Optimisation parameters when using exact solution.
        /// </summary>
        public MeanReversionStatisticsOptimisationParameters OptimisationParameters { get; set; }
    }

    /// <summary>
    /// Single-point Reversion Statistics calculation parameters.
    /// </summary>
    public abstract class MeanReversionStatisticsCalculationParametersBase : NestedPresentableObject
    {
        /// <summary>
        /// Initialize with default parameters.
        /// </summary>
        protected MeanReversionStatisticsCalculationParametersBase(ReturnMethod returnMethod, double? minAlpha, double? maxAlpha, double? maxSigma, double? alphaFixed)
        {
            ReturnMethod        = returnMethod;
            MinAlpha            = minAlpha;
            MaxAlpha            = maxAlpha;
            MaxSigma            = maxSigma;
            AlphaFixed          = alphaFixed;
        }

        /// <summary>
        /// Return method can be diff or log.
        /// </summary>
        public ReturnMethod ReturnMethod { get; set; }

        /// <summary>
        /// Lower bound on alpha.
        /// </summary>
        public double? MinAlpha { get; set; }

        /// <summary>
        /// Upper bound on alpha.
        /// </summary>
        public double? MaxAlpha { get; set; }

        /// <summary>
        /// Upper bound on sigma.
        /// </summary>
        public double? MaxSigma { get; set; }

        /// <summary>
        /// Flag to use a fixed alpha or not.
        /// </summary>
        public double? AlphaFixed { get; set; }
    }
}