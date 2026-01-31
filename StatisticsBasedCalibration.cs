// <summary>
// Base class for statistics-based calibrations.
// </summary>
using System;
using System.Collections.Generic;
using System.Diagnostics;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Statistics-based and benchmark-based Calibration base class.
    /// </summary>
    public abstract class StatisticsAndBetaBasedCalibration : StatisticsBasedCalibration
    {
        /// <summary>
        /// Tag name used for IPriceFactorGrouping.
        //  The definition of systemic driver may be given in a grouping file, under column "Benchmark".
        /// </summary>
        private const string DriverMappingTagName = "Benchmark";

        /// <summary>
        /// Set defaults.
        /// </summary>
        protected StatisticsAndBetaBasedCalibration()
        {
            Benchmark = string.Empty;
        }

        /// <summary>
        /// The driver ID used as benchmark. 
        /// If this is not set then the driver specified in the grouping file will be used.
        /// </summary>
        public string Benchmark { get; set; }

        /// <summary>
        /// Calculate beta given the archive index of the price factor and the benchmark price factor.
        /// </summary>
        protected double GetBeta(double[] nameValues, double[] nameValuesAtHorizon, double[] benchmarkValues, double[] benchmarkValuesAtHorizon, double nameVol, double benchmarkVol, ReturnMethod returnMethod)
        {
            // Calculate beta from statistics and correlations.
            double correlation;
            CorrelationCalculationHelper.GetCorrelation(nameValues, nameValuesAtHorizon, benchmarkValues, benchmarkValuesAtHorizon, returnMethod, out correlation);

            return benchmarkVol > 0.0 ? correlation * nameVol / benchmarkVol : 0.0;
        }

        /// <summary>
        /// Returns the driver factor ID with a list of matched point IDs.
        /// </summary>
        /// ///<remarks>
        ///  If no match found, return null.
        /// If the list of matched point IDs is null, we are matching all points IDs of price factor and benchmark price factor. 
        /// If the input field "Benchmark" is populated, the given benchmarket is used. Otherwise, this method goes into historical data 
        /// archive (e.g., ADA file), finds all price factor's corresponding points, and finds their matching pairs in the factor grouping 
        /// (e.g., Grouping file). 
        /// </remarks>
        protected List<GroupingFileMatch> GetDriverID(ICalibrationData calibrationData, IPriceFactor priceFactor, ErrorList errorList)
        {
            var groupingPoints = new List<GroupingFileMatch>();

            // If the benchmark is provided, we do not get it from the grouping file.
            if (!string.IsNullOrEmpty(Benchmark))
            {
                string driverFactorPoint;
                string driverFactorType;
                FactorID driverFactorID;
                StatisticsData.ParseScalarID(Benchmark, out driverFactorType, out driverFactorID, out driverFactorPoint);
                string matchedPoint = string.IsNullOrEmpty(driverFactorPoint) ? string.Empty : driverFactorPoint;

                groupingPoints.Add(
                    new GroupingFileMatch(
                        new FactorTypeAndIDAndPoint(priceFactor.TypeDotSubType(), priceFactor.GetID(), string.Empty),
                        new FactorTypeAndIDAndPoint(driverFactorType, driverFactorID, matchedPoint)));

                return groupingPoints;
            }

            // Reading from grouping file
            string priceFactorType = priceFactor.TypeDotSubType();
            FactorID priceFactorID = priceFactor.GetID();
            IPriceFactorGrouping priceFactorGrouping = calibrationData.Grouping;

            IHistoricalArchive archive = calibrationData.GetHistoricalMarketRates(
                new FactorTypeAndIDList { priceFactor }, calibrationData.ArchiveStartDate, calibrationData.ArchiveEndDate);

            foreach (ArchiveItemID item in archive.ItemIDs)
            {
                if (StatisticsData.GetTypeForFactor(item.fPriceFactorSubType) == priceFactorType && item.fPriceFactorID.CompareTo(priceFactorID) == 0)
                {
                    string riskFactorString = StatisticsData.BuildScalarID(priceFactorType, priceFactorID, item.fPointID);
                    string candidateDriver = priceFactorGrouping.GetTagValue(riskFactorString, DriverMappingTagName);

                    string candidateDriverFactorTypeName;
                    FactorID candidateDriverFactorID;
                    string candidateDriverPoint;
                    StatisticsData.ParseScalarID(candidateDriver, out candidateDriverFactorTypeName, out candidateDriverFactorID, out candidateDriverPoint);

                    // Both entries must be of the same type and the driver must have an ID.
                    if (string.IsNullOrEmpty(candidateDriverFactorID.ToCode()) || priceFactorType != candidateDriverFactorTypeName)
                        continue;

                    // Matched points may have different point IDs.
                    groupingPoints.Add(
                            new GroupingFileMatch(
                                new FactorTypeAndIDAndPoint(priceFactorType, priceFactorID, item.fPointID),
                                new FactorTypeAndIDAndPoint(candidateDriverFactorTypeName, candidateDriverFactorID, candidateDriverPoint)));
                }
            }

            // No match found.
            if (groupingPoints.Count == 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                    "Info about either systemic driver for {0} in grouping file or historical data for {0} is not available.",
                    priceFactor.GetKey()));
            }
            else
            {
                foreach (GroupingFileMatch match in groupingPoints)
                {
                    errorList.Add(ErrorLevel.Info,
                                  string.Format("Matched entry {0} to Benchmark {1}.", match.fResponse, match.fDriver));
                }
            }

            return groupingPoints;
        }

        /// <summary>
        /// Matching in a grouping file (price factor and a point is matched to another price factor and a point).
        /// </summary>
        protected struct GroupingFileMatch
        {
            /// <summary>
            /// Driver price factor and point.
            /// </summary>
            public readonly FactorTypeAndIDAndPoint fDriver;

            /// <summary>
            /// Response price factor and point.
            /// </summary>
            public FactorTypeAndIDAndPoint fResponse;

            /// <summary>
            /// Set the match.
            /// </summary>
            public GroupingFileMatch(FactorTypeAndIDAndPoint response, FactorTypeAndIDAndPoint driver)
            {
                fResponse = response;
                fDriver = driver;
            }
        }
    }

    /// <summary>
    /// Base class to calculate statistics from the historical archive or to read statistics from the statistics set. 
    /// </summary>
    public abstract class StatisticsBasedCalibration : PreComputedStatisticsBasedCalibration
    {
        /// <summary>
        /// Set defaults.
        /// </summary>
        protected StatisticsBasedCalibration()
        {
            Use_Pre_Computed_Statistics = YesNo.Yes;
            Data_Retrieval_Parameters   = new CalibratorDataRetrievalParameters<string>();
        }

        /// <summary>
        /// Use the statistics set if Yes; otherwise use archive market data.
        /// </summary>
        public YesNo Use_Pre_Computed_Statistics    { get; set; }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<string> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// Calculate the moment statistics (mean and vol) for a particular price factor of a model.
        /// </summary>
        /// <remarks>
        /// Return true if statistics successfully calculated.
        /// </remarks>
        protected bool GetMeanAndVolMomentStatistics(ICalibrationData calibrationData, PriceModel priceModel, ReturnMethod returnMethod, ErrorList errorList,
                                                       ref TimeSeriesStructure<string> timeSeriesStructure, out MomentStatistics[] momentStatistics,
                                                       out bool[] isValidStatistics)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(priceModel);
            bool isSuccess = GetMeanAndVolMomentStatistics(calibrationData, new FactorTypeAndID(priceFactor.TypeDotSubType(), priceFactor.GetID()),
                                                 returnMethod, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics);

            if (!isSuccess)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No moment statistics for {0}.", priceFactor.GetKey()));
                return false;
            }

            Debug.Assert(momentStatistics != null);
            Debug.Assert(momentStatistics.Length > 0);

            return true;
        }

        /// <summary>
        /// Calculate the mean-reversion statistics (alpha, sigma and theta) for a particular price factor.
        /// </summary>
        /// <remarks>
        /// Also retrieves the string ids and calculates an array of flag for valid statistics of each time series.
        /// Return true if statistics successfully calculated.
        /// </remarks>
        protected bool GetMeanReversionStatistics(ICalibrationData calibrationData, PriceModel model, MeanReversionStatisticsMultiPointCalculationParameters calculationParameters,
                                                  ErrorList errorList, ref TimeSeriesStructure<string> timeSeries, out MeanReversionStatistics[] meanReversionStatistics,
                                                  out bool[] isValidStatistics)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);
            if (!GetMeanReversionStatistics(calibrationData, new FactorTypeAndID(priceFactor.TypeDotSubType(), priceFactor.GetID()),
                                              calculationParameters, errorList, ref timeSeries, out meanReversionStatistics, out isValidStatistics))
            {
                errorList.Add(ErrorLevel.Error, string.Format("No mean-reversion statistics for {0}.", priceFactor.GetKey()));
                return false;
            }

            Debug.Assert(meanReversionStatistics != null);
            Debug.Assert(meanReversionStatistics.Length > 0);

            return true;
        }

        /// <summary>
        /// Calculate the moment statistics (mean and vol) for a particular price factor of a model.
        /// Return true if statistics successfully calculated.
        /// </summary>
        protected bool GetMeanAndVolMomentStatistics(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID, ReturnMethod returnMethod, ErrorList errorList,
                                                     ref TimeSeriesStructure<string> timeSeriesStructure, out MomentStatistics[] momentStatistics,
                                                     out bool[] isValidStatistics)
        {
            if (Use_Pre_Computed_Statistics == YesNo.Yes)
            {
                isValidStatistics = GetMomentStatisticsFromSet(calibrationData, factorTypeAndID, returnMethod, out momentStatistics);

                if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                                                "Calibration using pre-computed {0} Moment Statistics for {1} failed.", returnMethod,
                                                GetFactorTypeAndIDDescription(factorTypeAndID)));

                    isValidStatistics = GetMomentStatisticsFromArchive(calibrationData, factorTypeAndID, returnMethod, errorList, ref timeSeriesStructure, out momentStatistics);

                    if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                    {
                        errorList.Add(ErrorLevel.Error, string.Format(
                                          "Calibration using {0} Moment Statistics calculated from archive for {1} failed.", returnMethod,
                                          GetFactorTypeAndIDDescription(factorTypeAndID)));
                        return false;
                    }

                    errorList.Add(ErrorLevel.Warning, string.Format(
                                      "Calibration using {0} Moment Statistics calculated from archive for {1} instead.", returnMethod,
                                      GetFactorTypeAndIDDescription(factorTypeAndID)));
                    return true;
                }

                errorList.Add(ErrorLevel.Info, string.Format(
                                  "Calibration using pre-computed {0} Moment Statistics for {1} successful.", returnMethod,
                                  GetFactorTypeAndIDDescription(factorTypeAndID)));
                return true;
            }

            isValidStatistics = GetMomentStatisticsFromArchive(calibrationData, factorTypeAndID, returnMethod, errorList, ref timeSeriesStructure, out momentStatistics);

            if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                                         "Calibration using {0} Moment Statistics for {1} calculated from archive for {1} failed.", returnMethod,
                                         GetFactorTypeAndIDDescription(factorTypeAndID)));

                isValidStatistics = GetMomentStatisticsFromSet(calibrationData, factorTypeAndID, returnMethod, out momentStatistics);

                if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                                      "Calibration using pre-computed {0} Moment Statistics for {1} failed.", returnMethod,
                                      GetFactorTypeAndIDDescription(factorTypeAndID)));
                    return false;
                }

                errorList.Add(ErrorLevel.Warning, string.Format(
                                  "Calibration using pre-computed {0} Moment Statistics for {1} instead.", returnMethod,
                                  GetFactorTypeAndIDDescription(factorTypeAndID)));
                return true;
            }

            errorList.Add(ErrorLevel.Info, string.Format(
                              "Calibration using {0} Moment Statistics calculated from archive for {1} successful.", returnMethod,
                              GetFactorTypeAndIDDescription(factorTypeAndID)));
            return true;
        }

        /// <summary>
        /// Calculate the mean-reversion statistics (alpha, sigma and theta) for a particular price factor of a model.
        /// </summary>
        /// <remarks>
        /// Also retrieves the string ids and calculates an array of flag for valid statistics of each time series.
        /// Return true if statistics successfully calculated.
        /// </remarks>
        protected bool GetMeanReversionStatistics(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID, 
                                                  MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList, 
                                                  ref TimeSeriesStructure<string> timeSeries, out MeanReversionStatistics[] meanReversionStatistics,
                                                  out bool[] isValidStatistics)
        {
            if (Use_Pre_Computed_Statistics == YesNo.Yes)
            {
                isValidStatistics = GetMeanReversionStatisticsFromSet(calibrationData, factorTypeAndID, calculationParameters.ReturnMethod, out meanReversionStatistics);

                if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                                             "Calibration using pre-computed {0} Mean-Reversion Statistics for {1} failed.", calculationParameters.ReturnMethod,
                                             GetFactorTypeAndIDDescription(factorTypeAndID)));

                    isValidStatistics = GetMeanReversionStatisticsFromArchive(calibrationData, factorTypeAndID, calculationParameters, errorList, 
                                                                              ref timeSeries, out meanReversionStatistics);

                    if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                    {
                        errorList.Add(ErrorLevel.Error, string.Format(
                                          "Calibration using {0} Mean-Reversion Statistics calculated from archive for {1} failed.", calculationParameters.ReturnMethod,
                                          GetFactorTypeAndIDDescription(factorTypeAndID)));
                        return false;
                    }

                    errorList.Add(ErrorLevel.Warning, string.Format(
                                      "Calibration using {0} Mean-Reversion Statistics calculated from archive for {1} instead.", calculationParameters.ReturnMethod,
                                      GetFactorTypeAndIDDescription(factorTypeAndID)));
                    return true;
                }

                errorList.Add(ErrorLevel.Info, string.Format(
                                  "Calibration using pre-computed {0} Mean-Reversion Statistics for {1} successful.", calculationParameters.ReturnMethod,
                                  GetFactorTypeAndIDDescription(factorTypeAndID)));
                return true;
            }

            isValidStatistics = GetMeanReversionStatisticsFromArchive(calibrationData, factorTypeAndID, calculationParameters, errorList, 
                                                                      ref timeSeries, out meanReversionStatistics);

            if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                                         "Calibration using {0} Mean-Reversion Statistics for {1} calculated from archive for {1} failed.", calculationParameters.ReturnMethod,
                                         GetFactorTypeAndIDDescription(factorTypeAndID)));

                isValidStatistics = GetMeanReversionStatisticsFromSet(calibrationData, factorTypeAndID, calculationParameters.ReturnMethod, out meanReversionStatistics);

                if (CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                                      "Calibration using pre-computed {0} Mean-Reversion Statistics for {1} failed.", calculationParameters.ReturnMethod,
                                      GetFactorTypeAndIDDescription(factorTypeAndID)));
                    return false;
                }

                errorList.Add(ErrorLevel.Warning, string.Format(
                                  "Calibration using pre-computed {0} Mean-Reversion Statistics for {1} instead.", calculationParameters.ReturnMethod,
                                  GetFactorTypeAndIDDescription(factorTypeAndID)));
                return true;
            }

            errorList.Add(ErrorLevel.Info, string.Format(
                              "Calibration using {0} Mean-Reversion Statistics calculated from archive for {1} successful.", calculationParameters.ReturnMethod,
                              GetFactorTypeAndIDDescription(factorTypeAndID)));
            return true;
        }

        /// <summary>
        /// Returns true if moment statistics successfully calculated from archive data.
        /// </summary>
        protected bool[] GetMomentStatisticsFromArchive(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID, ReturnMethod returnMethod, ErrorList errorList,
                                                      ref TimeSeriesStructure<string> timeSeriesStructure, out MomentStatistics[] momentStatistics)
        {
            return MomentStatisticsCalculationHelper.ComputeMomentStatistics(calibrationData, factorTypeAndID, Data_Retrieval_Parameters, returnMethod,
                                                                             errorList, out momentStatistics, ref timeSeriesStructure);
        }

        /// <summary>
        /// Get drift or GBM drift from LOG statistics of spot price.
        /// </summary>
        protected bool TryGetDriftFromLogStatistics(ICalibrationData calibrationData, PriceModel model, ErrorList errorList, bool calcGbmDrift, out double drift)
        {
            drift = 0.0;

            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            IAssetPrice priceFactorToEvolve = (IAssetPrice)CalibrationHelper.GetPriceFactor(model);
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, new FactorTypeAndID(priceFactorToEvolve.GetType().Name, priceFactorToEvolve.GetID()), ReturnMethod.Log, errorList, ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                return false;
            }

            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics at first point for {0}.", priceFactorToEvolve.GetKey()));
                return false;
            }

            drift = calcGbmDrift ? momentStatistics[0].GetGbmDrift() : momentStatistics[0].Drift;
            return true;
        }

        /// <summary>
        /// Performs a shallow copy of the properties on this object to the target object.
        /// </summary>
        protected void ShallowCopyTo(StatisticsBasedCalibration target)
        {
            target.Data_Retrieval_Parameters   = this.Data_Retrieval_Parameters;
            target.Use_Pre_Computed_Statistics = this.Use_Pre_Computed_Statistics;
        }
        
        /// <summary>
        /// Returns true if mean-reversion statistics successfully calculated from archive data.
        /// </summary>
        /// <remarks> 
        /// Alpha is the mean-reversion speed, sigma is the mean-reversion volatility and longRunMean is the long run mean.
        /// We always use the approximate MLE method since the exact one might be time consuming and we are not interested in the correlation across tenors.
        /// </remarks>
        private bool[] GetMeanReversionStatisticsFromArchive(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID,
                                                             MeanReversionStatisticsMultiPointCalculationParameters calculationParameters, ErrorList errorList, 
                                                             ref TimeSeriesStructure<string> timeSeries,
                                                             out MeanReversionStatistics[] meanReversionStatistics)
        {
            CorrelatedMeanReversionStatisticsVectorAlpha correlatedStatistics;
            bool[] isSuccess = ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(calibrationData, factorTypeAndID, Data_Retrieval_Parameters, 
                                                                                                  calculationParameters, errorList, ref timeSeries, 
                                                                                                  out correlatedStatistics);

            if (isSuccess == null)
            {
                meanReversionStatistics = null;
                return null;
            }

            int length = correlatedStatistics.GetLength();
            meanReversionStatistics = new MeanReversionStatistics[length];
            for (int index = 0; index < length; ++index)
            {
                meanReversionStatistics[index] =
                    new MeanReversionStatistics(correlatedStatistics.PointsID[index],
                                                correlatedStatistics.MeanReversionSpeeds[index],
                                                correlatedStatistics.MeanReversionVols[index],
                                                correlatedStatistics.LongRunMeans[index]);
            }

            return isSuccess;
        }
        
        /// <summary>
        /// Return a description of the factor type and ID, e.g. InterestRateCreditSpread.EUR.SPREAD.
        /// </summary>
        private string GetFactorTypeAndIDDescription(FactorTypeAndID typeAndID)
        {
            return typeAndID.fFactorType + FactorID.CodeSeparator + typeAndID.fFactorID.ToCode();
        }
    }

    /// <summary>
    /// Base class to extract statistics from the statistics set.
    /// </summary>
    public abstract class PreComputedStatisticsBasedCalibration
    {
        /// <summary>
        /// Returns true for each time series extracted if moment statistics successfully read from the statistics set.
        /// </summary>
        protected bool[] GetMomentStatisticsFromSet(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID, ReturnMethod returnMethod, 
                                                    out MomentStatistics[] momentStatistics)
        {
            IStatistics statistics  = calibrationData.Statistics;

            List<int> archiveIndexes;
            List<string> pointsList;
            statistics.FindEntries(factorTypeAndID.fFactorType, factorTypeAndID.fFactorID, out archiveIndexes, out pointsList);

            int numTimeSeries = archiveIndexes.Count;

            if (numTimeSeries == 0)
            {
                momentStatistics = null;
                return null;
            }

            momentStatistics = new MomentStatistics[numTimeSeries];

            var isAllSuccess = new bool[numTimeSeries];

            // Get statistics for each index corresponding to a time series.
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                double mean;
                bool isMeanRead = StatisticsHelper.TryGetStatistic(statistics, archiveIndexes[timeSeriesIndex],
                                                                   BasePlugInStatisticsCalculations.DriftString, returnMethod, out mean);
                double vol;
                bool isVolRead = StatisticsHelper.TryGetStatistic(statistics, archiveIndexes[timeSeriesIndex],
                                                                  BasePlugInStatisticsCalculations.VolatilityString, returnMethod, out vol);

                isAllSuccess[timeSeriesIndex] = isVolRead && isMeanRead;

                if (isAllSuccess[timeSeriesIndex])
                    momentStatistics[timeSeriesIndex] = new MomentStatistics(pointsList[timeSeriesIndex], mean, vol);
            }

            return isAllSuccess;
        }

        /// <summary>
        /// Returns true for each time series extracted if mean-reversion statistics successfully read from the statistics set.
        /// </summary>
        protected bool[] GetMeanReversionStatisticsFromSet(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID, ReturnMethod returnMethod,
                                                           out MeanReversionStatistics[] meanReversionStatistics)
        {
            IStatistics statistics = calibrationData.Statistics;

            List<int> archiveIndexes;
            List<string> pointsList;
            statistics.FindEntries(factorTypeAndID.fFactorType, factorTypeAndID.fFactorID, out archiveIndexes, out pointsList);

            int numTimeSeries = archiveIndexes.Count;

            if (numTimeSeries == 0)
            {
                meanReversionStatistics = null;
                return null;
            }

            meanReversionStatistics = new MeanReversionStatistics[numTimeSeries];

            var isAllSuccess = new bool[numTimeSeries];

            // Get statistics for each index corresponding to a time series.
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                double alpha;
                bool isAlphaRead        = StatisticsHelper.TryGetStatistic(statistics, archiveIndexes[timeSeriesIndex],
                                                                    BasePlugInStatisticsCalculations.MeanRevSpeedString, returnMethod, out alpha);
                double sigma;
                bool isSigmaRead        = StatisticsHelper.TryGetStatistic(statistics, archiveIndexes[timeSeriesIndex],
                                                                    BasePlugInStatisticsCalculations.ReversionVolString, returnMethod, out sigma);
                double theta;
                bool isLongRunMeanRead  = StatisticsHelper.TryGetStatistic(statistics, archiveIndexes[timeSeriesIndex],
                                                                    BasePlugInStatisticsCalculations.LongRunMeanString, returnMethod, out theta);

                isAllSuccess[timeSeriesIndex] = isAlphaRead && isSigmaRead && isLongRunMeanRead;

                if (isAllSuccess[timeSeriesIndex])
                    meanReversionStatistics[timeSeriesIndex] = new MeanReversionStatistics(pointsList[timeSeriesIndex], alpha, sigma, theta);
            }

            return isAllSuccess;
        }

        /// <summary>
        /// Get the correlation matrix from the correlation statistics set.
        /// </summary>
        protected bool GetCorrelationMatrix(ICalibrationData calibrationData, PriceModel model, out SymmetricMatrix correlationMatrix)
        {
            IStatistics statistics  = calibrationData.Statistics;
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            List<int> archiveIndexes;
            List<string> pointsList;
            statistics.FindEntries(priceFactor.TypeDotSubType(), model.fID, out archiveIndexes, out pointsList);

            int numVars = archiveIndexes.Count;
            var stdev = new double[numVars];
            
            for (int index = 0; index < numVars; index++)
            {
                stdev[index] = 1.0;
            }

            correlationMatrix = statistics.CovarianceMatrix(archiveIndexes, stdev);

            return correlationMatrix != null;
        }

        /// <summary>
        /// Get the covariance matrix from the correlation statistics set.
        /// </summary>
        protected bool GetCovarianceMatrix(ICalibrationData calibrationData, PriceModel model, double[] stdevs, out SymmetricMatrix covarianceMatrix)
        {
            IStatistics statistics  = calibrationData.Statistics;
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            List<int> archiveIndexes;
            List<string> pointsList;
            statistics.FindEntries(priceFactor.TypeDotSubType(), model.fID, out archiveIndexes, out pointsList);

            covarianceMatrix = statistics.CovarianceMatrix(archiveIndexes, stdevs);

            return covarianceMatrix != null;
        }
    }

    /// <summary>
    /// Class representing a point, i.e., a factor type, a factor ID and a point ID
    /// </summary>
    public class FactorTypeAndIDAndPoint : FactorTypeAndID
    {
        /// <summary>
        /// Construct a new instance.
        /// </summary>
        public FactorTypeAndIDAndPoint(string factorType, FactorID factorID, string pointID)
            : base(factorType, factorID)
        {
            PointID = pointID;
        }

        /// <summary>
        /// Point ID.
        /// </summary>
        public string PointID { get; set; }

        /// <summary>
        /// Return the factor type and ID.
        /// </summary>
        public FactorTypeAndID GetFactorTypeAndID()
        {
            return new FactorTypeAndID(fFactorType, fFactorID);
        }

        /// <summary>
        /// Display description, e.g. InterestRateCreditSpread.EUR.SPREAD.3M.
        /// </summary>
        public override string ToString()
        {
            return fFactorType + FactorID.CodeSeparator + fFactorID.ToCode() +
                   (string.IsNullOrEmpty(PointID) ? String.Empty : FactorID.CodeSeparator + PointID);
        }
    }
}
