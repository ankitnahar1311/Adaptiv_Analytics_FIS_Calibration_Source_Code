using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Abstract base class for PCA calibration of a multi factor GBM model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get the pre-computed statistics. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file.
    /// </remarks>
    /// <typeparam name="TModelParameters">The type of model parameters that we are calibrating (either a <see cref="PriceModel"/> or a model parameters <see cref="PriceFactor"/>).</typeparam>
    public abstract class PCACalibration<TModelParameters>
        where TModelParameters : BaseFactor, IMultiFactorModelParameters
    {
        /// <summary>
        /// Group by all tag.
        /// </summary>
        /// <remarks>
        /// If a user is grouping by "ALL" then all factors of the target type are included in the group (e.g. all equities).
        /// </remarks>
        private const string GroupByAll = "ALL";

        /// <summary>
        /// The default tag name used for IPriceFactorGrouping.
        //  The definition of systemic driver may be given in a grouping file, under column "Benchmark".
        /// </summary>
        private const string DriverMappingDefaultTagName = "PCA GROUP";

        private int fNumPCAFactors;
        private Percentage fExplainedVariance;

        /// <summary>
        /// Default constructor.
        /// </summary>
        protected PCACalibration()
        {
            Group_By                  = GroupByAll;
            Calibrate_All_In_Group    = YesNo.No;
            Number_Of_PCA_Factors     = 3;
            Explained_Variance        = 1.0;
            Significance_Threshold    = CalcUtils.DefaultMinProcessWeight;
            Idiosyncratic_Factor      = YesNo.No;
            Data_Retrieval_Parameters = new CalibratorDataRetrievalParameters<string>();
        }

        /// <summary>
        /// If <see cref="YesNo.Yes"/> then use the statistics set.
        /// Otherwise use archive market data.
        /// </summary>
        public YesNo Use_Pre_Computed_Statistics
        {
            get;
            set;
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<string> Data_Retrieval_Parameters
        {
            get;
            set;
        }

        /// <summary>
        /// How to group. If ALL, it matches everything, otherwise it looks in the Grouping file for a matching tag. If left blank uses the "PCA GROUP" tag.
        /// </summary>
        public string Group_By
        {
            get;
            set;
        }

        /// <summary>
        /// If <see cref="YesNo.Yes"/> then calibrate the whole group.
        /// </summary>
        public YesNo Calibrate_All_In_Group
        {
            get;
            set;
        }

        /// <summary>
        /// Maximum number of factors for the PCA decomposition.
        /// </summary>
        public int Number_Of_PCA_Factors
        {
            get
            {
                return fNumPCAFactors;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Number of PCA factors must be strictly positive.");

                fNumPCAFactors = value;
            }
        }

        /// <summary>
        /// Required variance explanation.
        /// </summary>
        public Percentage Explained_Variance
        {
            get
            {
                return fExplainedVariance;
            }
            set
            {
                if (value < 0.0 || 1.0 < value)
                    throw new ArgumentException("The variance explained by the systematic drivers must be between 0% and 100%.");

                fExplainedVariance = value;
            }
        }

        /// <summary>
        /// Significance threshold.
        /// </summary>
        public double Significance_Threshold
        {
            get;
            set;
        }

        /// <summary>
        /// If <see cref="YesNo.Yes"/> then we include an idiosyncratic factor.
        /// </summary>
        public YesNo Idiosyncratic_Factor
        {
            get;
            set;
        }

        /// <summary>
        /// Get a <see cref="FactorTypeAndIDAndPoint"/> from an <see cref="IStatistics"/> object and statistics set index.
        /// </summary>
        protected static FactorTypeAndIDAndPoint GetFactorTypeAndIDAndPoint(IStatistics statistics, int idx)
        {
            string priceFactorType;
            FactorID id;
            string point;
            statistics.GetFactorDetails(idx, out priceFactorType, out id, out point);
            var currentFactorTypeIDPoint = new FactorTypeAndIDAndPoint(priceFactorType, id, point);
            return currentFactorTypeIDPoint;
        }

        /// <summary>
        /// Performs a shallow copy of the properties on this object to the target object.
        /// </summary>
        protected void ShallowCopyTo(PCACalibration<TModelParameters> target)
        {
            target.Calibrate_All_In_Group      = this.Calibrate_All_In_Group;
            target.Data_Retrieval_Parameters   = this.Data_Retrieval_Parameters;
            target.Explained_Variance          = this.Explained_Variance;
            target.Group_By                    = this.Group_By;
            target.Idiosyncratic_Factor        = this.Idiosyncratic_Factor;
            target.Number_Of_PCA_Factors       = this.Number_Of_PCA_Factors;
            target.Significance_Threshold      = this.Significance_Threshold;
            target.Use_Pre_Computed_Statistics = this.Use_Pre_Computed_Statistics;
        }

        /// <summary>
        /// Driver calibration routine.
        /// </summary>
        protected bool TryCalibrate(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            bool isSuccess = Use_Pre_Computed_Statistics == YesNo.Yes
                ? TryCalculateFromStatisticsSet(calibrationData, priceFactorToEvolve, modelParameters, errorList)
                : TryCalculateFromHistoricalArchive(calibrationData, priceFactorToEvolve, modelParameters, errorList);

            // Log success or failure, as appropriate.
            if (isSuccess)
            {
                errorList.Add(ErrorLevel.Info, string.Format("PCA Calibration of {0} successful.", modelParameters.GetKey()));
                return true;
            }

            errorList.Add(ErrorLevel.Error, string.Format("PCA Calibration of {0} failed.", modelParameters.GetKey()));
            return false;
        }

        /// <summary>
        /// Filter the list of indices to remove unwanted tenors.
        /// </summary>
        /// <remarks>
        /// This would have trivial implementation for <see cref="SpotPrice"/> models (e.g. equities), and non-trivial implementation for "curve" models (e.g. hazard rates).
        /// </remarks>
        protected virtual void RemoveInvalidTenorIndices(IStatistics statistics, ref List<int> currentIndices)
        {
            // Base implementation doesn't remove anything (typical for spot processes)
        }

        /// <summary>
        /// Filter the list of IDs to remove unwanted tenors.
        /// </summary>
        /// <remarks>
        /// This would have trivial implementation for <see cref="SpotPrice"/> models (e.g. equities), and non-trivial implementation for "curve" models (e.g. hazard rates).
        /// </remarks>
        protected virtual void RemoveInvalidTenorIDs(ref List<FactorTypeAndIDAndPoint> factorTypeAndIDAndPoint)
        {
            // Base implementation doesn't remove anything (typical for spot processes)
        }

        /// <summary>
        /// Set the drift and volatility of the price model for when <see cref="Use_Pre_Computed_Statistics"/> is set to <see cref="YesNo.Yes"/>.
        /// </summary>
        protected abstract void SetDriftAndVolStatistics(ICalibrationData calibrationData, TModelParameters modelParameters, double vol, int statisticsSetIndex);

        /// <summary>
        /// Set the drift and volatility of the price model for when <see cref="Use_Pre_Computed_Statistics"/> is set to <see cref="YesNo.No"/>.
        /// </summary>
        protected abstract void SetDriftAndVolHistorical(ICalibrationData calibrationData, TModelParameters modelParameters, double drift, double vol);

        /// <summary>
        /// For an <see cref="IModelCalibration"/>, this registers a <see cref="PriceModel"/> of the target type and with the supplied <paramref name="id"/> on the <paramref name="priceModelList"/>.
        /// For a model parameters <see cref="Bootstrapper"/> this simply creates a new instance of the model parameters <see cref="PriceFactor"/> with the supplied <paramref name="id"/>.
        /// </summary>
        protected abstract TModelParameters GetOrRegisterModelParameters(FactorID id, PriceModelList priceModelList);

        /// <summary>
        /// Method to store all calculated model parameters for an entire group.
        /// </summary>
        /// <remarks>
        /// Base implementation does not store model parameters, and will not be overridden by standard calibrations.
        /// Implementations of <see cref="IBootstrapperCalculateMultiple"/> will want to override this so that they can return all model parameters pricefactors when calibrating the entire group.
        /// </remarks>
        protected virtual void StoreModelParameters(TModelParameters currentModelParameters)
        {
            // Base implemetation does nothing.
        }

        /// <summary>
        /// Calculates the eigen-decomposition of a square matrix.
        /// </summary>
        /// <remarks>
        /// We make this virtual so that we can override it in //Dev/Adaptiv/Live/Analytics/Plugins/NativePCACalibration/.
        /// </remarks>
        protected virtual void CalculateEigenDecomposition(SymmetricMatrix covarianceMatrix, out Vector eigenValues, out Matrix eigenVectors, out bool isPositiveSemiDefinite)
        {
            Debug.Assert(covarianceMatrix.Rows == covarianceMatrix.Cols, "Covariance matrix should be square.");

            // Diagonalize the matrix and sort the eigenvalues in descending order.
            var eigenCalculateResult = covarianceMatrix.CalculateEigenVectorsAndValues();

            // Matrix is PSD iff all eigenvalues are non-negative.
            isPositiveSemiDefinite = -CalcUtils.SMALL <= eigenCalculateResult.EigenValues.MinElement();
            eigenValues = eigenCalculateResult.EigenValues;
            eigenVectors = eigenCalculateResult.EigenVectors;
        }

        /// <summary>
        /// Try get details of the PCA group to which the primary factor belongs as well the corresponding statistics set index for each group member.
        /// </summary>
        /// <remarks>
        /// May fail if our primary factor does not belong to a group, or if we are able to determine a group, but can't find statistics for the group members.
        /// </remarks>
        protected bool TryGetGroupNameAndStatisticsIndicesForGroup(ICalibrationData sd, FactorTypeAndID primaryFactorTypeAndId, ErrorList errorList, out int primaryFactorIndex, out List<int> groupIndices, out string groupName)
        {
            // Determine the group to which the primary factor belongs.
            groupName = GetGroupName(sd.Grouping, primaryFactorTypeAndId);
            if (string.IsNullOrWhiteSpace(groupName))
                errorList.Add(ErrorLevel.Warning, string.Format("No PCA group has been found for Group_By {0}", Group_By));

            primaryFactorIndex = -1;
            groupIndices = new List<int>();
            bool groupNamesExist = false;
            foreach (int index in sd.Statistics.AvailableIndexes())
            {
                // Get the details for the current entry in the statistics set.
                string currentFactorType;
                string currentPointID;
                FactorID currentFactorID;
                sd.Statistics.GetFactorDetails(index, out currentFactorType, out currentFactorID, out currentPointID);

                if (currentFactorType == primaryFactorTypeAndId.fFactorType)
                {
                    // Record which index corresponds to our primary pricefactor.
                    if (currentFactorID.Equals(primaryFactorTypeAndId.fFactorID))
                        primaryFactorIndex = index;

                    string pcaGroup = GetGroupName(sd.Grouping, new FactorTypeAndID(currentFactorType, currentFactorID));
                    if (pcaGroup != null)
                        groupNamesExist = true;

                    if (pcaGroup == groupName)
                        groupIndices.Add(index);
                }
            }

            if (!groupNamesExist && primaryFactorIndex >= 0)
            {
                // No entries in the statistics set belong to a group.
                // We have found the primary factor in the statistics set.
                // We fake a group containing all factors of the relevant type.
                groupName = "All" + primaryFactorTypeAndId.fFactorType;
            }

            // Remove any irrelevant tenors
            RemoveInvalidTenorIndices(sd.Statistics, ref groupIndices);

            if (groupIndices.Count == 0)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("No statistics found for members of PCA group {0} ", groupName));

                return false;
            }

            if (groupName == null)
            {
                // Some entries in the statistics set belong to a group (so we know that some groups exist).
                // However our primary factor does not belong to a group.
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("{0} {1} does not belong to a PCA group.", primaryFactorTypeAndId.fFactorType, primaryFactorTypeAndId.fFactorID));

                return false;
            }

            return true;
        }

        /// <summary>
        /// Try get details of the PCA group to which the primary factor belongs as well the corresponding factor ID's for the members of the group.
        /// </summary>
        /// <remarks>
        /// May fail if our primary factor does not belong to a group.
        /// </remarks>
        protected bool TryGetGroupNameAndArchiveFactorIDsForGroup(ICalibrationData calibrationData, FactorTypeAndID primaryFactorTypeAndIDAndPoint, ErrorList errorList,
                                                                  out List<FactorTypeAndIDAndPoint> validDriversFactorTypeAndIDAndPoint, out string primaryGroupName)
        {
            primaryGroupName = GetGroupName(calibrationData.Grouping, primaryFactorTypeAndIDAndPoint);
            if (string.IsNullOrEmpty(primaryGroupName))
            {
                // We were unable to find a group for the primary model parameters.
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("No PCA GROUP found for {0}.", primaryFactorTypeAndIDAndPoint));

                validDriversFactorTypeAndIDAndPoint = null;
                return false;
            }

             var factors = calibrationData.ListAvailableFactors(null);

            // Extract the entire archive.
             IHistoricalArchive archive = calibrationData.GetHistoricalMarketRates(factors, calibrationData.ArchiveStartDate, calibrationData.ArchiveEndDate);

            // Loop over each archive item, and store those that belong to the same group as our primary factor.
            bool foundTimeSeriesForPrimaryFactor = false;
            validDriversFactorTypeAndIDAndPoint = new List<FactorTypeAndIDAndPoint>();
            foreach (ArchiveItemID item in archive.ItemIDs)
            {
                var candidateFactorTypeAndIDAndPoint = new FactorTypeAndIDAndPoint(item.fPriceFactorSubType, item.fPriceFactorID, string.Empty);
                string candidateGroup = GetGroupName(calibrationData.Grouping, candidateFactorTypeAndIDAndPoint);
                if (String.Compare(candidateGroup, primaryGroupName, StringComparison.Ordinal) == 0)
                {
                    validDriversFactorTypeAndIDAndPoint.Add(new FactorTypeAndIDAndPoint(item.fPriceFactorSubType, item.fPriceFactorID, item.fPointID));
                }

                if (candidateFactorTypeAndIDAndPoint.fFactorID.Equals(primaryFactorTypeAndIDAndPoint.fFactorID) && candidateFactorTypeAndIDAndPoint.fFactorType == primaryFactorTypeAndIDAndPoint.fFactorType)
                    foundTimeSeriesForPrimaryFactor = true;
            }

            if (!foundTimeSeriesForPrimaryFactor)
            {
                // We were unable to find a time series for the primary factor in the historical archive.
                errorList.Add(ErrorLevel.Error, String.Format("Unable to find time series for the primary factor {0}", primaryFactorTypeAndIDAndPoint));
                validDriversFactorTypeAndIDAndPoint = null;
                return false;
            }

            if (validDriversFactorTypeAndIDAndPoint.Count == 0)
            {
                // Create a dummy group containing all factors of this type, e.g., AllEquityPrice.
                primaryGroupName = "All" + primaryFactorTypeAndIDAndPoint.fFactorType;

                foreach (ArchiveItemID item in archive.ItemIDs)
                {
                    // Find other price factors of the same type.
                    if (String.Compare(item.fPriceFactorSubType, primaryFactorTypeAndIDAndPoint.fFactorType, StringComparison.Ordinal) == 0)
                    {
                        validDriversFactorTypeAndIDAndPoint.Add(new FactorTypeAndIDAndPoint(item.fPriceFactorSubType, item.fPriceFactorID, item.fPointID));
                    }
                }
            }

            // Filter out tenors.
            RemoveInvalidTenorIDs(ref validDriversFactorTypeAndIDAndPoint);

            if (validDriversFactorTypeAndIDAndPoint.Count == 0)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("Unable to find valid drivers for group."));

                validDriversFactorTypeAndIDAndPoint = null;
                return false;
            }

            return true;
        }

        /// <summary>
        /// Separate the statistics indices into those with zero and non-zero vols.
        /// </summary>
        private static void SeparateZeroAndNonZeroVolIndices(ICalibrationData sd, List<int> allIndicesForGroup, List<int> nonZeroVolIndices, List<double> nonZeroVols, List<int> volZeroIndices)
        {
            foreach (int groupIndex in allIndicesForGroup)
            {
                double vol = StatisticsHelper.GetStatistic(sd.Statistics, groupIndex, BasePlugInStatisticsCalculations.VolatilityString, ReturnMethod.Log);
                if (vol > 0.0)
                {
                    nonZeroVolIndices.Add(groupIndex);
                    nonZeroVols.Add(vol);
                }
                else
                {
                    volZeroIndices.Add(groupIndex);
                }
            }
        }

        /// <summary>
        /// Returns the covariance matrix given the time series for each member in the group.
        /// </summary>
        private static SymmetricMatrix GetCovarianceMatrixFromTimeSeries(List<TimeSeriesStatisticsContainer> groupTimeSeries)
        {
            int numMembersInGroup = groupTimeSeries.Count;

            // Construct the underlying price sub-matrix with name and drivers.
            var covMatrix = SymmetricMatrix.Create(numMembersInGroup);
            for (int groupMemberIndex1 = 0; groupMemberIndex1 < numMembersInGroup; ++groupMemberIndex1)
            {
                TimeSeriesStructure<string> timeSeries1 = groupTimeSeries[groupMemberIndex1].fTimeSeries;

                for (int groupMemberIndex2 = groupMemberIndex1; groupMemberIndex2 < numMembersInGroup; ++groupMemberIndex2)
                {
                    TimeSeriesStructure<string> timeSeries2 = groupTimeSeries[groupMemberIndex2].fTimeSeries;

                    double returnsCovariance;
                    CorrelationCalculationHelper.GetCovariance(timeSeries1.Values[0], timeSeries1.ValuesAtHorizon[0], timeSeries2.Values[0], timeSeries2.ValuesAtHorizon[0], ReturnMethod.Log, timeSeries1.DeltaT, out returnsCovariance);
                    covMatrix[groupMemberIndex1, groupMemberIndex2] = returnsCovariance;
                }
            }

            return covMatrix;
        }

        /// <summary>
        /// Construct the covariance matrix from the statistics set.
        /// </summary>
        private static SymmetricMatrix GetCovarianceMatrixFromStatisticsSet(List<int> groupIndices, List<double> volatilities, IStatistics stats)
        {
            int size = groupIndices.Count;
            var covarianceMatrix = SymmetricMatrix.Create(size);
            for (int i = 0; i < size; ++i)
            {
                for (int j = i; j < size; ++j)
                {
                    double d = stats.Correlation(groupIndices[i], groupIndices[j]) * volatilities[i] * volatilities[j];
                    covarianceMatrix[i, j] = d;
                }
            }

            return covarianceMatrix;
        }

        /// <summary>
        /// Clear the weights of the supplied model and set the <see cref="IMultiFactorModelParameters.PCA_Factors"/> tag to <see cref="YesNo.Yes"/>.
        /// </summary>
        private static void ClearWeights(TModelParameters modelParameters)
        {
            modelParameters.Weights.Clear();
            modelParameters.Idiosyncratic_Weight = 0;

            // Driving factors are principal components
            modelParameters.PCA_Factors = YesNo.Yes;
        }

        /// <summary>
        /// Find the number of PCA factor required defined as min(<paramref name="maxNumPCAFactors"/>, numFromLevelOfVariance).
        /// </summary>
        /// <remarks>
        /// Where numFromLevelofVariance is the number of PCA factors needed
        /// to get <paramref name="explainedVarianceRequirement"/> of the variance explained by the systematic drivers.
        /// </remarks>
        private static int GetNumPCAFactors(bool isPositiveSemiDefinite, Vector eigenValues, double explainedVarianceRequirement, int maxNumPCAFactors)
        {
            if (!isPositiveSemiDefinite)
                return maxNumPCAFactors;

            // Keep adding factors until we have explained enough of the variance.
            double partialSum = 0.0;
            double sumEigenValues = eigenValues.Sum();
            int numRequiredFactors;
            for (numRequiredFactors = 1; numRequiredFactors <= maxNumPCAFactors; numRequiredFactors++)
            {
                partialSum += eigenValues[numRequiredFactors - 1];
                if (partialSum / sumEigenValues >= explainedVarianceRequirement)
                {
                    break;
                }
            }

            return Math.Min(maxNumPCAFactors, numRequiredFactors);
        }

        /// <summary>
        /// If priceFactor doesn't exist for the priceModel then set a default factor price.
        /// </summary>
        private static void SetPriceFactorForModel(MarketDataContainer marketData, TModelParameters primaryModelParameters, PriceModel currentParametersAsPriceModel, FactorID id)
        {
            IPriceFactor priceFactor = currentParametersAsPriceModel.GetPriceFactor();
            if (priceFactor == null)
            {
                var primaryModel = primaryModelParameters as PriceModel;
                Debug.Assert(primaryModel != null);
                priceFactor = (IPriceFactor)marketData.PriceFactors.Register(primaryModel.GetPriceFactor().GetType(), id);
                currentParametersAsPriceModel.SetPriceFactor(priceFactor);
            }

            // Need to flag calibration attempted so we don't retry a failing calibration for all of its group's members.
            currentParametersAsPriceModel.fCalibrationUsed = true;
        }

        /// <summary>
        /// Set parameters for those models with non-zero vol.
        /// </summary>
        private void SetNonZeroVolModelParameters(ICalibrationData sd, TModelParameters primaryModelParameters, List<int> nonZeroVolIndices, FactorTypeAndIDAndPoint primaryFactorTypeAndIDAndPoint, string groupName, bool isPositiveSemiDefinite, Matrix eigenVectors, Vector eigenValues, int numPCAFactors, List<double> nonZeroVols)
        {
            int numUnderlyings = nonZeroVolIndices.Count;
            for (int modelIndex = 0; modelIndex < numUnderlyings; ++modelIndex)
            {
                var currentFactorTypeIDPoint = GetFactorTypeAndIDAndPoint(sd.Statistics, nonZeroVolIndices[modelIndex]);
                TModelParameters currentModelParameters = LocateModelParameters(sd.MarketData, primaryModelParameters, primaryFactorTypeAndIDAndPoint, currentFactorTypeIDPoint, groupName);

                // do not add to PC curves when matrix isn't positive semi-definite
                if (currentModelParameters == null || !isPositiveSemiDefinite)
                    continue;

                // Set the PCA weights and drift and vol.
                ClearWeights(currentModelParameters);
                SetModelWeights(groupName, eigenVectors, eigenValues, numPCAFactors, numUnderlyings, modelIndex, nonZeroVols[modelIndex], currentModelParameters);
                SetDriftAndVolStatistics(sd, currentModelParameters, nonZeroVols[modelIndex], nonZeroVolIndices[modelIndex]);

                // Store this model in our list of calibrated models.
                StoreModelParameters(currentModelParameters);
            }
        }

        /// <summary>
        /// Set parameters for those models with zero vol.
        /// </summary>
        private void SetZeroVolModelParameters(ICalibrationData sd, TModelParameters primaryModelParameters, List<int> volZeroIndices, FactorTypeAndIDAndPoint primaryFactorTypeAndIDAndPoint, string groupName, bool isPositiveSemiDefinite)
        {
            foreach (int volZeroIndex in volZeroIndices)
            {
                var currentFactorTypeIDPoint = GetFactorTypeAndIDAndPoint(sd.Statistics, volZeroIndex);
                TModelParameters currentModelParameters = LocateModelParameters(sd.MarketData, primaryModelParameters, primaryFactorTypeAndIDAndPoint, currentFactorTypeIDPoint, groupName);

                if (currentModelParameters == null || !isPositiveSemiDefinite)
                {
                    // Don't set any values on the model parameters.
                    continue;
                }

                // Clear any weights on the model parameters and set the volatility to 0
                ClearWeights(primaryModelParameters);
                SetDriftAndVolStatistics(sd, primaryModelParameters, 0.0, volZeroIndex);

                // Store this model in our list of calibrated models.
                StoreModelParameters(currentModelParameters);
            }
        }

        /// <summary>
        /// Calculate the model parameters from the historical archive and the grouping file.
        /// </summary>
        private bool TryCalculateFromHistoricalArchive(ICalibrationData calibrationData, IPriceFactor priceFactorToEvolve, TModelParameters primaryModelParameters, ErrorList errorList)
        {
            // Type and ID of price factor.
            var primaryFactorTypeAndIDAndPoint = new FactorTypeAndIDAndPoint(priceFactorToEvolve.TypeDotSubType(), priceFactorToEvolve.GetID(), string.Empty);

            // (1) Determine which group the primary factor lies in.
            // (2) Also get a list of ID's from this group.
            List<FactorTypeAndIDAndPoint> groupIDs;
            string groupName;
            if (!TryGetGroupNameAndArchiveFactorIDsForGroup(calibrationData, primaryFactorTypeAndIDAndPoint, errorList, out groupIDs, out groupName))
            {
                // No group found for the primary factor, or no timeseries found for any members of the group.
                return false;
            }

            try
            {
                // Get time series for group members and separate into valid and invalid series.
                var groupTimeSeries = GetTimeSeriesForGroupMembers(calibrationData, errorList, groupIDs);
                List<TimeSeriesStatisticsContainer> validTimeSeries;
                List<TimeSeriesStatisticsContainer> invalidTimeSeries;
                SeparateValidAndInvalidTimeSeries(groupTimeSeries, priceFactorToEvolve, errorList, out validTimeSeries, out invalidTimeSeries);

                // If we don't have any valid time series we can't form a covariance matrix and do the eigen-decomposition.
                if (validTimeSeries.Count == 0)
                    throw new CalibrationException("A valid set of drivers is required.");

                // Calculate eigen-decomposition and get the number of PCA factors.
                SymmetricMatrix covarianceMatrix = GetCovarianceMatrixFromTimeSeries(validTimeSeries);
                Matrix eigenVectors;
                Vector eigenValues;
                bool isPositiveSemiDefinite;
                CalculateEigenDecomposition(covarianceMatrix, out eigenValues, out eigenVectors, out isPositiveSemiDefinite);
                int numPCAFactors   = GetNumPCAFactors(isPositiveSemiDefinite, eigenValues, Explained_Variance, Number_Of_PCA_Factors);

                // If calibration succeeds set models with no statistics to default values.
                // If calibration fails leave all models in original state but still flag that the calibration has been run.
                SetModelsWithInvalidStatistics(calibrationData, primaryModelParameters, invalidTimeSeries, primaryFactorTypeAndIDAndPoint, groupName, isPositiveSemiDefinite);

                // Models with non-zero volatilities.
                SetModelsWithValidStatistics(calibrationData, primaryModelParameters, validTimeSeries, primaryFactorTypeAndIDAndPoint, groupName, isPositiveSemiDefinite, eigenVectors, eigenValues, numPCAFactors);

                if (!isPositiveSemiDefinite)
                {
                    // Inform the user of a problem with the covariance matrix.
                    errorList.Add(ErrorLevel.Error, string.Format("Covariance matrix of group {0} is not positive-semidefinite.", groupName));
                    return false;
                }
            }
            catch (CalibrationException e)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Failed to calibrate group {0}.", groupName), e);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Calculate the model parameters using the statistics set and the grouping file.
        /// </summary>
        private bool TryCalculateFromStatisticsSet(ICalibrationData sd, IPriceFactor priceFactorToEvolve, TModelParameters primaryModelParameters, ErrorList errorList)
        {
            // Type and ID of price factor.
            var primaryFactorTypeAndIDAndPoint = new FactorTypeAndIDAndPoint(priceFactorToEvolve.GetType().Name, priceFactorToEvolve.GetID(), string.Empty);

            // (1) Determine which group the primary factor lies in.
            // (2) Also get a list of statistics indices from this group, as well as the statistics index for the primary factor.
            int primaryFactorIndex;
            List<int> allIndicesForGroup;
            string groupName;
            if (!TryGetGroupNameAndStatisticsIndicesForGroup(sd, primaryFactorTypeAndIDAndPoint, errorList, out primaryFactorIndex, out allIndicesForGroup, out groupName))
            {
                // No group found for the primary factor, or no statistics found for any members of the group.
                return false;
            }

            // If the factor is not the group where it is expected to be in. We
            // skip it. It wouldn't have had it's model populated anyway and we
            // would potentially be repeating an expensive calculation that wasn't
            // necessary (some clients have lots of equities in the AMD that do
            // not have statistics).
            if (primaryFactorIndex == -1)
            {
                string errorMessage = string.Format("{0} should be in group {1} but not present: is it missing from the statistics set? Calibration skipped.", primaryModelParameters.GetID(), groupName);
                errorList.Add(ErrorLevel.Error, errorMessage);
                return false;
            }

            try
            {
                // Separate the statistics indices into those with zero vols and those with non-zero vols.
                // Hopefully the number of zero vol indices is small.
                var volZeroIndices = new List<int>();
                var nonZeroVolIndices = new List<int>(allIndicesForGroup.Count);
                var nonZeroVols = new List<double>(allIndicesForGroup.Count);
                SeparateZeroAndNonZeroVolIndices(sd, allIndicesForGroup, nonZeroVolIndices, nonZeroVols, volZeroIndices);

                // If we don't have any non-zero vol indices we can't form covariance matrix and do the eigen-decomposition.
                if (nonZeroVols.Count == 0)
                    throw new CalibrationException("A valid set of drivers is required.");

                // Get the covariance matrix and its eigen-decomposition.
                var covarianceMatrix = GetCovarianceMatrixFromStatisticsSet(nonZeroVolIndices, nonZeroVols, sd.Statistics);
                Matrix eigenVectors;
                Vector eigenValues;
                bool isPositiveSemiDefinite;
                CalculateEigenDecomposition(covarianceMatrix, out eigenValues, out eigenVectors, out isPositiveSemiDefinite);
                int numPCAFactors = GetNumPCAFactors(isPositiveSemiDefinite, eigenValues, Explained_Variance, Number_Of_PCA_Factors);

                // If calibration succeeds set models for indices with zero vol to default values
                // If calibration fails leave all models in original state but still flag that the calibration has been run
                SetZeroVolModelParameters(sd, primaryModelParameters, volZeroIndices, primaryFactorTypeAndIDAndPoint, groupName, isPositiveSemiDefinite);

                // Set the PCA weights and drift/vol for the models with non-zero vol.
                SetNonZeroVolModelParameters(sd, primaryModelParameters, nonZeroVolIndices, primaryFactorTypeAndIDAndPoint, groupName, isPositiveSemiDefinite, eigenVectors, eigenValues, numPCAFactors, nonZeroVols);

                if (!isPositiveSemiDefinite)
                {
                    // Inform the user of a problem with the covariance matrix.
                    errorList.Add(ErrorLevel.Error, string.Format("Covariance matrix of group {0} is not positive-semidefinite.", groupName));
                    return false;
                }
            }
            catch (CalibrationException e)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Failed to calibrate group {0}.", groupName), e);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Return a collection of time series, hopefully one for each entry in <paramref name="factorTypeAndIDAndPoint"/>.
        /// </summary>
        /// <remarks>
        /// We may be unable to obtain a valid time series for some of the entries in <paramref name="factorTypeAndIDAndPoint"/>, in which case we log an error.
        /// </remarks>
        private List<Tuple<FactorTypeAndIDAndPoint, TimeSeriesStructure<string>>> GetTimeSeriesForGroupMembers(ICalibrationData calibrationData, ErrorList errorList, IEnumerable<FactorTypeAndIDAndPoint> factorTypeAndIDAndPoint)
        {
            var result = new List<Tuple<FactorTypeAndIDAndPoint, TimeSeriesStructure<string>>>();

            // Get drivers time series.
            foreach (FactorTypeAndIDAndPoint driverTypeAndID in factorTypeAndIDAndPoint)
            {
                TimeSeriesStructure<string> timeSeries;
                if (GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData,
                    driverTypeAndID.GetFactorTypeAndID(), Data_Retrieval_Parameters, errorList, out timeSeries))
                {
                    result.Add(new Tuple<FactorTypeAndIDAndPoint, TimeSeriesStructure<string>>(driverTypeAndID, timeSeries));
                }
            }

            return result;
        }

        /// <summary>
        /// Returns a list of valid and a list of invalid time series. 
        /// </summary>
        /// <remarks>
        /// An invalid timeseries is one which has zero volatility or no valid statistics.
        /// </remarks>
        private void SeparateValidAndInvalidTimeSeries(List<Tuple<FactorTypeAndIDAndPoint, TimeSeriesStructure<string>>> validDriversTimeSeries, IPriceFactor priceFactor, ErrorList errorList, out List<TimeSeriesStatisticsContainer> validTimeSeries, out List<TimeSeriesStatisticsContainer> invalidTimeSeries)
        {
            validTimeSeries = new List<TimeSeriesStatisticsContainer>();
            invalidTimeSeries = new List<TimeSeriesStatisticsContainer>();

            // Add driver volatilities - using first point.
            foreach (var tuple in validDriversTimeSeries)
            {
                FactorTypeAndIDAndPoint factorTypeAndIDAndPoint = tuple.Item1;
                TimeSeriesStructure<string> timeSeries = tuple.Item2;

                MomentStatistics[] momentStatistics;
                bool[] isValidStatistics = MomentStatisticsCalculationHelper.ComputeMomentStatistics(timeSeries, ReturnMethod.Log, errorList, out momentStatistics);

                int indexPoint = GetIndexOfPoint(factorTypeAndIDAndPoint, timeSeries, errorList);
                if (isValidStatistics == null || CalibrationHelper.AreAllStatisticsInvalid(isValidStatistics))
                {
                    // Exit calibration if invalid moment statistics available.
                    invalidTimeSeries.Add(new TimeSeriesStatisticsContainer(factorTypeAndIDAndPoint, timeSeries, momentStatistics[indexPoint]));
                    errorList.Add(ErrorLevel.Error, string.Format("No moment statistics for {0}.", priceFactor.GetKey()));
                }
                else if (momentStatistics[indexPoint].Vol < CalcUtils.TINY)
                {
                    // Exit calibration if volatility is zero.
                    invalidTimeSeries.Add(new TimeSeriesStatisticsContainer(factorTypeAndIDAndPoint, timeSeries, momentStatistics[indexPoint]));
                    errorList.Add(ErrorLevel.Error, string.Format("Invalid zero vol statistics for {0}.", priceFactor.GetKey()));
                }
                else
                {
                    validTimeSeries.Add(new TimeSeriesStatisticsContainer(factorTypeAndIDAndPoint, timeSeries, momentStatistics[indexPoint]));
                }
            }
        }

        /// <summary>
        /// Return the index of the point in the time series.
        /// Assume time series is not empty.
        /// </summary>
        private int GetIndexOfPoint(FactorTypeAndIDAndPoint factorTypeAndIDAndPoint, TimeSeriesStructure<string> timeSeries, ErrorList errorList)
        {
            if (string.IsNullOrEmpty(factorTypeAndIDAndPoint.PointID))
                return 0;

            for (int index = 0; index < timeSeries.NumTimeSeries; ++index)
            {
                if (String.Compare(timeSeries.TimeSeriesIDs[index], factorTypeAndIDAndPoint.PointID, StringComparison.Ordinal) == 0)
                    return index;
            }

            errorList.Add(ErrorLevel.Warning, string.Format("Could not find the time series for {0}. Using the first available time series for the price factor.", factorTypeAndIDAndPoint));
            return 0;
        }

        /// <summary>
        /// Return the price model of a driver.
        /// </summary>
        /// <remarks>
        /// If <paramref name="currentFactorTypeAndIDAndPoint"/> is that of the primary model we get the primary model back.
        /// For other members of the group this method will return null if we are not calibrating the whole group.
        /// </remarks>
        private TModelParameters LocateModelParameters(MarketDataContainer marketData, TModelParameters primaryModelParameters, FactorTypeAndIDAndPoint primaryFactorTypeAndIDAndPoint, FactorTypeAndIDAndPoint currentFactorTypeAndIDAndPoint, string groupName)
        {
            // Return primary model parameters.
            if (String.CompareOrdinal(primaryFactorTypeAndIDAndPoint.ToString(), currentFactorTypeAndIDAndPoint.ToString()) == 0)
                return primaryModelParameters;

            // If we ARE NOT calibrating the whole group then we only care about the primary model parameters, and can return null for other model parameters.
            if (Calibrate_All_In_Group == YesNo.No)
                return null;

            // If we ARE calibrating the whole group then we care about all model(parameters) in the group.
            FactorID id = currentFactorTypeAndIDAndPoint.fFactorID;
            string priceFactorType = currentFactorTypeAndIDAndPoint.fFactorType;

            if (priceFactorType != primaryFactorTypeAndIDAndPoint.fFactorType)
                throw new CalibrationException(string.Format("Group {0} has factors of the wrong type.", groupName));

            // We may be working with either a price model or a model parameters price factor.
            TModelParameters currentModelParameters = GetOrRegisterModelParameters(id, marketData.PriceModels);

            // For PriceModels we may need to set a default pricefactor.
            var currentParametersAsPriceModel = currentModelParameters as PriceModel;
            if (currentParametersAsPriceModel != null)
            {
                SetPriceFactorForModel(marketData, primaryModelParameters, currentParametersAsPriceModel, id);
            }

            return currentModelParameters;
        }

        /// <summary>
        /// Get the PCA group for id.
        /// </summary>
        private string GetGroupName(IPriceFactorGrouping grouping, FactorTypeAndID factorTypeAndID)
        {
            if (Group_By == GroupByAll)
                return GroupByAll;

            if (grouping == null)
                return null;

            string tagName = String.IsNullOrEmpty(Group_By) ? DriverMappingDefaultTagName : Group_By;
            string riskFactorString = StatisticsData.BuildScalarID(factorTypeAndID.fFactorType, factorTypeAndID.fFactorID, null);
            return grouping.GetTagValue(riskFactorString, tagName);
        }

        /// <summary>
        /// Set the model's PCA weights.
        /// </summary>
        private void SetModelWeights(string groupName, Matrix eigenVectors, Vector eigenValues, int numPCAFactors, int numUnderlyings, int targetModelIndex, double vol, TModelParameters model)
        {
            ClearWeights(model);
            Vector weights = Vector.Create(numUnderlyings);
            double sumSquared = 0.0;
            for (int otherModelIndex = 0; otherModelIndex < numUnderlyings; otherModelIndex++)
            {
                double eigenValue = eigenValues[otherModelIndex];
                if (eigenValue < 0)
                    break;

                weights[otherModelIndex] = eigenVectors[targetModelIndex, otherModelIndex] * Math.Sqrt(eigenValue) / vol;
                sumSquared += weights[otherModelIndex] * weights[otherModelIndex];
            }

            Vector principalComponentWeights = Vector.Create(numPCAFactors);
            principalComponentWeights.AssignElements(factorIndex => weights[factorIndex]);

            double epsilon;
            MultiFactorCalibrationHelper.ObtainRoundedNormalisedRegressionWeights(Significance_Threshold, sumSquared, Idiosyncratic_Factor, null, principalComponentWeights, out epsilon);
            for (int factorIndex = 0; factorIndex < numPCAFactors; factorIndex++)
            {
                model.Weights.Add(new FactorWeight
                {
                    Driver = groupName + (factorIndex + 1).ToString(CultureInfo.InvariantCulture),
                    Weight = weights[factorIndex]
                });
            }

            model.Idiosyncratic_Weight = epsilon;
            model.fDefaultUsed = false;
        }

        /// <summary>
        /// Set the model associated with the invalid time series (i.e., with zero volatilities).
        /// </summary>
        private void SetModelsWithInvalidStatistics(ICalibrationData calibrationData, TModelParameters primaryModelParameters, List<TimeSeriesStatisticsContainer> invalidTimeSeries,
                                                    FactorTypeAndIDAndPoint primaryFactorTypeAndIDAndPoint, string groupName, bool isPositiveSemiDefinite)
        {
            for (int timeSeriesIndex = 0; timeSeriesIndex < invalidTimeSeries.Count; ++timeSeriesIndex)
            {
                TModelParameters currentModelParameters = LocateModelParameters(calibrationData.MarketData, primaryModelParameters, primaryFactorTypeAndIDAndPoint, invalidTimeSeries[timeSeriesIndex].fID, groupName);
                if (currentModelParameters != null && isPositiveSemiDefinite)
                {
                    ClearWeights(currentModelParameters);
                    SetDriftAndVolHistorical(calibrationData, currentModelParameters, invalidTimeSeries[timeSeriesIndex].fMomentStatistics.Drift, CalcUtils.MinVolatility);
                    StoreModelParameters(currentModelParameters);
                }
            }
        }

        /// <summary>
        /// Set the parameters for a model associated with a valid time series.
        /// </summary>
        private void SetModelsWithValidStatistics(ICalibrationData calibrationData, TModelParameters primaryModelParameters, List<TimeSeriesStatisticsContainer> validTimeSeries, FactorTypeAndIDAndPoint primaryFactorTypeAndIDAndPoint,
                                                  string groupName, bool isPositiveSemiDefinite, Matrix eigenVectors, Vector eigenValues, int numPCAFactors)
        {
            int numValidTimeSeries = validTimeSeries.Count;
            for (int outerTimeSeriesIndex = 0; outerTimeSeriesIndex < numValidTimeSeries; ++outerTimeSeriesIndex)
            {
                // Get the model associated with type, ID and point.
                TModelParameters currentModelParameters = LocateModelParameters(calibrationData.MarketData, primaryModelParameters, primaryFactorTypeAndIDAndPoint, validTimeSeries[outerTimeSeriesIndex].fID, groupName);

                // Do not add to PC curves when matrix isn't positive semi-definite.
                if (currentModelParameters != null && isPositiveSemiDefinite)
                {
                    StoreModelParameters(currentModelParameters);

                    // Set model parameters.
                    double vol = validTimeSeries[outerTimeSeriesIndex].fMomentStatistics.Vol;
                    SetDriftAndVolHistorical(calibrationData, currentModelParameters, validTimeSeries[outerTimeSeriesIndex].fMomentStatistics.Drift, vol);
                    SetModelWeights(groupName, eigenVectors, eigenValues, numPCAFactors, numValidTimeSeries, outerTimeSeriesIndex, vol, currentModelParameters);
                }
            }
        }

        /// <summary>
        /// Container for statistics of a time series.
        /// </summary>
        private struct TimeSeriesStatisticsContainer
        {
            public readonly FactorTypeAndIDAndPoint fID;
            public readonly TimeSeriesStructure<string> fTimeSeries;
            public readonly MomentStatistics fMomentStatistics;

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public TimeSeriesStatisticsContainer(FactorTypeAndIDAndPoint id, TimeSeriesStructure<string> timeSeries, MomentStatistics momentStatistics)
            {
                fID               = id;
                fTimeSeries       = timeSeries;
                fMomentStatistics = momentStatistics;
            }
        }
    }

    /// <summary>
    /// Class for PCA calibration of a multi factor GBM asset price model.
    /// </summary>
    [DisplayName("Multi-Factor GBM Asset Price PCA Calibration")]
    public class MultiGBMAssetPricePCACalibration : PCACalibration<MultiGBMAssetPriceModel>, IModelCalibration
    {
        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var modelParameters = CalibrationHelper.Validate<MultiGBMAssetPriceModel>(calibrationData, priceModel, output);

            // Run the main calibration routine in the base class.
            IPriceFactor priceFactorToEvolve = CalibrationHelper.GetPriceFactor(priceModel);
            if (TryCalibrate(calibrationData, modelParameters, priceFactorToEvolve, errorList))
                GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// Set the vol and drift of the model.
        /// </summary>
        /// <param name="sd">Market Data and statistics set.</param>
        /// <param name="modelParameters">Price model you want to set.</param>
        /// <param name="vol">Volatility of the model.</param>
        /// <param name="index">Index of the model in the statistics set.</param>
        protected override void SetDriftAndVolStatistics(ICalibrationData sd, MultiGBMAssetPriceModel modelParameters, double vol, int index)
        {
            double driftStat = StatisticsHelper.GetGbmDrift(sd.Statistics, index, ReturnMethod.Log);
            modelParameters.Drift = driftStat;
            modelParameters.Vol = vol;
        }

        /// <inheritdoc />
        /// <remarks>Model is being used as a proxy to set by reference the priceModel parameters.</remarks>
        protected override void SetDriftAndVolHistorical(ICalibrationData calibrationData, MultiGBMAssetPriceModel modelParameters, double drift, double vol)
        {
            modelParameters.Drift = drift + 0.5 * vol * vol;
            modelParameters.Vol = vol;
        }

        /// <inheritdoc />
        protected override MultiGBMAssetPriceModel GetOrRegisterModelParameters(FactorID id, PriceModelList priceModelList)
        {
            return priceModelList.Register<MultiGBMAssetPriceModel>(id);
        }
    }

    /// <summary>
    /// Class for PCA calibration of a multi factor GBM credit model.
    /// </summary>
    [DisplayName("Multi-Factor GBM Hazard Rate PCA Calibration")]
    public class MultiGBMHazardRatePCACalibration : PCACalibration<MultiGBMHazardRateModel>, IModelCalibration
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public MultiGBMHazardRatePCACalibration()
        {
            Tenor = 5.0;
        }

        /// <summary>
        /// The target tenor to be used for calibration.
        /// </summary>
        public Period Tenor
        {
            get;
            set;
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMHazardRateModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var modelParameters = CalibrationHelper.Validate<MultiGBMHazardRateModel>(calibrationData, priceModel, output);

            // Run the main calibration routine in the base class.
            IPriceFactor priceFactorToEvolve = CalibrationHelper.GetPriceFactor(priceModel);
            if (TryCalibrate(calibrationData, modelParameters, priceFactorToEvolve, errorList))
                GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        /// <remarks>
        /// Set only the volatility since hazard rate model doesn't have Drift property.
        /// </remarks>
        protected override void SetDriftAndVolStatistics(ICalibrationData sd, MultiGBMHazardRateModel modelParameters, double vol, int index)
        {
            modelParameters.Vol = vol;
        }

        /// <inheritdoc />
        /// <remarks>
        /// Set only the volatility since hazard rate model doesn't have Drift property.
        /// </remarks>
        protected override void SetDriftAndVolHistorical(ICalibrationData calibrationData, MultiGBMHazardRateModel modelParameters, double drift, double vol)
        {
            modelParameters.Vol = vol;
        }

        /// <inheritdoc />
        /// <remarks> Return a list of points that match the user-define Parameter, Tenor.</remarks>
        protected override void RemoveInvalidTenorIDs(ref List<FactorTypeAndIDAndPoint> factorTypeAndIDAndPoint)
        {
            var filteredPoints = new List<FactorTypeAndIDAndPoint>();
            foreach (FactorTypeAndIDAndPoint factorTypeIdPoint in factorTypeAndIDAndPoint)
            {
                double tenor = CalcUtils.ParsePoint(factorTypeIdPoint.PointID);

                if ((Period)tenor == Tenor)
                    filteredPoints.Add(factorTypeIdPoint);
            }

            if (filteredPoints.Count > 0)
                factorTypeAndIDAndPoint = filteredPoints;
        }

        /// <summary>
        /// Filter the indices to keep only the indices for the specified tenor.
        /// </summary>
        protected override void RemoveInvalidTenorIndices(IStatistics statistics, ref List<int> indices)
        {
            var filteredPoints = new List<int>();
            foreach (int index in indices)
            {
                string priceFactorType;
                FactorID id;
                double tenor;
                StatisticsHelper.GetFactorDetails(statistics, index, out priceFactorType, out id, out tenor);
                if ((Period)tenor == Tenor)
                    filteredPoints.Add(index);
            }

            indices = filteredPoints;
        }

        /// <inheritdoc />
        protected override MultiGBMHazardRateModel GetOrRegisterModelParameters(FactorID id, PriceModelList priceModelList)
        {
            return priceModelList.Register<MultiGBMHazardRateModel>(id);
        }
    }

    /// <summary>
    /// PCA Calibration for the Multi GBM term-structure model.
    /// </summary>
    /// <remarks>
    /// We do not set the vol here, as the model obtains the vol from a separate <see cref="GBMTSImpliedParameters"/> pricefactor.
    /// </remarks>
    [DisplayName("Multi-Factor GBM TS Model PCA Calibration")]
    public class MultiGBMAssetPriceTSPCACalibration : PCACalibration<MultiGBMAssetPriceTSModel>, IModelCalibration
    {
        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public MultiGBMAssetPriceTSPCACalibration()
        {
            Set_Risk_Neutral_Drift = YesNo.No;
            Risk_Neutral_Drift     = YesNo.No;
            Vol_Calibration        = VolCalibration.Implied;
        }

        /// <summary>
        /// Indicates whether we wish to override the Risk_Neutral_Drift flag on the model.
        /// </summary>
        public YesNo Set_Risk_Neutral_Drift
        {
            get;
            set;
        }

        /// <summary>
        /// The value for the Risk_Neutral_Drift on the model (only used if Set_Risk_Neutral_Drift is "Yes").
        /// </summary>
        public YesNo Risk_Neutral_Drift
        {
            get;
            set;
        }

        /// <summary>
        /// If set to Historical, volatility will be set to a constant volatility, calculated from historical vols.
        /// </summary>
        public VolCalibration Vol_Calibration { get; set; }

        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceTSModelImplied);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<MultiGBMAssetPriceTSModel>(calibrationData, priceModel, output);

            // Run the main calibration routine in the base class.
            IPriceFactor priceFactorToEvolve = CalibrationHelper.GetPriceFactor(priceModel);
            if (TryCalibrate(calibrationData, model, priceFactorToEvolve, errorList))
                GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var multiGBMTSModel = (MultiGBMAssetPriceTSModel)priceModel;
            GBMAssetPriceTSCalibrationHelper.CalibrationOutput(multiGBMTSModel.Drift, multiGBMTSModel.Vol, output);
        }

        /// <inheritdoc />
        protected override void SetDriftAndVolStatistics(ICalibrationData calibrationData, MultiGBMAssetPriceTSModel modelParameters, double vol, int statisticsSetIndex)
        {
            // set the risk neutral flag of the model if requested
            if (Set_Risk_Neutral_Drift == YesNo.Yes)
                modelParameters.Risk_Neutral_Drift = Risk_Neutral_Drift;

            // Calibrate the vol.
            var priceFactorToEvolve = (IAssetPrice)CalibrationHelper.GetPriceFactor(modelParameters);

            if (Vol_Calibration == VolCalibration.Implied)
            {
                GBMAssetPriceTSCalibrationHelper.VolCalibrationFromPriceFactor((IAssetPriceWithVolatilityPriceFactor)priceFactorToEvolve, calibrationData.PriceFactors,
                                                                                modelParameters.Risk_Neutral_Drift, modelParameters, null);
            }
            else
            {
                // set the volatility to volatility from statistics set which is already calculated and passed to the method.
                modelParameters.Vol.AssignConstant(vol);
            }

            // Calibrate the drift.
            var drift = calibrationData.Statistics.GetStatistic(statisticsSetIndex, BasePlugInStatisticsCalculations.DriftString, ReturnMethod.Log);
            GBMAssetPriceTSCalibrationHelper.CreateDriftWithVolatilityCorrection(modelParameters, drift, modelParameters.GetIntegratedVol());
        }

        /// <inheritdoc />
        protected override void SetDriftAndVolHistorical(ICalibrationData calibrationData, MultiGBMAssetPriceTSModel modelParameters, double drift, double vol)
        {
            // Set the risk neutral flag of the model if requested.
            if (Set_Risk_Neutral_Drift == YesNo.Yes)
                modelParameters.Risk_Neutral_Drift = Risk_Neutral_Drift;

            // Set volatility.
            var priceFactorToEvolve = (IAssetPrice)CalibrationHelper.GetPriceFactor(modelParameters);

            if (Vol_Calibration == VolCalibration.Implied)
            {
                GBMAssetPriceTSCalibrationHelper.VolCalibrationFromPriceFactor((IAssetPriceWithVolatilityPriceFactor)priceFactorToEvolve, calibrationData.PriceFactors,
                                                                                modelParameters.Risk_Neutral_Drift, modelParameters, null);
            }
            else
            {
                // set the volatility to historical volatility which is already calculated and passed to the method.
                modelParameters.Vol.AssignConstant(vol);
            }

            // Set drift.
            GBMAssetPriceTSCalibrationHelper.CreateDriftWithVolatilityCorrection(modelParameters, drift, modelParameters.GetIntegratedVol());
        }

        /// <inheritdoc />
        protected override MultiGBMAssetPriceTSModel GetOrRegisterModelParameters(FactorID id, PriceModelList priceModelList)
        {
            return priceModelList.Register<MultiGBMAssetPriceTSModel>(id);
        }
    }

    /// <summary>
    /// PCA Calibration for the Multi GBM term-structure "Implied" model.
    /// </summary>
    /// <remarks>
    /// We do not set the vol here, as the model obtains the vol from a separate <see cref="GBMTSImpliedParameters"/> pricefactor.
    /// </remarks>
    [DisplayName("Multi-Factor GBM TS Implied Model PCA Calibration")]
    public class MultiGBMAssetPriceTSPCACalibrationImplied : PCACalibration<MultiGBMAssetPriceTSModelImplied>, IModelCalibration
    {
        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceTSModelImplied);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<MultiGBMAssetPriceTSModelImplied>(calibrationData, priceModel, output);

            // Run the main calibration routine in the base class.
            IPriceFactor priceFactorToEvolve = CalibrationHelper.GetPriceFactor(priceModel);
            if (TryCalibrate(calibrationData, model, priceFactorToEvolve, errorList))
                GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        protected override void SetDriftAndVolStatistics(ICalibrationData calibrationData, MultiGBMAssetPriceTSModelImplied modelParameters, double vol, int statisticsSetIndex)
        {
            // We do nothing.
            // The vol calibration is implied, so is handled by the bootstrapper.
            // We take Risk_Neutral_Drift = Yes for the implied model, so no drift required.
        }

        /// <inheritdoc />
        protected override void SetDriftAndVolHistorical(ICalibrationData calibrationData, MultiGBMAssetPriceTSModelImplied modelParameters, double drift, double vol)
        {
            // We do nothing.
            // The vol calibration is implied, so is handled by the bootstrapper.
            // We take Risk_Neutral_Drift = Yes for the implied model, so no drift required.
        }

        /// <inheritdoc />
        protected override MultiGBMAssetPriceTSModelImplied GetOrRegisterModelParameters(FactorID id, PriceModelList priceModelList)
        {
            return priceModelList.Register<MultiGBMAssetPriceTSModelImplied>(id);
        }
    }
}
