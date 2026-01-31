/// <author>
/// Alastair Wilkins, Nathalie Antunes, Phil Koop, Tom Judd, Ian Rootes
/// </author>
/// <owner>
/// Alastair Wilkins
/// </owner>
/// <summary>
/// Calibration that calculates risk factor correlations from the statistics set.
/// </summary>
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Interface that price factor models implement to describe the linear relationship 
    /// between risk factor returns and the returns on the factors in the statistics set.
    /// Used in the correlations calibration.
    /// </summary>
    public interface ICorrelatedNormalCoefficients
    {
        /// <summary>
        /// For each correlated Normal in modelCorrelatedNormals, adds a row to coefficientMatrix giving the
        /// coeeficients of that correlated normal with respect the values in the statistics set.
        /// </summary>
        /// <param name="statistics">Statistics set.</param>
        /// <param name="priceFactorType">Type of the underlying price factor.</param>
        /// <param name="modelCorrelatedNormals">The model's correlated normals.</param>
        /// <param name="modelCorrelatedNormalIndexes">Indexes of the model's correlated normals in the collection of all correlated normals.</param>
        /// <param name="coefficientMatrix">Dictionary mapping correlated normal index to matrix row, where matrix row is a dictionary mapping statistics set index to coefficient.</param>
        void CalculateCoefficients(IStatistics statistics, Type priceFactorType, RiskFactorList modelCorrelatedNormals, List<int> modelCorrelatedNormalIndexes, CoefficientMatrix coefficientMatrix);
    }

    /// <summary>
    /// Sparse matrix row: a mapping of statistics index to coefficient.
    /// </summary>
    public class CoefficientRow : Dictionary<int, double>
    {
        /// <summary>
        /// Create a new dictionary mapping statistics index to coefficient value.
        /// </summary>
        /// <param name="primaryRow">Index of first risk factor in group.</param>
        public CoefficientRow(int primaryRow)
        {
            PrimaryRow = primaryRow;
        }

        public int PrimaryRow
        {
            get; set;
        }
    }

    /// <summary>
    /// Sparse matrix: a mapping of risk factor index to sparse matrix row.
    /// </summary>
    public class CoefficientMatrix : Dictionary<int, CoefficientRow>
    {
    }

    /// <summary>
    /// Standard calculation of risk factor correlations from price factor correlations in correlations file
    /// </summary>
    [DisplayName("Correlations Calibration")]
    public class CorrelationsCalibration : IModelCalibration
    {
        /// <summary>
        /// Whether to use averaging or a single price factor for one-factor price models.
        /// </summary>
        public YesNo Use_Averaging
        {
            get; set;
        }

        /// <summary>
        /// Calibration calculation that calculates risk factor correlations.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // priceModel will be null
            CalibrationHelper.Validate(calibrationData, output);

            Calculate(calibrationData.MarketData, calibrationData.Statistics, errorList);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SparseMatrix);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model
        /// </summary>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Generic calculation of coefficients for single risk factor model.
        /// </summary>
        public static void SingleFactor(PriceModel priceModel, IStatistics statistics, Type priceFactorType, int correlatedNormalIndex, CoefficientMatrix coefficientMatrix, bool useAveraging)
        {
            List<int> indexes;
            List<string> points;

            if (!statistics.FindEntries(priceFactorType.Name, priceModel.fID, out indexes, out points))
            {
                // No statistics for priceFactor
                return;
            }

            // Create a new row for this correlatedNormal
            CoefficientRow row = new CoefficientRow(correlatedNormalIndex);
            coefficientMatrix[correlatedNormalIndex] = row;
            if (useAveraging)
            {
                foreach (int i in indexes)
                {
                    row[i] = 1;
                }
            }
            else
            {
                // use only the first price factor value
                row[indexes[0]] = 1;
            }
        }

        /// <summary>
        /// Calculate coefficients for a PCA model.
        /// </summary>
        public static void CalculatePCACoefficients(int dimensions, List<int> factorIndexes, List<int> modelCorrelatedNormalIndexes, Func<int, double> getEigenvalue, Func<int, int, double> getEigenvectorEntry, double[] volatilities, CoefficientMatrix coefficientMatrix)
        {
            int numFactors = Math.Min(dimensions, factorIndexes.Count);
            int primaryRow = -1;
            for (int factor = 0; factor < numFactors; ++factor)
            {
                if (primaryRow < 0)
                    primaryRow = modelCorrelatedNormalIndexes[factor];

                double overRootEigenvalue = getEigenvalue(factor);
                overRootEigenvalue = overRootEigenvalue > 0.0 ? 1.0 / Math.Sqrt(overRootEigenvalue) : 0.0;

                // Set the coefficients for the current PCA factor
                CoefficientRow row = new CoefficientRow(primaryRow);
                for (int i = 0; i < factorIndexes.Count; ++i)
                    row[factorIndexes[i]] = volatilities[i] * overRootEigenvalue * getEigenvectorEntry(factor, i);

                coefficientMatrix[modelCorrelatedNormalIndexes[factor]] = row;
            }
        }

        /// <summary>
        /// Generic linear calculation of risk factor correlations using a risk factor/price factor coefficient matrix (A), as per Theory Guide.
        /// </summary>
        private void Calculate(MarketDataContainer marketData, IStatistics statistics, ErrorList errors)
        {
            // Calculate the coefficient matrix (A)
            CoefficientMatrix coefficientMatrix = CalculateCoefficientMatrix(marketData, statistics, errors);

            // Calculate normalisation factors for the coefficient matrix
            CoefficientRow normaliseFactors = new CoefficientRow(0);
            foreach (int k in coefficientMatrix.Keys)
            {
                CoefficientRow row = coefficientMatrix[k];

                double normaliseFactor = 0.0;
                normaliseFactor = statistics.Correlation(row.Keys, row, row.Keys, row);

                normaliseFactor = normaliseFactor > 0.0 ? 1.0 / Math.Sqrt(normaliseFactor) : 0.0;
                normaliseFactors[k] = normaliseFactor;
            }

            // The risk factor correlation matrix is ACA^T, where C is the price factor correlation matrix
            foreach (int k in coefficientMatrix.Keys)
            {
                foreach (int l in coefficientMatrix.Keys)
                {
                    // The matrix is symmetric so quit when inner loops reaches the diagonal
                    if (l == k)
                    {
                        break;
                    }

                    double riskFactorCorrelation = 0.0;

                    // Correlations for PCA components against another PCA component of the same group are 0
                    if (coefficientMatrix[k].PrimaryRow != coefficientMatrix[l].PrimaryRow)
                    {
                        riskFactorCorrelation = statistics.Correlation(coefficientMatrix[k].Keys, coefficientMatrix[k], coefficientMatrix[l].Keys, coefficientMatrix[l]);

                        riskFactorCorrelation *= normaliseFactors[k] * normaliseFactors[l];

                        if (double.IsNaN(riskFactorCorrelation))
                            riskFactorCorrelation = 0.0;
                    }

                    // Any correlations with an absolute value calculated of less than TINY can be assumed to be zero.
                    // They are not actually calculated as zero due to errors at the limit of double precision.
                    // e.g. PCA1 of currency 1 against PCA2 of currency 2 where the time series for currencies 1 and 2 are identical.
                    if (Math.Abs(riskFactorCorrelation) > CalcUtils.TINY)
                        marketData.SetCorrelation(k, l, (float)riskFactorCorrelation);
                    else
                        marketData.SetCorrelation(k, l, 0.0F);
                }
            }
        }

        /// <summary>
        /// Calculate the correlated normal/price factor coefficient matrix from price models.
        /// </summary>
        private CoefficientMatrix CalculateCoefficientMatrix(MarketDataContainer marketData, IStatistics statistics, ErrorList errors)
        {
            CoefficientMatrix coefficientMatrix = new CoefficientMatrix();

            PriceModelList priceModels = marketData.PriceModels;
            var finder = CalibrationConfig.FindFactorsForModels(marketData.PriceFactors);
            foreach (PriceModel priceModel in priceModels)
            {
                priceModel.SetPriceFactor(finder.FindFactorForModel(priceModel));
                IPriceFactor modelPriceFactor = priceModel.GetPriceFactor();
                if (modelPriceFactor == null)
                    continue;

                // Get model's risk factors
                RiskFactorList modelRiskFactors = new RiskFactorList();
                priceModel.AttachProcesses(modelRiskFactors);   // Don't use GetRiskFactors() because we want the side-effects of attaching processes.
                CorrelatedNormals correlatedNormalDict = modelRiskFactors.CorrelatedNormals;
                if (correlatedNormalDict.Count == 0)
                    continue;

                // Indexes of risk factors in collection of all risk factors
                List<int> correlatedNormalDictIndexes = new List<int>(correlatedNormalDict.Count);

                // Have statistics for at least one of the model risk factors
                bool haveStatistics = false;

                foreach (CorrelatedNormal correlatedNormal in correlatedNormalDict)
                {
                    int correlatedNormalIndex = marketData.RiskFactors.CorrelatedNormals.GetOrAddIndex(correlatedNormal.Key);
                    correlatedNormalDictIndexes.Add(correlatedNormalIndex);

                    // Get index of risk factor in statistics set.
                    // Returns -1 if no statistics for this risk factor.
                    int statisticsIndex = statistics.GetIndex(correlatedNormal.Key, null);
                    if (statisticsIndex >= 0)
                    {
                        // There is one coefficient relating the risk factor to its entry in the statistics set
                        CoefficientRow row = new CoefficientRow(correlatedNormalIndex);
                        coefficientMatrix[correlatedNormalIndex] = row;
                        row[statisticsIndex] = 1.0;
                        haveStatistics = true;
                    }
                }

                if (haveStatistics)
                    continue;

                // Calculate risk factor/price factor coefficients when
                // the model implements IRiskFactorCoefficients or has a single risk factor.
                var correlatedNormalCoefficients = priceModel as ICorrelatedNormalCoefficients;
                if (correlatedNormalCoefficients != null)
                {
                    try
                    {
                        correlatedNormalCoefficients.CalculateCoefficients(statistics, modelPriceFactor.GetType(), modelRiskFactors, correlatedNormalDictIndexes, coefficientMatrix);
                    }
                    catch (CalibrationException exception)
                    {
                        errors.Add(exception);
                    }
                }
                else if (modelRiskFactors.Count == 1)
                {
                    SingleFactor(priceModel, statistics, modelPriceFactor.GetType(), correlatedNormalDictIndexes[0], coefficientMatrix, Use_Averaging == YesNo.Yes);
                }
            }

            return coefficientMatrix;
        }
    }
}