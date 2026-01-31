using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Drawing;
using System.Globalization;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// MLE calibration of PCA interest rate models without data cleaning or data filling.
    /// </summary>
    [DisplayName("PCA Interest Rate Calibration")]
    public class PCAInterestRateCalibration : IModelCalibration
    {
        // Any tenor with less than 50% data (compared to other tenors) will be ignored.
        private const double DataThresholdForValidTenor = 0.5;

        // Extracted market data.
        private TimeSeriesStructure<Period> fTtermStructure;

        /// <summary>
        /// Initializes a new instance of the PcaInterestRateMleCalibration class and sets the defaults.
        /// </summary>
        public PCAInterestRateCalibration()
        {
            // Use MLE by default.
            Calibration_Method = CalibrationMethod.MLE;
            MLE_Parameters = new PCAInterestRateMLEParameters();

            Number_Of_PCA_Factors = 3;
            Matrix_Type = PCSource.Correlation;
            Distribution_Type = ProbabilityDistribution.Lognormal;
        }

        /// <summary>
        /// Calibration method can be statistics or MLE.
        /// </summary>
        public CalibrationMethod Calibration_Method { get; set; }

        /// <summary>
        /// The number of PCA factors.
        /// </summary>
        public int Number_Of_PCA_Factors { get; set; }

        /// <summary>
        /// Model distribution type, i.e., normal or lognormal.
        /// </summary>
        public ProbabilityDistribution Distribution_Type { get; set; }

        /// <summary>
        /// User-defined property for whether to use covariance or correlation to calculate beta.
        /// </summary>
        public PCSource Matrix_Type { get; set; }

        /// <summary>
        /// User-defined property to set the drift of the model.
        /// </summary>
        public DriftMethod Rate_Drift_Model { get; set; }

        /// <summary>
        /// Parameters for the MLE method.
        /// </summary>
        public PCAInterestRateMLEParameters MLE_Parameters { get; set; }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="priceModel">The price model to calibrate.</param>
        /// <param name="output">Calibration output in an xml format.</param>
        /// <param name="errorList">List of infos, warnings and errors.</param>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            Debug.Assert(errorList != null, "Errorlist cannot be null.");

            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new PCAInterestRateStatisticsCalibration
                {
                    Number_Of_PCA_Factors = Number_Of_PCA_Factors,
                    Matrix_Type = Matrix_Type,
                    Rate_Drift_Model = Rate_Drift_Model,
                    Distribution_Type = Distribution_Type
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            // Validate the data against the model.
            PCAInterestRateModel pcaModel = CalibrationHelper.ValidateWithoutStatistics<PCAInterestRateModel>(calibrationData, priceModel, output);

            // Validate properties. Fail if any property not properly set.
            if (!ValidateProperties(errorList, calibrationData.ArchiveStartDate, calibrationData.ArchiveEndDate,
                priceModel.ToString()))
            {
                return;
            }

            pcaModel.Princ_Comp_Source = Matrix_Type;
            pcaModel.Distribution_Type = Distribution_Type;
            pcaModel.Rate_Drift_Model = Rate_Drift_Model;

            // Get the price factor for the model: we need the historical values corresponding to this factor.
            IPriceFactor priceFactor = pcaModel.GetPriceFactor();
            if (priceFactor == null)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Calibration of {0} requires the model's price factor to be present in the market data.", priceModel));
                return;
            }

            // Next, get the historical values corresponding to the response variables (the zero rates of various allTenors).
            if (!GenericDataRetrieval.TryGetHistoricalRateTermStructure(calibrationData, priceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out fTtermStructure))
            {
                return;
            }

            Period[] allHistoricalTenors = fTtermStructure.TimeSeriesIDs;
            double[][] allHistoricalZeroRates = fTtermStructure.Values;
            double[][] allHistoricalZeroRateAtHorizon = fTtermStructure.ValuesAtHorizon;

            int numTenors = fTtermStructure.NumTimeSeries;

            if (numTenors < Number_Of_PCA_Factors)
            {
                errorList.Add(ErrorLevel.Warning, string.Format(
                                  "Calibration of {0} requires number of all tenors to be greater or equal to the number of PCA factors. Setting to the number of tenors.", priceModel));

                Number_Of_PCA_Factors = numTenors;
            }

            // Select only a set of tenors for calibration based on tenor being greater or equal to the minimum tenor used for calibration
            // and sufficient historical data at each tenor.
            Period[] unselectedTenorsToExtrapolate = null;
            Period[] selectedTenors;
            double[][] selectedHistoricalZeroRates;
            double[][] selectedHistoricalZeroRatesAtHorizon;

            double? minTenorParameter = MLE_Parameters.GetMinTenor();
            if (minTenorParameter.HasValue)
            {
                List<int> selectedTenorIndices;
                Period[] invalidSelectedTenors;
                bool hadValidTenors = CalibrationHelper.TryGetSelectedTermStructureWithoutInvalidTenors(allHistoricalTenors,
                                                                                      allHistoricalZeroRates,
                                                                                      minTenorParameter.Value,
                                                                                      DataThresholdForValidTenor,
                                                                                      errorList,
                                                                                      out unselectedTenorsToExtrapolate,
                                                                                      out invalidSelectedTenors,
                                                                                      out selectedTenorIndices,
                                                                                      out selectedTenors,
                                                                                      out selectedHistoricalZeroRates);

                // No valid tenor. Failing calibration.
                if (!hadValidTenors)
                    return;

                selectedHistoricalZeroRatesAtHorizon = new double[selectedTenorIndices.Count][];
                for (int index = 0; index < selectedTenorIndices.Count; ++index)
                {
                    int selectedTenorIndex = selectedTenorIndices[index];
                    selectedHistoricalZeroRatesAtHorizon[index] = allHistoricalZeroRateAtHorizon[selectedTenorIndex];
                }
            }
            else
            {
                selectedTenors = allHistoricalTenors;
                selectedHistoricalZeroRates = allHistoricalZeroRates;
                selectedHistoricalZeroRatesAtHorizon = allHistoricalZeroRateAtHorizon;
            }

            // Output allTenors calibrating to.
            string tenorsToCalibrateStr = "{";
            int selectedLength = selectedTenors.Length;
            for (int index = 0; index < selectedLength; ++index)
            {
                tenorsToCalibrateStr += selectedTenors[index] + (index == selectedLength - 1 ? "}" : ", ");
            }

            errorList.Add(ErrorLevel.Info, string.Format("Calibrating to Tenors {0}.", tenorsToCalibrateStr));

            // Calibrate the unknown parameters (alpha, sigma, theta) for selected allTenors.
            double alpha;
            double[] selectedTheta;
            double[] selectedSigma;
            SymmetricMatrix selectedCovarianceMatrix;

            CalibratePcaMleAlphaSigmaTheta(selectedHistoricalZeroRates, selectedHistoricalZeroRatesAtHorizon,
                                           fTtermStructure.DeltaT, errorList,
                                           out alpha, out selectedSigma, out selectedTheta, out selectedCovarianceMatrix);

            bool isCorrelationMatrix = Matrix_Type == PCSource.Correlation;
            if (isCorrelationMatrix)
            {
                // Compute correlation matrix if required.
                for (int indexA = 0; indexA < selectedLength; ++indexA)
                {
                    for (int indexB = 0; indexB <= indexA; ++indexB)
                    {
                        selectedCovarianceMatrix[indexA, indexB] = selectedCovarianceMatrix[indexA, indexB] /
                                                                   (selectedSigma[indexA] * selectedSigma[indexB]);
                    }
                }
            }

            // Extrapolate model parameters
            double[] theta;
            double[] sigma;
            SymmetricMatrix covarianceMatrix;

            Period[] allValidTenors;
            if (!ExtrapolateModelParameters(unselectedTenorsToExtrapolate, selectedTenors, selectedSigma, selectedTheta, selectedCovarianceMatrix,
                                            errorList, isCorrelationMatrix, out allValidTenors, out sigma, out theta, out covarianceMatrix))
                return;

            int numValidTenors = allValidTenors.Length;

            // Set (alpha, sigma, theta) on the model.
            pcaModel.Yield_Volatility.Clear();
            pcaModel.Historical_Yield.Clear();
            for (int tenorIndex = 0; tenorIndex < numValidTenors; ++tenorIndex)
            {
                double tenor = allValidTenors[tenorIndex];
                pcaModel.Yield_Volatility[tenor] = sigma[tenorIndex];
                pcaModel.Historical_Yield[tenor] = theta[tenorIndex];
            }

            pcaModel.Reversion_Speed = alpha;

            // PCA with all historical data
            Matrix eigenVectors;
            Vector eigenValues;
            CalcUtils.ComputePcaCovarianceMethod(covarianceMatrix, out eigenValues, out eigenVectors);

            // Set eigenvectors and eigenvalues model parameters from a PCA with all historical data.
            pcaModel.Eigenvectors.Clear();
            for (int driverIndex = 0; driverIndex < Number_Of_PCA_Factors; ++driverIndex)
            {
                // Check that eigenvalue is not negative - since eigenvalues are sorted, any negative eigenvalue will be followed by negative eigenvalues.
                if (eigenValues[driverIndex] < 0)
                {
                    errorList.Add(ErrorLevel.Warning,
                                  string.Format(
                                      "Using only {0} PCA factors since lower PCA factors have negative eigenvalues (a possible shortcoming of the QL algorithm even on a PSD matrix).",
                                      driverIndex));
                    break;
                }

                var eigen = new Eigen<RateCurve> { Eigenvalue = eigenValues[driverIndex] };

                for (int tenorIndex = 0; tenorIndex < numValidTenors; ++tenorIndex)
                {
                    double tenor = allValidTenors[tenorIndex];
                    eigen.Eigenvector[tenor] = eigenVectors[tenorIndex, driverIndex];
                }

                pcaModel.Eigenvectors.Add(eigen);
            }

            // Write calibration data to an XML stream.
            GetUpdatedCurves(pcaModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(PCAInterestRateModel);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">The calibrated price model.</param>
        /// <param name="output">The calibration xml output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            PCAInterestRateModel pcaModel = (PCAInterestRateModel)priceModel;

            var graphData = new ResultSeriesList { fXAxis = { fLabel = "Time (Years)" }, fYAxis = { fLabel = "Value" } };
            graphData.Add(new ResultSeries("Historical Yield", LineStyle.Solid, 2, Color.Blue, pcaModel.Historical_Yield.Clone()));
            graphData.Add(new ResultSeries("Yield Volatility", LineStyle.Solid, 2, Color.Black, pcaModel.Yield_Volatility.Clone()));

            for (int k = 0; k < pcaModel.GetDimension(); ++k)
                graphData.Add(new ResultSeries("Eigenvector " + (k + 1), LineStyle.Dashed, pcaModel.Eigenvectors[k].Eigenvector.Clone()));

            graphData.ToXml(output);
            CalibrationHelper.WriteModelParametersToXml(priceModel, output);
        }

        /// <summary>
        /// Determine the optimal alpha, theta, sigma parameters under the MLE approach for the PCA interest rate model.
        /// </summary>
        /// <param name="timeSeries">The historical data a two-dimensional matrix.</param>
        /// <param name="timeSeriesAtHorizon">Two-dimensional timeseries at horizon.</param>
        /// <param name="horizon">Horizon of returns.</param>
        /// <param name="holidayCalendar">Holiday calendar.</param>
        /// <param name="errorList">List of infos, warnings and errors.</param>
        /// <param name="alpha">The mean-reversion speed.</param>
        /// <param name="sigma">The mean-reversion volatility.</param>
        /// <param name="theta">The mean-reversion level.</param>
        /// <param name="covarMatrix">The covariance matrix.</param>
        /// <returns>Return true if model parameters successfully computed.</returns>
        /// <remarks>We always use the exact MLE solution.</remarks>
        private void CalibratePcaMleAlphaSigmaTheta(double[][] timeSeries, double[][] timeSeriesAtHorizon,
                                                    double deltaT, ErrorList errorList,
                                                    out double alpha, out double[] sigma, out double[] theta, out SymmetricMatrix covarMatrix)
        {
            ReturnMethod returnMethod = Distribution_Type == ProbabilityDistribution.Lognormal ? ReturnMethod.Log : ReturnMethod.Diff;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(returnMethod,
                                                                                                   MLE_Parameters.Reversion_Speed_Lower_Bound,
                                                                                                   MLE_Parameters.GetMeanReversionSpeedUpperBound(),
                                                                                                   MLE_Parameters.GetYieldVolatilityUpperBound(),
                                                                                                   MLE_Parameters.GetFixedMeanReversionSpeed(),
                                                                                                   true)
            {
                OptimisationParameters = MLE_Parameters.Exact_Solution_Optimisation_Parameters
            };

            CorrelatedMeanReversionStatisticsScalarAlpha meanReversionStatistics;
            ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(timeSeries, timeSeriesAtHorizon, deltaT, calculationParameters,
                                                                               errorList, out meanReversionStatistics);

            alpha = meanReversionStatistics.MeanReversionSpeed;
            sigma = meanReversionStatistics.MeanReversionVols;
            theta = meanReversionStatistics.LongRunMeans;
            covarMatrix = meanReversionStatistics.CovarianceMatrix;
        }

        /// <summary>
        /// Extrapolate model parameters based on weighted average to calculate sigma and theta parameters for small tenors 
        /// (i.e., less than the minimum tenor used in calibration).
        /// We use weighted average rather than extrapolotion because of the unstable behaviour for small tenors.
        /// </summary>
        /// <param name="unselectedLowerTenors">Unselected tenors to extrapolate.</param>
        /// <param name="selectedTenors">The selected tenors for calibration.</param>
        /// <param name="selectedSigma">The selected mean-reversion volatilities.</param>
        /// <param name="selectedTheta">The selected mean-reversion levels.</param>
        /// <param name="selectedCovarMatrix">The selected covariance matrix.</param>
        /// <param name="errorList">List of infos, warnings and errors.</param>
        /// <param name="isCorrelationMatrix">Flag to use a covariance matrix or a correlation matrix.</param>
        /// <param name="allValidTenors">The valid tenors (selected and unselected).</param>
        /// <param name="sigma">The mean-reversion volatilities.</param>
        /// <param name="theta">The mean-reversion levels.</param>
        /// <param name="covarMatrix">The extrapolated covariance matrix.</param>
        /// <returns>Return true if extrapolation is successful.</returns>
        private bool ExtrapolateModelParameters(Period[] unselectedLowerTenors, Period[] selectedTenors, double[] selectedSigma,
                                                double[] selectedTheta, SymmetricMatrix selectedCovarMatrix, ErrorList errorList, bool isCorrelationMatrix,
                                                out Period[] allValidTenors, out double[] sigma, out double[] theta, out SymmetricMatrix covarMatrix)
        {
            if (unselectedLowerTenors == null || unselectedLowerTenors.Length == 0)
            {
                allValidTenors = selectedTenors;
                sigma = selectedSigma;
                theta = selectedTheta;
                covarMatrix = selectedCovarMatrix;
                return true;
            }

            int numSelectedTenors = selectedTenors.Length;
            int numLowerTenorToExtrapolate = unselectedLowerTenors.Length;
            allValidTenors = new Period[numSelectedTenors + numLowerTenorToExtrapolate];

            // Array of all tenors.
            for (int index = 0; index < numLowerTenorToExtrapolate; ++index)
                allValidTenors[index] = unselectedLowerTenors[index];

            for (int index = 0; index < numSelectedTenors; ++index)
                allValidTenors[numLowerTenorToExtrapolate + index] = selectedTenors[index];

            sigma = GetFlatLeftExtrapolation(selectedSigma, numLowerTenorToExtrapolate);
            theta = GetFlatLeftExtrapolation(selectedTheta, numLowerTenorToExtrapolate);

            return CorrelationExtrapolationHelper.ExtrapolateFlatLeftCovarianceMatrix(unselectedLowerTenors, selectedTenors, selectedCovarMatrix,
                                                                                      sigma, errorList, isCorrelationMatrix, out covarMatrix);
        }

        /// <summary>
        /// Extrapolate the curve backward.
        /// </summary>
        private double[] GetFlatLeftExtrapolation(double[] selectedTenors, int numLowerTenorsToLeftExtrapolate)
        {
            int numSelectedTenors = selectedTenors.Length;
            var extrapolatedValues = new double[numLowerTenorsToLeftExtrapolate + numSelectedTenors];

            // Flat left interpolation.
            for (int i = 0; i < numLowerTenorsToLeftExtrapolate; i++)
                extrapolatedValues[i] = selectedTenors[0];

            for (int i = 0; i < selectedTenors.Length; i++)
                extrapolatedValues[numLowerTenorsToLeftExtrapolate + i] = selectedTenors[i];

            return extrapolatedValues;
        }

        /// <summary>
        /// Properties validation. Provides user with a list of errors and warnings about all properties.
        /// </summary>
        private bool ValidateProperties(ErrorList errorList, DateTime archiveStartDate, DateTime archiveEndDate, string priceModelID)
        {
            bool validProperties = true;

            if (Number_Of_PCA_Factors <= 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("{0} calibration - the number of PCA factors must be greater or equal to 1.", priceModelID));
                validProperties = false;
            }

            var calibrationEndDate = MLE_Parameters.Data_Retrieval_Parameters.GetEndDate();
            if (calibrationEndDate.HasValue && calibrationEndDate < archiveStartDate)
            {
                errorList.Add(ErrorLevel.Error,
                    string.Format("The {0} calibration end date, {1}, must be on or after the archive start date, {2}",
                    priceModelID, calibrationEndDate.ToString(), archiveStartDate.ToString(TDate.DateFormat,
                        CultureInfo.InvariantCulture)));
                validProperties = false;
            }

            var calibrationStartDate = MLE_Parameters.Data_Retrieval_Parameters.GetStartDate();
            if (calibrationStartDate.HasValue && calibrationStartDate > archiveEndDate)
            {
                errorList.Add(ErrorLevel.Error,
                    string.Format("The {0} calibration start date, {1}, must be on or before the archive end date, {2}",
                    priceModelID, calibrationStartDate.ToString(), archiveEndDate.ToString(TDate.DateFormat,
                        CultureInfo.InvariantCulture)));
                validProperties = false;
            }

            if (!string.IsNullOrEmpty(MLE_Parameters.Min_Tenor) && new Period(MLE_Parameters.Min_Tenor) < 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                                  "The minimum tenor threshold for calibration of the {0} is not correctly set. It must be greater or equal 0D.", priceModelID));
                validProperties = false;
            }

            if (MLE_Parameters.Reversion_Speed_Upper_Bound.Exists && MLE_Parameters.Reversion_Speed_Lower_Bound >= MLE_Parameters.Reversion_Speed_Upper_Bound.Value)
            {
                errorList.Add(ErrorLevel.Error, string.Format(
                                  "Calibration of {0} - the lower bound on the reversion speed must be lower than the upper bound.", priceModelID));
                validProperties = false;
            }

            validProperties &= MLE_Parameters.Data_Retrieval_Parameters.ValidateProperties(errorList);

            return validProperties;
        }
    }

    /// <summary>
    /// Parameters for the MLE method.
    /// </summary>
    public class PCAInterestRateMLEParameters : NestedPresentableObject
    {
        private Period? fMinTenor;
        private string fMinTenorString;

        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public PCAInterestRateMLEParameters()
        {
            // Data cleaning (to remove outliers) used by default.
            Data_Retrieval_Parameters = new CalibratorDataRetrievalParameters<Period>
            {
                Data_Cleaning_Methods =
                {
                    new DataCleaningMethodWrapper<Period>
                    {
                        Method_Name = Property.DisplayName(typeof(DataCleaningZScoreMethod<Period>)),
                        Parameters  = new DataCleaningZScoreMethod<Period>()
                    }
                }
            };

            Reversion_Speed_Fixed = new DoubleNull();
            Reversion_Speed_Lower_Bound = 0.1;                                  // A lower reversion speed would be unstable when transforming a log-normal theta.
            Reversion_Speed_Upper_Bound = new DoubleNull(2.0, true);            // A higher reversion speed is unlikely and observed for very short time series.
            Yield_Volatility_Upper_Bound = new DoubleNull();
            Min_Tenor = "3M";                                 // The minimum tenor set to 3 months by default.

            // Use approximation by default.
            Exact_Solution_Optimisation_Parameters = new MeanReversionStatisticsOptimisationParameters();
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<Period> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// User-defined floor on zero-rate historical data. Optional.
        /// </summary>
        public string Min_Tenor
        {
            get
            {
                return fMinTenorString;
            }
            set
            {
                fMinTenorString = value;
                fMinTenor = string.IsNullOrEmpty(fMinTenorString) ? (Period?)null : new Period(fMinTenorString);
            }
        }

        /// <summary>
        /// User-defined property for an explicit alpha. Optional.
        /// </summary>
        public DoubleNull Reversion_Speed_Fixed { get; set; }

        /// <summary>
        /// User-defined property for a lower bound on alpha. We need to ensure a positive lower bound.
        /// </summary>
        public double Reversion_Speed_Lower_Bound { get; set; }

        /// <summary>
        /// User-defined property for an uppder bound on alpha. Optional.
        /// </summary>
        public DoubleNull Reversion_Speed_Upper_Bound { get; set; }

        /// <summary>
        /// User-defined property for a lower bound on alpha. Optional.
        /// </summary>
        public DoubleNull Yield_Volatility_Upper_Bound { get; set; }

        /// <summary>
        /// Optimisation parameters.
        /// </summary>
        public MeanReversionStatisticsOptimisationParameters Exact_Solution_Optimisation_Parameters { get; set; }

        /// <summary>
        /// Returns the min tenor for calibration.
        /// </summary>
        public Period? GetMinTenor()
        {
            return fMinTenor;
        }

        /// <summary>
        /// Return null if the value does not exist.
        /// </summary>
        public double? GetFixedMeanReversionSpeed()
        {
            return Reversion_Speed_Fixed.Exists ? Reversion_Speed_Fixed.Value : (double?)null;
        }

        /// <summary>
        /// Return null if the value does not exist.
        /// </summary>
        public double? GetMeanReversionSpeedUpperBound()
        {
            return Reversion_Speed_Upper_Bound.Exists ? Reversion_Speed_Upper_Bound.Value : (double?)null;
        }

        /// <summary>
        /// Return null if the value does not exist.
        /// </summary>
        public double? GetYieldVolatilityUpperBound()
        {
            return Yield_Volatility_Upper_Bound.Exists ? Yield_Volatility_Upper_Bound.Value : (double?)null;
        }
    }
}