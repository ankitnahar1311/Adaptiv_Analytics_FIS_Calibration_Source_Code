using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Stores the name of a driver.
    /// </summary>
    [Serializable]
    public class DriverName : DisplayableString
    {
    }

    /// <summary>
    /// Abstract class for regression calibration of a multi factor GBM model.
    /// </summary>
    /// <typeparam name="TModelParameters">The type of model (or model parameters object) we are calibrating.</typeparam>
    public abstract class RegressionCalibration<TModelParameters> : StatisticsBasedCalibration
        where TModelParameters : BaseFactor, IMultiFactorModelParameters
    {
        // Extracted drivers time series.
        protected readonly List<TimeSeriesStructure<string>> fDriverTimeSeries;

        // Extracted name time series.
        protected TimeSeriesStructure<string> fNameTimeSeries;

        /// <summary>
        /// Default constructor.
        /// </summary>
        protected RegressionCalibration()
        {
            Drivers                = new DisplayableList<DriverName>();
            Significance_Threshold = CalcUtils.DefaultMinProcessWeight;
            fDriverTimeSeries      = new List<TimeSeriesStructure<string>>();
        }

        /// <summary>
        /// List of driving factor IDs.
        /// </summary>
        public DisplayableList<DriverName> Drivers
        {
            get;
            set;
        }

        /// <summary>
        /// Value below which we round weights to zero.
        /// </summary>
        public double Significance_Threshold
        {
            get;
            set;
        }

        /// <summary>
        /// Worker calibration routine.
        /// </summary>
        public bool Calibrate(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            if (Use_Pre_Computed_Statistics == YesNo.No)
            {
                if (!TryGetNameAndDriversTimeSeries(calibrationData, priceFactorToEvolve, errorList))
                    return false;
            }

            return Calculate(calibrationData, modelParameters, priceFactorToEvolve, errorList);
        }

        /// <summary>
        /// Performs a shallow copy of the properties on this object to the target object.
        /// </summary>
        protected void ShallowCopyTo(RegressionCalibration<TModelParameters> target)
        {
            base.ShallowCopyTo(target);
            target.Drivers                = this.Drivers;
            target.Significance_Threshold = this.Significance_Threshold;
        }

        /// <summary>
        /// Returns whether the calibration uses an idiosyncratic factor.
        /// </summary>
        protected virtual YesNo GetIdiosyncraticFactorChoice()
        {
            return YesNo.Yes;
        }

        /// <summary>
        /// Calculate the correlation matrix of the name you are calibrating on,
        /// and the correlation matrix of the drivers using correlations from the statistics set.            
        /// </summary>
        /// <param name="calibrationData">Market Data and statistics set.</param>
        /// <param name="modelParameters">Multi-factor model (or model parameters object) to calibrate.</param>
        /// <param name="priceFactorToEvolve">The price factor that is being evolved.</param>
        /// <param name="errorList">List of info, warnings and errors.</param>
        /// <param name="nameDriverCorrelations">Correlations of name with drivers.</param>
        /// <param name="driverCorrelations">Correlation matrix of the drivers.</param>
        protected abstract void GetNameDriverStatisticsCorrelations(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList,
                                                                    out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations);

        /// <summary>
        /// Calculate the correlation matrix of the name you are calibrating on,
        /// and the correlation matrix of the drivers using correlations calculated from historical data.            
        /// </summary>
        /// <param name="calibrationData">Market Data and statistics set.</param>
        /// <param name="modelParameters">Multi-factor model (or model parameters object) to calibrate.</param>
        /// <param name="priceFactorToEvolve">The price factor that is being evolved.</param>
        /// <param name="errorList">List of info, warnings and errors.</param>
        /// <param name="nameDriverCorrelations">Correlations of name with drivers.</param>
        /// <param name="driverCorrelations">Correlation matrix of the drivers.</param>
        protected abstract bool GetNameDriverHistoricalCorrelations(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList, out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations);

        /// <summary>
        /// Calculate the model parameters.
        /// </summary>
        /// <param name="calibrationData">Statistics set.</param>
        /// <param name="modelParameters">Multi-factor model (or model parameters object) to calibrate.</param>
        /// <param name="priceFactorToEvolve">The price factor that is being evolved.</param>
        /// <param name="errorList">List of calibration warnings and errors.</param>
        private bool Calculate(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            int numberOfDrivers = Drivers.Count;
            if (numberOfDrivers == 0)
            {
                // No drivers specified, nothing to calibrate
                return false;
            }

            Vector nameDriverCorrelations;
            SymmetricMatrix driverCorrelations;

            // Calculate the correlation matrix of the name, the correlation matrix of the drivers
            if (Use_Pre_Computed_Statistics == YesNo.Yes)
                GetNameDriverStatisticsCorrelations(calibrationData, modelParameters, priceFactorToEvolve, errorList, out nameDriverCorrelations, out driverCorrelations);
            else if (!GetNameDriverHistoricalCorrelations(calibrationData, modelParameters, priceFactorToEvolve, errorList, out nameDriverCorrelations, out driverCorrelations))
                return false;

            // Solve the multivariate linear regression problem
            Vector beta;
            RegressionCalibrationHelper.ComputeParametersOfMultivariateRegression(out beta, driverCorrelations, nameDriverCorrelations, numberOfDrivers, errorList);
            if (beta == null)
            {
                return false;
            }

            double epsilon;
            MultiFactorCalibrationHelper.ObtainRoundedNormalisedRegressionWeights(Significance_Threshold, 1.0, GetIdiosyncraticFactorChoice(), driverCorrelations, beta, out epsilon);

            // set model parameters of the interface
            var model = (IMultiFactorModelParameters)modelParameters;

            model.Idiosyncratic_Weight = epsilon;

            model.Weights.Clear();
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                var factorWeight = new FactorWeight { Driver = Drivers[i], Weight = beta[i] };
                model.Weights.Add(factorWeight);
            }

            // This method gives normalised weights and factors are real world drivers
            model.PCA_Factors = YesNo.No;

            modelParameters.fDefaultUsed = false;
            return true;
        }

        /// <summary>
        /// Gets the name time series and drivers time series.
        /// </summary>
        private bool TryGetNameAndDriversTimeSeries(ICalibrationData calibrationData, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            // Get name time series.
            if (!GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData, priceFactorToEvolve,
                Data_Retrieval_Parameters, errorList, out fNameTimeSeries))
            {
                return false;
            }

            // Get drivers time series.
            foreach (DriverName driverName in Drivers)
            {
                var driverFactorID = new FactorID(driverName);
                var driverID = new FactorTypeAndID(priceFactorToEvolve.TypeDotSubType(), driverFactorID);
                TimeSeriesStructure<string> timeSeries;
                if (GenericDataRetrieval.TryGetHistoricalTimeSeriesStructure(calibrationData, driverID,
                    Data_Retrieval_Parameters, errorList, out timeSeries))
                {
                    fDriverTimeSeries.Add(timeSeries);
                }
            }

            return true;
        }
    }

    /// <summary>
    /// Base class for regression calibration of a multi factor GBM model on curve (Interest Rate and Hazard Rate).
    /// </summary>
    /// <typeparam name="TModelParameters">The type of model (or model parameters object) we are calibrating.</typeparam>
    public abstract class MultiGBMRateRegressionCalibration<TModelParameters> : RegressionCalibration<TModelParameters>
        where TModelParameters : BaseFactor, IMultiFactorModelParameters
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        protected MultiGBMRateRegressionCalibration()
        {
            Upper_Tenor = 10.0;
            Lower_Tenor = 1.0;
        }

        /// <summary>
        /// Lower tenor for hazard rate.
        /// </summary>  
        public Period Lower_Tenor
        {
            get;
            set;
        }

        /// <summary>
        /// Upper tenor for hazard rate.
        /// </summary>  
        public Period Upper_Tenor
        {
            get;
            set;
        }

        /// <inheritdoc />
        protected override void GetNameDriverStatisticsCorrelations(ICalibrationData sd, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList,
                                                                    out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations)
        {
            IStatistics statistics = sd.Statistics;
            Type factorType = priceFactorToEvolve.GetType();

            // Get indexes of name in statistics set
            List<int> indexes;
            List<double> tenors;
            StatisticsHelper.FindEntries(statistics, factorType.Name, modelParameters.fID, out indexes, out tenors);
            if (indexes.Count == 0)
            {
                throw new CalibrationException(string.Format("No statistics for {0}.", priceFactorToEvolve.GetKey()));
            }

            var nameIndexes = new List<int>();
            for (int k = 0; k < indexes.Count; ++k)
            {
                if (Lower_Tenor <= tenors[k] && tenors[k] <= Upper_Tenor)
                {
                    nameIndexes.Add(indexes[k]);
                }
            }

            if (nameIndexes.Count == 0)
            {
                throw new CalibrationException(string.Format("No statistics in range {0} to {1} for {2}.", Lower_Tenor, Upper_Tenor, priceFactorToEvolve.GetKey()));
            }

            // Calculate normalisation factor for name, = sum of entries in correlation matrix
            double nameNormalise = 0.0;
            for (int k = 1; k < nameIndexes.Count; ++k)
            {
                for (int l = 0; l < k; ++l)
                {
                    nameNormalise += statistics.Correlation(nameIndexes[k], nameIndexes[l]);
                }
            }

            nameNormalise = 2.0 * nameNormalise + nameIndexes.Count;
            nameNormalise = nameNormalise > CalcUtils.TINY ? 1.0 / Math.Sqrt(nameNormalise) : 0.0;

            int numberOfDrivers = Drivers.Count;

            // Get indexes of drivers in statistics set.
            // Define an array of lists of integers.
            List<int>[] driverIndexes = new List<int>[numberOfDrivers];
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                FactorID factorID = new FactorID(Drivers[i]);
                StatisticsHelper.FindEntries(statistics, factorType.Name, factorID, out indexes, out tenors);
                if (indexes.Count == 0)
                {
                    throw new CalibrationException(string.Format("No statistics for driver {0}.", Drivers[i]));
                }

                driverIndexes[i] = new List<int>();
                for (int k = 0; k < indexes.Count; ++k)
                {
                    if (Lower_Tenor <= tenors[k] && tenors[k] <= Upper_Tenor)
                    {
                        driverIndexes[i].Add(indexes[k]);
                    }
                }

                if (driverIndexes[i].Count == 0)
                {
                    throw new CalibrationException(string.Format("No statistics in range {0} to {1} for driver {2}.", Lower_Tenor, Upper_Tenor, Drivers[i]));
                }
            }

            // Calculate normalisation factors for drivers, = sum of entries in correlation matrix
            double[] driverNormalise = new double[numberOfDrivers];
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                driverNormalise[i] = 0.0;
                for (int k = 1; k < driverIndexes[i].Count; ++k)
                {
                    for (int l = 0; l < k; ++l)
                    {
                        driverNormalise[i] += statistics.Correlation(driverIndexes[i][k], driverIndexes[i][l]);
                    }
                }

                driverNormalise[i] = 2.0 * driverNormalise[i] + driverIndexes[i].Count;
                driverNormalise[i] = driverNormalise[i] > CalcUtils.TINY ? 1.0 / Math.Sqrt(driverNormalise[i]) : 0.0;
            }

            // Get correlations of name with drivers
            nameDriverCorrelations = Vector.Create(numberOfDrivers);
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                double correlation = 0.0;
                foreach (int index in nameIndexes)
                {
                    for (int l = 0; l < driverIndexes[i].Count; ++l)
                    {
                        correlation += statistics.Correlation(index, driverIndexes[i][l]);
                    }
                }

                correlation *= nameNormalise * driverNormalise[i];
                nameDriverCorrelations[i] = correlation;
            }

            // Get correlation matrix of drivers
            driverCorrelations = SymmetricMatrix.Create(numberOfDrivers);
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                driverCorrelations[i, i] = 1.0;
                for (int j = 0; j < i; ++j)
                {
                    double correlation = 0.0;
                    for (int k = 0; k < driverIndexes[i].Count; ++k)
                    {
                        for (int l = 0; l < driverIndexes[j].Count; ++l)
                        {
                            correlation += statistics.Correlation(driverIndexes[i][k], driverIndexes[j][l]);
                        }
                    }

                    correlation *= driverNormalise[i] * driverNormalise[j];
                    driverCorrelations[i, j] = correlation;
                }
            }

            SetVolatility(sd, modelParameters, errorList);
        }

        /// <inheritdoc />
        protected override bool GetNameDriverHistoricalCorrelations(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList, out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations)
        {
            nameDriverCorrelations = null;
            driverCorrelations = null;

            // Get correlation of name price factor.
            List<int> validNamePointIndices;
            SymmetricMatrix validNameCorrelationMatrix;
            double nameNormalisingFactor;

            if (!GetValidCorrelations(priceFactorToEvolve.GetID(), errorList, out validNamePointIndices, out validNameCorrelationMatrix, out nameNormalisingFactor))
                return false;

            // Get correlation of driver price factors.
            int numberOfDrivers = Drivers.Count;
            var validDriverPointIndices = new List<int>[numberOfDrivers];
            var validDriverCorrelationMatrices = new SymmetricMatrix[numberOfDrivers];
            var driverNormalisingFactors = new double[numberOfDrivers];
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                var factorID = new FactorID(Drivers[i]);

                if (!GetValidCorrelations(factorID, errorList, out validDriverPointIndices[i], out validDriverCorrelationMatrices[i], out driverNormalisingFactors[i]))
                    return false;
            }

            // Get correlations of name with drivers
            nameDriverCorrelations = Vector.Create(numberOfDrivers);
            for (int driverIndex = 0; driverIndex < numberOfDrivers; ++driverIndex)
            {
                double sumCorrelation = 0.0;
                foreach (int pointIndex in validNamePointIndices)
                {
                    for (int driverPointIndex = 0; driverPointIndex < validDriverPointIndices[driverIndex].Count; ++driverPointIndex)
                    {
                        double correl;
                        CorrelationCalculationHelper.GetCorrelation(
                            fNameTimeSeries.Values[pointIndex], 
                            fNameTimeSeries.ValuesAtHorizon[pointIndex],
                            fDriverTimeSeries[driverIndex].Values[validDriverPointIndices[driverIndex][driverPointIndex]],
                            fDriverTimeSeries[driverIndex].ValuesAtHorizon[validDriverPointIndices[driverIndex][driverPointIndex]], 
                            ReturnMethod.Log, out correl);

                        sumCorrelation += correl;
                    }
                }

                sumCorrelation *= nameNormalisingFactor * driverNormalisingFactors[driverIndex];
                nameDriverCorrelations[driverIndex] = sumCorrelation;
            }

            // Get correlation matrix of drivers.
            driverCorrelations = SymmetricMatrix.Create(numberOfDrivers);
            for (int rowIndex = 0; rowIndex < numberOfDrivers; ++rowIndex)
            {
                driverCorrelations[rowIndex, rowIndex] = 1.0;
                for (int colIndex = 0; colIndex < rowIndex; ++colIndex)
                {
                    double sumCorrelation = 0.0;
                    for (int k = 0; k < validDriverPointIndices[rowIndex].Count; ++k)
                    {
                        for (int l = 0; l < validDriverPointIndices[colIndex].Count; ++l)
                        {
                            double correl;
                            CorrelationCalculationHelper.GetCorrelation(
                                fDriverTimeSeries[rowIndex].Values[validDriverPointIndices[rowIndex][k]],
                                fDriverTimeSeries[rowIndex].ValuesAtHorizon[validDriverPointIndices[rowIndex][k]],
                                fDriverTimeSeries[colIndex].Values[validDriverPointIndices[colIndex][l]],
                                fDriverTimeSeries[colIndex].ValuesAtHorizon[validDriverPointIndices[colIndex][l]],
                                ReturnMethod.Log, out correl);

                            sumCorrelation += correl;
                        }
                    }

                    sumCorrelation *= driverNormalisingFactors[rowIndex] * driverNormalisingFactors[colIndex];
                    driverCorrelations[rowIndex, colIndex] = sumCorrelation;
                }
            }

            SetVolatility(calibrationData, modelParameters, errorList);

            return true;
        }

        /// <summary>
        /// Set the price model volatility.
        /// </summary>
        protected abstract void SetVolatility(ICalibrationData calibrationData, TModelParameters modelParameters, ErrorList errorList);

        /// <summary>
        /// Get a set of valid correlations given a time series.
        /// </summary>
        private bool GetValidCorrelations(FactorID factorID, ErrorList errorList, out List<int> validPointIndices, out SymmetricMatrix validCorrelations, out double normalisingFactor)
        {
            validPointIndices = null;
            validCorrelations = null;
            normalisingFactor = 1.0;

            SymmetricMatrix nameCorrel;
            bool isValid = CorrelationCalculationHelper.GetCorrelationMatrix(fNameTimeSeries, ReturnMethod.Log, out nameCorrel);

            if (!isValid)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error,
                                  "Calibration failed - correlation of name price factor cannot be calculated.");

                return false;
            }

            // Get valid points for name price factor.
            var points = new List<double>(fNameTimeSeries.TimeSeriesIDs.Length);
            foreach (string pointString in fNameTimeSeries.TimeSeriesIDs)
                points.Add(CalcUtils.ParsePoint(pointString));

            validPointIndices = new List<int>();
            for (int index = 0; index < points.Count; ++index)
            {
                if (Lower_Tenor <= points[index] && points[index] <= Upper_Tenor)
                    validPointIndices.Add(index);
            }

            if (validPointIndices.Count == 0)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error,
                                  string.Format("No statistics in range {0} to {1} for {2}.", Lower_Tenor, Upper_Tenor, factorID));

                return false;
            }

            validCorrelations = CorrelationCalculationHelper.GetSelectedMatrix(validPointIndices, nameCorrel);

            // Calculate normalisation factor for name, = sum of entries in correlation matrix
            normalisingFactor = 0.0;
            for (int rowIndex = 1; rowIndex < validCorrelations.Rows; ++rowIndex)
            {
                for (int colIndex = 0; colIndex < rowIndex; ++colIndex)
                {
                    normalisingFactor += validCorrelations[rowIndex, colIndex];
                }
            }

            normalisingFactor = 2.0 * normalisingFactor + validCorrelations.Rows;
            normalisingFactor = normalisingFactor > CalcUtils.TINY ? 1.0 / Math.Sqrt(normalisingFactor) : 0.0;

            return true;
        }
    }

    /// <summary>
    /// Class for regression calibration of a multi factor GBM model on Hazard Rate.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get the pre-computed statistics. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file.
    /// </remarks>
    [DisplayName("Multi-Factor GBM Hazard Rate Regression Calibration")]
    public class MultiGBMHazardRateRegressionCalibration : MultiGBMRateRegressionCalibration<MultiGBMHazardRateModel>, IModelCalibration
    {
        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMHazardRateModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var modelParameters = CalibrationHelper.Validate<MultiGBMHazardRateModel>(calibrationData, priceModel, output);
            if (!Calibrate(calibrationData, modelParameters, CalibrationHelper.GetPriceFactor(priceModel), errorList))
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        protected override void SetVolatility(ICalibrationData calibrationData, MultiGBMHazardRateModel modelParameters, ErrorList errorList)
        {
            // Moment statistics with mean with no volatility adjustment, 'BasePlugInStatisticsCalculations.VolatilityString'.
            MomentStatistics[] momentStatistics;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, modelParameters, ReturnMethod.Log, errorList, ref fNameTimeSeries, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            // Get average volatility across tenor points.
            int length = momentStatistics.Length;
            double volatility = 0.0;
            int validLength = 0;
            for (int index = 0; index < length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                volatility += momentStatistics[index].Vol;
                validLength++;
            }

            if (validLength > 0)
                modelParameters.Vol = volatility / validLength;
            else
                errorList.Add(ErrorLevel.Error, string.Format("Calibration of {0} failed.", modelParameters));
        }
    }

    /// <summary>
    /// Class for regression calibration of a multi-factor GOU model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get the pre-computed statistics. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file.
    /// </remarks>
    /// <typeparam name="TModel">The type of GOU model that we are calibrating.</typeparam>
    public class MultiGOURegressionCalibration<TModel> : MultiGBMRateRegressionCalibration<TModel>, IModelCalibration
        where TModel : PriceModel, IMultiGOUModel, IMultiFactorModelParameters
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public MultiGOURegressionCalibration()
        {
            Idiosyncratic_Factor = YesNo.No;
        }

        /// <summary>
        /// Flag indicating whether we wish to include an idiosyncratic factor.
        /// </summary>
        public YesNo Idiosyncratic_Factor
        {
            get;
            set;
        }

        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(TModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<TModel>(calibrationData, priceModel, output);
            if (!Calibrate(calibrationData, model, CalibrationHelper.GetPriceFactor(priceModel), errorList))
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        protected override void SetVolatility(ICalibrationData calibrationData, TModel model, ErrorList errorList)
        {
            // Get mean reversion statistics.
            MeanReversionStatistics[] meanReversionStatistics;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Log, null, null, null, null, false);
            bool isSuccess = GetMeanReversionStatistics(calibrationData, model, calculationParameters, errorList,
                                                        ref fNameTimeSeries, out meanReversionStatistics, out isValidStatistics);

            if (!isSuccess)
                return;

            MultiGOUModel underlyingModel = model.GetUnderlyingModel();
            underlyingModel.Vol.Clear();

            int length = meanReversionStatistics.Length;
            double alphaSum = 0.0;
            for (int index = 0; index < length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                alphaSum += meanReversionStatistics[index].Alpha;

                // Build volatility curve from all available tenors between fLowerTenor and fUpperTenor
                underlyingModel.Vol.X.Add(meanReversionStatistics[index].ParsePointIDToDouble());
                underlyingModel.Vol.Y.Add(meanReversionStatistics[index].Sigma);
            }

            // Mean reversion speed is set to the average of all available tenors between fLowerTenor and fUpperTenor
            underlyingModel.Mean_Reversion_Speed = alphaSum / length;
        }

        /// <inheritdoc />
        protected override YesNo GetIdiosyncraticFactorChoice()
        {
            return Idiosyncratic_Factor;
        }
    }

    /// <summary>
    /// Class for regression calibration of a multi-factor GOU model on Interest Rate.
    /// </summary>
    [DisplayName("Multi-Factor GOU Interest Rate Regression Calibration")]
    public class MultiGOUInterestRateRegressionCalibration : MultiGOURegressionCalibration<MultiGOUInterestRateModel>
    {
    }

    /// <summary>
    /// Class for regression calibration of a multi-factor GOU model on Hazard Rate.
    /// </summary>
    [DisplayName("Multi-Factor GOU Hazard Rate Regression Calibration")]
    public class MultiGOUHazardRateRegressionCalibration : MultiGOURegressionCalibration<MultiGOUHazardRateModel>
    {
    }

    /// <summary>
    /// Class for regression calibration of a multi factor GBM asset price model.
    /// </summary>
    /// <typeparam name="TModelParameters">The type of model (or model parameters object) we are calibrating.</typeparam>
    public abstract class MultiGBMAssetPriceRegressionCalibrationBase<TModelParameters> : RegressionCalibration<TModelParameters>
        where TModelParameters : BaseFactor, IMultiFactorModelParameters
    {
        /// <inheritdoc />
        protected override void GetNameDriverStatisticsCorrelations(ICalibrationData sd, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList,
                                                                    out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations)
        {
            IStatistics statistics = sd.Statistics;
            Type factorType = priceFactorToEvolve.GetType();

            // Get index of name in statistics set
            int nameIndex = statistics.GetIndex(factorType.Name, modelParameters.fID, null);
            if (nameIndex < 0)
            {
                throw new CalibrationException(string.Format("No statistics for {0}.", priceFactorToEvolve.GetKey()));
            }

            int numberOfDrivers = Drivers.Count;

            // Get indexes of drivers in statistics set.
            List<int> driverIndexes = new List<int>();
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                FactorID factorID = new FactorID(Drivers[i]);
                int index = statistics.GetIndex(factorType.Name, factorID, null);
                if (index < 0)
                {
                    throw new CalibrationException(string.Format("No statistics for driver {0}.", factorID));
                }

                driverIndexes.Add(index);
            }

            // Set the drift and volatility model
            SetModelDriftAndVol(sd, modelParameters, priceFactorToEvolve, errorList);

            // Get correlations of name with drivers
            nameDriverCorrelations = Vector.Create(numberOfDrivers);
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                nameDriverCorrelations[i] = statistics.Correlation(nameIndex, driverIndexes[i]);
            }

            // Get correlation matrix of drivers
            driverCorrelations = SymmetricMatrix.Create(numberOfDrivers);
            for (int i = 0; i < numberOfDrivers; ++i)
            {
                driverCorrelations[i, i] = 1.0;
                for (int j = 0; j < i; ++j)
                {
                    double correlation = statistics.Correlation(driverIndexes[i], driverIndexes[j]);
                    driverCorrelations[i, j] = correlation;
                }
            }
        }

        /// <inheritdoc />    
        protected override bool GetNameDriverHistoricalCorrelations(ICalibrationData calibrationData, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList, out Vector nameDriverCorrelations, out SymmetricMatrix driverCorrelations)
        {
            bool returnValue = true;

            // Set the drift and volatility of the model.
            SetModelDriftAndVol(calibrationData, modelParameters, priceFactorToEvolve, errorList);

            int numberOfDrivers = Drivers.Count;

            // Get correlations of name with drivers.
            nameDriverCorrelations = Vector.Create(numberOfDrivers);
            for (int driverIndex = 0; driverIndex < numberOfDrivers; ++driverIndex)
            {
                double correl;
                if (!CorrelationCalculationHelper.GetCorrelation(
                    fNameTimeSeries.Values[0],
                    fNameTimeSeries.ValuesAtHorizon[0],
                    fDriverTimeSeries[driverIndex].Values[0],
                    fDriverTimeSeries[driverIndex].ValuesAtHorizon[0],
                    ReturnMethod.Log, out correl))
                {
                    errorList.Add(ErrorLevel.Error, string.Format("Could not calculate correlation between {0} and driver {1}.", priceFactorToEvolve.GetID(), driverIndex));
                    returnValue = false;
                }
                    
                nameDriverCorrelations[driverIndex] = correl;
            }

            // Get correlation matrix of drivers.
            driverCorrelations = SymmetricMatrix.Create(numberOfDrivers);
            for (int rowIndex = 0; rowIndex < numberOfDrivers; ++rowIndex)
            {
                driverCorrelations[rowIndex, rowIndex] = 1.0;
                for (int colIndex = 0; colIndex < rowIndex; ++colIndex)
                {
                    double correlation;
                    if (!CorrelationCalculationHelper.GetCorrelation(
                        fDriverTimeSeries[rowIndex].Values[0],
                        fDriverTimeSeries[rowIndex].ValuesAtHorizon[0],
                        fDriverTimeSeries[colIndex].Values[0],
                        fDriverTimeSeries[colIndex].ValuesAtHorizon[0],
                        ReturnMethod.Log, out correlation))
                    {
                        errorList.Add(ErrorLevel.Error, string.Format("Could not calculate correlation between drivers {0} and {1}.", rowIndex, colIndex));
                        returnValue = false;
                    }

                    driverCorrelations[rowIndex, colIndex] = correlation;
                }
            }

            return returnValue;
        }

        /// <summary>
        /// Set the drift and vol of the model.
        /// </summary>
        /// <param name="sd">Market Data and Statistics set.</param>
        /// <param name="modelParameters">Model parameters you want to set.</param>
        /// <param name="priceFactorToEvolve">The pricefactor that is being evolved.</param>
        /// <param name="errorList">Error list of error message.</param>
        protected abstract void SetModelDriftAndVol(ICalibrationData sd, TModelParameters modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList);
    }

    /// <summary>
    /// Class for regression calibration of a multi factor GBM asset price model.
    /// </summary>
    [DisplayName("Multi-Factor GBM Asset Price Regression Calibration")]
    public class MultiGBMAssetPriceRegressionCalibration : MultiGBMAssetPriceRegressionCalibrationBase<MultiGBMAssetPriceModel>, IModelCalibration
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
            MultiGBMAssetPriceModel model = CalibrationHelper.Validate<MultiGBMAssetPriceModel>(calibrationData, priceModel, output);
            if (!Calibrate(calibrationData, model, CalibrationHelper.GetPriceFactor(priceModel), errorList))
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        protected override void SetModelDriftAndVol(ICalibrationData calibrationData, MultiGBMAssetPriceModel modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            var priceFactor = CalibrationHelper.GetPriceFactor(modelParameters);

            // Moment statistics using name time series.
            MomentStatistics[] momentStatistics;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, modelParameters, ReturnMethod.Log, errorList, ref fNameTimeSeries, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics at first point for {0}.", priceFactor.GetKey()));
                return;
            }

            // Set model parameters using earliest point.
            modelParameters.Drift = momentStatistics[0].GetGbmDrift();
            modelParameters.Vol = momentStatistics[0].Vol;
        }
    }

    /// <summary>
    /// Multi-Factor GBM Asset Price Beta Calibration with the benchmark read from the grouping file.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get the pre-computed statisticse. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file.
    /// </remarks>
    [DisplayName("Multi-Factor GBM Asset Price Beta Calibration")]
    public class MultiGBMAssetPriceBetaCalibration : StatisticsAndBetaBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public MultiGBMAssetPriceBetaCalibration()
        {
            Idiosyncratic_Factor = YesNo.No;
            Significance_Threshold = CalcUtils.DefaultMinProcessWeight;
        }

        /// <summary>
        /// Flag indicating whether we wish to include an idiosycratic factor.
        /// </summary>
        public YesNo Idiosyncratic_Factor
        {
            get;
            set;
        }

        /// <summary>
        /// Value below which we consider a weight to be zero.
        /// </summary>
        public double Significance_Threshold
        {
            get;
            set;
        }

        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceModel);
        }

        /// <summary>
        /// Calibrates the price model from the statistics files 
        /// using linear regression.
        /// </summary>
        /// <remarks>
        /// Due to a request from Services we keep this method virtual (in case they wish to modify the calibration).
        /// </remarks>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var modelParameters = CalibrationHelper.Validate<MultiGBMAssetPriceModel>(calibrationData, priceModel, output);

            bool isValid = (Use_Pre_Computed_Statistics == YesNo.Yes)
                               ? CalculateFromPreComputedStatistics(calibrationData, modelParameters, errorList)
                               : CalculateFromHistoricalStatistics(calibrationData, modelParameters, errorList);

            if (!isValid)
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// Calculate the model parameters using pre-computed statistics.
        /// </summary>
        /// <param name="sd">Statistics set.</param>
        /// <param name="modelParameters">Multi-factor model to calibrate.</param>
        /// <param name="errorList">Errorlist to output errors.</param>
        private bool CalculateFromPreComputedStatistics(ICalibrationData sd, MultiGBMAssetPriceModel modelParameters, ErrorList errorList)
        {
            Type factorType = CalibrationHelper.GetPriceFactor(modelParameters).GetType();

            IStatistics statistics = sd.Statistics;

            int nameIndex = statistics.GetIndex(factorType.Name, modelParameters.fID, null);
            if (nameIndex < 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No statistics for {0} {1}.", factorType.Name, modelParameters.fID));
                return false;
            }

            string driverString = statistics.GetEquityIndex(nameIndex);
            if (String.IsNullOrEmpty(driverString))
            {
                errorList.Add(ErrorLevel.Error, string.Format("No index for {0} {1}.", factorType.Name, modelParameters.fID));
                return false;
            }

            string codeSeparator = FactorID.CodeSeparator.ToString();
            if (!driverString.Contains(codeSeparator))
                driverString = factorType.Name + codeSeparator + driverString;

            string driverType, driverPoint;
            FactorID driverID;
            StatisticsData.ParseScalarID(driverString, out driverType, out driverID, out driverPoint);

            if (driverType != factorType.Name || !string.IsNullOrEmpty(driverPoint))
            {
                errorList.Add(ErrorLevel.Error, string.Format("Invalid index {0} for {1} {2}.", driverString, factorType.Name, modelParameters.fID));
                return false;
            }

            int driverIndex = statistics.GetIndex(driverType, driverID, driverPoint);
            if (driverIndex < 0)
            {
                errorList.Add(ErrorLevel.Error, string.Format("Statistics for the index {0} {1} is missing, calibration cannot assign the model volatility.", factorType.Name, driverID));
                return false;
            }

            // retrieve beta's name
            double beta = StatisticsHelper.GetStatistic(statistics, nameIndex, StatisticsData.StatisticBeta, ReturnMethod.Log); // beta always log method

            double volDriver;
            double driftDriver;
            bool hasDriverDrift = StatisticsHelper.TryGetGbmDrift(statistics, driverIndex, ReturnMethod.Log, out driftDriver);
            bool hasDriverVol = StatisticsHelper.TryGetGbmVol(statistics, driverIndex, ReturnMethod.Log, out volDriver);

            // retrieve the drift and volatility of response
            double volName;
            double drift;
            bool hasNameDrift = StatisticsHelper.TryGetGbmDrift(statistics, nameIndex, ReturnMethod.Log, out drift);
            bool hasNameVol = StatisticsHelper.TryGetGbmVol(statistics, nameIndex, ReturnMethod.Log, out volName);

            // TRADED asset without a vol
            if (!hasDriverVol)
            {
                throw new CalibrationException(string.Format("No vol for {0} {1} is supplied.", factorType.Name, driverID));
            }

            // Calibration fails if Idiosyncratic Factor = Yes and vol of asset is not provided.
            // Idiosyncratic Factor = No, vol of asset is not required.
            if (!hasNameVol && Idiosyncratic_Factor == YesNo.Yes)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No vol for {0} {1} is supplied.", factorType.Name, modelParameters.fID));
                return false;
            }

            if (!hasNameDrift)
            {
                // NON TRADED asset doesn't have a drift - use the driver's drift instead
                if (!hasDriverDrift)
                {
                    // driver doesn't have a drift
                    errorList.Add(ErrorLevel.Error, string.Format("No drift for name {0} {1} nor index {0} {2} is supplied.", factorType.Name, modelParameters.fID, driverID));
                    return false;
                }

                drift = driftDriver;
            }

            Vector systemicWeight = Vector.Create(1);
            double idiosyncraticWeight = 0.0;

            if (hasNameVol && Idiosyncratic_Factor == YesNo.Yes)
            {
                // traded asset price with idiosyncratic factor required
                modelParameters.Vol = volName;
                double correlation = beta * volDriver / volName;

                if (correlation > 1.0 || correlation < -1.0)
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                        "Correlation between {0} {1} and {0} {2} retrieved from beta and volatilities is not valid.",
                        factorType.Name, modelParameters.fID, driverID));
                    return false;
                }

                systemicWeight[0] = correlation;
                MultiFactorCalibrationHelper.ObtainRoundedNormalisedRegressionWeights(Significance_Threshold, 1.0, Idiosyncratic_Factor, null, systemicWeight, out idiosyncraticWeight);
            }
            else
            {
                // idiosyncratic flag doesn't apply to non traded asset price
                modelParameters.Vol = Math.Abs(beta) * volDriver;
                systemicWeight[0] = Math.Sign(beta);
            }

            // set model parameters
            modelParameters.Weights.Clear();
            var factorWeight = new FactorWeight
                                   {
                                           Driver = driverID.ToCode(),
                                           Weight = systemicWeight[0]
                                   };

            // Driver name should be based on unchanged ID (in case the ID contains a dot).
            modelParameters.Idiosyncratic_Weight = idiosyncraticWeight;
            modelParameters.Weights.Add(factorWeight);
            modelParameters.Drift = drift;

            // driving factor are real world driver and output weights are always normalised
            modelParameters.PCA_Factors = YesNo.No;

            modelParameters.fDefaultUsed = false;

            return true;
        }

        /// <summary>
        /// Calibration using historical market data.
        /// </summary>
        private bool CalculateFromHistoricalStatistics(ICalibrationData calibrationData, MultiGBMAssetPriceModel modelParameters, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(modelParameters);

            // Moment statistics with volatility adjustment for name price factor.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> nameTimeSeriesStructure = null;
            bool[] isNameValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, modelParameters, ReturnMethod.Log, errorList, ref nameTimeSeriesStructure, out momentStatistics, out isNameValidStatistics))
                return false;

            if (!isNameValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for name price factor at first point for {0}.", priceFactor.GetKey()));
                return false;
            }

            // Selecting first point.
            double volName = momentStatistics[0].Vol;
            double driftName = momentStatistics[0].GetGbmDrift();

            // Get groupings from file.
            List<GroupingFileMatch> groupFileMatches = GetDriverID(calibrationData, priceFactor, errorList);

            if (groupFileMatches == null || groupFileMatches.Count == 0)
            {
                errorList.Add(ErrorLevel.Error, "Driver not found in the grouping file.");
                return false;
            }

            // Selecting first point - as done in the first method.
            GroupingFileMatch match = groupFileMatches[0];
            var driverFactorTypeAndID = new FactorTypeAndID(priceFactor.TypeDotSubType(), match.fDriver.fFactorID);

            // Moment statistics with volatility adjustment for driver price factor.
            TimeSeriesStructure<string> driverTimeSeriesStructure = null;
            bool[] isDriverValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, driverFactorTypeAndID, ReturnMethod.Log, errorList,
                                               ref driverTimeSeriesStructure, out momentStatistics, out isDriverValidStatistics))
            {
                return false;
            }

            if (!isDriverValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics for drivers at first point for {0}.", priceFactor.GetKey()));
                return false;
            }

            // Selecting first point.
            double volDriver = momentStatistics[0].Vol;
            double driftDriver = momentStatistics[0].GetGbmDrift();

            // Beta always log method. Calculation based on log.
            double beta = GetBeta(nameTimeSeriesStructure.Values[0], nameTimeSeriesStructure.ValuesAtHorizon[0], driverTimeSeriesStructure.Values[0], driverTimeSeriesStructure.ValuesAtHorizon[0], volName, volDriver, ReturnMethod.Log);

            // TRADED asset without a vol
            if (!isDriverValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error,
                              string.Format("No vol for {0} {1} is supplied.", priceFactor.TypeDotSubType(),
                                            driverFactorTypeAndID.fFactorID));
                return false;
            }

            if (!isNameValidStatistics[0])
            {
                // NON TRADED asset doesn't have a drift - use the driver's drift instead
                if (!isDriverValidStatistics[0])
                {
                    // driver doesn't have a drift
                    errorList.Add(ErrorLevel.Error,
                                  string.Format("No drift for name {0} {1} nor index {0} {2} is supplied.",
                                                priceFactor.TypeDotSubType(), modelParameters.fID, driverFactorTypeAndID.fFactorID));
                    return false;
                }

                driftName = driftDriver;
            }

            Vector systemicWeight = Vector.Create(1);
            double idiosyncraticWeight = 0.0;

            if (isNameValidStatistics[0] && Idiosyncratic_Factor == YesNo.Yes)
            {
                // traded asset price with idiosyncratic factor required
                modelParameters.Vol = volName;
                double correlation = beta * volDriver / volName;

                if (correlation > 1.0 || correlation < -1.0)
                {
                    errorList.Add(ErrorLevel.Error, string.Format(
                        "Correlation between {0} {1} and {0} {2} retrieved from beta and volatilities is not valid.",
                        priceFactor.TypeDotSubType(), modelParameters.fID, driverFactorTypeAndID.fFactorID));

                    return false;
                }

                systemicWeight[0] = correlation;
                MultiFactorCalibrationHelper.ObtainRoundedNormalisedRegressionWeights(Significance_Threshold, 1.0, Idiosyncratic_Factor, null, systemicWeight, out idiosyncraticWeight);
            }
            else
            {
                // idiosyncratic flag doesn't apply to non traded asset price
                modelParameters.Vol = Math.Abs(beta) * volDriver;
                systemicWeight[0] = Math.Sign(beta);
            }

            // set model parameters
            modelParameters.Weights.Clear();
            var factorWeight = new FactorWeight
            {
                // Driver name should be based on unchanged ID (in case the ID contains a dot).
                Driver = driverFactorTypeAndID.fFactorID.ToCode(),
                Weight = systemicWeight[0]
            };

            modelParameters.Idiosyncratic_Weight = idiosyncraticWeight;
            modelParameters.Weights.Add(factorWeight);
            modelParameters.Drift = driftName;

            // driving factor are real world driver and output weights are always normalised
            modelParameters.PCA_Factors = YesNo.No;
            modelParameters.fDefaultUsed = false;

            return true;
        }
    }

    /// <summary>
    /// Regression Calibration for the Multi GBM term-structure "Implied" model.
    /// </summary>
    /// <remarks>
    /// We do not set the vol here, as the model obtains the vol from a separate <see cref="GBMTSImpliedParameters"/> pricefactor.
    /// </remarks>
    [DisplayName("Multi-Factor GBM TS Implied Model Calibration")]
    public class MultiGBMAssetPriceTSRegressionCalibrationImplied : MultiGBMAssetPriceRegressionCalibrationBase<MultiGBMAssetPriceTSModelImplied>, IModelCalibration
    {
        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceTSModelImplied);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<MultiGBMAssetPriceTSModelImplied>(calibrationData, priceModel, output);
            if (!Calibrate(calibrationData, model, CalibrationHelper.GetPriceFactor(priceModel), errorList))
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <inheritdoc />
        protected override void SetModelDriftAndVol(ICalibrationData calibrationData, MultiGBMAssetPriceTSModelImplied modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            // We do nothing.
            // The vol calibration is implied, so is handled by the bootstrapper.
            // We take Risk_Neutral_Drift = Yes for the implied model, so no drift required.
        }
    }

    /// <summary>
    /// Regression Calibration for the Multi GBM term-structure model.
    /// </summary>
    /// <remarks>
    /// We do not set the vol here, as the model obtains the vol from a separate <see cref="GBMTSImpliedParameters"/> pricefactor.
    /// </remarks>
    [DisplayName("Multi-Factor GBM TS Model Calibration")]
    public class MultiGBMAssetPriceTSRegressionCalibration : MultiGBMAssetPriceRegressionCalibrationBase<MultiGBMAssetPriceTSModel>, IModelCalibration
    {
        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public MultiGBMAssetPriceTSRegressionCalibration()
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
        /// If set to No, volatility will be set to a constant volatility, calculated from historical vols.
        /// </summary>
        public VolCalibration Vol_Calibration { get; set; }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(MultiGBMAssetPriceTSModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<MultiGBMAssetPriceTSModel>(calibrationData, priceModel, output);
            if (!Calibrate(calibrationData, model, CalibrationHelper.GetPriceFactor(priceModel), errorList))
                return;

            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var multiGBMTSModel = (MultiGBMAssetPriceTSModel)priceModel;
            GBMAssetPriceTSCalibrationHelper.CalibrationOutput(multiGBMTSModel.Drift, multiGBMTSModel.Vol, output);
        }

        /// <inheritdoc />
        protected override void SetModelDriftAndVol(ICalibrationData calibrationData, MultiGBMAssetPriceTSModel modelParameters, IPriceFactor priceFactorToEvolve, ErrorList errorList)
        {
            // Set the risk neutral flag on the model if requested
            if (Set_Risk_Neutral_Drift == YesNo.Yes)
                modelParameters.Risk_Neutral_Drift = Risk_Neutral_Drift;

            // Moment statistics with volatility adjustment.
            MomentStatistics[] momentStatistics;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, modelParameters, ReturnMethod.Log, errorList, ref fNameTimeSeries, out momentStatistics, out isValidStatistics))
            {
                return;
            }

            Debug.Assert(momentStatistics.Length == 1, "Asset prices are one dimensional.");

            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics at first point for {0}.", priceFactorToEvolve.GetKey()));
                return;
            }

            if (Vol_Calibration == VolCalibration.Implied)
            {
                GBMAssetPriceTSCalibrationHelper.VolCalibrationFromPriceFactor((IAssetPriceWithVolatilityPriceFactor)priceFactorToEvolve, calibrationData.PriceFactors,
                                                                                modelParameters.Risk_Neutral_Drift, modelParameters, null);
            }
            else
            {
                modelParameters.Vol.AssignConstant(momentStatistics[0].Vol);
            }

            // Set model parameters using earliest point.
            GBMAssetPriceTSCalibrationHelper.CreateDriftWithVolatilityCorrection(modelParameters, momentStatistics[0].Drift, modelParameters.GetIntegratedVol());
        }
    }
}
