using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Base class for calibration of PCA interest rate models.
    /// </summary>
    [DisplayName("PCA Interest Rate Calibration")]
    public class PCAInterestRateStatisticsCalibration
    {
        private int fNumPcaFactors;
        private PCSource fMatrixType;
        private DriftMethod fDriftMethod;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public PCAInterestRateStatisticsCalibration()
        {
            fNumPcaFactors = 3;
            fMatrixType = PCSource.Correlation;
            fDriftMethod = DriftMethod.Log_Drift_To_Blend;
            Distribution_Type = ProbabilityDistribution.Lognormal;
        }

        /// <summary>
        /// The number of PCA factors.
        /// </summary>
        public int Number_Of_PCA_Factors
        {
            get { return fNumPcaFactors; }
            set { fNumPcaFactors = Math.Max(value, 0); }
        }

        /// <summary>
        /// User-defined property for whether to use covariance or correlation to calculate beta.
        /// </summary>
        public PCSource Matrix_Type
        {
            get { return fMatrixType; }
            set { fMatrixType = value; }
        }

        /// <summary>
        /// User-defined property to set the drift of the model.
        /// </summary>
        public DriftMethod Rate_Drift_Model
        {
            get { return fDriftMethod; }
            set { fDriftMethod = value; }
        }

        /// <summary>
        /// Model distribution type, i.e., normal or lognormal.
        /// </summary>
        public ProbabilityDistribution Distribution_Type { get; set; }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            PCAInterestRateModel pcaModel = CalibrationHelper.Validate<PCAInterestRateModel>(calibrationData, priceModel, output);

            pcaModel.Princ_Comp_Source = fMatrixType;
            pcaModel.Rate_Drift_Model = fDriftMethod;
            pcaModel.Distribution_Type = Distribution_Type;

            CalibratePcaInterestRateModel(calibrationData.Statistics, pcaModel, errorList);

            GetUpdatedCurves(pcaModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        private void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            WriteCapVolCalibrationDataAsXML((PCAInterestRateModel)priceModel, output);
        }

        /// <summary>
        /// Writes calibration data to an XML stream.
        /// </summary>
        private static void WriteCapVolCalibrationDataAsXML(PCAInterestRateModel pcaModel, XmlWriter writer)
        {
            var graphData = new ResultSeriesList
            {
                fXAxis = { fLabel = "Time (Years)" },
                fYAxis = { fLabel = "Value" }
            };

            graphData.Add(new ResultSeries("Historical Yield", LineStyle.Solid, pcaModel.Historical_Yield.Clone()));
            graphData.Add(new ResultSeries("Yield Volatility", LineStyle.Solid, pcaModel.Yield_Volatility.Clone()));

            for (int k = 0; k < pcaModel.GetDimension(); ++k)
            {
                graphData.Add(new ResultSeries("Eigenvector " + (k + 1).ToString(CultureInfo.InvariantCulture), LineStyle.Solid,
                    pcaModel.Eigenvectors[k].Eigenvector.Clone()));
            }

            graphData.ToXml(writer);
        }

        /// <summary>
        /// Calculate the PCA model parameters from the statistics of log changes in zero rates.
        /// </summary>
        private void CalibratePcaInterestRateModel(IStatistics statistics, PCAInterestRateModel pcaModel, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(pcaModel);

            List<int> indexes;
            List<string> points;
            List<double> reversionVols, reversionSpeeds, longRunMeans;

            ReturnMethod returnMethod = pcaModel.Distribution_Type == ProbabilityDistribution.Lognormal ? ReturnMethod.Log : ReturnMethod.Diff;
            StatisticsHelper.GetMeanReversionStatistics(statistics, priceFactor.GetType(), pcaModel.fID, returnMethod, out indexes, out points, out reversionVols, out reversionSpeeds, out longRunMeans);

            // Number of points with non-zero volatility
            int nRates = indexes.Count;
            if (nRates < fNumPcaFactors)
            {
                errorList.Add(ErrorLevel.Error,
                              string.Format(
                                  "Cannot calibrate {0}-factor {1} because only {2} curve points in statistics ({3}.{4}) with non-zero volatility.",
                                  fNumPcaFactors, pcaModel.GetKey(), nRates, priceFactor.GetType().Name, pcaModel.fID));
                return;
            }

            errorList.Add(ErrorLevel.Info, string.Format(
                              "Calibration using pre-computed {0} Mean-Reversion Statistics for {1} successful.", returnMethod, pcaModel.GetKey()));

            // Convert point strings to zero rate tenors,
            // set volatility, calculate mean reversion speed and set long run mean curve.
            var terms = new double[nRates];

            pcaModel.Yield_Volatility.Clear();
            pcaModel.Historical_Yield.Clear();

            double sumReversionSpeed = 0.0;
            for (int i = 0; i < nRates; ++i)
            {
                terms[i] = CalcUtils.ParsePoint(points[i]);

                pcaModel.Historical_Yield[terms[i]] = longRunMeans[i];
                pcaModel.Yield_Volatility[terms[i]] = reversionVols[i];
                sumReversionSpeed += reversionSpeeds[i];
            }

            pcaModel.Reversion_Speed = sumReversionSpeed / nRates;

            // Construct the covariance or correlation matrix
            var correlationOrCovariance = SymmetricMatrix.Create(nRates);
            for (int i = 0; i < nRates; ++i)
            {
                for (int j = 0; j <= i; ++j)
                {
                    double c = j < i ? statistics.Correlation(indexes[i], indexes[j]) : 1.0;
                    if (fMatrixType == PCSource.Covariance)
                        c *= reversionVols[i] * reversionVols[j];

                    correlationOrCovariance[i, j] = c;
                }
            }

            // Eigenvector decomposition
            var eigenCalculationResults = correlationOrCovariance.CalculateEigenVectorsAndValues();

            // Set eigenvector/eigenvalue model parameters
            pcaModel.Eigenvectors.Clear();
            for (int k = 0; k < fNumPcaFactors; ++k)
            {
                pcaModel.Eigenvectors.Add(new Eigen<RateCurve>());
                pcaModel.Eigenvectors[k].Eigenvalue = eigenCalculationResults.EigenValues[k];

                for (int i = 0; i < nRates; ++i)
                {
                    pcaModel.Eigenvectors[k].Eigenvector[terms[i]] = eigenCalculationResults.EigenVectors[i, k];
                }
            }
        }
    }
}