using SunGard.Adaptiv.Analytics.Framework;
using System;
using System.Xml;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration for Clewlow-Strickland Multi-Factor Forward Price model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    public class CSMultiFactorForwardPriceCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Use a drift if Yes.
        /// </summary>
        public YesNo Use_Drift { get; set; }

        /// <summary>
        /// Calibrate multi-factored Clewlow-Strickland forward price model from historical forward prices
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            CSMultiFactorForwardPriceModel model = CalibrationHelper.Validate<CSMultiFactorForwardPriceModel>(calibrationData, priceModel, output);
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(priceModel);

            // Obtain covariance matrix and drift from file
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidMomentStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, priceModel, ReturnMethod.Log, errorList, ref timeSeriesStructure, out momentStatistics, out isValidMomentStatistics))
            {
                return;
            }

            int nPoints = momentStatistics.Length;

            var xValues = new DoubleArray();  // all fwd price curves have subset of same X values (e.g., X = [6m, 1y,...] = [0.5, 1,...])
            var driftsList = new DriftsList();
            var driftCurve = new RateCurve();
            var stdev = new double[nPoints];

            for (int index = 0; index < nPoints; index++)
            {
                if (!isValidMomentStatistics[index])
                    continue;

                xValues.Add(momentStatistics[index].ParsePointIDToDouble());

                stdev[index] = momentStatistics[index].Vol;
                driftCurve[xValues[index]] = Use_Drift == YesNo.Yes ? momentStatistics[index].GetGbmDrift() : 0.0;
            }

            var instantaneousDriftCurve = new InstantaneousDriftCurve { Drift = driftCurve };
            driftsList[0] = instantaneousDriftCurve;

            // Calculate covariance matrix from historical data or read the pre-computed matrix.
            SymmetricMatrix covarianceMatrix;
            bool isValidCovarianceMatrix = Use_Pre_Computed_Statistics == YesNo.Yes
                                           ? GetCovarianceMatrix(calibrationData, model, stdev, out covarianceMatrix)
                                           : CorrelationCalculationHelper.GetCovarianceMatrix(timeSeriesStructure, ReturnMethod.Log, out covarianceMatrix);

            if (!isValidCovarianceMatrix)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No covariance calculated for {0}.", priceFactor.GetKey()));
                return;
            }

            var eigenResult = covarianceMatrix.CalculateEigenVectorsAndValues();

            // Create new volatilities curves and populate the model
            int numFactors = model.NumOfFactors();
            var instantaneousVols = new VolatilitiesList(numFactors);

            for (int j = 0; j < numFactors; j++)
            {
                var volCurve = new RateCurve();
                if (j < nPoints)
                {
                    // Eigenvalues of covariance matrix should never be zero, but in case of estimation/numerical errors we floor it at 0.0.
                    for (int i = 0; i < nPoints; i++)
                        volCurve[xValues[i]] = Math.Sqrt(Math.Max(eigenResult.EigenValues[j], 0.0)) * eigenResult.EigenVectors[i, j];
                }
                else
                {
                    // If there are more PCA factors than historical variables, we assign the remaining PCA vols to zero. 
                    // In this case, entire volatility is captured, and we simply analyse the same covariance matrix in an orthogonal space.
                    foreach (float x in xValues)
                        volCurve[x] = 0.0;
                }

                // For display only; typically the first eigenvector all has the same sign, signifying a parallel shift. 
                // If unfortunately this sign comes out to be negative from PCA, we flip it to positive.
                if (j == 0 && volCurve[xValues[0]] < 0.0)
                {
                    foreach (float x in xValues)
                        volCurve[x] *= -1.0;
                }

                var instantaneousVolCurve = new InstantaneousVolCurve { Volatility = volCurve };
                instantaneousVols[j] = instantaneousVolCurve;
            }

            model.InstantaneousVols = instantaneousVols;
            model.InstantaneousDrift = driftsList;

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Get the updated curved for display in calibration
        /// </summary>
        /// <param name="priceModel"> Multi-factored CS forward price model</param>
        /// <param name="output">xml writer</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            CSMultiFactorForwardPriceModel model = (CSMultiFactorForwardPriceModel)priceModel;

            int numOfFactors = model.NumOfFactors();
            var graphData = new ResultSeriesList();

            graphData.fXAxis.fLabel = "Time to Maturity (Years)";
            graphData.fYAxis.fLabel = Use_Drift == YesNo.Yes ? "Volatility PCA Components and Drift (%)" : "Volatility PCA Components (%)";

            for (int i = 0; i < model.NumOfFactors(); i++)
            {
                var volCurve = model.InstantaneousVols[i].Volatility.Clone();
                volCurve.MultiplyBy(Percentage.OverPercentagePoint);
                graphData.Add(new ResultSeries(string.Format("Volatility PCA Component {0}", i), LineStyle.Solid, volCurve));
            }

            if (Use_Drift == YesNo.Yes)
            {
                var driftCurve = model.InstantaneousDrift[0].Drift.Clone();
                driftCurve.MultiplyBy(Percentage.OverPercentagePoint);
                graphData.Add(new ResultSeries("Drift", LineStyle.Dashed, driftCurve));
            }

            graphData.ToXml(output);
        }

        /// <summary>
        /// Get type of associated model
        /// </summary>
        /// <returns>CS forward price model type</returns>
        public Type ModelType()
        {
            return typeof(CSMultiFactorForwardPriceModel);
        }
    }
}