/// <author>
/// Alastair Wilkins
/// </author>
/// <owner>
/// Alastair Wilkins
/// </owner>
/// <summary>
/// Calibration of PCA volatility models.
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
    /// Generic base class for calibration of PCA volatility models.
    /// </summary>
    [Serializable]
    public abstract class PCAVolCalibration<TSurfaceClass, TPriceClass, TRiskClass> : StatisticsBasedCalibration
        where TSurfaceClass : VolSurface, new()
        where TPriceClass : class, IPriceFactor
        where TRiskClass : OrnsteinUhlenbeckProcess, new()
    {
        private int fNumberOfPcaFactors = 3;

        protected PCAVolCalibration()
        {
            Expiry_Grid    = new DisplayableList<Period>();
            Moneyness_Grid = new DisplayableList<DisplayableDouble>();
        }

        /// <summary>
        /// Number of factors for the PCA decomposition.
        /// </summary>
        public int Number_Of_PCA_Factors
        {
            get { return fNumberOfPcaFactors; } set { fNumberOfPcaFactors = Math.Max(value, 0); }
        }

        /// <summary>
        /// Expiry grid.
        /// </summary>
        public DisplayableList<Period> Expiry_Grid                  { get; set; }

        /// <summary>
        /// Money grid.
        /// </summary>
        public DisplayableList<DisplayableDouble> Moneyness_Grid    { get; set; }

        /// <summary>
        /// Returns the tenor grid.
        /// </summary>
        public virtual DisplayableList<Period> TenorGrid()
        {
            return null;
        }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            PCAVolModel<TSurfaceClass, TPriceClass, TRiskClass> pcaModel =
                CalibrationHelper.Validate<PCAVolModel<TSurfaceClass, TPriceClass, TRiskClass>>(calibrationData, priceModel, output);

            CalibratePcaVolModel(calibrationData, pcaModel, errorList);
            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// Calculate the PCA model parameters from the statistics of log changes in volatility.
        /// </summary>
        private void CalibratePcaVolModel(ICalibrationData calibrationData, PCAVolModel<TSurfaceClass, TPriceClass, TRiskClass> pcaModel, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(pcaModel);

            MeanReversionStatistics[] meanReversionStatistics;
            TimeSeriesStructure<string> timeSeries = null;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Log, null, null, null, null, false);
            bool isSuccess = GetMeanReversionStatistics(calibrationData, pcaModel,calculationParameters, errorList,
                                                        ref timeSeries, out meanReversionStatistics, out isValidStatistics);

            if (!isSuccess)
                return;

            int nPoints = meanReversionStatistics.Length;

            List<double> expiryGrid    = CalibrationHelper.ToDoubleList(Expiry_Grid, period => (double)period);
            List<double> tenorGrid     = CalibrationHelper.ToDoubleList(TenorGrid(), period => (double)period);
            List<double> moneynessGrid = CalibrationHelper.ToDoubleList(Moneyness_Grid, d => d.Value);

            var points          = new List<VolSurface.Point>(nPoints);
            var validIndexes    = new List<int>();
            for (int index = nPoints - 1; index >= 0; --index)
            {
                // Parse point string and construct point
                VolSurface.Point point = pcaModel.Historical_Surface.PointFromString(meanReversionStatistics[index].PointID);
                points.Add(point);

                if ((expiryGrid    == null || CalibrationHelper.IsInList(expiryGrid,    point.Expiry,    CalcUtils.WEEK_IN_YEARS)) &&
                    (tenorGrid     == null || CalibrationHelper.IsInList(tenorGrid,     point.Tenor,     CalcUtils.WEEK_IN_YEARS)) &&
                    (moneynessGrid == null || CalibrationHelper.IsInList(moneynessGrid, point.Moneyness, CalcUtils.SMALL)))
                {
                    // Surface point in statistics file is on (or close to) a point in the selected grid
                    validIndexes.Add(index);
                }
            }

            // Reverse list of points constructed in reverse order
            points.Reverse();
            validIndexes.Reverse();

            // Number of selected points with non-zero volatility
            if (nPoints < fNumberOfPcaFactors)
            {
                var withs = new List<string> {" non-zero volatility"};

                if (Expiry_Grid.Count > 0)
                    withs.Add(" expiry component in Expriy_Grid");

                if (TenorGrid() != null && TenorGrid().Count > 0)
                    withs.Add(" tenor component in Tenor_Grid");

                if (Moneyness_Grid.Count > 0)
                    withs.Add(" moneyness component in Moneyness_Grid");

                errorList.Add(ErrorLevel.Error, string.Format("Cannot calibrate {0}-factor {1} because only {2} surface points in statistics with{3}.", 
                    fNumberOfPcaFactors, pcaModel.GetKey(), nPoints, StringUtils.ArrayToString(withs.ToArray())));
                return;
            }

            // Calculate mean reversion speed and set long run mean surface
            double sumReversionSpeed = 0.0;
            pcaModel.Historical_Surface.Surface.Clear();
            int nValidPoints = 0;
            var reversionVols = new double[nPoints];
            for (int i = 0; i < nPoints; ++i)
            {
                reversionVols[i] = meanReversionStatistics[i].Sigma;

                // Use only valid indexes, i.e., points near a point on selected grid.
                if (!validIndexes.Contains(i))
                    continue;

                nValidPoints++;
                sumReversionSpeed += meanReversionStatistics[i].Alpha;
                pcaModel.Historical_Surface.SetValue(points[i], meanReversionStatistics[i].Theta);
            }

            pcaModel.Reversion_Speed = sumReversionSpeed / nValidPoints;

            // Construct the covariance matrix.
            // Calculate covariance matrix from historical data or read the pre-computed matrix.
            SymmetricMatrix correlationMatrix;
            isSuccess = Use_Pre_Computed_Statistics == YesNo.Yes
                                 ? GetCorrelationMatrix(calibrationData, pcaModel, out correlationMatrix)
                                 : CorrelationCalculationHelper.GetCorrelationMatrix(timeSeries, ReturnMethod.Log, out correlationMatrix);

            if (!isSuccess || correlationMatrix == null)
            {
                errorList.Add(ErrorLevel.Error, string.Format("No correlation matrix available for {0}.", priceFactor.GetKey()));
                return;
            }

            var covarianceMatrix = CorrelationCalculationHelper.GetCovarianceMatrix(correlationMatrix, reversionVols);

            // Create a covariance matrix for valid indexes.
            SymmetricMatrix validCovarianceMatrix;
            if (nValidPoints < nPoints)
            {
                validCovarianceMatrix = SymmetricMatrix.Create(nValidPoints);
                int validRowIndex = 0;
                for (int i = 0; i < nPoints; ++i)
                {
                    if (!validIndexes.Contains(i))
                        continue;

                    int validColumnIndex = 0;
                    for (int j = 0; j <= i; ++j)
                    {
                        if (!validIndexes.Contains(j))
                            continue;

                        validCovarianceMatrix[validRowIndex, validColumnIndex] = covarianceMatrix[i, j];
                        validColumnIndex++;
                    }

                    validRowIndex++;
                }
            }
            else
            {
                validCovarianceMatrix = covarianceMatrix;
            }

            // Eigenvector decomposition
            var eigenCalculateResult = validCovarianceMatrix.CalculateEigenVectorsAndValues();

            // Set eigenvector/eigenvalue model parameters
            pcaModel.Eigenvectors.Clear();
            for (int k = 0; k < fNumberOfPcaFactors; ++k)
            {
                pcaModel.Eigenvectors.Add(new Eigen<TSurfaceClass>());
                pcaModel.Eigenvectors[k].Eigenvalue = eigenCalculateResult.EigenValues[k];

                for (int index = 0; index < validIndexes.Count; index++)
                {
                    int validIndex = validIndexes[index];
                    pcaModel.Eigenvectors[k].Eigenvector.SetValue(points[validIndex], eigenCalculateResult.EigenVectors[index, k]);
                }
            }
        }
    }

    /// <summary>
    /// Calibration of PCA asset price volatility model.
    /// </summary>
    [DisplayName("PCA Asset Price Volatility Calibration")]
    public class PCAAssetPriceVolCalibration : PCAVolCalibration<AssetPriceVolSurface, ISpotProcessVol, AssetPriceVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(PCAAssetPriceVolModel);
        }
    }

    /// <summary>
    /// Calibration of PCA FX volatility model.
    /// </summary>
    [DisplayName("PCA FX Volatility Calibration")]
    public class PCAFXVolCalibration : PCAVolCalibration<AssetPriceVolSurface, ISpotProcessVol, AssetPriceVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(PCAFXVolModel);
        }
    }

    /// <summary>
    /// Base class for calibration of PCA interest volatility models.
    /// </summary>
    public abstract class PCAInterestVolCalibration<TRiskClass> : PCAVolCalibration<InterestVolSurface, IInterestVol, TRiskClass>
        where TRiskClass : OrnsteinUhlenbeckProcess, new()
    {
        protected PCAInterestVolCalibration()
        {
            Tenor_Grid = new DisplayableList<Period>();
        }

        public DisplayableList<Period> Tenor_Grid
        {
            get; set;
        }

        public override DisplayableList<Period> TenorGrid()
        {
            return Tenor_Grid;
        }
    }

    /// <summary>
    /// Calibration of PCA interest rate volatility model.
    /// </summary>
    [DisplayName("PCA Interest Rate Volatility Calibration")]
    public class PCAInterestRateVolCalibration : PCAInterestVolCalibration<InterestRateVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(PCAInterestRateVolModel);
        }
    }

    /// <summary>
    /// Calibration of PCA interest yield volatility model.
    /// </summary>
    [DisplayName("PCA Interest Yield Volatility Calibration")]
    public class PCAInterestYieldVolCalibration : PCAInterestVolCalibration<InterestYieldVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(PCAInterestYieldVolModel);
        }
    }
}