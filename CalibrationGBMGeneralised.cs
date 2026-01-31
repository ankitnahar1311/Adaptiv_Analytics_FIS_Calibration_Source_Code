using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    using HealingMethod = CorrelationMatrixHealer.CorrelationHealingMethod;

    /// <summary>
    /// Enum to determine which volatility calibration to use.
    /// </summary>
    public enum VolCalibration
    {
        Implied, Historical
    }

    /// <summary>
    /// Enum to determine which volatility calibration to use for FX Rate TS Extrapolation calibration.
    /// </summary>
    public enum VolCalibration2
    {
        Implied, Implied_Extrapolated, Historical
    }

    /// <summary>
    /// Enum to determine which drift calibration to use.
    /// </summary>
    public enum DriftCalibration
    {
        Historical_Drift, Historical_Interest_Rates
    }

    /// <summary>
    /// Interface for classes which store parameters for the GBMTS implied vol calibration.
    /// </summary>
    public interface IGBMTSModelImpliedParameters
    {
        /// <summary>
        /// Asset volatility.
        /// </summary>
        VolCurve Vol
        {
            get;
            set;
        }

        /// <summary>
        /// Quanto FX volatility.
        /// </summary>
        VolCurve Quanto_FX_Volatility
        {
            get;
            set;
        }

        /// <summary>
        /// Quanto FX correlation.
        /// </summary>
        double Quanto_FX_Correlation
        {
            get;
            set;
        }
    }

    /// <summary>
    /// Provides helping methods to all classes calibrating a GBM asset price TS model.
    /// </summary>
    public static class GBMAssetPriceTSCalibrationHelper
    {
        /// <summary>
        /// Name, type and currency properties of an asset price factor.
        /// </summary>
        public struct AssetPriceProperties
        {
            public string Name;
            public Type Type;
            public string ParentCurrency;
            public string DomesticCurrency;
            public string ForeignCurrency;
        }

        /// <summary>
        /// Update output curve for a model with term structure drift and vol.
        /// </summary>
        public static void CalibrationOutput(VolCurve drift, VolCurve vol, XmlWriter output)
        {
            ResultSeriesList graphData = new ResultSeriesList();
            var volatilityCurve = vol.Clone();
            volatilityCurve.MultiplyBy(Percentage.OverPercentagePoint);
            var driftCurve = drift.Clone();
            driftCurve.MultiplyBy(Percentage.OverPercentagePoint);

            graphData.fXAxis.fLabel = "Time (Years)";
            graphData.fYAxis.fLabel = "Volatility or Rate (%)";
            graphData.Add(new ResultSeries("Volatility", LineStyle.Solid, volatilityCurve));
            graphData.Add(new ResultSeries("Drift", LineStyle.Solid, driftCurve));
            graphData.ToXml(output);
        }

        /// <summary>
        /// Ensures that the supplied <see cref="SpotPrice" /> implements <see cref="IAssetPriceWithVolatilityPriceFactor"/>.
        /// </summary>
        public static IAssetPriceWithVolatilityPriceFactor ValidateAssetPrice(IPriceFactor assetPriceToBeEvolved, string nameOfUsingType)
        {
            var spotPriceAsISpotPriceWithVolatilityPriceFactor = assetPriceToBeEvolved as IAssetPriceWithVolatilityPriceFactor;
            if (spotPriceAsISpotPriceWithVolatilityPriceFactor == null)
                throw new AnalyticsException(string.Format("{0} cannot be used when the pricefactor being evolved is of type {1}", nameOfUsingType, assetPriceToBeEvolved.GetType().Name));

            return spotPriceAsISpotPriceWithVolatilityPriceFactor;
        }

        /// <summary>
        /// Register the appropriate volatility price factor (and possibly FXVols and correlations as well).
        /// </summary>
        /// <param name="priceFactorList">Register with this one.</param>
        /// <param name="priceFactors">Use these to calculate. As provided by ICalibrationData / IBootstrappingData.</param>
        /// <remarks>
        /// If an explicit vol ID has been supplied we use that. Otherwise:
        /// We first see if we can find a vol pricefactor, possibly with a different currency to the spot price.
        /// If we can't find any vol pricefactors we register one with the same currency as the spot price.
        /// We share a lot of logic with <see cref="VolCalibrationFromPriceFactor"/> unfortunately.
        /// </remarks>
        public static void RegisterVolatilityPriceFactor(IAssetPriceWithVolatilityPriceFactor priceFactor, PriceFactorList priceFactorList, PriceFactorList priceFactors,
                                                         YesNo riskNeutral, FactorID explicitVolID)
        {
            string baseCurrency = priceFactorList.BaseCcyCode;
            var asset = GetAssetPriceProperties(priceFactor, baseCurrency);

            if (RequireQuantoParameters(riskNeutral, baseCurrency, asset.DomesticCurrency, asset.ForeignCurrency))
            {
                FXVolHelper.Register(priceFactorList, asset.DomesticCurrency, baseCurrency);

                // For the quanto correction we require the correlation between 
                // the price of the asset in domestic currency (i.e. the asset price as quoted) and 
                // the price of domestic currency in base currency.
                if (asset.Type == typeof(FxRate))
                    CorrelationHelper.Register(priceFactorList, asset.Type, asset.ForeignCurrency, asset.DomesticCurrency, typeof(IFxRate), asset.DomesticCurrency, baseCurrency);
                else
                    CorrelationHelper.Register(priceFactorList, asset.Type, asset.Name, asset.DomesticCurrency, typeof(IFxRate), asset.DomesticCurrency, baseCurrency);
            }

            ISpotProcessVol volProcess;
            if (explicitVolID != null)
            {
                // The user supplied an explicit vol ID
                volProcess = priceFactor.RegisterVolatilityPriceFactorForCalibration(priceFactorList, explicitVolID);
            }
            else if (priceFactor.TryGetVolatilityPriceFactorForCalibration(priceFactors, asset.Name, asset.ParentCurrency, out volProcess))
            {
                // Found a volatility price factor in the price factor list for the underlying asset price or FX rate, register this price factor
                priceFactor.RegisterVolatilityPriceFactorForCalibration(priceFactorList, volProcess.GetID());
            }
            else
            {
                // Register a volatility price factor with the same currency as the asset price
                priceFactor.RegisterVolatilityPriceFactorForCalibration(priceFactorList, asset.Name, asset.ParentCurrency);
                return;
            }

            string requiredCurrency;
            AssetPriceProperties volAsset;
            if (NeedsTriangulation(asset, volProcess, out volAsset, out requiredCurrency))
            {
                // Register factors for the calculation of the volatility of the asset price in requiredCurrency.
                // Register the volAsset.DomesticCurrency / requiredCurrency FX volatility curve.
                FXVolHelper.Register(priceFactorList, volAsset.DomesticCurrency, requiredCurrency);

                // Register the correlation between asset price in volAsset.DomesticCurrency and price of volAsset.DomesticCurrency in requiredCurrency
                CorrelationHelper.Register(priceFactorList, volAsset.Type, volAsset.Name, volAsset.DomesticCurrency, typeof(IFxRate), volAsset.DomesticCurrency, requiredCurrency);
            }
        }

        /// <summary>
        /// Calibrate vol from the price factor.
        /// </summary>
        public static void VolCalibrationFromPriceFactor(IAssetPriceWithVolatilityPriceFactor priceFactor, PriceFactorList priceFactorList,
                                                         YesNo riskNeutral, IGBMTSModelImpliedParameters modelParameters, FactorID explicitVolID)
        {
            string baseCurrency = priceFactorList.BaseCcyCode;
            var asset = GetAssetPriceProperties(priceFactor, baseCurrency);

            if (RequireQuantoParameters(riskNeutral, baseCurrency, asset.DomesticCurrency, asset.ForeignCurrency))
            {
                IFXVol fxVol = FXVolHelper.Get(priceFactorList, asset.DomesticCurrency, baseCurrency).fFXVol;

                // For the quanto correction we require the correlation between 
                // the price of the asset in domestic currency (i.e. the asset price as quoted) and 
                // the price of domestic currency in base currency.
                CorrelationHelper correl;
                if (asset.Type == typeof(FxRate))
                    correl = CorrelationHelper.Get(priceFactorList, asset.Type, asset.ForeignCurrency, asset.DomesticCurrency, typeof(IFxRate), asset.DomesticCurrency, baseCurrency);
                else
                    correl = CorrelationHelper.Get(priceFactorList, asset.Type, asset.Name, asset.DomesticCurrency, typeof(IFxRate), asset.DomesticCurrency, baseCurrency);

                modelParameters.Quanto_FX_Volatility = new VolCurve() { Type = CurveType.Integrated };
                GetATMVolCurve(fxVol, modelParameters.Quanto_FX_Volatility);
                modelParameters.Quanto_FX_Correlation = correl.Value;
            }

            ISpotProcessVol volProcess;
            if (explicitVolID != null)
            {
                // The user supplied an explicit vol ID
                volProcess = priceFactor.GetVolatilityPriceFactorForCalibration(priceFactorList, explicitVolID);
            }
            else if (!priceFactor.TryGetVolatilityPriceFactorForCalibration(priceFactorList, asset.Name, asset.ParentCurrency, out volProcess))
            {
                // No volatility price factor in the price factor list for the underlying asset price or FX rate
                throw new CalibrationException(string.Format("Failed to obtain volatility price factor for asset price or FX Rate {0}.", asset.Name));
            }

            modelParameters.Vol.Clear();
            GetATMVolCurve(volProcess, modelParameters.Vol);

            string requiredCurrency;
            AssetPriceProperties volAsset;
            if (NeedsTriangulation(asset, volProcess, out volAsset, out requiredCurrency))
            {
                // Calculate the volatility of the asset price in requiredCurrency.
                // Get the volAsset.DomesticCurrency / requiredCurrency FX volatility curve.
                IFXVol fxVol = FXVolHelper.Get(priceFactorList, volAsset.DomesticCurrency, requiredCurrency).fFXVol;
                RateCurve fxVolCurve = new RateCurve();
                GetATMVolCurve(fxVol, fxVolCurve);

                // Get the correlation between asset price in volAsset.DomesticCurrency and price of volAsset.DomesticCurrency in requiredCurrency
                CorrelationHelper correl = CorrelationHelper.Get(priceFactorList, volAsset.Type, volAsset.Name, volAsset.DomesticCurrency, typeof(IFxRate), volAsset.DomesticCurrency, requiredCurrency);

                for (int i = 0; i < modelParameters.Vol.Count; ++i)
                {
                    double xi = modelParameters.Vol.X[i];
                    double vol = modelParameters.Vol.Y[i];
                    double fxv = fxVolCurve[xi];
                    double volSqr = fxv * fxv + vol * vol + 2.0 * fxv * vol * correl.Value;
                    modelParameters.Vol[xi] = Math.Sqrt(Math.Max(volSqr, 0.0));
                }
            }
        }

        /// <summary>
        /// Create drift curve from constant drift of log price and integrated volatility curve.
        /// </summary>
        public static void CreateDriftWithVolatilityCorrection(IModelWithDrift model, double drift, RateCurve integratedVol)
        {
            model.Drift.Clear();
            for (int i = 0; i < integratedVol.Count; ++i)
            {
                var t = integratedVol.X[i];
                var vol = integratedVol.Y[i];
                model.Drift[t] = drift + 0.5 * vol * vol;
            }
        }

        /// <summary>
        /// Extract the ATM vol curve from the vol surface.
        /// </summary>
        public static void GetATMVolCurve(ISpotProcessVol volatilityPriceFactor, RateCurve atmVolCurve)
        {
            #region Validation

            if (volatilityPriceFactor == null)
                throw new ArgumentNullException("volatilityPriceFactor", "volatilityPriceFactor is null.");
            if (volatilityPriceFactor != null && volatilityPriceFactor.GetSurface().Dimension != 2)
                throw new ArgumentNullException("volatilityPriceFactor", "volatilityPriceFactor.Surface.Dimension is not 2.");

            #endregion

            atmVolCurve.Clear();

            Surface surface = volatilityPriceFactor.GetSurface();
            if (surface.Count > 0)
            {
                // Create the coordinate the we will use to get values from the implied vol surface.
                var volSurfaceCoordinate = new double[2];

                // Set the ATM moneyness to 1.
                const double ATMMoneyness = 1.0;

                // Loop over each smile curve and get the ATM point on the smile curve
                for (int i = 0; i < surface.X.Count; ++i)
                {
                    Surface smileCurve = surface.Y[i];
                    if (smileCurve.Curve != null)
                    {
                        double tExpiry = surface.X[i];
                        volSurfaceCoordinate[0] = ATMMoneyness;
                        volSurfaceCoordinate[1] = tExpiry;

                        // Query the surface at the specified time-to-expiry and moneyness of 1.
                        double atmVol = surface.GetValueVarArgs(volSurfaceCoordinate);
                        atmVolCurve[tExpiry] = atmVol;
                    }
                }
            }
            else
            {
                atmVolCurve.Assign(volatilityPriceFactor.GetAtmVol());
            }
        }

        /// <summary>
        /// Returns true when riskNeutral is Yes (not setting Risk_Neutral_Drift or setting Risk_Neutral_Drift to Yes), 
        /// and domestic currency is not base currency, 
        /// and foreign currency is not base currency (for FX rate price factors). 
        /// </summary>
        private static bool RequireQuantoParameters(YesNo riskNeutral, string baseCurrency, string domesticCurrency, string foreignCurrency)
        {
            return riskNeutral == YesNo.Yes && domesticCurrency != baseCurrency && (foreignCurrency == null || foreignCurrency != baseCurrency);
        }

        /// <summary>
        /// Get the asset price properties.
        /// </summary>
        private static AssetPriceProperties GetAssetPriceProperties(IAssetPriceWithVolatilityPriceFactor priceFactor, string baseCurrency)
        {
            var app = new AssetPriceProperties();

            app.Type = priceFactor.GetType();
            app.Name = priceFactor.GetID().ToCode();

            app.ParentCurrency = priceFactor.ParentCurrency();
            if (string.IsNullOrEmpty(app.ParentCurrency))
                app.ParentCurrency = baseCurrency;

            app.DomesticCurrency = priceFactor.DomesticCurrency();
            if (string.IsNullOrEmpty(app.DomesticCurrency))
                app.DomesticCurrency = baseCurrency;

            app.ForeignCurrency = null;
            if (app.Type == typeof(FxRate))
            {
                app.ForeignCurrency = ((IFxRate)priceFactor).ForeignCurrency();
                if (string.IsNullOrEmpty(app.ForeignCurrency))
                    app.ForeignCurrency = app.Name;
            }

            return app;
        }

        /// <summary>
        /// Return true when the given asset price and the underlying asset of the volatility price factor have 
        /// a different currency (equity and commodity) or are not the same currency pair (FX).
        /// </summary>
        private static bool NeedsTriangulation(AssetPriceProperties asset, ISpotProcessVol volProcess, out AssetPriceProperties volAsset, out string requiredCurrency)
        {
            volAsset = new AssetPriceProperties() { Type = asset.Type };
            requiredCurrency = string.Empty;

            var assetVolProcess = volProcess as IAssetPriceVol;
            if (assetVolProcess != null)
            {
                string volCurrency = assetVolProcess.GetCurrency();
                requiredCurrency = asset.DomesticCurrency;
                if (volCurrency == requiredCurrency)
                    return false;  // no triangulation required

                volAsset.Name = assetVolProcess.GetAssetPriceId();
                volAsset.DomesticCurrency = volCurrency;
                volAsset.ParentCurrency = volCurrency;
                return true;
            }

            var fxVolProcess = volProcess as IFXVol;
            if (fxVolProcess != null)
            {
                string volCurrency1 = fxVolProcess.Currency1();
                string volCurrency2 = fxVolProcess.Currency2();

                // The currency pair {volCurrency1, volCurrency2} must contain at least one element of  
                // the currency pair {asset.ForeignCurrency, asset.DomesticCurrency}.
                if (volCurrency1 != asset.ForeignCurrency && volCurrency1 != asset.DomesticCurrency &&
                    volCurrency2 != asset.ForeignCurrency && volCurrency2 != asset.DomesticCurrency)
                    throw new CalibrationException(string.Format("Cannot calibrate FX rate {0} / {1} using volatility of {2} / {3}", asset.ForeignCurrency, asset.DomesticCurrency, volCurrency1, volCurrency2));

                // If the currency pairs are the same then no triangulation is required and we return false.
                // If the currency pairs have precisely one common currency then triangulation is required and 
                // we return true and the details of the underlying asset for the given fx volatility price factor.
                // In this case, requiredCurrency is the element of 
                // {asset.ForeignCurrency, asset.DomesticCurrency} that is not in {volCurrency1, volCurrency2}.
                // There are four cases of the currency pairs having precisely one common currency.
                if (volCurrency1 == asset.ForeignCurrency && volCurrency2 != asset.DomesticCurrency)
                {
                    requiredCurrency = asset.DomesticCurrency;
                    volAsset.ForeignCurrency = volCurrency1;
                    volAsset.DomesticCurrency = volCurrency2;
                }
                else if (volCurrency2 == asset.ForeignCurrency && volCurrency1 != asset.DomesticCurrency)
                {
                    requiredCurrency = asset.DomesticCurrency;
                    volAsset.ForeignCurrency = volCurrency2;
                    volAsset.DomesticCurrency = volCurrency1;
                }
                else if (volCurrency1 == asset.DomesticCurrency && volCurrency2 != asset.ForeignCurrency)
                {
                    requiredCurrency = asset.ForeignCurrency;
                    volAsset.ForeignCurrency = volCurrency1;
                    volAsset.DomesticCurrency = volCurrency2;
                }
                else if (volCurrency2 == asset.DomesticCurrency && volCurrency1 != asset.ForeignCurrency)
                {
                    requiredCurrency = asset.ForeignCurrency;
                    volAsset.ForeignCurrency = volCurrency2;
                    volAsset.DomesticCurrency = volCurrency1;
                }

                if (!string.IsNullOrEmpty(requiredCurrency))
                {
                    volAsset.Name = volAsset.ForeignCurrency;
                    volAsset.ParentCurrency = volAsset.DomesticCurrency;
                    return true;
                }
            }

            return false;  // no triangulation required
        }
    }

    /// <summary>
    /// Helper class for common multi factor calibration methods.
    /// </summary>
    public static class MultiFactorCalibrationHelper
    {
        /// <summary>
        /// Outputs the normalised systemic weights rounded to zero (if less than the significance threshold)
        /// or one (if there is only one systemic weight above the significance threshold such that the
        /// idiosynchratic factor's weight would be less than the significance threshold.
        /// </summary>
        /// <remarks>If Matrix is null, then correlation matrix is treated as a diagonal matrix with uncorrelated risk factors.</remarks>
        /// <remarks>If weights are normalised to 1.0, sumsquared must be set to 1.0.</remarks>
        public static void ObtainRoundedNormalisedRegressionWeights(double significanceThreshold, double sumSquared, YesNo idiosynchraticFactor, SymmetricMatrix driverCorrelations, Vector beta, out double epsilon)
        {
            double systematicWeightedSum = 0.0;

            for (int i = 0; i < beta.Count; ++i)
            {
                if (Math.Abs(beta[i]) < significanceThreshold)
                    beta[i] = 0.0;
                else
                {
                    if (sumSquared - CalcUtils.Sqr(beta[i]) < significanceThreshold) // in this case, all other risk factor weights are rounded to zero so this weight needs to be rounded up to sqrt(sumSquared).
                        beta[i] = Math.Sqrt(sumSquared) * Math.Sign(beta[i]);

                    systematicWeightedSum += CalcUtils.Sqr(beta[i]);
                }
            }

            // Add cross terms if systematic drivers are correlated
            if (driverCorrelations != null)
            {
                for (int i = 0; i < beta.Count; ++i)
                {
                    if (beta[i] == 0.0)
                        continue;

                    for (int j = 0; j < i; ++j)
                    {
                        systematicWeightedSum += 2.0 * driverCorrelations[i, j] * beta[i] * beta[j];
                    }
                }
            }

            epsilon = Math.Sqrt(Math.Max(sumSquared - systematicWeightedSum, 0.0));
            if (idiosynchraticFactor == YesNo.No || (epsilon < significanceThreshold && sumSquared != systematicWeightedSum))
            {
                epsilon = 0.0;
                Debug.Assert(systematicWeightedSum != 0.0, "Normalisation coefficient is null, cannot rescale the weight list");
                beta.MultiplyBy(Math.Sqrt(sumSquared / systematicWeightedSum));
            }
        }
    }

    /// <summary>
    /// Helper class for Regression Calibration.
    /// </summary>
    public static class RegressionCalibrationHelper
    {
        private const int MaximumNumberOfDriversForAlternatingProjectionHealingMethod = 200;

        /// <summary>
        /// Compute the parameters of a multivariate linear regression 
        /// driverCorrelations * beta = responseDriverCorrelations using Singular Value Decomposition.
        /// </summary>
        /// <param name="beta">Returns a vector beta as the solution to the linear regression.</param>
        /// <param name="driverCorrelations">Correlation between the drivers.</param>
        /// <param name="responseDriverCorrelations">Correlation between the response and the drivers.</param>
        /// <param name="numDrivers">Number of drivers.</param>
        /// <param name="errorList">List of calibration warnings and errors.</param>
        public static void ComputeParametersOfMultivariateRegression(out Vector beta, SymmetricMatrix driverCorrelations,
                Vector responseDriverCorrelations, int numDrivers, ErrorList errorList)
        {
            // Compute the correlation matrix by putting together the correlation of drivers
            // matrix and correlation between the response and drivers time series
            var fullCorrelation = SymmetricMatrix.Create(numDrivers + 1);
            for (int i = 0; i < numDrivers; ++i)
            {
                for (int j = i + 1; j < numDrivers; ++j)
                {
                    double rhoij = driverCorrelations[i, j];
                    fullCorrelation[i, j] = rhoij;
                }

                double gammai = responseDriverCorrelations[i];
                fullCorrelation[numDrivers, i] = gammai;
                fullCorrelation[i, i] = 1.0;
            }

            fullCorrelation[numDrivers, numDrivers] = 1.0;

            // Check that the full matrix is positive semi-definite
            bool isMatrixPositiveSemidefinite = driverCorrelations.IsPositiveSemiDefiniteSymmetric();

            // Attempt to heal matrix
            if (!isMatrixPositiveSemidefinite)
            {
                HealingMethod healingMethod = numDrivers < MaximumNumberOfDriversForAlternatingProjectionHealingMethod
                                        ? HealingMethod.Alternating_Projections : HealingMethod.Eigenvalue_Raising;
                var healedMatrix = CorrelationMatrixHealer.Heal(fullCorrelation, healingMethod, errorList);

                // Extract driver correlation matrix and response-driver correlation vector
                for (int i = 0; i < numDrivers; ++i)
                {
                    for (int j = i + 1; j < numDrivers; ++j)
                    {
                        driverCorrelations[i, j] = healedMatrix[i, j];
                    }

                    responseDriverCorrelations[i] = healedMatrix[numDrivers, i];
                    driverCorrelations[i, i] = 1.0;
                }
            }

            // Solve the matrix equation driverCorrelations * beta = responseDriverCorrelations.
            var svdResult = driverCorrelations.ToMatrix().CalculateSVD();
            beta = Matrix.SolveLinearSystemUsingSVD(svdResult, responseDriverCorrelations);
        }
    }

    /// <summary>
    /// Calibration base class of generalised single GBM model.
    /// </summary>
    public abstract class GBMAssetPriceTSCalibrationBase : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        protected GBMAssetPriceTSCalibrationBase()
        {
            Set_Risk_Neutral_Drift = YesNo.No;
            Risk_Neutral_Drift = YesNo.No;
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
        /// Calibrates the generalised single GBM model from the statistics set
        /// regardless of the risk neutral flag.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var model = CalibrationHelper.Validate<GBMAssetPriceTSModel>(calibrationData, priceModel, output);

            // set the risk neutral flag of the model if requested
            if (Set_Risk_Neutral_Drift == YesNo.Yes)
                model.Risk_Neutral_Drift = Risk_Neutral_Drift;

            Calculate(calibrationData, model, errorList);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public virtual void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var gbmTSModel = (GBMAssetPriceTSModel)priceModel;
            GBMAssetPriceTSCalibrationHelper.CalibrationOutput(gbmTSModel.Drift, gbmTSModel.Vol, output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        /// <returns>Type of price model.</returns>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(GBMAssetPriceTSModel);
        }

        protected abstract void Calculate(ICalibrationData sd, GBMAssetPriceTSModel model, ErrorList errorList);

        public void HistoricalVolCalibration(ICalibrationData sd, GBMAssetPriceTSModel model, ErrorList errorList,
            IAssetPriceWithVolatilityPriceFactor priceFactor)
        {
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            var priceFactorToEvolve = (IAssetPrice)priceFactor;
            bool[] isValidStatistics;
            if (GetMeanAndVolMomentStatistics(sd,
                new FactorTypeAndID(priceFactorToEvolve.GetType().Name, priceFactorToEvolve.GetID()), ReturnMethod.Log, errorList,
                ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                Debug.Assert(momentStatistics.Length == 1, "Asset prices are one dimensional.");

                if (isValidStatistics[0])
                {
                    model.Vol.AssignConstant(momentStatistics[0].Vol);
                }
                else
                {
                    errorList.Add(ErrorLevel.Error, "No valid historical volatility found.");
                }
            }
        }
    }

    /// <summary>
    /// Standard calibration for single factor GBM model with vol and drift term structure.
    /// </summary>
    [DisplayName("GBM Asset Price Term Structure Calibration")]
    public class GBMAssetPriceTSCalibration : GBMAssetPriceTSCalibrationBase
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public GBMAssetPriceTSCalibration()
        {
            Vol_Calibration = VolCalibration.Implied;
        }

        /// <summary>
        /// If set to No, volatility will be set to a constant volatility, calculated from historical vols.
        /// </summary>
        public VolCalibration Vol_Calibration { get; set; }

        /// <summary>
        /// Calibrates GBMAssetPriceTSModel from statistics set.
        /// </summary>
        protected override void Calculate(ICalibrationData sd, GBMAssetPriceTSModel model, ErrorList errorList)
        {
            var priceFactor = (IAssetPriceWithVolatilityPriceFactor)CalibrationHelper.GetPriceFactor(model);

            if (Vol_Calibration == VolCalibration.Implied)
            {
                GBMAssetPriceTSCalibrationHelper.VolCalibrationFromPriceFactor(priceFactor, sd.PriceFactors, model.Risk_Neutral_Drift, model, null);
            }
            else
            {
                HistoricalVolCalibration(sd, model, errorList, priceFactor);
            }

            // Calculate drift when not setting Risk_Neutral_Drift, or setting Risk_Neutral_Drift to No.
            if (Set_Risk_Neutral_Drift == YesNo.No || Risk_Neutral_Drift == YesNo.No)
            {
                double drift;
                if (TryGetDriftFromLogStatistics(sd, model, errorList, false, out drift))
                    GBMAssetPriceTSCalibrationHelper.CreateDriftWithVolatilityCorrection(model, drift, model.GetIntegratedVol());
            }
        }
    }

    /// <summary>
    /// Calibration of GBM asset price term structure model from archived market data for FX rates.
    /// </summary>
    [DisplayName("FX Rate Term Structure Extrapolation Calibration")]
    public class FXRateTSExtrapolationCalibration : GBMAssetPriceTSCalibrationBase
    {
        protected Period fBetaPeriod = new Period("3Y");
        protected Period fAlphaPeriod = new Period("1M");
        protected double fBetaTol = CalcUtils.SMALL;

        public FXRateTSExtrapolationCalibration()
        {
            Calibration_Date = DateTime.Now;
        }

        public VolCalibration2 Vol_Calibration
        {
            get;
            set;
        }

        public DriftCalibration Drift_Calibration
        {
            get;
            set;
        }

        public Period Calibration_Period
        {
            get { return fBetaPeriod; }
            set { SetBetaPeriod(value); }
        }

        public Period Averaging_Period
        {
            get { return fAlphaPeriod; }
            set { SetAlphaPeriod(value); }
        }

        public double Convergence_Tolerance
        {
            get { return fBetaTol; }
            set { SetBetaTol(value); }
        }

        public TDate Calibration_Date
        {
            get;
            set;
        }

        /// <summary>
        /// Calibrate the volatility from the Market Data Archive.
        /// </summary>
        public void VolCalibrationFromMDA(IHistoricalArchive mda, GBMAssetPriceTSModel model, FactorID fxVolID)
        {
            FXVol currentVol = new FXVol();
            currentVol.fID = fxVolID;
            // We need to set a value on the surface to ensure the GetValue call in GetImpliedVolsFromMDA uses the 2D surface representation
            currentVol.Surface.SetValue(new double[] { 1.0, 1.0 }, CalcUtils.MinVolatility);

            string[] volPointIDs = mda.GetPointIDs(FXVol.TypeDotSubType(typeof(FXVol), FXVol.Vanilla), fxVolID);
            if (volPointIDs.Length == 0)
                throw new CalibrationException(String.Format("No entries in market data archive for FXVol {0}.", fxVolID));

            List<double> terms = new List<double>();
            foreach (string point in volPointIDs)
            {
                string[] parts = point.Split(',');
                terms.Add(GetExactTerm(parts[1]));
            }
            terms.Sort();

            // impliedVols indexed as [date][term]
            List<List<double>> impliedVols = GetImpliedVolsFromMDA(mda, fxVolID, volPointIDs, terms, currentVol);

            VolCurve volCurve = new VolCurve();
            AddAverageVolsToModel(impliedVols, terms, volCurve);

            List<double> alphas = new List<double>();
            double beta = 0.1;
            double oldBeta = double.MaxValue;

            int numDatePoints = impliedVols.Count;
            int numIterations = 0;

            while (Math.Abs(beta - oldBeta) > fBetaTol)
            {
                oldBeta = beta;
                alphas.Clear();

                for (int i = 0; i < numDatePoints; ++i)
                {
                    alphas.Add(SolveForAlpha(i, impliedVols, beta, terms));
                }

                beta = CalculateBeta(impliedVols, alphas, terms);
                numIterations++;

                if (numIterations > 10000)
                    throw new CalibrationException(string.Format("Calibration of {0} failed to converge.", model.GetKey()));
            }

            double alpha = CalculateAverageAlpha(alphas);

            AddExtrapolatedVolsToModel(alpha, beta, volCurve);

            model.Vol = volCurve;
        }

        /// <summary>
        /// Calibrate the drift from the Market Data Archive.
        /// </summary>
        public void DriftCalibrationFromMDA(IHistoricalArchive mda, GBMAssetPriceTSModel model, ICalibrationData sd, FactorID foreignInterestRateId, FactorID domesticInterestRateId, bool isSpotInverted)
        {
            InterestRate foreignInterestRate = new InterestRate();
            InterestRate domesticInterestRate = new InterestRate();

            int index = sd.MarketData.PriceFactors.IndexOf(typeof(InterestRate), foreignInterestRateId);
            if ((index >= 0 && sd.MarketData.PriceFactors[index].GetPriceModel() != null && sd.MarketData.PriceFactors[index].GetPointIDs().Length > 0)
                || (index >= 0 && sd.MarketData.PriceFactors[index].GetPriceModel() == null))
            {
                foreignInterestRate = (InterestRate)sd.MarketData.PriceFactors[index].Clone();
            }
            else if (index >= 0)
            {
                // If there exists an InterestRate with no data, then remove it.
                // Then Get the required InterestRate which will load data on demand.
                // Replace the old InterestRate back so as not to change the price factors during the calibration.
                InterestRate temp = (InterestRate)sd.MarketData.PriceFactors[index];
                sd.MarketData.PriceFactors.Remove(temp);
                foreignInterestRate = (InterestRate)sd.MarketData.PriceFactors.Get<InterestRate>(foreignInterestRateId).Clone();
                sd.MarketData.PriceFactors.Add(temp);
            }
            else
            {
                foreignInterestRate = (InterestRate)sd.MarketData.PriceFactors.Get<InterestRate>(foreignInterestRateId).Clone();
            }

            index = sd.MarketData.PriceFactors.IndexOf(typeof(InterestRate), domesticInterestRateId);
            if ((index >= 0 && sd.MarketData.PriceFactors[index].GetPriceModel() != null && sd.MarketData.PriceFactors[index].GetPointIDs().Length > 0)
                || (index >= 0 && sd.MarketData.PriceFactors[index].GetPriceModel() == null))
            {
                domesticInterestRate = (InterestRate)sd.MarketData.PriceFactors[index].Clone();
            }
            else if (index >= 0)
            {
                // If there exists an InterestRate with no data, then remove it.
                // Then Get the required InterestRate which will load data on demand.
                // Replace the old InterestRate back so as not to change the price factors during the calibration.
                InterestRate temp = (InterestRate)sd.MarketData.PriceFactors[index];
                sd.MarketData.PriceFactors.Remove(temp);
                domesticInterestRate = (InterestRate)sd.MarketData.PriceFactors.Get<InterestRate>(domesticInterestRateId).Clone();
                sd.MarketData.PriceFactors.Add(temp);
            }
            else
            {
                domesticInterestRate = (InterestRate)sd.MarketData.PriceFactors.Get<InterestRate>(domesticInterestRateId).Clone();
            }

            string[] foreignPointIDs = mda.GetPointIDs(InterestRate.TypeDotSubType(typeof(InterestRate), InterestRate.Vanilla), foreignInterestRateId);
            if (foreignPointIDs.Length == 0)
                throw new CalibrationException(String.Format("No entries in market data archive for InterestRate {0}.", foreignInterestRateId));

            string[] domesticPointIDs = mda.GetPointIDs(InterestRate.TypeDotSubType(typeof(InterestRate), InterestRate.Vanilla), domesticInterestRateId);
            if (domesticPointIDs.Length == 0)
                throw new CalibrationException(String.Format("No entries in market data archive for InterestRate {0}.", domesticInterestRateId));

            List<double> forwardRateTerms = (isSpotInverted ? domesticPointIDs : foreignPointIDs).Select(GetExactTerm).ToList();
            forwardRateTerms.Sort();

            List<List<double>> forwardRates = GetForwardRatesFromMDA(foreignInterestRateId, domesticInterestRateId, mda, foreignPointIDs, domesticPointIDs, foreignInterestRate, domesticInterestRate, forwardRateTerms);

            AddAverageDriftToModel(forwardRates, forwardRateTerms, model);
        }

        /// <summary>
        /// Calibrates GBMAssetPriceTSModel for FX rates using either archived market data or statistics.
        /// </summary>
        protected override void Calculate(ICalibrationData sd, GBMAssetPriceTSModel model, ErrorList errorList)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);
            Type factorType = priceFactor.GetType();

            if (factorType != typeof(FxRate))
                throw new CalibrationException(string.Format("{0} cannot calibrate model of {1}.", GetType().Name, factorType.Name));

            var fxRate = sd.PriceFactors.GetInterface<IFxRate>(priceFactor.GetID());
            var foreignCurrency = fxRate.ForeignCurrency();
            var domesticCurrency = fxRate.DomesticCurrency();
            var isSpotInverted = fxRate.SpotIsInverted();

            FactorTypeAndIDList factors = new FactorTypeAndIDList();
            FactorID fxVolID = null;
            FactorID foreignInterestRateId = null;
            FactorID domesticInterestRateId = null;

            if (Drift_Calibration == DriftCalibration.Historical_Interest_Rates)
            {
                var parentFxRate = sd.PriceFactors.GetInterface<IFxRate>(fxRate.ParentCurrency());

                if (isSpotInverted)
                {
                    foreignInterestRateId = parentFxRate.FactorRate().GetID();
                    domesticInterestRateId = fxRate.FactorRate().GetID();
                }
                else
                {
                    foreignInterestRateId = fxRate.FactorRate().GetID();
                    domesticInterestRateId = parentFxRate.FactorRate().GetID();
                }

                string typeName = typeof(InterestRate).Name;
                factors.Add(new FactorTypeAndID(typeName, foreignInterestRateId));
                factors.Add(new FactorTypeAndID(typeName, domesticInterestRateId));
            }

            if (Vol_Calibration == VolCalibration2.Implied_Extrapolated)
            {
                // If invert is true the moneyness axis is inverted, but we only use moneyness = 1
                bool invert;
                fxVolID = FXVol.MakeCanonicalID(foreignCurrency, domesticCurrency, out invert);
                factors.Add(new FactorTypeAndID(typeof(FXVol).Name, fxVolID));
            }

            IHistoricalArchive mda = null;
            if (factors.Count > 0)
            {
                DateTime startDate = DateAdjuster.Subtract(DateTime.FromOADate(Calibration_Date), Period.ValueToTerm(fBetaPeriod));
                mda = sd.GetHistoricalMarketRates(factors, startDate, Calibration_Date);
            }

            var priceFactoAsAssetPriceWithVol = (IAssetPriceWithVolatilityPriceFactor)priceFactor;

            if (Vol_Calibration == VolCalibration2.Implied_Extrapolated)
            {
                VolCalibrationFromMDA(mda, model, fxVolID);
            }
            else if (Vol_Calibration == VolCalibration2.Implied)
            {
                GBMAssetPriceTSCalibrationHelper.VolCalibrationFromPriceFactor(priceFactoAsAssetPriceWithVol, sd.PriceFactors, model.Risk_Neutral_Drift, model, null);
            }
            else
            {
                // Historical
                HistoricalVolCalibration(sd, model, errorList, priceFactoAsAssetPriceWithVol);
            }

            if (Drift_Calibration == DriftCalibration.Historical_Interest_Rates)
            {
                DriftCalibrationFromMDA(mda, model, sd, foreignInterestRateId, domesticInterestRateId, isSpotInverted);
            }
            else if (Set_Risk_Neutral_Drift == YesNo.No || Risk_Neutral_Drift == YesNo.No)
            {
                double drift;
                if (TryGetDriftFromLogStatistics(sd, model, errorList, false, out drift))
                    GBMAssetPriceTSCalibrationHelper.CreateDriftWithVolatilityCorrection(model, drift, model.GetIntegratedVol());
            }
        }

        private static double SolveForAlpha(int date, List<List<double>> impliedVols, double beta, List<double> terms)
        {
            // using bisection

            double x_lower = 0.0;
            double x_upper = 1.0;
            double x_mid = 0.0;
            double f_upper = 0.0;
            double f_mid = 0.0;
            int numIterations = 0;

            f_upper = Function(x_upper, date, impliedVols, beta, terms);

            while (f_upper < 0.0)
            {
                x_upper *= 2.0;
                f_upper = Function(x_upper, date, impliedVols, beta, terms);
            }

            while (Math.Abs(x_lower - x_upper) > 1E-4) // alpha tolerance
            {
                x_mid = (x_lower + x_upper) / 2.0;
                f_mid = Function(x_mid, date, impliedVols, beta, terms);

                if (f_mid < 0.0)
                    x_lower = x_mid;
                else
                    x_upper = x_mid;

                numIterations++;

                if (numIterations > 10000)
                    throw new CalibrationException("The calibration has failed to converge.");
            }
            x_mid = (x_lower + x_upper) / 2.0;

            if (x_mid == 0.0)
                throw new CalibrationException("The calibration has failed to converge.");

            return x_mid;
        }

        private static double Function(double x_n, int date, List<List<double>> impliedVols, double beta, List<double> terms)
        {
            double function = 0.0;

            for (int j = 0; j < terms.Count; ++j)
            {
                double impliedVol = impliedVols[date][j];
                double term = terms[j];
                double exponential = Math.Exp(-1.0 * x_n * term);

                if (term == 0.0)
                    continue;

                function += ((beta / term) * (1 - exponential) - impliedVol * impliedVol) * exponential;
            }
            return function;
        }

        private static double CalculateBeta(List<List<double>> impliedVols, List<double> alphas, List<double> terms)
        {
            double numerator = 0.0;
            double denominator = 0.0;

            for (int i = 0; i < alphas.Count; ++i)
            {
                for (int j = 0; j < terms.Count; ++j)
                {
                    double alpha = alphas[i];
                    double term = terms[j];
                    double impliedVol = impliedVols[i][j];

                    if (term == 0.0)
                        continue;

                    double exponential = Math.Exp(-1.0 * alpha * term);
                    double oneMinusExpOverTerm = (1.0 - exponential) / term;

                    numerator += oneMinusExpOverTerm * impliedVol * impliedVol;
                    denominator += oneMinusExpOverTerm * oneMinusExpOverTerm;
                }
            }

            return numerator / denominator;
        }

        private static void AddExtrapolatedVolsToModel(double alpha, double beta, RateCurve volCurve)
        {
            int numVolPoints = volCurve.X.Count;
            if (numVolPoints <= 0)
                return;

            double lastTime = volCurve.X[numVolPoints - 1];
            double lastVol = volCurve.Y[numVolPoints - 1];

            List<double> termsToAdd = new List<double>();

            int year = (int)Math.Ceiling(lastTime + 1);

            while (year < 10)
            {
                termsToAdd.Add(year);
                ++year;
            }
            while (year <= 30)
            {
                termsToAdd.Add(year);
                year += 5;
            }
            termsToAdd.Add(100);

            double delta = 0.0;

            foreach (int term in termsToAdd)
            {
                delta = Math.Exp(-0.1 * (term - (lastTime)));

                double value = beta * (1.0 - Math.Exp(-1.0 * alpha * term)) + delta * (lastVol * lastVol * lastTime - beta * (1.0 - Math.Exp(-1.0 * alpha * lastTime)));

                if (Double.IsNaN(value))
                    continue;

                numVolPoints = volCurve.X.Count;

                double vol = Math.Sqrt(value / term);
                double previousVol = volCurve.Y[numVolPoints - 1];
                double previousTerm = volCurve.X[numVolPoints - 1];

                double varianceCondition = 3.0 * term * vol - 2.0 * term * previousVol - previousTerm * vol;

                if (varianceCondition >= 0.0)
                {
                    value = value / term;
                    value = Math.Sqrt(value);

                    volCurve.X.Add(term);
                    volCurve.Y.Add(value);
                }
            }
        }

        private static void AddAverageDriftToModel(List<List<double>> forwardRates, List<double> forwardRateTerms, GBMAssetPriceTSModel model)
        {
            model.Drift.Clear();

            int numTerms = forwardRateTerms.Count;
            int numPointsToAverage = forwardRates.Count;

            for (int i = 0; i < numTerms; ++i)
            {
                double term = forwardRateTerms[i];
                double driftTotal = 0.0;

                if (term == 0.0)
                    continue;

                for (int j = 0; j < numPointsToAverage; ++j)
                {
                    double drift = Math.Log(forwardRates[j][i]);

                    driftTotal += drift;
                }

                model.Drift.X.Add(term);
                model.Drift.Y.Add((driftTotal / (numPointsToAverage * term)));
            }
        }

        private void SetBetaPeriod(Period value)
        {
            if (value > fAlphaPeriod)
                fBetaPeriod = value;
            else
                throw new ApplicationException("The Beta period must be longer than the Alpha period.");
        }

        private void SetAlphaPeriod(Period value)
        {
            if (value < fBetaPeriod)
                fAlphaPeriod = value;
            else
                throw new ApplicationException("The Alpha period must be shorter than the Beta period.");
        }

        private void SetBetaTol(double value)
        {
            if (value > 1E-11)
                fBetaTol = value;
            else
                throw new ApplicationException("The convergence tolerance must be greater than 1E-11");
        }

        private List<List<double>> GetImpliedVolsFromMDA(IHistoricalArchive mda, FactorID factorID, string[] volPointIDs, List<double> terms, FXVol currentVol)
        {
            List<List<double>> impliedVols = new List<List<double>>();

            DateTime[] availableDates = mda.HistoricalDates;
            DateTime endDate = DateTime.FromOADate(Math.Min(availableDates[availableDates.Length - 1].ToOADate(), Calibration_Date));
            DateTime startDate = DateAdjuster.Subtract(endDate, Period.ValueToTerm(fBetaPeriod));

            FXVol volSurface = currentVol;
            volSurface.fParent = new PriceFactorList();
            using (var cache = Vector.Cache(1))
            {
                Vector one = cache.Get(1.0);
                Vector vec = cache.Get();

                foreach (DateTime date in availableDates)
                {
                    if (date < startDate)
                        continue;

                    ScenarioAssetPriceVol deterministicModel = new ScenarioAssetPriceVol(1, 2);

                    deterministicModel.SetPriceFactor(volSurface);
                    deterministicModel.fParent = volSurface.fParent;

                    foreach (string pointID in volPointIDs)
                    {
                        int index = mda.IndexOfItem(FXVol.TypeDotSubType(typeof(FXVol), FXVol.Vanilla), factorID, pointID);
                        double vol = 0.0F;
                        bool succeeded = mda.TryGetNonMissingValue(index, date, true, out vol);

                        if (!succeeded)
                        {
                            if (index < 0)
                                throw new CalibrationException(String.Format("No value in market data archive for FXVol {0} {1}.", factorID, pointID));

                            vol = CalcUtils.MinVolatility;
                        }

                        vec.Assign(vol);
                        deterministicModel.SetPointValue(0.0F, pointID, vec);
                    }

                    volSurface.SetPriceModel(deterministicModel);

                    List<double> impliedVolsByTerm = new List<double>();

                    foreach (double term in terms)
                    {
                        volSurface.GetValue(vec, 0.0, one, one, term);
                        impliedVolsByTerm.Add(vec[0]);
                    }

                    impliedVols.Add(impliedVolsByTerm);
                }
            }

            return impliedVols;
        }

        private double CalculateAverageAlpha(List<double> alphas)
        {
            int numPointsToAverage = (int)Math.Ceiling(5.0 * Period.ValueToTerm(fAlphaPeriod).ToDays() / 7.0); // business days adjustment
            int numPoints = alphas.Count;
            double alphaTotal = 0.0;

            for (int i = numPoints - 1; i > numPoints - numPointsToAverage - 1 && i >= 0; i--)
            {
                alphaTotal += alphas[i];
            }

            return alphaTotal / Math.Min(numPointsToAverage, numPoints);
        }

        private void AddAverageVolsToModel(List<List<double>> impliedVols, List<double> terms, RateCurve volCurve)
        {
            volCurve.Clear();

            int numTerms = terms.Count;
            int numPointsToAverage = (int)Math.Ceiling(5.0 * Period.ValueToTerm(fAlphaPeriod).ToDays() / 7.0); // business days adjustment
            int numPoints = impliedVols.Count;

            for (int i = 0; i < numTerms; ++i)
            {
                double impliedVolSquareTotal = 0.0;
                double term = terms[i];

                for (int j = numPoints - 1; j > numPoints - numPointsToAverage - 1 && j >= 0; j--)
                {
                    impliedVolSquareTotal += CalcUtils.Sqr(impliedVols[j][i]);
                }

                volCurve.X.Add(term);
                volCurve.Y.Add(Math.Sqrt(impliedVolSquareTotal / Math.Min(numPointsToAverage, numPoints)));
            }
        }

        private double GetExactTerm(string term)
        {
            double archiveTerm = CalcUtils.ParsePoint(term);
            int numDays = Convert.ToInt32(archiveTerm * CalcUtils.DAYS_IN_YEAR);
            return numDays * CalcUtils.DAY_IN_YEARS;
        }

        private List<List<double>> GetForwardRatesFromMDA(FactorID foreignFactorID, FactorID domesticFactorID, IHistoricalArchive mda, string[] foreignIRPointIDs, string[] domesticIRPointIDs, InterestRate currentForeignRate, InterestRate currentDomesticRate, List<double> forwardRateTerms)
        {
            List<List<double>> forwardRates = new List<List<double>>();

            DateTime[] availableDates = mda.HistoricalDates;
            DateTime endDate = DateTime.FromOADate(Math.Min(availableDates[availableDates.Length - 1].ToOADate(), Calibration_Date));
            DateTime startDate = DateAdjuster.Subtract(endDate, Period.ValueToTerm(fAlphaPeriod));

            InterestRate foreignInterestRate = currentForeignRate;
            foreignInterestRate.fParent = new PriceFactorList();

            InterestRate domesticInterestRate = currentDomesticRate;
            domesticInterestRate.fParent = new PriceFactorList();

            foreach (DateTime date in availableDates)
            {
                if (date < startDate)
                    continue;

                List<double> forwardRatesByTerm = new List<double>();

                SetScenarioInterestRateValues(foreignFactorID, foreignInterestRate, mda, date, foreignIRPointIDs);
                SetScenarioInterestRateValues(domesticFactorID, domesticInterestRate, mda, date, domesticIRPointIDs);

                foreach (double term in forwardRateTerms)
                {
                    using (var cache = Vector.Cache(1))
                    {
                        Vector domesticRatePoint = cache.Get();
                        Vector forwardRatePoint = cache.Get();

                        domesticInterestRate.GetValue(domesticRatePoint, 0, term);
                        foreignInterestRate.GetValue(forwardRatePoint, 0, term);

                        double forwardRate = forwardRatePoint[0] / domesticRatePoint[0];
                        forwardRatesByTerm.Add(forwardRate);
                    }
                }

                forwardRates.Add(forwardRatesByTerm);
            }

            return forwardRates;
        }

        private static void SetScenarioInterestRateValues(FactorID factorID, InterestRate interestRate, IHistoricalArchive mda, DateTime date, IEnumerable<string> pointIDs)
        {
            ScenarioInterestRate scenarioInterestRate = new ScenarioInterestRate(1, 2);
            scenarioInterestRate.SetPriceFactor(interestRate);
            interestRate.SetPriceModel(scenarioInterestRate);
            Vector vector = Vector.Create(1);

            foreach (string pointID in pointIDs)
            {
                int index = mda.IndexOfItem(InterestRate.TypeDotSubType(typeof(InterestRate), InterestRate.Vanilla), factorID, pointID);
                double rate;
                bool succeeded = mda.TryGetNonMissingValue(index, date, true, out rate);

                if (!succeeded)
                {
                    if (index < 0)
                        throw new CalibrationException(String.Format("No value in market data archive for InterestRate {0} {1}.", factorID, pointID));

                    rate = CalcUtils.MinInterestRate;
                }

                vector.Assign(rate);
                scenarioInterestRate.SetPointValue(0, pointID, vector);
            }
        }
    }
}
