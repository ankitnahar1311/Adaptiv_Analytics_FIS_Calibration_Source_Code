using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;
using SunGard.Adaptiv.Core;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Interface for classes that store Heston model parameters (excluding the Drift and Risk_Neutral_Drift flag).
    /// </summary>
    public interface IHestonModelParameters
    {
        /// <summary>
        /// Mean reversion speed for instantaneous variance.
        /// </summary>
        double Variance_Reversion_Speed
        {
            get;
            set;
        }

        /// <summary>
        /// Mean reversion level for instantaneous variance.
        /// </summary>
        double Variance_Reversion_Level
        {
            get;
            set;
        }

        /// <summary>
        /// Volatility for instantaneous variance.
        /// </summary>
        double Variance_Vol
        {
            get;
            set;
        }

        /// <summary>
        /// Risk premium for instantaneous variance.
        /// </summary>
        double Variance_Risk_Premium
        {
            get;
            set;
        }

        /// <summary>
        /// Correlation between the asset price and variance processes.
        /// </summary>
        double Spot_Variance_Correlation
        {
            get;
            set;
        }

        /// <summary>
        /// Initial Variance.
        /// </summary>
        double Initial_Variance
        {
            get;
            set;
        }
    }

    /// <summary>
    /// Base class for Heston Calibration/Bootstrapping.
    /// </summary>
    /// <typeparam name="TModelParameters">The type of model parameters (could be a model or a model parameters object).</typeparam>
    [Serializable]
    public abstract class HestonImpliedVolatilityCalibrationBase<TModelParameters>
        where TModelParameters : BaseFactor, IHestonModelParameters
    {
        protected CalibrationParameters fCalibrationParameters;

        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        protected HestonImpliedVolatilityCalibrationBase()
        {
            // Default calibration settings
            Lower_Time_To_Expiry = new Period("3M");
            Upper_Time_To_Expiry = new Period("10y");
            Lower_Moneyness = 0.5;
            Upper_Moneyness = 2.0;
            Use_Last_Parameters_As_Initial_Guesses = YesNo.Yes;
            Use_Penalty = YesNo.No;

            // Default parameter limits
            Variance_Reversion_Speed_Lower_Bound  = 0.01;
            Variance_Reversion_Speed_Upper_Bound  = 10.0;
            Variance_Reversion_Level_Lower_Bound  = 0.005;
            Variance_Reversion_Level_Upper_Bound  = 1.0;
            Variance_Vol_Lower_Bound              = 0.005;
            Variance_Vol_Upper_Bound              = 2.0;
            Spot_Variance_Correlation_Lower_Bound = -1.0;
            Spot_Variance_Correlation_Upper_Bound = 1.0;
            Initial_Variance_Lower_Bound          = 0.005;
            Initial_Variance_Upper_Bound          = 1.0;

            Fix_Initial_Variance = YesNo.No;
            Optimiser_Fractional_Tolerance = CalcUtils.SMALL;
        }

        /// <summary>
        /// Lower time to expiry for implied vol set.
        /// </summary>
        public Period Lower_Time_To_Expiry { get; set; }

        /// <summary>
        /// Upper time to expiry for implied vol set.
        /// </summary>
        public Period Upper_Time_To_Expiry { get; set; }

        /// <summary>
        /// Lower moneyness for implied vol set.
        /// </summary>
        public double Lower_Moneyness { get; set; }

        /// <summary>
        /// Upper moneyness for implied vol set.
        /// </summary>
        public double Upper_Moneyness { get; set; }

        /// <summary>
        /// Use last parameters as initial guesses flag.
        /// </summary>
        public YesNo Use_Last_Parameters_As_Initial_Guesses { get; set; }

        /// <summary>
        /// Use penalty function flag.
        /// </summary>
        public YesNo Use_Penalty { get; set; }

        /// <summary>
        /// Lower bound for Variance Reversion Speed.
        /// </summary>
        public double Variance_Reversion_Speed_Lower_Bound { get; set; }

        /// <summary>
        /// Upper bound for Variance Reversion Speed.
        /// </summary>
        public double Variance_Reversion_Speed_Upper_Bound { get; set; }

        /// <summary>
        /// Lower bound for Variance Reversion Level.
        /// </summary>
        public double Variance_Reversion_Level_Lower_Bound { get; set; }

        /// <summary>
        /// Upper bound for Variance Reversion Level.
        /// </summary>
        public double Variance_Reversion_Level_Upper_Bound { get; set; }

        /// <summary>
        /// Lower bound for Variance Vol.
        /// </summary>
        public double Variance_Vol_Lower_Bound { get; set; }

        /// <summary>
        /// Upper bound for Variance Vol.
        /// </summary>
        public double Variance_Vol_Upper_Bound { get; set; }

        /// <summary>
        /// Lower bound for Spot_Variance_Correlation.
        /// </summary>
        public double Spot_Variance_Correlation_Lower_Bound { get; set; }

        /// <summary>
        /// Upper bound for Spot_Variance_Correlation.
        /// </summary>
        public double Spot_Variance_Correlation_Upper_Bound { get; set; }

        /// <summary>
        /// Lower bound for Initial Variance.
        /// </summary>
        public double Initial_Variance_Lower_Bound { get; set; }

        /// <summary>
        /// Upper bound for Initial Variance.
        /// </summary>
        public double Initial_Variance_Upper_Bound { get; set; }

        /// <summary>
        /// Flag to fix the initial variance or let it be part of the optimization.
        /// </summary>
        public YesNo Fix_Initial_Variance { get; set; }

        /// <summary>
        /// Fractional tolerance used in the optimizer.
        /// </summary>
        public double Optimiser_Fractional_Tolerance { get; set; }

        /// <summary>
        /// Calibrate the Heston model to implied volatilities.
        /// </summary>
        public bool TryCalculate(PriceFactorList priceFactors, ICalibrationData calibrationData, IAssetPrice priceFactorToEvolve, HestonAssetPriceModel modelForCalibration, TModelParameters modelParameters, ErrorList errorList)
        {
            // Calibration needs access to the asset vol price factor and built spot price
            AssetPriceAndVol assetPriceAndVol = GetSpotPriceAndVolatilityPriceFactor(priceFactors, priceFactorToEvolve);
            modelForCalibration.SetPriceFactor(assetPriceAndVol.BuiltAssetPrice);

            // We create a temporary model to make use of the CalculateUndiscountedDomesticEuropeanOptionPrice() routine.
            // We calibrate the model, and then copy the Property values to our parameters pricefactor.
            ValidateParametersBound(errorList);
            if (errorList.MaxLevel() >= ErrorLevel.Error)
            {
                return false;
            }

            // Set the drift on the model.
            SetDrift(calibrationData, modelForCalibration, errorList);

            // Dimension of the optimization problem
            const int OptimizationDimension = 5;

            // Optimizer maximum iteration
            const int MaxIterations = 2000;

            // Initialise the calibration parameters object.
            fCalibrationParameters = new CalibrationParameters();
            fCalibrationParameters.Initialise(this, modelForCalibration, errorList, OptimizationDimension, assetPriceAndVol.VolPriceFactor, assetPriceAndVol.Inverted);

            // Initialize necessary downhill simplex objects.
            int valuationCount;
            double distance;
            double[] modelMin = new double[OptimizationDimension];
            double[,] simplex = new double[OptimizationDimension + 1, OptimizationDimension];
            double[] xMin = new double[OptimizationDimension];

            fCalibrationParameters.TransformInitialParametersOnR(xMin);
            fCalibrationParameters.PopulateSimplex(simplex, xMin);

            bool success = Minimizer<CalibrationParameters>.MinimizeDownhillSimplex(HestonOptionPriceDistance, fCalibrationParameters, simplex, MaxIterations, Optimiser_Fractional_Tolerance, out xMin, out distance, out valuationCount);
            if (!success)
            {
                errorList.Add(ErrorLevel.Error, "Exceeded number of iterations in downhill simplex method for calibration of " + modelParameters.GetKey() + ".");
                return false;
            }

            fCalibrationParameters.InverseTransformParameters(xMin, modelMin);

            // Set model properties
            SetModelParameters(modelMin, modelParameters);
            return true;
        }

        /// <summary>
        /// Registers the spot and volatility pricefactors.
        /// </summary>
        protected static void RegisterAssetAndVolatilityPriceFactor(PriceFactorList priceFactorList, IAssetPrice assetPrice)
        {
            var equityPrice = assetPrice as IEquityPrice;
            if (equityPrice != null)
            {
                RegisterAssetAndVolatilityPriceFactor<IEquityPrice, IEquityPriceVol>(priceFactorList, equityPrice);
                return;
            }

            var commodityPrice = assetPrice as ICommodityPrice;
            if (commodityPrice != null)
            {
                RegisterAssetAndVolatilityPriceFactor<ICommodityPrice, ICommodityPriceVol>(priceFactorList, commodityPrice);
                return;
            }

            var fxRate = assetPrice as IFxRate;
            if (fxRate != null)
            {
                RegisterFXAndVolatilityPriceFactor(priceFactorList, fxRate);
                return;
            }

            throw new ArgumentException(string.Format("Heston calibration: Asset price should be an EquityPrice, CommodityPrice or FxRate. Actual type: {0}", assetPrice.GetType().Name), "assetPrice");
        }

        /// <summary>
        /// Gets the built spot price and volatility pricefactor.
        /// </summary>
        protected static AssetPriceAndVol GetSpotPriceAndVolatilityPriceFactor(PriceFactorList priceFactorList, IAssetPrice assetPrice)
        {
            var equityPrice = assetPrice as IEquityPrice;
            if (equityPrice != null)
            {
                return GetAssetAndVolatilityPriceFactor<IEquityPrice, IEquityPriceVol>(priceFactorList, equityPrice);
            }

            var commodityPrice = assetPrice as ICommodityPrice;
            if (commodityPrice != null)
            {
                return GetAssetAndVolatilityPriceFactor<ICommodityPrice, ICommodityPriceVol>(priceFactorList, commodityPrice);
            }

            var fxRate = assetPrice as IFxRate;
            if (fxRate != null)
            {
                return GetFXAndVolatilityPriceFactor(priceFactorList, fxRate);
            }

            throw new ArgumentException(string.Format("Heston calibration: Spot price should be an EquityPrice, CommodityPrice or FxRate. Actual type: {0}", assetPrice.GetType().Name), "assetPrice");
        }

        /// <summary>
        /// Set the drift on the model.
        /// </summary>
        protected abstract void SetDrift(ICalibrationData calibrationData, HestonAssetPriceModel model, ErrorList errorList);
        
        /// <summary>
        /// Registers the spot and volatility pricefactors.
        /// </summary>
        /// <remarks>
        /// For Equities and Commodities.
        /// </remarks>
        /// <typeparam name="TAsset">The type of asset price.</typeparam>
        /// <typeparam name="TVol">The type of volatility pricefactor.</typeparam>
        private static void RegisterAssetAndVolatilityPriceFactor<TAsset, TVol>(PriceFactorList priceFactorList, TAsset assetPrice)
            where TAsset : class, IAssetPrice
            where TVol : class, ISpotProcessVol
        {
            var builtAssetPrice =  priceFactorList.RegisterInterface<TAsset>(assetPrice.GetID()[0]);
            priceFactorList.RegisterInterface<TVol>(new[] { builtAssetPrice.GetID()[0], builtAssetPrice.ParentCurrency() });
        }

        /// <summary>
        /// Gets the spot and volatility pricefactors.
        /// </summary>
        /// <remarks>
        /// For Equities and Commodities.
        /// </remarks>
        /// <typeparam name="TAsset">The type of asset price.</typeparam>
        /// <typeparam name="TVol">The type of volatility pricefactor.</typeparam>
        private static AssetPriceAndVol GetAssetAndVolatilityPriceFactor<TAsset, TVol>(PriceFactorList priceFactorList, TAsset assetPrice)
            where TAsset : class, IAssetPrice
            where TVol : class, ISpotProcessVol
        {
            var builtAssetPrice = priceFactorList.GetInterface<TAsset>(assetPrice.GetID()[0]);
            var spotPriceVol = priceFactorList.GetInterface<TVol>(new[] { builtAssetPrice.GetID()[0], builtAssetPrice.ParentCurrency() });
            return new AssetPriceAndVol(builtAssetPrice, spotPriceVol);
        }

        /// <summary>
        /// Registers the FxRate and FXVol pricefactors.
        /// </summary>
        private static void RegisterFXAndVolatilityPriceFactor(PriceFactorList priceFactorList, IFxRate spotPrice)
        {
            var builtSpotPrice = priceFactorList.RegisterInterface<IFxRate>(spotPrice.GetID()[0]);
            bool inverted;
            priceFactorList.RegisterInterface<IFXVol>(FXVol.MakeCanonicalID(builtSpotPrice.GetID()[0], builtSpotPrice.ParentCurrency(), out inverted));
        }

        /// <summary>
        /// Gets the FXRate and FXVol pricefactors.
        /// </summary>
        private static AssetPriceAndVol GetFXAndVolatilityPriceFactor(PriceFactorList priceFactorList, IFxRate spotPrice)
        {
            var builtSpotPrice = priceFactorList.GetInterface<IFxRate>(spotPrice.GetID()[0]);

            // Retrieve the FX Vol in the canonical form for currency pair defined by the FxRate ID and modelling currency
            bool inverted;
            var fxVol = priceFactorList.GetInterface<IFXVol>(FXVol.MakeCanonicalID(builtSpotPrice.GetID()[0], builtSpotPrice.ParentCurrency(), out inverted));
            return new AssetPriceAndVol(builtSpotPrice, fxVol, inverted);
        }

        /// <summary>
        /// Set model parameters from an array of parameters.
        /// </summary>
        private static void SetModelParameters(double[] parameters, IHestonModelParameters modelParameters)
        {
            modelParameters.Variance_Reversion_Speed  = parameters[0];
            modelParameters.Variance_Reversion_Level  = parameters[1];
            modelParameters.Variance_Vol              = parameters[2];
            modelParameters.Spot_Variance_Correlation = parameters[3];
            modelParameters.Initial_Variance          = parameters[4];
        }

        /// <summary>
        /// Objective function for Heston model.
        /// </summary>
        /// <remarks>
        /// The function is the square difference of out of the money option between the current model parameters 
        /// and observed option price weighted by the option's vega.
        /// The objective function can have a penalty equal to the Euclidian norm squared of the difference between the current 
        /// solution and the initial solution, this setting forces the method to converge to a solution close to the initial guess.
        /// This setting is best used in conjunction with Use_Last_Parameters_As_Initial_Guesses = Yes.
        /// </remarks>
        private static bool HestonOptionPriceDistance(double[] xParameters, CalibrationParameters calibrationParameters, out double distance)
        {
            HestonAssetPriceModel hestonModel = calibrationParameters.HestonModel;

            double[] inverseParams = new double[xParameters.Length];

            // Inverse transform to retrieve the current solution on the original parameter space
            calibrationParameters.InverseTransformParameters(xParameters, inverseParams);
            SetModelParameters(inverseParams, calibrationParameters.HestonModel);

            distance = 0.0;

            if (calibrationParameters.AddPenalty == YesNo.Yes)
                distance += CalcUtils.NormSqr(inverseParams, calibrationParameters.InitialParameters);

            using (var cache = Vector.Cache(calibrationParameters.NumberOfScenarios))
            {
                Vector spotModel = cache.Get();
                Vector strikeModel = cache.Get(1.0);
                Vector sign = cache.Get();
                Vector forward = cache.Get();
                Vector modelOptionPrice = cache.Get();

                Vector observedOptionPrice = cache.Get();
                Vector d1 = cache.Get();
                Vector d1Indicator = cache.Get();

                Vector strikeConstant = cache.Get();
                Vector stdDev = cache.Get();

                Vector fitError = cache.GetClear();
                Vector vega = cache.Get();

                foreach (ObservedImpliedVol vol in calibrationParameters.ObservedImpliedVols)
                {
                    spotModel.Assign(vol.Moneyness);
                    strikeModel.Assign(1.0);

                    // Calculate the out of the money option price under Heston for the set of parameters
                    hestonModel.CalculateUndiscountedDomesticEuropeanOptionPrice(modelOptionPrice, forward, 0.0, spotModel, strikeModel, vol.Time_To_Expiry, sign, false);

                    // Calculate option price under Black Scholes for the observed implied vol
                    stdDev.Assign(vol.Observed_Implied_Vol * Math.Sqrt(vol.Time_To_Expiry));
                    strikeConstant.Assign(-sign * strikeModel);
                    PricingFunctions.BlackFunction(observedOptionPrice, d1Indicator, d1, sign, strikeConstant, sign, forward, strikeModel, stdDev);

                    // vega = forward * CalcUtils.nd(d1) * sqrt(tau)
                    vega.Assign(forward * CalcUtils.nd(d1) * Math.Sqrt(vol.Time_To_Expiry));
                    vega.AssignConditional(vega > 0.01, vega, 0.01);

                    // Weigh the option price difference with vega
                    fitError.Add(VectorMath.Sqr((modelOptionPrice - observedOptionPrice) / vega));
                }

                distance += fitError[0];
            }

            return true;
        }

        /// <summary>
        /// Validate parameter limits.
        /// </summary>
        private void ValidateParametersBound(ErrorList errorList)
        {
            if (Lower_Moneyness > Upper_Moneyness)
                errorList.Add(ErrorLevel.Error, string.Format("Upper_Moneyness must be greater than Lower_Moneyness for {0}.", GetType().Name));

            if (Lower_Moneyness < 0.0)
                errorList.Add(ErrorLevel.Error, string.Format("Lower_Moneyness must be positive for {0}.", GetType().Name));

            if (Lower_Time_To_Expiry > Upper_Time_To_Expiry)
                errorList.Add(ErrorLevel.Error, string.Format("Upper_Time_To_Expiry must be greater than Lower_Time_To_Expiry for {0}.", GetType().Name));

            if (Variance_Reversion_Speed_Lower_Bound < 0 || Variance_Reversion_Speed_Upper_Bound < Variance_Reversion_Speed_Lower_Bound)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Reversion_Speed_Lower_Bound cannot be negative or greater than Variance_Reversion_Speed_Upper_Bound for {0}.", GetType().Name));

            if (Variance_Reversion_Speed_Upper_Bound <= 0.0)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Reversion_Speed_Upper_Bound must be positive for {0}.", GetType().Name));

            if (Variance_Reversion_Level_Lower_Bound < 0.0 || Variance_Reversion_Level_Upper_Bound < Variance_Reversion_Level_Lower_Bound)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Reversion_Level_Lower_Bound cannot be negative or greater than Variance_Reversion_Level_Upper_Bound for {0}.", GetType().Name));

            if (Variance_Reversion_Level_Upper_Bound < 0.0)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Reversion_Level_Upper_Bound must be positive for {0}.", GetType().Name));

            if (Variance_Vol_Lower_Bound < 0.0 || Variance_Vol_Upper_Bound < Variance_Vol_Lower_Bound)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Vol_Lower_Bound cannot be negative or greater than Variance_Vol_Upper_Bound for {0}.", GetType().Name));

            if (Variance_Vol_Upper_Bound < 0.0)
                errorList.Add(ErrorLevel.Error, string.Format("Variance_Vol_Upper_Bound must be positive for {0}.", GetType().Name));

            if (Spot_Variance_Correlation_Upper_Bound < Spot_Variance_Correlation_Lower_Bound)
                errorList.Add(ErrorLevel.Error, string.Format("Spot_Variance_Correlation_Upper_Bound must be greater than or equal to Correlation_Lower_Bound for {0}.", GetType().Name));

            if (Spot_Variance_Correlation_Lower_Bound < -1.0)
                errorList.Add(ErrorLevel.Error, string.Format("Spot_Vol_Correlation_Lower_Bound must be greater than or equal to -1.0 for {0}.", GetType().Name));

            if (Spot_Variance_Correlation_Upper_Bound > 1.0)
                errorList.Add(ErrorLevel.Error, string.Format("Spot_Vol_Correlation_Upper_Bound must be less than or equal to 1.0 for {0}.", GetType().Name));

            if (Initial_Variance_Lower_Bound < 0.0 || Initial_Variance_Upper_Bound < Initial_Variance_Lower_Bound)
                errorList.Add(ErrorLevel.Error, string.Format("Initial_Variance_Lower_Bound cannot be negative or greater than Initial_Variance_Upper_Bound for {0}.", GetType().Name));

            if (Initial_Variance_Upper_Bound < 0.0)
                errorList.Add(ErrorLevel.Error, string.Format("Initial_Variance_Upper_Bound must be positive for {0}.", GetType().Name));
        }

        /// <summary>
        /// Implied Vol set.
        /// </summary>
        [Serializable]
        protected struct ObservedImpliedVol
        {
            /// <summary>
            /// Initialies a new instance.
            /// </summary>
            public ObservedImpliedVol(double timeToExpiry, double moneyness, double observedImpliedVol)
                : this()
            {
                Time_To_Expiry       = timeToExpiry;
                Moneyness            = moneyness;
                Observed_Implied_Vol = observedImpliedVol;
            }

            /// <summary>
            /// Implied vol time to expiry.
            /// </summary>
            public double Time_To_Expiry
            {
                get;
                private set;
            }

            /// <summary>
            /// Implied vol moneyness.
            /// </summary>
            public double Moneyness
            {
                get;
                private set;
            }

            /// <summary>
            /// Observed implied vol.
            /// </summary>
            public double Observed_Implied_Vol
            {
                get;
                private set;
            }
        }

        /// <summary>
        /// Holds the built spot price and volatility pricefactor (and the "inverted" flag - for FxRates only).
        /// </summary>
        protected struct AssetPriceAndVol
        {
            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public AssetPriceAndVol(IAssetPrice assetPrice, ISpotProcessVol spotProcessVol, bool inverted = false)
                : this()
            {
                Ensure.ArgumentNotNull(assetPrice, "assetPrice");
                Ensure.ArgumentNotNull(spotProcessVol, "spotProcessVol");
                
                BuiltAssetPrice = assetPrice;
                VolPriceFactor = spotProcessVol;
                Inverted       = inverted;
            }

            /// <summary>
            /// The built spot price.
            /// </summary>
            public IAssetPrice BuiltAssetPrice
            {
                get;
                private set;
            }

            /// <summary>
            /// The volatility price factor.
            /// </summary>
            public ISpotProcessVol VolPriceFactor
            {
                get;
                private set;
            }

            /// <summary>
            /// Flag indicating whether the FXVol surface is inverted.
            /// </summary>
            public bool Inverted
            {
                get;
                private set;
            }
        }

        /// <summary>
        /// Calibration parameters passed to the optimization algorithm.
        /// </summary>
        [Serializable]
        protected class CalibrationParameters
        {
            private double fShortTermATMVariance = 0.0;

            /// <summary>
            /// Default parameters.
            /// </summary>
            public CalibrationParameters()
            {
                NumberOfScenarios = -1;
                ObservedImpliedVols = new List<ObservedImpliedVol>();
            }

            /// <summary>
            /// List of observed implied vol used in the calibration.
            /// </summary>
            public List<ObservedImpliedVol> ObservedImpliedVols { get; private set; }

            /// <summary>
            /// Reference to price model to calibrate.
            /// </summary>
            public HestonAssetPriceModel HestonModel { get; private set; }

            /// <summary>
            /// Initial parameters.
            /// </summary>
            public double[] InitialParameters { get; private set; }

            /// <summary>
            /// Constants used to transform parameters from their original search space on to the real axis.
            /// </summary>
            public double[,] TransformationConstants { get; private set; }

            /// <summary>
            /// Number of scenarios.
            /// </summary>
            public int NumberOfScenarios { get; private set; }

            /// <summary>
            /// Add Penalty flag.
            /// </summary>
            public YesNo AddPenalty { get; private set; }

            /// <summary>
            /// Initialise calibration parameters.
            /// </summary>
            public void Initialise(HestonImpliedVolatilityCalibrationBase<TModelParameters> calibration, HestonAssetPriceModel priceModel, ErrorList errorList, int dim, ISpotProcessVol volPriceFactor, bool inverted)
            {
                AddPenalty = calibration.Use_Penalty;

                NumberOfScenarios = priceModel.GetPriceFactor().GetParent().NumScenarios;

                HestonModel = priceModel;
                HestonModel.InitialiseCurrentState(NumberOfScenarios);

                ReadImpliedVolatility(volPriceFactor, calibration.Lower_Time_To_Expiry, calibration.Upper_Time_To_Expiry, calibration.Lower_Moneyness, calibration.Upper_Moneyness, inverted);

                SetSearchSpace(calibration, dim);

                ValidateInitialParameters(calibration, priceModel, errorList);
            }

            /// <summary>
            /// Transform initial set of model parameters into the transform parameters used in the optimization method.
            /// </summary>
            public void TransformInitialParametersOnR(double[] transformedParameters)
            {
                for (int i = 0; i < InitialParameters.Length; ++i)
                {
                    if (TransformationConstants[i, 0] == 0.0)
                        transformedParameters[i] = InitialParameters[i];
                    else
                        transformedParameters[i] = Math.Tan(Math.PI / TransformationConstants[i, 0] * (InitialParameters[i] - TransformationConstants[i, 1]));
                }
            }

            /// <summary>
            /// Retrieve the model parameters by inverse transformation.
            /// </summary>
            public void InverseTransformParameters(double[] parameters, double[] inverseTransformedParameters)
            {
                for (int i = 0; i < parameters.Length; ++i)
                {
                    if (TransformationConstants[i, 0] == 0.0)
                        inverseTransformedParameters[i] = InitialParameters[i];
                    else
                        inverseTransformedParameters[i] = Math.Atan(parameters[i]) * TransformationConstants[i, 0] / Math.PI + TransformationConstants[i, 1];
                }
            }

            /// <summary>
            /// Populate the given simplex with one vertex at the given parameter vector.
            /// </summary>
            /// <param name="simplex">The simplex to populate</param>
            /// <param name="x">The initial guesses of the parameter values</param>
            public void PopulateSimplex(double[,] simplex, double[] x)
            {
                const double Scale = 0.5;

                // Populate each vertex
                for (int v = 0; v < simplex.GetLength(0); v++)
                {
                    // Populate each dimension of the current vertex
                    for (int d = 0; d < simplex.GetLength(1); d++)
                    {
                        simplex[v, d] = x[d];
                    }

                    // Vertexes are offset from x by Scale in one dimension.
                    // If TransformationConstants[d, 0] == 0.0 then all simple[v, d] are the same.
                    if (v < simplex.GetLength(1) && TransformationConstants[v, 0] != 0.0)
                        simplex[v, v] += Scale;
                }
            }

            /// <summary>
            /// Calculate transformation parameters to map the search space imposed by lower and upper limits on model parameters
            /// on to R^dim where dim is the dimension of the optimization problem.
            /// Set initial guesses for model parameters.
            /// </summary>
            private void SetSearchSpace(HestonImpliedVolatilityCalibrationBase<TModelParameters> calibration, int dim)
            {
                TransformationConstants = new double[dim, 2];

                TransformationConstants[0, 0] = calibration.Variance_Reversion_Speed_Upper_Bound - calibration.Variance_Reversion_Speed_Lower_Bound;
                TransformationConstants[0, 1] = 0.5 * (calibration.Variance_Reversion_Speed_Upper_Bound + calibration.Variance_Reversion_Speed_Lower_Bound);

                TransformationConstants[1, 0] = calibration.Variance_Reversion_Level_Upper_Bound - calibration.Variance_Reversion_Level_Lower_Bound;
                TransformationConstants[1, 1] = 0.5 * (calibration.Variance_Reversion_Level_Upper_Bound + calibration.Variance_Reversion_Level_Lower_Bound);

                TransformationConstants[2, 0] = calibration.Variance_Vol_Upper_Bound - calibration.Variance_Vol_Lower_Bound;
                TransformationConstants[2, 1] = 0.5 * (calibration.Variance_Vol_Upper_Bound + calibration.Variance_Vol_Lower_Bound);

                TransformationConstants[3, 0] = calibration.Spot_Variance_Correlation_Upper_Bound - calibration.Spot_Variance_Correlation_Lower_Bound;
                TransformationConstants[3, 1] = 0.5 * (calibration.Spot_Variance_Correlation_Upper_Bound + calibration.Spot_Variance_Correlation_Lower_Bound);

                if (calibration.Fix_Initial_Variance == YesNo.No)
                {
                    TransformationConstants[4, 0] = calibration.Initial_Variance_Upper_Bound - calibration.Initial_Variance_Lower_Bound;
                    TransformationConstants[4, 1] = 0.5 * (calibration.Initial_Variance_Upper_Bound + calibration.Initial_Variance_Lower_Bound);
                }
                else
                {
                    // force the search interval to be [fShortTermATMVariance, fShortTermATMVariance]
                    TransformationConstants[4, 0] = 0.0;
                    TransformationConstants[4, 1] = fShortTermATMVariance;
                }

                InitialParameters = new double[dim];

                for (int i = 0; i < dim; ++i)
                {
                    // Set the initial default guesses to the mid point in the interval
                    InitialParameters[i] = TransformationConstants[i, 1];
                }
            }

            /// <summary>
            /// Validate initial parameters.
            /// </summary>
            /// <remarks>
            /// If we use last model parameters and the parameter is not within the search interval then use the middle of the interval instead.
            /// </remarks>
            private void ValidateInitialParameters(HestonImpliedVolatilityCalibrationBase<TModelParameters> calibration, HestonAssetPriceModel priceModel, ErrorList errorList)
            {
                if (calibration.Use_Last_Parameters_As_Initial_Guesses == YesNo.No)
                    return;

                if (TransformationConstants[0, 0] != 0.0)
                {
                    if (priceModel.Variance_Reversion_Speed >= calibration.Variance_Reversion_Speed_Upper_Bound || priceModel.Variance_Reversion_Speed <= calibration.Variance_Reversion_Speed_Lower_Bound)
                    {
                        errorList.Add(ErrorLevel.Info, string.Format("Initial guess for Mean_Reversion_Level is set the mid point of the search interval because Use_Last_Parameters_As_Initial_Guesses is set to Yes and Mean_Reversion_Level is outside the search interval for {0}.", priceModel.GetKey()));
                    }
                    else
                    {
                        InitialParameters[0] = priceModel.Variance_Reversion_Speed;
                    }
                }

                if (TransformationConstants[1, 0] != 0.0)
                {
                    if (priceModel.Variance_Reversion_Level >= calibration.Variance_Reversion_Level_Upper_Bound || priceModel.Variance_Reversion_Level <= calibration.Variance_Reversion_Level_Lower_Bound)
                    {
                        errorList.Add(ErrorLevel.Info, string.Format("Initial guess for Reversion_Level is set the mid point of the search interval because Use_Last_Parameters_As_Initial_Guesses is set to Yes and Reversion_Level is outside the search interval for {0}.", priceModel.GetKey()));
                    }
                    else
                    {
                        InitialParameters[1] = priceModel.Variance_Reversion_Level;
                    }
                }

                if (TransformationConstants[2, 0] != 0.0)
                {
                    if (priceModel.Variance_Vol >= calibration.Variance_Vol_Upper_Bound || priceModel.Variance_Vol <= calibration.Variance_Vol_Lower_Bound)
                    {
                        errorList.Add(ErrorLevel.Info, string.Format("Initial guess for Vol Of Vol is set the mid point of the search interval because Use_Last_Parameters_As_Initial_Guesses is set to Yes and Vol Of Vol is outside the search interval for {0}.", priceModel.GetKey()));
                    }
                    else
                    {
                        InitialParameters[2] = priceModel.Variance_Vol;
                    }
                }

                if (TransformationConstants[3, 0] != 0.0)
                {
                    if (priceModel.Spot_Variance_Correlation >= calibration.Spot_Variance_Correlation_Upper_Bound || priceModel.Spot_Variance_Correlation <= calibration.Spot_Variance_Correlation_Lower_Bound)
                    {
                        errorList.Add(ErrorLevel.Info, string.Format("Initial guess for Spot Vol Correlation is set the mid point of the search interval because Use_Last_Parameters_As_Initial_Guesses is set to Yes and Spot Vol Correlation is outside the search interval for {0}.", priceModel.GetKey()));
                    }
                    else
                    {
                        InitialParameters[3] = priceModel.Spot_Variance_Correlation;
                    }
                }

                if (TransformationConstants[4, 0] != 0.0)
                {
                    if (priceModel.Initial_Variance >= calibration.Initial_Variance_Upper_Bound || priceModel.Initial_Variance <= calibration.Initial_Variance_Lower_Bound)
                    {
                        errorList.Add(ErrorLevel.Info, string.Format("Initial guess for Initial Variance is set the mid point of the search interval because Use_Last_Parameters_As_Initial_Guesses is set to Yes and Initial Variance is outside the search interval for {0}.", priceModel.GetKey()));
                    }
                    else
                    {
                        InitialParameters[4] = priceModel.Initial_Variance;
                    }
                }
            }

            /// <summary>
            /// Read the implied volatility from the price factor vol surface between lower time to expiry 
            /// and upper time to expiry and between lower moneyness and upper moneyness.
            /// </summary>
            /// <remarks>
            /// The calibration routine for FX is set up for the currency pair defined on the price factor 
            /// i.e. where the foreign currency is the FX rate ID and the domestic currency is the modelling currency.
            /// When it is not in the canonical form, the calibration properties are assumed to be defined for this 
            /// currency pair and the moneyness is adjusted to use the canonical form.
            /// </remarks>
            private void ReadImpliedVolatility(ISpotProcessVol volPriceFactor, Period lowerTimeToExpiry, Period upperTimeToExpiry, double lowerMoneyness, double upperMoneyness, bool inverted)
            {
                Surface impliedVolSurface = volPriceFactor.GetSurface();

                double canonicalLowerMoneyness = inverted ? (1.0 / upperMoneyness) : lowerMoneyness;
                double canonicalUpperMoneyness = inverted ? (1.0 / lowerMoneyness) : upperMoneyness;

                // Create the coordinate that we will use to get values from the implied vol surface.
                var volSurfaceCoordinate = new double[2];

                for (int timeIndex = 0; timeIndex < impliedVolSurface.X.Count; ++timeIndex)
                {
                    double timeToExpiry = impliedVolSurface.X[timeIndex];

                    if (timeToExpiry < lowerTimeToExpiry)
                        continue;

                    if (timeToExpiry > upperTimeToExpiry)
                        break;

                    DoubleCurve smile = impliedVolSurface.Y[timeIndex].Curve;

                    for (int moneynessIndex = 0; moneynessIndex < smile.X.Count; ++moneynessIndex)
                    {
                        double moneyness = smile.X[moneynessIndex];
                        volSurfaceCoordinate[0] = moneyness;
                        volSurfaceCoordinate[1] = timeToExpiry;

                        double impliedVol = impliedVolSurface.GetValueVarArgs(volSurfaceCoordinate);

                        if (moneyness > canonicalUpperMoneyness)
                            break;

                        if (moneyness >= canonicalLowerMoneyness && CalcUtils.MinVolatility <= impliedVol && impliedVol <= CalcUtils.MaxVolatility)
                        {
                            ObservedImpliedVols.Add(new ObservedImpliedVol(timeToExpiry, inverted ? 1.0 / moneyness : moneyness, impliedVol));
                        }
                    }
                }

                if (ObservedImpliedVols.Count == 0)
                    throw new AnalyticsException("Calibration of " + HestonModel.GetKey() + " to an empty set of implied volatilities. No values were found on " + volPriceFactor.GetKey() + " for moneyness between " + canonicalLowerMoneyness + " and " + canonicalUpperMoneyness + " and for time to expiry between " + lowerTimeToExpiry + " and " + upperTimeToExpiry + ".");

                // index of the shortest time to expiry
                int firstExpiryIndex = impliedVolSurface.X.IndexOf(ObservedImpliedVols[0].Time_To_Expiry);

                // Index of the closest to ATM vol
                int atmIndex = impliedVolSurface.Y[firstExpiryIndex].Curve.X.Locate(1.0f);

                // Take the first and only point in the smile if necessary
                atmIndex = atmIndex > 0 ? atmIndex : 0;

                double shortestTimeToExpiry = impliedVolSurface.X[firstExpiryIndex];
                double moneynessForShortestTimeToExpiry = impliedVolSurface.Y[firstExpiryIndex].Curve.X[atmIndex];

                // Store the variance for this point.
                volSurfaceCoordinate[0] = moneynessForShortestTimeToExpiry;
                volSurfaceCoordinate[1] = shortestTimeToExpiry;
                fShortTermATMVariance = CalcUtils.Sqr(impliedVolSurface.GetValueVarArgs(volSurfaceCoordinate));
            }
        }
    }

    /// <summary>
    /// Heston calibration to implied volatilities.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get any statistics required from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate any statistics required from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks> 
    public class HestonImpliedVolatilityCalibration : HestonImpliedVolatilityCalibrationBase<HestonAssetPriceModel>, IModelCalibration
    {
        private readonly HestonDriftCalibrator fDriftCalibrator = new HestonDriftCalibrator();

        /// <summary>
        /// Use the statistics set if Yes; otherwise use archive market data.
        /// </summary>
        public YesNo Use_Pre_Computed_Statistics
        {
            get
            {
                return fDriftCalibrator.Use_Pre_Computed_Statistics;
            }
            set
            {
                fDriftCalibrator.Use_Pre_Computed_Statistics = value;
            }
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<string> Data_Retrieval_Parameters
        {
            get
            {
                return fDriftCalibrator.Data_Retrieval_Parameters;
            }
            set
            {
                fDriftCalibrator.Data_Retrieval_Parameters = value;
            }
        }

        /// <summary>
        /// Set risk neutral flag.
        /// </summary>
        public YesNo Set_Risk_Neutral_Drift
        {
            get;
            set;
        }

        /// <summary>
        /// Risk neutral flag.
        /// </summary>
        public YesNo Risk_Neutral_Drift
        {
            get;
            set;
        }

        /// <inheritdoc />
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(HestonAssetPriceModel);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            var assetPrice = priceModel.GetPriceFactor() as IAssetPrice;
            if (assetPrice == null)
            {
                errorList.Add(ErrorLevel.Error, string.Format("{0} cannot calibrate model of {1} because no price factor is attached.", GetType().Name, priceModel.GetKey()));
                return;
            }

            var hestonModel = (HestonAssetPriceModel)priceModel;
            TryCalculate(calibrationData.PriceFactors, calibrationData, assetPrice, hestonModel, hestonModel, errorList);
            GetUpdatedCurves(priceModel, output);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var graphData = new ResultSeriesList();

            List<ObservedImpliedVol> impliedVols = fCalibrationParameters.ObservedImpliedVols;

            var marketCurves = new Dictionary<string, Curve>();
            var modelCurves = new Dictionary<string, Curve>();

            using (var cache = Vector.Cache(fCalibrationParameters.NumberOfScenarios))
            {
                Vector spotModel = cache.Get();
                Vector strikeModel = cache.Get();
                Vector sign = cache.Get();
                Vector forward = cache.Get();
                Vector modelOptionPrice = cache.Get();
                Vector modelVolatility = cache.Get();

                foreach (var vol in impliedVols)
                {
                    string expiry = (new Period(vol.Time_To_Expiry)).ToString();

                    if (!marketCurves.ContainsKey(expiry))
                    {
                        marketCurves.Add(expiry, new Curve());
                        modelCurves.Add(expiry, new Curve());
                    }

                    double moneyness = vol.Moneyness;

                    spotModel.Assign(vol.Moneyness);
                    strikeModel.Assign(1.0);

                    // Calculate the option price under Heston for the set of parameters
                    fCalibrationParameters.HestonModel.CalculateUndiscountedDomesticEuropeanOptionPrice(modelOptionPrice, forward, 0.0, spotModel, strikeModel, vol.Time_To_Expiry, sign, false);
                    PricingFunctions.ComputeVolFromUndiscountedEuropeanOptionPrice(modelVolatility, modelOptionPrice, forward, strikeModel, vol.Time_To_Expiry, sign);

                    marketCurves[expiry][moneyness] = vol.Observed_Implied_Vol * Percentage.OverPercentagePoint;
                    modelCurves[expiry][moneyness] = modelVolatility[0] * Percentage.OverPercentagePoint;
                }
            }

            graphData.fXAxis.fLabel = "Moneyness";
            graphData.fYAxis.fLabel = "Value (%)";

            // Add a market implied smile and a model implied smile for all time to expiry used in the calibration routine
            foreach (KeyValuePair<string, Curve> item in marketCurves)
            {
                string label = item.Key;
                graphData.Add(new ResultSeries(string.Format("{0} {1} Market Volatility", priceModel.fID, label), LineStyle.Solid, marketCurves[label]));
                graphData.Add(new ResultSeries(string.Format("{0} {1} Model Implied Volatility", priceModel.fID, label), LineStyle.Solid, modelCurves[label]));
            }

            graphData.ToXml(output);
            CalibrationHelper.WriteModelParametersToXml(priceModel, output);
        }

        /// <inheritdoc />
        protected override void SetDrift(ICalibrationData calibrationData, HestonAssetPriceModel model, ErrorList errorList)
        {
            // Set the risk neutral flag of the model if requested
            if (Set_Risk_Neutral_Drift == YesNo.Yes)
                model.Risk_Neutral_Drift = Risk_Neutral_Drift;

            // Calulate drift when not setting Risk_Neutral_Drift, or setting Risk_Neutral_Drift to No
            if (Set_Risk_Neutral_Drift == YesNo.No || Risk_Neutral_Drift == YesNo.No)
            {
                fDriftCalibrator.CalibrateDrift(calibrationData, model, errorList);
            }
        }
    }

    /// <summary>
    /// Concrete implementation of <see cref="StatisticsBasedCalibration"/> to expose the Drift calibration.
    /// </summary>
    /// <remarks>
    /// This is required as the <see cref="HestonImpliedVolatilityCalibration"/> inherits from <see cref="HestonImpliedVolatilityCalibrationBase{TModelParameters}"/>,
    /// and not from <see cref="StatisticsBasedCalibration"/>.
    /// </remarks>
    public class HestonDriftCalibrator : StatisticsBasedCalibration
    {
        /// <summary>
        /// Calibrates the drift for the Heston model.
        /// </summary>
        public void CalibrateDrift(ICalibrationData calibrationData, HestonAssetPriceModel model, ErrorList errorList)
        {
            model.Drift.Clear();
            double drift;
            if (TryGetDriftFromLogStatistics(calibrationData, model, errorList, true, out drift))
                model.Drift[1.0] = drift;
        }
    }
}
