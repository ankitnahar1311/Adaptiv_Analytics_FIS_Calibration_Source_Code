using System;
using System.Diagnostics;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Collection of utility methods for <see cref="HW2FHistoricalVolCalibration"/>.
    /// </summary>
    public static class HW2FVolCalibrationUtilities
    {
        /// <summary>
        /// Update the supplied 2-factor Hull-White model with the calculated parameters.
        /// </summary>
        public static void UpdateModel(HW2FVolParameters hwParams, HullWhite2FactorInterestRateModel hw2fModel)
        {
            Debug.Assert(hw2fModel != null, "hw2fModel should not be null.");

            // (0) Set non-relevant properties to defaults.
            hw2fModel.Quanto_FX_Correlation_1 = 0.0;
            hw2fModel.Quanto_FX_Correlation_2 = 0.0;
            hw2fModel.Quanto_FX_Volatility.Clear();

            // (1) Set the mean reversion.
            hw2fModel.Alpha_1 = hwParams.Alpha_1.Value;
            hw2fModel.Alpha_2 = hwParams.Alpha_2.Value;

            // (2) Set the volatility.
            Action<DoubleCurve, SafeVolatility> updateSigma = delegate(DoubleCurve curve, SafeVolatility sigma)
            {
                curve.Clear();
                curve.AppendPoint(0.0, sigma.Value);
            };

            updateSigma(hw2fModel.Sigma_1, hwParams.Sigma_1);
            updateSigma(hw2fModel.Sigma_2, hwParams.Sigma_2);

            // (3) Set the correlation.
            hw2fModel.Correlation = hwParams.Rho.Value;
        }

        /// <summary>
        /// Reorder the factors (if necessary) so that the first factor has the larger volatility.
        /// </summary>
        public static HW2FVolParameters OrderFactors(HW2FVolParameters hwParams)
        {
            if (hwParams.Sigma_1.Value >= hwParams.Sigma_2.Value)
                return hwParams;
            else
                return new HW2FVolParameters(hwParams.Alpha_2, hwParams.Alpha_1, hwParams.Sigma_2, hwParams.Sigma_1, hwParams.Rho);
        }
        
        /// <summary>
        /// Calculate B.
        /// </summary>
        public static void GetB(Vector tenors, SafeReversionRate alpha, Vector vout)
        {
            Debug.Assert(tenors != null, "The tenors vector should not be null.");
            Debug.Assert(vout != null, "vout should not be null.");

            vout.Assign((1.0 - VectorMath.Exp(-alpha.Value * tenors)) / alpha.Value);
        }

        /// <summary>
        /// Calculates the variance of a sum of 2 random variables.
        /// </summary>
        public static void GetVarianceOf2RandomVariables(Vector sigma1, Vector sigma2, SafeCorrelation rho, Vector vout)
        {
            Debug.Assert(sigma1 != null, "sigma1 should not be null.");
            Debug.Assert(sigma2 != null, "sigma2 should not be null.");
            Debug.Assert(vout != null, "vout should not be null.");

            vout.Assign(VectorMath.Sqr(sigma1) + VectorMath.Sqr(sigma2) + 2.0 * sigma1 * sigma2 * rho.Value);
        }

        /// <summary>
        /// Calculate the predicted vols from the tenors and HW2F parameters.
        /// </summary>
        public static void GetPredictedVols(Vector tenors, HW2FVolParameters hwParams, Vector vout)
        {
            Debug.Assert(tenors != null, "The tenors vector should not be null.");
            Debug.Assert(vout != null, "vout should not be null.");

            using (var cache = Vector.CacheLike(tenors))
            {
                Vector b1 = cache.Get();
                Vector b2 = cache.Get();
                GetB(tenors, hwParams.Alpha_1, b1);
                b1.MultiplyBy(1.0 / tenors);
                GetB(tenors, hwParams.Alpha_2, b2);
                b2.MultiplyBy(1.0 / tenors);

                Vector scaledSigma1 = cache.Get(hwParams.Sigma_1.Value * b1);
                Vector scaledSigma2 = cache.Get(hwParams.Sigma_2.Value * b2);

                GetVarianceOf2RandomVariables(scaledSigma1, scaledSigma2, hwParams.Rho, vout);
                vout.AssignSqrt(vout);
            }
        }

        /// <summary>
        /// Get the sum-of-squared errors between the observed and predicted volatilities.
        /// </summary>
        public static double GetObjectiveValue(HW2FVolParameters hwParams, Vector tenors, Vector observedVols, Vector weights)
        {
            Debug.Assert(tenors != null, "The tenors vector should not be null.");
            Debug.Assert(observedVols != null, "The observedVols vector should not be null.");

            using (var cache = Vector.Cache(observedVols.Count))
            {
                Vector theoreticalVols = cache.Get();
                GetPredictedVols(tenors, hwParams, theoreticalVols);
                Vector squaredDifferences = cache.Get(theoreticalVols - observedVols);
                squaredDifferences.MultiplyBy(squaredDifferences);
                squaredDifferences.MultiplyBy(weights);
                return squaredDifferences.Sum();
            }
        }

        /// <summary>
        /// Calculate B-dash.
        /// </summary>
        public static void GetBDash(Vector tenors, SafeReversionRate alpha, Vector vout)
        {
            Debug.Assert(tenors != null, "The tenors vector should not be null.");
            Debug.Assert(vout != null, "vout should not be null.");

            using (var cache = Vector.CacheLike(tenors))
            {
                Vector exp = cache.Get(VectorMath.Exp(-alpha.Value * tenors));
                vout.Assign((alpha.Value * tenors * exp - (1.0 - exp)) / (alpha.Value * alpha.Value));
            }
        }

        /// <summary>
        /// Calculate the gradient analytically.
        /// </summary>
        public static void GetGradient(Vector tenors, HW2FVolParameters hwParams, Vector observedVols, Vector weights, double[] vout)
        {
            using (var cache = Vector.CacheLike(tenors))
            {
                double sigma1 = hwParams.Sigma_1.Value;
                double sigma2 = hwParams.Sigma_2.Value;
                double rho = hwParams.Rho.Value;

                Vector vol = cache.Get();
                GetPredictedVols(tenors, hwParams, vol);
                Vector weightedfMinusY = cache.Get(vol - observedVols);
                weightedfMinusY.MultiplyBy(weights);

                Vector b1 = cache.Get();
                Vector derivB1 = cache.Get();
                GetB(tenors, hwParams.Alpha_1, b1);
                b1.Assign(b1 / tenors);
                GetBDash(tenors, hwParams.Alpha_1, derivB1);
                derivB1.Assign(derivB1 / tenors);

                Vector b2 = cache.Get();
                Vector derivB2 = cache.Get();
                GetB(tenors, hwParams.Alpha_2, b2);
                b2.Assign(b2 / tenors);
                GetBDash(tenors, hwParams.Alpha_2, derivB2);
                derivB2.Assign(derivB2 / tenors);

                Vector dVardAlpha1 = cache.Get(2.0 * sigma1 * sigma1 * b1 * derivB1 + 2.0 * rho * sigma1 * sigma2 * b2 * derivB1);
                Vector dVardAlpha2 = cache.Get(2.0 * sigma2 * sigma2 * b2 * derivB2 + 2.0 * rho * sigma1 * sigma2 * b1 * derivB2);
                Vector dVardSigma1 = cache.Get(2.0 * sigma1 * b1 * b1 + 2.0 * rho * sigma2 * b1 * b2);
                Vector dVardSigma2 = cache.Get(2.0 * sigma2 * b2 * b2 + 2.0 * rho * sigma1 * b1 * b2);
                Vector dVardRho = cache.Get(2.0 * sigma1 * sigma2 * b1 * b2);

                Vector dVoldAlpha1 = cache.Get(dVardAlpha1 / (2.0 * vol));
                Vector dVoldAlpha2 = cache.Get(dVardAlpha2 / (2.0 * vol));
                Vector dVoldSigma1 = cache.Get(dVardSigma1 / (2.0 * vol));
                Vector dVoldSigma2 = cache.Get(dVardSigma2 / (2.0 * vol));
                Vector dVoldRho = cache.Get(dVardRho / (2.0 * vol));

                vout[0] = 2.0 * weightedfMinusY.Dot(dVoldAlpha1);
                vout[1] = 2.0 * weightedfMinusY.Dot(dVoldAlpha2);
                vout[2] = 2.0 * weightedfMinusY.Dot(dVoldSigma1);
                vout[3] = 2.0 * weightedfMinusY.Dot(dVoldSigma2);
                vout[4] = 2.0 * weightedfMinusY.Dot(dVoldRho);
            }
        }

        /// <summary>
        /// Wrapper for correlation values to handle validation.
        /// </summary>
        /// <remarks>
        /// Ensures that the value stored is in [<see cref="CalcUtils.MinCorrelation"/>, <see cref="CalcUtils.MaxCorrelation"/>].
        /// </remarks>
        [Serializable]
        public struct SafeCorrelation
        {
            private readonly double fValue;

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public SafeCorrelation(double correlation)
            {
                if (correlation < CalcUtils.MinCorrelation || CalcUtils.MaxCorrelation < correlation)
                    throw new ArgumentException(string.Format("The correlation should lie in the range [{0}, {1}] (inclusive).", CalcUtils.MinVolatility, CalcUtils.MaxVolatility));

                fValue = correlation;
            }

            /// <summary>
            /// Get the underlying value.
            /// </summary>
            public double Value
            {
                get
                {
                    return fValue;
                }
            }
        }

        /// <summary>
        /// Wrapper for reversion rate values to handle validation.
        /// </summary>
        /// <remarks>
        /// Ensures that the value stored is in [<see cref="CalcUtils.MinReversionRate"/>, <see cref="double.PositiveInfinity"/>].
        /// </remarks>
        public struct SafeReversionRate
        {
            private readonly double fValue;

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public SafeReversionRate(double reversionRate)
            {
                if (reversionRate < CalcUtils.MinReversionRate)
                    throw new ArgumentException(string.Format("The mean reversion rate cannot be less than {0}", CalcUtils.MinReversionRate));

                fValue = reversionRate;
            }

            /// <summary>
            /// Get the underlying value.
            /// </summary>
            public double Value
            {
                get
                {
                    return fValue;
                }
            }
        }

        /// <summary>
        /// Wrapper for volatility values to handle validation.
        /// </summary>
        /// <remarks>
        /// Ensures that the value stored is in [<see cref="CalcUtils.MinVolatility"/>, <see cref="CalcUtils.MaxVolatility"/>].
        /// </remarks>
        public struct SafeVolatility
        {
            private readonly double fValue;

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public SafeVolatility(double volatility)
            {
                if (volatility < CalcUtils.MinVolatility || CalcUtils.MaxVolatility < volatility)
                    throw new ArgumentException(string.Format("The volatility should lie in the range [{0}, {1}] (inclusive).", CalcUtils.MinVolatility, CalcUtils.MaxVolatility));

                fValue = volatility;
            }

            /// <summary>
            /// Get the underlying value.
            /// </summary>
            public double Value
            {
                get
                {
                    return fValue;
                }
            }
        }
    }
}
