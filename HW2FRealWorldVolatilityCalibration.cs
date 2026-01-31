using System;
using System.ComponentModel;
using System.Drawing;
using System.Globalization;
using System.Xml;

using SunGard.Adaptiv.Analytics.Algorithms;
using SunGard.Adaptiv.Analytics.Algorithms.Optimisation;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibrates a Hull-White 2-Factor model by performing a least squares fit on volatilities.
    /// See F.C. Park (2004), "Implementing Interest Rate Models: A Practical Guide".
    /// </summary>
    /// <remarks>
    /// Since the objective function has multiple local minima, we need to use
    /// global optimization methods.
    /// </remarks>
    [DisplayName("Hull-White 2 Factor Historical Volatility Calibration")]
    public sealed class HW2FHistoricalVolCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        private const int MinimumNumberOfTenors = 6;

        private SqrtPositiveEpsilon fScaleFactor;
        private Vector fTenors;
        private Vector fObservedVols;
        private Vector fTenorWeights;
        private int fNumObjectiveEvaluations;
        private int fNumGradientEvaluations;

        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public HW2FHistoricalVolCalibration()
        {
            Multi_Start_Options            = new HW2FVolMultiStartCalibration.Options();
            Local_Optimizer_Options        = new DisplayableQuasiNewtonOptions();
            Variables                      = new HW2FVariables();
            Tenor_Weights                  = new TenorWeights();
            Differential_Evolution_Options = new HW2FVolDifferentialEvolutionCalibration.Options();
        }

        /// <summary>
        /// Enumeration of optimisation methods for the <see cref="HW2FHistoricalVolCalibration"/>.
        /// </summary>
        public enum HW2FVolOptimisationMethod
        {
            Multi_Start,
            Stepped_Rho,
            Differential_Evolution
        }

        /// <summary>
        /// The variable bounds and initial values.
        /// </summary>
        public HW2FVariables Variables
        {
            get;
            set;
        }

        /// <summary>
        /// The optimization method to use.
        /// </summary>
        public HW2FVolOptimisationMethod Method
        {
            get;
            set;
        }

        /// <summary>
        /// Options for the local optimiser.
        /// </summary>
        public DisplayableQuasiNewtonOptions Local_Optimizer_Options
        {
            get;
            set;
        }

        /// <summary>
        /// Options for the multi-start optimization method.
        /// </summary>
        public HW2FVolMultiStartCalibration.Options Multi_Start_Options
        {
            get;
            set;
        }

        /// <summary>
        /// Options for the differential evolution method.
        /// </summary>
        public HW2FVolDifferentialEvolutionCalibration.Options Differential_Evolution_Options
        {
            get;
            set;
        }

        /// <summary>
        /// List of weights to apply to the tenors.
        /// </summary>
        public TenorWeights Tenor_Weights
        {
            get;
            set;
        }

        /// <inheritdoc />
        public Type ModelType()
        {
            return typeof(HullWhite2FactorInterestRateModel);
        }

        /// <inheritdoc />
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            var hw2fModel = (HullWhite2FactorInterestRateModel)priceModel;
            Vector theoreticalVolVector = Vector.Create(fTenors.Count);
            var hwParams = new HW2FVolParameters(hw2fModel.Alpha_1, hw2fModel.Alpha_2, hw2fModel.Sigma_1.Y[0], hw2fModel.Sigma_2.Y[0], hw2fModel.Correlation);
            HW2FVolCalibrationUtilities.GetPredictedVols(fTenors, hwParams, theoreticalVolVector);

            var observedVolCurve = new Curve();
            var theoreticalVolCurve = new Curve();
            for (int tenorIndex = 0; tenorIndex < fTenors.Count; tenorIndex++)
            {
                observedVolCurve.AppendPoint((float)fTenors[tenorIndex], (float)fObservedVols[tenorIndex]);
                theoreticalVolCurve.AppendPoint((float)fTenors[tenorIndex], (float)theoreticalVolVector[tenorIndex]);
            }

            var graphs = new ResultSeriesList
            {
                new ResultSeries("Historical vols", LineStyle.Dashed, 4, Color.Black, observedVolCurve),
                new ResultSeries("Fitted vols", LineStyle.Solid, 2, Color.Red, theoreticalVolCurve)
            };

            graphs.fXAxis.fLabel = "Tenor (years)";
            graphs.fYAxis.fLabel = "Volatility";
            graphs.ToXml(output);
        }

        /// <inheritdoc />
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            fNumObjectiveEvaluations = 0;
            fNumGradientEvaluations = 0;

            var hw2fModel = CalibrationHelper.Validate<HullWhite2FactorInterestRateModel>(calibrationData, priceModel, output);

            BoundedOptVariable[] boundedVariables;
            if (!Variables.TryGetBoundedVariables(errorList, out boundedVariables))
            {
                // The variables configured by the user are not valid.
                return;
            }

            DoubleCurve observedVolCurve;
            if (!TryGetObservedVols(calibrationData, priceModel, errorList, out observedVolCurve))
            {
                // We were unable to obtain any observed vols.
                return;
            }

            // Store the tenors and observed vols as vectors.
            fTenors       = Vector.Create(observedVolCurve.Count);
            fObservedVols = Vector.Create(observedVolCurve.Count);
            fTenorWeights = Vector.Create(observedVolCurve.Count);
            fTenors.AssignElements(i => observedVolCurve.X[i]);
            fObservedVols.AssignElements(i => observedVolCurve.Y[i]);
            if (!TryPopulateTenorWeights(fTenorWeights))
            {
                errorList.Add(ErrorLevel.Error, string.Format("If tenor weights are supplied then the number supplied '{0}' should match the number of observed tenors {1}.",
                                                                    Tenor_Weights.Count, fTenors.Count));
                return;
            }

            // Use the observed vols to automatically scale the objective function.
            fScaleFactor = 1.0 / Math.Max(fObservedVols.Mean(), CalcUtils.SMALL);

            double[] finalPoint;
            double finalObjectiveValue;
            switch (Method)
            {
                case HW2FVolOptimisationMethod.Multi_Start:
                    {
                        var multistartCalibration = new HW2FVolMultiStartCalibration(GetObjectiveValue, GetGradient, boundedVariables, GetLocalOptimiser(), Multi_Start_Options);
                        var optimiserResult = multistartCalibration.Optimise();
                        finalPoint = optimiserResult.FinalPoint;
                        finalObjectiveValue = optimiserResult.FinalObjective;

                        if (!optimiserResult.DidConverge)
                            errorList.Add(ErrorLevel.Warning, string.Format("Calibration of {0} did not converge.", priceModel.Description()));
                        
                        break;
                    }

                case HW2FVolOptimisationMethod.Stepped_Rho:
                    {
                        var steppedRhoCalibration = new HW2FVolSteppedRhoCalibration(GetObjectiveValue, GetGradient, boundedVariables, GetLocalOptimiser());
                        var optimiserResult = steppedRhoCalibration.Optimise();
                        finalPoint = optimiserResult.FinalPoint;
                        finalObjectiveValue = optimiserResult.FinalObjective;

                        if (!optimiserResult.DidConverge)
                            errorList.Add(ErrorLevel.Warning, string.Format("Calibration of {0} did not converge.", priceModel.Description()));

                        break;
                    }

                case HW2FVolOptimisationMethod.Differential_Evolution:
                    {
                        var optimiserResult = HW2FVolDifferentialEvolutionCalibration.Optimise(GetObjectiveValue, boundedVariables, Differential_Evolution_Options);
                        finalPoint = optimiserResult.FinalPoint;
                        finalObjectiveValue = optimiserResult.FinalObjective;
                        errorList.Add(ErrorLevel.Info, optimiserResult.Message);
                        break;
                    }

                default:
                    throw new InvalidOperationException("Unrecognized HW2FVol calibration method.");
            }

            errorList.Add(ErrorLevel.Info, string.Format(CultureInfo.InvariantCulture, "Final objective value: {0}", finalObjectiveValue));
            errorList.Add(ErrorLevel.Info, string.Format("Total num objective evaluations: {0}", fNumObjectiveEvaluations));
            errorList.Add(ErrorLevel.Info, string.Format("Total num gradient evaluations: {0}", fNumGradientEvaluations));

            var finalParameters = HW2FVolParameters.FromArray(finalPoint);
            finalParameters = HW2FVolCalibrationUtilities.OrderFactors(finalParameters);
            HW2FVolCalibrationUtilities.UpdateModel(finalParameters, hw2fModel);
            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Try to populate the tenor weights. Return success/failure code.
        /// </summary>
        private bool TryPopulateTenorWeights(Vector tenorWeights)
        {
            if (Tenor_Weights.Count == 0)
            {
                tenorWeights.Assign(1.0);
                return true;
            }
            else
            {
                if (Tenor_Weights.Count != fTenors.Count)
                {
                    return false;
                }
                else
                {
                    tenorWeights.AssignElements(i => Tenor_Weights[i].Weight);
                    return true;
                }
            }
        }

        /// <summary>
        /// Returns a newly constructed local optimiser.
        /// </summary>
        private BFGSOptimiser GetLocalOptimiser()
        {
            var localOptimiser = new BFGSOptimiser
            {
                Logger = NullOptimisationLogger.Singleton,
                Options = new QuasiNewtonOptions
                {
                    MaxMajorIterations = Local_Optimizer_Options.GetMaxMajorIterations(),
                    GradientTolerance = Local_Optimizer_Options.GetGradientTolerance(),
                    LineSearcher = QuasiNewtonOptions.LineSearcherType.Armijo
                }
            };

            return localOptimiser;
        }

        /// <summary>
        /// Try to get the curve of observed vols (either by reading it directly from the statistics or calculating it from the archive).
        /// </summary>
        private bool TryGetObservedVols(ICalibrationData calibrationData, PriceModel priceModel, ErrorList errorList, out DoubleCurve observedVolCurve)
        {
            // Reads historical vols from RiskMetrics into a Curve.
            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, priceModel, ReturnMethod.Diff, errorList,
                                               ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
            {
                observedVolCurve = null;
                errorList.Add(ErrorLevel.Error, string.Format("Unable to obtain vol statistics for '{0}'.", priceModel.GetPriceFactor()));
                return false;
            }

            var historicalVolCurve = new DoubleCurve();
            for (int index = 0; index < momentStatistics.Length; ++index)
            {
                if (!isValidStatistics[index])
                    continue;

                double term = momentStatistics[index].ParsePointIDToDouble();
                historicalVolCurve[term] = momentStatistics[index].Vol;
            }

            if (historicalVolCurve.Count == 0)
            {
                observedVolCurve = null;
                errorList.Add(ErrorLevel.Error, string.Format("No valid volatilities found for '{0}'.", priceModel.GetPriceFactor()));
                return false;
            }
            else
            {
                if (historicalVolCurve.Count < MinimumNumberOfTenors)
                    errorList.Add(ErrorLevel.Warning, string.Format("Should have observed volatilities for at least 6 tenors. Actual number of tenors: '{0}'.", historicalVolCurve.Count));

                observedVolCurve = historicalVolCurve;
                return true;
            }
        }

        /// <summary>
        /// Calculates the objective value for the supplied coordinate.
        /// </summary>
        private double GetObjectiveValue(double[] x)
        {
            fNumObjectiveEvaluations++;
            HW2FVolParameters hwParams = HW2FVolParameters.FromArray(x);
            return fScaleFactor * HW2FVolCalibrationUtilities.GetObjectiveValue(hwParams, fTenors, fObservedVols, fTenorWeights);
        }

        /// <summary>
        /// Calculates the gradient for the supplied coordinate.
        /// </summary>
        private void GetGradient(double[] x, double[] vout)
        {
            fNumGradientEvaluations++;
            HW2FVolParameters hwParams = HW2FVolParameters.FromArray(x);
            HW2FVolCalibrationUtilities.GetGradient(fTenors, hwParams, fObservedVols, fTenorWeights, vout);
            for (int i = 0; i < vout.Length; i++)
            {
                vout[i] *= fScaleFactor;
            }
        }

        /// <summary>
        /// List of tenor weights for the <see cref="HW2FHistoricalVolCalibration"/>.
        /// </summary>
        public sealed class TenorWeights : DisplayableList<TenorWeight>
        {
        }

        /// <summary>
        /// Small class to store and validate tenor weights.
        /// </summary>
        public sealed class TenorWeight : NestedPresentableObject
        {
            private double fWeight = 1.0;

            /// <summary>
            /// The weight to apply in the objective function for this particular tenor.
            /// </summary>
            public double Weight
            {
                get
                {
                    return fWeight;
                }
                set
                {
                    if (value < 0.0)
                        throw new ArgumentOutOfRangeException("value", value, "The tenor weight cannot be negative.");

                    fWeight = value;
                }
            }
        }
    }
}
