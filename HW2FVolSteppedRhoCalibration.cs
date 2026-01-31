using System;
using System.Collections.Generic;
using System.Linq;

using SunGard.Adaptiv.Analytics.Algorithms.Optimisation;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Stepped-rho approach for the HW2F Historical Vol Calibration.
    /// </summary>
    /// <remarks>
    /// We fix rho at various values and runs an optimization for each value of rho.
    /// We then "fine-tune" our result by running a local optimization over all variables (including rho) starting at the best fixed-rho result.
    /// </remarks>
    public class HW2FVolSteppedRhoCalibration
    {
        private const int NumParentParams = 5;
        private const int NumChildParams = 4;

        private readonly BoundedOptVariable[]       fParentBoundedVariables;
        private readonly Func<double[], double>     fParentObjectiveFunction;
        private readonly Action<double[], double[]> fParentGradientFunction;
        private readonly BFGSOptimiser              fLocalOptimiser;
        private readonly double[]                   fFullXArray = new double[NumParentParams];
        private readonly double[]                   fFullVoutArray = new double[NumParentParams];

        private double fCurrentRho = double.NaN;

        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public HW2FVolSteppedRhoCalibration(Func<double[], double> objectiveFunction, Action<double[], double[]> gradientFunction, BoundedOptVariable[] parentBoundedVariables, BFGSOptimiser localOptimiser)
        {
            fParentObjectiveFunction = objectiveFunction;
            fParentGradientFunction  = gradientFunction;

            fParentBoundedVariables  = parentBoundedVariables;
            fLocalOptimiser          = localOptimiser;
        }

        /// <summary>
        /// Run the optimisation and return the result.
        /// </summary>
        public OptimiserResult Optimise()
        {
            // Run a general optimization from the best fixed-rho result.
            var fineTuningProblem = new BoxConstrainedProblem(fParentObjectiveFunction, fParentGradientFunction, GetNewStartingPoint(GetBestResult()));
            return fLocalOptimiser.OptimiseProblem(fineTuningProblem);
        }

        /// <summary>
        /// Runs an optimization for each fixed value of rho and returns the best result.
        /// </summary>
        private HW2FVolParameters GetBestResult()
        {
            var results = new List<Tuple<double, HW2FVolParameters>>();
            foreach (double rho in GetRhoValues())
            {
                fCurrentRho         = rho;
                var problem         = new BoxConstrainedProblem(GetFixedRhoObjectiveValue, GetFixedRhoGradient, GetChildBoundedVariables());
                var optimizerResult = fLocalOptimiser.OptimiseProblem(problem);

                results.Add(new Tuple<double, HW2FVolParameters>(optimizerResult.FinalObjective, GetParametersFromSolution(optimizerResult.FinalPoint)));
            }

            var sortedResults = results.OrderBy(result => result.Item1);
            return sortedResults.First().Item2;
        }

        /// <summary>
        /// Constructs a new starting point from the supplied <see cref="HW2FVolParameters"/>.
        /// </summary>
        private BoundedOptVariable[] GetNewStartingPoint(HW2FVolParameters hwParams)
        {
            double[] arrayParams = hwParams.ToArray();
            return fParentBoundedVariables.Select((variable, index) => new BoundedOptVariable(variable.Bounds.LowerBound, arrayParams[index], variable.Bounds.UpperBound)).ToArray();
        }

        /// <summary>
        /// Returns the objective value for a fixed value of rho.
        /// </summary>
        /// <remarks>
        /// Note that the full objective function has 5 parameters.
        /// Since we are fixing rho we take 4 parameters here and add rho (which is fixed) separately.
        /// </remarks>
        private double GetFixedRhoObjectiveValue(double[] x)
        {
            Array.Copy(x, fFullXArray, x.Length);
            fFullXArray[4] = fCurrentRho;

            return fParentObjectiveFunction(fFullXArray);
        }

        /// <summary>
        /// Calculates the gradient for a fixed value of rho.
        /// </summary>
        /// <remarks>
        /// Note that the full objective function has 5 parameters.
        /// Since we are fixing rho we take 4 parameters here and add rho (which is fixed) separately.
        /// </remarks>
        private void GetFixedRhoGradient(double[] x, double[] vout)
        {
            Array.Copy(x, fFullXArray, x.Length);
            fFullXArray[4] = fCurrentRho;

            // Get the 5-d gradient:
            fParentGradientFunction(fFullXArray, fFullVoutArray);

            // Read off the first 4 entries for the subproblem:
            Array.Copy(fFullVoutArray, vout, vout.Length);
        }

        /// <summary>
        /// Generates a discrete list of correlation values.
        /// </summary>
        private IEnumerable<double> GetRhoValues()
        {
            var result = new List<double>();
            for (int i = 0; i < 10; i++)
            {
                result.Add(-1.0 + 0.001 * i);
            }

            for (int i = 0; i < 9; i++)
            {
                result.Add(-0.99 + 0.01 * i);
            }

            for (int i = 0; i < 10; i++)
            {
                result.Add(-0.90 + 0.1 * i);
            }

            return result;
        }

        /// <summary>
        /// Convert the supplied array of length 4 into an instance of <see cref="HW2FVolParameters"/>.
        /// </summary>
        private HW2FVolParameters GetParametersFromSolution(double[] x)
        {
            return HW2FVolParameters.FromArray(new double[] { x[0], x[1], x[2], x[3], fCurrentRho });
        }

        /// <summary>
        /// Get all the bounded variables, except rho (which is last).
        /// </summary>
        private BoundedOptVariable[] GetChildBoundedVariables()
        {
            return fParentBoundedVariables.Take(NumChildParams).ToArray();
        }
    }
}
