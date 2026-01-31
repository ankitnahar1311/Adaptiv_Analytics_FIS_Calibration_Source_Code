using System;
using System.Collections.Generic;
using System.Linq;

using SunGard.Adaptiv.Analytics.Algorithms;
using SunGard.Adaptiv.Analytics.Algorithms.Optimisation;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Multi-start approach to HW2F Historical Vol Calibration.
    /// </summary>
    /// <remarks>
    /// We begin with a population of starting points in 5-D box.
    /// We run a local optimization for each starting point, and return the result that has the best objective value.
    /// Note: The initial point supplied by the user is included in the population of starting points.
    /// </remarks>
    public class HW2FVolMultiStartCalibration
    {
        private readonly BoundedOptVariable[]       fParentBoundedVariables;
        private readonly Func<double[], double>     fParentObjectiveFunction;
        private readonly Action<double[], double[]> fParentGradientFunction;
        private readonly BFGSOptimiser              fLocalOptimiser;
        private readonly Options                    fOptions;
        private readonly RandomNumberGenerator      fRandom;

        /// <summary>
        /// Initializes a new instance.
        /// </summary>
        public HW2FVolMultiStartCalibration(Func<double[], double> objectiveFunction, Action<double[], double[]> gradientFunction, BoundedOptVariable[] parentBoundedVariables, BFGSOptimiser localOptimiser, Options options)
        {
            fParentObjectiveFunction = objectiveFunction;
            fParentGradientFunction = gradientFunction;

            fParentBoundedVariables  = parentBoundedVariables;
            fLocalOptimiser          = localOptimiser;
            fOptions                 = options;
            
            fRandom = new RandomNumberGenerator();
            fRandom.SetSeed(CalcUtils.ReplaceZeroWithRandom(fOptions.Seed));
        }

        /// <summary>
        /// Run the optimisation and return the result.
        /// </summary>
        public OptimiserResult Optimise()
        {
            var results = new List<Tuple<double, OptimiserResult>>();
            var startingPoints = GetStartingPoints();
            foreach (var startingPoint in startingPoints)
            {
                var problem = new BoxConstrainedProblem(fParentObjectiveFunction, fParentGradientFunction, startingPoint);
                var optimiserResult = fLocalOptimiser.OptimiseProblem(problem);
                results.Add(new Tuple<double, OptimiserResult>(optimiserResult.FinalObjective, optimiserResult));
            }

            var sortedResults = results.OrderBy(result => result.Item1);
            var bestResult = sortedResults.First().Item2;
            int numObjectiveEvaluations = results.Sum(item => item.Item2.NumObjectiveEvaluations);
            int numGradientEvaluations = results.Sum(item => item.Item2.NumGradientEvaluations);

            fLocalOptimiser.Logger.LogError("Total num objective evaluations: {0}", numObjectiveEvaluations);
            fLocalOptimiser.Logger.LogError("Total Num gradient evaluations: {0}", numGradientEvaluations);

            return bestResult;
        }
        
        /// <summary>
        /// Get the collection of random starting points (and the initial point provided by the user).
        /// </summary>
        public List<BoundedOptVariable[]> GetStartingPoints()
        {
            var result = new List<BoundedOptVariable[]>();

            // Include the starting point provided by the user.
            result.Add(fParentBoundedVariables);

            // Add random population.
            for (int i = 0; i < fOptions.Initial_Population_Size; i++)
            {
                result.Add(GetRandomPoint(fParentBoundedVariables));
            }

            return result;
        }

        /// <summary>
        /// Get a random point in the N-D box defined by the variable lower and upper bounds.
        /// </summary>
        private BoundedOptVariable[] GetRandomPoint(IEnumerable<BoundedOptVariable> boundedVariables)
        {
            return boundedVariables.Select(GetUniformVariateInRange).ToArray();
        }

        /// <summary>
        /// Returns a uniform variate between the lower and upper bounds of the variable.
        /// </summary>
        private BoundedOptVariable GetUniformVariateInRange(BoundedOptVariable boundedVariable)
        {
            double randomInitialValue = boundedVariable.Bounds.LowerBound + (boundedVariable.Bounds.UpperBound - boundedVariable.Bounds.LowerBound) * fRandom.NextDouble();
            return new BoundedOptVariable(boundedVariable.Name, boundedVariable.Bounds.LowerBound, randomInitialValue, boundedVariable.Bounds.UpperBound);
        }

        /// <summary>
        /// Collection of options for a multi-start optimization for the HW2F real-world vol calibration.
        /// </summary>
        public class Options : NestedPresentableObject
        {
            private StrictlyPosInt fPopulationSize;

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public Options()
            {
                fPopulationSize = 1000;
                Seed = 0;
            }

            /// <summary>
            /// The number of starting points to generate.
            /// </summary>
            public int Initial_Population_Size
            {
                get
                {
                    return fPopulationSize;
                }
                set
                {
                    fPopulationSize = value;
                }
            }

            /// <summary>
            /// The seed to use in the random number generator.
            /// </summary>
            public int Seed
            {
                get;
                set;
            }
        }
    }
}
