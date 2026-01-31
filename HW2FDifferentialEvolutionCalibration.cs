using System;

using SunGard.Adaptiv.Analytics.Algorithms;
using SunGard.Adaptiv.Analytics.Algorithms.Optimisation;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Differential Evolution approach to the Hull-White Two-Factor Historical Vol Calibration.
    /// </summary>
    public static class HW2FVolDifferentialEvolutionCalibration
    {
        /// <summary>
        /// Run the optimisation and return the result.
        /// </summary>
        public static DifferentialEvolutionOptimiser.OptimiserResult Optimise(Func<double[], double> objectiveFunction, BoundedOptVariable[] parentBoundedVariables, Options options)
        {
            var problem = new BoxConstrainedProblem(objectiveFunction, parentBoundedVariables);
            return DifferentialEvolutionOptimiser.Optimise(problem, options.GetAlgorithmsOptionsObject());
        }

        /// <summary>
        /// Collection of options for the <see cref="DifferentialEvolutionOptimiser"/>.
        /// </summary>
        public class Options : NestedPresentableObject
        {
            private readonly DifferentialEvolutionOptimiser.Options fOptions = new DifferentialEvolutionOptimiser.Options();

            /// <summary>
            /// Initializes a new instance.
            /// </summary>
            public Options()
            {
                Number_Of_Agents              = 100;
                Mean_Differential_Weight      = 0.5;
                Cross_Over_Probability        = 0.9;
                Max_Generations               = 1000;
                Age_No_Substitution_Stop      = 3;
                Age_No_Min_Change_Stop        = 20;
                Relative_Objective_Gap        = 1E-6;
                Dithering                     = DifferentialEvolutionOptimiser.DitheringMethod.Every_Agent;
                Initial_Population_Generation = DifferentialEvolutionOptimiser.InitialPopulationGenerationMethod.Halton;
                Constraint_Handling           = DifferentialEvolutionOptimiser.ConstraintHandlingMethod.Reflect;
                Seed                          = 1;
            }

            /// <summary>
            /// The size of the population.
            /// </summary>
            public int Number_Of_Agents
            {
                get
                {
                    return fOptions.NumAgents;
                }
                set
                {
                    fOptions.NumAgents = value;
                }
            }

            /// <summary>
            /// Differential weight to use (either with or without dithering).
            /// </summary>
            public double Mean_Differential_Weight
            {
                get
                {
                    return fOptions.DifferentialWeight;
                }
                set
                {
                    fOptions.DifferentialWeight = value;
                }
            }

            /// <summary>
            /// The probability of applying the DE update to a specific coordinate in a point.
            /// </summary>
            public double Cross_Over_Probability
            {
                get
                {
                    return fOptions.CrossOverProbability;
                }
                set
                {
                    fOptions.CrossOverProbability = value;
                }
            }

            /// <summary>
            /// The maximum number of generations.
            /// </summary>
            public int Max_Generations
            {
                get
                {
                    return fOptions.MaxGenerations;
                }
                set
                {
                    fOptions.MaxGenerations = value;
                }
            }

            /// <summary>
            /// Maximum number of generations with no substitutions.
            /// </summary>
            public int Age_No_Substitution_Stop
            {
                get
                {
                    return fOptions.AgeNoSubStop;
                }
                set
                {
                    fOptions.AgeNoSubStop = value;
                }
            }

            /// <summary>
            /// Maximum number of generations with no change in the best objective value.
            /// </summary>
            public int Age_No_Min_Change_Stop
            {
                get
                {
                    return fOptions.AgeNoMinChangeStop;
                }
                set
                {
                    fOptions.AgeNoMinChangeStop = value;
                }
            }

            /// <summary>
            /// Terminate the optimization if the relative difference between the best and worst objective values is less than this threshold.
            /// </summary>
            public double Relative_Objective_Gap
            {
                get
                {
                    return fOptions.RelativeObjectiveGap;
                }
                set
                {
                    fOptions.RelativeObjectiveGap = value;
                }
            }

            /// <summary>
            /// The method to use for handling the box constraints.
            /// </summary>
            public DifferentialEvolutionOptimiser.ConstraintHandlingMethod Constraint_Handling
            {
                get
                {
                    return fOptions.ConstraintHandling;
                }
                set
                {
                    fOptions.ConstraintHandling = value;
                }
            }

            /// <summary>
            /// The dithering method for the differential weight.
            /// </summary>
            public DifferentialEvolutionOptimiser.DitheringMethod Dithering
            {
                get
                {
                    return fOptions.Dithering;
                }
                set
                {
                    fOptions.Dithering = value;
                }
            }

            /// <summary>
            /// The method for generating the population of initial points.
            /// </summary>
            public DifferentialEvolutionOptimiser.InitialPopulationGenerationMethod Initial_Population_Generation
            {
                get
                {
                    return fOptions.InitialPopulationGeneration;
                }
                set
                {
                    fOptions.InitialPopulationGeneration = value;
                }
            }

            /// <summary>
            /// The random seed for generating the population of initial points.
            /// </summary>
            public int Seed
            {
                get
                {
                    return fOptions.Seed;
                }
                set
                {
                    fOptions.Seed = value;
                }
            }

            /// <summary>
            /// Get the actual algorithms options object.
            /// </summary>
            public DifferentialEvolutionOptimiser.Options GetAlgorithmsOptionsObject()
            {
                return fOptions;
            }
        }
    }
}
