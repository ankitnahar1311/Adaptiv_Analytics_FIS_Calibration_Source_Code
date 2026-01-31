// ------------------------------------------------------------------------------------
// Owner: Perukrishnen Vytelingum
// Helper class to simulate a model.
// ------------------------------------------------------------------------------------

using System;
using SunGard.Adaptiv.Analytics.Framework;
using SunGard.Adaptiv.Core;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Helper class for anything related to simulating a model.
    /// </summary>
    public static class ModelSimulationHelper
    {
        /// <summary>
        /// Returns correlated random numbers given a list of flags where samples are required.
        /// ValidSample is used to avoid computing correlated random numbers at non-required timestamps. Set ValidSample to number to compute for timestamps.
        /// </summary>
        /// <param name="requiredCorrelatedSample">An array of flag to determine on which date to generate correlated random samples across all time series.</param>
        /// <param name="correlMatrix">Correlation matrix of all time series.</param>
        /// <param name="numTimeSeries">The number of time series.</param>
        /// <param name="lengthTimeSeries">The length of the time series.</param>
        /// <param name="seed">The random seed.</param>
        /// <param name="correlatedSample">Output a two-dimensional array of correlated samples.</param>
        /// <returns>Return true if correlated samples successfully generated.</returns>
        /// <remarks>
        /// For rows where correlated normals are not required, we set the value to NaN. 
        /// This goes for rows with partially missing data where requiredCorrelatedSample is false since we are looking at filling in 
        /// completely missing rows.
        /// </remarks>
        public static bool GetCorrelatedRandomNumbers(bool[] requiredCorrelatedSample, SymmetricMatrix correlMatrix, int numTimeSeries, int lengthTimeSeries, int seed, 
            out double[][] correlatedSample)
        {
            var random = seed == 0 ? new Random() : new Random(seed);

            if (requiredCorrelatedSample.Length < lengthTimeSeries)
                requiredCorrelatedSample = null;

            var uncorrelatedRandomNumbers = new double[numTimeSeries][];
            correlatedSample = new double[numTimeSeries][];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                uncorrelatedRandomNumbers[index] = new double[lengthTimeSeries];
                correlatedSample[index] = new double[lengthTimeSeries];

                for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
                    correlatedSample[index][indexLength] = double.NaN;
            }

            // Generate random numbers
            for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
            {
                if (requiredCorrelatedSample == null || requiredCorrelatedSample[indexLength])
                {
                    for (int indexTimeSeries = 0; indexTimeSeries < numTimeSeries; ++indexTimeSeries)
                    {
                        uncorrelatedRandomNumbers[indexTimeSeries][indexLength] = GetNormalValue(random);
                    }
                }
            }

            if (correlMatrix == null || correlMatrix.Cols == 1)
            {
                correlatedSample = uncorrelatedRandomNumbers;
                return true;
            }

            // Heal matrix if not positive definite. Should never happen in practice.
            Matrix choleskyDecomposition = Matrix.Create(numTimeSeries, numTimeSeries);
            if (!correlMatrix.TryCalculateCholeskyDecomposition(choleskyDecomposition))
            {
                var healedMatrix = CorrelationMatrixHealer.Heal(correlMatrix, CorrelationMatrixHealer.CorrelationHealingMethod.Eigenvalue_Raising, null);
                healedMatrix.TryCalculateCholeskyDecomposition(choleskyDecomposition);
            }

            // Only compute the correlated samples where required
            for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
            {
                if (requiredCorrelatedSample == null || requiredCorrelatedSample[indexLength])
                {
                    for (int index = 0; index < numTimeSeries; ++index)
                    {
                        double sum = 0.0;
                        for (int indexMm = 0; indexMm < numTimeSeries; ++indexMm)
                        {
                            sum += uncorrelatedRandomNumbers[indexMm][indexLength] * choleskyDecomposition[index, indexMm];
                        }

                        correlatedSample[index][indexLength] = sum;
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Returns correlated random numbers given a list of flags where samples are required.
        /// </summary>
        /// <param name="requiredCorrelatedSample">A two-dimensional array of flags where first dimension represents length of timeseries and second dimension represents number of columns in timeseries.
        /// This is used to determine for which points on which dates to generate correlated random samples.</param>
        /// <param name="correlationMatrix">Correlation matrix of all time series.</param>
        /// <param name="seed">The random seed.</param>
        /// <returns>A two-dimensional array of correlated samples.</returns>
        /// <remarks>
        /// For rows where correlated normals are not required, we set the value to NaN. 
        /// This goes for rows with partially missing data where requiredCorrelatedSample is false since we are looking at filling in 
        /// completely missing rows.
        /// </remarks>
        public static double[][] GetCorrelatedRandomNumbers(bool[,] requiredCorrelatedSample, SymmetricMatrix correlationMatrix, int seed)
        {
            Ensure.ArgumentNotNull(requiredCorrelatedSample, "requiredCorrelatedSample");
            Ensure.ArgumentNotNull(correlationMatrix, "correlationMatrix");

            if (requiredCorrelatedSample.Length == 0)
            {
                throw new ArgumentException("The correlated random numbers cannot be generated as required correlated sample flags aren't provided.", "requiredCorrelatedSample");
            }

            int lengthTimeSeries = requiredCorrelatedSample.GetLength(0);
            int numTimeSeries = requiredCorrelatedSample.GetLength(1);

            if (correlationMatrix.Cols != numTimeSeries)
            {
                throw new ArgumentException("Number of columns in the correlation matrix must be same as number of columns in the time series.", "correlationMatrix");
            }

            var random = seed == 0 ? new Random() : new Random(seed);

            var uncorrelatedRandomNumbers = new double[numTimeSeries][];
            var correlatedSample = new double[numTimeSeries][];
            for (int index = 0; index < numTimeSeries; ++index)
            {
                uncorrelatedRandomNumbers[index] = new double[lengthTimeSeries];
                correlatedSample[index] = new double[lengthTimeSeries];

                for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
                    correlatedSample[index][indexLength] = double.NaN;
            }

            // Generate random numbers
            for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
            {
                for (int indexTimeSeries = 0; indexTimeSeries < numTimeSeries; ++indexTimeSeries)
                {
                    if (requiredCorrelatedSample[indexLength, indexTimeSeries])
                    {
                        uncorrelatedRandomNumbers[indexTimeSeries][indexLength] = GetNormalValue(random);
                    }
                }
            }

            // Heal matrix if not positive definite. Should never happen in practice.
            Matrix choleskyDecomposition = Matrix.Create(numTimeSeries, numTimeSeries);
            if (!correlationMatrix.TryCalculateCholeskyDecomposition(choleskyDecomposition))
            {
                var healedMatrix = CorrelationMatrixHealer.Heal(correlationMatrix, CorrelationMatrixHealer.CorrelationHealingMethod.Eigenvalue_Raising, null);
                healedMatrix.TryCalculateCholeskyDecomposition(choleskyDecomposition);
            }

            // Only compute the correlated samples where required
            for (int indexLength = 0; indexLength < lengthTimeSeries; ++indexLength)
            {
                for (int index = 0; index < numTimeSeries; ++index)
                {
                    if (requiredCorrelatedSample[indexLength, index])
                    {
                        double sum = 0.0;
                        for (int indexMm = 0; indexMm < numTimeSeries; ++indexMm)
                        {
                            sum += uncorrelatedRandomNumbers[indexMm][indexLength] * choleskyDecomposition[index, indexMm];
                        }

                        correlatedSample[index][indexLength] = sum;
                    }
                }
            }

            return correlatedSample;
        }

        /// <summary>
        /// Returns a random number drawn from a normal distribution with mean 0 and variance 1.
        /// </summary>
        private static double GetNormalValue(Random randomSeed)
        {
            double u1;
            double s2;

            do
            {
                u1 = 2 * randomSeed.NextDouble() - 1;
                double u2 = 2 * randomSeed.NextDouble() - 1;
                s2 = CalcUtils.Sqr(u1) + CalcUtils.Sqr(u2);
            }
            while ((s2 >= 1) && (Math.Abs(s2 - 0) > CalcUtils.TINY));

            double root = Math.Sqrt(-2 * Math.Log(s2) / s2);
            return u1 * root;
        }
    }
}
