using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    public interface IMarketModelParameters
    {
        Period Last_Tenor
        {
            get; set;
        }

        LiborRateList Rates
        {
            get; set;
        }

        CurveList Betas
        {
            get; set;
        }
    }

    public abstract class MarketModelCalibrationBase
    {
        /// <summary>
        /// Set up and run a downhill simplex minimization to fit the given parameter vector using
        /// the given distance function.
        /// </summary>
        /// <param name="x">The vector of parameters to fit</param>
        /// <param name="distanceFunction">The goodness of fit function</param>
        /// <param name="tolerance">Tolerance to be achieved in the function value</param>
        /// <param name="parameters">Auxiliary data needed by the distance function</param>
        /// <param name="distance">The distance of the returned fit as reported by the supplied function</param>
        /// <param name="iterations">Maximum number of iterations to achieve convergence.</param>
        /// <returns>False on error</returns>
        public bool DownhillSimplex<T>(double[] x, FunctionDoubleArrayToDouble<T> distanceFunction, int iterations, double tolerance, T parameters, out double distance)
        {
            // restart the minimization a few times if there is no convergence
            const int NumRestarts = 4;

            distance = 0;
            double[,] simplex = CreateSimplex(x.Length);

            bool result = false;
            for (int i = 0; !result && (i < NumRestarts); i++)
            {
                PopulateSimplex(simplex, x);
                int valuationCount; // not used currently
                double[] xMin;
                bool success = Minimizer<T>.MinimizeDownhillSimplex(distanceFunction, parameters, simplex, (i + 1) * iterations, tolerance, out xMin, out distance, out valuationCount);
                if (success)
                {
                    result = true;
                }

                for (int j = 0; j < xMin.Length; j++)
                    x[j] = xMin[j];
            }

            return result;
        }

        /// <summary>
        /// Create a simplex for N parameters
        /// </summary>
        /// <param name="N">The number of parameters</param>
        /// <returns>The simplex</returns>
        protected static double[,] CreateSimplex(int N)
        {
            // The simplex is an N-dimensional triangle in N-space defined by N+1 vertexes
            return new double[N + 1, N];
        }

        /// <summary>
        /// Populate the given simplex with one vertex at the given parameter vector
        /// </summary>
        /// <param name="simplex">The simplex to populate</param>
        /// <param name="x">The initial guesses of the parameter values</param>
        protected static void PopulateSimplex(double[,] simplex, double[] x)
        {
            const double scale = 0.01;

            // Populate each vertex
            for (int v = 0; v < simplex.GetLength(0); v++)
            {
                // Populate each dimension of the current vertex
                for (int d = 0; d < simplex.GetLength(1); d++)
                {
                    simplex[v, d] = x[d];
                }

                // All vertexes except one are offset from x by scale in one dimension
                if (v < simplex.GetLength(1))
                    simplex[v, v] += scale;
            }
        }
    }

    /// <summary>
    /// Calibrate the Market Interest Rate Model to cap volatity quotations. The general
    /// plan is to fit the instantaneous volatility function as closely as possible to
    /// the following functional form in time t and rate expiry T:
    /// 
    ///         vol(t,T) = f(T-t)g(t)K_T
    ///         
    /// It is economically desirable for the volatility function to be as time-homogeneous
    /// as possible, with as little correction from the rate-specific constant K as possible.
    /// For this reason, f and g are fit independently, with f fit before g. The adjustment
    /// factor K is left until the end, and if the earlier fit is good will be close to 1.
    /// 
    /// This approach is based on "Modern Pricing of Interest-Rate Derivatives" by Riccardo
    /// Rebonato.
    /// </summary>
    [DisplayName("Market Model Cap Volatility Calibration")]
    public class MarketModelCapVolCalibration : MarketModelCapVolCalibrationGeneric<MarketInterestRateModel>
    {
    }

    public abstract class MarketModelCapVolCalibrationGeneric<M> : MarketModelCalibrationBase, IModelCalibration
        where M : PriceModel, IMarketModelParameters
    {
        [NonSerialized]
        protected InputParameters fInputParameters = new InputParameters();
        [NonSerialized]
        protected HomogeneousParameters fHomogeneousParameters = new HomogeneousParameters();
        [NonSerialized]
        protected TemporalParameters fTemporalParameters = new TemporalParameters();

        protected MarketModelCapVolCalibrationGeneric()
        {
            Target_Tenor = 0.25;
            Max_Horizon = 15;
            Max_Iterations = 200;
            Tolerance = 0.0001;
            Random_Seed = 0;
            Temporal_Functions = 3;
            Discretization = 1.0 / 12.0;

            fHomogeneousParameters.Input = fInputParameters;
            fTemporalParameters.Input = fInputParameters;
            fTemporalParameters.Homogeneous = fHomogeneousParameters;
        }

        public Period Target_Tenor
        {
            get; set;
        }

        public Period Max_Horizon
        {
            get; set;
        }

        public int Max_Iterations
        {
            get; set;
        }

        public double Tolerance
        {
            get; set;
        }

        public int Random_Seed
        {
            get; set;
        }

        public int Temporal_Functions
        {
            get; set;
        }

        public Period Discretization
        {
            get; set;
        }

        /// <summary>
        /// Validates input parameters and returns the priceModel as a M model.
        /// </summary>
        public virtual M Validate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output)
        {
            M model = CalibrationHelper.Validate<M>(calibrationData, priceModel, output);

            // Do not allow calibration of spread curves
            if (model.fID.Count > 1)
                throw new CalibrationException(string.Format("{0} cannot calibrate model of {1}.", GetType().Name, CalibrationHelper.SafeSubType(priceModel)));

            // Calibration only suitable for scripting models and Interest Rate price factors.
            IPriceFactor priceFactor = model.GetPriceFactor();
            if (model.GetType() == typeof(MarketInterestRateModel) && priceFactor != null && priceFactor.GetType() != typeof(InterestRate))
                throw new CalibrationException(string.Format("{0} cannot calibrate model of {1}.", GetType().Name, priceFactor.GetType().Name));

            return model;
        }

        /// <summary>
        /// Calibrates the Market model using the current caplet vols.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            M model = Validate(calibrationData, priceModel, output);

            // Set up the input parameters (target of the calibration)
            fInputParameters.Initialize(this, priceModel, calibrationData.MarketData, calibrationData.PriceFactors);

            // Fit the time-homogeneous volatility parameters
            fHomogeneousParameters.Calibrate(this);

            // Fit the time-varying volatility parameters
            fTemporalParameters.Calibrate(this);

            // Populate the model
            model.Last_Tenor = Target_Tenor;
            model.Rates.Clear();
            for (int i = 0; i < fInputParameters.Expiries.Count; i++)
            {
                double expiry = fInputParameters.Expiries[i];
                LiborRate rate = new LiborRate();
                model.Rates.Add(rate);
                rate.Volatility = new RateCurve();
                rate.Volatility[0d] = fTemporalParameters.Volatility(0, expiry);
                for (int tIdx = 0; tIdx <= fInputParameters.Discretization.Count; tIdx++)
                {
                    double t = Math.Min(expiry, fInputParameters.Discretization[tIdx]);
                    rate.Volatility[t] = fTemporalParameters.Volatility(t, expiry);

                    if (t >= expiry)
                        break;
                }

                // Adjust the variance of each expiry to exactly match the market price
                double K = Math.Sqrt(fInputParameters.BlacksVariance(expiry) / fTemporalParameters.Variance(0, expiry, i));
                rate.Volatility.MultiplyBy(K);
            }

            // 1-factor model has constant beta = 1
            model.Betas = new CurveList();
            RateCurve beta1 = new RateCurve();
            beta1[0] = 1;
            model.Betas.Add(beta1);

            WriteCalibrationDataAsXML(output);
        }

        /// <summary>
        /// Rebuild output curve after manual changes to price model
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for</param>
        /// <param name="output">Revised output</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            WriteCalibrationDataAsXML(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(M);
        }

        /// <summary>
        /// Write calibration data to an XML stream.
        /// </summary>
        /// <param name="model">The calibrated model</param>
        /// <param name="writer">The output XML stream</param>
        protected void WriteCalibrationDataAsXML(XmlWriter writer)
        {
            try
            {
                var marketVolatility = fInputParameters.TargetCurve.Clone();
                marketVolatility.MultiplyBy(Percentage.OverPercentagePoint);

                var homogeneousVolatility = new DoubleCurve();
                var parametricVolatility = new DoubleCurve();
                var instantaneousVolatility = new DoubleCurve();
                for (int i = 0; i < fInputParameters.Expiries.Count; i++)
                {
                    double t = fInputParameters.Expiries[i];
                    double v = fHomogeneousParameters.Variance(0, t, t);
                    homogeneousVolatility[t] = Math.Sqrt(v / t);

                    v = fTemporalParameters.Variance(0, t, i);
                    parametricVolatility[t] = Math.Sqrt(v / t);

                    instantaneousVolatility[t] = fHomogeneousParameters.Volatility(t);
                }

                homogeneousVolatility.MultiplyBy(Percentage.OverPercentagePoint);
                parametricVolatility.MultiplyBy(Percentage.OverPercentagePoint);
                instantaneousVolatility.MultiplyBy(Percentage.OverPercentagePoint);

                ResultSeriesList graphData = new ResultSeriesList();
                graphData.fXAxis.fLabel = "Time (Years)";
                graphData.fYAxis.fLabel = "Volatility or Rate (%)";
                graphData.Add(new ResultSeries("Market Volatility", LineStyle.Solid, marketVolatility));
                graphData.Add(new ResultSeries("Total RMS Parametric Volatility", LineStyle.Solid, parametricVolatility));
                graphData.Add(new ResultSeries("Homogeneous RMS Parametric Volatility", LineStyle.Solid, homogeneousVolatility));
                graphData.Add(new ResultSeries("Instantaneous Parametric Volatility (in T-t)", LineStyle.Solid, instantaneousVolatility));

                graphData.ToXml(writer);
            }
            catch (XmlException e)
            {
                throw new CalibrationException("Error writing calibration data", e);
            }
        }

        /// <summary>
        /// Parameters relating to fitting the time-homogeneous functional form f(T-t).
        /// The form of this function is assumed to be:
        ///               f(t) = [a + b(T-t)] exp[-c(T-t)] + d
        /// </summary>
        public class HomogeneousParameters
        {
            protected double[] fX = new double[4];
            protected bool fRecalculate = true;
            protected RateCurve fVariance = new RateCurve();

            public InputParameters Input
            {
                get; set;
            }

            public double A
            {
                get { return fX[0]; }
                set { fRecalculate = true; fX[0] = value; }
            }

            public double B
            {
                get { return fX[1]; }
                set { fRecalculate = true; fX[1] = value; }
            }

            public double C
            {
                get { return fX[2]; }
                set { fRecalculate = true; fX[2] = value; }
            }

            public double D
            {
                get { return fX[3]; }
                set { fRecalculate = true; fX[3] = value; }
            }

            public double Distance
            {
                get; set;
            }

            /// <summary>
            /// Returns true if the function parameters are valid, false otherwise
            /// </summary>
            public bool IsValid()
            {
                return (C <= 0 || D <= 0 || (A + D) <= 0) ? false : true;
            }

            /// <summary>
            /// Returns the instantaneous volatility as a function of T - t.
            /// </summary>
            /// <param name="dT">Expiry T minus time t</param>
            public double Volatility(double dT)
            {
                return (A + B * dT) * Math.Exp(-C * dT) + D;
            }

            /// <summary>
            /// Returns the cumulative variance as a function of T - t from t=t1 to t=t2.
            /// </summary>
            /// <param name="t1">The start of the integration period</param>
            /// <param name="t2">The end of the integration period</param>
            /// <param name="T">The expiry of the relevant interest rate</param>
            public double Variance(double t1, double t2, double T)
            {
                Debug.Assert(t1 < t2, "t1 must be less than t2");
                Debug.Assert(t2 <= T, "t2 must be less than or equal to expiry");

                if (fRecalculate)
                    RecalculateVariance();
                fRecalculate = false;

                // The variance curve is an increasing function of dT i.e. decreasing in t
                return fVariance[T - t1] - fVariance[T - t2];
            }

            /// <summary>
            /// Fit the time-homogeneous volatility function parameters a, b, c, and d.
            /// </summary>
            /// <param name="calibration">The calling calibration object.</param>
            public void Calibrate(MarketModelCapVolCalibrationGeneric<M> calibration)
            {
                // get close to a solution of the right form using heuristics about the shape
                // of f(T-t)
                A = 0;
                B = 0;
                C = 0;
                D = 0;
                EstimateD();
                EstimateA();
                EstimateBandC();

                double distance = 0;
                if (!calibration.DownhillSimplex(fX, DistanceFunction, calibration.fInputParameters.Iterations, calibration.fInputParameters.Tolerance, this, out distance))
                    throw new CalibrationException("Homogeneous minimization fit failed.");

                Distance = distance;
            }

            /// <summary>
            /// Distance function for fitting time-homogeneous volatility function parameters a,
            /// b, c and d.
            /// </summary>
            /// <param name="x">parameter vector for which to calculate the distance.</param>
            /// <param name="homogeneousParameters">Homogeneous parameters.</param>
            /// <param name="distance">The distance from the ideal fit.</param>
            /// <returns>bool indicating success.</returns>
            /// <exception cref="ArgumentException"></exception>
            protected static bool DistanceFunction(double[] x, HomogeneousParameters homogeneousParameters, out double distance)
            {
                distance = 0;

                homogeneousParameters.A = x[0];
                homogeneousParameters.B = x[1];
                homogeneousParameters.C = x[2];
                homogeneousParameters.D = x[3];
                if (!homogeneousParameters.IsValid())
                    return false;

                distance = homogeneousParameters.CalculateDistance();

                return true;
            }

            /// <summary>
            /// Calculate and store the time-homogeneous variance as a function of T-t.
            /// </summary>
            protected void RecalculateVariance()
            {
                const int NumSteps = 10;

                fVariance.Clear();
                fVariance[0] = 0;
                double variance = 0;
                double oldT = 0;
                for (int i = 0; i < Input.Discretization.Count; i++)
                {
                    double newT = Input.Discretization[i];
                    double dV = 0;
                    double stepSize = (newT - oldT) / NumSteps;
                    for (int j = 0; j < NumSteps; j++)
                    {
                        double x = oldT + (j + 0.5) * stepSize;
                        double fx = Volatility(x);
                        dV += fx * fx;
                    }
                    dV *= stepSize;
                    variance += dV;
                    fVariance[newT] = variance;
                    oldT = newT;
                }
                fRecalculate = false;
            }

            /// <summary>
            /// Estimate an initial value for A from the target Black's volatilty using the
            /// fact that instantaneous and average volatilities converge as dT -> zero.
            /// </summary>
            protected void EstimateA()
            {
                double dT = Input.TargetCurve.X[0];
                A = (Input.TargetCurve.Y[0] - D) * Math.Exp(C * dT) - B * dT;
            }

            /// <summary>
            /// Estimate initial values for B and C using the facts that the first derivative
            /// of the volatility function at the origin is m = b-ca and the location of the
            /// extremum of the volatility function is k = (b-ca)/cb. These two equations with two
            /// two unknowns yield a quadratic function for c by substituting b:
            ///                   a c^2 + m c - m / k = 0
            /// </summary>
            protected void EstimateBandC()
            {
                DoubleCurve curve = Input.TargetCurve;
                double slope = curve.Gradient(0);
                for (int i = 1; slope == 0 && i < curve.Count - 1; i++)
                    slope = curve.Gradient(i);

                if (Math.Abs(slope) < 1.0E-6)
                {
                    B = 0;
                    C = 1.0E8;
                }
                else
                {
                    double xMax = 0;
                    double yMax = curve[0];
                    for (int i = 1; i < curve.Count && yMax <= curve.Y[i]; i++)
                    {
                        xMax = curve.X[i];
                        yMax = curve.Y[i];
                    }

                    if (xMax > 0)
                    {
                        double a = A;
                        double b = slope;
                        double c = -slope / xMax;

                        if (a == 0)
                        {
                            C = c / -b;
                        }
                        else
                        {
                            double sign_b = (b >= 0) ? +1 : -1;
                            double square = b * b - 4 * a * c;
                            if (square > 0)
                            {
                                double q = -0.5 * (b + sign_b * Math.Sqrt(square));
                                double x1 = q / a;
                                double x2 = c / q;

                                // C must be positive; if both roots are positive, choose the smallest
                                if (x1 <= 0)
                                    C = x2;
                                else if (x2 <= 0)
                                    C = x1;
                                else
                                    C = x1 < x2 ? x1 : x2;
                            }
                            else
                                C = c / -b;
                        }
                    }
                    else
                    {
                        C = -slope / A;
                    }

                    // B is determined by slope at origin = b - ca
                    B = slope + A * C;
                }
            }

            /// <summary>
            /// Estimate an initial value for D using the fact that the long-run volatility
            /// converges to D.
            /// </summary>
            protected void EstimateD()
            {
                double dT = Input.TargetCurve.LastX();
                D = Input.TargetCurve.LastY() - (A + B * dT) * Math.Exp(-C * dT);
            }

            /// <summary>
            /// Returns the sum of the squared differences between the target volatility 
            /// and the integrated homogeneous volatility at each rate expiry
            /// </summary>
            protected double CalculateDistance()
            {
                double distance = 0;
                foreach (float expiry in Input.Expiries)
                {
                    double difference = Input.TargetCurve[expiry] - Math.Sqrt(Variance(0, expiry, expiry) / expiry);
                    distance += difference * difference;
                }
                return distance;
            }
        }

        /// <summary>
        /// Parameters relating to fitting the time-varying functional form g(t).
        /// The form of this function is assumed to be:
        /// 
        ///            g(t) = 1 + sum_1^N e_i sin[(t pi i)/maturity + k_i]exp(-a t)
        ///  
        /// For some number of basis functions N. 
        /// 
        /// The parameters of g(t) are stored in an array of doubles in the order [e_1, k_1, e_2,
        /// k_2, ... , e_N, k_N, a].
        /// </summary>
        public class TemporalParameters
        {
            protected RateCurve[] fVariance = null;

            public TemporalParameters()
            {
                Basis = new BasisFunctions();
            }

            public InputParameters Input
            {
                get; set;
            }

            public HomogeneousParameters Homogeneous
            {
                get; set;
            }

            public BasisFunctions Basis
            {
                get; set;
            }

            public double Distance
            {
                get; set;
            }

            /// <summary>
            /// Fit the time-varying function parameters.
            /// </summary>
            /// <param name="calibration">The calling calibration object</param>
            public void Calibrate(MarketModelCapVolCalibrationGeneric<M> calibration)
            {
                // allocate the variance cache
                fVariance = new RateCurve[Input.Expiries.Count];
                for (int i = 0; i < fVariance.Length; i++)
                    fVariance[i] = new RateCurve();

                // initialize the basis function parameters
                Basis.Initialize(Math.Max(1, calibration.Temporal_Functions), Input.LastExpiry, Input.Seed);

                // copy the function parameters out
                double[] x = new double[Basis.Length];
                for (int i = 0; i < x.Length; i++)
                    x[i] = Basis[i];

                double distance = 0;
                if (!calibration.DownhillSimplex(x, DistanceFunction, calibration.fInputParameters.Iterations, calibration.fInputParameters.Tolerance, calibration.fTemporalParameters, out distance))
                    throw new CalibrationException("Temporal minimization fit failed.");

                // copy the function parameters back in
                for (int i = 0; i < x.Length; i++)
                    Basis[i] = x[i];

                Distance = distance;
            }

            /// <summary>
            /// Returns the instantaneous volatility at t f(T-t)g(t) for the rate associated
            /// with expiry time T.
            /// </summary>
            /// <param name="t">The time</param>
            /// <param name="expiry">The expiry of the associated rate</param>
            public double Volatility(double t, double expiry)
            {
                Debug.Assert(t <= expiry, "t must be less than or equal to expiry");

                return Homogeneous.Volatility(expiry - t) * Basis.G(t);
            }

            /// <summary>
            /// Returns the cumulative variance as a function of time t and expiry T
            /// from t=t1 to t=t2.
            /// </summary>
            /// <param name="t1">The start of the integration period</param>
            /// <param name="t2">The end of the integration period</param>
            /// <param name="rateIdx">The expiry index of the associated rate</param>
            public double Variance(double t1, double t2, int rateIdx)
            {
                Debug.Assert(t1 < t2, "t1 must be less than t2");
                Debug.Assert(t2 <= Input.Expiries[rateIdx], "t2 must be less than or equal to expiry");

                if (Basis.Recalculate)
                    RecalculateVariance();

                // Each rate-specific variance curve is an increasing function of time
                return fVariance[rateIdx][t2] - fVariance[rateIdx][t1];
            }

            /// <summary>
            /// Distance function for fitting the time-varying volatility function parameter vector E
            /// </summary>
            /// <param name="x">Parameter vector for which to calculate the distance.</param>
            /// <param name="temporalParameters">Temporal parameters.</param>
            /// <param name="distance">The distance from the ideal fit.</param>
            /// <returns>bool indicating success.</returns>
            /// <exception cref="ArgumentException"></exception>
            protected static bool DistanceFunction(double[] x, TemporalParameters temporalParameters, out double distance)
            {
                distance = 0;

                for (int i = 0; i < x.Length; i++)
                    temporalParameters.Basis[i] = x[i];

                // The distance is the sum of the squared integrated volatility differences
                for (int i = 0; i < temporalParameters.Input.Expiries.Count; i++)
                {
                    double expiry = temporalParameters.Input.Expiries[i];
                    double difference = temporalParameters.Input.TargetCurve[expiry] - Math.Sqrt(temporalParameters.Variance(0, expiry, i) / expiry);
                    distance += difference * difference;
                }

                return true;
            }

            protected void RecalculateVariance()
            {
                for (int rateIdx = 0; rateIdx < Input.Expiries.Count; rateIdx++)
                    RecalculateVariance(rateIdx);

                Basis.Recalculate = false;
            }

            protected void RecalculateVariance(int rateIdx)
            {
                const int NumSteps = 10;

                fVariance[rateIdx].Clear();
                fVariance[rateIdx][0] = 0;
                double expiry = Input.Expiries[rateIdx];
                double oldT = 0;
                double variance = 0;
                for (int i = 0; i < Input.Discretization.Count && oldT <= expiry; i++)
                {
                    double newT = Math.Min(Input.Discretization[i], expiry);
                    double dV = 0;
                    double stepSize = (newT - oldT) / NumSteps;
                    for (int j = 0; j < NumSteps; j++)
                    {
                        double x = oldT + (j + 0.5) * stepSize;
                        double fx = Volatility(x, expiry);
                        dV += fx * fx;
                    }
                    dV *= stepSize;
                    variance += dV;
                    fVariance[rateIdx][newT] = variance;
                    oldT = newT;
                }
            }

            /// <summary>
            /// This class manages the internal representation of the temporal function parameters.
            /// </summary>
            public class BasisFunctions
            {
                protected double[] fE; // parameters of g(t)

                public int Dimension
                {
                    get { return (fE.Length - 1) / 2; }
                    set { fE = new double[2 * value + 1]; }
                }

                public double Last
                {
                    get { return fE[fE.Length - 1]; }
                    set { fE[fE.Length - 1] = value; }
                }

                public int Length
                {
                    get { return fE.Length; }
                }

                public double Maturity
                {
                    get; set;
                }

                public bool Recalculate
                {
                    get; set;
                }

                public double this[int i]
                {
                    get { return fE[i]; }
                    set { Recalculate = true; fE[i] = value; }
                }

                public void Initialize(int dimension, double maturity, int seed)
                {
                    RandomNumberGenerator generator = new RandomNumberGenerator();
                    generator.SetSeed(seed);
                    Dimension = dimension;
                    Maturity = maturity;
                    for (int i = 0; i < Dimension; i++)
                    {
                        fE[2 * i] = generator.NextDouble();
                        fE[2 * i + 1] = generator.NextDouble();
                    }
                    Last = 0.1;
                }

                /// <summary>
                /// Calculate the value of g(t).
                /// </summary>
                /// <param name="t">The time at which to calculate g(t)</param>
                public double G(double t)
                {
                    double fT = 0;
                    for (int b = 0; b < Dimension; b++)
                    {
                        fT += fE[2 * b] * Math.Sin(t * Math.PI * (b + 1) / Maturity + fE[2 * b + 1]);
                    }
                    fT *= Math.Exp(-Last * t);
                    return 1 + fT;
                }
            }
        }

        /// <summary>
        /// Parameters defining the target of the calibration and controlling its process.
        /// </summary>
        public class InputParameters
        {
            public InputParameters()
            {
                TargetCurve = new RateCurve();
                Expiries = new FloatArray();
                Discretization = new FloatArray();
            }

            public RateCurve TargetCurve
            {
                get; set;
            }

            public FloatArray Expiries
            {
                get; set;
            }

            public int Iterations
            {
                get; set;
            }

            public double Tolerance
            {
                get; set;
            }

            public int Seed
            {
                get; set;
            }

            public FloatArray Discretization
            {
                get; set;
            }

            public double LastExpiry
            {
                get { return Expiries[Expiries.Count - 1]; }
            }

            /// <summary>
            /// Initialize the input parameters from the calibration, price model, and market data.
            /// </summary>
            /// <param name="calibration">The calling calibration</param>
            /// <param name="priceModel">The price model to be calibrated</param>
            /// <param name="marketData">The market data containing the target volatilities</param>
            public void Initialize(MarketModelCapVolCalibrationGeneric<M> calibration, PriceModel priceModel, MarketDataContainer marketData, PriceFactorList priceFactors)
            {
                TargetCurve.Clear();
                Expiries.Clear();
                Discretization.Clear();

                Iterations = calibration.Max_Iterations;
                Tolerance = calibration.Tolerance;
                Seed = CalcUtils.ReplaceZeroWithRandom(calibration.Random_Seed);

                // construct the rate expiry array and the target volatility curve
                string factorID = priceModel.fID.ToCode();
                IInterestRate df;

                if (priceModel.fID.Count > 1)
                {
                    List<ScriptAssetPriceBase> assets;
                    List<ScriptRateBase> interestRates;
                    string payoffCurrency;
                    string discountCurrencyID;
                    bool populated;
                    ScriptBaseFactor.PopulateScriptFactorElements(priceFactors, out assets, out interestRates, out payoffCurrency, out discountCurrencyID, priceModel.fID, out populated);
                    FactorID interestRateID = new FactorID(interestRates[0].Rate);
                    df = priceFactors.GetInterface<IInterestRate>(interestRateID);
                }
                else
                {
                    df = priceFactors.GetInterface<IInterestRate>(factorID);
                }

                IInterestRateVol vol = marketData.PriceFactors.GetInterface<IInterestRateVol>(df.GetCurrency());
                if (vol.GetSurface().X.Count <= 0)
                    throw new CalibrationException("Interest rate volatility surface is empty.");
                double expiry = calibration.Target_Tenor;
                while (expiry <= calibration.Max_Horizon)
                {
                    Expiries.Add((float)expiry);

                    double volatility = vol.GetUnflooredFactorVol(0.05, 0.05, expiry, calibration.Target_Tenor);
                    volatility = Math.Max(volatility, CalcUtils.MinVolatility); // floor at min
                    volatility = Math.Min(volatility, CalcUtils.MaxVolatility); // cap at max

                    TargetCurve[expiry] = volatility;
                    expiry += calibration.Target_Tenor;
                }

                // construct the discretization array ensuring that it ends at the last expiry
                calibration.Discretization = Math.Min(calibration.Discretization, calibration.Target_Tenor);
                double dT = calibration.Discretization;
                double time = 0;
                do
                {
                    time += dT;
                    time = Math.Min(calibration.Max_Horizon, time);
                    Discretization.Add((float)time);
                } while (time < calibration.Max_Horizon);
            }

            /// <summary>
            /// Return the total variance to time T corresponding to the Black's volatility
            /// quoted by the market.
            /// </summary>
            public double BlacksVariance(double T)
            {
                double vol = TargetCurve[T];
                return vol * vol * T;
            }
        }
    }

    /// <summary>
    /// Calibrate the factor loadings of the Market Interest Rate model to  constant-maturity
    /// rate correlations. The correlations of the fixed-maturity rates in
    /// the model during simulation are determined by their remaining time to maturity - the
    /// rate correlation functions are required to be completely homogeneous.
    /// </summary>
    public abstract class MarketModelCorrelationCalibrationBase : MarketModelCalibrationBase
    {
        protected const double MaxCorr = 1;
        protected const double MinCorr = 0.01;

        protected int fFactors = 0;
        [NonSerialized]
        protected double[] fExpiries = null;
        [NonSerialized]
        protected double[] fTenors = null;
        [NonSerialized]
        protected Matrix fTarget = null;
        [NonSerialized]
        protected Matrix fBetas = null;
        [NonSerialized]
        protected double[] fThetas = null;

        protected MarketModelCorrelationCalibrationBase()
        {
            Factors = 3;
            Max_Iterations = 10000;
            Tolerance = 1E-7;
        }

        public int Factors
        {
            get { return fFactors; }
            set { fFactors = Math.Max(1, value); }
        }

        public int Max_Iterations
        {
            get; set;
        }

        public double Tolerance
        {
            get; set;
        }

        public int Random_Seed
        {
            get; set;
        }

        /// <summary>
        /// Find the best fitting set of "angles" theta to the target correlation matrix and
        /// store the betas (model driver vectors) that they imply.
        /// </summary>
        public void FitThetas()
        {
            if (Factors < 2)
            {
                for (int i = 0; i < fExpiries.Length; i++)
                    fBetas[i, 0] = 1;

                return;
            }

            RandomNumberGenerator generator = new RandomNumberGenerator();
            generator.SetSeed(CalcUtils.ReplaceZeroWithRandom(Random_Seed));
            for (int i = 0; i < fThetas.Length; i++)
                fThetas[i] = 0.5 * Math.PI * generator.NextDouble();

            double distance = 0;
            if (!DownhillSimplex(fThetas, DistanceFunction, Max_Iterations, Tolerance, this, out distance))
            {
                // tolerate failures here instead of throwing an exception - the downhill simplex is too unreliable in this application
                //throw new CalibrationException(typeof(M), "Correlation", "Theta minimization fit failed.");
            }
        }

        /// <summary>
        /// Calculate the current distance for the downhill simplex minimizer. This function
        /// causes changes to the betas.
        /// </summary>
        /// <param name="x">Set of "angles" theta that control the betas</param>
        /// <param name="parameters">The calibration object</param>
        /// <param name="distance">The current distance of the model from the target</param>
        /// <returns></returns>
        protected static bool DistanceFunction(double[] x, MarketModelCorrelationCalibrationBase parameters, out double distance)
        {
            parameters.SetThetas(x);
            distance = parameters.CalculateDistance();

            return true;
        }

        /// <summary>
        /// Write calibration data to an XML stream.
        /// </summary>
        /// <param name="writer">The output XML stream</param>
        protected void WriteCalibrationDataAsXML(XmlWriter writer)
        {
            try
            {
                const float scale = (float)Percentage.OverPercentagePoint;
                Curve targetCorrelation = new Curve();
                Curve modelCorrelation = new Curve();
                //Curve beta1 = new Curve();
                //Curve beta2 = new Curve();
                //Curve beta3 = new Curve();
                int iBase = fExpiries.Length / 2;
                for (int i = 0; i < fExpiries.Length; i++)
                {
                    targetCorrelation[fExpiries[i]] = fTarget[iBase, i];
                    double rho = 0;
                    for (int j = 0; j < Factors; j++)
                        rho += fBetas[i, j] * fBetas[iBase, j];
                    modelCorrelation[fExpiries[i]] = rho;
                    //beta1[fExpiries[i]] = fBetas[i, 0];
                    //beta2[fExpiries[i]] = fBetas[i, 1];
                    //beta3[fExpiries[i]] = fBetas[i, 2];
                }
                targetCorrelation.MultiplyBy(scale);
                modelCorrelation.MultiplyBy(scale);
                //beta1.MultiplyBy(scale);
                //beta2.MultiplyBy(scale);
                //beta3.MultiplyBy(scale);

                ResultSeriesList graphData = new ResultSeriesList();
                graphData.fXAxis.fLabel = "Rate Expiry (Years)";
                graphData.fYAxis.fLabel = "Correlation to " + fExpiries[iBase].ToString(CultureInfo.InvariantCulture) + "Year rate";
                //graphData.YAxis.Label = "Model driver value";
                graphData.Add(new ResultSeries("Target", LineStyle.Solid, targetCorrelation));
                graphData.Add(new ResultSeries("Model", LineStyle.Solid, modelCorrelation));
                //graphData.Add(new ResultSeries("Beta1", LineStyle.Solid, beta1));
                //graphData.Add(new ResultSeries("Beta2", LineStyle.Solid, beta2));
                //graphData.Add(new ResultSeries("Beta3", LineStyle.Solid, beta3));
                graphData.ToXml(writer);
            }
            catch (XmlException e)
            {
                throw new CalibrationException("Error writing calibration data.", e);
            }
        }

        /// <summary>
        /// Initialize internal state from current model state and calibration parameters.
        /// </summary>
        /// <param name="model">The market model to be calibrated</param>
        protected void Initialize(IMarketModelParameters model)
        {
            int nRates = model.Rates.Count;
            if (nRates < 1)
                throw new CalibrationException("There are no rates in the input market model - calibrate the volatility first.");

            Factors = Math.Min(Factors, nRates);
            fTarget = Matrix.Create(nRates, nRates);
            fBetas = Matrix.Create(nRates, Factors);

            fExpiries = new double[nRates];
            fTenors = new double[nRates];
            double t1 = 0;
            double t2 = model.Rates[0].GetReset();
            fExpiries[0] = t2;
            for (int i = 1; i < nRates; i++)
            {
                t2 = model.Rates[i].GetReset();
                fExpiries[i] = t2;
                fTenors[i - 1] = t2 - t1;
                t1 = t2;
            }
            fTenors[nRates - 1] = model.Last_Tenor;

            fThetas = new double[nRates * (Factors - 1)];
        }

        /// <summary>
        /// Build the target correlation matrix.
        /// </summary>
        protected abstract void SetTarget(IStatistics statistics, PriceModel priceModel);

        /// <summary>
        /// Orthogonalize the model driver vectors (betas).
        /// </summary>
        protected void Orthogonalize()
        {
            int nRates = fExpiries.Length;

            // construct the model correlation matrix
            var rho = SymmetricMatrix.Create(nRates);
            for (int i = 0; i < nRates; i++)
            {
                for (int j = i; j < nRates; j++)
                {
                    rho[i, j] = ModelCorrelation(i, j);
                }
            }

            // Orthogonalize the model correlation matrix, and sort by eigenvalue
            var eigenCalculateResult = rho.CalculateEigenVectorsAndValues();

            // Calculate the square root of the eigenvalues, accounting for the possibility
            // that small eigenvalues might be negative after finite-precision calculation
            var flooredRootEigenValues = eigenCalculateResult.EigenValues;
            flooredRootEigenValues.Assign(VectorMath.Sqrt(VectorMath.Max(flooredRootEigenValues, 0.0)));

            // Assign the orthogonal basis to beta after normalizing to account for numerical errors
            var eigenVectors = eigenCalculateResult.EigenVectors;
            for (int rIdx = 0; rIdx < nRates; rIdx++)
            {
                double normalize = 0;
                for (int fIdx = 0; fIdx < Factors; fIdx++)
                {
                    double value = eigenVectors[rIdx, fIdx] * flooredRootEigenValues[fIdx];
                    fBetas[rIdx, fIdx] = value;
                    normalize += value * value;
                }

                if (normalize > 0.0)
                {
                    normalize = Math.Sqrt(1 / normalize);

                    for (int fIdx = 0; fIdx < Factors; fIdx++)
                    {
                        fBetas[rIdx, fIdx] *= normalize;
                    }
                }
            }
        }

        /// <summary>
        /// Copy the given thetas to the model; set the betas that the imply.
        /// </summary>
        /// <param name="x">Thetas to copy.</param>
        protected void SetThetas(double[] x)
        {
            for (int i = 0; i < x.Length; i++)
                fThetas[i] = x[i];

            for (int i = 0; i < fExpiries.Length; i++)
                SetBeta(i);
        }

        /// <summary>
        /// Retrieve the theta value corresponding to the given rate and theta indexes.
        /// </summary>
        /// <param name="rateIdx">Index of the forward rate</param>
        /// <param name="thetaIdx">Index of the theta for the given rate</param>
        /// <returns></returns>
        protected double GetTheta(int rateIdx, int thetaIdx)
        {
            int offset = rateIdx * (Factors - 1);
            return fThetas[offset + thetaIdx];
        }

        /// <summary>
        /// Set the values of the beta vectors (model drivers) at the given rate index position.
        /// </summary>
        /// <param name="rateIdx">Position in the beta vectors to set.</param>
        protected void SetBeta(int rateIdx)
        {
            if (Factors <= 1)
                fBetas[rateIdx, 1] = 1;
            else
            {
                double sinProduct = 1;
                for (int i = 0; i < Factors - 1; i++)
                {
                    double theta = GetTheta(rateIdx, i);
                    fBetas[rateIdx, i] = sinProduct * Math.Cos(theta);
                    sinProduct *= Math.Sin(theta);
                }
                fBetas[rateIdx, Factors - 1] = sinProduct;
            }
        }

        /// <summary>
        /// Calculate the correlation implied by the model between the given rates.
        /// </summary>
        /// <param name="rate1Idx">Index of the first rate</param>
        /// <param name="rate2Idx">Index of the second rate</param>
        /// <returns>The model correlation</returns>
        protected double ModelCorrelation(int rate1Idx, int rate2Idx)
        {
            double rho = 0;
            for (int k = 0; k < Factors; k++)
                rho += fBetas[rate1Idx, k] * fBetas[rate2Idx, k];

            return rho;
        }

        /// <summary>
        /// Calculates the distance between the model and target correlation matrixes.
        /// </summary>
        /// <returns>The distance</returns>
        protected double CalculateDistance()
        {
            double distance = 0;
            for (int i = 0; i < fExpiries.Length; i++)
            {
                for (int j = i + 1; j < fExpiries.Length; j++)
                {
                    double rho = ModelCorrelation(i, j);

                    double difference = fTarget[i, j] - rho;
                    distance += difference * difference;
                }
            }

            return distance;
        }
    }

    public abstract class MarketModelCorrelationCalibrationGeneric<M> : MarketModelCorrelationCalibrationBase, IModelCalibration
        where M : PriceModel, IMarketModelParameters
    {
        /// <summary>
        /// Calibrates the Market model using the historical statistics
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            M model = CalibrationHelper.Validate<M>(calibrationData, priceModel, output);

            // Calibration only suitable for scripting models and Interest Rate price factors.
            IPriceFactor priceFactor = model.GetPriceFactor();
            if (model.GetType() == typeof(MarketInterestRateModel) && priceFactor != null && priceFactor.GetType() != typeof(InterestRate))
                throw new CalibrationException(string.Format("{0} cannot calibrate model of {1}.", GetType().Name, priceFactor.GetType().Name));

            Initialize(model);
            SetTarget(calibrationData.Statistics, priceModel);
            FitThetas();
            Orthogonalize();

            // populate the model
            model.Betas.Clear();
            for (int bIdx = 0; bIdx < Factors; bIdx++)
            {
                RateCurve beta = new RateCurve();
                for (int eIdx = 0; eIdx < fExpiries.Length; eIdx++)
                    beta[fExpiries[eIdx]] = fBetas[eIdx, bIdx];

                model.Betas.Add(beta);
            }

            WriteCalibrationDataAsXML(output);
        }

        /// <summary>
        /// Rebuild output curve after manual changes to price model
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for</param>
        /// <param name="output">Revised output</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            M model = priceModel as M;
            WriteCalibrationDataAsXML(output);
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(M);
        }
    }

    [DisplayName("Market Model Parametric Correlation Calibration")]
    public class MarketModelCorrelationCalibrationParametric : MarketModelCorrelationCalibrationParametricGeneric<MarketInterestRateModel>
    {
    }

    public abstract class MarketModelCorrelationCalibrationParametricGeneric<M> : MarketModelCorrelationCalibrationGeneric<M>
        where M : PriceModel, IMarketModelParameters
    {
        protected Percentage fLongCorr = MinCorr;

        protected MarketModelCorrelationCalibrationParametricGeneric()
        {
            Long_Correlation = 0.65;
            Beta = 0.25;
            Gamma = 0.5;
        }

        public Percentage Long_Correlation
        {
            get { return fLongCorr; }
            set { fLongCorr = Math.Min(MaxCorr, Math.Max(MinCorr, value)); }
        }

        public double Beta
        {
            get; set;
        }

        public double Gamma
        {
            get; set;
        }

        /// <summary>
        /// Initialize the target correlation matrix from a correlation surface implied
        /// by the function
        /// 
        ///     rho(t, Ti, Tj) = long_rho + (1 - long_rho) exp(-beta abs((Ti-t)^gamma - (Tj-t)^gamma))
        /// </summary>
        protected override void SetTarget(IStatistics statistics, PriceModel priceModel)
        {
            for (int i = 0; i < fExpiries.Length; i++)
            {
                fTarget[i, i] = 1;
                double ti = Math.Pow(fExpiries[i], Gamma);
                for (int j = i + 1; j < fExpiries.Length; j++)
                {
                    double tj = Math.Pow(fExpiries[j], Gamma);
                    fTarget[i, j] = Long_Correlation + (1 - Long_Correlation) * Math.Exp(-Beta * Math.Abs(ti - tj));
                    fTarget[j, i] = fTarget[i, j];
                }
            }
        }
    }

    [DisplayName("Market Model Statistical Correlation Calibration")]
    public class MarketModelCorrelationCalibrationStatistical : MarketModelCorrelationCalibrationStatisticalGeneric<MarketInterestRateModel>
    {
    }

    public abstract class MarketModelCorrelationCalibrationStatisticalGeneric<M> : MarketModelCorrelationCalibrationGeneric<M>
        where M : PriceModel, IMarketModelParameters
    {
        /// <summary>
        /// Initialize the target correlation surface from a statistics set.
        /// </summary>
        protected override void SetTarget(IStatistics statistics, PriceModel priceModel)
        {
            List<int> indexes;
            List<string> points;

            statistics.FindEntries(MarketInterestRateModel.FwdRateName, priceModel.fID, out indexes, out points);

            ForwardRateDictionary forwardRates = new ForwardRateDictionary();
            forwardRates.SetEntries(indexes, points);

            for (int i = 0; i < fExpiries.Length; i++)
            {
                fTarget[i, i] = 1;
                int idxI = forwardRates.GetIndex(fExpiries[i], fTenors[i]);
                for (int j = i + 1; j < fExpiries.Length; j++)
                {
                    int idxJ = forwardRates.GetIndex(fExpiries[j], fTenors[j]);
                    fTarget[i, j] = statistics.Correlation(idxI, idxJ);
                    fTarget[j, i] = fTarget[i, j];
                }
            }
        }
    }
}
