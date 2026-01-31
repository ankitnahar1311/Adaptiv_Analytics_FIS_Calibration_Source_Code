/// <author>
/// Justin Chan
/// </author>
/// <owner>
/// Phil Koop
/// </owner>
/// <summary>
/// Calibration of one-factor volatility models using Maximum Likelihood Estimation techniques.
/// </summary>
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Generic base class for calibration of one-factor volatility models.
    /// </summary>
    /// <typeparam name="TVolSurface">A VolSurface object.</typeparam>
    /// <typeparam name="TPriceFactor">A Price Factor.</typeparam>
    /// <typeparam name="TRiskFactor">An Ornstein Uhlenbeck Process.</typeparam>
    [Serializable]
    public abstract class SimpleVolMLECalibration<TVolSurface, TPriceFactor, TRiskFactor>
        where TVolSurface : VolSurface, new()
        where TPriceFactor : class, IPriceFactor
        where TRiskFactor : OrnsteinUhlenbeckProcess, new()
    {
        // Selection window on volatility surface to be used in calibration
        protected readonly SurfaceWindow<TVolSurface> fVolatilityWindow; 

        protected TimeSeriesStructure<VolSurface.Point> fVolSurfaceStructure;
        protected List<TDate> fHistoricalDates;
        protected double[] fProcessValues;
        protected double[] fProcessValuesAtHorizon;
        private bool[] fIsValidCoordinate; // fIsValidCoordinate[i] true if fCoordinates lies within the window

        /// <summary>
        /// Initializes a new instance of the SimpleVolMLECalibration class.
        /// </summary>
        protected SimpleVolMLECalibration()
        {
            fVolatilityWindow = new SurfaceWindow<TVolSurface>(VolatilitySurfaceDimension());

            // Use MLE by default.
            Calibration_Method  = CalibrationMethod.MLE;
            MLE_Parameters      = new SimpleVolMLEParameters<TVolSurface>(fVolatilityWindow);
        }

        /// <summary>
        /// Calibration method can be statistics or MLE.
        /// </summary>
        public CalibrationMethod Calibration_Method                 { get; set; }

        /// <summary>
        /// Blending Method for model calibration.
        /// </summary>
        public BlendMethod Blend_Method                             { get; set; }

        /// <summary>
        /// Parameters for the MLE method.
        /// </summary>
        public SimpleVolMLEParameters<TVolSurface> MLE_Parameters   { get; set; }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public abstract int VolatilitySurfaceDimension();

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Step (1): Validation
            fVolatilityWindow.Validate();
            SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor> model =
                CalibrationHelper.Validate<SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor>>(calibrationData, priceModel, output);

            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            if (model.Theta().Surface.Dimension != VolatilitySurfaceDimension() ||
                model.InitialSurface().Dimension != VolatilitySurfaceDimension())
                throw new CalibrationException(
                    string.Format("This calibrator assumes volatility surface for {0} to be {1}-dimensional", priceFactor.GetKey(), VolatilitySurfaceDimension()));

            if (!GenericDataRetrieval.TryGetHistoricalVolSurface<TVolSurface>(calibrationData, priceFactor,
                MLE_Parameters.Data_Retrieval_Parameters, errorList, out fVolSurfaceStructure))
            {
                return;
            }

            if (MLE_Parameters.Detrended_Time_Series == YesNo.Yes)
            {
                double[] drifts;
                IHolidayCalendar calendar = CalibrationHelper.GetHolidayCalendar(calibrationData, errorList, MLE_Parameters.Data_Retrieval_Parameters.Calendar);
                if (!TimeSeriesAnalysis.TryDetrendExponential(fVolSurfaceStructure, errorList, calendar, out drifts))
                {
                    return;
                }
            }

            fHistoricalDates = new List<TDate>(fVolSurfaceStructure.Dates);
            var coordinates  = new List<VolSurface.Point>(fVolSurfaceStructure.TimeSeriesIDs);

            var validCoordinates = fVolatilityWindow.GetPointsInWindow(coordinates);
            if (!validCoordinates.Any())
                throw new CalibrationException(string.Format("No historical data for {0} in selected window.", priceFactor.GetKey()));

            fIsValidCoordinate = coordinates.Select(point => validCoordinates.Contains(point)).ToArray();

            // Step(4): calculate observations of average process
            fProcessValues          = new double[fVolSurfaceStructure.Length];
            fProcessValuesAtHorizon = new double[fVolSurfaceStructure.Length];
            for (int i = 0; i < fVolSurfaceStructure.Length; i++)
            {
                double sum          = 0.0;
                double sumAtHorizon = 0.0;
                for (int j = 0; j < fVolSurfaceStructure.NumTimeSeries; j++)
                {
                    if (!fIsValidCoordinate[j])
                        continue;

                    sum          += Math.Log(fVolSurfaceStructure.Values[j][i]);
                    sumAtHorizon += Math.Log(fVolSurfaceStructure.ValuesAtHorizon[j][i]);
                }

                fProcessValues[i]          = sum / validCoordinates.Count;
                fProcessValuesAtHorizon[i] = sumAtHorizon / validCoordinates.Count;
            }

            // Step(5): calculate alpha and sigma of average process 
            var alphaFixed = MLE_Parameters.Condition_On_Reversion_Speed == YesNo.Yes ? (double?)MLE_Parameters.Reversion_Speed : null;
            MeanReversionStatistics mrs;
            if (!ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(fProcessValues, fProcessValuesAtHorizon, ReturnMethod.Diff,
                                                                                    fVolSurfaceStructure.DeltaT, alphaFixed, null, true, errorList,
                                                                                    out mrs))
            {
                errorList.Add(ErrorLevel.Error, "Calibration Error");
                return;
            }

            double alpha = mrs.Alpha;
            double sigma = mrs.Sigma;

            // Step(6): Assigning parameters to model
            TVolSurface longRunMeanSurface = model.Theta(); 
            longRunMeanSurface.Surface.Clear();
            for (int i = 0; i < fVolSurfaceStructure.NumTimeSeries; i++)
            {
                // Calculate theta of point using alpha of average process
                if (!ReversionStatisticsCalculationHelper.ComputeMleReversionStatistics(fVolSurfaceStructure.Values[i], fVolSurfaceStructure.ValuesAtHorizon[i],
                                                                                        ReturnMethod.Log, fVolSurfaceStructure.DeltaT, alpha, null, false, errorList,
                                                                                        out mrs))
                {
                    errorList.Add(ErrorLevel.Error, "Calibration Error");
                    return;
                }

                mrs.Sigma = sigma;
                mrs.TransformToThetaOfGouProcess();
                longRunMeanSurface.SetValue(coordinates[i], mrs.Theta);
            }

            model.Blend_Method    = Blend_Method;
            model.Reversion_Speed = alpha;
            model.Vol_of_Vol      = sigma;

            // Step(7): Special cases: optionally fitting the long run mean surface to current vol surface
            if (MLE_Parameters.Use_Current_Smile == YesNo.Yes)
                FitLongRunMeanSurfaceToCurrentSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleVolMLECalibration<TVolSurface, TPriceFactor, TRiskFactor>))
                GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        public virtual void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            if (fVolSurfaceStructure == null || fVolSurfaceStructure.Length < StatisticsCalculationHelper.MinTimeSeriesLength)
                return;

            SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor> model = (SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor>)priceModel;
            var graphData = new ResultSeriesList
                                {
                                    fXAxis = { fLabel = "Date" },
                                    fYAxis = { fLabel = "Historical Volatilities(%)" }
                                };

            graphData.fXAxis.fAxisType = AxisType.Date;
            
            var coordinates = fVolSurfaceStructure.TimeSeriesIDs;
            var dates       = fVolSurfaceStructure.Dates;
            var values      = fVolSurfaceStructure.Values;

            for (int i = 0; i < coordinates.Length; i++)
            {
                if (!fIsValidCoordinate[i])
                    continue;

                // For each point on the grid surface, plot a curve
                var historicalData = new Curve();
                for (int j = 0; j < dates.Length; j++)
                    historicalData[(float)dates[j]] = values[i][j] * Percentage.OverPercentagePoint;
                
                graphData.Add(fVolatilityWindow.fDimension == 2
                                  ? new ResultSeries(string.Format(CultureInfo.InvariantCulture, "({0}, {1})",
                                                      coordinates[i].Moneyness, coordinates[i].Expiry), 
                                                      LineStyle.Dashed, historicalData)
                                  : new ResultSeries(string.Format(CultureInfo.InvariantCulture, "({0}, {1}, {2})",
                                                      coordinates[i].Moneyness, coordinates[i].Expiry, coordinates[i].Tenor), 
                                                      LineStyle.Dashed, historicalData));
            }

            // Graphical aide to visualise results; plotting one standard deviation range around the average of all long run means.
            TVolSurface longRunMeanSurface = model.Theta();

            var iter                        = new SurfaceIter(longRunMeanSurface.Surface);
            var args                        = new double[longRunMeanSurface.Surface.Dimension];
            double sumVolLongRunMean        = 0.0;
            int numPointsLongRunMean        = 0;
            var longRunOneStdevUpperRange   = new Curve();
            var longRunOneStdevLowerRange   = new Curve();

            while (!iter.AtEnd)
            {
                double currentVol = 0.0;
                iter.GetPoint(ref args, ref currentVol);

                sumVolLongRunMean += currentVol;
                numPointsLongRunMean++;

                iter.Next();
            }

            double avgLongRunMean = sumVolLongRunMean / numPointsLongRunMean;

            if (model.Reversion_Speed > 0.0)
            {
                foreach (TDate date in dates)
                {
                    longRunOneStdevUpperRange[date] = avgLongRunMean * Math.Exp(model.Vol_of_Vol / Math.Sqrt(2.0 * model.Reversion_Speed)) * Percentage.OverPercentagePoint;
                    longRunOneStdevLowerRange[date] = avgLongRunMean * Math.Exp(-model.Vol_of_Vol / Math.Sqrt(2.0 * model.Reversion_Speed)) * Percentage.OverPercentagePoint;
                }

                graphData.Add(new ResultSeries(string.Format("Long Run One Standard Deviation Upper Range"), LineStyle.Solid, longRunOneStdevUpperRange));
                graphData.Add(new ResultSeries(string.Format("Long Run One Standard Deviation Lower Range"), LineStyle.Solid, longRunOneStdevLowerRange));
            }

            graphData.ToXml(output);
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected abstract void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<TVolSurface, TPriceFactor, TRiskFactor> model, bool useSinglePoint, double[] singlePoint);
    }

    /// <summary>
    /// N-dimensional window of volatiilty surface, e.g., to be used as selection window of calibration on historical surfaces.
    /// </summary>
    /// <remarks>
    /// The bounds can degenerate into a single point. If the values of all bounds are zero, then it is assumed that entire surface is to be used.
    /// </remarks>
    /// <typeparam name="TVolSurface">A VolSurface object.</typeparam>
    public class SurfaceWindow<TVolSurface> where TVolSurface : VolSurface, new()
    {
        /// <summary>
        /// Returns the number of dimensions of the window; getter only.
        /// </summary>
        public readonly int fDimension;

        private readonly DoubleArray fUpperBound;
        private readonly DoubleArray fLowerBound;

        /// <summary>
        /// Initializes a new instance of the SurfaceWindow class.
        /// </summary>
        public SurfaceWindow(int dim)
        {
            if (dim > 0)
                fDimension = dim;
            else
                throw new AnalyticsException(
                    string.Format("Surface selection window dimension {0} cannot be non-positive.", dim));

            fUpperBound = new DoubleArray(fDimension, fDimension);
            fLowerBound = new DoubleArray(fDimension, fDimension);
        }

        /// <summary>
        /// Validates the bounds of the surface window.
        /// </summary>
        public void Validate()
        {
            for (int i = 0; i < fDimension; i++)
                if (fUpperBound[i] < fLowerBound[i])
                    throw new AnalyticsException(
                        string.Format(
                            "Invalid surface selection window: upper bound {0} cannot be smaller than lower bound {1}.",
                            fUpperBound[i], fLowerBound[i]));
        }

        /// <summary>
        /// Get value of upper bound of a given dimension.
        /// </summary>
        public double GetUpperBound(int dim)
        {
            if (dim < fUpperBound.Count)
                return fUpperBound[dim];

            return float.NaN;
        }

        /// <summary>
        /// Set value of upper bound of a given dimension.
        /// </summary>
        public void SetUpperBound(int dim, double value)
        {
            if (dim < fUpperBound.Count)
                fUpperBound[dim] = value;
        }

        /// <summary>
        /// Get value of lower bound of a given dimension.
        /// </summary>
        public double GetLowerBound(int dim)
        {
            if (dim < fLowerBound.Count)
                return fLowerBound[dim];

            return double.NaN;
        }

        /// <summary>
        /// Set value of lower bound of a given dimension.
        /// </summary>
        public void SetLowerBound(int dim, double value)
        {
            if (dim < fLowerBound.Count)
                fLowerBound[dim] = value;
        }

        /// <summary>
        /// Returns true if no particular window on surface is selected.
        /// </summary>
        public bool UseEntireSurface()
        {
            foreach (var bound in fUpperBound)
                if (Math.Abs(bound - 0.0) > CalcUtils.TINY)
                    return false;

            foreach (var bound in fLowerBound)
                if (Math.Abs(bound - 0.0) > CalcUtils.TINY)
                    return false;

            return true;
        }

        /// <summary>
        /// Returns true if selection window is single point.
        /// </summary>
        public bool UseSinglePoint()
        {
            if (UseEntireSurface())
                return false;

            for (int i = 0; i < fUpperBound.Count; i++)
                if (Math.Abs(fUpperBound[i] - fLowerBound[i]) > CalcUtils.TINY)
                    return false;

            return true;
        }

        /// <summary>
        /// Returns the single point, assuming the whole selection window is single point.
        /// </summary>
        public double[] GetSinglePointArgs()
        {
            return fUpperBound.ToArray();
        }

        /// <summary>
        /// Given a list of points, return the subset of which falls within the window.
        /// </summary>
        public List<VolSurface.Point> GetPointsInWindow(List<VolSurface.Point> points)
        {
            if (UseEntireSurface() || points.Count == 0)
                return points;

            var validPoints = new List<VolSurface.Point>();
            var surface = new TVolSurface();

            if (UseSinglePoint())
            {
                validPoints.Add(ClosestPoint(surface.PointFromArgs(GetSinglePointArgs()), points));
            }
            else
            {
                foreach (VolSurface.Point point in points)
                    if (IsPointInWindow(point))
                        validPoints.Add(point);
            }

            return validPoints;
        }

        /// <summary>
        /// Returns true if the given point is in the selection window.
        /// </summary>
        private bool IsPointInWindow(VolSurface.Point gridPoint)
        {
            for (int i = 0; i < fDimension; i++)
            {
                switch (i)
                {
                    case 0:
                        if (gridPoint.Moneyness < fLowerBound[i] || fUpperBound[i] < gridPoint.Moneyness)
                            return false;
                        break;
                    case 1:
                        if (gridPoint.Expiry < fLowerBound[i] || fUpperBound[i] < gridPoint.Expiry)
                            return false;
                        break;
                    case 2:
                        if (gridPoint.Tenor < fLowerBound[i] || fUpperBound[i] < gridPoint.Tenor)
                            return false;
                        break;
                    default:
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Returns the closest point (among a given list of candidates) to a given point.
        /// </summary>
        /// <param name="singlePoint">The single point to which we wish to find the closest point.</param>
        /// <param name="points">The list of points among which we wish to find hte closest point.</param>
        private VolSurface.Point ClosestPoint(VolSurface.Point singlePoint, List<VolSurface.Point> points)
        {
            if (points.Count == 0)
                throw new CalibrationException(string.Format("Given list of vol surface points cannot be empty."));

            VolSurface.Point closestPoint = points[0];
            double minDistance = double.MaxValue;
            foreach (VolSurface.Point point in points)
            {
                double distance = VolSurface.Point.Distance(singlePoint, point);
                if (!(distance < minDistance))
                    continue;

                minDistance = distance;
                closestPoint = point;
            }

            return closestPoint;
        }
    }

    /// <summary>
    /// Parameters for the MLE method.
    /// </summary>
    /// <typeparam name="TVolSurface">A VolSurface object.</typeparam>
    public class SimpleVolMLEParameters<TVolSurface> : NestedPresentableObject where TVolSurface : VolSurface, new()
    {
        private readonly SurfaceWindow<TVolSurface> fVolatilityWindow;
        private double fAlpha;

        /// <summary>
        /// Initialise parameters.
        /// </summary>
        public SimpleVolMLEParameters(SurfaceWindow<TVolSurface> surfaceWindow)
        {
            Data_Retrieval_Parameters       = new CalibratorDataRetrievalParameters<VolSurface.Point>();
            Detrended_Time_Series           = YesNo.Yes;
            Condition_On_Reversion_Speed    = YesNo.No;
            Use_Current_Smile               = YesNo.No;

            fVolatilityWindow               = surfaceWindow;
        }

        /// <summary>
        /// Parameters for historical data retrieval.
        /// </summary>
        public CalibratorDataRetrievalParameters<VolSurface.Point> Data_Retrieval_Parameters { get; set; }

        /// <summary>
        /// Optionally remove the timeseries' drift before MLE calibration; this is a calibration setting.
        /// </summary>
        public YesNo Detrended_Time_Series { get; set; }

        /// <summary>
        /// Optionally using conditional Reversion Speed (i.e., Alpha); output parameters will be conditional on the given Alpha.
        /// </summary>
        public YesNo Condition_On_Reversion_Speed { get; set; }

        /// <summary>
        /// Optional user-defined Reversion_Speed (i.e., Alpha) for conditional MLE estimation; used in conjunction with Condition_On_Reversion_Speed.
        /// </summary>
        public double Reversion_Speed
        {
            get { return fAlpha; }
            set { fAlpha = Math.Abs(value) >= CalcUtils.MinReversionRate ? value : 0.0; }
        }

        /// <summary>
        /// Optional user-defined upper bound of moneyness of the window of historical volatility surfaces to be used for calibration.
        /// </summary>
        public double Moneyness_Upper_Bound
        {
            get
            {
                return fVolatilityWindow.GetUpperBound(0);
            }
            set
            {
                if (value >= 0.0)
                    fVolatilityWindow.SetUpperBound(0, value);
                else
                    throw new AnalyticsException(string.Format("Moneyness upper bound {0} cannot be negative", value));
            }
        }

        /// <summary>
        /// Optional user-defined lower bound of moneyness of the window of historical volatility surfaces to be used for calibration.
        /// </summary>
        public double Moneyness_Lower_Bound
        {
            get
            {
                return fVolatilityWindow.GetLowerBound(0);
            }
            set
            {
                if (value >= 0.0)
                    fVolatilityWindow.SetLowerBound(0, value);
                else
                    throw new AnalyticsException(string.Format("Moneyness lower bound {0} cannot be negative", value));
            }
        }

        /// <summary>
        /// Optional user-defined upper bound of expiry of the window of historical volatility surfaces to be used for calibration.
        /// </summary>
        public Period Expiry_Upper_Bound
        {
            get
            {
                return fVolatilityWindow.GetUpperBound(1);
            }
            set
            {
                if (value >= 0.0)
                    fVolatilityWindow.SetUpperBound(1, (float)value);
                else
                    throw new AnalyticsException(string.Format("Expiry upper bound {0} cannot be negative", value));
            }
        }

        /// <summary>
        /// Optional user-defined lower bound of expiry of the window of historical volatility surfaces to be used for calibration.
        /// </summary>
        public Period Expiry_Lower_Bound
        {
            get
            {
                return fVolatilityWindow.GetLowerBound(1);
            }
            set
            {
                if (value >= 0.0)
                    fVolatilityWindow.SetLowerBound(1, (float)value);
                else
                    throw new AnalyticsException(string.Format("Expiry lower bound {0} cannot be negative", value));
            }
        }

        /// <summary>
        /// Optional user-defined choice for using the current shape of the smile of volatility for the shape of the long run mean surface (i.e., Theta surface of the model).
        /// </summary>
        public YesNo Use_Current_Smile { get; set; }
    }

    /// <summary>
    /// Model Calibrator for Simple Asset Price Vol Model.
    /// </summary>
    [DisplayName("Simple Asset Price Volatility Calibration")]
    public class SimpleAssetPriceVolCalibration : SimpleVolMLECalibration<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimpleAssetPriceVolStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleAssetPriceVolCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public virtual Type ModelType()
        {
            return typeof(SimpleAssetPriceVolModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 2;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess> model, 
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }

    /// <summary>
    /// Model Calibrator for Simple FX Vol Model.
    /// </summary>
    [DisplayName("Simple FX Volatility Calibration")]
    public class SimpleFXVolCalibration : SimpleVolMLECalibration<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimpleFxVolStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleFXVolCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SimpleFXVolModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 2;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess> model, 
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }

    /// <summary>
    /// Calibration of the SimplePriceIndexVolatilityModel.
    /// </summary>
    [DisplayName("Simple Price Index Volatility Calibration")]
    public class SimplePriceIndexVolatilityCalibration : SimpleVolMLECalibration<PriceIndexVolSurface, IPriceIndexVolatility, PriceIndexVolatilityProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>MLE
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimplePriceIndexVolatilityStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimplePriceIndexVolatilityCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SimplePriceIndexVolatilityModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 3;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<PriceIndexVolSurface, IPriceIndexVolatility, PriceIndexVolatilityProcess> model, 
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }

    /// <summary>
    /// Calibration of the SimpleForwardPriceVolatilityModel.
    /// </summary>
    [DisplayName("Simple Forward Price Volatility Calibration")]
    public class SimpleForwardPriceVolatilityCalibration : SimpleVolMLECalibration<ForwardPriceVolSurface, IForwardPriceVol, ForwardPriceVolProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>MLE
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimpleForwardPriceVolatilityStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleForwardPriceVolatilityCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SimpleForwardPriceVolatilityModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 3;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<ForwardPriceVolSurface, IForwardPriceVol, ForwardPriceVolProcess> model,
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }

    /// <summary>
    /// Calibration of the SimpleInterestRateVolModel.
    /// </summary>
    [DisplayName("Simple Interest Rate Volatility Calibration")]
    public class SimpleInterestRateVolCalibration : SimpleVolMLECalibration<InterestVolSurface, IInterestVol, InterestRateVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>MLE
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimpleInterestRateVolStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleInterestRateVolCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SimpleInterestRateVolModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 3;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<InterestVolSurface, IInterestVol, InterestRateVolOUProcess> model,
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }

    /// <summary>
    /// Calibration of the SimpleInterestYieldVolModel.
    /// </summary>
    [DisplayName("Simple Interest Yield Volatility Calibration")]
    public class SimpleInterestYieldVolCalibration : SimpleVolMLECalibration<InterestVolSurface, IInterestVol, InterestYieldVolOUProcess>, IModelCalibration
    {
        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>MLE
        public override void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            // Calibrated using statistics-based method if selected.
            if (Calibration_Method == CalibrationMethod.Pre_Computed_Statistics)
            {
                var statisticsCalibration = new SimpleInterestYieldVolStatisticsStatisticsCalibration
                {
                    Blend_Method = Blend_Method
                };

                statisticsCalibration.Calibrate(calibrationData, priceModel, output, errorList);
                return;
            }

            base.Calibrate(calibrationData, priceModel, output, errorList);

            // If a child calls base.Calibrate(...), the base's graphing function may not be appropriate anymore.
            if (GetType() == typeof(SimpleInterestYieldVolCalibration))
                base.GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Type of price factor model to calibrate.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(SimpleInterestYieldVolModel);
        }

        /// <summary>
        /// Number of dimensions on the volatilty surface.
        /// </summary>
        /// <remarks>
        /// Any concrete child needs to specify this upfront because the user has to specify upper/lower bounds of volatility window used in calibration.
        /// This base class assumes at least the surface is 2D. If a 3D surface is to be used, additional public properties (e.g., Tenor_Lower_Bound)
        /// needs to be added to the child.
        /// </remarks>
        public override int VolatilitySurfaceDimension()
        {
            return 3;
        }

        /// <summary>
        /// Fits the smileness of the long run mean surface to current volatility surface.
        /// </summary>
        protected override void FitLongRunMeanSurfaceToCurrentSurface(SimpleVolModel<InterestVolSurface, IInterestVol, InterestYieldVolOUProcess> model,
            bool useSinglePoint, double[] singlePoint)
        {
            CalibrationHelper.CalculateLongRunMeanSurface(model, fVolatilityWindow.UseSinglePoint(), fVolatilityWindow.GetSinglePointArgs());
        }
    }
}