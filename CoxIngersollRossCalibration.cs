using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Calibration for CIR Interest Rate model.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Cox-Ingersoll-Ross Interest Rate Calibration")]
    public class CIRInterestRateCalibration : StatisticsBasedCalibration, IModelCalibration
    {
        /// <summary>
        /// Calibrates the CIR Interest Rate Price model from the statistics.
        /// </summary>
        /// <remarks>
        /// Kappa and Theta are set to the average of statistics of each tenor.
        /// Sigma is a function of the reversion vol of the stat set and model's parameters.
        /// </remarks>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output,
                                      ErrorList errorList)
        {
            CIRInterestRateModel model = CalibrationHelper.Validate<CIRInterestRateModel>(calibrationData, priceModel, output);

            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(priceModel);

            // Retrieve the reversion statistics using diff return.
            MeanReversionStatistics[] meanReversionStatistics;
            TimeSeriesStructure<string> timeSeries = null;
            bool[] isValidStatistics;
            var calculationParameters = new MeanReversionStatisticsMultiPointCalculationParameters(ReturnMethod.Diff, null, null, null, null, false);
            if (!GetMeanReversionStatistics(calibrationData, priceModel, calculationParameters, errorList,
                                            ref timeSeries, out meanReversionStatistics, out isValidStatistics))
            {
                return;
            }

            if (!isValidStatistics[0])
            {
                errorList.Add(ErrorLevel.Error, string.Format("No valid statistics at first point for {0}.", priceFactor.GetKey()));
                return;
            }

            // Retrieve the data with the smallest term available in the statistics assuming it is sorted.
            model.Reversion_Speed = meanReversionStatistics[0].Alpha;
            model.Reversion_Level = meanReversionStatistics[0].Theta;

            double term = meanReversionStatistics[0].ParsePointIDToDouble();

            // Adjust the vol from the stat to account for the model square root type diffusion.
            double gamma = Math.Exp(-model.Reversion_Speed * term);
            double oneMinusGamma = 1.0 - gamma;
            using (var cache = Vector.Cache(1))
            {
                Vector initialShortRate = cache.Get();
                model.ComputeInitialShortRate(initialShortRate);
                Vector vol = cache.Get(meanReversionStatistics[0].Sigma * VectorMath.Sqrt(term * model.Reversion_Speed / 
                    ((initialShortRate * gamma + 0.5 * model.Reversion_Level * oneMinusGamma) * oneMinusGamma)));
                model.Vol = vol[0];

            }

            SetModelParameters(model);

            GetUpdatedCurves(priceModel, output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="priceModel">Price model to calculate curve(s) for.</param>
        /// <param name="output">Revised output.</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            CIRInterestRateModel model = (CIRInterestRateModel)priceModel;

            var modelCurve = BuildModelCurve(model);
            try
            {
                var observedCurve = ((IRateBase)model.GetPriceFactor()).Curve.Clone();
                observedCurve.MultiplyBy(Percentage.OverPercentagePoint);
                modelCurve.MultiplyBy(Percentage.OverPercentagePoint);

                var graphData = new ResultSeriesList
                {
                    fXAxis = { fLabel = "Time (Years)" },
                    fYAxis = { fLabel = "Zero Rate (%)" }
                };
                graphData.Add(new ResultSeries("Observed", LineStyle.Solid, observedCurve));
                graphData.Add(new ResultSeries("Model Implied", LineStyle.Solid, modelCurve));
                graphData.ToXml(output);

                CalibrationHelper.WriteModelParametersToXml(priceModel, output);
            }
            catch (XmlException e)
            {
                throw new CalibrationException(string.Format("Failed to calibrate {0}.", priceModel.GetKey()), e);
            }
        }

        /// <summary>
        /// The price model which this calibration is for.
        /// </summary>
        [DebuggerStepThrough]
        public Type ModelType()
        {
            return typeof(CIRInterestRateModel);
        }

        /// <summary>
        /// Set the deterministic shift parameters if needed.
        /// </summary>
        /// <remarks>
        /// This calibration doesn't apply to deterministic shift model.
        /// </remarks>
        protected virtual void SetModelParameters(CIRInterestRateModel model)
        {
            model.Use_Deterministic_Shift = YesNo.No;
        }

        /// <summary>
        /// Build the Zero Rate curve implied by the model at time 0.
        /// </summary>
        private static DoubleCurve BuildModelCurve(CIRInterestRateModel model)
        {
            var priceFactor = (IRateBase)model.GetPriceFactor();
            using (var cache = Vector.Cache(1))
            {
                Vector rate = cache.Get(model.GetInitialProcessValue());
                Vector df = cache.Get();

                var resultCurve = priceFactor.Curve.Clone();
                foreach (double x in priceFactor.Curve.X)
                {
                    if (x <= 0.0)
                        continue;

                    model.DiscountFactorValue(df, 0.0, x, model.Reversion_Speed, model.Reversion_Level, model.Vol, rate);
                    resultCurve[x] = -Math.Log(df[0]) / x;
                }

                return resultCurve;
            }
        }
    }

    /// <summary>
    /// Calibration for CIR Interest Rate model with deterministic shift.
    /// </summary>
    /// <remarks>
    /// If Use_Pre_Computed_Statistics is "Yes", we try to get statistics from the set and if those are not available, 
    /// we try to calculate the statistics from the archive file. 
    /// If Use_Pre_Computed_Statistics is "No", we try to calculate the statistics from the archive file and if those 
    /// are not available, we try to get statistics from the set.
    /// </remarks>
    [DisplayName("Cox-Ingersoll-Ross Deterministic Shift Interest Rate Calibration")]
    public class CIRDeterministicShiftInterestRateCalibration : CIRInterestRateCalibration
    {
        /// <summary>
        /// Set the deterministic shift parameters.
        /// </summary>
        protected override void SetModelParameters(CIRInterestRateModel model)
        {
            using (var cache = Vector.Cache(VectorCurveHelpers.TargetVectorSize(model.GetPriceFactor())))
            {
                model.Use_Deterministic_Shift = YesNo.Yes;

                Vector initialRate = cache.Get();
                model.ComputeInitialShortRate(initialRate);
                Vector initialProcessValue = cache.Get(0.5 * initialRate); //arbitrary value between 0 and r0
                model.Initial_Process_Value = initialProcessValue[0];
            }
        }
    }
}