/// <author>
/// Alastair Wilkins
/// </author>
/// <owner>
/// Alastair Wilkins
/// </owner>
/// <summary>
/// Calibration of one-factor volatility models.
/// </summary>
using System;
using System.Collections.Generic;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Generic base class for calibration of one-factor volatility models.
    /// </summary>
    [Serializable]
    public abstract class SimpleVolStatisticsCalibration<SurfaceClass, PriceClass, RiskClass>
        where SurfaceClass : VolSurface, new()
        where PriceClass : class, IPriceFactor
        where RiskClass : OrnsteinUhlenbeckProcess, new()
    {
        /// <summary>
        /// Blend method.
        /// </summary>
        public BlendMethod Blend_Method { get; set; }

        /// <summary>
        /// Calibrates the price factor model.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            SimpleVolModel<SurfaceClass, PriceClass, RiskClass> model =
                CalibrationHelper.Validate<SimpleVolModel<SurfaceClass, PriceClass, RiskClass>>(calibrationData, priceModel, output);

            CalibrateSimpleVolModel(calibrationData.Statistics, model);
            GetUpdatedCurves(output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        private void GetUpdatedCurves(XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// Calculate the PCA model parameters from the statistics of log changes in volatility.
        /// </summary>
        private void CalibrateSimpleVolModel(IStatistics statistics, SimpleVolModel<SurfaceClass, PriceClass, RiskClass> model)
        {
            IPriceFactor priceFactor = CalibrationHelper.GetPriceFactor(model);

            List<int> indexes;
            List<string> pointStrings;
            List<double> reversionVols, reversionSpeeds, longRunMeans;
            StatisticsHelper.GetMeanReversionStatistics(statistics, priceFactor.GetType(), model.fID, ReturnMethod.Log, out indexes, out pointStrings, out reversionVols, out reversionSpeeds, out longRunMeans);

            // Number of points with non-zero volatility
            int nPoints = indexes.Count;
            if (nPoints == 0)
                throw new CalibrationException(string.Format("No statistics for {0}.", priceFactor.GetKey()));

            SurfaceClass theta = model.Theta();
            double atmMoneyness = theta.AtmMoneyness();

            // Parse point strings, calculate vol of vol, mean reversion speed and set long run mean surface
            var points = new VolSurface.Point[nPoints];
            double sumVolOfVols = 0.0;
            double sumReversionSpeed = 0.0;
            theta.Surface.Clear();
            bool allAtm = true;

            for (int i = 0; i < nPoints; ++i)
            {
                points[i] = model.Theta().PointFromString(pointStrings[i]);
                sumReversionSpeed += reversionSpeeds[i];
                sumVolOfVols += reversionVols[i];
                theta.SetValue(points[i], longRunMeans[i]);

                if (allAtm && points[i].Moneyness != atmMoneyness)
                    allAtm = false;
            }

            model.Vol_of_Vol = sumVolOfVols / nPoints;
            model.Reversion_Speed = sumReversionSpeed / nPoints;
            model.Blend_Method = Blend_Method;

            if (allAtm)
                CalibrationHelper.CalculateLongRunMeanSurface(model, false, null);
        }
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimpleInterestRateVolStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<InterestVolSurface, IInterestVol, InterestRateVolOUProcess>
    {
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimpleInterestYieldVolStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<InterestVolSurface, IInterestVol, InterestYieldVolOUProcess>
    {
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimpleAssetPriceVolStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess>
    {
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimpleFxVolStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<AssetPriceVolSurface, ISpotProcessVol, OrnsteinUhlenbeckProcess>
    {
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimpleForwardPriceVolatilityStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<ForwardPriceVolSurface, IForwardPriceVol, ForwardPriceVolProcess>
    {
    }

    /// <summary>
    /// Statistics-based calibration.
    /// </summary>
    public class SimplePriceIndexVolatilityStatisticsStatisticsCalibration : SimpleVolStatisticsCalibration<PriceIndexVolSurface, IPriceIndexVolatility, PriceIndexVolatilityProcess>
    {
    }
}
