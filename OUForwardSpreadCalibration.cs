/// <author>
/// Philip Koop
/// </author>
/// <owner>
/// Philip Koop, Al Wilkins
/// </owner>
/// <summary>
/// Statistical calibrations for OU forward spread models.
/// </summary>
using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Base class for statistical calibratiohn of OU models.
    /// </summary>
    public abstract class OUForwardSpreadStatisticsCalibrationBase
    {
        /// <summary>
        /// Base date used in calculation.
        /// </summary>
        public TDate Base_Date          { get; set; }

        /// <summary>
        /// Last rollover data.
        /// </summary>
        public TDate Last_Rollover_Date { get; set; }

        /// <summary>
        /// Calibrates the price model from the statistics files.
        /// </summary>
        public void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            ForwardPriceModel model = CalibrationHelper.Validate<ForwardPriceModel>(calibrationData, priceModel, output);

            Calculate(calibrationData.Statistics, model, errorList);
            GetUpdatedCurves(output);
        }

        /// <summary>
        /// Rebuilds output curve after manual changes to price model.
        /// </summary>
        /// <param name="output">Revised output.</param>
        private void GetUpdatedCurves(XmlWriter output)
        {
            var graphData = new ResultSeriesList();
            graphData.ToXml(output);
        }

        /// <summary>
        /// Set the volatility and reversion speed from statistics.
        /// </summary>
        /// <remarks>
        /// We must convert the relative time of the statistics to absolute date format.
        /// </remarks>
        protected abstract void SetStatisticalParameters(List<int> indexes, List<double> points, IStatistics statistics, ForwardPrice price, ForwardPriceModel model);

        /// <summary>
        /// Calculate the model parameters.
        /// </summary>
        /// <param name="statistics">statistics set.</param>
        /// <param name="priceModel">model to calibrate.</param>
        /// <param name="errorList">List of errors, warnings and infos.</param>
        private void Calculate(IStatistics statistics, ForwardPriceModel priceModel, ErrorList errorList)
        {
            ForwardPrice price = (ForwardPrice)CalibrationHelper.GetPriceFactor(priceModel);

            if (price.Curve.Count < 1)
                throw new CalibrationException("No points in spread curve to be modeled.");

            // Fetch the statistics for the forward spread to be modeled
            string       priceType = price.GetType().Name;
            List<int>    indexes;
            List<double> points;
            StatisticsHelper.FindEntries(statistics, priceType, priceModel.fID, out indexes, out points);
            if (!indexes.Any())
            {
                errorList.Add(ErrorLevel.Error, string.Format("No statistics for {0}.", price.GetKey()));
                return;
            }

            // Set volatilty and reversion speed from statistics
            SetStatisticalParameters(indexes, points, statistics, price, priceModel);
        }
    }

    /// <summary>
    /// Calibrate the OU Forward Spread TS model from statistics.
    /// </summary>
    public class OUForwardSpreadTSStatisticsCalibration : OUForwardSpreadStatisticsCalibrationBase
    {
        /// <summary>
        /// Set the volatility and reversion speed from statistics.
        /// </summary>
        /// <remarks>
        /// We must convert the relative time of the statistics to absolute date format.
        /// </remarks>
        protected override void SetStatisticalParameters(List<int> indexes, List<double> points, IStatistics statistics, ForwardPrice price, ForwardPriceModel model)
        {
            var volatility     = new PriceCurve();
            var reversionSpeed = new PriceCurve();
            for (int i =  0; i < indexes.Count; ++ i)
            {
                TDate pointDate = Base_Date + CalcUtils.YearsToDays(points[i]);
                int   index     = indexes[i];

                volatility[pointDate]     = StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.ReversionVolString, ReturnMethod.Diff);
                reversionSpeed[pointDate] = StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.MeanRevSpeedString, ReturnMethod.Diff);
            }

            var forwardSpreadModel = (OUForwardSpreadTSModel)model;

            forwardSpreadModel.Volatility      = volatility;
            forwardSpreadModel.Reversion_Speed = reversionSpeed;

            // Set mean reversion level from the price factor curve; that curve is already in absolute date format
            var reversionLevel = new PriceCurve();
            reversionLevel.Assign(price.Curve);
            forwardSpreadModel.Reversion_Level = reversionLevel;

            // Set the last rollover date from the calibration parameter
            forwardSpreadModel.Last_Rollover_Date = Last_Rollover_Date;
        }
    }

    /// <summary>
    /// Calibrate the OU Forward Spread model from statistics.
    /// </summary>
    public class OUForwardSpreadStatisticsCalibration : OUForwardSpreadStatisticsCalibrationBase
    {
        /// <summary>
        /// Set the volatility and reversion speed from statistics.
        /// </summary>
        /// <remarks>
        /// The reversion speed is set to the simple average of the point reversion speeds. We must convert the relative time of the
        /// statistics to absolute date format.
        /// </remarks>
        protected override void SetStatisticalParameters(List<int> indexes, List<double> points, IStatistics statistics, ForwardPrice price, ForwardPriceModel model)
        {
            var volatility = new PriceCurve();
            double sumOfReversionSpeeds = 0.0;
            for (int i = 0; i < indexes.Count; ++i)
            {
                TDate pointDate = Base_Date + CalcUtils.YearsToDays(points[i]);
                int index = indexes[i];

                volatility[pointDate] = StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.ReversionVolString, ReturnMethod.Diff);
                sumOfReversionSpeeds += StatisticsHelper.GetStatistic(statistics, index, BasePlugInStatisticsCalculations.MeanRevSpeedString, ReturnMethod.Diff);
            }

            var forwardSpreadModel = (OUForwardSpreadModel)model;

            forwardSpreadModel.Volatility       = volatility;
            forwardSpreadModel.Reversion_Speed  = sumOfReversionSpeeds / indexes.Count;

            // Set mean reversion level from the price factor curve; that curve is already in absolute date format
            var reversionLevel = new PriceCurve();
            reversionLevel.Assign(price.Curve);
            forwardSpreadModel.Reversion_Level = reversionLevel;

            // Set the last rollover date from the calibration parameter
            forwardSpreadModel.Last_Rollover_Date = Last_Rollover_Date;
        }
    }
}