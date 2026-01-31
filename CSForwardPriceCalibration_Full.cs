using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Xml;
using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Pair of doubles representing a tenor and the historical volatility for that tenor. 
    /// </summary>
    public struct TenorAndHistoricalVol
    {
        public readonly double fTenor;
        public readonly double fHistoricalVol;

        /// <summary>
        /// Construct a new instance from parameters.
        /// </summary>
        public TenorAndHistoricalVol(double tenor, double historicalVol)
        {
            fTenor = tenor;
            fHistoricalVol = historicalVol;
        }
    }

    [DisplayName("Clewlow-Strickland Forward Price Full Calibration")]
    public class CSForwardPriceCalibration_Full : StatisticsBasedCalibration, IModelCalibration
    {
        protected int fNumSims = 100;

        public int Number_Of_Iterations
        {
            get { return fNumSims; } set { fNumSims = Math.Max(value, 1); }
        }

        /// <summary>
        /// Core method for the calibration. This function takes the historical data and price model, performs analysis of
        /// the historical data and sets the price model parameters Alpha and Sigma to those which produce a variance of
        /// returns closest to those predicted by the model.
        /// </summary>
        public virtual void Calibrate(ICalibrationData calibrationData, PriceModel priceModel, XmlWriter output, ErrorList errorList)
        {
            CSForwardPriceModel model = CalibrationHelper.Validate<CSForwardPriceModel>(calibrationData, priceModel, output);

            MomentStatistics[] momentStatistics;
            TimeSeriesStructure<string> timeSeriesStructure = null;
            bool[] isValidStatistics;
            if (!GetMeanAndVolMomentStatistics(calibrationData, priceModel, ReturnMethod.Log, errorList,
                                               ref timeSeriesStructure, out momentStatistics, out isValidStatistics))
                return;

            var validMomentStatistics = momentStatistics.Where((x, i) => isValidStatistics[i]);
            if (!validMomentStatistics.Any())
                return;

            var parameters = validMomentStatistics.Select(x => new TenorAndHistoricalVol(x.ParsePointIDToDouble(), x.Vol)).ToList();

            // perform calibration
            double alpha, sigma;
            GetCalibratedParameters(parameters, out alpha, out sigma);

            model.Alpha = alpha;
            model.Sigma = sigma;

            WriteCalibrationDataToXml(model, output);
        }

        /// <summary>
        /// Get the updated curved for display in calibration
        /// </summary>
        /// <param name="priceModel">CS forward price model</param>
        /// <param name="output">xml writer</param>
        public void GetUpdatedCurves(PriceModel priceModel, XmlWriter output)
        {
            WriteCalibrationDataToXml(priceModel as CSForwardPriceModel, output);
        }

        /// <summary>
        /// Get type of associated model
        /// </summary>
        /// <returns>CS forward price model type</returns>
        public Type ModelType()
        {
            return typeof(CSForwardPriceModel);
        }

        /// <summary>
        /// Merit function for historical calibration of the CS forward price model.
        /// </summary>
        /// <param name="alpha">Alpha; one of the parameters of the CS model</param>
        /// <param name="sigma">Sigma; one of the parameters of the CS model</param>
        /// <param name="parameters">List of tenors and historical volatilities for those tenors.</param>
        /// <param name="meritFunction">Merit function, which should be minimised by varying alpha & sigma in order to best fit the model to the data</param>
        /// <returns>bool to indicate success of calculation</returns>
        private static bool MeritFunction(double alpha, double sigma, List<TenorAndHistoricalVol> parameters, out double meritFunction)
        {
            meritFunction = 0.0;

            foreach (var parameter in parameters)
            {
                double modelVolatility = sigma * Math.Exp(-alpha * parameter.fTenor);
                meritFunction += CalcUtils.Sqr(modelVolatility - parameter.fHistoricalVol);
            }

            return true;
        }

        /// <summary>
        /// Writes calibration data to XML so it can be utilised.
        /// </summary>
        /// <param name="priceModel"></param>
        /// <param name="output"></param>
        private static void WriteCalibrationDataToXml(CSForwardPriceModel priceModel, XmlWriter output)
        {
            ResultSeriesList graphData = new ResultSeriesList();
            graphData.fXAxis.fLabel = "Time in Years";
            graphData.fYAxis.fLabel = "Value";
            graphData.ToXml(output);

            output.WriteStartElement("CSForwardPriceCalibration");
            output.WriteStartElement("Parameters");
            output.WriteStartElement("Alpha");
            output.WriteValue(priceModel.Alpha.ToString(CultureInfo.InvariantCulture));
            output.WriteEndElement();
            output.WriteStartElement("Sigma");
            output.WriteValue(priceModel.Sigma.ToString(CultureInfo.InvariantCulture));
            output.WriteEndElement();
            output.WriteEndElement();
            output.WriteEndElement();
        }

        /// <summary>
        /// Gets the calibrated alpha and sigma values.
        /// </summary>
        /// <param name="parameters">Merit function parameters to be used for the minimization function.</param>
        /// <param name="sigma">A double to which sigma will be asigned.</param>
        /// <param name="alpha">A double to which alpha will be asigned.</param>
        /// <returns>A boolean indicating success.</returns>
        private bool GetCalibratedParameters(List<TenorAndHistoricalVol> parameters, out double alpha, out double sigma)
        {
            // simplex is just a triangle in (alpha, sigma) space
            double[,] simplex = new double[3, 2];

            const double scale = 0.001; //scale chosen for magnitude of real alpha, sigma values

            // first vertex at (scale, scale)
            simplex[0, 0] = scale;
            simplex[0, 1] = scale;
            // second vertex at (2 * scale, scale)
            simplex[1, 0] = 2 * scale;
            simplex[1, 1] = scale;
            // third vertex at (scale, 2 * scale)
            simplex[2, 0] = scale;
            simplex[2, 1] = 2 * scale;

            int valuationCount; // not used currently 
            double minValue; // not used currently - could be used to output how well the calibration worked if required
            return Minimizer<List<TenorAndHistoricalVol>>.MinimizeDownhillSimplex(MeritFunction, parameters, simplex, fNumSims, out alpha, out sigma, out minValue, out valuationCount);
        }
    }
}