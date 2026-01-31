// <summary>
// Clean data retrieval with pre-calibration methods (for data cleaning and filling).
// </summary>
using System;
using System.Collections.Generic;
using System.Text;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Static data retrieval routines for calibrations use.
    /// </summary>
    public static class GenericDataRetrieval
    {
        /// <summary>
        /// Return the historical times. The calendar is extracted from calibration data given an explicit Calendar ID in data retrieval parameters.
        /// </summary>
        /// <remarks>
        /// We select all points for the given factor type and ID.
        /// </remarks>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive and the calendar.</param>
        /// <param name="factorTypeAndID">Price factors to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="idConverter">Function to convert archive ID to point.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="timeSeriesStructure">Output a time series structure with generic time series string ids.</param>
        /// <returns>Return true if timeseries structure is successfully extracted.</returns>
        /// <typeparam name="TPoint">Point type.</typeparam>
        public static bool TryGetHistoricalTimeSeriesStructure<TPoint>(ICalibrationData calibrationData,
            FactorTypeAndID factorTypeAndID, DataRetrievalParametersWithCalendar<TPoint> dataRetrievalParameters,
            Func<string, TPoint> idConverter, ErrorList errorList, out TimeSeriesStructure<TPoint> timeSeriesStructure)
        {
            var calendar = dataRetrievalParameters.GetHolidayCalendar(calibrationData, errorList);
            return TryGetHistoricalTimeSeriesStructure(calibrationData, calendar, factorTypeAndID,
                dataRetrievalParameters, s => true, idConverter, errorList, out timeSeriesStructure);
        }

        /// <summary>
        /// Return the historical times. We default to archive calendar given no explicit calendar.
        /// </summary>
        /// <param name="historicalData">IHistoricalDataRetrieval object containing historical data archive.</param>
        /// <param name="factorTypeAndID">Price factors to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="idConverter">Function to convert archive ID to point.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="timeSeriesStructure">Output a time series structure with generic time series string ids.</param>
        /// <returns>Return true if timeseries structure is successfully extracted.</returns>
        /// <typeparam name="TPoint">Point type.</typeparam>
        public static bool TryGetHistoricalTimeSeriesStructure<TPoint>(IHistoricalDataRetrieval historicalData,
            FactorTypeAndID factorTypeAndID, DataRetrievalParametersBase<TPoint> dataRetrievalParameters, Predicate<ExtendedArchiveItemID> pointSelector,
            Func<string, TPoint> idConverter, ErrorList errorList, out TimeSeriesStructure<TPoint> timeSeriesStructure)
        {
            var calendar = historicalData.GetArchiveCalendar();
            return TryGetHistoricalTimeSeriesStructure(historicalData, calendar, factorTypeAndID,
                dataRetrievalParameters, pointSelector, idConverter, errorList, out timeSeriesStructure);
        }

        /// <summary>
        /// Return the historical times, and the different time series associated with the list of factors at those times. 
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="factorTypeAndID">Price factors to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="timeSeriesStructure">Output a time series structure with generic time series string ids.</param>
        /// <returns>Return true if timeseries structure is successfully extracted.</returns>
        public static bool TryGetHistoricalTimeSeriesStructure(ICalibrationData calibrationData, FactorTypeAndID factorTypeAndID,
            DataRetrievalParametersWithCalendar<string> dataRetrievalParameters, ErrorList errorList,
            out TimeSeriesStructure<string> timeSeriesStructure)
        {
            return TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID, dataRetrievalParameters,
                pointString => pointString,  errorList, out timeSeriesStructure);
        }

        /// <summary>
        /// Data retrieval for a set of time series with string ids, e.g., { EquityPrice.BT, EquityPrice.VOD, EquityPrice.GE }.
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="priceFactor">Price factor to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="timeSeriesStructure">Output a time series structure for a different assets.</param>
        /// <returns>Return true if term structure is successfully extracted.</returns>
        public static bool TryGetHistoricalTimeSeriesStructure(ICalibrationData calibrationData, IPriceFactor priceFactor,
            DataRetrievalParametersWithCalendar<string> dataRetrievalParameters, ErrorList errorList,
            out TimeSeriesStructure<string> timeSeriesStructure)
        {
            var factorTypeAndID = new FactorTypeAndID(priceFactor.TypeDotSubType(), priceFactor.GetID());
            return TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID, dataRetrievalParameters, errorList,
                out timeSeriesStructure);
        }

        /// <summary>
        /// Return the historical times, and the corresponding zero rates over all tenors at those times. 
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="priceFactor">Price factor to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="termStructure">Output a term structure over a set of tenors.</param>
        /// <returns>Return true if term structure is successfully extracted.</returns>
        public static bool TryGetHistoricalRateTermStructure(ICalibrationData calibrationData, IPriceFactor priceFactor,
            DataRetrievalParametersWithCalendar<Period> dataRetrievalParameters, ErrorList errorList,
            out TimeSeriesStructure<Period> timeSeriesStructure)
        {
            var factorTypeAndID = new FactorTypeAndID(priceFactor.TypeDotSubType(), priceFactor.GetID());
            return TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID, dataRetrievalParameters,
                        pointString => CalcUtils.ParsePoint(pointString), errorList, out timeSeriesStructure);
        }

        /// <summary>
        /// Return the historical times, and the corresponding zero rates over all points for a Vol Surface, e.g., EquityPriceVol.
        /// </summary>
        /// <param name="calibrationData">ICalibrationData object containing historical data archive.</param>
        /// <param name="priceFactor">Single asset to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="volStructure">Output a volatilty surface structure.</param>
        /// <returns>Return true if time series is successfully extracted.</returns>
        /// <typeparam name="TVol">TVol can be any VolSurface, e.g., AssetPriceVolSurface, ForwardPriceVolSurface, InterestRateVolSurface and PriceIndexVolSurface. </typeparam>
        public static bool TryGetHistoricalVolSurface<TVol>(ICalibrationData calibrationData, IPriceFactor priceFactor,
            DataRetrievalParametersWithCalendar<VolSurface.Point> dataRetrievalParameters, ErrorList errorList,
            out TimeSeriesStructure<VolSurface.Point> timeSeriesStructure) where TVol : VolSurface, new()
        {
            var factorTypeAndID = new FactorTypeAndID(priceFactor.TypeDotSubType(), priceFactor.GetID());

            // Dummy vol surface needed for the idConverter
            VolSurface volSurface = new TVol();
            
            return TryGetHistoricalTimeSeriesStructure(calibrationData, factorTypeAndID,
                dataRetrievalParameters, pointString => volSurface.PointFromString(pointString), errorList,
                out timeSeriesStructure);
        }
  
        /// <summary>
        /// Return the historical times, and the different time series associated with the list of factors at those times. 
        /// </summary>
        /// <param name="historicalData">IHistoricalDataRetrieval object containing historical data archive.</param>
        /// <param name="calendar">Calendar used to calculate returns.</param>
        /// <param name="factorTypeAndID">Price factors to extract.</param>
        /// <param name="dataRetrievalParameters"><see cref="DataRetrievalParametersWithCalendar{T}">Data retrieval parameters</see> for more information.</param>
        /// <param name="idConverter">Function to convert archive ID to point.</param>
        /// <param name="errorList">List of infos, warnings and errors. Can be null to turn off diagnostics.</param>
        /// <param name="timeSeriesStructure">Output a time series structure with generic time series string ids.</param>
        /// <returns>Return true if timeseries structure is successfully extracted.</returns>
        /// <typeparam name="TPoint">Point type.</typeparam>
        internal static bool TryGetHistoricalTimeSeriesStructure<TPoint>(IHistoricalDataRetrieval historicalData, IHolidayCalendar calendar,
            FactorTypeAndID factorTypeAndID, DataRetrievalParametersBase<TPoint> dataRetrievalParameters, Predicate<ExtendedArchiveItemID> pointSelector,
            Func<string, TPoint> idConverter, ErrorList errorList, out TimeSeriesStructure<TPoint> timeSeriesStructure)
        {
            Requires.Assert(historicalData != null);
            Requires.Assert(factorTypeAndID != null);
            Requires.Assert(dataRetrievalParameters != null);

            // Extract market data
            TDate[] archiveDates;
            double[][] archiveData;
            List<ExtendedArchiveItemID> archiveItems;
            bool isSuccess = TryExtractMarketData(historicalData, factorTypeAndID, dataRetrievalParameters, pointSelector, out archiveItems,
                out archiveDates, out archiveData);

            var priceFactorName = factorTypeAndID.fFactorType + FactorID.CodeSeparator + factorTypeAndID.fFactorID.ToCode();

            if (!isSuccess)
            {
                if (errorList != null)
                    errorList.Add(ErrorLevel.Error, string.Format("Historical data not available for {0}.", priceFactorName));

                timeSeriesStructure = null;
                return false;
            }

            // Extract time series points.
            int numTimeSeries = archiveData.Length;
            var timeSeriesIDs = new TPoint[numTimeSeries];
            for (int timeSeriesIndex = 0; timeSeriesIndex < numTimeSeries; ++timeSeriesIndex)
            {
                timeSeriesIDs[timeSeriesIndex] = idConverter(archiveItems[timeSeriesIndex].ItemID.fPointID);
            }

            // Create clean term structure.
            timeSeriesStructure = new TimeSeriesStructure<TPoint>(archiveDates, timeSeriesIDs, archiveData, calendar, dataRetrievalParameters, errorList);

            if (timeSeriesStructure.Length < StatisticsCalculationHelper.MinTimeSeriesLength)
            {
                if (errorList != null)
                {
                    errorList.Add(ErrorLevel.Error,
                        string.Format("Not enough historical data for {0}; there must be at least 3 observations.",
                        priceFactorName));
                }

                return false;
            }

            // Data retrieval diagnostic.
            if (errorList != null)
                DataRetrievalDiagnostics(errorList, dataRetrievalParameters.Diagnostics_Error_Level, priceFactorName, timeSeriesStructure);

            return true;
        }

        /// <summary>
        /// Extract market data give a list of price factors.
        /// </summary>
        /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
        internal static bool TryExtractMarketData<T>(IHistoricalDataRetrieval calibrationData, 
            FactorTypeAndID factorTypeAndID, DataRetrievalParametersBase<T> dataRetrievalParameters, 
            Predicate<ExtendedArchiveItemID> pointSelector, out List<ExtendedArchiveItemID> archiveItems, 
            out TDate[] historicalDates, out double[][] historicalValues)
        {
            Requires.Assert(calibrationData != null);

            var lengthData = dataRetrievalParameters.GetSamplePeriod();

            var userStartDate = dataRetrievalParameters.GetStartDate().HasValue ? dataRetrievalParameters.GetStartDate().Value : 0;
            var userEndDate = dataRetrievalParameters.GetEndDate().HasValue ? dataRetrievalParameters.GetEndDate().Value : 0;

            var dates = StatisticsSettings.GenerateDatesWithinBounds(userStartDate, userEndDate, lengthData, calibrationData.ArchiveStartDate, calibrationData.ArchiveEndDate);

            var startDate = dates.StartDate;
            var endDate = dates.EndDate;

            // Determine the number of tenors stored in the archive.
            IHistoricalArchive archive = calibrationData.GetHistoricalMarketRates(new FactorTypeAndIDList { factorTypeAndID },
                startDate, endDate);

            ThrowIfArchiveDatesInvalid(archive);

            var allLocalArchiveItems = CleanArchiveGenerator.GetArchiveItems(archive, factorTypeAndID);
            var selectedArchiveItems = allLocalArchiveItems.FindAll(pointSelector);
            var numTimeSeries = selectedArchiveItems.Count;
            if (numTimeSeries == 0)
            {
                archiveItems = null;
                historicalDates = null;
                historicalValues = null;
                return false;
            }

            var historicalValuesList = new List<double[]>();
            var historicalDatesList  = new List<TDate>();

            // Iterate over the historical dates in the archive, and store each date and the zero rates for all tenors at that date.
            foreach (DateTime date in archive.HistoricalDates)
            {
                // Start date for data retrieval.
                if (date < startDate)
                    continue;

                // End date for data retrieval.
                if (date > endDate)
                    break;

                int dateIndex = archive.IndexOfDate(date);

                // Get zero-rates across all tenors.
                var values = new double[numTimeSeries];
                for (int tenorIndex = 0; tenorIndex < numTimeSeries; ++tenorIndex)
                {
                    values[tenorIndex] = archive[selectedArchiveItems[tenorIndex].Index, dateIndex];
                }

                historicalValuesList.Add(values);
                historicalDatesList.Add(date);
            }

            // Validate market data
            if (historicalDatesList.Count < StatisticsCalculationHelper.MinTimeSeriesLength)
            {
                archiveItems = null;
                historicalDates = null;
                historicalValues = null;
                return false;
            }

            archiveItems = selectedArchiveItems;
            historicalDates = historicalDatesList.ToArray();
            historicalValues = CleanArchiveGenerator.GetHistoricalValuesArray(historicalValuesList);
            return true;
        }

        /// <summary>
        /// Throws if the archive doesn't contain any data or if the archive dates are not in the strictly increasing order.
        /// </summary>
        internal static void ThrowIfArchiveDatesInvalid(IHistoricalArchive historicalArchive)
        {
            var archiveDates = historicalArchive.HistoricalDates;
            if (archiveDates.Length == 0)
                throw new CalibrationException("The archive doesn't contain any data.");

            DateTime previousArchiveDate = DateTime.MinValue;
            foreach (DateTime date in archiveDates)
            {
                if (date <= previousArchiveDate)
                    throw new CalibrationException("The archive dates must be stricly increasing.");

                previousArchiveDate = date;
            }
        }

        /// <summary>
        /// Data retrieval diagnostics.
        /// </summary>
        /// <param name="errorList">Error list to output info, warning and error messages.</param>
        /// <param name="userDefinedMinErrorLevel">A minimum error level that user wants to see for the diagnostics.</param>
        /// <param name="priceFactorName">Price factor name.</param>
        /// <param name="timeSeriesStructure">Time series structure.</param>
        /// <typeparam name="T">T is the type of time series ids. </typeparam>
        private static void DataRetrievalDiagnostics<T>(ErrorList errorList, ErrorLevel userDefinedMinErrorLevel,
            string priceFactorName, TimeSeriesStructure<T> timeSeriesStructure)
        {
            if (userDefinedMinErrorLevel == ErrorLevel.None || ErrorLevel.Info < userDefinedMinErrorLevel)
                return;

            errorList.Add(ErrorLevel.Info, string.Format("Data retrieval for {0} complete with {1} sampled dates.",
                                            priceFactorName, timeSeriesStructure.Length));

            var firstDate = timeSeriesStructure.Dates[0];
            var lastDate = timeSeriesStructure.Dates[timeSeriesStructure.Length - 1];
            Term term;
            Term.DateDiff(firstDate, lastDate, out term);

            int missingPercentage   = (int)Math.Round(100.0 * timeSeriesStructure.RatioMissingData);

            var stringBuilder       = new StringBuilder(string.Format("Historical data over {0} extracted to calibrate the model", term));

            if (missingPercentage > 0)
                stringBuilder.Append(string.Format(" and {0}% of the data during that period missing", missingPercentage));

            stringBuilder.Append(".");

            errorList.Add(ErrorLevel.Info, stringBuilder.ToString());
        }
    }
}
