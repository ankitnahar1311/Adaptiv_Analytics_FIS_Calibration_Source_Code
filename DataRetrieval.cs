using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing.Design;
using System.Globalization;
using System.Linq;

using SunGard.Adaptiv.Analytics.Framework;

namespace SunGard.Adaptiv.Analytics.Models
{
    /// <summary>
    /// Historical data retrieval parameters (without a calendar and with a horizon) for bootstrapping.
    /// </summary>
    public class BootstapperDataRetrievalParameters : DataRetrievalParametersBase<string>
    {
        private readonly HorizonContainer fHorizon;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public BootstapperDataRetrievalParameters()
        {
            fHorizon = new HorizonContainer();
        }

        /// <summary>
        /// Copy constructor.
        /// </summary>
        public BootstapperDataRetrievalParameters(BootstapperDataRetrievalParameters other) : base(other)
        {
            fHorizon = new HorizonContainer(other.fHorizon);
        }

        /// <summary>
        /// Horizon, i.e., duration of returns. If left blank, we use Frequency for non-overlapping returns.
        /// </summary>
        public string Horizon
        {
            get
            {
                return fHorizon.TermString;
            }
            set
            {
                fHorizon.TermString = value;
            }
        }

        /// <summary>
        /// Returns the either the horizon, or the frequency if no horizon is specified.
        /// </summary>
        public override Term? GetHorizon()
        {
            return fHorizon.TermValue ?? Frequency;
        }
    }

    /// <summary>
    /// Historical data retrieval parameters (with a calendar and with a horizon) for calibration.
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    public class CalibratorDataRetrievalParameters<T> : DataRetrievalParametersWithCalendar<T>
    {
        private readonly HorizonContainer fHorizon;

        /// <summary>
        /// Default constructor.
        /// </summary>
        public CalibratorDataRetrievalParameters()
        {
            fHorizon = new HorizonContainer();
        }

        /// <summary>
        /// Default constructor.
        /// </summary>
        public CalibratorDataRetrievalParameters(CalibratorDataRetrievalParameters<T> other) : base(other)
        {
            fHorizon = new HorizonContainer(other.fHorizon);
        }

        /// <summary>
        /// Horizon, i.e., duration of returns. If left blank, we use Frequency for non-overlapping returns.
        /// </summary>
        public string Horizon
        {
            get
            {
                return fHorizon.TermString;
            }
            set
            {
                fHorizon.TermString = value;
            }
        }

        /// <summary>
        /// Returns the either the horizon, or the frequency if no horizon is specified.
        /// </summary>
        public override Term? GetHorizon()
        {
            return fHorizon.TermValue ?? Frequency;
        }
    }

    /// <summary>
    /// Historical data retrieval parameters (with a calendar and without a horizon).
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    public class DataRetrievalParametersWithCalendar<T> : DataRetrievalParametersBase<T>
    {
        // Holiday calendar used for data retrieval - null at initialisation.
        private IHolidayCalendar fHolidayCalendar;

        // False at initialisation.
        private bool fHolidayCalendarInitialised;

        /// <summary>
        /// Default constructor with pre-calibration methods turned off by default.
        /// </summary>
        public DataRetrievalParametersWithCalendar()
        {
            Calendar = string.Empty;
        }

        /// <summary>
        /// Copy constructor.
        /// </summary>
        public DataRetrievalParametersWithCalendar(DataRetrievalParametersWithCalendar<T> other) : base(other)
        {
            fHolidayCalendarInitialised = false;
            Calendar = other.Calendar;
        }

        /// <summary>
        /// ID of Holiday Calendar.
        /// </summary>
        [LocationStringsFormAttribute]
        [Editor(typeof(ModalEditor), typeof(UITypeEditor))]
        public string Calendar { get; set; }

        /// <summary>
        /// Get the holiday calendar - create one if uninitialized.
        /// </summary>
        public IHolidayCalendar GetHolidayCalendar(ICalibrationData calibrationData, ErrorList errorList)
        {
            if (fHolidayCalendarInitialised == false)
            {
                fHolidayCalendar = CalibrationHelper.GetHolidayCalendar(calibrationData, errorList, Calendar);
                fHolidayCalendarInitialised = true;
            }

            return fHolidayCalendar;
        }
    }

    /// <summary>
    /// Historical data retrieval parameters (without a calendar and without a horizon).
    /// </summary>
    /// <typeparam name="T">T is the type of id of time series, e.g., Period, string or Vol.Point. </typeparam>
    public abstract class DataRetrievalParametersBase<T> : NestedPresentableObject
    {
        // Data cleaning method - null at initialisation.
        private List<IDataCleaningMethod<T>> fDataCleaningMethodList;

        private int? fBusinessDaysInYear = CalcUtils.BUSINESS_DAYS_IN_YEAR;
        private Term fFrequency = new Term(1);
        private Term? fSamplePeriod;

        /// <summary>
        /// Default constructor with pre-calibration methods turned off by default.
        /// </summary>
        protected DataRetrievalParametersBase()
        {
            Business_Days_In_Year   = "260";
            Data_Cleaning_Methods   = new DisplayableList<DataCleaningMethodWrapper<T>>();
            Diagnostics_Error_Level = ErrorLevel.Info;
            Sample_Period = string.Empty;
        }

        /// <summary>
        /// Default constructor with pre-calibration methods turned off by default.
        /// </summary>
        protected DataRetrievalParametersBase(DataRetrievalParametersBase<T> other)
        {
            Start_Date              = other.Start_Date;
            End_Date                = other.End_Date;
            Sample_Period           = other.Sample_Period;
            Frequency               = other.Frequency;
            Business_Days_In_Year   = other.Business_Days_In_Year;
            Diagnostics_Error_Level = other.Diagnostics_Error_Level;
            Data_Cleaning_Methods   = new DisplayableList<DataCleaningMethodWrapper<T>>(other.Data_Cleaning_Methods);
        }

        /// <summary>
        /// Uder-defined start date for calibration. Optional.
        /// </summary>
        public TDate Start_Date { get; set; }

        /// <summary>
        /// Uder-defined start date for calibration. Optional.
        /// </summary>
        public TDate End_Date { get; set; }

        /// <summary>
        /// Length of market data.
        /// </summary>
        public string Sample_Period
        {
            get
            {
                return fSamplePeriod.HasValue ? fSamplePeriod.Value.ToString() : string.Empty;
            }
            set
            {
                if (string.IsNullOrWhiteSpace(value))
                    fSamplePeriod = null;
                else
                    fSamplePeriod = new Term(value);
            }
        }

        /// <summary>
        /// Frequency, i.e., period between observations.
        /// </summary>
        public Term Frequency
        {
            get
            {
                return fFrequency;
            }
            set
            {
                if (!value.IsZero())
                    fFrequency = value;
            }
        }

        /// <summary>
        /// Number of businsess days in year
        /// </summary>
        /// <remarks>
        /// Will be calculated from historical data if empty.
        /// </remarks>
        public string Business_Days_In_Year
        {
            get
            {
                return fBusinessDaysInYear.HasValue ? fBusinessDaysInYear.Value.ToString(CultureInfo.InvariantCulture) : string.Empty;
            }
            set
            {
                int d;
                if (!string.IsNullOrEmpty(value) && int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out d) && d > 0)
                    fBusinessDaysInYear = d;
                else
                    fBusinessDaysInYear = null;
            }
        }

        /// <summary>
        /// Gets or sets whether data retrieval diagnostics level.
        /// </summary>
        public ErrorLevel Diagnostics_Error_Level { get; set; }

        /// <summary>
        /// List of data cleaning methods with their parameters.
        /// </summary>
        public DisplayableList<DataCleaningMethodWrapper<T>> Data_Cleaning_Methods { get; set; }

        /// <summary>
        /// Returns the sample period of the market data.
        /// </summary>
        public Term? GetSamplePeriod()
        {
            return fSamplePeriod;
        }

        /// <summary>
        /// Returns the horizon.
        /// </summary>
        public virtual Term? GetHorizon()
        {
            return null;
        }

        /// <summary>
        /// Returns business days in year.
        /// </summary>
        public int? GetBusinessDaysInYear()
        {
            return fBusinessDaysInYear;
        }

        /// <summary>
        /// Returns the start date. An empty field is given as a TDate 0. We transform it to null then.
        /// </summary>
        public TDate? GetStartDate()
        {
            return Start_Date == 0 ? (TDate?)null : Start_Date;
        }

        /// <summary>
        /// Returns the start date. An empty field is given as a TDate 0. We transform it to null then.
        /// </summary>
        public TDate? GetEndDate()
        {
            return End_Date == 0 ? (TDate?)null : End_Date;
        }

        /// <summary>
        /// Returns a list of data cleaning methods.
        /// </summary>
        public List<IDataCleaningMethod<T>> GetDataCleaningMethods()
        {
            if (fDataCleaningMethodList == null)
            {
                fDataCleaningMethodList = Data_Cleaning_Methods.Select(m => m.Parameters).ToList();
            }

            return fDataCleaningMethodList;
        }

        /// <summary>
        /// Properties validation. Provides user with a list of errors and warnings about all properties.
        /// </summary>
        public bool ValidateProperties(ErrorList errorList)
        {
            bool validProperties = true;

            // Validate data cleaning parameters.
            foreach (var dataCleaningMethod in Data_Cleaning_Methods)
            {
                validProperties = validProperties && dataCleaningMethod.Parameters.ValidateParameters(errorList);
            }

            if (Start_Date > 0.0 && End_Date > 0.0 && Start_Date > End_Date)
            {
                errorList.Add(ErrorLevel.Error, "The calibration end date must be on or after the calibration start date.");
                validProperties = false;
            }

            return validProperties;
        }
    }

    /// <summary>
    /// Container for the horizon.
    /// </summary>
    internal class HorizonContainer
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public HorizonContainer()
        {
            TermString = string.Empty;
        }

        /// <summary>
        /// Copy constructor.
        /// </summary>
        public HorizonContainer(HorizonContainer other)
        {
            TermString = other.TermString;
        }

        /// <summary>
        /// Horizon, i.e., duration of returns. If left blank, we use Frequency for non-overlapping returns.
        /// </summary>
        public string TermString
        {
            get
            {
                return TermValue.HasValue ? TermValue.Value.ToString() : string.Empty;
            }
            set
            {
                TermValue = string.IsNullOrEmpty(value) ? (Term?)null : new Term(value);
            }
        }

        /// <summary>
        /// Return the horizon.
        /// </summary>
        public Term? TermValue { get; private set; }
    }
}