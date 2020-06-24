#include "Calendar.h"
#include "Date.h"
#include "ConstantNumbers.h"

CCalendar::CCalendar()
{
}
CCalendar::~CCalendar()
{
}
CCalendar::CCalendar(string HOLIDAYSNAME)
{
	holidaysname = HOLIDAYSNAME;

	static string holidaysset[] = { "SEB","NYB","LNB","LNNYB","TKB","LNTKB","TGT","CME","CMELNB","NYSEB" };

	static int numholidaysset[] = { NumSEHolidays, NumNYHolidays, NumLNHolidays, NumLNNYHolidays, NumTKHolidays, NumLNTKHolidays, NumTARGETHolidays, NumCMEHolidays, NumCMELNHolidays, NumNYSEHolidays };

	int i = 0;

	if (holidaysname == "NYB")
	{
		numholidays = numholidaysset[1];

	}
	else if (holidaysname == "LNB")
	{
		numholidays = numholidaysset[2];

	}
	else if (holidaysname == "LNNYB")
	{
		numholidays = numholidaysset[3];

	}
	else if (holidaysname == "TKB")
	{
		numholidays = numholidaysset[4];

	}
	else if (holidaysname == "LNTKB")
	{
		numholidays = numholidaysset[5];

	}
	else if (holidaysname == "TGT")
	{
		numholidays = numholidaysset[6];

	}
	else if (holidaysname == "CME")
	{
		numholidays = numholidaysset[7];

	}
	else if (holidaysname == "CMELNB")
	{
		numholidays = numholidaysset[8];

	}
	else if (holidaysname == "NYSEB")
	{
		numholidays = numholidaysset[9];

	}
	else
	{
		numholidays = numholidaysset[0];

	}

	if (holidaysname == "NYB")
	{
		holidays = vector<CDate>(NYHolidays, NYHolidays + numholidays);
	}
	else if (holidaysname == "LNB")
	{
		holidays = vector<CDate>(LNHolidays, LNHolidays + numholidays);
	}
	else if (holidaysname == "LNNYB")
	{
		holidays = vector<CDate>(LNNYHolidays, LNNYHolidays + numholidays);
	}
	else if (holidaysname == "TKB")
	{
		holidays = vector<CDate>(TKHolidays, TKHolidays + numholidays);
	}
	else if (holidaysname == "LNTKB")
	{
		holidays = vector<CDate>(LNTKHolidays, LNTKHolidays + numholidays);
	}
	else if (holidaysname == "TGT")
	{
		holidays = vector<CDate>(TGTHolidays, TGTHolidays + numholidays);
	}
	else if (holidaysname == "CME")
	{
		holidays = vector<CDate>(CMEHolidays, CMEHolidays + numholidays);
	}
	else if (holidaysname == "CMELNB")
	{
		holidays = vector<CDate>(CMELNHolidays, CMELNHolidays + numholidays);
	}
	else if (holidaysname == "NYSEB")
	{
		holidays = vector<CDate>(NYSEHolidays, NYSEHolidays + numholidays);
	}
	else
	{
		holidays = vector<CDate>(SEHolidays, SEHolidays + numholidays);
	}
}

CDate CCalendar::operator[](int i)
{
	return holidays[i];
}

int CCalendar::get_numholidays()
{
	return numholidays;
}

bool is_holiday(CDate date, CCalendar holidays)
{
	int i, weekdayflag = get_weekday(date);
	bool flag = false;

	if (weekdayflag == 0 || weekdayflag == 6)
	{
		flag = true;
	}
	else
	{
		for (i=0;i<holidays.get_numholidays();i++)
		{
			if (date == holidays[i])
			{
				flag = true;
				break;
			}
		}
	}

	return flag;
}

void ShiftBusDate(CDate date, CCalendar holidays, int numshift, CDate& newdate)
{
	newdate = date;
	int sgn = sign(numshift);
	if (sgn!=0)
	{
		int i;
		for (i=0; i<abs(numshift); i++)
		{
			newdate = DateAdd('d', sgn, newdate);
			while (is_holiday(newdate, holidays))
			{
				newdate = DateAdd('d', sgn, newdate);
			}
		}
	}
}

void ConvReflDate(CDate date, CCalendar holidays, string conv, CDate& newdate)
{
	newdate = date;
	if (conv =="MODFOL")
	{
		while (is_holiday(newdate,holidays))
		{
			newdate = DateAdd('d', 1, newdate);
		}
		if (newdate[1]!=date[1])
		{
			newdate = DateAdd('d', -1, newdate);
			while (is_holiday(newdate, holidays))
			{
				newdate = DateAdd('d', -1, newdate);
			}
		}
	}
	else if (conv == "FOL")
	{
		while (is_holiday(newdate, holidays))
		{
			newdate = DateAdd('d', 1, newdate);
		}
	}
	else if (conv == "MODPREC")
	{
		while (is_holiday(newdate, holidays))
		{
			newdate = DateAdd('d', -1, newdate);
		}
		if (newdate[1] != date[1])
		{
			newdate = DateAdd('d', +1, newdate);
			while (is_holiday(newdate, holidays))
			{
				newdate = DateAdd('d', +1, newdate);
			}
		}
	}
	else if (conv == "PREC")
	{
		while (is_holiday(newdate, holidays))
		{
			newdate = DateAdd('d', -1, newdate);
		}

	}
}


