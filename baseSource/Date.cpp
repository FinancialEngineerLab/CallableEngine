#define MONTH_JAN "JAN"
#define MONTH_FEB "FEB"
#define MONTH_MAR "MAR"
#define MONTH_APR "APR"
#define MONTH_MAY "MAY"
#define MONTH_JUN "JUN"
#define MONTH_JUL "JUL"
#define MONTH_AUG "AUG"
#define MONTH_SEP "SEP"
#define MONTH_OCT "OCT"
#define MONTH_NOV "NOV"
#define MONTH_DEC "DEC"

#define SUN 0
#define MON 1
#define TUE 2
#define WED 3
#define THU 4
#define FRI 5
#define SAT 6

#include "Date.h"


CDate::CDate()
{
}

CDate::CDate(int Y, int M, int D)
{
	y = Y;
	m = M;
	d = D;

	ymd = vector<int>(3);
	
	ymd[0] = y;
	ymd[1] = m;
	ymd[2] = d;
	
	yyyymmdd = set_yyyymmdd(ymd);

	totaldays = get_totaldays(y, m, d);

	t = double(totaldays);
}

CDate::CDate(vector<int> YMD)
{
	ymd = YMD;

	y = ymd[0];
	m = ymd[1];
	d = ymd[2];

	yyyymmdd = set_yyyymmdd(ymd);

	totaldays = get_totaldays(y, m, d);

	t = double(totaldays);
}

CDate::CDate(string YYYYMMDD)
{
	yyyymmdd = YYYYMMDD;
	ymd = get_yyyymmdd(yyyymmdd);
	
	y = ymd[0];
	m = ymd[1];
	d = ymd[2];

	totaldays = get_totaldays(y, m, d);

	t = double(totaldays);
}

CDate::CDate(int TOTALDAYS)
{
	totaldays = TOTALDAYS;
	t = double(totaldays);
	ymd = totaldays2YYYYMMDD(totaldays);

	y = ymd[0];
	m = ymd[1];
	d = ymd[2];

	yyyymmdd = set_yyyymmdd(ymd);
}

string CDate::get_str_date()
{
	return yyyymmdd;
}

vector<int> CDate::get_vec_date()
{
	return ymd;
}

int CDate::get_int_date()
{
	return totaldays;
}

int CDate::get_year()
{
	return y;
}

int CDate::get_month()
{
	return m;
}

int CDate::get_day()
{
	return d;
}

void CDate::set_year(int year)
{
	y = year;
}

void CDate::set_month(int month)
{
	m = month;
}

void CDate::set_day(int day)
{
	d = day;
}

double CDate::get_dbl_date()
{
	return t;
}

int CDate::operator[](int i)
{
	if (i==0)
	{
		return y;
	}
	else if (i==1)
	{
		return m;
	}
	else
	{
		return d;
	}
}

bool CDate::operator==(const CDate& date)
{
	return (yyyymmdd == date.yyyymmdd);
}

bool CDate::operator!=(const CDate& date)
{
	return (yyyymmdd != date.yyyymmdd);
}

bool CDate::operator<(const CDate& date)
{
	return (t < date.t);
}

bool CDate::operator>(const CDate& date)
{
	return (t > date.t);
}

bool CDate::operator<=(const CDate& date)
{
	return (t <= date.t);
}

bool CDate::operator>=(const CDate& date)
{
	return (t >= date.t);
}


double CDate::operator-(const CDate& date) const
{
	return t - date.t;
}

CDate::~CDate()
{
}

double cvg(CDate t_st, CDate t_end, string dcb)
{
	double ans = 0.0;
	
	if (dcb == "30/365")
	{
		ans = (fmax(30.0 - t_st[2], 0.0) + fmin(t_end[2], 30) + 360.0 * (t_end[0] - t_st[0]) + 30.0 * (t_end[1]-t_st[1] - 1)) / 360.0;
	}
	else if (dcb == "ACT/360")
	{
		ans = (t_end - t_st) / 360.0;
	}
	else if (dcb == "30/365NL")
	{
		int i, count = 0, y1 = t_st.get_year(), y2 = t_end.get_year();

		for (i=y1; i<=y2; i++)
		{
			if (is_leapyear(i))
			{
				CDate date(i, 2, 29);

				//if (t_st < date && date < t_end)
				//{
				//	count = count + 1;
				//}

				if (t_st < date && date < t_end)
				{
					count = count + 1;
				}
			}
		}

		ans = (t_end - t_st - count) / 365.0;
	}
	else if(dcb == "ACT/ACTKRW")
	{
		int y1 = t_st.get_year(), y2 = t_end.get_year();
		double den = 0.0;

		if (y1==y2)
		{
			den = (is_leapyear(y1) ? 366.0 : 365.0);
		}
		else
		{
			CDate date1(y1, 2, 28), date2(y2, 2, 28);
			den = ((is_leapyear(y1) && t_st <= date1) || (is_leapyear(y2) && t_end >= date2) ? 366.0 : 365.0);
		}
		ans = (t_end - t_st) / den;

	}
	else
	{
		ans = (t_end - t_st) / 365.0;
	}

	return ans;
}

int is_leapyear(int year)
{
	return !(year % 400) || (year % 100) && !(year % 4);
}

int get_monthdays(int year, int month)
{
	static int mdays[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int md = mdays[month - 1];
	if (is_leapyear(year) && month ==2)
	{
		md++;
	}
	return md;
}

int get_totaldays(int year, int month, int day)
{
	int y1 = year - 1;
	int total = 365 * y1;
	
	total += (y1 / 400 - y1 / 100 + y1 / 4);

	for (int m = 1; m < month; m++)
	{
		total += get_monthdays(year, m);
	}
	total += day;
	return total;
}

string set_yyyymmdd(vector<int> YYYYMMDD)
{
	char yyyy[65], mm[65], dd[65];

	string yyyymmdd;

	_itoa_s(YYYYMMDD[0], yyyy, 65, 10);
	_itoa_s(YYYYMMDD[1], mm, 65, 10);
	_itoa_s(YYYYMMDD[2], dd, 65, 10);

	yyyymmdd.append(yyyy);
	yyyymmdd.append("-");

	if (YYYYMMDD[1]<10)
	{
		yyyymmdd.append("0");
	}

	yyyymmdd.append(mm);
	yyyymmdd.append("-");

	if (YYYYMMDD[2] < 10)
	{
		yyyymmdd.append("0");
	}

	yyyymmdd.append(dd);

	return yyyymmdd;

}

vector<int> get_yyyymmdd(string yyyymmdd)
{
	char yyyy[4], mm[2], dd[2];

	int ny, nm, nd;

	vector<int> YYYYMMDD(3);

	ny = int(yyyymmdd._Copy_s(yyyy, 4, 4));
	nm = int(yyyymmdd._Copy_s(mm, 2,2,5));
	nd = int(yyyymmdd._Copy_s(dd, 2,2,8));

	YYYYMMDD[0] = atoi(yyyy);
	YYYYMMDD[1] = atoi(mm);
	YYYYMMDD[2] = atoi(dd);

	return YYYYMMDD; 

}

vector<int> totaldays2YYYYMMDD(int totaldays)
{
	vector<int> YYYYMMDD(3);

	int tempdays = totaldays, YYYY = 1, MM = 1, DD, leapflag = is_leapyear(YYYY);

	for ( ; tempdays > 365 + leapflag; )
	{
		tempdays = tempdays - 365 - leapflag;
		YYYY = YYYY + 1;
		leapflag = is_leapyear(YYYY);
	}

	YYYYMMDD[0] = YYYY;
	int mdays[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	
	if (is_leapyear(YYYY))
	{
		mdays[1] = mdays[1] + 1;
	}

	for (;tempdays > mdays[MM-1];)
	{
		tempdays = tempdays - mdays[MM - 1];
		MM = MM + 1;
	}

	YYYYMMDD[1] = MM;
	DD = tempdays;
	YYYYMMDD[2] = DD;

	return YYYYMMDD;
}

int get_weekday(CDate date)
{
	return (get_totaldays(date[0], date[1], date[2])) % 7;
}

CDate DateAdd(char ymd, int numymd, CDate date)
{
	CDate newdate = date;

	if (ymd == 'y')
	{
		newdate.set_year(newdate[0] + numymd);
		if (date[1] == 2 && date[2] == 29)
		{
			if (!(is_leapyear(newdate[0])))
			{
				newdate.set_day(28);
			}
		}
		newdate = CDate(date[0], date[1], date[2]);
	}
	else if (ymd == 'm')
	{
		int adj = (sign(sign(date[1] + numymd) - 1) - 1) / 2;
		int q, d;

		q = (date[1]+numymd) / 12;

		newdate.set_month(date[1] + numymd - (q + adj) * 12);

		if (newdate[1]<1)
		{
			q = q - 1;
		}

		newdate.set_month(date[1] + numymd - (q + adj) * 12);
		newdate.set_year(newdate[0]+q+adj);
		
		d = get_monthdays(newdate[0], newdate[1]);

		if (date[2]>d)
		{
			newdate.set_day(d);
		}

		newdate = CDate(newdate[0], newdate[1], newdate[2]);


	}
	else
	{
		newdate = CDate(get_totaldays(date[0], date[1], date[2]) + numymd);
		//newdate = get_totaldays(date[0], date[1], date[2]);

	}

	return newdate;
}

int get_month_from_MMMYY(string str)
{
	basic_string<char>::iterator str_Iter;
	str_Iter = str.erase(str.end() - 2, str.end());

	if (str==MONTH_JAN)
	{
		return 1;
	}
	else if (str == MONTH_FEB)
	{
		return 2;
	}
	else if (str == MONTH_MAR)
	{
		return 3;
	}
	else if (str == MONTH_APR)
	{
		return 4;
	}
	else if (str == MONTH_MAY)
	{
		return 5;
	}
	else if (str == MONTH_JUN)
	{
		return 6;
	}
	else if (str == MONTH_JUL)
	{
		return 7;
	}
	else if (str == MONTH_AUG)
	{
		return 8;
	}
	else if (str == MONTH_SEP)
	{
		return 9;
	}
	else if (str == MONTH_OCT)
	{
		return 10;
	}
	else if (str == MONTH_NOV)
	{
		return 11;
	}
	else
	{
		return 12;
	}
}

int get_year_from_MMMYY(string str)
{
	basic_string<char>::iterator str_Iter;
	str_Iter = str.erase(str.begin(), str.begin() + 3);

	int yy = atoi(&(str[0]));

	return 2000 + yy;

}














