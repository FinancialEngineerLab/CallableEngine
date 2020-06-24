#include "Schedule.h"
#include "Date.h"

void paymentschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate *holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	vector<CDate>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, num_schedule, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		tempdate = startdate;
		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m', (i - 1) * frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);

		for (i=0;i<num_schedule-1;i++)
		{
			ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[num_schedule - 1]);
		}
	}
	else
	{
		tempdate = maturitydate;
		while (tempdate > startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);

		for (i = 0; i < num_schedule ; i++)
		{
			ConvReflDate(DateAdd('m',(i-num_schedule+1)*frq, maturitydate),holidays, numholiday, conv, schedule[i]);
		}
	}
};

void paymentschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	paymentschedule(startdate, maturitydate,	holidays, numholiday, stub, direction, conv, freq, schedule);
};

void calculationstartschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<CDate> &schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, num_schedule, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		tempdate = startdate;
		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m', (i - 1) * frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);
		schedule[0] = startdate;

		for (i = 1; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', i * frq, startdate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
	else
	{
		tempdate = maturitydate;
		while (tempdate > startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);
		schedule[0] = startdate;

		for (i = 1; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', (i-num_schedule) * frq, maturitydate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
};

void calculationstartschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	calculationstartschedule(startdate,	maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, schedule);

};

void calculationendschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<CDate>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, num_schedule, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		tempdate = startdate;
		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m', (i - 1) * frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);

		for (i = 0; i < num_schedule-1; i++)
		{
			schedule[i] = DateAdd('m', (i+1) * frq, startdate);
		}

		schedule[num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i = 0; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
	else
	{
		tempdate = maturitydate;
		while (tempdate > startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule = vector<CDate>(num_schedule);

		for (i = 0; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', (i - num_schedule+1) * frq, maturitydate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 0; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
};

void calculationendschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	calculationendschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, schedule);
};

void fixingschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	string fixlag,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(fixlag);
	int i, d = ymd[2];
	
	if (setin == "ARREAR")
	{
		calculationendschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			schedule
		);
	}
	else
	{
		calculationstartschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			schedule
		);

		int num_schedule = int(schedule.size());

		for (i=0;i<num_schedule;i++)
		{
			ShiftBusDate(schedule[i],fixingholidays,numfixingholiday,d,schedule[i]);
		};
	}
};

void fixingschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	string fixlag,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixingschedule
	(
		startdate,
		maturitydate,
		holidays,
		numholiday,
		fixingholidays,
		numfixingholiday,	
		stub,
		direction,
		conv,
		freq,
		adjflag,
		setin,
		fixlag,
		schedule
	);
};

void fixingschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	int fixlag,
	vector<CDate>& schedule
)
{
	int i;

	if (setin == "ARREAR")
	{
		calculationendschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			schedule
		);
	}
	else
	{
		calculationstartschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			schedule
		);

		int num_schedule = int(schedule.size());

		for (i = 0; i < num_schedule; i++)
		{
			ShiftBusDate(schedule[i], fixingholidays, numfixingholiday, fixlag, schedule[i]);
		};
	}
};

void fixingschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	int fixlag,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixingschedule
	(
		startdate,
		maturitydate,
		holidays,
		numholiday,
		fixingholidays,
		numfixingholiday,
		stub,
		direction,
		conv,
		freq,
		adjflag,
		setin,
		fixlag,
		schedule
	);
};

void depositeschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string adjflag,
	string payin,
	int spotlag,
	CDate& calstartdate,
	CDate& calenddate,
	CDate& paydate
)
{
	calstartdate = startdate;
	ConvReflDate(maturitydate, holidays, numholiday, conv, calenddate);
	ConvReflDate(calenddate, holidays, numholiday, conv, paydate);
};

void depositeschedule
(
	CDate today,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string adjflag,
	string payin,
	int spotlag,
	CDate& calstartdate,
	CDate& calenddate,
	CDate& paydate
)
{
	if (tenor == "o/n")
	{
		calstartdate = today;
		ShiftBusDate(calstartdate, holidays, numholiday, 1, calenddate);
		paydate = calenddate;

		if (payin != "ARREAR")
		{
			paydate = calstartdate;
		}
	}
	else if (tenor == "t/n")
	{
		ShiftBusDate(today, holidays, numholiday, 1, calstartdate);
		ShiftBusDate(calstartdate, holidays, numholiday, 1, calenddate);

		paydate = calenddate;

		if (payin != "ARREAR")
		{
			paydate = calstartdate;
		}
	}
	else
	{
		ShiftBusDate(today, holidays, numholiday, spotlag, calstartdate);

		if (tenor == "s/n")
		{
			ShiftBusDate(calstartdate, holidays, numholiday, 1, calenddate);

			paydate = calenddate;

			if (payin != "ARREAR")
			{
				paydate = calstartdate;
			}
		}
		else
		{
			vector<int> ymd = YMD2I(tenor);
			int d = ymd[2], m = ymd[2] * Nummonthayear + ymd[1];

			if (d>0)
			{
				calenddate = DateAdd('d', d, calstartdate);
				ConvReflDate(calenddate, holidays, numholiday, conv, paydate);

				if (adjflag == "ADJUSTED")
				{
					calenddate = paydate;
				}

				if (payin != "ARREAR")
				{
					paydate = calstartdate;
				}
			}
			else
			{
				calenddate = DateAdd('m', m, calstartdate);
				ConvReflDate(calenddate, holidays, numholiday, conv, paydate);

				if (adjflag == "ADJUSTED")
				{
					calenddate = paydate;
				}

				if (payin != "ARREAR")
				{
					paydate = calstartdate;
				}
			}

		}
	}
};

void fixedlegcashflowschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<vector<CDate>>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, num_schedule, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		tempdate = startdate;
		
		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m', (i - 1) * frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule[0] = vector<CDate>(num_schedule);
		schedule[1] = vector<CDate>(num_schedule);
		schedule[2] = vector<CDate>(num_schedule);

		schedule[1][0] = startdate;
		schedule[2][num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i=0;i<num_schedule-1;i++)
			{
				ConvReflDate(DateAdd('m', (i + 1) * frq, startdate), holidays, numholiday, conv, schedule[0][i]);

				schedule[1][i + 1] = schedule[0][i];
				schedule[2][i] = schedule[0][i];
			}

			ConvReflDate(schedule[2][num_schedule - 1], holidays, numholiday, conv, schedule[2][num_schedule - 1]);
		}
		else
		{
			for (i=0; i<num_schedule-1;i++)
			{
				schedule[1][i + 1] = DateAdd('m', (i + 1) * frq, startdate);
				
				ConvReflDate(schedule[1][i + 1], holidays, numholiday, conv, schedule[0][i]);
				schedule[2][i] = schedule[1][i + 1];
			}
		}
		ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[0][num_schedule - 1]);
	}
	else
	{
		tempdate = maturitydate;

		while (tempdate < startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}

		schedule[0] = vector<CDate>(num_schedule);
		schedule[1] = vector<CDate>(num_schedule);
		schedule[2] = vector<CDate>(num_schedule);

		schedule[1][0] = startdate;
		schedule[2][num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(DateAdd('m', (i-num_schedule) * frq, maturitydate), holidays, numholiday, conv, schedule[0][i-1]);

				schedule[1][i] = schedule[0][i-1];
				schedule[2][i-1] = schedule[0][i-1];
			}

			ConvReflDate(schedule[2][num_schedule - 1], holidays, numholiday, conv, schedule[2][num_schedule - 1]);
		}
		else
		{
			for (i = 1; i < num_schedule; i++)
			{
				schedule[1][i] = DateAdd('m', (i -num_schedule) * frq, maturitydate);

				ConvReflDate(schedule[1][i], holidays, numholiday, conv, schedule[0][i]);
				schedule[2][i] = schedule[1][i + 1];
			}
		}
		ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[0][num_schedule - 1]);
	}
};

void fixedlegcashflowschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<vector<CDate>>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixedlegcashflowschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, schedule);
};

void ShiftBusDate
(
	CDate date,
	CDate* holidays,
	int numholiday,
	int numshift,
	CDate& newdate
)
{
	newdate = date;

	int sgn = sign(numshift);

	if (sgn!=0)
	{
		int i;
		for (i=0;i<abs(numshift);i++)
		{
			newdate = DateAdd('d', sgn, newdate);
			while (is_holiday(newdate,holidays, numholiday))
			{
				newdate = DateAdd('d', sgn, newdate);
			}
		}
	}
};

void ConvReflDate
(
	CDate date,
	CDate* holidays,
	int numholiday,
	string conv,
	CDate& newdate
)
{
	newdate = date;
	if (conv == "MODFOL")
	{
		while (is_holiday(newdate, holidays, numholiday))
		{
			newdate = DateAdd('d',1,newdate);
		}

		if (newdate[1] != date[1])
		{
			newdate = DateAdd('d', -1, newdate);
			while (is_holiday(newdate, holidays, numholiday))
			{
				newdate = DateAdd('d', 1, newdate);
			}
		}
	}
	else if (conv == "FOL")
	{
		while (is_holiday(newdate, holidays, numholiday))
		{
			newdate = DateAdd('d', 1, newdate);
		}

	}
	else if (conv == "MODPREC")
	{
		while (is_holiday(newdate, holidays, numholiday))
		{
			newdate = DateAdd('d', -1, newdate);
		}

		if (newdate[1] != date[1])
		{
			newdate = DateAdd('d', +1, newdate);
			while (is_holiday(newdate, holidays, numholiday))
			{
				newdate = DateAdd('d', +1, newdate);
			}
		}
	}
	else if (conv == "PREC")
	{
		while (is_holiday(newdate, holidays, numholiday))
		{
			newdate = DateAdd('d', -1, newdate);
		}
	}
};

bool is_holiday
(
	CDate date,
	CDate* holidays,
	int numholiday
)
{
	int i, weekdayflag = get_weekday(date);
	bool flag = false;

	if (weekdayflag == 0 || weekdayflag == 6)
	{
		flag = true;
	}
	else
	{
		for (i=0; i<numholiday; i++)
		{
			if (date == holidays[i])
			{
				flag = true;
				break;
			}
		}
	}

	return flag;
};

int findnumschedule
(
	CDate startdate,
	CDate maturitydate,
	string stub,
	string direction,
	string freq
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, num_schedule, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		tempdate = startdate;

		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m', (i - 1) * frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule>1)
		{
			num_schedule = num_schedule - 1;
		}
	}
	else
	{
		tempdate = maturitydate;

		while (tempdate > startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}
	}

	return num_schedule;
};

int findnumschedule
(
	CDate startdate,
	string tenor,
	string stub,
	string direction,
	string freq
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);
	int num_schedule = findnumschedule(startdate, maturitydate, stub, direction, freq);

	return num_schedule;

};

int findnumschedule
(
	CDate startdate,
	CDate maturitydate,
	string stub,
	string direction,
	int frq
)
{
	CDate tempdate;
	int i = 1, num_schedule;

	if (direction == "FORWARD")
	{
		tempdate = startdate;

		while (tempdate < maturitydate)
		{
			i = i + 1;
			tempdate = DateAdd('m',(i-1)*frq, startdate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}
	}
	else
	{
		tempdate = maturitydate;

		while (tempdate < startdate)
		{
			i = i + 1;
			tempdate = DateAdd('m', -(i - 1) * frq, maturitydate);
		}

		num_schedule = i - 1;

		if (stub == "LONG" && num_schedule > 1)
		{
			num_schedule = num_schedule - 1;
		}
	}

	return num_schedule;
};

int findnumschedule
(
	CDate startdate,
	string tenor,
	string stub,
	string direction,
	int frq

)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);
	int num_schedule = findnumschedule(startdate, maturitydate, stub, direction, frq);

	return num_schedule;

};

void paymentschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(freq);
	int i=1, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		for (i=0; i<num_schedule-1;i++)
		{
			ConvReflDate(DateAdd('m', (i + 1) * frq, startdate), holidays, numholiday, conv, schedule[i]);
		}
		ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[num_schedule-1]);
	}
	else
	{
		for (i = 0; i < num_schedule; i++)
		{
			ConvReflDate(DateAdd('m', (i-num_schedule+1) * frq, maturitydate), holidays, numholiday, conv, schedule[i]);
		}
	}
};

void paymentschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	paymentschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, num_schedule, schedule);
};


void calculationstartschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		schedule[0] = startdate;

		for (i = 1; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', i * frq, startdate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
	else
	{
		schedule[0] = startdate;

		for (i = 1; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', (i - num_schedule) * frq, maturitydate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
};

void calculationstartschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	calculationstartschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, num_schedule, schedule);
};

void calculationendschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		for (i = 0; i < num_schedule - 1; i++)
		{
			schedule[i] = DateAdd('m', (i + 1) * frq, startdate);
		}

		schedule[num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i = 0; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
	else
	{
		for (i = 0; i < num_schedule; i++)
		{
			schedule[i] = DateAdd('m', (i - num_schedule + 1) * frq, maturitydate);
		}

		if (adjflag == "ADJUSTED")
		{
			for (i = 0; i < num_schedule; i++)
			{
				ConvReflDate(schedule[i], holidays, numholiday, conv, schedule[i]);
			}
		}
	}
};

void calculationendschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	calculationendschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, num_schedule, schedule);
};

void fixingschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	string fixlag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(fixlag);
	int i, d = ymd[2];

	if (setin == "ARREAR")
	{
		calculationendschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			num_schedule,
			schedule
		);
	}
	else
	{
		calculationstartschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			num_schedule,
			schedule
		);

		for (i = 0; i < num_schedule; i++)
		{
			ShiftBusDate(schedule[i], fixingholidays, numfixingholiday, d, schedule[i]);
		};
	}
};

void fixingschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	string fixlag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixingschedule
	(
		startdate,
		maturitydate,
		holidays,
		numholiday,
		fixingholidays,
		numfixingholiday,
		stub,
		direction,
		conv,
		freq,
		adjflag,
		setin,
		fixlag,
		num_schedule,
		schedule
	);
};


void fixingschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	int fixlag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	int i;

	if (setin == "ARREAR")
	{
		calculationendschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			num_schedule,
			schedule
		);
	}
	else
	{
		calculationstartschedule
		(
			startdate,
			maturitydate,
			holidays,
			numholiday,
			stub,
			direction,
			conv,
			freq,
			adjflag,
			num_schedule,
			schedule
		);

		for (i = 0; i < num_schedule; i++)
		{
			ShiftBusDate(schedule[i], fixingholidays, numfixingholiday, fixlag, schedule[i]);
		};
	}
};

void fixingschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	CDate* fixingholidays,
	int numfixingholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	string setin,
	int fixlag,
	int num_schedule,
	vector<CDate>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixingschedule
	(
		startdate,
		maturitydate,
		holidays,
		numholiday,
		fixingholidays,
		numfixingholiday,
		stub,
		direction,
		conv,
		freq,
		adjflag,
		setin,
		fixlag,
		num_schedule,
		schedule
	);
};


void fixedlegcashflowschedule
(
	CDate startdate,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<vector<CDate>>& schedule
)
{
	CDate tempdate;
	vector<int> ymd = YMD2I(freq);
	int i = 1, frq = ymd[0] * Nummonthayear + ymd[1];

	if (direction == "FORWARD")
	{
		schedule[1][0] = startdate;
		schedule[2][num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i = 0; i < num_schedule - 1; i++)
			{
				ConvReflDate(DateAdd('m', (i + 1) * frq, startdate), holidays, numholiday, conv, schedule[0][i]);

				schedule[1][i + 1] = schedule[0][i];
				schedule[2][i] = schedule[0][i];
			}

			ConvReflDate(schedule[2][num_schedule - 1], holidays, numholiday, conv, schedule[2][num_schedule - 1]);
		}
		else
		{
			for (i = 0; i < num_schedule - 1; i++)
			{
				schedule[1][i + 1] = DateAdd('m', (i + 1) * frq, startdate);

				ConvReflDate(schedule[1][i + 1], holidays, numholiday, conv, schedule[0][i]);
				schedule[2][i] = schedule[1][i + 1];
			}
		}
		ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[0][num_schedule - 1]);
	}
	else
	{
		schedule[1][0] = startdate;
		schedule[2][num_schedule - 1] = maturitydate;

		if (adjflag == "ADJUSTED")
		{
			for (i = 1; i < num_schedule; i++)
			{
				ConvReflDate(DateAdd('m', (i - num_schedule) * frq, maturitydate), holidays, numholiday, conv, schedule[0][i - 1]);

				schedule[1][i] = schedule[0][i - 1];
				schedule[2][i - 1] = schedule[0][i - 1];
			}

			ConvReflDate(schedule[2][num_schedule - 1], holidays, numholiday, conv, schedule[2][num_schedule - 1]);
		}
		else
		{
			for (i = 1; i < num_schedule; i++)
			{
				schedule[1][i] = DateAdd('m', (i - num_schedule) * frq, maturitydate);

				ConvReflDate(schedule[1][i], holidays, numholiday, conv, schedule[0][i-1]);
				schedule[2][i-1] = schedule[1][i];
			}
		}
		ConvReflDate(maturitydate, holidays, numholiday, conv, schedule[0][num_schedule - 1]);
	}
};

void fixedlegcashflowschedule
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	int num_schedule,
	vector<vector<CDate>>& schedule
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];
	CDate maturitydate = DateAdd('m', m, startdate);

	fixedlegcashflowschedule(startdate, maturitydate, holidays, numholiday, stub, direction, conv, freq, adjflag, num_schedule,schedule);
};

void findmaturity
(
	CDate startdate,
	int tenor,
	CDate* holidays,
	int num_holiday,
	string conv,
	CDate& maturity
)
{
	ConvReflDate(DateAdd('m', tenor, startdate), holidays, num_holiday, conv, maturity);
};

void findmaturity
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int num_holiday,
	string conv,
	CDate& maturity
)
{
	vector<int> ymd = YMD2I(tenor);
	int m = ymd[0] * Nummonthayear + ymd[1];

	findmaturity(startdate,m,holidays, num_holiday, conv, maturity);
};

CDate* Holidays(string holidaysname)
{
	if (holidaysname == "NYB")
	{
		return NYHolidays;
	}
	else if (holidaysname == "LNB")
	{
		return LNHolidays;
	}
	else if (holidaysname == "LNNYB")
	{
		return LNNYHolidays;
	}
	else if (holidaysname == "TKB")
	{
		return TKHolidays;
	}
	else if (holidaysname == "LNTKB")
	{
		return LNTKHolidays;
	}
	else if (holidaysname == "TGT")
	{
		return TGTHolidays;
	}
	else
	{
		return SEHolidays;
	}
};

int NumHolidays(string holidaysname)
{
	if (holidaysname == "NYB")
	{
		return NumNYHolidays;
	}
	else if (holidaysname == "LNB")
	{
		return NumLNHolidays;
	}
	else if (holidaysname == "LNNYB")
	{
		return NumLNNYHolidays;
	}
	else if (holidaysname == "TKB")
	{
		return NumTKHolidays;
	}
	else if (holidaysname == "LNTKB")
	{
		return NumLNTKHolidays;
	}
	else if (holidaysname == "TGT")
	{
		return NumTARGETHolidays;
	}
	else
	{
		return NumSEHolidays;
	}
};

void avgfixedschedule
(
	string avg_freq,
	double avg_interval,
	string avg_weekday,
	string avg_sched_from,
	string avg_fixing_setin,
	string avg_fixing_lag,
	string avg_fixing_conv,

	vector<CDate> _couponleg_calc_startdate,
	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,
	vector<CDate> _couponleg_paydate,

	vector<CDate*> index1_holidayss,
	vector<int> num_index1holidayss,

	vector<bool> _ra_flag,

	vector<vector<CDate>> couponindex_fixingdate,
	vector<CDate> date,

	vector<vector<CDate>>& couponindex_avg_fixeddate,
	vector<vector<bool>>& avg_fixedflag,

	vector<bool>& avg_flag,
	vector<int>& avg_fixeddenom
)
{
	int i, j, k, jj, jjj, num_remained_couponleg_cf = int(_couponleg_paydate.size()), numshift = YMD2I(avg_fixing_lag)[2], count = 0;

	CDate tempdate, fixingtempdate;

	if (avg_freq == MONTHLY)
	{
		if (avg_sched_from == "CALLSTART")
		{
		}
		else if (avg_sched_from == "CALLEND")
		{
		}
		else
		{
			for (i=0; i<num_remained_couponleg_cf; i++)
			{
				if (_ra_flag[i])
				{
					if (avg_fixing_setin == "ARREAR")
					{
					}
					else
					{
						count = 0;

						ConvReflDate(DateAdd('m',count, couponleg_calc_shiftedstartdate[i]),index1_holidayss[i], num_index1holidayss[i],"PREC",tempdate);

						couponindex_avg_fixeddate[i].push_back(tempdate);

						while (tempdate > couponleg_calc_shiftedstartdate[i])
						{
							count = count - int(avg_interval);

							ConvReflDate(DateAdd('m', count, couponleg_calc_shiftedstartdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

							couponindex_avg_fixeddate[i].insert(couponindex_avg_fixeddate[i].begin(),tempdate);
						}
						
						/*
						count = 0;

						ConvReflDate(DateAdd('m', count, couponleg_calc_shiftedstartdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

						ShiftBusDate(tempdate, index1_holidayss[i], num_index1holidayss[i],numshift,fixingtempdate);
						
						couponindex_avg_fixeddate[i].push_back(fixingtempdate);

						while (tempdate > _couponleg_calc_startdate[i])
						{
							count = count - int(avg_interval);

							ConvReflDate(DateAdd('m', count, _couponleg_calc_startdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

							ShiftBusDate(tempdate, index1_holidayss[i], num_index1holidayss[i], numshift, fixingtempdate);

							couponindex_avg_fixeddate[i].insert(couponindex_avg_fixeddate[i].begin(), fixingtempdate);
						}
						*/

						/*
						count = -int(avg_interval);

						ConvReflDate(DateAdd('m', count, _couponleg_paydate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

						couponindex_avg_fixeddate[i].push_back(tempdate);

						while (tempdate > _couponleg_calc_startdate[i])
						{
							count = count - int(avg_interval);

							ConvReflDate(DateAdd('m', count, _couponleg_paydate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

							couponindex_avg_fixeddate[i].insert(couponindex_avg_fixeddate[i].begin(), tempdate);
						}
						*/

						avg_fixeddenom[i] = int(couponindex_avg_fixeddate[i].size());
						jj = 0; jjj = 0;

						avg_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

						for (j=0;j<avg_fixeddenom[i];j++)
						{
							ShiftBusDate(couponindex_avg_fixeddate[i][j],index1_holidayss[i], num_index1holidayss[i],numshift,couponindex_avg_fixeddate[i][j]);
							
							for (k = jj; k < int(couponindex_fixingdate[i].size()); k++)
							{
								if (couponindex_avg_fixeddate[i][j] == couponindex_fixingdate[i][k])
								{
									avg_fixedflag[i][k] = true;
									jj = k;
									break;
								}
							}

							for (k = jjj; k < int(date.size()); k++)
							{
								if (couponindex_avg_fixeddate[i][j]==date[k])
								{
									avg_flag[k] = true;
									jjj = k;
									break;
								}
							}
						}

					}
				}
			}
		}
	}
	else if (avg_freq == WEEKLY)
	{
		CDate tmpstrtdate;
		int startday, targetday;

		if (avg_weekday == "MON")
		{
			targetday = 1;
		}
		else if (avg_weekday == "TUE")
		{
			targetday = 2;
		}
		else if (avg_weekday == "WED")
		{
			targetday = 3;
		}
		else if (avg_weekday == "THU")
		{
			targetday = 4;
		}
		else if (avg_weekday == "FRI")
		{
			targetday = 5;
		}
		else if (avg_weekday == "SAT")
		{
			targetday = 6;
		}
		else
		{
			targetday = 0;
		}

		for (i=0;i<num_remained_couponleg_cf;i++)
		{
			if (_ra_flag[i])
			{
				tmpstrtdate = couponleg_calc_shiftedstartdate[i];
				startday = get_weekday(tmpstrtdate);

				while (startday !=targetday)
				{
					tmpstrtdate = DateAdd('d', 1, tmpstrtdate);
					startday = get_weekday(tmpstrtdate);
				}

				count = 0;

			}
			
			ConvReflDate(tmpstrtdate, index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

			couponindex_avg_fixeddate[i].push_back(tempdate);

			while (tempdate < couponleg_calc_shiftedenddate[i])
			{
				count = count + 7 * int(avg_interval);

				ConvReflDate(DateAdd('d',count,tmpstrtdate), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

				couponindex_avg_fixeddate[i].push_back(tempdate);
			}

			avg_fixeddenom[i] = int(couponindex_avg_fixeddate[i].size());

			jj = 0; jjj = 0;

			avg_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

			for (j=0;j<avg_fixeddenom[i];j++)
			{
				ShiftBusDate(couponindex_avg_fixeddate[i][j],index1_holidayss[i],num_index1holidayss[i],numshift,couponindex_avg_fixeddate[i][j]);

				for (k=jj;k<int(couponindex_fixingdate[i].size());k++)
				{
					if (couponindex_avg_fixeddate[i][j]==couponindex_fixingdate[i][k])
					{
						avg_fixedflag[i][k] = true;
						jj = k;
						break;
					}
				}

				for (k=jjj;k<int(date.size());k++)
				{
					if (couponindex_avg_fixeddate[i][j]==date[k])
					{
						avg_flag[k] = true;
						jjj = k;
						break;
					}
				}
			}
		}
	}
	else
	{
	}
};

void StartEndSchedule
(
	string obsv_freq,
	int obsv_interval,
	string obsv_fixing_mode,
	int num_remained_couponleg_cf,
	
	vector<CDate*> index1_holidayss,
	
	vector<int> num_index1holidayss,
	vector<bool> _ra_flag,
	
	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,

	vector<vector<CDate>> couponindex_fixingdate,
	vector<CDate> date,

	vector<vector<CDate>>& adv_date,
	vector<vector<CDate>>& arr_date,

	vector<vector<bool>>& adv_fixedflag,
	vector<vector<bool>>& arr_fixedflag,

	vector<bool>& adv_flag,
	vector<bool>& arr_flag
)
{
	int i, j, jj, jjj, jjjj, jjjjj, k, count = 0;
	CDate tempdate;

	for (i=0;i<num_remained_couponleg_cf;i++)
	{
		if (obsv_freq=="MONTHLY")
		{
			if (_ra_flag[i] && obsv_interval==1)
			{
				count = 0;

				ConvReflDate(DateAdd('m', count, couponleg_calc_shiftedstartdate[i]),index1_holidayss[i],num_index1holidayss[i],"PREC",tempdate);

				adv_date[i].push_back(tempdate);

				ConvReflDate(DateAdd('m', count+obsv_interval, couponleg_calc_shiftedstartdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

				arr_date[i].push_back(tempdate);

				while (tempdate <= couponleg_calc_shiftedenddate[i])
				{
					count = count + obsv_interval;

					ConvReflDate(DateAdd('m',(-12+count),couponleg_calc_shiftedenddate[i]),index1_holidayss[i],num_index1holidayss[i],"PREC",tempdate);

					adv_date[i].push_back(tempdate);

					ConvReflDate(DateAdd('m', (-12 + count+obsv_interval), couponleg_calc_shiftedenddate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

					arr_date[i].push_back(tempdate);
				}

				adv_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				arr_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				jj = 0; jjj = 0; jjjj = 0; jjjjj = 0;

				for (j=0;j<int(adv_date[i].size());j++)
				{
					for (k=jj; k<int(couponindex_fixingdate[i].size()); k++)
					{
						if (adv_date[i][j]==couponindex_fixingdate[i][k])
						{
							adv_fixedflag[i][k] = true;
							jj = k;
							break;
						}
					}

					for (k=jjj; k<int(couponindex_fixingdate[i].size()); k++)
					{
						if (arr_date[i][j] == couponindex_fixingdate[i][k])
						{
							arr_fixedflag[i][k] = true;
							jjj = k;
							break;
						}

					}

					for (k = jjjj; k<int(date.size()); k++)
					{
						if (adv_date[i][j] == date[k])
						{
							adv_flag[k] = true;
							jjjj = k;
							break;
						}

					}

					for (k = jjjjj; k<int(date.size()); k++)
					{
						if (arr_date[i][j] == date[k])
						{
							arr_flag[k] = true;
							jjjjj = k;
							break;
						}
					}
				}
			}
			else if (_ra_flag[i] && obsv_interval > 1)
			{
				adv_date[i].push_back(couponleg_calc_shiftedstartdate[i]);

				arr_date[i].push_back(couponleg_calc_shiftedenddate[i]);

				/*
				ConvReflDate(DateAdd('m', count, couponleg_calc_shiftedstartdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

				adv_date[i].push_back(couponleg_calc_shiftedstartdate[i]);

				ConvReflDate(DateAdd('m', count+obsv_interval, couponleg_calc_shiftedstartdate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

				arr_date[i].push_back(couponleg_calc_shiftedstartdate[i]);
				*/

				//count = count + obsv_interval;

				/*
				while (tempdate <= couponleg_calc_shiftedenddate[i])
				{
					count = count + obsv_interval;

					ConvReflDate(DateAdd('m', (-12 + count), couponleg_calc_shiftedenddate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

					ConvReflDate(DateAdd('m', (count), couponleg_calc_shiftedenddate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

					adv_date[i].push_back(tempdate);

					ConvReflDate(DateAdd('m', (-12 + count+obsv_interval), couponleg_calc_shiftedenddate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

					ConvReflDate(DateAdd('m', (count+obsv_interval), couponleg_calc_shiftedenddate[i]), index1_holidayss[i], num_index1holidayss[i], "PREC", tempdate);

					arr_date[i].push_back(tempdate);
				}
				*/

				adv_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				arr_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				jj = 0; jjj = 0; jjjj = 0; jjjjj = 0;

				for (j=0;j<int(adv_date[i].size());j++)
				{
					for (k=jj;k<int(couponindex_fixingdate[i].size());k++)
					{
						if (adv_date[i][j]==couponindex_fixingdate[i][k])
						{
							adv_fixedflag[i][k] = true;
							jj = k;
							break;
						}
					}

					for (k = jjj; k<int(couponindex_fixingdate[i].size()); k++)
					{
						if (arr_date[i][j] == couponindex_fixingdate[i][k])
						{
							arr_fixedflag[i][k] = true;
							jjj = k;
							break;
						}
					}

					for (k = jjjj; k < int(date.size()); k++)
					{
						if (adv_date[i][j] == date[k])
						{
							adv_flag[k] = true;
							jjjj = k;
							break;
						}
					}

					for (k = jjjjj; k < int(date.size()); k++)
					{
						if (arr_date[i][j] == date[k])
						{
							arr_flag[k] = true;
							jjjjj = k;
							break;
						}
					}

				}

			}
		}
	}
};

void MinMaxSchedule
(
	string obsv_freq,
	int obsv_interval,
	string obsv_fixing_mode,
	int num_remained_couponleg_cf,

	vector<CDate*> index1_holidayss,

	vector<int> num_index1holidayss,
	vector<bool> _ra_flag,

	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,

	vector<vector<CDate>> couponindex_fixingdate,
	vector<CDate> date,

	vector<vector<CDate>>& adv_date,
	vector<vector<CDate>>& arr_date,

	vector<vector<bool>>& adv_fixedflag,
	vector<vector<bool>>& arr_fixedflag,

	vector<bool>& adv_flag,
	vector<bool>& arr_flag
)
{
	int i, j, jj, jjj, jjjj, jjjjj, k, count = 0;
	CDate tempdate;

	for (i = 0; i < num_remained_couponleg_cf; i++)
	{
		if (obsv_freq == "WEEKLY")
		{
			if (_ra_flag[i])
			{
				count = 0;

				ConvReflDate
				(
					DateAdd('d', count, couponleg_calc_shiftedstartdate[i]),
					index1_holidayss[i],
					num_index1holidayss[i],
					"PREC",
					tempdate
				);

				adv_date[i].push_back(tempdate);

				ConvReflDate
				(
					DateAdd('d', count + obsv_interval, couponleg_calc_shiftedstartdate[i]),
					index1_holidayss[i],
					num_index1holidayss[i],
					"PREC",
					tempdate
				);

				arr_date[i].push_back(tempdate);

				while (tempdate <= couponleg_calc_shiftedenddate[i])
				{
					count = count + obsv_interval;

					ConvReflDate
					(
						DateAdd('d', (-365 + count),couponleg_calc_shiftedenddate[i]),
						index1_holidayss[i],
						num_index1holidayss[i], 
						"PREC", 
						tempdate
					);

					adv_date[i].push_back(tempdate);

					ConvReflDate
					(
						DateAdd('d', (-365 + count + obsv_interval), couponleg_calc_shiftedenddate[i]),
						index1_holidayss[i], 
						num_index1holidayss[i], 
						"PREC", 
						tempdate
					);

					arr_date[i].push_back(tempdate);
				}

				adv_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				arr_fixedflag[i] = vector<bool>(int(couponindex_fixingdate[i].size()), false);

				jj = 0; jjj = 0; jjjj = 0; jjjjj = 0;

				for (j = 0; j < int(adv_date[i].size()); j++)
				{
					for (k = jj; k < int(couponindex_fixingdate[i].size()); k++)
					{
						if (adv_date[i][j] == couponindex_fixingdate[i][k])
						{
							adv_fixedflag[i][k] = true;
							jj = k;
							break;
						}
					}

					for (k = jjj; k < int(couponindex_fixingdate[i].size()); k++)
					{
						if (arr_date[i][j] == couponindex_fixingdate[i][k])
						{
							arr_fixedflag[i][k] = true;
							jjj = k;
							break;
						}

					}

					for (k = jjjj; k < int(date.size()); k++)
					{
						if (adv_date[i][j] == date[k])
						{
							adv_flag[k] = true;
							jjjj = k;
							break;
						}

					}

					for (k = jjjjj; k < int(date.size()); k++)
					{
						if (arr_date[i][j] == date[k])
						{
							arr_flag[k] = true;
							jjjjj = k;
							break;
						}
					}
				}
			}
		}
	}
};

/*
void avgfixingschedule
(
	string avg_freq,
	double avg_interval,
	string avg_weekday,
	string avg_sched_from,
	string avg_fixing_setin,
	string avg_fixing_lag,
	string avg_fixing_conv,

	vector<CDate> _couponleg_calc_startdate,
	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,
	vector<CDate> _couponleg_paydate,

	vector<CDate*> index1_holidayss,
	vector<int> num_index1holidayss,

	vector<bool> _ra_flag,

	vector<vector<CDate>> couponindex_fixingdate,
	vector<vector<CDate>> couponindex_avg_fixingdate,
	
	vector<vector<bool>>& avg_flag,
	vector<int>& avg_denom
)
{
	int i, j, num_remained_couponleg_cf = int(_couponleg_paydate.size()), numshift = YMD2I(avg_fixing_lag)[2], count = 0;

	CDate tempdate;

	if (avg_freq == MONTHLY)
	{
		if (avg_sched_from == "CALLSTART")
		{
		}
		else if (avg_sched_from == "CALLEND")
		{
		}
		else
		{
			for (i = 0; i < num_remained_couponleg_cf; i++)
			{
				if (avg_fixing_setin == "ARREAR")
				{
				}
				else
				{
					count = -int(avg_interval);

					ConvReflDate
					(
						DateAdd('m', count, _couponleg_paydate[i]),
						index1_holidayss[i], 
						num_index1holidayss[i], 
						"PREC", 
						tempdate
					);

					while (tempdate > couponleg_calc_shiftedstartdate[i])
					{
						count = count - int(avg_interval);

						ConvReflDate
						(
							DateAdd('m', count, _couponleg_paydate[i]),
							index1_holidayss[i],
							num_index1holidayss[i],
							"PREC", 
							tempdate
						);

						couponindex_avg_fixingdate[i].insert(couponindex_avg_fixingdate[i].begin(), tempdate);
					}

					avg_denom[i] = int(couponindex_fixingdate[i].size());

					for (j=0; j < avg_denom[i];j++)
					{
						ShiftBusDate(couponindex_fixingdate[i][j],index1_holidayss[i],num_index1holidayss[i],numshift,couponindex_avg_fixingdate[i][j]);
					}
				}		
			}
		}
	}
	else
	{
	}
};

*/

CDate find_nth_day_of_yyyymm
(
	int yyyy,
	int mm,
	int nth,
	int day,
	CDate* holidays,
	int numholiday
)
{
	CDate firstdate = CDate(yyyy, mm, 1);
	int tmpday = get_weekday(firstdate);
	CDate targetfirstdate, targetnthdate, targetnextdate;

	if (tmpday <=day)
	{
		targetfirstdate = DateAdd('d', day - tmpday, firstdate);
	}
	else
	{
		targetfirstdate = DateAdd('d', day + 7 - tmpday, firstdate);
	}

	targetnthdate = targetfirstdate;
	targetnextdate = targetnthdate;
	int count = 0;

	while (count < nth)
	{
		targetnthdate = targetnextdate;
		
		if (is_holiday(targetnthdate, holidays, numholiday))
		{
			targetnthdate = DateAdd('d', 7, targetnthdate);
		}
		else
		{
			count = count + 1;
			targetnextdate = DateAdd('d', 7, targetnthdate);
		}
	}

	return targetnthdate;
};

CDate find_nth_day_of_MMMYY
(
	string MMMYY,
	int nth,
	int day,
	CDate* holidays,
	int numholiday
)
{
	return find_nth_day_of_yyyymm(get_year_from_MMMYY(MMMYY),get_month_from_MMMYY(MMMYY),nth,day,holidays,numholiday);
};

void shortfuturesschedule
(
	string mmmyy,
	CDate* holiday,
	int numholiday,
	int mm,
	int nth,
	int day,
	string conv,
	int fixinglag,
	CDate& calstartdate,
	CDate& calenddate,
	CDate& lasttradedate
)
{
	calstartdate = find_nth_day_of_MMMYY(mmmyy, nth, day, holiday, numholiday);
	calenddate = DateAdd('m', mm, calstartdate);

};

void fxswapschedule
(
	CDate today,
	string swappointtenor,
	CDate* fxswapholiday,
	int numfxswapholiday,
	CDate* fxswapstartholiday,
	int numfxswapstartholiday,
	string fxswapconv,
	int fxswaplag,
	CDate& calstartdate,
	CDate& calenddate
)
{
	if (swappointtenor == "o/n")
	{
		calstartdate = today;
		ShiftBusDate(calstartdate, fxswapholiday, numfxswapholiday, 1, calenddate);
	}
	else if(swappointtenor == "t/n")
	{
		ShiftBusDate(today, fxswapstartholiday, numfxswapstartholiday, 1, calstartdate);
		ShiftBusDate(calstartdate,fxswapholiday,numfxswapholiday,1,calenddate);
	}
	else
	{
		ShiftBusDate(today,fxswapstartholiday,numfxswapstartholiday,fxswaplag,calstartdate);
		vector<int> ymd = YMD2I(swappointtenor);
		int d = ymd[2], m = ymd[0] * Nummonthayear + ymd[1];

		if (d > 0)
		{
			calenddate = DateAdd('d', d, calstartdate);
			ConvReflDate(calenddate, fxswapholiday, numfxswapholiday, fxswapconv, calenddate);
		}
		else
		{
			calenddate = DateAdd('m', m, calstartdate);
			ConvReflDate(calenddate, fxswapholiday, numfxswapholiday, fxswapconv, calenddate);
		}
	}
}
void MinMaxSchedule
(
	string obsv_freq,
	int obsv_interval,
	string obsv_fixing_mode,
	int num_remained_couponleg_cf,
	vector<CDate*> index1_holidays,
	vector<int> num_index1holidays,
	vector<bool> _ra_flag,
	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,
	vector<vector<CDate>>& adv_date,
	vector<vector<CDate>>& arr_date,
	vector<vector<bool>>& adv_fixedflag,
	vector<vector<bool>>& arr_fixedflag,
	vector<bool>& adv_flag,
	vector<bool>& arr_flag)
{
}
;





