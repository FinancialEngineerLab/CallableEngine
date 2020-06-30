#include "Instrument.h"

CInstrument::CInstrument()
{

}

CInstrument::CInstrument
(
	CCurrency CRCY,
	string TYPE,
	CDate* HOLIDAY,
	int NUMHOLIDAY,
	string STUB,
	string DIRECTION,
	string CONV,
	string ADJFLAG,
	string PAYIN,
	string SETIN,
	string BASIS,
	int SPOTLAG,
	CDate* FIXINGHOLIDAY,
	int NUMFIXINGHOLIDAY
)
{
	crcy = CRCY,
	type = TYPE;
	holiday = HOLIDAY;
	numholiday = NUMHOLIDAY;
	stub = STUB;
	direction = DIRECTION;
	conv = CONV;
	adjflag = ADJFLAG;
	payin = PAYIN;
	setin = SETIN;
	basis = BASIS;
	spotlag = SPOTLAG;
	fixingholiday = FIXINGHOLIDAY;
	numfixingholiday = NUMFIXINGHOLIDAY;
}

CInstrument::CInstrument
(
	CCurrency CRCY,
	string TYPE,
	CDate* HOLIDAY,

	int NUMHOLIDAY,
	CDate* STARTHOLIDAY,
	int NUMSTARTHOLIDAY,

	string STUB,
	string DIRECTION,
	string CONV,
	string ADJFLAG,
	string PAYIN,
	string SETIN,
	string BASIS,

	int SPOTLAG,
	CDate* FIXINGHOLIDAY,
	int NUMFIXINGHOLIDAY
)
{
	crcy = CRCY;
	type = TYPE;
	holiday = HOLIDAY;
	numholiday = NUMHOLIDAY;
	startholiday = STARTHOLIDAY;
	numstartholiday = NUMSTARTHOLIDAY;
	stub = STUB;
	direction = DIRECTION;
	conv = CONV;
	adjflag = ADJFLAG;
	payin = PAYIN;
	setin = SETIN;
	basis = BASIS;
	spotlag = SPOTLAG;
	fixingholiday = FIXINGHOLIDAY;
	numfixingholiday = NUMFIXINGHOLIDAY;

}

CInstrument::CInstrument(CCurrency CRCY, string TYPE)
{
	crcy = CRCY;
	type = TYPE;
	string symbol = crcy.get_symbol();

	if (symbol == "USD")
	{
		holiday = NYHolidays;
		numholiday = NumNYHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/360";
		spotlag = 2;
		fixingholiday = LNHolidays;
		numfixingholiday = NumLNHolidays;

		if (type == "SWAP")
		{
			holiday = LNNYHolidays;
			numholiday = NumLNNYHolidays;
		}
		else if (type == "SFUT")
		{
			holiday = CMELNHolidays;
			numholiday = NumCMELNHolidays;
			futmm = 3;
			futnth = 3;
			futday = 3;
		}
		else if (type == "BASIS")
		{
			holiday = LNNYHolidays;
			numholiday = NumLNNYHolidays;
		}
	}
	else if (symbol=="JPY")
	{
		holiday = TKHolidays;
		numholiday = NumTKHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/360";
		spotlag = 2;
		fixingholiday = TKHolidays;
		numfixingholiday = NumTKHolidays;

		if (type == "SWAP")
		{
			holiday = LNTKHolidays;
			numholiday = NumTKHolidays;
			basis = "ACT/365NL";
		}
	}
	else if (symbol=="EUR")
	{
		holiday = TGTHolidays;
		numholiday = NumTARGETHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/360";
		spotlag = 2;
		fixingholiday = TGTHolidays;
		numfixingholiday = NumTARGETHolidays;

		if (type == "SWAP")
		{
			basis = "30/360";
		}
	}
	else if (symbol=="USDKRW")
	{
		holiday = NYSEHolidays;
		numholiday = NumNYSEHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/365";
		spotlag = 2;
		fixingholiday = LNHolidays;
		numfixingholiday = NumLNHolidays;

		if (type == "FX")
		{
			startholiday = SEHolidays;
			numstartholiday = NumSEHolidays;
			fixingholiday = SEHolidays;
			numfixingholiday = NumSEHolidays;
		}

	}
	else if (symbol == "USDJPY")
	{
		holiday = NYTKHolidays;
		numholiday = NumNYTKHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/365";
		spotlag = 2;
		fixingholiday = LNHolidays;
		numfixingholiday = NumLNHolidays;

		if (type == "FX")
		{
			startholiday = TKHolidays;
			numstartholiday = NumTKHolidays;
		}
		else if (type == "BASIS")
		{
			holiday = LNTKHolidays;
			numholiday = NumLNNYTKHolidays;

		}
		else if (type == "SWAP")
		{
			holiday = LNTKHolidays;
			numholiday = NumLNTKHolidays;
			basis = "ACT/365NL";
		}
	}
	else if (symbol == "EURUSD")
	{
		holiday = NYTGTHolidays;
		numholiday = NumNYTGTHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/365";
		spotlag = 2;
		fixingholiday = TGTHolidays;
		numfixingholiday = NumTARGETHolidays;

		if (type == "FX")
		{
			startholiday = TGTHolidays;
			numstartholiday = NumTARGETHolidays;
		}
		else if (type == "BASIS")
		{
			holiday = LNNYTGTHolidays;
			numholiday = NumLNNYTGTHolidays;

		}
		else if (type == "SWAP")
		{
			holiday = TGTHolidays;
			numholiday = NumTARGETHolidays;
			basis = "30/360";
		}

	}
	else
	{
		holiday = SEHolidays;
		numholiday = NumSEHolidays;
		stub = "SHORT";
		direction = "BACKWARD";
		conv = "MODFOL";
		adjflag = "ADJUSTED";
		payin = "ARREAR";
		setin = "ADVANCE";
		basis = "ACT/365";
		spotlag = 1;
		fixingholiday = SEHolidays;
		numfixingholiday = NumSEHolidays;

		if (type == "BOND")
		{
			holiday = KRWHolidays;
			numholiday = NumKRWHolidays;
			basis = "ACT/ACTKRW";
			bondquotedigit = 2;
		}

		if (type == "KTBSWAP")
		{
			adjflag = "UNADJUSTED";
		}

	}

}

CInstrument::~CInstrument()
{

}

CCurrency CInstrument::get_currency()
{
	return crcy;
}

string CInstrument::get_type()
{
	return type;
}

CDate* CInstrument::get_holiday()
{
	return holiday;
}

int CInstrument::get_numholiday()
{
	return numholiday;
}

string CInstrument::get_stub()
{
	return stub;
}

string CInstrument::get_direction()
{
	return direction;
}

string CInstrument::get_conv()
{
	return conv;
}

string CInstrument::get_adjflag()
{
	return adjflag;
}

string CInstrument::get_payin()
{
	return payin;
}

string CInstrument::get_setin()
{
	return setin;
}

string CInstrument::get_basis()
{
	return basis;
}

int CInstrument::get_spotlag()
{
	return spotlag;
}

int CInstrument::get_futmm()
{
	return futmm;
}

int CInstrument::get_futnth()
{
	return futnth;
}

int CInstrument::get_futday()
{
	return futday;
}

int CInstrument::get_bondquotedigit()
{
	return bondquotedigit;
}

CDate* CInstrument::get_fixingholiday()
{
	return fixingholiday;
}

int CInstrument::get_numfixingholiday()
{
	return numfixingholiday;
}

CDate* CInstrument::get_startholiday()
{
	return fixingholiday;
}

int CInstrument::get_numstartholiday()
{
	return numfixingholiday;
}

void depositeamount
(
	CCurrency crcy,
	double notional,
	CDate today,
	string tenor,
	double rate,
	double& amt
)
{
	CDate calstartdate, calenddate, paydate;
	CInstrument depo(crcy, "DEPO");
	
	CDate* holiday = depo.get_holiday();
	int numholiday = depo.get_numholiday();

	string basis = depo.get_basis(), stub = depo.get_stub(), direction = depo.get_direction(), conv = depo.get_conv(), adjflag = depo.get_adjflag(), payin = depo.get_payin(), setin = depo.get_setin();

	int spotflag = depo.get_spotlag();

	depositeschedule(today, tenor, holiday, numholiday, stub, direction, conv, adjflag, payin, spotflag, calstartdate, calenddate, paydate);

	amt = notional * (1.0 + rate * cvg(calstartdate, calenddate, basis));
}

void depositeamount(double notional, CDate calstartdate, CDate calenddate, double rate, string basis, double& amt)
{
	amt = notional * (1.0 + rate * cvg(calstartdate, calenddate, basis));
}

void shortfuturesamount(double notional, CDate calstartdate, CDate calenddate, double price, string basis, double& amt)
{
	amt = notional * (1.0 + (1.0 - price) * cvg(calstartdate, calenddate, basis));
}


