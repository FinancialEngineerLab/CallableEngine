#pragma once

#ifndef _INSTRUMENT_H_
#define _INSTRUMENT_H_

#include <string>
#include <stdlib.h>
#include "Currency.h"
#include "Schedule.h"

using namespace std;

class CInstrument
{
private:
	CCurrency crcy;
	string type;
	CDate* holiday;
	
	int numholiday;
	
	string stub;
	string direction;
	string conv;
	string adjflag;
	string payin;
	string basis;
	string setin;

	int spotlag;
	int futmm;
	int futnth;
	int futday;
	int bondquotedigit;

	CDate* fixingholiday;
	int numfixingholiday;
	CDate* startholiday;
	int numstartholiday;

public:
	CInstrument();
	~CInstrument();

	CInstrument
	(
		CCurrency crcy,
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
	);

	CInstrument
	(
		CCurrency crcy,
		string TYPE,
		CDate* HOLIDAY,

		int NUMHOLIDAY,

		CDate* startholiday,
		int numstartholiday,

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
	);

	CInstrument(CCurrency CRCY, string TYPE);

	CCurrency get_currency();
	string get_type();

	CDate* get_holiday();
	int get_numholiday();
	CDate* get_startholiday();
	int get_numstartholiday();

	string get_stub();
	string get_direction();
	string get_conv();
	string get_adjflag();
	string get_payin();
	string get_basis();
	string get_setin();
	
	int get_spotlag();
	int get_futmm();
	int get_futnth();
	int get_futday();
	int get_bondquotedigit();

	CDate* get_fixingholiday();
	int get_numfixingholiday();

};

void depositeamount(CCurrency crcy, double notional, CDate today, string tenor, double rate, double& amt);

void depositeamount(double notional, CDate calstartdate, CDate calenddate, double rate, string basis, double& amt);

void shortfuturesamount(double notional, CDate calstartdate, CDate calenddate, double price, string basis, double& amt);

#endif 