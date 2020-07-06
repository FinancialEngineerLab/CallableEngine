#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MyUtility.h"
#include "Interpolation.h"
#include "InterestRateCurve.h"
#include "ConstantNumbers.h"
#include "Instrument.h"
#include "PlainPricer.h"

// void main
void main
(
	int argc,
	char** argv
)
{
	cout << "My Simple Callable Swap Pricer" << endl;

	int
		i,
		j,
		k,
		l;

	char *label = new char[100];

	ifstream fin;
	fin.open(argv[1]);

	string
		str_temp
		, couponleg_conv
		, couponleg_freq
		, couponleg_adjflag
		, couponleg_dcb
		, couponleg_stub
		, couponleg_direction
		, couponleg_holiday
		, fundingleg_conv
		, fundingleg_freq
		, fundingleg_adjflag
		, fundingleg_dcb
		, fundingleg_stub
		, fundingleg_direction
		, fundingleg_holiday

		, fundinglegindex_type
		, fundinglegindex_tenor
		, fundinglegindexfixed_freq
		, fundinglegindexfloating_freq

		, fixing_setin
		, fixing_holiday
		, fixing_lag

		, call_freq
		, call_conv
		, call_stub
		, call_direction
		, call_holiday

		, parameterschedfrom
		, calibratoin_type;

	double
		notional
		, couponleg_fixedrate
		, couponleg_margin
		, fundingleg_margin
		, fundingleg_indexmulti
		, callfee
		, input_a
		, input_sig
		, deltabumpingsize
		, gammabumpingsize
		, swaptionvegabumpingsize
		, dT
		, dx
		, quantile;

	int
		int_tmp
		, couponleg_frq
		, fundingleg_frq
		, fundinglegindx_tenor
		, fix_lag
		, call_frq
		, callperiod
		, num_mktrate
		, num_swaptionoption
		, num_swaptionswap
		, num_couponleg_cf
		, num_fundingleg_cf
		, num_calldate
		, num_parameter
		, num_fixing_history;

	bool
		coupleg_cfgen_flag
		, fundingleg_cfgen_flag
		, callable_flag
		, call_schedgen_flag
		, parameter_schedgen_flag
		, calibration_flag
		, deltaflag
		, gammaflag
		, swaptionvegaflag;

	vector<int> ymd;

	fin >> label; fin >> str_temp; CDate today(str_temp);
	fin >> label; fin >> str_temp; CDate start_date(str_temp);
	fin >> label; fin >> str_temp; CDate maturity(str_temp);
	fin >> label; fin >> str_temp; CDate settlement_date(str_temp);

	fin >> label; fin >> notional;
	fin >> label; fin >> str_temp; CCurrency crcy(str_temp);

	fin >> label; fin >> couponleg_fixedrate;
	fin >> label; fin >> couponleg_margin;
	fin >> label; fin >> couponleg_freq; ymd = YMD2I(couponleg_freq); couponleg_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponleg_conv;
	fin >> label; fin >> couponleg_adjflag;
	fin >> label; fin >> couponleg_dcb;
	fin >> label; fin >> couponleg_stub;
	fin >> label; fin >> couponleg_direction;
	fin >> label; fin >> couponleg_holiday;
	fin >> label; fin >> coupleg_cfgen_flag;

	fin >> label; fin >> fundingleg_margin;
	fin >> label; fin >> fundingleg_freq; ymd = YMD2I(fundingleg_freq); fundingleg_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> fundingleg_indexmulti;
	fin >> label; fin >> fundingleg_conv;
	fin >> label; fin >> fundingleg_adjflag;
	fin >> label; fin >> fundingleg_dcb;
	fin >> label; fin >> fundingleg_stub;
	fin >> label; fin >> fundingleg_direction;
	fin >> label; fin >> fundingleg_holiday;
	fin >> label; fin >> fundingleg_cfgen_flag;

	fin >> label; fin >> str_temp; CCurrency fundinglegindex_crcy(str_temp);
	fin >> label; fin >> 
		fundinglegindex_tenor; 
		ymd = YMD2I(fundinglegindex_tenor);
		fundinglegindx_tenor = ymd[0] * Nummonthayear + ymd[1];

	fin >> label; fin >> fundinglegindex_type;
	fin >> label; fin >> fundinglegindexfixed_freq;
	fin >> label; fin >> fundinglegindexfloating_freq;
	
	fin >> label; fin >> fixing_setin;
	fin >> label; fin >> fixing_lag;
	ymd = YMD2I(fixing_lag);
	fix_lag = ymd[2];
	
	fin >> label; fin >> fixing_holiday;

	fin >> label; fin >> callable_flag;

	fin >> label; fin >> str_temp;
	CDate callstart_date(str_temp);
	fin >> label; fin >> str_temp;
	CDate callend_date(str_temp);

	fin >> label; fin >> call_freq;
	ymd = YMD2I(call_freq); 
	call_frq = ymd[0] * Nummonthayear + ymd[1];

	fin >> label; fin >> call_conv;
	fin >> label; fin >> call_stub;
	fin >> label; fin >> call_direction;
	fin >> label; fin >> call_holiday;

	fin >> label; fin >> str_temp;
	ymd = YMD2I(str_temp);
	callperiod = ymd[2];

	fin >> label; fin >> str_temp;
	CDate call_date(str_temp);

	fin >> label; fin >> callfee;
	fin >> label; fin >> call_schedgen_flag;

	fin >> label; fin >> input_a;
	fin >> label; fin >> input_sig;
	fin >> label; fin >> parameterschedfrom;
	fin >> label; fin >> parameter_schedgen_flag; 

	fin >> label; fin >> calibration_flag;
	fin >> label; fin >> calibratoin_type;

	fin >> label; fin >> deltaflag;
	fin >> label; fin >> deltabumpingsize;

	fin >> label; fin >> gammaflag;
	fin >> label; fin >> gammabumpingsize;

	fin >> label; fin >> swaptionvegaflag;
	fin >> label; fin >> swaptionvegabumpingsize;

	fin >> label; fin >> num_mktrate;

	vector<string>
		tenor(num_mktrate)
		, type(num_mktrate)
		, fixedrateleg_freq(num_mktrate)
		, floatingrateleg_freq(num_mktrate);

	vector<double> mkt_rate(num_mktrate);

	for (i=0; i<num_mktrate; i++)
	{
		fin >> tenor[i];
		fin >> type[i];
		fin >> fixedrateleg_freq[i];
		fin >> floatingrateleg_freq[i];
		fin >> mkt_rate[i];
	}

	fin >> label; fin >> num_swaptionoption;

	vector<string> swaptionoptiontenor(num_swaptionoption);

	for (i=0; i<num_swaptionoption; i++)
	{
		fin >> swaptionoptiontenor[i];
	}

	fin >> label; fin >> num_swaptionswap;

	vector<string> swaptionswaptenor(num_swaptionswap);

	for (i = 0; i < num_swaptionswap; i++)
	{
		fin >> swaptionswaptenor[i];
	}

	vector<vector<double>> mktswaptionvol(num_swaptionswap);

	fin >> label;

	for (i = 0; i<num_swaptionswap; i++)
	{
		mktswaptionvol[i] = vector<double>(num_swaptionoption);

		for (j = 0; j<num_swaptionoption; j++)
		{
			fin >> mktswaptionvol[i][j];
		}
	}

	fin >> label; fin >> dT;
	fin >> label; fin >> dx;
	fin >> label; fin >> quantile;

	CDate * couponleg_holidays = Holidays(couponleg_holiday)
		, * fundingleg_holidays = Holidays(fundingleg_holiday)
		, * fixing_holidays = Holidays(fixing_holiday)
		, * call_holidays = Holidays(call_holiday);
	
	int num_couponlegholidays = NumHolidays(couponleg_holiday)
		, num_fundinglegholidays = NumHolidays(fundingleg_holiday)
		, num_callholidays = NumHolidays(call_holiday)
		, num_fixingholidays = NumHolidays(fixing_holiday);

	vector<CDate>
		couponleg_calc_startdate
		, couponleg_calc_enddate
		, couponleg_paydate
		, fundingleg_calc_startdate
		, fundingleg_calc_enddate
		, fundingleg_paydate
		, fixing_date
		, fixingindex_maturity
		, calleffective_date
		, callnotice_date
		, parameter_date;

	vector<double>
		couponleg_notional
		, couponleg_couponrate
		, couponleg_spread
		, couponlegindex1_mult
		, rahigh_bdry
		, ralow_bdry
		, rahigh_bdryoh
		, ralow_bdryoh
		, fundingleg_notional
		, fundingleg_mult
		, fundingleg_spread
		, callexe_fee
		, a
		, sig
		, imsi_sig;

	vector<bool>
		ra_flag
		, rahigh_bdryin_flag
		, ralow_bdryin_flag
		, callability_flag;

	vector<int>
		couponlegindx1_tenors
		, couponlegindex1fixed_frqs
		, rahigh_bdryoh_flag
		, ralow_bdryoh_flag
		, num_couponlegindex1holidays;

	vector<CCurrency> couponlegindex1_crcys;

	vector<string>
		couponlegindex1_tenors
		, couponlegindex1_types
		, couponlegindex1fixed_freqs
		, couponlegindex1floating_freqs;

	vector<CDate*> couponlegindex1_holidays;

	if (coupleg_cfgen_flag)
	{
		num_couponleg_cf = findnumschedule(start_date, maturity, couponleg_stub, couponleg_direction, couponleg_frq);

		couponleg_calc_startdate = vector<CDate>(num_couponleg_cf);

		calculationstartschedule
		(
			start_date
			, maturity
			, couponleg_holidays
			, num_couponlegholidays
			, couponleg_stub
			, fundingleg_direction
			, couponleg_conv
			, couponleg_freq
			, couponleg_adjflag
			, num_couponleg_cf
			, couponleg_calc_startdate
		);

		couponleg_calc_enddate = vector<CDate>(num_couponleg_cf);

		calculationendschedule
		(
			start_date
			, maturity
			, couponleg_holidays
			, num_couponlegholidays
			, couponleg_stub
			, fundingleg_direction
			, couponleg_conv
			, couponleg_freq
			, couponleg_adjflag
			, num_couponleg_cf
			, couponleg_calc_enddate
		);

		couponleg_paydate = vector<CDate>(num_couponleg_cf);

		paymentschedule
		(
			start_date
			, maturity
			, couponleg_holidays
			, num_couponlegholidays
			, couponleg_stub
			, fundingleg_direction
			, couponleg_conv
			, couponleg_freq
			, num_couponleg_cf
			, couponleg_paydate
		);

		couponleg_notional = vector<double>(num_couponleg_cf);
		couponleg_couponrate = vector<double>(num_couponleg_cf);
		couponleg_spread = vector<double>(num_couponleg_cf);
		
	}
	else
	{
		fin >> label; fin >> num_couponleg_cf;

		couponleg_calc_startdate = vector<CDate>(num_couponleg_cf);
		couponleg_calc_enddate = vector<CDate>(num_couponleg_cf);
		couponleg_paydate = vector<CDate>(num_couponleg_cf);

		couponleg_notional = vector<double>(num_couponleg_cf);
		couponleg_couponrate = vector<double>(num_couponleg_cf);
		couponleg_spread = vector<double>(num_couponleg_cf);

		for (i=0; i<num_couponleg_cf; i++)
		{
			fin >> str_temp; couponleg_calc_startdate[i] = CDate(str_temp);
			fin >> str_temp; couponleg_calc_enddate[i] = CDate(str_temp);
			fin >> str_temp; couponleg_paydate[i] = CDate(str_temp);

			fin >> couponleg_notional[i];
			fin >> couponleg_couponrate[i];
			fin >> couponleg_spread[i];

			couponleg_couponrate[i] = couponleg_couponrate[i] + couponleg_spread[i];
		}

	}

	if (fundingleg_cfgen_flag)
	{
		num_fundingleg_cf = findnumschedule(start_date, maturity, fundingleg_stub, fundingleg_direction, fundingleg_frq);

		fundingleg_calc_startdate = vector<CDate>(num_fundingleg_cf);

		calculationstartschedule
		(
			start_date
			, maturity
			, fundingleg_holidays
			, num_fundinglegholidays
			, fundingleg_stub
			, fundingleg_direction
			, fundingleg_conv
			, fundingleg_freq
			, fundingleg_adjflag
			, num_fundingleg_cf
			, fundingleg_calc_startdate
		);

		fundingleg_calc_enddate = vector<CDate>(num_fundingleg_cf);

		calculationendschedule
		(
			start_date
			, maturity
			, fundingleg_holidays
			, num_fundinglegholidays
			, fundingleg_stub
			, fundingleg_direction
			, fundingleg_conv
			, fundingleg_freq
			, fundingleg_adjflag
			, num_fundingleg_cf
			, fundingleg_calc_enddate
		);

		fundingleg_paydate = vector<CDate>(num_fundinglegholidays);

		paymentschedule
		(
			start_date
			, maturity
			, fundingleg_holidays
			, num_fundinglegholidays
			, fundingleg_stub
			, fundingleg_direction
			, fundingleg_conv
			, fundingleg_freq
			, num_fundingleg_cf
			, fundingleg_paydate			
		);

		fixing_date = vector<CDate>(num_fundinglegholidays);

		fixingschedule
		(
			start_date
			, maturity
			, fundingleg_holidays
			, num_fundinglegholidays
			, fixing_holidays
			, num_fixingholidays
			, fundingleg_stub
			, fundingleg_direction
			, fundingleg_conv
			, fundingleg_freq
			, fundingleg_adjflag
			, fixing_setin
			, fixing_lag
			, num_fundingleg_cf
			, fixing_date			
		);

		fundingleg_notional = vector<double>(num_fundingleg_cf);
		fundingleg_mult = vector<double>(num_fundingleg_cf);
		fundingleg_spread = vector<double>(num_fundingleg_cf);

		for (i=0; i<num_fundingleg_cf; i++)
		{
			fundingleg_notional[i] = notional;
			fundingleg_mult[i] = fundingleg_indexmulti;
			fundingleg_spread[i] = fundingleg_margin;
		}

	}
	else
	{
		fin >> label; fin >> num_fundingleg_cf;
		fixing_date = vector<CDate>(num_fundingleg_cf);
		fundingleg_calc_startdate = vector<CDate>(num_fundingleg_cf);
		fundingleg_calc_enddate = vector<CDate>(num_fundingleg_cf);
		fundingleg_paydate = vector<CDate>(num_fundingleg_cf);
		fundingleg_notional = vector<double>(num_fundingleg_cf);
		fundingleg_mult = vector<double>(num_fundingleg_cf);
		fundingleg_spread = vector<double>(num_fundingleg_cf);

		for (i=0; i<num_fundingleg_cf; i++)
		{
			fin >> str_temp; fixing_date[i] = CDate(str_temp);
			fin >> str_temp; fundingleg_calc_startdate[i] = CDate(str_temp);
			fin >> str_temp; fundingleg_calc_enddate[i] = CDate(str_temp);
			fin >> str_temp; fundingleg_paydate[i] = CDate(str_temp);
			
			fin >> fundingleg_notional[i];
			fin >> fundingleg_mult[i];
			fin >> fundingleg_spread[i];
		}
	}

	if (call_schedgen_flag)
	{
		num_calldate = findnumschedule(callstart_date, callend_date, call_stub, call_direction, call_frq);
		calleffective_date = vector<CDate>(num_calldate);
		callnotice_date = vector<CDate>(num_calldate);
		callexe_fee = vector<double>(num_calldate);
		callability_flag = vector<bool>(num_calldate);

		if(num_calldate > 0)
		{
			calculationstartschedule
			(
				callstart_date
				, callend_date
				, call_holidays
				, num_callholidays
				, call_stub
				, call_direction
				, call_conv
				, call_freq
				, "ADJUSTED"
				, num_calldate
				,calleffective_date
			);
		}

		for (i=0; i<num_calldate; i++)
		{
			ShiftBusDate
			(
				calleffective_date[i]
				, call_holidays
				, num_callholidays
				, callperiod
				, callnotice_date[i]
			);

			callexe_fee[i] = callfee;
			callability_flag[i] = callable_flag;
		}
	}
	else
	{
		fin >> label; fin >> num_calldate;
		calleffective_date = vector<CDate>(num_calldate);
		callnotice_date = vector<CDate>(num_calldate);
		callexe_fee = vector<double>(num_calldate);
		callability_flag = vector<bool>(num_calldate);

		for (i=0; i<num_calldate; i++)
		{
			fin >> str_temp; callnotice_date[i] = CDate(str_temp);
			fin >> str_temp; calleffective_date[i] = CDate(str_temp);
			fin >> callexe_fee[i];
			fin >> int_tmp; 

			if(int_tmp==1)
			{
				callability_flag[i] = 1;
			}
			else
			{
				callability_flag[i] = 0;
			}
		}
	}

	vector<int> aftercalldate_fundingstarti(num_calldate,0), aftercalldate_couponstarti(num_calldate,0);

	for (i=0; i<num_calldate; i++)
	{
		while (callnotice_date[i] > fundingleg_calc_startdate[aftercalldate_fundingstarti[i]])
		{
			aftercalldate_fundingstarti[i] = aftercalldate_fundingstarti[i]+1;

			if (aftercalldate_fundingstarti[i]>=num_fundingleg_cf-1)
			{
				break;
			}
		}

		while (callnotice_date[i] > couponleg_calc_startdate[aftercalldate_couponstarti[i]])
		{
			aftercalldate_couponstarti[i] = aftercalldate_couponstarti[i] + 1;

			if (aftercalldate_couponstarti[i] >= num_couponleg_cf-1)
			{
				break;
			}
		}
	}

	int callstarti = 0;

	if (num_calldate > 0)
	{
		while(today > callnotice_date[callstarti])
		{
			callstarti = callstarti + 1;
				
			if(callstarti >= num_calldate - 1)
			{
				break;
			}
		}
	}

	if (parameter_schedgen_flag)
	{

	}
	else
	{
		fin >> label; fin >> num_parameter;
		parameter_date = vector<CDate>(num_parameter);

		a = vector<double>(num_parameter);
		sig = vector<double>(num_parameter);

		for (i=0; i<num_parameter; i++)
		{
			fin >> str_temp; parameter_date[i] = CDate(str_temp);
			fin >> a[i];
			fin >> sig[i];
		}
	}

	fin >> label; fin >> num_fixing_history;
	vector<CDate> fixinghistory_date(num_fixing_history);
	vector<double> fixinghistory_rate(num_fixing_history);

	for(i=0; i<num_fixing_history; i++)
	{
		fin >> str_temp; fixinghistory_date[i] = CDate(str_temp);
		fin >> fixinghistory_rate[i];
	}

	vector<CDate> swaptionoptionmaturity_date(num_swaptionoption);

	for (i=0; i<num_swaptionoption; i++)
	{
		findmaturity
		(
			today
			, swaptionoptiontenor[i]
			, couponleg_holidays
			, num_couponlegholidays
			, "FOL"
			, swaptionoptionmaturity_date[i]
		);
	}

	string dcb = "ACT/365";
	vector<CDate> mtrty;
	vector<double> df, zero; 

	zerocurve
	(
		crcy
		, today
		, tenor
		, type
		, fixedrateleg_freq
		, floatingrateleg_freq
		, mkt_rate
		, mtrty
		, df
		, zero
	);

	vector<double> mtrty_t(int(mtrty.size()));

	for (i=0; i<int(mtrty.size()); i++)
	{
		mtrty_t[i] = cvg(today, mtrty[i], dcb);
	}

	CInterpolation zc(mtrty, zero), zct(mtrty_t, zero);

	fixingindex_maturity = vector<CDate>(num_fundingleg_cf);

	CInstrument fixingindex(fundinglegindex_crcy, fundinglegindex_type);

	string fixingindex_conv = fixingindex.get_conv();
	string fixingindex_dcb = fixingindex.get_basis();

	CDate *fixingindex_holidays = fixingindex.get_holiday();

	int num_fixingindexholidays = fixingindex.get_numholiday();

	for (i=0; i<num_fundingleg_cf; i++)
	{
		ShiftBusDate
		(
			fixing_date[i]
			, fixing_holidays
			, num_fixingholidays
			, -fix_lag
			, fixingindex_maturity[i]
		);

		findmaturity
		(
			fixingindex_maturity[i]
			, fundinglegindex_tenor
			, fixingindex_holidays
			, num_fixingindexholidays
			, fixingindex_conv
			, fixingindex_maturity[i]
		);
	}
	
	int couponleg_cf_starti = 0;

	while (today >= couponleg_paydate[couponleg_cf_starti])
	{
		couponleg_cf_starti = couponleg_cf_starti + 1;

		if(couponleg_cf_starti >= num_couponleg_cf -1)
		{
			break;
		}
	}
			
	int fundingleg_cf_starti = 0;

	while (today >= fundingleg_paydate[fundingleg_cf_starti])
	{
		fundingleg_cf_starti = fundingleg_cf_starti + 1;

		if(fundingleg_cf_starti >= num_fundingleg_cf - 1)
		{
			break;
		}
	}
			
	for (i = fundingleg_cf_starti; i < num_fundingleg_cf; i++)
	{
		if (today <= fixing_date[i])
		{
			break;
		}		
	}

	int fundinglegnonfixed_cf_starti = i;

	double ondf = 1.0, 
	settlement_date_df = exp(-zc(settlement_date)*cvg(today, settlement_date, dcb));

	if(fundinglegnonfixed_cf_starti < num_fundingleg_cf)
	{
		if(today == fixing_date[i]
		&& fixinghistory_rate[i] > MinFixingRate)
		{
			fundinglegnonfixed_cf_starti = fundinglegnonfixed_cf_starti + 1;
			CDate ondate;
			CInstrument temp(crcy, "DEPO");

			ShiftBusDate
			(today
			, temp.get_holiday()
			, temp.get_numholiday()
			, 1
			, ondate);
					
			ondf = exp(-zc(ondate)*cvg(today, ondate, dcb));
		}
	}

	settlement_date_df = settlement_date_df / ondf;

	int num_remained_fundingleg_cf = num_fundingleg_cf - fundinglegnonfixed_cf_starti;

	vector<double> _fundingleg_notional(num_remained_fundingleg_cf)
	, _fundingleg_mult(num_remained_fundingleg_cf)
	, _fundingleg_spread(num_remained_fundingleg_cf)
	, fundinglegcalcstartT(num_remained_fundingleg_cf)
	, fundinglegtau(num_remained_fundingleg_cf)
	, fundinglegT(num_remained_fundingleg_cf);

	vector<CDate> _fundingleg_calc_startdate(num_remained_fundingleg_cf)
	, _fundingleg_paydate(num_remained_fundingleg_cf);

	CDate fundingleg_calc_startdate0 = fundingleg_calc_startdate[num_fundingleg_cf - 1];

	if(fundinglegnonfixed_cf_starti < num_remained_fundingleg_cf)
	{
		fundingleg_calc_startdate0 = fundingleg_calc_startdate[fundinglegnonfixed_cf_starti];
	}

	double fundinglegcalcstartT0 = cvg(today, fundingleg_calc_startdate0, dcb);

	for (i=0; i<num_remained_fundingleg_cf; i++)
	{
		_fundingleg_notional[i] = fundingleg_notional[i+fundinglegnonfixed_cf_starti];
		_fundingleg_spread[i] = fundingleg_spread[i+fundinglegnonfixed_cf_starti];

		_fundingleg_calc_startdate[i] = fundingleg_calc_startdate[i+fundinglegnonfixed_cf_starti];

		_fundingleg_paydate[i] = fundingleg_paydate[i+fundinglegnonfixed_cf_starti];

		_fundingleg_mult[i] = fundingleg_mult[i+fundinglegnonfixed_cf_starti];

		fundinglegtau[i] = cvg(fundingleg_calc_startdate[i+fundinglegnonfixed_cf_starti]
							, fundingleg_calc_enddate[i+fundinglegnonfixed_cf_starti]
							, fundingleg_dcb);
				
		fundinglegcalcstartT[i] = cvg(today
									, fundingleg_calc_startdate[i+fundinglegnonfixed_cf_starti]
									, dcb);
				
		fundinglegT[i] = cvg(today
							, fundingleg_paydate[i+fundinglegnonfixed_cf_starti]
							, dcb);
	}
				
	int num_remained_couponleg_cf = num_couponleg_cf - couponleg_cf_starti
	, couponindexfixinghistory_starti = 0
	, couponleg_nonfixedcf_starti = couponleg_cf_starti;

	vector<vector<vector<CDate>>> targetswaptionswap_date(num_remained_couponleg_cf);

	vector<double> _couponleg_notional(num_remained_couponleg_cf)
		, _couponleg_couponrate(num_remained_couponleg_cf)
		, _couponleg_spread(num_remained_couponleg_cf)
		, _couponlegindex1_mult(num_remained_couponleg_cf)
		, couponlegtau(num_remained_couponleg_cf)
		, couponlegT(num_remained_couponleg_cf);

	vector<CDate> _couponleg_calc_startdate(num_remained_couponleg_cf)
		, _couponleg_paydate(num_remained_couponleg_cf);

	int max_nweek = 4
	, max_nmonth = 11
	, max_nquater = 100
	, max_nday = 10000;

	for (i=0; i<num_remained_couponleg_cf; i++)
	{
		_couponleg_notional[i] = couponleg_notional[i+couponleg_cf_starti];
		_couponleg_couponrate[i] = couponleg_couponrate[i+couponleg_cf_starti];
		_couponleg_spread[i] = couponleg_spread[i+couponleg_cf_starti];
		_couponleg_calc_startdate[i] = couponleg_calc_startdate[i+couponleg_cf_starti];
		_couponleg_paydate[i] = couponleg_paydate[i+couponleg_cf_starti];

		couponlegtau[i] = cvg(couponleg_calc_startdate[i+couponleg_cf_starti]
							, couponleg_calc_enddate[i+couponleg_cf_starti]
							, couponleg_dcb);

		couponlegT[i] = cvg(today
							, couponleg_paydate[i+couponleg_cf_starti]
							, dcb);		
	}

	imsi_sig = vector<double>(num_parameter);

	imsi_sig[0] = sig[0];

	for (i = 0; i< num_parameter-1; i++)
	{
		imsi_sig[i+1] = sig[i];
	}

	double floatinglegprice
		, fixedlegprice
		, atmswaprate
		, swapprice;

	MyPlainSwapPrice(crcy
		, today
		, settlement_date
		, notional
		, zc
		, couponleg_calc_startdate
		, couponleg_calc_enddate
		, couponleg_paydate
		, couponleg_dcb
		, couponleg_notional
		, couponleg_couponrate
		, fundingleg_calc_startdate
		, fundingleg_calc_enddate
		, fundingleg_paydate
		, fundingleg_dcb
		, fundingleg_notional
		, fundingleg_mult
		, fundingleg_spread
		, fixing_date
		, fixingindex_maturity
		, fixingindex_dcb
		, fixinghistory_date
		, fixinghistory_rate
		, false
		, floatinglegprice
		, fixedlegprice
		, atmswaprate
		, swapprice);

	// swap end

	// bermudan option start

	double tmpv = 0.0;

	if (num_calldate > 0)
	{
		double maxT
		, tmpdt
		, tmpAt
		, tmpBt
		, lambda1 = 0.5
		, lambda2 = 0.5;

		int Nx, maxNx, ii;

		string cont_dcb = "ACT/365";

		vector<double> couponleg_payt(num_couponleg_cf)
		, fundingleg_payt(num_fundingleg_cf)
		, callnotice_t(num_calldate)
		, parameter_t(num_parameter);

		for (i = 0; i < num_couponleg_cf; i++)
		{
			couponleg_payt[i] = cvg(today, couponleg_paydate[i], cont_dcb);
		}

		for (i = 0; i < num_calldate; i++)
		{
			callnotice_t[i] = cvg(today, callnotice_date[i], cont_dcb);
		}

		for (i = 0; i < num_fundingleg_cf; i++)
		{
			fundingleg_payt[i] = cvg(today, fundingleg_paydate[i], cont_dcb);
		}


		for (i = 0; i < num_parameter; i++)
		{
			parameter_t[i] = cvg(today, parameter_date[i], cont_dcb);
		}

		maxT = max(max(couponleg_payt[num_couponleg_cf-1]
						, fundingleg_payt[num_fundingleg_cf-1])
						, callnotice_t[num_calldate-1]);

		double stdr = quantile*sqrt(0.5*sig[num_parameter-1]*sig[num_parameter-1]/a[num_parameter-1]
									* (1.0-exp(-2.0*a[num_parameter-1]*maxT)));

		for (i = num_parameter - 1; i > 0; i--)
		{
			stdr = max(stdr, quantile*sqrt(0.5*sig[i-1]*sig[i-1]/a[i-1]
							*(1.0-exp(-2.0*a[i-1]*parameter_t[i])))
						);
		}

		Nx = int(ceil(stdr/dx));
		maxNx = 2 * Nx;

		vector<double> x(maxNx+1);

		for(i=1; i<=Nx; i++)
		{
			x[Nx+i] = i*dx;
			x[Nx-i] = -i*dx;
		}

		x[Nx] = 0.0;

		int callstarti = num_calldate;

		for(i=0; i<num_calldate; i++)
		{
			if
			(callnotice_t[i] > 0.0
			&& callability_flag[i])
			{
				callstarti = i;
				break;
			}
		}

		int callendi = -1;

		for(i=num_calldate-1; i>=callstarti; i--)
		{
			if(callability_flag[i])
			{
				callendi = i;
				break;
			}
		}

		int tempi;

		vector<double> dt, Vts, Uts, PMt;

		vector<int> Nt, aftercall_fund_starti, aftercall_coupon_starti;

		vector<vector<double>> phi, t;

		vector<vector<double>> tT, V_tT, PM_tT, xcoef_tT;

		vector<double> Vnew(maxNx+1), VV(maxNx+1);

		for ( i = callstarti; i <= callendi; i++ )
		{
			if (i == callstarti)
			{
				Nt.push_back(2*(int(ceil(callnotice_t[i]/(2.0*dT)))));
				dt.push_back(callnotice_t[i]/double(Nt[i-callstarti]));
			}
			else
			{
				Nt.push_back(2*(int(ceil((callnotice_t[i]-callnotice_t[i-1])/(2.0*dT)))));
				dt.push_back((callnotice_t[i]-callnotice_t[i-1])/double(Nt[i-callstarti]));
			}

			phi.push_back(vector<double>(Nt[i-callstarti]));
			t.push_back(vector<double>(Nt[i-callstarti]));

			for (tempi = 0; tempi < num_remained_fundingleg_cf; tempi++)
			{
				if (callnotice_date[i] <= _fundingleg_calc_startdate[tempi])
				{
					aftercall_fund_starti.push_back(tempi);
					break;
				}
			}
			for (tempi = 0; tempi < num_remained_couponleg_cf; tempi++)
			{
				if (callnotice_date[i] <= couponleg_calc_startdate[tempi])
				{
					aftercall_coupon_starti.push_back(tempi);
					break;
				}			
			}
		}

		vector<vector<double>> Vold(callendi-callstarti+1)
		, Vcoup(callendi-callstarti+1)
		, Vfund(callendi-callstarti+1);

		for (k = 0; k <= callendi - callstarti; k++)
		{
			cout << k << endl;

			int tmpfundingleg_cf_starti = 0;

			while (callnotice_date[k+callstarti] >= _fundingleg_paydate[tmpfundingleg_cf_starti])
			{
				tmpfundingleg_cf_starti = tmpfundingleg_cf_starti + 1;

				if (tmpfundingleg_cf_starti >= num_fundingleg_cf-1)
				{
					break;
				}
			}

			CDate tmpfundingleg_calc_startdate0 = _fundingleg_calc_startdate[aftercall_fund_starti[k]];

			double tmpfundinglegcalcstartT0 = cvg(today
												, tmpfundingleg_calc_startdate0
												, dcb);

			Vold[k] = vector<double>(maxNx + 1);
			Vcoup[k] = vector<double>(maxNx + 1);
			Vfund[k] = vector<double>(maxNx + 1);

			PMt.push_back(exp(-zct(callnotice_t[k+callstarti])*callnotice_t[k+callstarti]));

			Vts.push_back	(
							0.5
							* pow(	sig[k+callstarti]
									* 	(	1.0
											- exp(	-a[k+callstarti]
													* callnotice_t[k+callstarti]
												)
										)
									/ a[k+callstarti]
									,2.0
								)
							);

			Uts.push_back
			(	0.25
				* pow
					(
						sig[k+callstarti]
						,2.0
					)
				* 	(
						1.0
						- exp
							(
								-2.0
								* a[k+callstarti]
								* callnotice_t[k+callstarti]
							)
					)
				/ a[k+callstarti]
			);

			vector<CDate> tmpmtrty;
			vector<double> tmpdf, tmpzero;

			zerocurve
			(
				crcy
				, callnotice_date[k+callstarti]
				, tenor
				, type
				, fixedrateleg_freq
				, floatingrateleg_freq
				, mkt_rate
				, tmpmtrty
				, tmpdf
				, tmpzero
			);

			vector<double>
				tmpmtrty_t(int(tmpmtrty.size()))
				, tmpmtrty_tT(int(tmpmtrty.size()))
				, tmpVtT(int(tmpmtrty.size()))
				, tmpPMtT(int(tmpmtrty.size()))
				, tmpxcoef(int(tmpmtrty.size()));

			for (i = 0; i < int(tmpmtrty.size()); i++)
			{
				tmpmtrty_t[i] = cvg(today, tmpmtrty[i], dcb);
				tmpmtrty_tT[i] = tmpmtrty_t[i] - callnotice_t[k+callstarti];
				tmpPMtT[i] = exp(-zct(tmpmtrty_t[i])*tmpmtrty_t[i])/PMt[k]; 
				tmpxcoef[i] = - (
								1.0
								- exp
									(
										-a[k+callstarti]
										* tmpmtrty_tT[i]
									)
								)
								/ a[k+callstarti]; // B(t,T)
				tmpVtT[i] = pow(tmpxcoef[i],2.0); 
			}
						
			tT.push_back(tmpmtrty_tT);
			V_tT.push_back(tmpVtT);
			PM_tT.push_back(tmpPMtT);
			xcoef_tT.push_back(tmpxcoef);

			CDate ondate; 
			CInstrument tempinst(crcy, "DEPO");
			ShiftBusDate(callnotice_date[k+callstarti], tempinst.get_holiday(), tempinst.get_numholiday(), 1, ondate);
			double ont = cvg(callnotice_date[k+callstarti], ondate, dcb);

			for (i = 0; i <=maxNx; i++)
			{
				for (l = 0; l < num_mktrate; l++)
				{
					tmpzero[l] = -	(
									log(PM_tT[k][l])
									- V_tT[k][l]
									* Uts[k]
									- xcoef_tT[k][l]
									* 	(
									-Vts[k]
										-x[i]
										)
									)
									/ tT[k][l]; 
				}

				CInterpolation tmpzct(tT[k], tmpzero);
				CInterpolation tmpzc(tmpmtrty, tmpzero);
				double tmpondf = exp(-tmpzct(ont)*ont);

				PlainSwapFundinglegNonFixedPricexy
				(
					callnotice_date[k+callstarti]
					, callnotice_t[k+callstarti]
					, tmpondf
					, tmpzc
					, aftercall_fund_starti[k]
					, num_remained_fundingleg_cf
					, _fundingleg_notional
					, _fundingleg_mult
					, _fundingleg_spread
					, tmpfundingleg_calc_startdate0
					, tmpfundinglegcalcstartT0
					, _fundingleg_paydate
					, fundinglegtau
					, fundinglegT
					, Vfund[k][i]
				);

				NonFixedPricexy
				(
					callnotice_t[k+callstarti]
					, aftercall_coupon_starti[k]
					, tmpondf
					, tmpzc
					, num_remained_couponleg_cf
					, _couponleg_notional
					, _couponleg_couponrate
					, _couponleg_paydate
					, couponlegtau
					, couponlegT
					, Vcoup[k][i]
				);

				Vold[k][i] = Vcoup[k][i] - Vfund[k][i];							
			}
		}

		for (k = callendi-callstarti; k >= 0; k--)
		{
			tmpdt = dt[k];

			if (k == callendi - callstarti)
			{
				for (ii = 0; ii <= maxNx; ii++)
				{
					VV[ii] = max(Vold[k][ii],0.0);
				}
			}
			else
			{
				for (ii = 0; ii <= maxNx; ii++)
				{
					VV[ii] = max(Vold[k][ii],VV[ii]);
				}
			}
						
			for (l = Nt[k]-1; l >= 0; l--)
			{
				t[k][l] = callnotice_t[k+callstarti] - dt[k]*(Nt[k] - l);

				phi[k][l] = (
								pow
								(
									sig[k+callstarti]
									/ a[k+callstarti]
									, 2.0
								)
							)
							*
							(
								0.5
								- exp
									(
										-a[k+callstarti]
										* t[k][l]
									)
								+ 0.5
								* exp
									(
										-2.0
										* a[k+callstarti]
										* t[k][l]
									)
							);
				vector<double> A(maxNx+1)
				, B(maxNx+1)
				, C(maxNx+1)
				, xl(maxNx+1)
				, xd(maxNx+1)
				, xu(maxNx+1);

				tmpAt = -0.5
						* sig[k+callstarti]
						* sig[k+callstarti]
						/ (dx*dx);
				tmpBt = 1.0 / tmpdt - 2.0*tmpAt+phi[k][l];
				vector<double> ww(maxNx+1), zz(maxNx+1), xx(maxNx+1);
				for (ii = 1; ii < maxNx; ii++)
				{
					A[ii] = tmpAt - (ii-Nx)*dx*a[k+callstarti]*0.5/dx;
					B[ii] = tmpBt - (ii-Nx)*dx;
					C[ii] = tmpAt + (ii-Nx)*dx*a[k+callstarti]*0.5/dx;
				}

				B[1] = B[1] + 2.0*C[1];
				A[1] = A[1] - C[1];
				B[maxNx-1] = 2.0*A[maxNx-1]+B[maxNx-1];
				C[maxNx-1] = C[maxNx-1]-A[maxNx-1];
				xu = C;
				xd[maxNx-1] = B[maxNx-1];
				ww[maxNx-1] = VV[maxNx-1]/tmpdt;

				for (ii = maxNx-2; ii >= 1; ii--)
				{
					xl[ii] = A[ii]/xd[ii+1];
					xd[ii] = B[ii]-xl[ii]*xu[ii+1];
					ww[ii] = VV[ii]/tmpdt - xl[ii]*ww[ii+1];
				}
				xx[1] = ww[1] / xd[2];
				for (ii = 2; ii < maxNx; ii++)
				{
					xx[ii] = (ww[ii] - xu[ii]*xx[ii-1])/xd[ii];
				}
				
				xx[0] = 2.0*xx[1] - xx[2];
				xx[maxNx] = 2.0*xx[maxNx-1]-xx[maxNx-2];
				VV = xx;				
			}	
		}
		tmpv = VV[Nx];
	}
	// Bermudan Option End

	ofstream  fout(argv[2]);

	vector<double> value(6, 0.0);

	value[0] = floatinglegprice;

	value[1] = -fixedlegprice;

	value[2] = swapprice;

	value[3] = tmpv;

	value[4] = value[2] + value[3];

	value[5] =	- value[1]
				- value[3]
				+ _couponleg_notional[num_remained_couponleg_cf - 1]
				/ ondf
				/ settlement_date_df
				* exp(-zc(maturity) * cvg(today, maturity, dcb));

	fout << "Price_" << endl;
	fout << setprecision(15) << value[0] << endl;

	fout << "Price_" << endl;
	fout << setprecision(15) << value[1] << endl;

	fout << "Price_" << endl;
	fout << setprecision(15) << value[2] << endl;
	
	fout << "Price_" << endl;
	fout << setprecision(15) << value[3] << endl;

	fout << "Price_" << endl;	
	fout << setprecision(15) << value[4] << endl;

	fout << "Price_Note" << endl;
	fout << setprecision(15) << value[5] << endl;

	delete[]label;

}





