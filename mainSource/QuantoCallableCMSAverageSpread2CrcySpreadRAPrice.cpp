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

/**
 *  2017/02/13 made by yeojin
 *  [history]
 *   1. use the HW2F function to generate both USD30Y and USD2Y index, / avg spread
 *   2. tempprice=tempprice+로 되어 있지않고 위치가 잘못 지정되어 있는 부분 수정
 */
void main(int argc, char** argv)
{
	cout << "QuantoCallableCMSAverageSpread2CrcySpreadRAPrice Module - 20180124" << endl;

	int i, j, k, l;
	string dcb = "ACT/365";

	char* label = new char[100];

	ifstream fin;
	fin.open(argv[1]);

	string str_temp, couponleg_conv, couponleg_freq, couponleg_adjflag, couponleg_dcb, couponleg_stub, couponleg_direction, couponleg_holiday, fundingleg_conv, fundingleg_freq, fundingleg_adjflag, fundingleg_dcb, fundingleg_stub, fundingleg_direction, fundingleg_holiday, fundinglegindex_type, fundinglegindex_tenor, fundinglegindexfixed_freq, fundinglegindexfloating_freq, fixing_setin, fixing_holiday, fixing_lag, couponlegindex1_tenor, couponlegindex1_type, couponlegindex1fixed_freq, couponlegindex1floating_freq, couponlegindex1_holiday, couponlegindex2_tenor, couponlegindex2_type, couponlegindex2fixed_freq, couponlegindex2floating_freq, couponlegindex2_holiday;
	string couponlegrngchkindex_tenor, couponlegrngchkindex_type, couponlegrngchkindexfixed_freq, couponlegrngchkindexfloating_freq, couponlegrngchkindex_holiday, couponlegrngchkindex2_tenor, couponlegrngchkindex2_type, couponlegrngchkindex2fixed_freq, couponlegrngchkindex2floating_freq, couponlegrngchkindex2_holiday, call_freq, call_conv, call_stub, call_direction, call_holiday, paramtermshcedfrom, calibration_type, raoptm_period1, raoptm_period2, raoptm_freq1, raoptm_freq2, avg_freq, avg_weekday, avg_sched_from, avg_fixing_setin, avg_fixing_lag, avg_fixing_conv;
	double notional, fxspot, couponleg_fixedrate, couponleg_margin, fundingleg_margin, fundingleg_indexmulti, couponlegindex1_multi, couponlegindex2_multi, couponlegrngchkindex_multi, couponlegrngchkindex2_multi, caprate, floorrate, rahighbdry, ralowbdry, rahighbdryoh, ralowbdryoh, centervol, volinterval, callfee, input_corr, nu2, a1, sig1, a2, b2, sig2, eta2, CX2x2, CX2y2, gamx1x2, gamx1y2, rho2, deltabumpingsize, gammabumpingsize, crossgammabumpingsize, swaptionvegabumpingsize, capvegabumpingsize, dT, dx, dy, quantile, npath, avg_interval;
	int int_tmp, couponleg_frq, fundingleg_frq, fundinglegindx_tenor, couponlegindx1_tenor, couponlegindex1fixed_frq, couponlegindx2_tenor, couponlegindex2fixed_frq, couponlegrngchkindx_tenor, couponlegrngchkindexfixed_frq, couponlegrngchkindx2_tenor, couponlegrngchkindex2fixed_frq, fix_lag, rahighbdryoh_flag, ralowbdryoh_flag, nshift, nlockout, call_frq, callperiod, num_atmcapvolstrikes, num_mktrate, num_indexmktrate, num_rngchkindexmktrate, num_swaptionoption, num_swaptionswap, num_mktcapvolstrikes, num_capkeytenor, num_terminalcorrtenor, raoptm_prd1, raoptm_prd2, raoptm_frq1, raoptm_frq2, num_couponleg_cf, num_fundingleg_cf, num_calldate, num_paramterm, num_fixing_history, num_couponindexfixing_history;
	bool coupleg_cfgen_flag, fundingleg_cfgen_flag, rahighbdry_flag, ralowbdry_flag, rafirstdatein_flag, ralastdatein_flag, callable_flag, volsurfgen_flag, call_schedgen_flag, paramterm_schedgen_flag, calibration_flag, deltaflag, gammaflag, crossgammaflag, swaptionvegaflag, capvegaflag, flag_temp;
	vector<int> ymd;

	fin >> label;
	fin >> str_temp;
	CDate today(str_temp);
	fin >> label; fin >> str_temp; CDate start_date(str_temp);
	fin >> label; fin >> str_temp; CDate maturity(str_temp);
	fin >> label; fin >> str_temp; CDate settlement_date(str_temp);
	fin >> label; fin >> notional;
	fin >> label; fin >> str_temp; CCurrency crcy(str_temp);
	fin >> label; fin >> str_temp; CCurrency stlmcrcy(str_temp);
	fin >> label; fin >> fxspot;

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
	fin >> label; fin >> fundinglegindex_tenor; ymd = YMD2I(fundinglegindex_tenor); fundinglegindx_tenor = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> fundinglegindex_type;
	fin >> label; fin >> fundinglegindexfixed_freq;
	fin >> label; fin >> fundinglegindexfloating_freq;
	fin >> label; fin >> fixing_setin;
	fin >> label; fin >> fixing_lag; ymd = YMD2I(fixing_lag); fix_lag = ymd[2];
	fin >> label; fin >> fixing_holiday;

	fin >> label; fin >> str_temp; CCurrency couponlegindex1_crcy(str_temp);
	fin >> label; fin >> couponlegindex1_tenor; ymd = YMD2I(couponlegindex1_tenor); couponlegindx1_tenor = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegindex1_type;
	fin >> label; fin >> couponlegindex1fixed_freq; ymd = YMD2I(couponlegindex1fixed_freq); couponlegindex1fixed_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegindex1floating_freq;
	fin >> label; fin >> couponlegindex1_holiday; CDate* couponlegindex1_holidays = Holidays(couponlegindex1_holiday); int num_couponlegindex1holidays = NumHolidays(couponlegindex1_holiday);
	fin >> label; fin >> couponlegindex1_multi;

	fin >> label; fin >> str_temp; CCurrency couponlegindex2_crcy(str_temp);
	fin >> label; fin >> couponlegindex2_tenor; ymd = YMD2I(couponlegindex2_tenor); couponlegindx2_tenor = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegindex2_type;
	fin >> label; fin >> couponlegindex2fixed_freq; ymd = YMD2I(couponlegindex2fixed_freq); couponlegindex2fixed_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegindex2floating_freq;
	fin >> label; fin >> couponlegindex2_holiday; CDate* couponlegindex2_holidays = Holidays(couponlegindex2_holiday); int num_couponlegindex2holidays = NumHolidays(couponlegindex2_holiday);
	fin >> label; fin >> couponlegindex2_multi;

	fin >> label; fin >> caprate;			//not global					
	fin >> label; fin >> floorrate;			//not global					

	fin >> label; fin >> str_temp; CCurrency couponlegrngchkindex_crcy(str_temp);
	fin >> label; fin >> couponlegrngchkindex_tenor; ymd = YMD2I(couponlegrngchkindex_tenor); couponlegrngchkindx_tenor = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegrngchkindex_type;
	fin >> label; fin >> couponlegrngchkindexfixed_freq; ymd = YMD2I(couponlegrngchkindexfixed_freq); couponlegrngchkindexfixed_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegrngchkindexfloating_freq;
	fin >> label; fin >> couponlegrngchkindex_holiday; CDate* couponlegrngchkindex_holidays = Holidays(couponlegrngchkindex_holiday); int num_couponlegrngchkindexholidays = NumHolidays(couponlegrngchkindex_holiday);
	fin >> label; fin >> couponlegrngchkindex_multi;

	fin >> label; fin >> str_temp; CCurrency couponlegrngchkindex2_crcy(str_temp);
	fin >> label; fin >> couponlegrngchkindex2_tenor; ymd = YMD2I(couponlegrngchkindex2_tenor); couponlegrngchkindx2_tenor = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegrngchkindex2_type;
	fin >> label; fin >> couponlegrngchkindex2fixed_freq; ymd = YMD2I(couponlegrngchkindex2fixed_freq); couponlegrngchkindex2fixed_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> couponlegrngchkindex2floating_freq;
	fin >> label; fin >> couponlegrngchkindex2_holiday; CDate* couponlegrngchkindex2_holidays = Holidays(couponlegrngchkindex2_holiday); int num_couponlegrngchkindex2holidays = NumHolidays(couponlegrngchkindex2_holiday);
	fin >> label; fin >> couponlegrngchkindex2_multi;

	fin >> label; fin >> rahighbdry;
	fin >> label; fin >> rahighbdry_flag;
	fin >> label; fin >> ralowbdry;
	fin >> label; fin >> ralowbdry_flag;
	fin >> label; fin >> rahighbdryoh;
	fin >> label; fin >> rahighbdryoh_flag;
	fin >> label; fin >> ralowbdryoh;
	fin >> label; fin >> ralowbdryoh_flag;

	fin >> label; fin >> avg_freq;
	fin >> label; fin >> avg_interval;
	fin >> label; fin >> avg_weekday;
	fin >> label; fin >> avg_sched_from;
	fin >> label; fin >> avg_fixing_setin;
	fin >> label; fin >> avg_fixing_lag;
	fin >> label; fin >> avg_fixing_conv;

	fin >> label; fin >> str_temp; ymd = YMD2I(str_temp); nshift = ymd[2];
	fin >> label; fin >> str_temp; ymd = YMD2I(str_temp); nlockout = ymd[2];
	fin >> label; fin >> rafirstdatein_flag;
	fin >> label; fin >> ralastdatein_flag;

	fin >> label; fin >> callable_flag;
	fin >> label; fin >> str_temp; CDate callstart_date(str_temp);
	fin >> label; fin >> str_temp; CDate callend_date(str_temp);
	fin >> label; fin >> call_freq; ymd = YMD2I(call_freq); call_frq = ymd[0] * Nummonthayear + ymd[1];
	fin >> label; fin >> call_conv;
	fin >> label; fin >> call_stub;
	fin >> label; fin >> call_direction;
	fin >> label; fin >> call_holiday;
	fin >> label; fin >> str_temp; ymd = YMD2I(str_temp); callperiod = ymd[2];
	fin >> label; fin >> str_temp; CDate call_date(str_temp);
	fin >> label; fin >> callfee;
	fin >> label; fin >> call_schedgen_flag;

	fin >> label; fin >> volsurfgen_flag;
	fin >> label; fin >> centervol;
	fin >> label; fin >> volinterval;
	fin >> label; fin >> num_atmcapvolstrikes;

	fin >> label; fin >> input_corr;
	fin >> label; fin >> nu2;
	fin >> label; fin >> a1;
	fin >> label; fin >> sig1;
	fin >> label; fin >> a2;
	fin >> label; fin >> b2;
	fin >> label; fin >> sig2;
	fin >> label; fin >> eta2;
	fin >> label; fin >> CX2x2;
	fin >> label; fin >> CX2y2;
	fin >> label; fin >> gamx1x2;
	fin >> label; fin >> gamx1y2;
	fin >> label; fin >> rho2;
	fin >> label; fin >> paramtermshcedfrom;
	fin >> label; fin >> paramterm_schedgen_flag;
	fin >> label; fin >> calibration_flag;
	fin >> label; fin >> calibration_type;

	fin >> label; fin >> deltaflag;
	fin >> label; fin >> deltabumpingsize;
	fin >> label; fin >> gammaflag;
	fin >> label; fin >> gammabumpingsize;
	fin >> label; fin >> crossgammaflag;
	fin >> label; fin >> crossgammabumpingsize;
	fin >> label; fin >> swaptionvegaflag;
	fin >> label; fin >> swaptionvegabumpingsize;
	fin >> label; fin >> capvegaflag;
	fin >> label; fin >> capvegabumpingsize;

	fin >> label; fin >> num_mktrate;
	vector<string> tenor(num_mktrate), type(num_mktrate), fixedrateleg_freq(num_mktrate), floatingrateleg_freq(num_mktrate);
	vector<double> mkt_rate(num_mktrate);
	for (i = 0; i < num_mktrate; i++) { fin >> tenor[i]; fin >> type[i]; fin >> fixedrateleg_freq[i]; fin >> floatingrateleg_freq[i]; fin >> mkt_rate[i]; }

	fin >> label; fin >> num_indexmktrate;
	vector<string> indextenor(num_indexmktrate), indextype(num_indexmktrate), indexfixedrateleg_freq(num_indexmktrate), indexfloatingrateleg_freq(num_indexmktrate);
	vector<double> indexmkt_rate(num_indexmktrate);
	for (i = 0; i < num_indexmktrate; i++) { fin >> indextenor[i]; fin >> indextype[i]; fin >> indexfixedrateleg_freq[i]; fin >> indexfloatingrateleg_freq[i]; fin >> indexmkt_rate[i]; }

	num_rngchkindexmktrate = num_mktrate;
	vector<string> rngchkindextenor = tenor, rngchkindextype = type, rngchkindexfixedrateleg_freq = fixedrateleg_freq, rngchkindexfloatingrateleg_freq = floatingrateleg_freq;
	vector<double> rngchkindexmkt_rate = mkt_rate;


	fin >> label; fin >> num_swaptionoption;
	vector<string> swaptionoptiontenor(num_swaptionoption);
	for (i = 0; i < num_swaptionoption; i++) fin >> swaptionoptiontenor[i];
	fin >> label; fin >> num_swaptionswap;
	vector<string> swaptionswaptenor(num_swaptionswap);
	for (i = 0; i < num_swaptionswap; i++) fin >> swaptionswaptenor[i];
	vector<vector<double>> mktswaptionvol(num_swaptionswap);
	fin >> label;
	for (i = 0; i < num_swaptionswap; i++)
	{
		mktswaptionvol[i] = vector<double>(num_swaptionoption);
		for (j = 0; j < num_swaptionoption; j++) fin >> mktswaptionvol[i][j];
	}

	fin >> label; fin >> num_capkeytenor;
	vector<string> capkeytenor(num_capkeytenor), capfreq(num_capkeytenor);
	vector<double> capatmvol(num_capkeytenor);
	for (i = 0; i < num_capkeytenor; i++) { fin >> capkeytenor[i]; fin >> capfreq[i]; fin >> capatmvol[i]; }
	fin >> label; fin >> num_mktcapvolstrikes;
	vector<double> mktcapstrikes(num_mktcapvolstrikes);
	for (i = 0; i < num_mktcapvolstrikes; i++) fin >> mktcapstrikes[i];
	vector<vector<double>> mktcapvol(num_capkeytenor);
	fin >> label;
	for (i = 0; i < num_capkeytenor; i++)
	{
		mktcapvol[i] = vector<double>(num_mktcapvolstrikes);
		for (j = 0; j < num_mktcapvolstrikes; j++) fin >> mktcapvol[i][j];
	}

	fin >> label; fin >> num_terminalcorrtenor;
	vector<string> terminalcorrtenor(num_terminalcorrtenor);
	vector<vector<double>> terminalcorr(num_terminalcorrtenor);
	for (i = 0; i < num_terminalcorrtenor; i++)
	{
		fin >> terminalcorrtenor[i];
		terminalcorr[i] = vector<double>(num_terminalcorrtenor);
		for (j = 0; j < num_terminalcorrtenor; j++) { fin >> label; fin >> terminalcorr[i][j]; }
	}

	fin >> label; fin >> raoptm_period1; ymd = YMD2I(raoptm_period1); raoptm_prd1 = ymd[0] * Nummonthayear + ymd[1]; fin >> raoptm_freq1; ymd = YMD2I(raoptm_freq1); raoptm_frq1 = ymd[2];
	fin >> label; fin >> raoptm_period2; ymd = YMD2I(raoptm_period2); raoptm_prd2 = ymd[0] * Nummonthayear + ymd[1]; fin >> raoptm_freq2; ymd = YMD2I(raoptm_freq2); raoptm_frq2 = ymd[2];

	fin >> label; fin >> dT;
	fin >> label; fin >> dx;
	fin >> label; fin >> dy;
	fin >> label; fin >> quantile;
	fin >> label; fin >> npath;

	CDate* couponleg_holidays = Holidays(couponleg_holiday), * fundingleg_holidays = Holidays(fundingleg_holiday), * fixing_holidays = Holidays(fixing_holiday), * call_holidays = Holidays(call_holiday);
	int num_couponlegholidays = NumHolidays(couponleg_holiday), num_fundinglegholidays = NumHolidays(fundingleg_holiday), num_callhoidays = NumHolidays(call_holiday), num_fixingholidays = NumHolidays(fixing_holiday);
	vector<CDate> couponleg_calc_startdate, couponleg_calc_enddate, couponleg_paydate, fundingleg_calc_startdate, fundingleg_calc_enddate, fundingleg_paydate, fixing_date, fixingindex_maturity, calleffective_date, callnotice_date, paramterm_date;
	vector<double> couponleg_notional, couponleg_couponrate, couponleg_spread, couponlegindex1_mult, couponlegindex2_mult, couponlegrngchkindex_mult, couponlegrngchkindex2_mult, caprates, floorrates, rahigh_bdry, ralow_bdry, rahigh_bdryoh, ralow_bdryoh, fundingleg_notional, fundingleg_mult, fundingleg_spread, callexe_fee, nu2s, a1s, sig1s, a2s, b2s, sig2s, eta2s, CX2x2s, CX2y2s, gamx1x2s, gamx1y2s, rho2s;
	vector<bool> ra_flag, rahigh_bdryin_flag, ralow_bdryin_flag, callability_flag;
	vector<int> couponlegindx1_tenors, couponlegindx2_tenors, couponlegrngchkindx_tenors, couponlegrngchkindx2_tenors, couponlegindex1fixed_frqs, couponlegindex2fixed_frqs, couponlegrngchkindexfixed_frqs, couponlegrngchkindex2fixed_frqs, rahigh_bdryoh_flag, ralow_bdryoh_flag, num_couponlegindex1holidayss, num_couponlegindex2holidayss, num_couponlegrngchkindexholidayss, num_couponlegrngchkindex2holidayss;
	vector<CCurrency> couponlegindex1_crcys, couponlegindex2_crcys, couponlegrngchkindex_crcys, couponlegrngchkindex2_crcys;
	vector<string> couponlegindex1_tenors, couponlegindex2_tenors, couponlegrngchkindex_tenors, couponlegrngchkindex2_tenors, couponlegindex1_types, couponlegindex2_types, couponlegrngchkindex_types, couponlegrngchkindex2_types, couponlegindex1fixed_freqs, couponlegindex2fixed_freqs, couponlegrngchkindexfixed_freqs, couponlegrngchkindex2fixed_freqs, couponlegindex1floating_freqs, couponlegrngchkindexfloating_freqs, couponlegrngchkindex2floating_freqs, couponlegindex2floating_freqs;
	vector<CDate*> couponlegindex1_holidayss, couponlegindex2_holidayss, couponlegrngchkindex_holidayss, couponlegrngchkindex2_holidayss;
	if (coupleg_cfgen_flag)
	{
		num_couponleg_cf = findnumschedule(start_date, maturity, couponleg_stub, couponleg_direction, couponleg_frq);
		couponleg_calc_startdate = vector<CDate>(num_couponleg_cf);
		calculationstartschedule(start_date, maturity, couponleg_holidays, num_couponlegholidays, couponleg_stub, fundingleg_direction, couponleg_conv, couponleg_freq, couponleg_adjflag, num_couponleg_cf, couponleg_calc_startdate);
		couponleg_calc_enddate = vector<CDate>(num_couponleg_cf);
		calculationendschedule(start_date, maturity, couponleg_holidays, num_couponlegholidays, couponleg_stub, fundingleg_direction, couponleg_conv, couponleg_freq, couponleg_adjflag, num_couponleg_cf, couponleg_calc_enddate);
		couponleg_paydate = vector<CDate>(num_couponleg_cf);
		paymentschedule(start_date, maturity, couponleg_holidays, num_couponlegholidays, couponleg_stub, fundingleg_direction, couponleg_conv, couponleg_freq, num_couponleg_cf, couponleg_paydate);
		couponleg_notional = vector<double>(num_couponleg_cf);
		couponleg_couponrate = vector<double>(num_couponleg_cf);
		couponleg_spread = vector<double>(num_couponleg_cf);
		ra_flag = vector<bool>(num_couponleg_cf, true);
		couponlegindex1_crcys = vector<CCurrency>(num_couponleg_cf, couponlegindex1_crcy);
		couponlegindex1_tenors = vector<string>(num_couponleg_cf, couponlegindex1_tenor);
		couponlegindx1_tenors = vector<int>(num_couponleg_cf, couponlegindx1_tenor);
		couponlegindex1_types = vector<string>(num_couponleg_cf, couponlegindex1_type);
		couponlegindex1fixed_freqs = vector<string>(num_couponleg_cf, couponlegindex1fixed_freq);
		couponlegindex1fixed_frqs = vector<int>(num_couponleg_cf, couponlegindex1fixed_frq);
		couponlegindex1floating_freqs = vector<string>(num_couponleg_cf, couponlegindex1floating_freq);
		couponlegindex1_holidayss = vector<CDate*>(num_couponleg_cf, couponlegindex1_holidays);
		num_couponlegindex1holidayss = vector<int>(num_couponleg_cf, num_couponlegindex1holidays);
		couponlegindex1_mult = vector<double>(num_couponleg_cf);
		couponlegindex2_crcys = vector<CCurrency>(num_couponleg_cf, couponlegindex2_crcy);
		couponlegindex2_tenors = vector<string>(num_couponleg_cf, couponlegindex2_tenor);
		couponlegindx2_tenors = vector<int>(num_couponleg_cf, couponlegindx2_tenor);
		couponlegindex2_types = vector<string>(num_couponleg_cf, couponlegindex2_type);
		couponlegindex2fixed_freqs = vector<string>(num_couponleg_cf, couponlegindex2fixed_freq);
		couponlegindex2fixed_frqs = vector<int>(num_couponleg_cf, couponlegindex2fixed_frq);
		couponlegindex2floating_freqs = vector<string>(num_couponleg_cf, couponlegindex2floating_freq);
		couponlegindex2_holidayss = vector<CDate*>(num_couponleg_cf, couponlegindex2_holidays);
		num_couponlegindex2holidayss = vector<int>(num_couponleg_cf, num_couponlegindex2holidays);
		couponlegindex2_mult = vector<double>(num_couponleg_cf);

		couponlegrngchkindex_crcys = vector<CCurrency>(num_couponleg_cf, couponlegrngchkindex_crcy);
		couponlegrngchkindex_tenors = vector<string>(num_couponleg_cf, couponlegrngchkindex_tenor);
		couponlegrngchkindx_tenors = vector<int>(num_couponleg_cf, couponlegrngchkindx_tenor);
		couponlegrngchkindex_types = vector<string>(num_couponleg_cf, couponlegrngchkindex_type);
		couponlegrngchkindexfixed_freqs = vector<string>(num_couponleg_cf, couponlegrngchkindexfixed_freq);
		couponlegrngchkindexfixed_frqs = vector<int>(num_couponleg_cf, couponlegrngchkindexfixed_frq);
		couponlegrngchkindexfloating_freqs = vector<string>(num_couponleg_cf, couponlegrngchkindexfloating_freq);
		couponlegrngchkindex_holidayss = vector<CDate*>(num_couponleg_cf, couponlegrngchkindex_holidays);
		num_couponlegrngchkindexholidayss = vector<int>(num_couponleg_cf, num_couponlegrngchkindexholidays);
		couponlegrngchkindex_mult = vector<double>(num_couponleg_cf);

		couponlegrngchkindex2_crcys = vector<CCurrency>(num_couponleg_cf, couponlegrngchkindex2_crcy);
		couponlegrngchkindex2_tenors = vector<string>(num_couponleg_cf, couponlegrngchkindex2_tenor);
		couponlegrngchkindx2_tenors = vector<int>(num_couponleg_cf, couponlegrngchkindx2_tenor);
		couponlegrngchkindex2_types = vector<string>(num_couponleg_cf, couponlegrngchkindex2_type);
		couponlegrngchkindex2fixed_freqs = vector<string>(num_couponleg_cf, couponlegrngchkindex2fixed_freq);
		couponlegrngchkindex2fixed_frqs = vector<int>(num_couponleg_cf, couponlegrngchkindex2fixed_frq);
		couponlegrngchkindex2floating_freqs = vector<string>(num_couponleg_cf, couponlegrngchkindex2floating_freq);
		couponlegrngchkindex2_holidayss = vector<CDate*>(num_couponleg_cf, couponlegrngchkindex2_holidays);
		num_couponlegrngchkindex2holidayss = vector<int>(num_couponleg_cf, num_couponlegrngchkindex2holidays);
		couponlegrngchkindex2_mult = vector<double>(num_couponleg_cf);

		caprates = vector<double>(num_couponleg_cf, caprate);
		floorrates = vector<double>(num_couponleg_cf, floorrate);

		rahigh_bdry = vector<double>(num_couponleg_cf);
		rahigh_bdryin_flag = vector<bool>(num_couponleg_cf);
		rahigh_bdryoh = vector<double>(num_couponleg_cf);
		rahigh_bdryoh_flag = vector<int>(num_couponleg_cf);
		ralow_bdry = vector<double>(num_couponleg_cf);
		ralow_bdryin_flag = vector<bool>(num_couponleg_cf);
		ralow_bdryoh = vector<double>(num_couponleg_cf);
		ralow_bdryoh_flag = vector<int>(num_couponleg_cf);
		for (i = 0; i < num_couponleg_cf; i++)
		{
			couponleg_notional[i] = notional;
			couponleg_couponrate[i] = couponleg_fixedrate;
			couponleg_spread[i] = couponleg_margin;
			couponlegindex1_mult[i] = couponlegindex1_multi;
			couponlegindex2_mult[i] = couponlegindex2_multi;
			rahigh_bdry[i] = rahighbdry;
			rahigh_bdryin_flag[i] = rahighbdry_flag;
			rahigh_bdryoh[i] = rahighbdryoh;
			rahigh_bdryoh_flag[i] = rahighbdryoh_flag;
			ralow_bdry[i] = ralowbdry;
			ralow_bdryin_flag[i] = ralowbdry_flag;
			ralow_bdryoh[i] = ralowbdryoh;
			ralow_bdryoh_flag[i] = ralowbdryoh_flag;
		}
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
		ra_flag = vector<bool>(num_couponleg_cf);
		couponlegindex1_crcys = vector<CCurrency>(num_couponleg_cf);
		couponlegindex1_tenors = vector<string>(num_couponleg_cf);
		couponlegindx1_tenors = vector<int>(num_couponleg_cf);
		couponlegindex1_types = vector<string>(num_couponleg_cf);
		couponlegindex1fixed_freqs = vector<string>(num_couponleg_cf);
		couponlegindex1fixed_frqs = vector<int>(num_couponleg_cf);
		couponlegindex1floating_freqs = vector<string>(num_couponleg_cf);
		couponlegindex1_holidayss = vector<CDate*>(num_couponleg_cf);
		num_couponlegindex1holidayss = vector<int>(num_couponleg_cf);

		couponlegindex1_mult = vector<double>(num_couponleg_cf);
		couponlegindex2_crcys = vector<CCurrency>(num_couponleg_cf);
		couponlegindex2_tenors = vector<string>(num_couponleg_cf);
		couponlegindx2_tenors = vector<int>(num_couponleg_cf);
		couponlegindex2_types = vector<string>(num_couponleg_cf);
		couponlegindex2fixed_freqs = vector<string>(num_couponleg_cf);
		couponlegindex2fixed_frqs = vector<int>(num_couponleg_cf);
		couponlegindex2floating_freqs = vector<string>(num_couponleg_cf);
		couponlegindex2_holidayss = vector<CDate*>(num_couponleg_cf);
		num_couponlegindex2holidayss = vector<int>(num_couponleg_cf);
		couponlegindex2_mult = vector<double>(num_couponleg_cf);

		couponlegrngchkindex_crcys = vector<CCurrency>(num_couponleg_cf);
		couponlegrngchkindex_tenors = vector<string>(num_couponleg_cf);
		couponlegrngchkindx_tenors = vector<int>(num_couponleg_cf);
		couponlegrngchkindex_types = vector<string>(num_couponleg_cf);
		couponlegrngchkindexfixed_freqs = vector<string>(num_couponleg_cf);
		couponlegrngchkindexfixed_frqs = vector<int>(num_couponleg_cf);
		couponlegrngchkindexfloating_freqs = vector<string>(num_couponleg_cf);
		couponlegrngchkindex_holidayss = vector<CDate*>(num_couponleg_cf);
		num_couponlegrngchkindexholidayss = vector<int>(num_couponleg_cf);
		couponlegrngchkindex_mult = vector<double>(num_couponleg_cf);

		couponlegrngchkindex2_crcys = vector<CCurrency>(num_couponleg_cf);
		couponlegrngchkindex2_tenors = vector<string>(num_couponleg_cf);
		couponlegrngchkindx2_tenors = vector<int>(num_couponleg_cf);
		couponlegrngchkindex2_types = vector<string>(num_couponleg_cf);
		couponlegrngchkindex2fixed_freqs = vector<string>(num_couponleg_cf);
		couponlegrngchkindex2fixed_frqs = vector<int>(num_couponleg_cf);
		couponlegrngchkindex2floating_freqs = vector<string>(num_couponleg_cf);
		couponlegrngchkindex2_holidayss = vector<CDate*>(num_couponleg_cf);
		num_couponlegrngchkindex2holidayss = vector<int>(num_couponleg_cf);
		couponlegrngchkindex2_mult = vector<double>(num_couponleg_cf);

		caprates = vector<double>(num_couponleg_cf);
		floorrates = vector<double>(num_couponleg_cf);

		rahigh_bdry = vector<double>(num_couponleg_cf);
		rahigh_bdryin_flag = vector<bool>(num_couponleg_cf);
		rahigh_bdryoh = vector<double>(num_couponleg_cf);
		rahigh_bdryoh_flag = vector<int>(num_couponleg_cf);
		ralow_bdry = vector<double>(num_couponleg_cf);
		ralow_bdryin_flag = vector<bool>(num_couponleg_cf);
		ralow_bdryoh = vector<double>(num_couponleg_cf);
		ralow_bdryoh_flag = vector<int>(num_couponleg_cf);
		for (i = 0; i < num_couponleg_cf; i++)
		{
			fin >> str_temp; couponleg_calc_startdate[i] = CDate(str_temp);
			fin >> str_temp; couponleg_calc_enddate[i] = CDate(str_temp);
			fin >> str_temp; couponleg_paydate[i] = CDate(str_temp);
			fin >> couponleg_notional[i];
			fin >> couponleg_couponrate[i];
			fin >> couponleg_spread[i];
			fin >> flag_temp; ra_flag[i] = flag_temp;

			fin >> str_temp; couponlegindex1_crcys[i] = CCurrency(str_temp);
			fin >> couponlegindex1_tenors[i]; ymd = YMD2I(couponlegindex1_tenors[i]); couponlegindx1_tenors[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegindex1_types[i];
			fin >> couponlegindex1fixed_freqs[i]; ymd = YMD2I(couponlegindex1fixed_freqs[i]); couponlegindex1fixed_frqs[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegindex1floating_freqs[i];
			fin >> str_temp; couponlegindex1_holidayss[i] = Holidays(str_temp); num_couponlegindex1holidayss[i] = NumHolidays(str_temp);
			fin >> couponlegindex1_mult[i];

			fin >> str_temp; couponlegindex2_crcys[i] = CCurrency(str_temp);
			fin >> couponlegindex2_tenors[i]; ymd = YMD2I(couponlegindex2_tenors[i]); couponlegindx2_tenors[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegindex2_types[i];
			fin >> couponlegindex2fixed_freqs[i]; ymd = YMD2I(couponlegindex2fixed_freqs[i]); couponlegindex2fixed_frqs[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegindex2floating_freqs[i];
			fin >> str_temp; couponlegindex2_holidayss[i] = Holidays(str_temp); num_couponlegindex2holidayss[i] = NumHolidays(str_temp);
			fin >> couponlegindex2_mult[i];

			fin >> caprates[i];
			fin >> floorrates[i];

			fin >> str_temp; couponlegrngchkindex_crcys[i] = CCurrency(str_temp);
			fin >> couponlegrngchkindex_tenors[i]; ymd = YMD2I(couponlegrngchkindex_tenors[i]); couponlegrngchkindx_tenors[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegrngchkindex_types[i];
			fin >> couponlegrngchkindexfixed_freqs[i]; ymd = YMD2I(couponlegrngchkindexfixed_freqs[i]); couponlegrngchkindexfixed_frqs[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegrngchkindexfloating_freqs[i];
			fin >> str_temp; couponlegrngchkindex_holidayss[i] = Holidays(str_temp); num_couponlegrngchkindexholidayss[i] = NumHolidays(str_temp);
			fin >> couponlegrngchkindex_mult[i];

			fin >> str_temp; couponlegrngchkindex2_crcys[i] = CCurrency(str_temp);
			fin >> couponlegrngchkindex2_tenors[i]; ymd = YMD2I(couponlegrngchkindex2_tenors[i]); couponlegrngchkindx2_tenors[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegrngchkindex2_types[i];
			fin >> couponlegrngchkindex2fixed_freqs[i]; ymd = YMD2I(couponlegrngchkindex2fixed_freqs[i]); couponlegrngchkindex2fixed_frqs[i] = ymd[0] * Nummonthayear + ymd[1];
			fin >> couponlegrngchkindex2floating_freqs[i];
			fin >> str_temp; couponlegrngchkindex2_holidayss[i] = Holidays(str_temp); num_couponlegrngchkindex2holidayss[i] = NumHolidays(str_temp);
			fin >> couponlegrngchkindex2_mult[i];

			fin >> rahigh_bdry[i];
			fin >> int_tmp; if (int_tmp == 1) rahigh_bdryin_flag[i] = 1; else rahigh_bdryin_flag[i] = 0;
			fin >> rahigh_bdryoh[i];
			fin >> rahigh_bdryoh_flag[i];
			fin >> ralow_bdry[i];
			fin >> int_tmp; if (int_tmp == 1) ralow_bdryin_flag[i] = 1; else ralow_bdryin_flag[i] = 0;
			fin >> ralow_bdryoh[i];
			fin >> ralow_bdryoh_flag[i];
		}
	}

	if (fundingleg_cfgen_flag)
	{
		num_fundingleg_cf = findnumschedule(start_date, maturity, fundingleg_stub, fundingleg_direction, fundingleg_frq);
		fundingleg_calc_startdate = vector<CDate>(num_fundingleg_cf);
		calculationstartschedule(start_date, maturity, fundingleg_holidays, num_fundinglegholidays, fundingleg_stub, fundingleg_direction, fundingleg_conv, fundingleg_freq, fundingleg_adjflag, num_fundingleg_cf, fundingleg_calc_startdate);
		fundingleg_calc_enddate = vector<CDate>(num_fundingleg_cf);
		calculationendschedule(start_date, maturity, fundingleg_holidays, num_fundinglegholidays, fundingleg_stub, fundingleg_direction, fundingleg_conv, fundingleg_freq, fundingleg_adjflag, num_fundingleg_cf, fundingleg_calc_enddate);
		fundingleg_paydate = vector<CDate>(num_fundingleg_cf);
		paymentschedule(start_date, maturity, fundingleg_holidays, num_fundinglegholidays, fundingleg_stub, fundingleg_direction, fundingleg_conv, fundingleg_freq, num_fundingleg_cf, fundingleg_paydate);
		fixing_date = vector<CDate>(num_fundingleg_cf);
		fixingschedule(start_date, maturity, fundingleg_holidays, num_fundinglegholidays, fixing_holidays, num_fixingholidays, fundingleg_stub, fundingleg_direction, fundingleg_conv, fundingleg_freq, fundingleg_adjflag, fixing_setin, fixing_lag, num_fundingleg_cf, fixing_date);
		fundingleg_notional = vector<double>(num_fundingleg_cf);
		fundingleg_mult = vector<double>(num_fundingleg_cf);
		fundingleg_spread = vector<double>(num_fundingleg_cf);
		for (i = 0; i < num_fundingleg_cf; i++)
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
		for (i = 0; i < num_fundingleg_cf; i++)
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
		if (num_calldate > 0) calculationstartschedule(callstart_date, callend_date, call_holidays, num_callhoidays, call_stub, call_direction, call_conv, call_freq, "ADJUSTED", num_calldate, calleffective_date);
		for (i = 0; i < num_calldate; i++)
		{
			ShiftBusDate(calleffective_date[i], call_holidays, num_callhoidays, callperiod, callnotice_date[i]);
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
		for (i = 0; i < num_calldate; i++)
		{
			fin >> str_temp; callnotice_date[i] = CDate(str_temp);
			fin >> str_temp; calleffective_date[i] = CDate(str_temp);
			fin >> callexe_fee[i];
			fin >> int_tmp; if (int_tmp == 1) callability_flag[i] = 1; else callability_flag[i] = 0;
		}
	}
	int callstarti = num_calldate;
	if (num_calldate > 0)
	{
		callstarti = 0;
		while (today > callnotice_date[callstarti])
		{
			callstarti = callstarti + 1;
			if (callstarti >= num_calldate) break;
		}
	}
	double b2sadj = 3.42, eta2sadj = 0.0077, rhosadj = -0.99;
	if (paramterm_schedgen_flag)
	{
	}
	else
	{
		fin >> label; fin >> num_paramterm;
		paramterm_date = vector<CDate>(num_paramterm);
		nu2s = vector<double>(num_paramterm);
		a1s = vector<double>(num_paramterm);
		sig1s = vector<double>(num_paramterm);
		a2s = vector<double>(num_paramterm);
		b2s = vector<double>(num_paramterm);
		sig2s = vector<double>(num_paramterm);
		eta2s = vector<double>(num_paramterm);
		CX2x2s = vector<double>(num_paramterm);
		CX2y2s = vector<double>(num_paramterm);

		gamx1x2s = vector<double>(num_paramterm);
		gamx1y2s = vector<double>(num_paramterm);
		rho2s = vector<double>(num_paramterm);

		for (i = 0; i < num_paramterm; i++)
		{
			fin >> str_temp; paramterm_date[i] = CDate(str_temp);
			fin >> nu2s[i];
			fin >> a1s[i];
			fin >> sig1s[i];
			fin >> a2s[i];
			fin >> b2s[i];
			fin >> sig2s[i];
			fin >> eta2s[i];
			fin >> CX2x2s[i];
			fin >> CX2y2s[i];
			fin >> gamx1x2s[i];
			fin >> gamx1y2s[i];
			fin >> rho2s[i];
		}
	}

	fin >> label; fin >> num_fixing_history;
	vector<CDate> fixinghistory_date(num_fixing_history);
	vector<double> fixinghistory_rate(num_fixing_history);
	for (i = 0; i < num_fixing_history; i++)
	{
		fin >> str_temp; fixinghistory_date[i] = CDate(str_temp);
		fin >> fixinghistory_rate[i];
	}

	fin >> label; fin >> num_couponindexfixing_history;
	vector<CDate> couponindexfixinghistory_date(num_couponindexfixing_history);
	vector<double> couponindex1fixinghistory_rate(num_couponindexfixing_history), couponindex2fixinghistory_rate(num_couponindexfixing_history), couponrngchkindexfixinghistory_rate(num_couponindexfixing_history), couponrngchkindex2fixinghistory_rate(num_couponindexfixing_history);
	for (i = 0; i < num_couponindexfixing_history; i++)
	{
		fin >> str_temp; couponindexfixinghistory_date[i] = CDate(str_temp);
		fin >> couponindex1fixinghistory_rate[i];
		fin >> couponindex2fixinghistory_rate[i];
		fin >> couponrngchkindexfixinghistory_rate[i];
		fin >> couponrngchkindex2fixinghistory_rate[i];
	}

	vector<CDate> swaptionoptionmaturity_date(num_swaptionoption);
	for (i = 0; i < num_swaptionoption; i++) findmaturity(today, swaptionoptiontenor[i], couponleg_holidays, num_couponlegholidays, "FOL", swaptionoptionmaturity_date[i]);

	vector<CDate> mtrty, index1mtrty, index2mtrty, rngchkindexmtrty;
	vector<double> df, zero, index1df, index1zero, index2df, index2zero, rngchkindexdf, rngchkindexzero, rngchkindex2df, rngchkindex2zero, mtrtyt, index1mtrtyt, index2mtrtyt, rngchkindexmtrtyt, rngchkindex2mtrtyt;

	zerocurve(crcy, today, tenor, type, fixedrateleg_freq, floatingrateleg_freq, mkt_rate, mtrty, df, zero);
	CInterpolation zc(mtrty, zero);
	for (i = 0; i<int(mtrty.size()); i++) mtrtyt.push_back(cvg(today, mtrty[i], dcb));
	CInterpolation zct(mtrtyt, zero);

	//CInterpolation index1zc=zc;								
	//CInterpolation index1zct=zct;								
	//index1df=df;								
	//index1zero=zero;								
	//index1mtrtyt=mtrtyt;								

	CInterpolation rngchkindexzc = zc;
	CInterpolation rngchkindexzct = zct;
	rngchkindexdf = df;
	rngchkindexzero = zero;
	rngchkindexmtrtyt = mtrtyt;

	zerocurve(couponlegindex2_crcy, today, indextenor, indextype, indexfixedrateleg_freq, indexfloatingrateleg_freq, indexmkt_rate, index2mtrty, index2df, index2zero);
	CInterpolation index2zc(index2mtrty, index2zero);
	for (i = 0; i<int(index2mtrty.size()); i++) index2mtrtyt.push_back(cvg(today, index2mtrty[i], dcb));
	CInterpolation index2zct(index2mtrtyt, index2zero);

	CInterpolation index1zc = index2zc;
	CInterpolation index1zct = index2zct;
	index1df = index2df;
	index1zero = index2zero;
	index1mtrtyt = index2mtrtyt;

	CInterpolation rngchkindex2zc = index2zc;
	CInterpolation rngchkindex2zct = index2zct;
	rngchkindex2df = index2df;
	rngchkindex2zero = index2zero;
	rngchkindex2mtrtyt = index2mtrtyt;

	/*
	zerocurve(couponlegrngchkindex_crcy,today,rngchkindextenor,rngchkindextype,rngchkindexfixedrateleg_freq,rngchkindexfloatingrateleg_freq,rngchkindexmkt_rate,rngchkindexmtrty,rngchkindexdf,rngchkindexzero);
	CInterpolation rngchkindexzc(rngchkindexmtrty,rngchkindexzero);
	for(i=0;i<int(rngchkindexmtrty.size());i++) rngchkindexmtrtyt.push_back(cvg(today,rngchkindexmtrty[i],dcb));
	CInterpolation rngchkindexzct(rngchkindexmtrtyt,rngchkindexzero);
	*/

	fixingindex_maturity = vector<CDate>(num_fundingleg_cf);
	CInstrument fixingindex(fundinglegindex_crcy, fundinglegindex_type);
	string fixingindex_conv = fixingindex.get_conv();
	string fixingindex_dcb = fixingindex.get_basis();
	CDate* fixingindex_holidays = fixingindex.get_holiday();
	int num_fixingindexholidays = fixingindex.get_numholiday();
	for (i = 0; i < num_fundingleg_cf; i++)
	{
		ShiftBusDate(fixing_date[i], fixing_holidays, num_fixingholidays, -fix_lag, fixingindex_maturity[i]);
		findmaturity(fixingindex_maturity[i], fundinglegindx_tenor, fixingindex_holidays, num_fixingindexholidays, fixingindex_conv, fixingindex_maturity[i]);
	}

	CDate ra_period1_maturity = DateAdd('m', raoptm_prd1, today), ra_period2_maturity = DateAdd('m', raoptm_prd2, today);


	int couponleg_cf_starti = 0;
	while (today >= couponleg_paydate[couponleg_cf_starti])
	{
		couponleg_cf_starti = couponleg_cf_starti + 1;
		if (couponleg_cf_starti >= num_couponleg_cf - 1) break;
	}

	int fundingleg_cf_starti = 0;
	while (today >= fundingleg_paydate[fundingleg_cf_starti])
	{
		fundingleg_cf_starti = fundingleg_cf_starti + 1;
		if (fundingleg_cf_starti >= num_fundingleg_cf - 1) break;
	}

	for (i = fundingleg_cf_starti; i < num_fundingleg_cf; i++) if (today <= fixing_date[i]) break;

	int fundinglegnonfixed_cf_starti = i;
	double ondf = 1.0, settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb));
	if (fundinglegnonfixed_cf_starti < num_fundingleg_cf)
	{
		if (today == fixing_date[i] && fixinghistory_rate[i] > MinFixingRate)
		{
			fundinglegnonfixed_cf_starti = fundinglegnonfixed_cf_starti + 1;
			CDate ondate;
			CInstrument temp(crcy, "DEPO");
			ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
			ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
		}
	}
	settlement_date_df = settlement_date_df / ondf;

	int num_remained_fundingleg_cf = num_fundingleg_cf - fundinglegnonfixed_cf_starti;
	vector<double> _fundingleg_notional(num_remained_fundingleg_cf), _fundingleg_mult(num_remained_fundingleg_cf), _fundingleg_spread(num_remained_fundingleg_cf), fundinglegcalcstartT(num_remained_fundingleg_cf), fundinglegtau(num_remained_fundingleg_cf), fundinglegT(num_remained_fundingleg_cf);;
	vector<CDate> _fixing_date(num_remained_fundingleg_cf), _fundingleg_calc_startdate(num_remained_fundingleg_cf), _fundingleg_paydate(num_remained_fundingleg_cf);
	CDate fundingleg_calc_startdate0 = fundingleg_calc_startdate[num_fundingleg_cf - 1];
	if (fundinglegnonfixed_cf_starti < num_fundingleg_cf) fundingleg_calc_startdate0 = fundingleg_calc_startdate[fundinglegnonfixed_cf_starti];
	double fundinglegcalcstartT0 = cvg(today, fundingleg_calc_startdate0, dcb);
	for (i = 0; i < num_remained_fundingleg_cf; i++)
	{
		_fundingleg_notional[i] = fundingleg_notional[i + fundinglegnonfixed_cf_starti];
		_fundingleg_spread[i] = fundingleg_spread[i + fundinglegnonfixed_cf_starti];
		_fixing_date[i] = fixing_date[i + fundinglegnonfixed_cf_starti];
		_fundingleg_calc_startdate[i] = fundingleg_calc_startdate[i + fundinglegnonfixed_cf_starti];
		_fundingleg_paydate[i] = fundingleg_paydate[i + fundinglegnonfixed_cf_starti];
		_fundingleg_mult[i] = fundingleg_mult[i + fundinglegnonfixed_cf_starti];
		fundinglegtau[i] = cvg(fundingleg_calc_startdate[i + fundinglegnonfixed_cf_starti], fundingleg_calc_enddate[i + fundinglegnonfixed_cf_starti], fundingleg_dcb);
		fundinglegT[i] = cvg(today, fundingleg_paydate[i + fundinglegnonfixed_cf_starti], dcb);
	}

	int num_remained_couponleg_cf = num_couponleg_cf - couponleg_cf_starti, couponindexfixinghistory_starti = 0, couponleg_nonfixedcf_starti = couponleg_cf_starti, _couponindexfixinghistory_starti = 0;
	vector<CInstrument> index1(num_couponleg_cf), index2(num_couponleg_cf), rngchkindex(num_couponleg_cf), rngchkindex2(num_couponleg_cf);
	vector<string> index1_dcb(num_couponleg_cf), index1_stub(num_couponleg_cf), index1_direction(num_couponleg_cf), index1_conv(num_couponleg_cf), index1_adjflag(num_couponleg_cf), index1_tenors(num_couponleg_cf), index1fixed_freqs(num_couponleg_cf), index2_dcb(num_couponleg_cf), index2_stub(num_couponleg_cf), index2_direction(num_couponleg_cf), index2_conv(num_couponleg_cf), index2_adjflag(num_couponleg_cf), index2_tenors(num_couponleg_cf), index2fixed_freqs(num_couponleg_cf),
		rngchkindex_dcb(num_couponleg_cf), rngchkindex_stub(num_couponleg_cf), rngchkindex_direction(num_couponleg_cf), rngchkindex_conv(num_couponleg_cf), rngchkindex_adjflag(num_couponleg_cf), rngchkindex_tenors(num_couponleg_cf), rngchkindexfixed_freqs(num_couponleg_cf),
		rngchkindex2_dcb(num_couponleg_cf), rngchkindex2_stub(num_couponleg_cf), rngchkindex2_direction(num_couponleg_cf), rngchkindex2_conv(num_couponleg_cf), rngchkindex2_adjflag(num_couponleg_cf), rngchkindex2_tenors(num_couponleg_cf), rngchkindex2fixed_freqs(num_couponleg_cf);
	//////////////////////////									
	vector<CInstrument> _index1(num_remained_couponleg_cf), _index2(num_remained_couponleg_cf), _rngchkindex(num_remained_couponleg_cf), _rngchkindex2(num_remained_couponleg_cf);
	vector<string> _index1_dcb(num_remained_couponleg_cf), _index1_stub(num_remained_couponleg_cf), _index1_direction(num_remained_couponleg_cf), _index1_conv(num_remained_couponleg_cf), _index1_adjflag(num_remained_couponleg_cf), _index1_tenors(num_remained_couponleg_cf), _index1fixed_freqs(num_remained_couponleg_cf), _index2_dcb(num_remained_couponleg_cf), _index2_stub(num_remained_couponleg_cf), _index2_direction(num_remained_couponleg_cf), _index2_conv(num_remained_couponleg_cf), _index2_adjflag(num_remained_couponleg_cf), _index2_tenors(num_remained_couponleg_cf), _index2fixed_freqs(num_remained_couponleg_cf),
		_rngchkindex_dcb(num_remained_couponleg_cf), _rngchkindex_stub(num_remained_couponleg_cf), _rngchkindex_direction(num_remained_couponleg_cf), _rngchkindex_conv(num_remained_couponleg_cf), _rngchkindex_adjflag(num_remained_couponleg_cf), _rngchkindex_tenors(num_remained_couponleg_cf), _rngchkindexfixed_freqs(num_remained_couponleg_cf),
		_rngchkindex2_dcb(num_remained_couponleg_cf), _rngchkindex2_stub(num_remained_couponleg_cf), _rngchkindex2_direction(num_remained_couponleg_cf), _rngchkindex2_conv(num_remained_couponleg_cf), _rngchkindex2_adjflag(num_remained_couponleg_cf), _rngchkindex2_tenors(num_remained_couponleg_cf), _rngchkindex2fixed_freqs(num_remained_couponleg_cf);;
	////////////////////////									


	vector<int> index1_spotlag(num_couponleg_cf), index2_spotlag(num_couponleg_cf), rngchkindex_spotlag(num_couponleg_cf), rngchkindex2_spotlag(num_couponleg_cf), num_fixing(num_couponleg_cf), num_lockout(num_couponleg_cf), num_realfixing(num_couponleg_cf), num_fixed(num_couponleg_cf, 0);
	vector<CDate> couponleg_calc_shiftedstartdate(num_couponleg_cf), couponleg_calc_shiftedenddate(num_couponleg_cf), lockout_startdate(num_couponleg_cf);
	vector<vector<CDate>> couponindex_fixingdate(num_couponleg_cf);
	vector<vector<double>> fixedrate1(num_couponleg_cf), fixedrate2(num_couponleg_cf), rngchkfixedrate(num_couponleg_cf), rngchkfixedrate2(num_couponleg_cf);

	//////////////////////////////////////									
	vector<int> _index1_spotlag(num_remained_couponleg_cf), _index2_spotlag(num_remained_couponleg_cf), _rngchkindex_spotlag(num_remained_couponleg_cf), _rngchkindex2_spotlag(num_remained_couponleg_cf), _num_fixing(num_remained_couponleg_cf), _num_lockout(num_remained_couponleg_cf), _num_realfixing(num_remained_couponleg_cf), _num_fixed(num_remained_couponleg_cf, 0);
	vector<CDate> _couponleg_calc_startdate(num_remained_couponleg_cf), _couponleg_calc_shiftedstartdate(num_remained_couponleg_cf), _couponleg_calc_shiftedenddate(num_remained_couponleg_cf), _lockout_startdate(num_remained_couponleg_cf);
	vector<vector<CDate>> _couponindex_fixingdate(num_remained_couponleg_cf);
	vector<vector<double>> _fixedrate1(num_remained_couponleg_cf), _fixedrate2(num_remained_couponleg_cf), _rngchkfixedrate(num_remained_couponleg_cf), _rngchkfixedrate2(num_remained_couponleg_cf);
	//////////////////////////////////////									


	vector<double> couponlegtau(num_couponleg_cf), couponlegT(num_couponleg_cf);

	vector<double> _couponleg_notional(num_remained_couponleg_cf), _couponleg_couponrate(num_remained_couponleg_cf), _couponleg_spread(num_remained_couponleg_cf), _couponlegindex1_mult(num_remained_couponleg_cf), _couponlegindex2_mult(num_remained_couponleg_cf), _rngchkindex_mult(num_remained_couponleg_cf), _rngchkindex2_mult(num_remained_couponleg_cf), _caprates(num_remained_couponleg_cf), _floorrates(num_remained_couponleg_cf), _rahigh_bdry(num_remained_couponleg_cf), _ralow_bdry(num_remained_couponleg_cf), _rahigh_bdryoh(num_remained_couponleg_cf), _ralow_bdryoh(num_remained_couponleg_cf), _couponlegtau(num_remained_couponleg_cf), _couponlegT(num_remained_couponleg_cf);
	vector<bool> _ra_flag(num_remained_couponleg_cf), _rahigh_bdryin_flag(num_remained_couponleg_cf), _ralow_bdryin_flag(num_remained_couponleg_cf);
	vector<CDate> _couponleg_paydate(num_remained_couponleg_cf);

	vector<vector<double>> couponindex_fixingT(num_couponleg_cf);
	vector<vector<double>> _couponindex_fixingT(num_remained_couponleg_cf);

	vector<double> hk(num_couponleg_cf), hs(num_couponleg_cf), lk(num_couponleg_cf), ls(num_couponleg_cf);

	vector<double> _hk(num_remained_couponleg_cf), _hs(num_remained_couponleg_cf), _lk(num_remained_couponleg_cf), _ls(num_remained_couponleg_cf);


	vector<CDate*> index1_holidayss(num_couponleg_cf), index2_holidayss(num_couponleg_cf), rngchkindex_holidayss(num_couponleg_cf), rngchkindex2_holidayss(num_couponleg_cf);
	vector<int> num_index1_cf(num_couponleg_cf), num_index2_cf(num_couponleg_cf), num_rngchkindex_cf(num_couponleg_cf), num_rngchkindex2_cf(num_couponleg_cf), num_index1holidayss(num_couponleg_cf), num_index2holidayss(num_couponleg_cf), num_rngchkindexholidayss(num_couponleg_cf), num_rngchkindex2holidayss(num_couponleg_cf);
	vector<vector<vector<CDate>>> index1_cf(num_couponleg_cf), index2_cf(num_couponleg_cf), rngchkindex_cf(num_couponleg_cf), rngchkindex2_cf(num_couponleg_cf);

	vector<CDate*> _index1_holidayss(num_remained_couponleg_cf), _index2_holidayss(num_remained_couponleg_cf), _rngchkindex_holidayss(num_remained_couponleg_cf), _rngchkindex2_holidayss(num_remained_couponleg_cf);
	vector<int> _num_index1_cf(num_remained_couponleg_cf), _num_index2_cf(num_remained_couponleg_cf), _num_rngchkindex_cf(num_remained_couponleg_cf), _num_rngchkindex2_cf(num_remained_couponleg_cf), _num_index1holidayss(num_remained_couponleg_cf), _num_index2holidayss(num_remained_couponleg_cf), _num_rngchkindexholidayss(num_remained_couponleg_cf), _num_rngchkindex2holidayss(num_remained_couponleg_cf);
	vector<vector<vector<CDate>>> _index1_cf(num_remained_couponleg_cf), _index2_cf(num_remained_couponleg_cf), _rngchkindex_cf(num_remained_couponleg_cf), _rngchkindex2_cf(num_remained_couponleg_cf);

	CDate first_simulate_date = DateAdd('d', 1, couponleg_paydate[num_couponleg_cf - 1]);

	int couponleg_shifted_cf_starti = num_couponleg_cf, _couponleg_shifted_cf_starti = num_remained_couponleg_cf;
	for (i = 0; i < num_couponleg_cf; i++)
	{
		couponlegtau[i] = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
		couponlegT[i] = cvg(today, couponleg_paydate[i], dcb);

		index1[i] = CInstrument(couponlegindex1_crcys[i], couponlegindex1_types[i]);
		index1_dcb[i] = index1[i].get_basis();
		index1_stub[i] = index1[i].get_stub();
		index1_direction[i] = index1[i].get_direction();
		index1_conv[i] = index1[i].get_conv();
		index1_adjflag[i] = index1[i].get_adjflag();
		index1_spotlag[i] = index1[i].get_spotlag();
		index1_holidayss[i] = couponlegindex1_holidayss[i];
		index1_tenors[i] = couponlegindex1_tenors[i];
		index1fixed_freqs[i] = couponlegindex1fixed_freqs[i];
		num_index1holidayss[i] = num_couponlegindex1holidayss[i];

		index2[i] = CInstrument(couponlegindex2_crcys[i], couponlegindex2_types[i]);
		index2_dcb[i] = index2[i].get_basis();
		index2_stub[i] = index2[i].get_stub();
		index2_direction[i] = index2[i].get_direction();
		index2_conv[i] = index2[i].get_conv();
		index2_adjflag[i] = index2[i].get_adjflag();
		index2_spotlag[i] = index2[i].get_spotlag();
		index2_holidayss[i] = couponlegindex2_holidayss[i];
		index2_tenors[i] = couponlegindex2_tenors[i];
		index2fixed_freqs[i] = couponlegindex2fixed_freqs[i];
		num_index2holidayss[i] = num_couponlegindex2holidayss[i];

		rngchkindex[i] = CInstrument(couponlegrngchkindex_crcys[i], couponlegrngchkindex_types[i]);
		rngchkindex_dcb[i] = rngchkindex[i].get_basis();
		rngchkindex_stub[i] = rngchkindex[i].get_stub();
		rngchkindex_direction[i] = rngchkindex[i].get_direction();
		rngchkindex_conv[i] = rngchkindex[i].get_conv();
		rngchkindex_adjflag[i] = rngchkindex[i].get_adjflag();
		rngchkindex_spotlag[i] = rngchkindex[i].get_spotlag();
		rngchkindex_holidayss[i] = couponlegrngchkindex_holidayss[i];
		rngchkindex_tenors[i] = couponlegrngchkindex_tenors[i];
		rngchkindexfixed_freqs[i] = couponlegrngchkindexfixed_freqs[i];
		num_rngchkindexholidayss[i] = num_couponlegrngchkindexholidayss[i];

		rngchkindex2[i] = CInstrument(couponlegrngchkindex2_crcys[i], couponlegrngchkindex2_types[i]);
		rngchkindex2_dcb[i] = rngchkindex2[i].get_basis();
		rngchkindex2_stub[i] = rngchkindex2[i].get_stub();
		rngchkindex2_direction[i] = rngchkindex2[i].get_direction();
		rngchkindex2_conv[i] = rngchkindex2[i].get_conv();
		rngchkindex2_adjflag[i] = rngchkindex2[i].get_adjflag();
		rngchkindex2_spotlag[i] = rngchkindex2[i].get_spotlag();
		rngchkindex2_holidayss[i] = couponlegrngchkindex2_holidayss[i];
		rngchkindex2_tenors[i] = couponlegrngchkindex2_tenors[i];
		rngchkindex2fixed_freqs[i] = couponlegrngchkindex2fixed_freqs[i];
		num_rngchkindex2holidayss[i] = num_couponlegrngchkindex2holidayss[i];

		num_index1_cf[i] = findnumschedule(today, index1_tenors[i], index1_stub[i], index1_direction[i], index1fixed_freqs[i]);
		num_index2_cf[i] = findnumschedule(today, index2_tenors[i], index2_stub[i], index2_direction[i], index2fixed_freqs[i]);
		num_rngchkindex_cf[i] = findnumschedule(today, rngchkindex_tenors[i], rngchkindex_stub[i], rngchkindex_direction[i], rngchkindexfixed_freqs[i]);
		num_rngchkindex2_cf[i] = findnumschedule(today, rngchkindex2_tenors[i], rngchkindex2_stub[i], rngchkindex2_direction[i], rngchkindex2fixed_freqs[i]);


		index1_cf[i] = vector<vector<CDate>>(3);
		index2_cf[i] = vector<vector<CDate>>(3);
		rngchkindex_cf[i] = vector<vector<CDate>>(3);
		rngchkindex2_cf[i] = vector<vector<CDate>>(3);

		for (j = 0; j < 3; j++)
		{
			index1_cf[i][j] = vector<CDate>(num_index1_cf[i]);
			index2_cf[i][j] = vector<CDate>(num_index2_cf[i]);
			rngchkindex_cf[i][j] = vector<CDate>(num_rngchkindex_cf[i]);
			rngchkindex2_cf[i][j] = vector<CDate>(num_rngchkindex2_cf[i]);
		}

		if (ra_flag[i])
		{
			couponleg_shifted_cf_starti = min(couponleg_shifted_cf_starti, i);
			couponleg_calc_shiftedstartdate[i] = couponleg_calc_startdate[i];
			couponleg_calc_shiftedenddate[i] = couponleg_calc_enddate[i];
			ShiftBusDate(couponleg_calc_shiftedstartdate[i], couponleg_holidays, num_couponlegholidays, nshift, couponleg_calc_shiftedstartdate[i]);
			ShiftBusDate(couponleg_calc_shiftedenddate[i], couponleg_holidays, num_couponlegholidays, nshift, couponleg_calc_shiftedenddate[i]);
			lockout_startdate[i] = couponleg_calc_shiftedenddate[i];
			ShiftBusDate(lockout_startdate[i], couponleg_holidays, num_couponlegholidays, nlockout, lockout_startdate[i]);
			num_fixing[i] = int(couponleg_calc_shiftedenddate[i] - couponleg_calc_shiftedstartdate[i]) + (ralastdatein_flag - 1) + rafirstdatein_flag;
			num_lockout[i] = int(couponleg_calc_shiftedenddate[i] - lockout_startdate[i]) + (ralastdatein_flag - 1);
			num_realfixing[i] = num_fixing[i] - num_lockout[i];
			couponindex_fixingdate[i] = vector<CDate>(num_fixing[i]);
			if (num_couponindexfixing_history > 0 && couponindexfixinghistory_starti < num_couponindexfixing_history)
			{
				while (couponindexfixinghistory_date[couponindexfixinghistory_starti] < couponleg_calc_shiftedstartdate[i])
				{
					couponindexfixinghistory_starti = couponindexfixinghistory_starti + 1;
					if (couponindexfixinghistory_starti >= num_couponindexfixing_history) break;
				}

				j = couponindexfixinghistory_starti;

				while (j < min(num_couponindexfixing_history, couponindexfixinghistory_starti + num_fixing[i]))
				{
					couponindex_fixingdate[i][j - couponindexfixinghistory_starti] = couponindexfixinghistory_date[j];
					if (is_holiday(couponindexfixinghistory_date[j], couponlegindex1_holidayss[i], num_couponlegindex1holidayss[i]))
					{
						if (int(fixedrate1[i].size()) > 0) fixedrate1[i].push_back(fixedrate1[i][int(fixedrate1[i].size()) - 1]);
					}
					else fixedrate1[i].push_back(couponindex1fixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegindex2_holidayss[i], num_couponlegindex2holidayss[i]))
					{
						if (int(fixedrate2[i].size()) > 0) fixedrate2[i].push_back(fixedrate2[i][int(fixedrate2[i].size()) - 1]);
					}
					else fixedrate2[i].push_back(couponindex2fixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegrngchkindex_holidayss[i], num_couponlegrngchkindexholidayss[i]))
					{
						if (int(rngchkfixedrate[i].size()) > 0) rngchkfixedrate[i].push_back(rngchkfixedrate[i][int(rngchkfixedrate[i].size()) - 1]);
					}
					else rngchkfixedrate[i].push_back(couponrngchkindexfixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegrngchkindex2_holidayss[i], num_couponlegrngchkindex2holidayss[i]))
					{
						if (int(rngchkfixedrate2[i].size()) > 0) rngchkfixedrate2[i].push_back(rngchkfixedrate2[i][int(rngchkfixedrate2[i].size()) - 1]);
					}
					else rngchkfixedrate2[i].push_back(couponrngchkindex2fixinghistory_rate[j]);

					j = j + 1;
				}
				num_fixed[i] = int(fixedrate1[i].size());
				for (j = num_realfixing[i]; j < num_fixed[i]; j++)
				{
					fixedrate1[i][j] = fixedrate1[i][j - 1];
					fixedrate2[i][j] = fixedrate2[i][j - 1];
					rngchkfixedrate[i][j] = rngchkfixedrate[i][j - 1];
					rngchkfixedrate2[i][j] = rngchkfixedrate2[i][j - 1];
				}
				if (num_fixed[i] > num_realfixing[i])
				{
					for (j = num_fixed[i]; j < num_fixing[i]; j++)
					{
						couponindex_fixingdate[i][j] = DateAdd('d', 1, couponindex_fixingdate[i][j - 1]);
						fixedrate1[i].push_back(fixedrate1[i][j - 1]);
						fixedrate2[i].push_back(fixedrate2[i][j - 1]);
						rngchkfixedrate[i].push_back(rngchkfixedrate[i][j - 1]);
						rngchkfixedrate2[i].push_back(rngchkfixedrate2[i][j - 1]);
					}
					num_fixed[i] = int(fixedrate1[i].size());
				}
				if (num_fixed[i] > 0) { if (fixedrate1[i][num_fixed[i] - 1] <= MinFixingRate) num_fixed[i] = num_fixed[i] - 1; }
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////									
	for (i = 0; i < num_remained_couponleg_cf; i++)
	{
		_couponleg_notional[i] = couponleg_notional[i + couponleg_cf_starti];
		_couponleg_couponrate[i] = couponleg_couponrate[i + couponleg_cf_starti];
		_couponleg_spread[i] = couponleg_spread[i + couponleg_cf_starti];
		_couponleg_calc_startdate[i] = couponleg_calc_startdate[i + couponleg_cf_starti];
		_couponleg_paydate[i] = couponleg_paydate[i + couponleg_cf_starti];
		_couponlegindex1_mult[i] = couponlegindex1_mult[i + couponleg_cf_starti];
		_couponlegindex2_mult[i] = couponlegindex2_mult[i + couponleg_cf_starti];

		_caprates[i] = caprates[i + couponleg_cf_starti];
		_floorrates[i] = floorrates[i + couponleg_cf_starti];
		_ra_flag[i] = ra_flag[i + couponleg_cf_starti];
		_rahigh_bdry[i] = rahigh_bdry[i + couponleg_cf_starti];
		_ralow_bdry[i] = ralow_bdry[i + couponleg_cf_starti];
		_rahigh_bdryoh[i] = rahigh_bdryoh[i + couponleg_cf_starti];
		_ralow_bdryoh[i] = ralow_bdryoh[i + couponleg_cf_starti];
		_rahigh_bdryin_flag[i] = rahigh_bdryin_flag[i + couponleg_cf_starti];
		_ralow_bdryin_flag[i] = ralow_bdryin_flag[i + couponleg_cf_starti];
		_couponlegtau[i] = couponlegtau[i + couponleg_cf_starti];
		_couponlegT[i] = couponlegT[i + couponleg_cf_starti];

		_index1[i] = index1[i + couponleg_cf_starti];
		_index1_dcb[i] = index1_dcb[i + couponleg_cf_starti];
		_index1_stub[i] = index1_stub[i + couponleg_cf_starti];
		_index1_direction[i] = index1_direction[i + couponleg_cf_starti];
		_index1_conv[i] = index1_conv[i + couponleg_cf_starti];
		_index1_adjflag[i] = index1_adjflag[i + couponleg_cf_starti];
		_index1_spotlag[i] = index1_spotlag[i + couponleg_cf_starti];
		_index1_holidayss[i] = index1_holidayss[i + couponleg_cf_starti];
		_index1_tenors[i] = index1_tenors[i + couponleg_cf_starti];
		_index1fixed_freqs[i] = index1fixed_freqs[i + couponleg_cf_starti];
		_num_index1holidayss[i] = num_index1holidayss[i + couponleg_cf_starti];

		_index2[i] = index2[i + couponleg_cf_starti];
		_index2_dcb[i] = index2_dcb[i + couponleg_cf_starti];
		_index2_stub[i] = index2_stub[i + couponleg_cf_starti];
		_index2_direction[i] = index2_direction[i + couponleg_cf_starti];
		_index2_conv[i] = index2_conv[i + couponleg_cf_starti];
		_index2_adjflag[i] = index2_adjflag[i + couponleg_cf_starti];
		_index2_spotlag[i] = index2_spotlag[i + couponleg_cf_starti];
		_index2_holidayss[i] = index2_holidayss[i + couponleg_cf_starti];
		_index2_tenors[i] = index2_tenors[i + couponleg_cf_starti];
		_index2fixed_freqs[i] = index2fixed_freqs[i + couponleg_cf_starti];
		_num_index2holidayss[i] = num_index2holidayss[i + couponleg_cf_starti];

		_rngchkindex[i] = rngchkindex[i + couponleg_cf_starti];
		_rngchkindex_dcb[i] = rngchkindex_dcb[i + couponleg_cf_starti];
		_rngchkindex_stub[i] = rngchkindex_stub[i + couponleg_cf_starti];
		_rngchkindex_direction[i] = rngchkindex_direction[i + couponleg_cf_starti];
		_rngchkindex_conv[i] = rngchkindex_conv[i + couponleg_cf_starti];
		_rngchkindex_adjflag[i] = rngchkindex_adjflag[i + couponleg_cf_starti];
		_rngchkindex_spotlag[i] = rngchkindex_spotlag[i + couponleg_cf_starti];
		_rngchkindex_holidayss[i] = rngchkindex_holidayss[i + couponleg_cf_starti];
		_rngchkindex_tenors[i] = rngchkindex_tenors[i + couponleg_cf_starti];
		_rngchkindexfixed_freqs[i] = rngchkindexfixed_freqs[i + couponleg_cf_starti];
		_num_rngchkindexholidayss[i] = num_rngchkindexholidayss[i + couponleg_cf_starti];
		_rngchkindex_mult[i] = couponlegrngchkindex_mult[i + couponleg_cf_starti];

		_rngchkindex2[i] = rngchkindex2[i + couponleg_cf_starti];
		_rngchkindex2_dcb[i] = rngchkindex2_dcb[i + couponleg_cf_starti];
		_rngchkindex2_stub[i] = rngchkindex2_stub[i + couponleg_cf_starti];
		_rngchkindex2_direction[i] = rngchkindex2_direction[i + couponleg_cf_starti];
		_rngchkindex2_conv[i] = rngchkindex2_conv[i + couponleg_cf_starti];
		_rngchkindex2_adjflag[i] = rngchkindex2_adjflag[i + couponleg_cf_starti];
		_rngchkindex2_spotlag[i] = rngchkindex2_spotlag[i + couponleg_cf_starti];
		_rngchkindex2_holidayss[i] = rngchkindex2_holidayss[i + couponleg_cf_starti];
		_rngchkindex2_tenors[i] = rngchkindex2_tenors[i + couponleg_cf_starti];
		_rngchkindex2fixed_freqs[i] = rngchkindex2fixed_freqs[i + couponleg_cf_starti];
		_num_rngchkindex2holidayss[i] = num_rngchkindex2holidayss[i + couponleg_cf_starti];
		_rngchkindex2_mult[i] = couponlegrngchkindex2_mult[i + couponleg_cf_starti];

		_num_index1_cf[i] = num_index1_cf[i + couponleg_cf_starti];
		_num_index2_cf[i] = num_index2_cf[i + couponleg_cf_starti];
		_num_rngchkindex_cf[i] = num_rngchkindex_cf[i + couponleg_cf_starti];
		_num_rngchkindex2_cf[i] = num_rngchkindex2_cf[i + couponleg_cf_starti];

		_index1_cf[i] = index1_cf[i + couponleg_cf_starti];
		_index2_cf[i] = index2_cf[i + couponleg_cf_starti];
		_rngchkindex_cf[i] = rngchkindex_cf[i + couponleg_cf_starti];
		_rngchkindex2_cf[i] = rngchkindex2_cf[i + couponleg_cf_starti];

		if (_ra_flag[i])
		{
			_couponleg_shifted_cf_starti = min(_couponleg_shifted_cf_starti, i);
			_couponleg_calc_shiftedstartdate[i] = couponleg_calc_shiftedstartdate[i + couponleg_cf_starti];
			_couponleg_calc_shiftedenddate[i] = couponleg_calc_shiftedenddate[i + couponleg_cf_starti];
			_lockout_startdate[i] = lockout_startdate[i + couponleg_cf_starti];
			_num_fixing[i] = num_fixing[i + couponleg_cf_starti];
			_num_lockout[i] = num_lockout[i + couponleg_cf_starti];
			_num_realfixing[i] = num_realfixing[i + couponleg_cf_starti];
			_couponindex_fixingdate[i] = vector<CDate>(_num_fixing[i]);
			if (num_couponindexfixing_history > 0 && _couponindexfixinghistory_starti < num_couponindexfixing_history)
			{
				while (couponindexfixinghistory_date[_couponindexfixinghistory_starti] < _couponleg_calc_shiftedstartdate[i])
				{
					_couponindexfixinghistory_starti = _couponindexfixinghistory_starti + 1;
					if (_couponindexfixinghistory_starti >= num_couponindexfixing_history) break;
				}

				j = _couponindexfixinghistory_starti;

				while (j < min(num_couponindexfixing_history, _couponindexfixinghistory_starti + _num_fixing[i]))
				{
					_couponindex_fixingdate[i][j - _couponindexfixinghistory_starti] = couponindexfixinghistory_date[j];

					if (is_holiday(couponindexfixinghistory_date[j], couponlegindex1_holidayss[i], num_couponlegindex1holidayss[i]))
					{
						if (int(_fixedrate1[i].size()) > 0) _fixedrate1[i].push_back(_fixedrate1[i][int(_fixedrate1[i].size()) - 1]);
					}
					else _fixedrate1[i].push_back(couponindex1fixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegindex2_holidayss[i], num_couponlegindex2holidayss[i]))
					{
						if (int(_fixedrate2[i].size()) > 0) _fixedrate2[i].push_back(_fixedrate2[i][int(_fixedrate2[i].size()) - 1]);
					}
					else _fixedrate2[i].push_back(couponindex2fixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegrngchkindex_holidayss[i], num_couponlegrngchkindexholidayss[i]))
					{
						if (int(_rngchkfixedrate[i].size()) > 0) _rngchkfixedrate[i].push_back(_rngchkfixedrate[i][int(_rngchkfixedrate[i].size()) - 1]);
					}
					else _rngchkfixedrate[i].push_back(couponrngchkindexfixinghistory_rate[j]);

					if (is_holiday(couponindexfixinghistory_date[j], couponlegrngchkindex2_holidayss[i], num_couponlegrngchkindex2holidayss[i]))
					{
						if (int(_rngchkfixedrate2[i].size()) > 0) _rngchkfixedrate2[i].push_back(_rngchkfixedrate2[i][int(_rngchkfixedrate2[i].size()) - 1]);
					}
					else _rngchkfixedrate2[i].push_back(couponrngchkindex2fixinghistory_rate[j]);

					j = j + 1;
				}
				_num_fixed[i] = int(_fixedrate1[i].size());
				for (j = _num_realfixing[i]; j < _num_fixed[i]; j++)
				{
					_fixedrate1[i][j] = _fixedrate1[i][j - 1];
					_fixedrate2[i][j] = _fixedrate2[i][j - 1];
					_rngchkfixedrate[i][j] = _rngchkfixedrate[i][j - 1];
					_rngchkfixedrate2[i][j] = _rngchkfixedrate2[i][j - 1];
				}
				if (_num_fixed[i] > _num_realfixing[i])
				{
					for (j = _num_fixed[i]; j < _num_fixing[i]; j++)
					{
						_couponindex_fixingdate[i][j] = DateAdd('d', 1, _couponindex_fixingdate[i][j - 1]);
						_fixedrate1[i].push_back(fixedrate1[i][j - 1]);
						_fixedrate2[i].push_back(fixedrate2[i][j - 1]);
						_rngchkfixedrate[i].push_back(rngchkfixedrate[i][j - 1]);
						_rngchkfixedrate2[i].push_back(rngchkfixedrate2[i][j - 1]);
					}
					_num_fixed[i] = int(_fixedrate1[i].size());
				}
				if (_num_fixed[i] > 0) { if (_fixedrate1[i][_num_fixed[i] - 1] <= MinFixingRate) _num_fixed[i] = _num_fixed[i] - 1; }
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////									
	double
		floatinglegfixedrpice,
		floatinglegnonfixedprice = 0.0,
		floatinglegprice = 0.0,
		fixedlegfixedprice = 0.0,
		fixedlegaccrualfixedprice = 0.0,
		fixedlegnonfixedprice = 0.0;

	PlainSwapFundinglegFixedPrice
	(
		today,
		ondf,
		settlement_date_df,
		zc,
		fixinghistory_rate,
		num_fundingleg_cf,
		fundingleg_cf_starti, fixing_date,
		fundingleg_calc_startdate,
		fundingleg_calc_enddate,
		fundingleg_paydate,
		fundingleg_dcb,
		fundingleg_notional,
		fundingleg_mult,
		fundingleg_spread,
		floatinglegfixedrpice
	);

	PlainSwapFundinglegNonFixedPrice
	(
		today,
		ondf,
		settlement_date_df,
		zc,
		num_remained_fundingleg_cf,
		_fundingleg_notional,
		_fundingleg_mult,
		_fundingleg_spread,
		fundingleg_calc_startdate0,
		fundinglegcalcstartT0,
		_fundingleg_paydate,
		fundinglegtau,
		fundinglegT,
		floatinglegnonfixedprice
	);

	PlainSwapCouponlegFixedPrice
	(
		today,
		ondf,
		settlement_date_df,
		zc,
		num_remained_couponleg_cf,
		_ra_flag,
		_couponleg_notional,
		_couponleg_couponrate,
		_couponleg_paydate,
		couponlegtau,
		couponlegT,
		fixedlegfixedprice
	);

	vector<double>
		ratecount(num_couponleg_cf, 0.0),
		fixedavg(num_couponleg_cf, 0.0),
		_ratecount(num_remained_couponleg_cf, 0.0),
		_fixedavg(num_remained_couponleg_cf, 0.0);

	//QuantoLeveragedAverageFixedAccrualPrice2(today,ondf,settlement_date_df,zc,num_remained_couponleg_cf,_num_fixing,_num_fixed,_fixedrate1,_fixedrate2,_rngchkfixedrate, _rngchkfixedrate2, _rngchkindex_mult, _rngchkindex2_mult, _ra_flag,_couponleg_notional,_couponleg_spread,_couponlegindex1_mult,_couponlegindex2_mult,_caprates,_floorrates,_rahigh_bdry,_rahigh_bdryin_flag,_ralow_bdry,_ralow_bdryin_flag,_couponleg_paydate,_couponlegtau,_couponlegT,_ratecount,_fixedavg);								

	int num_past_couponleg_cf = num_couponleg_cf - num_remained_couponleg_cf;

	QuantoAvgSpreadLiborRAFixedPrice2
	(
		today,
		ondf,
		settlement_date_df,
		zc,
		num_past_couponleg_cf,
		num_fixing,
		num_fixed,
		couponleg_couponrate,
		fixedrate1,
		fixedrate2,
		rngchkfixedrate,
		rngchkfixedrate2,
		couponlegrngchkindex_mult,
		couponlegrngchkindex2_mult,
		ra_flag,
		couponleg_notional,
		couponleg_spread,
		couponlegindex1_mult,
		couponlegindex2_mult,
		caprates,
		floorrates,
		rahigh_bdry,
		rahigh_bdryin_flag,
		ralow_bdry,
		ralow_bdryin_flag,
		couponleg_paydate,
		couponlegtau,
		couponlegT,
		ratecount,
		fixedavg,
		fixedlegaccrualfixedprice
	);

	if (num_calldate > 0)
	{
		if (callstarti < num_calldate)
		{
			if (first_simulate_date > callnotice_date[callstarti])
			{
				first_simulate_date = callnotice_date[callstarti];
			}
		}
	}

	int ra_starti = 0;

	if (num_remained_couponleg_cf > 0)
	{
		while (!_ra_flag[ra_starti])
		{
			ra_starti = ra_starti + 1;

			if (ra_starti >= num_remained_couponleg_cf - 1)
			{
				break;
			}
		}
	}
	if (_ra_flag[ra_starti])
	{
		if (today >= _couponleg_calc_shiftedstartdate[ra_starti])
		{
			first_simulate_date = today;
		}
	}

	int
		num_rates = max(int(couponleg_paydate[num_couponleg_cf - 1] - first_simulate_date + 1), 0),
		accrualstarti = num_rates;

	vector<double>
		dt(num_rates, Daily),
		t(num_rates),
		r1(num_rates),
		r2(num_rates),
		fM1(num_rates),
		fM2(num_rates),
		alpha1(num_rates),
		alpha2(num_rates),
		t_firststartT(num_rates),
		t_lastpayT(num_rates),
		index1t_firststartT(num_rates),
		index2t_firststartT(num_rates),
		rngchkindext_firststartT(num_rates),
		rngchkindex2t_firststartT(num_rates),
		ema1dt(num_rates), ema2dt(num_rates),
		emb2dt(num_rates), qadjx2(num_rates),
		qadjy2(num_rates),
		PMt_firststartT(num_rates),
		PMt_lastpayT(num_rates),
		PMindex1t_firststartT(num_rates),
		PMindex2t_firststartT(num_rates),
		PMrngchkindext_firststartT(num_rates),
		PMrngchkindex2t_firststartT(num_rates),
		Bt_firststartT(num_rates),
		Bt_lastpayT(num_rates),
		Bindex1a2t_firststartT(num_rates),
		Bindex1b2t_firststartT(num_rates),
		Bindex2a2t_firststartT(num_rates),
		Bindex2b2t_firststartT(num_rates),
		Brngchkindexa1t_firststartT(num_rates),
		Brngchkindex2a2t_firststartT(num_rates),
		Brngchkindex2b2t_firststartT(num_rates),
		Vt(num_rates),
		Vindex1t(num_rates),
		Vindex1_firststartT(num_rates),
		Vindex1t_firststartT(num_rates),
		Vindex2t(num_rates),
		Vindex2_firststartT(num_rates),
		Vindex2t_firststartT(num_rates),
		Vrngchkindext(num_rates),
		Vrngchkindex2t(num_rates),
		Vrngchkindex_firststartT(num_rates),
		Vrngchkindex2_firststartT(num_rates),
		Vrngchkindext_firststartT(num_rates),
		Vrngchkindex2t_firststartT(num_rates),
		coeft_firststartT(num_rates),
		coeft_lastpayT(num_rates),
		coefindex1t_firststartT(num_rates),
		coefindex2t_firststartT(num_rates),
		coefrngchkindext_firststartT(num_rates),
		coefrngchkindex2t_firststartT(num_rates);

	vector<CDate>
		date(num_rates),
		fund_firststartdate(num_rates),
		fund_lastpaydate(num_rates),
		index1_firststartdate(num_rates),
		index1_lastpaydate(num_rates),
		index2_firststartdate(num_rates),
		index2_lastpaydate(num_rates),
		rngchkindex_firststartdate(num_rates),
		rngchkindex_lastpaydate(num_rates),
		rngchkindex2_firststartdate(num_rates),
		rngchkindex2_lastpaydate(num_rates);

	vector<vector<CDate>>
		index1_paydate(num_rates),
		index2_paydate(num_rates),
		rngchkindex_paydate(num_rates),
		rngchkindex2_paydate(num_rates);

	vector<double>
		fund_firststartT(num_rates),
		fund_lastpayT(num_rates),
		fund_tau(num_rates),
		index1_firststartT(num_rates),
		index2_firststartT(num_rates),
		rngchkindex_firststartT(num_rates),
		rngchkindex2_firststartT(num_rates);

	vector<vector<double>>
		index1_payT(num_rates),
		index2_payT(num_rates),
		rngchkindex_payT(num_rates),
		rngchkindex2_payT(num_rates),
		index1_tau(num_rates),
		index2_tau(num_rates),
		rngchkindex_tau(num_rates),
		rngchkindex2_tau(num_rates),
		index1t_payT(num_rates),
		index2t_payT(num_rates),
		rngchkindext_payT(num_rates),
		rngchkindex2t_payT(num_rates),
		PMindex1t_payT(num_rates),
		PMindex2t_payT(num_rates),
		PMrngchkindext_payT(num_rates),
		PMrngchkindex2t_payT(num_rates),
		Bindex1a2t_payT(num_rates),
		Bindex1b2t_payT(num_rates),
		Bindex2a2t_payT(num_rates),
		Bindex2b2t_payT(num_rates),
		Brngchkindexa1t_payT(num_rates),
		Brngchkindex2a2t_payT(num_rates),
		Brngchkindex2b2t_payT(num_rates),
		Vindex1_payT(num_rates),
		Vindex1t_payT(num_rates),
		Vindex2_payT(num_rates),
		Vindex2t_payT(num_rates),
		Vrngchkindex_payT(num_rates),
		Vrngchkindext_payT(num_rates),
		Vrngchkindex2_payT(num_rates),
		Vrngchkindex2t_payT(num_rates),
		coefindex1t_payT(num_rates),
		coefindex2t_payT(num_rates),
		coefrngchkindext_payT(num_rates),
		coefrngchkindex2t_payT(num_rates);

	int num_factor = 3;			//check					

	vector<double>
		stdevx1(num_rates),
		stdevx2(num_rates),
		stdevy2(num_rates),
		corrx1x2(num_rates),
		corrx1y2(num_rates), corrx2y2(num_rates);

	vector<vector<vector<double>>>
		l_mtrx(num_rates),
		corr(num_rates);

	double deltat = Daily;

	int param_indx = 0;

	while (first_simulate_date > paramterm_date[param_indx])
	{
		param_indx = param_indx + 1;
		if (param_indx >= num_paramterm - 1) break;
	}

	int coupleg_indx = 0;

	for (i = 0; i < num_remained_couponleg_cf; i++)
	{
		if (_ra_flag[i])
		{
			if (first_simulate_date <= _couponleg_calc_shiftedenddate[i])
			{
				coupleg_indx = i;
				break;
			}
		}
		else
		{
			if (first_simulate_date <= couponleg_calc_enddate[i])
			{
				coupleg_indx = i;
				break;
			}
		}
	}

	date[0] = first_simulate_date;

	vector<int>
		i_call,
		i_fixing,
		i_fund_calc_start,
		i_fund_pay,
		i_coup_calc_shifted_start,
		i_coup_calc_shifted_end,
		i_coup_pay;

	int
		count_call = callstarti,
		count_fixing = 0,
		count_fund_calc_start = 0,
		count_fund_pay = 0,
		count_coup_calc_shifted_start = couponleg_shifted_cf_starti,
		count_coup_calc_shifted_end = couponleg_shifted_cf_starti,
		count_coup_pay = couponleg_shifted_cf_starti;

	vector<double> fund_pay_df, coup_pay_df;

	if (callstarti < num_calldate)
	{
		if (date[0] == callnotice_date[count_call])
		{
			i_call.push_back(0);
			if (count_call < num_calldate - 1) count_call = count_call + 1;
		}
	}

	if (num_remained_fundingleg_cf > 0)
	{
		if (date[0] > _fundingleg_calc_startdate[count_fund_calc_start] && count_fund_calc_start < num_remained_fundingleg_cf - 1)
		{
			while (date[0] > _fundingleg_calc_startdate[count_fund_calc_start])
			{
				count_fund_calc_start = count_fund_calc_start + 1;
				if (count_fund_calc_start >= num_remained_fundingleg_cf - 1) break;
			}
		}
		if (date[0] == _fundingleg_calc_startdate[count_fund_calc_start])
		{
			i_fund_calc_start.push_back(0);
			if (count_fund_calc_start < num_remained_fundingleg_cf - 1) count_fund_calc_start = count_fund_calc_start + 1;
		}
		if (date[0] > _fixing_date[count_fixing] && count_fixing < num_remained_fundingleg_cf - 1)
		{
			while (date[0] > _fixing_date[count_fixing])
			{
				count_fixing = count_fixing + 1;
				if (count_fixing >= num_remained_fundingleg_cf - 1) break;
			}
		}
		if (date[0] == _fixing_date[count_fixing])
		{
			i_fixing.push_back(0);
			if (count_fixing < num_remained_fundingleg_cf - 1) count_fixing = count_fixing + 1;
		}
		while (date[0] > _fundingleg_paydate[count_fund_pay])
		{
			count_fund_pay = count_fund_pay + 1;
			if (count_fund_pay >= num_remained_fundingleg_cf - 1) break;
		}
		if (date[0] == _fundingleg_paydate[count_fund_pay])
		{
			fund_pay_df.push_back(exp(-zc(_fundingleg_paydate[count_fund_pay]) * fundinglegT[count_fund_pay]));
			i_fund_pay.push_back(0);
			if (count_fund_pay < num_remained_fundingleg_cf - 1) count_fund_pay = count_fund_pay + 1;
		}
	}
	if (date[0] > _couponleg_calc_shiftedstartdate[count_coup_calc_shifted_start] && count_coup_calc_shifted_start < num_remained_couponleg_cf - 1)
	{
		while (date[0] > _couponleg_calc_shiftedstartdate[count_coup_calc_shifted_start])
		{
			count_coup_calc_shifted_start = count_coup_calc_shifted_start + 1;
			if (count_coup_calc_shifted_start >= num_remained_couponleg_cf - 1) break;
		}
	}
	if (date[0] == _couponleg_calc_shiftedstartdate[count_coup_calc_shifted_start])
	{
		i_coup_calc_shifted_start.push_back(0);
		if (count_coup_calc_shifted_start < num_remained_couponleg_cf - 1) count_coup_calc_shifted_start = count_coup_calc_shifted_start + 1;
	}
	if (date[0] > _couponleg_calc_shiftedenddate[count_coup_calc_shifted_end] && count_coup_calc_shifted_end < num_remained_couponleg_cf - 1)
	{
		while (date[0] > _couponleg_calc_shiftedenddate[count_coup_calc_shifted_end])
		{
			count_coup_calc_shifted_end = count_coup_calc_shifted_end + 1;
			if (count_coup_calc_shifted_end >= num_remained_couponleg_cf - 1) break;
		}
	}
	if (date[0] == _couponleg_calc_shiftedenddate[count_coup_calc_shifted_end])
	{
		i_coup_calc_shifted_end.push_back(0);
		if (count_coup_calc_shifted_end < num_remained_couponleg_cf - 1) count_coup_calc_shifted_end = count_coup_calc_shifted_end + 1;
	}
	while (date[0] > _couponleg_paydate[count_coup_pay])
	{
		count_coup_pay = count_coup_pay + 1;
		if (count_coup_pay >= num_remained_couponleg_cf - 1) break;
	}
	if (date[0] == _couponleg_paydate[count_coup_pay])
	{
		coup_pay_df.push_back(exp(-zc(_couponleg_paydate[count_coup_pay]) * couponlegT[count_coup_pay]));
		i_coup_pay.push_back(0);
		if (count_coup_pay < num_remained_couponleg_cf - 1) count_coup_pay = count_coup_pay + 1;
	}
	if (date[0] >= _couponleg_calc_shiftedstartdate[couponleg_shifted_cf_starti]) accrualstarti = min(0, accrualstarti);
	ShiftBusDate(date[0], fixing_holidays, num_fixingholidays, -fix_lag, fund_firststartdate[0]);
	findmaturity(fund_firststartdate[0], fundinglegindx_tenor, fixing_holidays, num_fixingholidays, fixingindex_conv, fund_lastpaydate[0]);
	t[0] = cvg(today, first_simulate_date, dcb);
	fund_firststartT[0] = cvg(today, fund_firststartdate[0], dcb);
	fund_lastpayT[0] = cvg(today, fund_lastpaydate[0], dcb);
	fund_tau[0] = cvg(fund_firststartdate[0], fund_lastpaydate[0], fundingleg_dcb);
	t_firststartT[0] = fund_firststartT[0] - t[0];
	t_lastpayT[0] = fund_lastpayT[0] - t[0];
	PMt_firststartT[0] = exp(zct(t[0]) * t[0] - zct(fund_firststartT[0]) * fund_firststartT[0]);
	PMt_lastpayT[0] = exp(zct(t[0]) * t[0] - zct(fund_lastpayT[0]) * fund_lastpayT[0]);

	ShiftBusDate(date[0], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_spotlag[coupleg_indx], index1_firststartdate[0]);
	findmaturity(index1_firststartdate[0], index1_tenors[coupleg_indx], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_conv[coupleg_indx], index1_lastpaydate[0]);
	index1_firststartT[0] = cvg(today, index1_firststartdate[0], dcb);
	index1t_firststartT[0] = index1_firststartT[0] - t[0];
	PMindex1t_firststartT[0] = exp(index1zct(t[0]) * t[0] - index1zct(index1_firststartT[0]) * index1_firststartT[0]);
	fixedlegcashflowschedule(index1_firststartdate[0], index1_tenors[coupleg_indx], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_stub[coupleg_indx], index1_direction[coupleg_indx], index1_conv[coupleg_indx], index1fixed_freqs[coupleg_indx], index1_adjflag[coupleg_indx], index1_cf[coupleg_indx]);
	index1_paydate[0] = index1_cf[coupleg_indx][0];
	for (i = 0; i < num_index1_cf[coupleg_indx]; i++)
	{
		index1_payT[0].push_back(cvg(today, index1_paydate[0][i], dcb));
		index1_tau[0].push_back(cvg(index1_cf[coupleg_indx][1][i], index1_cf[coupleg_indx][2][i], index1_dcb[coupleg_indx]));
		index1t_payT[0].push_back(index1_payT[0][i] - t[0]);
		PMindex1t_payT[0].push_back(exp(index1zct(t[0]) * t[0] - index1zct(index1_payT[0][i]) * index1_payT[0][i]));
	}

	ShiftBusDate(date[0], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_spotlag[coupleg_indx], index2_firststartdate[0]);
	findmaturity(index2_firststartdate[0], index2_tenors[coupleg_indx], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_conv[coupleg_indx], index2_lastpaydate[0]);
	index2_firststartT[0] = cvg(today, index2_firststartdate[0], dcb);
	index2t_firststartT[0] = index2_firststartT[0] - t[0];
	PMindex2t_firststartT[0] = exp(index2zct(t[0]) * t[0] - index2zct(index2_firststartT[0]) * index2_firststartT[0]);
	fixedlegcashflowschedule(index2_firststartdate[0], index2_tenors[coupleg_indx], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_stub[coupleg_indx], index2_direction[coupleg_indx], index2_conv[coupleg_indx], index2fixed_freqs[coupleg_indx], index2_adjflag[coupleg_indx], index2_cf[coupleg_indx]);
	index2_paydate[0] = index2_cf[coupleg_indx][0];
	for (i = 0; i < num_index2_cf[coupleg_indx]; i++)
	{
		index2_payT[0].push_back(cvg(today, index2_paydate[0][i], dcb));
		index2_tau[0].push_back(cvg(index2_cf[coupleg_indx][1][i], index2_cf[coupleg_indx][2][i], index2_dcb[coupleg_indx]));
		index2t_payT[0].push_back(index2_payT[0][i] - t[0]);
		PMindex2t_payT[0].push_back(exp(index2zct(t[0]) * t[0] - index2zct(index2_payT[0][i]) * index2_payT[0][i]));
	}

	ShiftBusDate(date[0], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_spotlag[coupleg_indx], rngchkindex_firststartdate[0]);
	findmaturity(rngchkindex_firststartdate[0], rngchkindex_tenors[coupleg_indx], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_conv[coupleg_indx], rngchkindex_lastpaydate[0]);
	rngchkindex_firststartT[0] = cvg(today, rngchkindex_firststartdate[0], dcb);
	rngchkindext_firststartT[0] = rngchkindex_firststartT[0] - t[0];
	PMrngchkindext_firststartT[0] = exp(rngchkindexzct(t[0]) * t[0] - rngchkindexzct(rngchkindex_firststartT[0]) * rngchkindex_firststartT[0]);
	fixedlegcashflowschedule(rngchkindex_firststartdate[0], rngchkindex_tenors[coupleg_indx], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_stub[coupleg_indx], rngchkindex_direction[coupleg_indx], rngchkindex_conv[coupleg_indx], rngchkindexfixed_freqs[coupleg_indx], rngchkindex_adjflag[coupleg_indx], rngchkindex_cf[coupleg_indx]);
	rngchkindex_paydate[0] = rngchkindex_cf[coupleg_indx][0];
	for (i = 0; i < num_rngchkindex_cf[coupleg_indx]; i++)
	{
		rngchkindex_payT[0].push_back(cvg(today, rngchkindex_paydate[0][i], dcb));
		rngchkindex_tau[0].push_back(cvg(rngchkindex_cf[coupleg_indx][1][i], rngchkindex_cf[coupleg_indx][2][i], rngchkindex_dcb[coupleg_indx]));
		rngchkindext_payT[0].push_back(rngchkindex_payT[0][i] - t[0]);
		PMrngchkindext_payT[0].push_back(exp(rngchkindexzct(t[0]) * t[0] - rngchkindexzct(rngchkindex_payT[0][i]) * rngchkindex_payT[0][i]));
	}

	//added on 2017/02/16								
	ShiftBusDate(date[0], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_spotlag[coupleg_indx], rngchkindex2_firststartdate[0]);
	findmaturity(rngchkindex2_firststartdate[0], rngchkindex2_tenors[coupleg_indx], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_conv[coupleg_indx], rngchkindex2_lastpaydate[0]);
	rngchkindex2_firststartT[0] = cvg(today, rngchkindex2_firststartdate[0], dcb);
	rngchkindex2t_firststartT[0] = rngchkindex2_firststartT[0] - t[0];
	PMrngchkindex2t_firststartT[0] = exp(rngchkindex2zct(t[0]) * t[0] - rngchkindex2zct(rngchkindex2_firststartT[0]) * rngchkindex2_firststartT[0]);
	fixedlegcashflowschedule(rngchkindex2_firststartdate[0], rngchkindex2_tenors[coupleg_indx], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_stub[coupleg_indx], rngchkindex2_direction[coupleg_indx], rngchkindex2_conv[coupleg_indx], rngchkindex2fixed_freqs[coupleg_indx], rngchkindex2_adjflag[coupleg_indx], rngchkindex2_cf[coupleg_indx]);
	rngchkindex2_paydate[0] = rngchkindex2_cf[coupleg_indx][0];
	for (i = 0; i < num_rngchkindex2_cf[coupleg_indx]; i++)
	{
		rngchkindex2_payT[0].push_back(cvg(today, rngchkindex2_paydate[0][i], dcb));
		rngchkindex2_tau[0].push_back(cvg(rngchkindex2_cf[coupleg_indx][1][i], rngchkindex2_cf[coupleg_indx][2][i], rngchkindex2_dcb[coupleg_indx]));
		rngchkindex2t_payT[0].push_back(rngchkindex2_payT[0][i] - t[0]);
		PMrngchkindex2t_payT[0].push_back(exp(rngchkindex2zct(t[0]) * t[0] - rngchkindex2zct(rngchkindex2_payT[0][i]) * rngchkindex2_payT[0][i]));
	}

	double
		stdevx1_0
		, stdevx2_0
		, stdevy2_0
		, corrx1x2_0 = 0.0
		, corrx1y2_0 = 0.0
		, corrx2y2_0 = 0.0
		, ema1dt_0
		, ema2dt_0
		, emb2dt_0
		, qadjx2_0
		, qadjy2_0;

	vector<vector<double>> corr_0(num_factor);
	vector<vector<double>> l_mtrx_0;
	stdevx1_0 = sig1s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a1s[param_indx] * t[0])) / a1s[param_indx]);
	stdevx2_0 = sig2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a2s[param_indx] * t[0])) / a2s[param_indx]);
	stdevy2_0 = eta2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * b2s[param_indx] * t[0])) / b2s[param_indx]);
	if (stdevx1_0 * stdevx2_0 > 0) corrx1x2_0 = sig1s[param_indx] * sig2s[param_indx] * gamx1x2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + a2s[param_indx]) * t[0])) / (a1s[param_indx] + a2s[param_indx]) / (stdevx1_0 * stdevx2_0);
	if (stdevx1_0 * stdevy2_0 > 0) corrx1y2_0 = sig1s[param_indx] * eta2s[param_indx] * gamx1y2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + b2s[param_indx]) * t[0])) / (a1s[param_indx] + b2s[param_indx]) / (stdevx1_0 * stdevy2_0);
	if (stdevx2_0 * stdevy2_0 > 0) corrx2y2_0 = sig2s[param_indx] * eta2s[param_indx] * rho2s[param_indx] * (1.0 - exp(-(a2s[param_indx] + b2s[param_indx]) * t[0])) / (a2s[param_indx] + b2s[param_indx]) / (stdevx2_0 * stdevy2_0);

	for (j = 0; j < num_factor; j++) corr_0[j] = vector<double>(num_factor, 1.0);
	corr_0[0][1] = corrx1x2_0;
	corr_0[0][2] = corrx1y2_0;
	corr_0[1][2] = corrx2y2_0;

	corr_0[1][0] = corr_0[0][1];
	corr_0[2][0] = corr_0[0][2];
	corr_0[2][1] = corr_0[1][2];

	l_mtrx_0 = choldc(num_factor, corr_0);

	ema1dt_0 = exp(-a1s[param_indx] * t[0]);
	ema2dt_0 = exp(-a2s[param_indx] * t[0]);
	emb2dt_0 = exp(-b2s[param_indx] * t[0]);

	qadjx2_0 = -sig2s[param_indx] * nu2s[param_indx] * CX2x2s[param_indx] * (1.0 - ema2dt_0) / a2s[param_indx];
	qadjy2_0 = -eta2s[param_indx] * nu2s[param_indx] * CX2y2s[param_indx] * (1.0 - emb2dt_0) / b2s[param_indx];

	stdevx1[0] = sig1s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a1s[param_indx] * dt[0])) / a1s[param_indx]);
	stdevx2[0] = sig2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a2s[param_indx] * dt[0])) / a2s[param_indx]);
	stdevy2[0] = eta2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * b2s[param_indx] * dt[0])) / b2s[param_indx]);

	corrx1x2[0] = sig1s[param_indx] * sig2s[param_indx] * gamx1x2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + a2s[param_indx]) * dt[0])) / (a1s[param_indx] + a2s[param_indx]) / (stdevx1[0] * stdevx2[0]);
	corrx1y2[0] = sig1s[param_indx] * eta2s[param_indx] * gamx1y2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + b2s[param_indx]) * dt[0])) / (a1s[param_indx] + b2s[param_indx]) / (stdevx1[0] * stdevy2[0]);
	corrx2y2[0] = sig2s[param_indx] * eta2s[param_indx] * rho2s[param_indx] * (1.0 - exp(-(a2s[param_indx] + b2s[param_indx]) * dt[0])) / (a2s[param_indx] + b2s[param_indx]) / (stdevx2[0] * stdevy2[0]);

	corr[0] = vector<vector<double>>(num_factor);

	for (j = 0; j < num_factor; j++) corr[0][j] = vector<double>(num_factor, 1.0);
	corr[0][0][1] = corrx1x2[0];
	corr[0][0][2] = corrx1y2[0];
	corr[0][1][2] = corrx2y2[0];
	corr[0][1][0] = corr[0][0][1];
	corr[0][2][0] = corr[0][0][2];
	corr[0][2][1] = corr[0][1][2];

	l_mtrx[0] = choldc(num_factor, corr[0]);

	ema1dt[0] = exp(-a1s[param_indx] * dt[0]);
	ema2dt[0] = exp(-a2s[param_indx] * dt[0]);
	emb2dt[0] = exp(-b2s[param_indx] * dt[0]);

	qadjx2[0] = -sig2s[param_indx] * nu2s[param_indx] * CX2x2s[param_indx] * (1.0 - ema2dt[0]) / a2s[param_indx];
	qadjy2[0] = -eta2s[param_indx] * nu2s[param_indx] * CX2y2s[param_indx] * (1.0 - emb2dt[0]) / b2s[param_indx];

	Bt_firststartT[0] = (1.0 - exp(-a1s[param_indx] * t_firststartT[0])) / a1s[param_indx];
	Bt_lastpayT[0] = (1.0 - exp(-a1s[param_indx] * t_lastpayT[0])) / a1s[param_indx];

	//USD 30Y - 2F								
	Bindex1a2t_firststartT[0] = (1.0 - exp(-a2s[param_indx] * index1t_firststartT[0])) / a2s[param_indx];
	Bindex1b2t_firststartT[0] = (1.0 - exp(-b2s[param_indx] * index1t_firststartT[0])) / b2s[param_indx];
	for (i = 0; i < num_index1_cf[coupleg_indx]; i++)
	{
		Bindex1a2t_payT[0].push_back((1.0 - exp(-a2s[param_indx] * index1t_payT[0][i])) / a2s[param_indx]);
		Bindex1b2t_payT[0].push_back((1.0 - exp(-b2s[param_indx] * index1t_payT[0][i])) / b2s[param_indx]);
	}

	//USD 2Y - 2F								
	Bindex2a2t_firststartT[0] = (1.0 - exp(-a2s[param_indx] * index2t_firststartT[0])) / a2s[param_indx];
	Bindex2b2t_firststartT[0] = (1.0 - exp(-b2s[param_indx] * index2t_firststartT[0])) / b2s[param_indx];
	for (i = 0; i < num_index2_cf[coupleg_indx]; i++)
	{
		Bindex2a2t_payT[0].push_back((1.0 - exp(-a2s[param_indx] * index2t_payT[0][i])) / a2s[param_indx]);
		Bindex2b2t_payT[0].push_back((1.0 - exp(-b2s[param_indx] * index2t_payT[0][i])) / b2s[param_indx]);
	}

	//KRW 10Y - 1F								
	Brngchkindexa1t_firststartT[0] = (1.0 - exp(-a1s[param_indx] * rngchkindext_firststartT[0])) / a1s[param_indx];
	for (i = 0; i < num_rngchkindex_cf[coupleg_indx]; i++)
	{
		Brngchkindexa1t_payT[0].push_back((1.0 - exp(-a1s[param_indx] * rngchkindext_payT[0][i])) / a1s[param_indx]);
	}

	//USD 10Y - 2F								
	//added on 2017/02/16								
	Brngchkindex2a2t_firststartT[0] = (1.0 - exp(-a2s[param_indx] * rngchkindex2t_firststartT[0])) / a2s[param_indx];
	Brngchkindex2b2t_firststartT[0] = (1.0 - exp(-b2s[param_indx] * rngchkindex2t_firststartT[0])) / b2s[param_indx];
	for (i = 0; i < num_rngchkindex2_cf[coupleg_indx]; i++)
	{
		Brngchkindex2a2t_payT[0].push_back((1.0 - exp(-a2s[param_indx] * rngchkindex2t_payT[0][i])) / a2s[param_indx]);
		Brngchkindex2b2t_payT[0].push_back((1.0 - exp(-b2s[param_indx] * rngchkindex2t_payT[0][i])) / b2s[param_indx]);
	}

	//Vindex1t(num_rates),Vindex2t(num_rates),								
	coeft_firststartT[0] = PMt_firststartT[0] * exp(-Bt_firststartT[0] * (0.5 * pow(sig1s[param_indx] * (1.0 - exp(-a1s[param_indx] * t[0])) / a1s[param_indx], 2.0) + 0.25 * pow(sig1s[param_indx], 2.0) * (1.0 - exp(-2.0 * a1s[param_indx] * t[0])) / a1s[param_indx] * Bt_firststartT[0]));
	coeft_lastpayT[0] = PMt_lastpayT[0] * exp(-Bt_lastpayT[0] * (0.5 * pow(sig1s[param_indx] * (1.0 - exp(-a1s[param_indx] * t[0])) / a1s[param_indx], 2.0) + 0.25 * pow(sig1s[param_indx], 2.0) * (1.0 - exp(-2.0 * a1s[param_indx] * t[0])) / a1s[param_indx] * Bt_lastpayT[0]));

	Vindex1t[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (t[0] + (2.0 * exp(-a2s[param_indx] * t[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * t[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (t[0] + (2.0 * exp(-b2s[param_indx] * t[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * t[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (t[0] + (exp(-a2s[param_indx] * t[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * t[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * t[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	Vindex1_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1_firststartT[0] + (2.0 * exp(-a2s[param_indx] * index1_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1_firststartT[0] + (2.0 * exp(-b2s[param_indx] * index1_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1_firststartT[0] + (exp(-a2s[param_indx] * index1_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	Vindex1t_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1t_firststartT[0] + (2.0 * exp(-a2s[param_indx] * index1t_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1t_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1t_firststartT[0] + (2.0 * exp(-b2s[param_indx] * index1t_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1t_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1t_firststartT[0] + (exp(-a2s[param_indx] * index1t_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1t_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1t_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	coefindex1t_firststartT[0] = PMindex1t_firststartT[0] * exp(0.5 * (Vindex1t_firststartT[0] - Vindex1_firststartT[0] + Vindex1t[0]) - Bindex1a2t_firststartT[0] * qadjx2_0 - Bindex1b2t_firststartT[0] * qadjy2_0);

	for (j = 0; j < num_index1_cf[coupleg_indx]; j++)
	{
		Vindex1_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1_payT[0][j] + (2.0 * exp(-a2s[param_indx] * index1_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1_payT[0][j] + (2.0 * exp(-b2s[param_indx] * index1_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1_payT[0][j] + (exp(-a2s[param_indx] * index1_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		Vindex1t_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1t_payT[0][j] + (2.0 * exp(-a2s[param_indx] * index1t_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1t_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1t_payT[0][j] + (2.0 * exp(-b2s[param_indx] * index1t_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1t_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1t_payT[0][j] + (exp(-a2s[param_indx] * index1t_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1t_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1t_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		coefindex1t_payT[0].push_back(PMindex1t_payT[0][j] * exp(0.5 * (Vindex1t_payT[0][j] - Vindex1_payT[0][j] + Vindex1t[0]) - Bindex1a2t_payT[0][j] * qadjx2_0 - Bindex1b2t_payT[0][j] * qadjy2_0));
	}

	Vindex2t[0] = Vindex1t[0];//=pow(sig2s[param_indx]/a2s[param_indx],2.0)*(t[0]+(2.0*exp(-a2s[param_indx]*t[0])-0.5*exp(-2.0*a2s[param_indx]*t[0])-1.5)/a2s[param_indx]);								

	Vindex2_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2_firststartT[0] + (2.0 * exp(-a2s[param_indx] * index2_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2_firststartT[0] + (2.0 * exp(-b2s[param_indx] * index2_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2_firststartT[0] + (exp(-a2s[param_indx] * index2_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	Vindex2t_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2t_firststartT[0] + (2.0 * exp(-a2s[param_indx] * index2t_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2t_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2t_firststartT[0] + (2.0 * exp(-b2s[param_indx] * index2t_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2t_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2t_firststartT[0] + (exp(-a2s[param_indx] * index2t_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2t_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2t_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	coefindex2t_firststartT[0] = PMindex2t_firststartT[0] * exp(0.5 * (Vindex2t_firststartT[0] - Vindex2_firststartT[0] + Vindex2t[0]) - Bindex2a2t_firststartT[0] * qadjx2_0 - Bindex2b2t_firststartT[0] * qadjy2_0);

	for (j = 0; j < num_index2_cf[coupleg_indx]; j++)
	{
		Vindex2_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2_payT[0][j] + (2.0 * exp(-a2s[param_indx] * index2_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2_payT[0][j] + (2.0 * exp(-b2s[param_indx] * index2_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2_payT[0][j] + (exp(-a2s[param_indx] * index2_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		Vindex2t_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2t_payT[0][j] + (2.0 * exp(-a2s[param_indx] * index2t_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2t_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2t_payT[0][j] + (2.0 * exp(-b2s[param_indx] * index2t_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2t_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2t_payT[0][j] + (exp(-a2s[param_indx] * index2t_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2t_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2t_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		coefindex2t_payT[0].push_back(PMindex2t_payT[0][j] * exp(0.5 * (Vindex2t_payT[0][j] - Vindex2_payT[0][j] + Vindex2t[0]) - Bindex2a2t_payT[0][j] * qadjx2_0 - Bindex2b2t_payT[0][j] * qadjy2_0));
	}


	Vrngchkindext[0] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (t[0] + (2.0 * exp(-a1s[param_indx] * t[0]) - 0.5 * exp(-2.0 * a1s[param_indx] * t[0]) - 1.5) / a1s[param_indx]);

	Vrngchkindex_firststartT[0] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindex_firststartT[0] + (2.0 * exp(-a1s[param_indx] * rngchkindex_firststartT[0]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindex_firststartT[0]) - 1.5) / a1s[param_indx]);

	Vrngchkindext_firststartT[0] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindext_firststartT[0] + (2.0 * exp(-a1s[param_indx] * rngchkindext_firststartT[0]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindext_firststartT[0]) - 1.5) / a1s[param_indx]);

	coefrngchkindext_firststartT[0] = PMrngchkindext_firststartT[0] * exp(0.5 * (Vrngchkindext_firststartT[0] - Vrngchkindex_firststartT[0] + Vrngchkindext[0]));

	for (j = 0; j < num_rngchkindex_cf[coupleg_indx]; j++)
	{
		Vrngchkindex_payT[0].push_back(pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindex_payT[0][j] + (2.0 * exp(-a1s[param_indx] * rngchkindex_payT[0][j]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindex_payT[0][j]) - 1.5) / a1s[param_indx]));
		Vrngchkindext_payT[0].push_back(pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindext_payT[0][j] + (2.0 * exp(-a1s[param_indx] * rngchkindext_payT[0][j]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindext_payT[0][j]) - 1.5) / a1s[param_indx]));
		coefrngchkindext_payT[0].push_back(PMrngchkindext_payT[0][j] * exp(0.5 * (Vrngchkindext_payT[0][j] - Vrngchkindex_payT[0][j] + Vrngchkindext[0])));
	}


	Vrngchkindex2t[0] = Vindex2t[0];

	Vrngchkindex2_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2_firststartT[0] + (2.0 * exp(-a2s[param_indx] * rngchkindex2_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2_firststartT[0] + (2.0 * exp(-b2s[param_indx] * rngchkindex2_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2_firststartT[0] + (exp(-a2s[param_indx] * rngchkindex2_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	Vrngchkindex2t_firststartT[0] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2t_firststartT[0] + (2.0 * exp(-a2s[param_indx] * rngchkindex2t_firststartT[0]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2t_firststartT[0]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2t_firststartT[0] + (2.0 * exp(-b2s[param_indx] * rngchkindex2t_firststartT[0]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2t_firststartT[0]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2t_firststartT[0] + (exp(-a2s[param_indx] * rngchkindex2t_firststartT[0]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2t_firststartT[0]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2t_firststartT[0]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));

	coefrngchkindex2t_firststartT[0] = PMrngchkindex2t_firststartT[0] * exp(0.5 * (Vrngchkindex2t_firststartT[0] - Vrngchkindex2_firststartT[0] + Vrngchkindex2t[0]) - Brngchkindex2a2t_firststartT[0] * qadjx2_0 - Bindex2b2t_firststartT[0] * qadjy2_0);

	for (j = 0; j < num_rngchkindex2_cf[coupleg_indx]; j++)
	{
		Vrngchkindex2_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2_payT[0][j] + (2.0 * exp(-a2s[param_indx] * rngchkindex2_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2_payT[0][j] + (2.0 * exp(-b2s[param_indx] * rngchkindex2_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2_payT[0][j] + (exp(-a2s[param_indx] * rngchkindex2_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		Vrngchkindex2t_payT[0].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2t_payT[0][j] + (2.0 * exp(-a2s[param_indx] * rngchkindex2t_payT[0][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2t_payT[0][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2t_payT[0][j] + (2.0 * exp(-b2s[param_indx] * rngchkindex2t_payT[0][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2t_payT[0][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2t_payT[0][j] + (exp(-a2s[param_indx] * rngchkindex2t_payT[0][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2t_payT[0][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2t_payT[0][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
		coefrngchkindex2t_payT[0].push_back(PMrngchkindex2t_payT[0][j] * exp(0.5 * (Vrngchkindex2t_payT[0][j] - Vrngchkindex2_payT[0][j] + Vrngchkindex2t[0]) - Brngchkindex2a2t_payT[0][j] * qadjx2_0 - Brngchkindex2b2t_payT[0][j] * qadjy2_0));
	}

	for (i = 1; i < num_rates; i++)
	{
		date[i] = DateAdd('d', i, first_simulate_date);

		if (callstarti < num_calldate)
		{
			if (date[i] == callnotice_date[count_call])
			{
				i_call.push_back(i);
				if (count_call < num_calldate - 1) count_call = count_call + 1;
			}
		}
		if (num_remained_fundingleg_cf > 0)
		{
			if (date[i] == _fundingleg_calc_startdate[count_fund_calc_start])
			{
				i_fund_calc_start.push_back(i);
				if (count_fund_calc_start < num_remained_fundingleg_cf - 1) count_fund_calc_start = count_fund_calc_start + 1;
			}
			if (date[i] == _fixing_date[count_fixing])
			{
				i_fixing.push_back(i);
				if (count_fixing < num_remained_fundingleg_cf - 1) count_fixing = count_fixing + 1;
			}
			if (date[i] == _fundingleg_paydate[count_fund_pay])
			{
				fund_pay_df.push_back(exp(-zc(_fundingleg_paydate[count_fund_pay]) * fundinglegT[count_fund_pay]));
				i_fund_pay.push_back(i);
				if (count_fund_pay < num_remained_fundingleg_cf - 1) count_fund_pay = count_fund_pay + 1;
			}
		}
		if (date[i] == _couponleg_calc_shiftedstartdate[count_coup_calc_shifted_start])
		{
			i_coup_calc_shifted_start.push_back(i);
			if (count_coup_calc_shifted_start < num_remained_couponleg_cf - 1) count_coup_calc_shifted_start = count_coup_calc_shifted_start + 1;
		}
		if (date[i] == _couponleg_calc_shiftedenddate[count_coup_calc_shifted_end])
		{
			i_coup_calc_shifted_end.push_back(i);
			if (count_coup_calc_shifted_end < num_remained_couponleg_cf - 1) count_coup_calc_shifted_end = count_coup_calc_shifted_end + 1;
		}
		if (date[i] == _couponleg_paydate[count_coup_pay])
		{
			coup_pay_df.push_back(exp(-zc(_couponleg_paydate[count_coup_pay]) * couponlegT[count_coup_pay]));
			i_coup_pay.push_back(i);
			if (count_coup_pay < num_remained_couponleg_cf - 1) count_coup_pay = count_coup_pay + 1;
		}
		if (date[i] >= _couponleg_calc_shiftedstartdate[couponleg_shifted_cf_starti]) accrualstarti = min(i, accrualstarti);


		ShiftBusDate(date[i], fixing_holidays, num_fixingholidays, -fix_lag, fund_firststartdate[i]);
		findmaturity(fund_firststartdate[i], fundinglegindx_tenor, fixing_holidays, num_fixingholidays, fixingindex_conv, fund_lastpaydate[i]);

		t[i] = cvg(today, date[i], dcb);
		fund_firststartT[i] = cvg(today, fund_firststartdate[i], dcb);
		fund_lastpayT[i] = cvg(today, fund_lastpaydate[i], dcb);
		fund_tau[i] = cvg(fund_firststartdate[i], fund_lastpaydate[i], fundingleg_dcb);
		t_firststartT[i] = fund_firststartT[i] - t[i];
		t_lastpayT[i] = fund_lastpayT[i] - t[i];
		PMt_firststartT[i] = exp(zct(t[i]) * t[i] - zct(fund_firststartT[i]) * fund_firststartT[i]);
		PMt_lastpayT[i] = exp(zct(t[i]) * t[i] - zct(fund_lastpayT[i]) * fund_lastpayT[i]);

		if (date[i] > _couponleg_calc_shiftedenddate[coupleg_indx] && coupleg_indx < num_remained_couponleg_cf - 1) coupleg_indx = coupleg_indx + 1;

		ShiftBusDate(date[i], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_spotlag[coupleg_indx], index1_firststartdate[i]);
		findmaturity(index1_firststartdate[i], index1_tenors[coupleg_indx], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_conv[coupleg_indx], index1_lastpaydate[i]);
		index1_firststartT[i] = cvg(today, index1_firststartdate[i], dcb);
		index1t_firststartT[i] = index1_firststartT[i] - t[i];
		PMindex1t_firststartT[i] = exp(index1zct(t[i]) * t[i] - index1zct(index1_firststartT[i]) * index1_firststartT[i]);
		fixedlegcashflowschedule(index1_firststartdate[i], index1_tenors[coupleg_indx], index1_holidayss[coupleg_indx], num_index1holidayss[coupleg_indx], index1_stub[coupleg_indx], index1_direction[coupleg_indx], index1_conv[coupleg_indx], index1fixed_freqs[coupleg_indx], index1_adjflag[coupleg_indx], index1_cf[coupleg_indx]);
		index1_paydate[i] = index1_cf[coupleg_indx][0];
		for (j = 0; j < num_index1_cf[coupleg_indx]; j++)
		{
			index1_payT[i].push_back(cvg(today, index1_paydate[i][j], dcb));
			index1_tau[i].push_back(cvg(index1_cf[coupleg_indx][1][j], index1_cf[coupleg_indx][2][j], index1_dcb[coupleg_indx]));
			index1t_payT[i].push_back(index1_payT[i][j] - t[i]);
			PMindex1t_payT[i].push_back(exp(index1zct(t[i]) * t[i] - index1zct(index1_payT[i][j]) * index1_payT[i][j]));
		}

		ShiftBusDate(date[i], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_spotlag[coupleg_indx], index2_firststartdate[i]);
		findmaturity(index2_firststartdate[i], index2_tenors[coupleg_indx], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_conv[coupleg_indx], index2_lastpaydate[i]);
		index2_firststartT[i] = cvg(today, index2_firststartdate[i], dcb);
		index2t_firststartT[i] = index2_firststartT[i] - t[i];
		PMindex2t_firststartT[i] = exp(index2zct(t[i]) * t[i] - index2zct(index2_firststartT[i]) * index2_firststartT[i]);
		fixedlegcashflowschedule(index2_firststartdate[i], index2_tenors[coupleg_indx], index2_holidayss[coupleg_indx], num_index2holidayss[coupleg_indx], index2_stub[coupleg_indx], index2_direction[coupleg_indx], index2_conv[coupleg_indx], index2fixed_freqs[coupleg_indx], index2_adjflag[coupleg_indx], index2_cf[coupleg_indx]);
		index2_paydate[i] = index2_cf[coupleg_indx][0];
		for (j = 0; j < num_index2_cf[coupleg_indx]; j++)
		{
			index2_payT[i].push_back(cvg(today, index2_paydate[i][j], dcb));
			index2_tau[i].push_back(cvg(index2_cf[coupleg_indx][1][j], index2_cf[coupleg_indx][2][j], index2_dcb[coupleg_indx]));
			index2t_payT[i].push_back(index2_payT[i][j] - t[i]);
			PMindex2t_payT[i].push_back(exp(index2zct(t[i]) * t[i] - index2zct(index2_payT[i][j]) * index2_payT[i][j]));
		}

		ShiftBusDate(date[i], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_spotlag[coupleg_indx], rngchkindex_firststartdate[i]);
		findmaturity(rngchkindex_firststartdate[i], rngchkindex_tenors[coupleg_indx], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_conv[coupleg_indx], rngchkindex_lastpaydate[i]);
		rngchkindex_firststartT[i] = cvg(today, rngchkindex_firststartdate[i], dcb);
		rngchkindext_firststartT[i] = rngchkindex_firststartT[i] - t[i];
		PMrngchkindext_firststartT[i] = exp(rngchkindexzct(t[i]) * t[i] - rngchkindexzct(rngchkindex_firststartT[i]) * rngchkindex_firststartT[i]);
		fixedlegcashflowschedule(rngchkindex_firststartdate[i], rngchkindex_tenors[coupleg_indx], rngchkindex_holidayss[coupleg_indx], num_rngchkindexholidayss[coupleg_indx], rngchkindex_stub[coupleg_indx], rngchkindex_direction[coupleg_indx], rngchkindex_conv[coupleg_indx], rngchkindexfixed_freqs[coupleg_indx], rngchkindex_adjflag[coupleg_indx], rngchkindex_cf[coupleg_indx]);
		rngchkindex_paydate[i] = rngchkindex_cf[coupleg_indx][0];
		for (j = 0; j < num_rngchkindex_cf[coupleg_indx]; j++)
		{
			rngchkindex_payT[i].push_back(cvg(today, rngchkindex_paydate[i][j], dcb));
			rngchkindex_tau[i].push_back(cvg(rngchkindex_cf[coupleg_indx][1][j], rngchkindex_cf[coupleg_indx][2][j], rngchkindex_dcb[coupleg_indx]));
			rngchkindext_payT[i].push_back(rngchkindex_payT[i][j] - t[i]);
			PMrngchkindext_payT[i].push_back(exp(rngchkindexzct(t[i]) * t[i] - rngchkindexzct(rngchkindex_payT[i][j]) * rngchkindex_payT[i][j]));
		}

		ShiftBusDate(date[i], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_spotlag[coupleg_indx], rngchkindex2_firststartdate[i]);
		findmaturity(rngchkindex2_firststartdate[i], rngchkindex2_tenors[coupleg_indx], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_conv[coupleg_indx], rngchkindex2_lastpaydate[i]);
		rngchkindex2_firststartT[i] = cvg(today, rngchkindex2_firststartdate[i], dcb);
		rngchkindex2t_firststartT[i] = rngchkindex2_firststartT[i] - t[i];
		PMrngchkindex2t_firststartT[i] = exp(rngchkindex2zct(t[i]) * t[i] - rngchkindex2zct(rngchkindex2_firststartT[i]) * rngchkindex2_firststartT[i]);
		fixedlegcashflowschedule(rngchkindex2_firststartdate[i], rngchkindex2_tenors[coupleg_indx], rngchkindex2_holidayss[coupleg_indx], num_rngchkindex2holidayss[coupleg_indx], rngchkindex2_stub[coupleg_indx], rngchkindex2_direction[coupleg_indx], rngchkindex2_conv[coupleg_indx], rngchkindex2fixed_freqs[coupleg_indx], rngchkindex2_adjflag[coupleg_indx], rngchkindex2_cf[coupleg_indx]);
		rngchkindex2_paydate[i] = rngchkindex2_cf[coupleg_indx][0];
		for (j = 0; j < num_rngchkindex2_cf[coupleg_indx]; j++)
		{
			rngchkindex2_payT[i].push_back(cvg(today, rngchkindex2_paydate[i][j], dcb));
			rngchkindex2_tau[i].push_back(cvg(rngchkindex2_cf[coupleg_indx][1][j], rngchkindex2_cf[coupleg_indx][2][j], rngchkindex2_dcb[coupleg_indx]));
			rngchkindex2t_payT[i].push_back(rngchkindex2_payT[i][j] - t[i]);
			PMrngchkindex2t_payT[i].push_back(exp(rngchkindex2zct(t[i]) * t[i] - rngchkindex2zct(rngchkindex2_payT[i][j]) * rngchkindex2_payT[i][j]));
		}

		stdevx1[i] = stdevx1[i - 1];
		stdevx2[i] = stdevx2[i - 1];
		stdevy2[i] = stdevy2[i - 1];

		corrx1x2[i] = corrx1x2[i - 1];
		corrx1y2[i] = corrx1y2[i - 1];
		corrx2y2[i] = corrx2y2[i - 1];

		corr[i] = corr[i - 1];
		l_mtrx[i] = l_mtrx[i - 1];

		ema1dt[i] = ema1dt[i - 1];
		ema2dt[i] = ema2dt[i - 1];
		emb2dt[i] = emb2dt[i - 1];

		qadjx2[i] = qadjx2[i - 1];
		qadjy2[i] = qadjy2[i - 1];

		if (_ra_flag[coupleg_indx])
		{
			if (date[i] > _couponleg_calc_shiftedenddate[coupleg_indx] && coupleg_indx < num_remained_couponleg_cf - 1)
			{
				coupleg_indx = coupleg_indx + 1;
			}
		}
		else
		{
			if (date[i] > couponleg_calc_enddate[coupleg_indx] && coupleg_indx < num_remained_couponleg_cf - 1)
			{
				coupleg_indx = coupleg_indx + 1;
			}
		}

		if (date[i] > paramterm_date[param_indx])
		{
			if (param_indx < num_paramterm - 1)
			{
				param_indx = param_indx + 1;
				stdevx1[i] = sig1s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a1s[param_indx] * dt[i])) / a1s[param_indx]);
				stdevx2[i] = sig2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * a2s[param_indx] * dt[i])) / a2s[param_indx]);
				stdevy2[i] = eta2s[param_indx] / M_SQRT2 * sqrt((1.0 - exp(-2.0 * b2s[param_indx] * dt[i])) / b2s[param_indx]);

				corrx1x2[i] = sig1s[param_indx] * sig2s[param_indx] * gamx1x2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + a2s[param_indx]) * dt[i])) / (a1s[param_indx] + a2s[param_indx]) / (stdevx1[i] * stdevx2[i]);
				corrx1y2[i] = sig1s[param_indx] * eta2s[param_indx] * gamx1y2s[param_indx] * (1.0 - exp(-(a1s[param_indx] + b2s[param_indx]) * dt[i])) / (a1s[param_indx] + b2s[param_indx]) / (stdevx1[i] * stdevy2[i]);
				corrx2y2[i] = sig2s[param_indx] * eta2s[param_indx] * rho2s[param_indx] * (1.0 - exp(-(a2s[param_indx] + b2s[param_indx]) * dt[i])) / (a2s[param_indx] + b2s[param_indx]) / (stdevx2[i] * stdevy2[i]);

				corr[i] = vector<vector<double>>(num_factor);
				for (j = 0; j < num_factor; j++) corr[i][j] = vector<double>(num_factor, 1.0);
				corr[i][0][1] = corrx1x2[i];
				corr[i][0][2] = corrx1y2[i];
				corr[i][1][2] = corrx2y2[i];

				corr[i][1][0] = corr[i][0][1];
				corr[i][2][0] = corr[i][0][2];
				corr[i][2][1] = corr[i][1][2];

				l_mtrx[i] = choldc(num_factor, corr[i]);

				ema1dt[i] = exp(-a1s[param_indx] * dt[i]);
				ema2dt[i] = exp(-a2s[param_indx] * dt[i]);
				emb2dt[i] = exp(-b2s[param_indx] * dt[i]);

				qadjx2[i] = -sig2s[param_indx] * nu2s[param_indx] * CX2x2s[param_indx] * (1.0 - ema2dt[i]) / a2s[param_indx];
				qadjy2[i] = -eta2s[param_indx] * nu2s[param_indx] * CX2y2s[param_indx] * (1.0 - emb2dt[i]) / b2s[param_indx];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////									
		Bt_firststartT[i] = (1.0 - exp(-a1s[param_indx] * t_firststartT[i])) / a1s[param_indx];
		Bt_lastpayT[i] = (1.0 - exp(-a1s[param_indx] * t_lastpayT[i])) / a1s[param_indx];

		Vindex1_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1_firststartT[i] + (2.0 * exp(-a2s[param_indx] * index1_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1_firststartT[i] + (2.0 * exp(-b2s[param_indx] * index1_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1_firststartT[i] + (exp(-a2s[param_indx] * index1_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		for (j = 0; j < num_index1_cf[coupleg_indx]; j++) Vindex1_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1_payT[i][j] + (2.0 * exp(-a2s[param_indx] * index1_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1_payT[i][j] + (2.0 * exp(-b2s[param_indx] * index1_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1_payT[i][j] + (exp(-a2s[param_indx] * index1_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));

		Vindex2_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2_firststartT[i] + (2.0 * exp(-a2s[param_indx] * index2_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2_firststartT[i] + (2.0 * exp(-b2s[param_indx] * index2_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2_firststartT[i] + (exp(-a2s[param_indx] * index2_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		for (j = 0; j < num_index2_cf[coupleg_indx]; j++) Vindex2_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2_payT[i][j] + (2.0 * exp(-a2s[param_indx] * index2_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2_payT[i][j] + (2.0 * exp(-b2s[param_indx] * index2_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2_payT[i][j] + (exp(-a2s[param_indx] * index2_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));

		Vrngchkindex_firststartT[i] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindex_firststartT[i] + (2.0 * exp(-a1s[param_indx] * rngchkindex_firststartT[i]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindex_firststartT[i]) - 1.5) / a1s[param_indx]);
		for (j = 0; j < num_rngchkindex_cf[coupleg_indx]; j++) Vrngchkindex_payT[i].push_back(pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindex_payT[i][j] + (2.0 * exp(-a1s[param_indx] * rngchkindex_payT[i][j]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindex_payT[i][j]) - 1.5) / a1s[param_indx]));

		Vrngchkindex2_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2_firststartT[i] + (2.0 * exp(-a2s[param_indx] * rngchkindex2_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2_firststartT[i] + (2.0 * exp(-b2s[param_indx] * rngchkindex2_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2_firststartT[i] + (exp(-a2s[param_indx] * rngchkindex2_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		for (j = 0; j < num_rngchkindex2_cf[coupleg_indx]; j++) Vrngchkindex2_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2_payT[i][j] + (2.0 * exp(-a2s[param_indx] * rngchkindex2_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2_payT[i][j] + (2.0 * exp(-b2s[param_indx] * rngchkindex2_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2_payT[i][j] + (exp(-a2s[param_indx] * rngchkindex2_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));

		///////////////////////////////////////////////////////////////////////////////////									

		Bindex1a2t_firststartT[i] = (1.0 - exp(-a2s[param_indx] * index1t_firststartT[i])) / a2s[param_indx];
		Bindex1b2t_firststartT[i] = (1.0 - exp(-b2s[param_indx] * index1t_firststartT[i])) / b2s[param_indx];
		for (j = 0; j < num_index1_cf[coupleg_indx]; j++)
		{
			Bindex1a2t_payT[i].push_back((1.0 - exp(-a2s[param_indx] * index1t_payT[i][j])) / a2s[param_indx]);
			Bindex1b2t_payT[i].push_back((1.0 - exp(-b2s[param_indx] * index1t_payT[i][j])) / b2s[param_indx]);
		}

		Bindex2a2t_firststartT[i] = (1.0 - exp(-a2s[param_indx] * index2t_firststartT[i])) / a2s[param_indx];
		Bindex2b2t_firststartT[i] = (1.0 - exp(-b2s[param_indx] * index2t_firststartT[i])) / b2s[param_indx];
		for (j = 0; j < num_index2_cf[coupleg_indx]; j++)
		{
			Bindex2a2t_payT[i].push_back((1.0 - exp(-a2s[param_indx] * index2t_payT[i][j])) / a2s[param_indx]);
			Bindex2b2t_payT[i].push_back((1.0 - exp(-b2s[param_indx] * index2t_payT[i][j])) / b2s[param_indx]);
		}

		Brngchkindexa1t_firststartT[i] = (1.0 - exp(-a1s[param_indx] * rngchkindext_firststartT[i])) / a1s[param_indx];
		for (j = 0; j < num_rngchkindex_cf[coupleg_indx]; j++)
		{
			Brngchkindexa1t_payT[i].push_back((1.0 - exp(-a1s[param_indx] * rngchkindext_payT[i][j])) / a1s[param_indx]);
		}

		Brngchkindex2a2t_firststartT[i] = (1.0 - exp(-a2s[param_indx] * rngchkindex2t_firststartT[i])) / a2s[param_indx];
		Brngchkindex2b2t_firststartT[i] = (1.0 - exp(-b2s[param_indx] * rngchkindex2t_firststartT[i])) / b2s[param_indx];
		for (j = 0; j < num_rngchkindex2_cf[coupleg_indx]; j++)
		{
			Brngchkindex2a2t_payT[i].push_back((1.0 - exp(-a2s[param_indx] * rngchkindex2t_payT[i][j])) / a2s[param_indx]);
			Brngchkindex2b2t_payT[i].push_back((1.0 - exp(-b2s[param_indx] * rngchkindex2t_payT[i][j])) / b2s[param_indx]);
		}

		coeft_firststartT[i] = PMt_firststartT[i] * exp(-Bt_firststartT[i] * (0.5 * pow(sig1s[param_indx] * (1.0 - exp(-a1s[param_indx] * t[i])) / a1s[param_indx], 2.0) + 0.25 * pow(sig1s[param_indx], 2.0) * (1.0 - exp(-2.0 * a1s[param_indx] * t[i])) / a1s[param_indx] * Bt_firststartT[i]));
		coeft_lastpayT[i] = PMt_lastpayT[i] * exp(-Bt_lastpayT[i] * (0.5 * pow(sig1s[param_indx] * (1.0 - exp(-a1s[param_indx] * t[i])) / a1s[param_indx], 2.0) + 0.25 * pow(sig1s[param_indx], 2.0) * (1.0 - exp(-2.0 * a1s[param_indx] * t[i])) / a1s[param_indx] * Bt_lastpayT[i]));

		//Vt[i]=pow(sig2s[param_indx]/a2s[param_indx],2.0)*(t[i]+(2.0*exp(-a2s[param_indx]*t[i])-0.5*exp(-2.0*a2s[param_indx]*t[i])-1.5)/a2s[param_indx])+pow(eta2s[param_indx]/b2s[param_indx],2.0)*(t[i]+(2.0*exp(-b2s[param_indx]*t[i])-0.5*exp(-2.0*b2s[param_indx]*t[i])-1.5)/b2s[param_indx])+2.0*rho2s[param_indx]*sig2s[param_indx]*eta2s[param_indx]/(a2s[param_indx]*b2s[param_indx])*(t[i]+(exp(-a2s[param_indx]*t[i])-1.0)/a2s[param_indx]+(exp(-b2s[param_indx]*t[i])-1.0)/b2s[param_indx]-(exp(-(a2s[param_indx]+b2s[param_indx])*t[i])-1.0)/(a2s[param_indx]+b2s[param_indx]));							

		Vindex1t[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (t[i] + (2.0 * exp(-a2s[param_indx] * t[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * t[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (t[i] + (2.0 * exp(-b2s[param_indx] * t[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * t[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (t[i] + (exp(-a2s[param_indx] * t[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * t[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * t[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		Vindex1t_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1t_firststartT[i] + (2.0 * exp(-a2s[param_indx] * index1t_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1t_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1t_firststartT[i] + (2.0 * exp(-b2s[param_indx] * index1t_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1t_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1t_firststartT[i] + (exp(-a2s[param_indx] * index1t_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1t_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1t_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		coefindex1t_firststartT[i] = PMindex1t_firststartT[i] * exp(0.5 * (Vindex1t_firststartT[i] - Vindex1_firststartT[i] + Vindex1t[i]) - Bindex1a2t_firststartT[i] * qadjx2[i] - Bindex1b2t_firststartT[i] * qadjy2[i]);
		for (j = 0; j < num_index1_cf[coupleg_indx]; j++)
		{
			Vindex1t_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index1t_payT[i][j] + (2.0 * exp(-a2s[param_indx] * index1t_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index1t_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index1t_payT[i][j] + (2.0 * exp(-b2s[param_indx] * index1t_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index1t_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index1t_payT[i][j] + (exp(-a2s[param_indx] * index1t_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index1t_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index1t_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
			coefindex1t_payT[i].push_back(PMindex1t_payT[i][j] * exp(0.5 * (Vindex1t_payT[i][j] - Vindex1_payT[i][j] + Vindex1t[i]) - Bindex1a2t_payT[i][j] * qadjx2[i] - Bindex1b2t_payT[i][j] * qadjy2[i]));
		}

		Vindex2t[i] = Vindex1t[i];//pow(sig2s[param_indx]/a2s[param_indx],2.0)*(t[i]+(2.0*exp(-a2s[param_indx]*t[i])-0.5*exp(-2.0*a2s[param_indx]*t[i])-1.5)/a2s[param_indx]);							
		Vindex2t_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2t_firststartT[i] + (2.0 * exp(-a2s[param_indx] * index2t_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2t_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2t_firststartT[i] + (2.0 * exp(-b2s[param_indx] * index2t_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2t_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2t_firststartT[i] + (exp(-a2s[param_indx] * index2t_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2t_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2t_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		coefindex2t_firststartT[i] = PMindex2t_firststartT[i] * exp(0.5 * (Vindex2t_firststartT[i] - Vindex2_firststartT[i] + Vindex2t[i]) - Bindex2a2t_firststartT[i] * qadjx2[i] - Bindex2b2t_firststartT[i] * qadjy2[i]);
		for (j = 0; j < num_index2_cf[coupleg_indx]; j++)
		{
			Vindex2t_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (index2t_payT[i][j] + (2.0 * exp(-a2s[param_indx] * index2t_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * index2t_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (index2t_payT[i][j] + (2.0 * exp(-b2s[param_indx] * index2t_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * index2t_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (index2t_payT[i][j] + (exp(-a2s[param_indx] * index2t_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * index2t_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * index2t_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
			coefindex2t_payT[i].push_back(PMindex2t_payT[i][j] * exp(0.5 * (Vindex2t_payT[i][j] - Vindex2_payT[i][j] + Vindex2t[i]) - Bindex2a2t_payT[i][j] * qadjx2[i] - Bindex2b2t_payT[i][j] * qadjy2[i]));
		}

		Vrngchkindext[i] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (t[i] + (2.0 * exp(-a1s[param_indx] * t[i]) - 0.5 * exp(-2.0 * a1s[param_indx] * t[i]) - 1.5) / a1s[param_indx]);
		Vrngchkindext_firststartT[i] = pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindext_firststartT[i] + (2.0 * exp(-a1s[param_indx] * rngchkindext_firststartT[i]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindext_firststartT[i]) - 1.5) / a1s[param_indx]);
		coefrngchkindext_firststartT[i] = PMrngchkindext_firststartT[i] * exp(0.5 * (Vrngchkindext_firststartT[i] - Vrngchkindex_firststartT[i] + Vrngchkindext[i]));
		for (j = 0; j < num_rngchkindex_cf[coupleg_indx]; j++)
		{
			Vrngchkindext_payT[i].push_back(pow(sig1s[param_indx] / a1s[param_indx], 2.0) * (rngchkindext_payT[i][j] + (2.0 * exp(-a1s[param_indx] * rngchkindext_payT[i][j]) - 0.5 * exp(-2.0 * a1s[param_indx] * rngchkindext_payT[i][j]) - 1.5) / a1s[param_indx]));
			coefrngchkindext_payT[i].push_back(PMrngchkindext_payT[i][j] * exp(0.5 * (Vrngchkindext_payT[i][j] - Vrngchkindex_payT[i][j] + Vrngchkindext[i])));
		}

		Vrngchkindex2t[i] = Vindex2t[i];
		Vrngchkindex2t_firststartT[i] = pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2t_firststartT[i] + (2.0 * exp(-a2s[param_indx] * rngchkindex2t_firststartT[i]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2t_firststartT[i]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2t_firststartT[i] + (2.0 * exp(-b2s[param_indx] * rngchkindex2t_firststartT[i]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2t_firststartT[i]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2t_firststartT[i] + (exp(-a2s[param_indx] * rngchkindex2t_firststartT[i]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2t_firststartT[i]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2t_firststartT[i]) - 1.0) / (a2s[param_indx] + b2s[param_indx]));
		coefrngchkindex2t_firststartT[i] = PMrngchkindex2t_firststartT[i] * exp(0.5 * (Vrngchkindex2t_firststartT[i] - Vrngchkindex2_firststartT[i] + Vrngchkindex2t[i]) - Brngchkindex2a2t_firststartT[i] * qadjx2[i] - Brngchkindex2b2t_firststartT[i] * qadjy2[i]);
		for (j = 0; j < num_rngchkindex2_cf[coupleg_indx]; j++)
		{
			Vrngchkindex2t_payT[i].push_back(pow(sig2s[param_indx] / a2s[param_indx], 2.0) * (rngchkindex2t_payT[i][j] + (2.0 * exp(-a2s[param_indx] * rngchkindex2t_payT[i][j]) - 0.5 * exp(-2.0 * a2s[param_indx] * rngchkindex2t_payT[i][j]) - 1.5) / a2s[param_indx]) + pow(eta2s[param_indx] / b2s[param_indx], 2.0) * (rngchkindex2t_payT[i][j] + (2.0 * exp(-b2s[param_indx] * rngchkindex2t_payT[i][j]) - 0.5 * exp(-2.0 * b2s[param_indx] * rngchkindex2t_payT[i][j]) - 1.5) / b2s[param_indx]) + 2.0 * rho2s[param_indx] * sig2s[param_indx] * eta2s[param_indx] / (a2s[param_indx] * b2s[param_indx]) * (rngchkindex2t_payT[i][j] + (exp(-a2s[param_indx] * rngchkindex2t_payT[i][j]) - 1.0) / a2s[param_indx] + (exp(-b2s[param_indx] * rngchkindex2t_payT[i][j]) - 1.0) / b2s[param_indx] - (exp(-(a2s[param_indx] + b2s[param_indx]) * rngchkindex2t_payT[i][j]) - 1.0) / (a2s[param_indx] + b2s[param_indx])));
			coefrngchkindex2t_payT[i].push_back(PMrngchkindex2t_payT[i][j] * exp(0.5 * (Vrngchkindex2t_payT[i][j] - Vrngchkindex2_payT[i][j] + Vrngchkindex2t[i]) - Brngchkindex2a2t_payT[i][j] * qadjx2[i] - Brngchkindex2b2t_payT[i][j] * qadjy2[i]));
		}
	}

	int seed = 1;
	srand(seed);
	int Npath = int(npath);
	vector<vector<double>> r(Npath), index1_r(Npath), index2_r(Npath), rngchkindex_r(Npath), rngchkindex2_r(Npath);
	vector<double> ex(num_factor), xv(num_factor);
	double index1_level, index2_level, rngchkindex_level, rngchkindex2_level;
	for (i = 0; i < npath; i++)
	{
		r[i] = vector<double>(num_rates);
		index1_r[i] = vector<double>(num_rates);
		index2_r[i] = vector<double>(num_rates);
		rngchkindex_r[i] = vector<double>(num_rates);
		rngchkindex2_r[i] = vector<double>(num_rates);

		index1_level = 0.0;
		index2_level = 0.0;
		rngchkindex_level = 0.0;
		rngchkindex2_level = 0.0;

		for (k = 0; k < num_factor; k++) ex[k] = snrnd();

		xv[0] = stdevx1_0 * ex[0];
		xv[1] = stdevx1_0 * l_mtrx_0[1][0] * ex[0] + stdevx2_0 * l_mtrx_0[1][1] * ex[1];
		xv[2] = stdevx1_0 * l_mtrx_0[2][0] * ex[0] + stdevx2_0 * l_mtrx_0[2][1] * ex[1] + stdevy2_0 * l_mtrx_0[2][2] * ex[2];

		r[i][0] = (coeft_firststartT[0] * exp(-Bt_firststartT[0] * xv[0]) / (coeft_lastpayT[0] * exp(-Bt_lastpayT[0] * xv[0])) - 1.0) / fund_tau[0];

		level_HW2F(num_index1_cf[coupleg_indx], index1_tau[0], coefindex1t_payT[0], Bindex1a2t_payT[0], Bindex1b2t_payT[0], xv[1], xv[2], index1_level);
		index1_r[i][0] = (coefindex1t_firststartT[0] * exp(-Bindex1a2t_firststartT[0] * xv[1] - Bindex1b2t_firststartT[0] * xv[2]) - coefindex1t_payT[0][num_index1_cf[coupleg_indx] - 1] * exp(-Bindex1a2t_payT[0][num_index1_cf[coupleg_indx] - 1] * xv[1] - Bindex1b2t_payT[0][num_index1_cf[coupleg_indx] - 1] * xv[2])) / index1_level;

		level_HW2F(num_index2_cf[coupleg_indx], index2_tau[0], coefindex2t_payT[0], Bindex2a2t_payT[0], Bindex2b2t_payT[0], xv[1], xv[2], index2_level);
		index2_r[i][0] = (coefindex2t_firststartT[0] * exp(-Bindex2a2t_firststartT[0] * xv[1] - Bindex2b2t_firststartT[0] * xv[2]) - coefindex2t_payT[0][num_index2_cf[coupleg_indx] - 1] * exp(-Bindex2a2t_payT[0][num_index2_cf[coupleg_indx] - 1] * xv[1] - Bindex2b2t_payT[0][num_index2_cf[coupleg_indx] - 1] * xv[2])) / index2_level;

		level_HW1F(num_rngchkindex_cf[coupleg_indx], rngchkindex_tau[0], coefrngchkindext_payT[0], Brngchkindexa1t_payT[0], xv[0], rngchkindex_level);
		rngchkindex_r[i][0] = (coefrngchkindext_firststartT[0] * exp(-Brngchkindexa1t_firststartT[0] * xv[0]) - coefrngchkindext_payT[0][num_rngchkindex_cf[coupleg_indx] - 1] * exp(-Brngchkindexa1t_payT[0][num_rngchkindex_cf[coupleg_indx] - 1] * xv[0])) / rngchkindex_level;

		level_HW2F(num_rngchkindex2_cf[coupleg_indx], rngchkindex2_tau[0], coefrngchkindex2t_payT[0], Brngchkindex2a2t_payT[0], Brngchkindex2b2t_payT[0], xv[1], xv[2], rngchkindex2_level);
		rngchkindex2_r[i][0] = (coefrngchkindex2t_firststartT[0] * exp(-Brngchkindex2a2t_firststartT[0] * xv[1] - Brngchkindex2b2t_firststartT[0] * xv[2]) - coefrngchkindex2t_payT[0][num_rngchkindex2_cf[coupleg_indx] - 1] * exp(-Brngchkindex2a2t_payT[0][num_rngchkindex2_cf[coupleg_indx] - 1] * xv[1] - Brngchkindex2b2t_payT[0][num_rngchkindex2_cf[coupleg_indx] - 1] * xv[2])) / rngchkindex2_level;

		alpha1[0] = fM1[0];

		for (j = 1; j < num_rates; j++)
		{
			index1_level = 0.0;
			index2_level = 0.0;
			rngchkindex_level = 0.0;
			rngchkindex2_level = 0.0;
			for (k = 0; k < num_factor; k++) ex[k] = snrnd();

			xv[0] = xv[0] * ema1dt[j] + stdevx1[j] * ex[0];
			xv[1] = xv[1] * ema2dt[j] + stdevx1[j] * l_mtrx[j][1][0] * ex[0] + stdevx2[j] * l_mtrx[j][1][1] * ex[1];
			xv[2] = xv[2] * emb2dt[j] + stdevx1[j] * l_mtrx[j][2][0] * ex[0] + stdevx2[j] * l_mtrx[j][2][1] * ex[1] + stdevy2[j] * l_mtrx[j][2][2] * ex[2];

			r[i][j] = (coeft_firststartT[j] * exp(-Bt_firststartT[j] * xv[0]) / (coeft_lastpayT[j] * exp(-Bt_lastpayT[j] * xv[0])) - 1.0) / fund_tau[j];

			level_HW2F(num_index1_cf[coupleg_indx], index1_tau[j], coefindex1t_payT[j], Bindex1a2t_payT[j], Bindex1b2t_payT[j], xv[1], xv[2], index1_level);
			index1_r[i][j] = (coefindex1t_firststartT[j] * exp(-Bindex1a2t_firststartT[j] * xv[1] - Bindex1b2t_firststartT[j] * xv[2]) - coefindex1t_payT[j][num_index1_cf[coupleg_indx] - 1] * exp(-Bindex1a2t_payT[j][num_index1_cf[coupleg_indx] - 1] * xv[1] - Bindex1b2t_payT[j][num_index1_cf[coupleg_indx] - 1] * xv[2])) / index1_level;

			level_HW2F(num_index2_cf[coupleg_indx], index2_tau[j], coefindex2t_payT[j], Bindex2a2t_payT[j], Bindex2b2t_payT[j], xv[1], xv[2], index2_level);
			index2_r[i][j] = (coefindex2t_firststartT[j] * exp(-Bindex2a2t_firststartT[j] * xv[1] - Bindex2b2t_firststartT[j] * xv[2]) - coefindex2t_payT[j][num_index2_cf[coupleg_indx] - 1] * exp(-Bindex2a2t_payT[j][num_index2_cf[coupleg_indx] - 1] * xv[1] - Bindex2b2t_payT[j][num_index2_cf[coupleg_indx] - 1] * xv[2])) / index2_level;

			level_HW1F(num_rngchkindex_cf[coupleg_indx], rngchkindex_tau[j], coefrngchkindext_payT[j], Brngchkindexa1t_payT[j], xv[0], rngchkindex_level);
			rngchkindex_r[i][j] = (coefrngchkindext_firststartT[j] * exp(-Brngchkindexa1t_firststartT[j] * xv[0]) - coefrngchkindext_payT[j][num_rngchkindex_cf[coupleg_indx] - 1] * exp(-Brngchkindexa1t_payT[j][num_rngchkindex_cf[coupleg_indx] - 1] * xv[0])) / rngchkindex_level;

			level_HW2F(num_rngchkindex2_cf[coupleg_indx], rngchkindex2_tau[j], coefrngchkindex2t_payT[j], Brngchkindex2a2t_payT[j], Brngchkindex2b2t_payT[j], xv[1], xv[2], rngchkindex2_level);
			rngchkindex2_r[i][j] = (coefrngchkindex2t_firststartT[j] * exp(-Brngchkindex2a2t_firststartT[j] * xv[1] - Brngchkindex2b2t_firststartT[j] * xv[2]) - coefrngchkindex2t_payT[j][num_rngchkindex2_cf[coupleg_indx] - 1] * exp(-Brngchkindex2a2t_payT[j][num_rngchkindex2_cf[coupleg_indx] - 1] * xv[1] - Brngchkindex2b2t_payT[j][num_rngchkindex2_cf[coupleg_indx] - 1] * xv[2])) / rngchkindex2_level;
		}
	}

	vector<int> aftercall_fund_starti(num_calldate - callstarti, 0), aftercall_coup_starti(num_calldate - callstarti, 0);
	vector<bool> _callability_flag(num_calldate - callstarti);
	for (i = callstarti; i < num_calldate; i++)
	{
		_callability_flag[i - callstarti] = callability_flag[i];
		while (calleffective_date[i] > _fundingleg_calc_startdate[aftercall_fund_starti[i - callstarti]])
		{
			aftercall_fund_starti[i - callstarti] = aftercall_fund_starti[i - callstarti] + 1;
			if (aftercall_fund_starti[i - callstarti] >= num_remained_fundingleg_cf - 1) break;
		}
		while (calleffective_date[i] > _couponleg_calc_startdate[aftercall_coup_starti[i - callstarti]])
		{
			aftercall_coup_starti[i - callstarti] = aftercall_coup_starti[i - callstarti] + 1;
			if (aftercall_coup_starti[i - callstarti] >= num_remained_couponleg_cf - 1) break;
		}
		if (_fundingleg_calc_startdate[aftercall_fund_starti[i - callstarti]] > _couponleg_calc_startdate[aftercall_coup_starti[i - callstarti]])
		{
			while (_fundingleg_calc_startdate[aftercall_fund_starti[i - callstarti]] > _couponleg_calc_startdate[aftercall_coup_starti[i - callstarti]])
			{
				aftercall_coup_starti[i - callstarti] = aftercall_coup_starti[i - callstarti] + 1;
				if (aftercall_coup_starti[i - callstarti] >= num_remained_couponleg_cf - 1) break;
			}
		}
		else if (_fundingleg_calc_startdate[aftercall_fund_starti[i - callstarti]] < _couponleg_calc_startdate[aftercall_coup_starti[i - callstarti]])
		{
			while (_fundingleg_calc_startdate[aftercall_fund_starti[i - callstarti]] < _couponleg_calc_startdate[aftercall_coup_starti[i - callstarti]])
			{
				aftercall_fund_starti[i - callstarti] = aftercall_fund_starti[i - callstarti] + 1;
				if (aftercall_fund_starti[i - callstarti] >= num_remained_fundingleg_cf - 1) break;
			}
		}
	}

	//double adj0=0.000;//20110124~20110331,20111230								
	//double adj0=0.03;//20110630 //adj0=0.007 ; 20130319								
	double adjc = 1.0;
	//double adj0=0.03;								
	double adj0 = 0.0;
	vector<vector<double>> aftercall_fund_price(num_calldate - callstarti), aftercall_coup_price(num_calldate - callstarti), hX(num_calldate - callstarti), uptocall_coup(num_calldate - callstarti + 1), call_coup(num_calldate - callstarti);
	vector<vector<vector<double>>> X(num_calldate - callstarti);
	uptocall_coup[0] = vector<double>(Npath, 0.0);
	for (k = 0; k < num_calldate - callstarti; k++)
	{
		aftercall_fund_price[k] = vector<double>(Npath, 0.0);
		aftercall_coup_price[k] = vector<double>(Npath, 0.0);
		hX[k] = vector<double>(Npath, 0.0);
		X[k] = vector<vector<double>>(Npath);
		for (i = 0; i < npath; i++) X[k][i] = vector<double>(num_factor);
		uptocall_coup[k + 1] = vector<double>(Npath, 0.0);
		call_coup[k] = vector<double>(Npath, 0.0);
	}
	double tmp_fund_price = 0.0, tmp_coup_price = 0.0, tmp_ratein_count = 0.0, aftercall_coup_price_avg = 0.0;
	int fixing_indx_adj = max(0, num_remained_fundingleg_cf - int(i_fixing.size())), fundpay_indx_adj = max(0, num_remained_fundingleg_cf - int(i_fund_pay.size())), start_indx_adj = max(0, num_remained_couponleg_cf - int(i_coup_calc_shifted_start.size())), end_indx_adj = max(0, num_remained_couponleg_cf - int(i_coup_calc_shifted_end.size())), pay_indx_adj = max(0, num_remained_couponleg_cf - int(i_coup_pay.size()));

	if (num_calldate - callstarti > 0)
	{
		for (i = 0; i < npath; i++)
		{
			for (k = aftercall_fund_starti[num_calldate - callstarti - 1]; k < num_remained_fundingleg_cf; k++) aftercall_fund_price[num_calldate - callstarti - 1][i] = aftercall_fund_price[num_calldate - callstarti - 1][i] + (r[i][i_fixing[k - fixing_indx_adj]] + _fundingleg_spread[k]) * _fundingleg_notional[k] * fund_pay_df[k - fundpay_indx_adj] * fundinglegtau[k];
			for (k = aftercall_coup_starti[num_calldate - callstarti - 1]; k < num_remained_couponleg_cf; k++)
			{
				tmp_ratein_count = 0.0;
				ratecount[k] = 0.0;
				for (l = i_coup_calc_shifted_start[k - start_indx_adj]; l < i_coup_calc_shifted_end[k - end_indx_adj]; l++)
				{
					if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[k], _ralow_bdry[k], _rahigh_bdryin_flag[k], _ralow_bdryin_flag[k])) tmp_ratein_count = tmp_ratein_count + 1.0;
					ratecount[k] = ratecount[k] + _couponlegindex1_mult[k] * index1_r[i][l] + _couponlegindex2_mult[k] * index2_r[i][l] + _couponleg_spread[k] + adj0;

				}
				//call_coup[num_calldate-callstarti-1][i]=min(max(ratecount[k]/double(num_fixed[k]+num_fixing[k]),_floorrates[k]),_caprates[k])*_couponleg_notional[k]*couponlegtau[k]*tmp_ratein_count/double(num_fixed[k]+num_fixing[k]);					
				call_coup[num_calldate - callstarti - 1][i] = min(max(ratecount[k] / double(num_fixing[k]), _floorrates[k]), _caprates[k]) * _couponleg_notional[k] * coup_pay_df[k - pay_indx_adj] * couponlegtau[k] * tmp_ratein_count / double(num_fixing[k]);
				aftercall_coup_price[num_calldate - callstarti - 1][i] = aftercall_coup_price[num_calldate - callstarti - 1][i] + call_coup[num_calldate - callstarti - 1][i];
			}
			aftercall_coup_price_avg = aftercall_coup_price_avg + aftercall_coup_price[num_calldate - callstarti - 1][i];
			hX[num_calldate - callstarti - 1][i] = aftercall_coup_price[num_calldate - callstarti - 1][i] - aftercall_fund_price[num_calldate - callstarti - 1][i];
			//hX[num_calldate-callstarti-1][i]=aftercall_coup_price[num_calldate-callstarti-1][i]*coup_pay_df[num_remained_couponleg_cf-1-pay_indx_adj]-aftercall_fund_price[num_calldate-callstarti-1][i];						
			hX[num_calldate - callstarti - 1][i] = hX[num_calldate - callstarti - 1][i] * adjc;
			X[num_calldate - callstarti - 1][i][0] = r[i][i_call[num_calldate - callstarti - 1]];
			X[num_calldate - callstarti - 1][i][1] = index1_r[i][i_call[num_calldate - callstarti - 1]];
			X[num_calldate - callstarti - 1][i][2] = index2_r[i][i_call[num_calldate - callstarti - 1]];
		}
	}
	for (j = num_calldate - 2; j >= callstarti; j--)
	{
		for (i = 0; i < npath; i++)
		{

			aftercall_fund_price[j - callstarti][i] = aftercall_fund_price[j - callstarti + 1][i] + aftercall_fund_price[j - callstarti][i];
			for (k = aftercall_fund_starti[j - callstarti]; k < aftercall_fund_starti[j - callstarti + 1]; k++) aftercall_fund_price[j - callstarti][i] = aftercall_fund_price[j - callstarti][i] + (r[i][i_fixing[k - fixing_indx_adj]] + _fundingleg_spread[k]) * _fundingleg_notional[k] * fund_pay_df[k - fundpay_indx_adj] * fundinglegtau[k];
			for (k = aftercall_coup_starti[j - callstarti]; k < aftercall_coup_starti[j - callstarti + 1]; k++)
			{

				tmp_ratein_count = 0.0;
				ratecount[k] = 0.0;
				if (k - start_indx_adj >= 0) {		//반기인 경우 indx 처리시 에러 발생 accrued에서 이미 해줬으므로 빼고,, 향후에 추가 검토			
					for (l = i_coup_calc_shifted_start[k - start_indx_adj]; l < i_coup_calc_shifted_end[k - end_indx_adj]; l++)
					{
						if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[k], _ralow_bdry[k], _rahigh_bdryin_flag[k], _ralow_bdryin_flag[k])) tmp_ratein_count = tmp_ratein_count + 1.0;
						ratecount[k] = ratecount[k] + _couponlegindex1_mult[k] * index1_r[i][l] + _couponlegindex2_mult[k] * index2_r[i][l] + _couponleg_spread[k] + adj0;

					}
					//aftercall_coup_price[j-callstarti][i]=aftercall_coup_price[j-callstarti][i]+min(max(ratecount[k]/double(num_fixed[k]+num_fixing[k]),_floorrates[k]),_caprates[k])*_couponleg_notional[k]*coup_pay_df[k-pay_indx_adj]*couponlegtau[k]*tmp_ratein_count/double(num_fixed[k]+num_fixing[k]);				
					//call_coup[j-callstarti][i]=min(max(ratecount[k]/double(num_fixed[k]+num_fixing[k]),_floorrates[k]),_caprates[k])*_couponleg_notional[k]*couponlegtau[k]*tmp_ratein_count/double(num_fixed[k]+num_fixing[k]);				
					call_coup[j - callstarti][i] = min(max(ratecount[k] / double(num_fixing[k]), _floorrates[k]), _caprates[k]) * _couponleg_notional[k] * coup_pay_df[k - pay_indx_adj] * couponlegtau[k] * tmp_ratein_count / double(num_fixing[k]);
					aftercall_coup_price[j - callstarti][i] = aftercall_coup_price[j - callstarti][i] + call_coup[j - callstarti][i];
				}
			}
			aftercall_coup_price_avg = aftercall_coup_price_avg + aftercall_coup_price[j - callstarti][i];
			aftercall_coup_price[j - callstarti][i] = aftercall_coup_price[j - callstarti][i] + aftercall_coup_price[j - callstarti + 1][i];

			hX[j - callstarti][i] = aftercall_coup_price[j - callstarti][i] - aftercall_fund_price[j - callstarti][i];
			//hX[j-callstarti][i]=aftercall_coup_price[j-callstarti][i]*coup_pay_df[num_remained_couponleg_cf-1-pay_indx_adj]-aftercall_fund_price[j-callstarti][i];						
			hX[j - callstarti][i] = hX[j - callstarti][i] * adjc;
			X[j - callstarti][i][0] = r[i][i_call[j - callstarti]];
			X[j - callstarti][i][1] = index1_r[i][i_call[j - callstarti]];
			X[j - callstarti][i][2] = index2_r[i][i_call[j - callstarti]];
		}
	}
	aftercall_coup_price_avg = aftercall_coup_price_avg / npath;
	//aftercall_coup_price_avg=aftercall_coup_price_avg*coup_pay_df[num_remained_couponleg_cf-1-pay_indx_adj]/npath;								

	//double adj=0;//20110124~20110331,20111230								
	//double adj=0.03;//20110630								
	double adj = 0.00;
	double fixedrateincount = 0.0, tempprice = 0.0;
	if (i_coup_calc_shifted_end.size() > 0)
	{
		if (num_calldate > 0 && callstarti > 0)
		{
			if (today > callnotice_date[callstarti - 1])
			{
				for (i = 0; i < npath; i++)
				{
					fixedrateincount = 0.0;
					ratecount[end_indx_adj] = 0.0;
					for (l = accrualstarti; l < i_coup_calc_shifted_end[0]; l++)
					{
						if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj], _ralow_bdry[end_indx_adj], _rahigh_bdryin_flag[end_indx_adj], _ralow_bdryin_flag[end_indx_adj])) fixedrateincount = fixedrateincount + 1.0;
						ratecount[end_indx_adj] = ratecount[end_indx_adj] + _couponlegindex1_mult[end_indx_adj] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj] * index2_r[i][l] + _couponleg_spread[end_indx_adj] + adj;
					}
					tempprice = tempprice + min(max(ratecount[end_indx_adj] / double(num_fixing[end_indx_adj]), _floorrates[end_indx_adj]), _caprates[end_indx_adj]) * fixedrateincount / double(num_fixing[end_indx_adj]);
				}
				if (coup_pay_df.size() > 0) fixedlegnonfixedprice = tempprice * _couponleg_notional[end_indx_adj] * coup_pay_df[0] * couponlegtau[end_indx_adj] / npath;
				//if(coup_pay_df.size()>0) fixedlegnonfixedprice=tempprice*_couponleg_notional[end_indx_adj]*couponlegtau[end_indx_adj]/npath;					
				if (_couponleg_calc_startdate.size() > 1 && i_coup_calc_shifted_end.size() > 1)
				{
					if (today < _couponleg_calc_shiftedstartdate[1] && _couponleg_calc_shiftedstartdate[1] < callnotice_date[callstarti])
					{
						for (k = 0; k < aftercall_coup_starti[0]; k++)
						{
							//////////////////////////////////////////////									
							for (i = 0; i < npath; i++)
							{
								fixedrateincount = 0.0;
								ratecount[end_indx_adj + 1 + k] = 0.0;
								for (l = i_coup_calc_shifted_end[0 + k]; l < i_coup_calc_shifted_end[1 + k]; l++)
								{
									if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj + 1], _ralow_bdry[end_indx_adj + 1], _rahigh_bdryin_flag[end_indx_adj + 1], _ralow_bdryin_flag[end_indx_adj + 1])) fixedrateincount = fixedrateincount + 1.0;
									ratecount[end_indx_adj + 1 + k] = ratecount[end_indx_adj + 1 + k] + _couponlegindex1_mult[end_indx_adj + 1 + k] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj + 1 + k] * index2_r[i][l] + _couponleg_spread[end_indx_adj + 1 + k] + adj;
								}
								tempprice = tempprice + min(max(ratecount[end_indx_adj + 1 + k] / double(num_fixing[end_indx_adj + 1 + k]), _floorrates[end_indx_adj + 1 + k]), _caprates[end_indx_adj + 1 + k]) * fixedrateincount / double(num_fixing[end_indx_adj + 1 + k]);
							}
							if (coup_pay_df.size() > 1) fixedlegnonfixedprice = fixedlegnonfixedprice + tempprice * _couponleg_notional[end_indx_adj + 1 + k] * coup_pay_df[1 + k] * couponlegtau[end_indx_adj + 1 + k] / npath;
							//if(coup_pay_df.size()>1) fixedlegnonfixedprice=fixedlegnonfixedprice+tempprice*_couponleg_notional[end_indx_adj+1+k]*couponlegtau[end_indx_adj+1+k]/npath;		
							tempprice = 0.0;
						}
					}
				}
			}
		}
		else if (num_calldate > 0 && callstarti <= 0)
		{
			for (i = 0; i < npath; i++)
			{
				fixedrateincount = 0.0;
				ratecount[end_indx_adj] = 0.0;
				for (l = accrualstarti; l < i_coup_calc_shifted_end[0]; l++)
				{
					if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj], _ralow_bdry[end_indx_adj], _rahigh_bdryin_flag[end_indx_adj], _ralow_bdryin_flag[end_indx_adj])) fixedrateincount = fixedrateincount + 1.0;
					ratecount[end_indx_adj] = ratecount[end_indx_adj] + _couponlegindex1_mult[end_indx_adj] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj] * index2_r[i][l] + _couponleg_spread[end_indx_adj] + adj;
				}
				tempprice = tempprice + min(max(ratecount[end_indx_adj] / double(num_fixing[end_indx_adj]), _floorrates[end_indx_adj]), _caprates[end_indx_adj]) * fixedrateincount / double(num_fixing[end_indx_adj]);
			}
			if (coup_pay_df.size() > 0) fixedlegnonfixedprice = tempprice * _couponleg_notional[end_indx_adj] * coup_pay_df[0] * couponlegtau[end_indx_adj] / npath;
			//if(coup_pay_df.size()>0) fixedlegnonfixedprice=tempprice*_couponleg_notional[end_indx_adj]*couponlegtau[end_indx_adj]/npath;						
			if (_couponleg_calc_startdate.size() > 1 && i_coup_calc_shifted_end.size() > 1)
			{
				if (today < _couponleg_calc_shiftedstartdate[1] && _couponleg_calc_shiftedstartdate[1] < callnotice_date[callstarti])
				{
					for (k = 0; k < aftercall_coup_starti[0] && _ra_flag[end_indx_adj + 1 + k]; k++)
					{
						//////////////////////////////////////////////									
						for (i = 0; i < npath; i++)
						{
							fixedrateincount = 0.0;
							ratecount[end_indx_adj + 1 + k] = 0.0;
							for (l = i_coup_calc_shifted_end[0 + k]; l < i_coup_calc_shifted_end[1 + k]; l++)
							{
								if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj + 1], _ralow_bdry[end_indx_adj + 1], _rahigh_bdryin_flag[end_indx_adj + 1], _ralow_bdryin_flag[end_indx_adj + 1])) fixedrateincount = fixedrateincount + 1.0;
								ratecount[end_indx_adj + 1 + k] = ratecount[end_indx_adj + 1 + k] + _couponlegindex1_mult[end_indx_adj + 1 + k] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj + 1 + k] * index2_r[i][l] + _couponleg_spread[end_indx_adj + 1 + k] + adj;
							}
							tempprice = tempprice + min(max(ratecount[end_indx_adj + 1 + k] / double(num_fixing[end_indx_adj + 1 + k]), _floorrates[end_indx_adj + 1 + k]), _caprates[end_indx_adj + 1 + k]) * fixedrateincount / double(num_fixing[end_indx_adj + 1 + k]);
						}
						if (coup_pay_df.size() > 1) fixedlegnonfixedprice = fixedlegnonfixedprice + tempprice * _couponleg_notional[end_indx_adj + 1 + k] * coup_pay_df[1 + k] * couponlegtau[end_indx_adj + 1 + k] / npath;
						//if(coup_pay_df.size()>1) fixedlegnonfixedprice=fixedlegnonfixedprice+tempprice*_couponleg_notional[end_indx_adj+1+k]*couponlegtau[end_indx_adj+1+k]/npath;			
						tempprice = 0.0;
					}
				}
			}
		}
		else
		{
			for (i = 0; i < npath; i++)
			{
				fixedrateincount = 0.0;
				ratecount[end_indx_adj] = 0.0;
				for (l = accrualstarti; l < i_coup_calc_shifted_end[0]; l++)
				{
					if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj], _ralow_bdry[end_indx_adj], _rahigh_bdryin_flag[end_indx_adj], _ralow_bdryin_flag[end_indx_adj])) fixedrateincount = fixedrateincount + 1.0;
					ratecount[end_indx_adj] = ratecount[end_indx_adj] + _couponlegindex1_mult[end_indx_adj] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj] * index2_r[i][l] + _couponleg_spread[end_indx_adj] + adj;
				}
				tempprice = tempprice + min(max(ratecount[end_indx_adj] / double(num_fixing[end_indx_adj]), _floorrates[end_indx_adj]), _caprates[end_indx_adj]) * fixedrateincount / double(num_fixing[end_indx_adj]);
			}
			if (coup_pay_df.size() > 0) fixedlegnonfixedprice = tempprice * _couponleg_notional[end_indx_adj] * coup_pay_df[0] * couponlegtau[end_indx_adj] / npath;
			//if(coup_pay_df.size()>0) fixedlegnonfixedprice=tempprice*_couponleg_notional[end_indx_adj]*couponlegtau[end_indx_adj]/npath;						
			if (_couponleg_calc_startdate.size() > 1 && i_coup_calc_shifted_end.size() > 1)
			{
				if (today < _couponleg_calc_shiftedstartdate[1])
				{
					for (k = 0; k<int(i_coup_calc_shifted_end.size()) - 1; k++)
					{
						//////////////////////////////////////////////									
						for (i = 0; i < npath; i++)
						{
							fixedrateincount = 0.0;
							ratecount[end_indx_adj + 1 + k] = 0.0;
							for (l = i_coup_calc_shifted_end[0 + k]; l < i_coup_calc_shifted_end[1 + k]; l++)
							{
								if (rangeincheck(_rngchkindex_mult[k] * rngchkindex_r[i][l] + _rngchkindex2_mult[k] * rngchkindex2_r[i][l], _rahigh_bdry[end_indx_adj + 1], _ralow_bdry[end_indx_adj + 1], _rahigh_bdryin_flag[end_indx_adj + 1], _ralow_bdryin_flag[end_indx_adj + 1])) fixedrateincount = fixedrateincount + 1.0;
								ratecount[end_indx_adj + 1 + k] = ratecount[end_indx_adj + 1 + k] + _couponlegindex1_mult[end_indx_adj + 1 + k] * index1_r[i][l] + _couponlegindex2_mult[end_indx_adj + 1 + k] * index2_r[i][l] + _couponleg_spread[end_indx_adj + 1 + k] + adj;
							}
							tempprice = tempprice + min(max(ratecount[end_indx_adj + 1 + k] / double(num_fixing[end_indx_adj + 1 + k]), _floorrates[end_indx_adj + 1 + k]), _caprates[end_indx_adj + 1 + k]) * fixedrateincount / double(num_fixing[end_indx_adj + 1 + k]);
						}
						if (coup_pay_df.size() > 1) fixedlegnonfixedprice = fixedlegnonfixedprice + tempprice * _couponleg_notional[end_indx_adj + 1 + k] * coup_pay_df[1 + k] * couponlegtau[end_indx_adj + 1 + k] / npath;
						//if(coup_pay_df.size()>1) fixedlegnonfixedprice=fixedlegnonfixedprice+tempprice*_couponleg_notional[end_indx_adj+1+k]*couponlegtau[end_indx_adj+1+k]/npath;			
//////////////////////////////////////////////////////////									
						tempprice = 0.0;
					}
				}
			}
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////									
	/*	if(i_coup_calc_shifted_end.size()>0)
		{
			if(num_calldate>0 && callstarti>0)
			{
				if(today>callnotice_date[callstarti-1])
				{
					for(i=0;i<npath;i++)
					{
						uptocall_coup[0][i]=fixedlegnonfixedprice;
					}
				}
			}
			else if(num_calldate>0 && callstarti<=0)
			{
				//if(today>callnotice_date[callstarti-1])
				{
					for(i=0;i<npath;i++)
					{
						uptocall_coup[0][i]=fixedlegnonfixedprice;
					}
				}
			}
		}
		for(j=callstarti;j<=num_calldate-2;j++)
		{
			for(i=0;i<npath;i++)
			{
				//uptocall_coup[j-callstarti+1][j]=uptocall_coup[j-callstarti][j]+call_coup[j-callstarti][i];
				uptocall_coup[j-callstarti+1][i]=uptocall_coup[j-callstarti][i]+call_coup[j-callstarti][i];
			}
		}
		for(j=callstarti;j<=num_calldate-1;j++)
		{
			for(i=0;i<npath;i++)
			{
				//hX[j-callstarti][i]=hX[j-callstarti][i]+uptocall_coup[j-callstarti][j]*coup_pay_df[aftercall_coup_starti[j-callstarti]-1-pay_indx_adj];
				hX[j-callstarti][i]=hX[j-callstarti][i]+uptocall_coup[j-callstarti][i]*coup_pay_df[aftercall_coup_starti[j-callstarti]-1-pay_indx_adj];
			}
		}*/
		/////////////////////////////////////////////////////////////////////////////////////////////////////////									
	fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
	//fixedlegnonfixedprice=fixedlegnonfixedprice/ondf/settlement_date_df*coup_pay_df[num_remained_couponleg_cf-1-pay_indx_adj];								


	int
		nx = num_factor,
		deg = 2;

	double sum = 0.0;

	LSRegression
	(
		X,
		hX,
		_callability_flag,
		num_calldate - callstarti,
		Npath,
		deg,
		nx,
		sum
	);

	////////////////CMS Spread Range Accrual Swap   End/////////////////////////////////////////////////////////////////////////////									

	ofstream fout(argv[2]);

	vector<double> value(5, 0.0);

	value[1] = sum;

	value[3]
		#NAME ?
		+floatinglegnonfixedprice;

	value[4]
		= -(
			fixedlegfixedprice
			#NAME ?
			#NAME ?
			);

	if (aftercall_coup_price.size() > 0)
	{
		value[4]
			= value[4]
			#NAME ?
			/ ondf
			/ settlement_date_df;
	}

	value[2]
		= value[3]
		+ value[4];

	value[0]
		= value[1]
		+ value[2];

	fout << "Price_Structured" << endl;
	fout << setprecision(15) << value[0] << endl;

	fout << "Price_Bermudan" << endl;
	fout << setprecision(15) << value[1] << endl;

	fout << "Price_Swap" << endl;
	fout << setprecision(15) << value[2] << endl;

	fout << "Price_FundingLeg" << endl;
	fout << setprecision(15) << value[3] << endl;

	fout << "Price_CouponLeg" << endl;
	fout << setprecision(15) << value[4] << endl;

	delete[]label, couponleg_holidays, fundingleg_holidays, fixing_holidays, call_holidays;

}
