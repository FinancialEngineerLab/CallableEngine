#pragma once
#include "Date.h"
#include "ConstantNumbers.h"

using namespace std;

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schedule
);

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
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

void depositeschedule
(
	CDate today,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string adjflag,
	string payin,
	int spotlag,
	CDate &calstartdate,
	CDate& calenddate,
	CDate& paydate
);

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
);

void fixedlegcashflowschedule
(
	CDate today,
	CDate maturitydate,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<vector<CDate>> &schedule
);

void fixedlegcashflowschedule
(
	CDate today,
	string tenor,
	CDate* holidays,
	int numholiday,
	string stub,
	string direction,
	string conv,
	string freq,
	string adjflag,
	vector<vector<CDate>>& schedule
);

void ShiftBusDate
(
	CDate date,
	CDate *holidays,
	int numholiday,
	int numshift,
	CDate &newdate
);

void ConvReflDate
(
	CDate date,
	CDate* holidays,
	int numholiday,
	string conv,
	CDate &newdate
);

bool is_holiday
(
	CDate date,
	CDate *holidays,
	int numholiday
);

int findnumschedule
(
	CDate startdate,
	CDate maturitydate,
	string stub,
	string direction,
	string freq
);

int findnumschedule
(
	CDate startdate,
	string tenor,
	string stub,
	string direction,
	string freq
);

int findnumschedule
(
	CDate startdate,
	CDate maturitydate,
	string stub,
	string direction,
	int frq
);

int findnumschedule
(
	CDate startdate,
	string tenor,
	string stub,
	string direction,
	int frq
);

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
	int num_schedule,
	vector<CDate> &schdule
);

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
	vector<CDate>& schdule
);

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
);

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
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<CDate>& schdule
);

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
	vector<vector<CDate>>& schdule
);

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
	vector<vector<CDate>>& schdule
);

void findmaturity
(
	CDate startdate,
	int tenor,
	CDate *holidays,
	int num_holiday,
	string conv,
	CDate &maturity
);

void findmaturity
(
	CDate startdate,
	string tenor,
	CDate* holidays,
	int num_holiday,
	string conv,
	CDate& maturity
);

CDate* Holidays(string holidaysname);

int NumHolidays(string holidaysname);

void avgfixedschedule
(
	string avg_freq,
	double avg_interval,
	string avg_weekday,
	string avg_sched_from,
	string avg_fixing_setin,
	string avg_fixing_lag,
	string avg_fixing_conv,
	vector<CDate> _coupleleg_calc_startdate,
	vector<CDate> couponleg_calc_shiftedstartdate,
	vector<CDate> couponleg_calc_shiftedenddate,
	vector<CDate> _couponleg_paydate,
	vector<CDate*> index1_holidays,
	vector<int> num_index1holidays,
	vector<bool> _ra_flag,
	vector<vector<CDate>> couponindex_fixingdate,
	vector<CDate> date,
	vector<vector<CDate>> &couponindex_avg_fixeddate,
	vector<vector<bool>> &avg_fixedflag,
	vector<bool> &avg_flag,
	vector<int> &avg_fixeddenom
);

void StartEndSchedule
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
	vector<vector<CDate>> couponindex_fixingdate,
	vector<CDate> date,
	vector<vector<CDate>> &adv_date,
	vector<vector<CDate>>& arr_date,
	vector<vector<bool>>& adv_fixedflag,
	vector<vector<bool>>& arr_fixedflag,
	vector<bool>& adv_flag,
	vector<bool>& arr_flag
);

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
	vector<CDate> couponleg_calc_enddate,
	vector<CDate> _couponleg_paydate,
	vector<CDate*> index1_holidays,
	vector<int> num_index1holidays,
	vector<bool> _ra_flag,
	vector<vector<CDate>> couponindex_fixingdate,
	vector<vector<CDate>> &couponindex_avg_fixeddate,
	vector<vector<bool>> &avg_fixedflag,
	vector<int> &avg_fixeddenom
);
*/

CDate find_nth_day_of_yyyymm
(
	int yyyy,
	int mm,
	int nth,
	int day,
	CDate *holidays,
	int numholiday
);

CDate find_nth_day_of_MMMYY
(
	string MMMYY,
	int nth,
	int day,
	CDate* holidays,
	int numholiday
);

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
	CDate &calstartdate,
	CDate &calenddate,
	CDate &lasttradedate
);

void fxswapschedule
(
	CDate today,
	string swappointtenor,
	CDate *fxswapholiday,
	int numfxswapholiday,
	CDate *fxswapstartholiday,
	int numfxswapstartholiday,
	string fxswapconv,
	int fxswaplag,
	CDate &calstartdate,
	CDate &calenddate
);

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
	vector<bool>& arr_flag
);


