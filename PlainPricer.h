#pragma once

#include <vector>
#include "Interpolation.h"
#include "InterestRateCurve.h"
#include "ConstantNumbers.h"
#include "MyUtility.h"

using namespace std;

void PlainSwapPrice
(
    CCurrency crcy
    , CDate today
    , CDate settlement_date
    , double notional
    , CInterpolation zc
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void MyPlainSwapPrice
(
    CCurrency crcy
    , CDate today
    , CDate settlement_date
    , double notional
    , CInterpolation zc
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void PlainCRSPrice
(
    double fxspot
    , CCurrency crcy1
    , CCurrency crcy2
    , CDate today
    , CDate start_date
    , CDate settlement_date
    , CInterpolation zc1
    , CInterpolation zc2
    , CCurrency coupcrcy
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double &floatinglegprice
    , double &fixedlegprice
    , double &atmswaprate
    , double &swapprice
);

void PlainCRS2FixedPrice
(
    double fxspot
    , CCurrency crcy1
    , CCurrency crcy2
    , CDate today
    , CDate start_date
    , CDate settlement_date
    , CInterpolation zc1
    , CInterpolation zc2
    , CCurrency coupcrcy
    , vector<CDate> couponleg1_calc_startdate
    , vector<CDate> couponleg1_calc_enddate
    , vector<CDate> couponleg1_paydate
    , string couponleg1_dcb
    , vector<double> couponleg1_notional
    , vector<double> couponleg1_couponrate
    , vector<CDate> couponleg2_calc_startdate
    , vector<CDate> couponleg2_calc_enddate
    , vector<CDate> couponleg2_paydate
    , string couponleg2_dcb
    , vector<double> couponleg2_notional
    , vector<double> couponleg2_mult
    , vector<double> couponleg2_couponrate
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& couponleg2price
    , double& couponleg1price
    , double& atmswaprate
    , double& swapprice
);

void PlainCRSPrice_org
(
    double fxspot
    , CCurrency crcy
    , CCurrency crcy1
    , CDate today
    , CDate start_date
    , CDate settlement_date
    , CInterpolation zc
    , CInterpolation cszc
    , CInterpolation basiszc
    , CCurrency coupcrcy
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void PlainZeroCouponSwapPrice
(
    CCurrency crcy
    , CDate today
    , CDate settlement_date
    , double notional
    , CInterpolation zc
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void PlainZeroCouponCompoundSwapPrice
(
    CCurrency crcy
    , CDate today
    , CDate settlement_date
    , double notional
    , CInterpolation zc
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void MyPlainZeroCouponCompoundSwapPrice
(
    CCurrency crcy
    , CDate today
    , CDate settlement_date
    , double notional
    , CInterpolation zc
    , vector<CDate> couponleg_calc_startdate
    , vector<CDate> couponleg_calc_enddate
    , vector<CDate> couponleg_paydate
    , string couponleg_dcb
    , vector<double> couponleg_notional
    , vector<double> couponleg_couponrate
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , vector<CDate> fixing_date
    , vector<CDate> fixingindex_maturity
    , string fixingindex_dcb
    , vector<CDate> fixinghistory_date
    , vector<double> fixinghistory_rate
    , bool estflag
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
);

void PlainSwapPrice(CCurrency crcy, CDate today, CDate settlement_date, double notional, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double couponleg_couponrate, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, double fixinghistory_rate, bool estflag, double& floatinglegprice, double& fixedlegprice, double& atmswaprate, double& swapprice);

void bondpricezc(CDate today, CDate spotdate, CInterpolation zc, double coupon, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double& bondd, double& bondc);

double findbondytmzc(CDate today, CDate spotdate, CInterpolation zc, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double bondp);

void bondpriceytm(CDate today, CDate settledate, double coupon, double ytm, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double& bondd, double& bondc);

double findbondytm(CDate today, CDate spotdate, CDate settledate, double coupon, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double bondp);

void CapFloorForwardVolMatrix(CCurrency crcy, CDate today, CDate* indexholiday, int num_indexholiday, CInterpolation zc, int num_strike, double centervol, double volinterval, vector<string> tenor, vector<string> freq, vector<double> capatmvol, vector<double>& capstrike, vector<CDate>& cpltmtrty, vector<vector<double>>& forwardvol);

void CapFloorForwardVolMatrix(CCurrency crcy, CDate today, CDate* indexholiday, int num_indexholiday, CInterpolation zc, vector<string> tenor, vector<string> freq, vector<vector<double>> capvolskew, vector<double> capstrike, vector<CDate>& cpltmtrty, vector<vector<double>>& forwardvol);

void PlainCapFloorDigitalPrice(CCurrency crcy, CDate today, CDate settlement_date, double notional, double strike, double inputvol, string callputflag, bool firstcashflowincludeflag, double firstfixingrate, CInterpolation zc, CInterpolation vol, vector<double> fundingleg_notional, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, vector<CDate> fundinglegindex_calc_startdate, string fundingleg_dcb, vector<double> fundingleg_mult, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, double& atmrate, double& flatvol, double& price);

void PlainCapFloorDigitalPrice2(CCurrency crcy, CDate today, CDate settlement_date, double notional, double strike, double inputvol, string callputflag, bool firstcashflowincludeflag, double firstfixingrate, CInterpolation zc, CInterpolation vol, vector<double> fundingleg_notional, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, vector<CDate> fundinglegindex_calc_startdate, string fundingleg_dcb, vector<double> fundingleg_mult, vector<CDate> fixing_date, vector<CDate> fixingindex_calc_startdate, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, double& atmrate, double& flatvol, double& price);

void CapPrice_Bl(CCurrency crcy, CDate today, CDate settlement_date, CInterpolation zc, CInterpolation vol, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<CDate> fixing_date, vector<CDate> fixingindex_calc_startdate, vector<CDate> fixingindex_maturity, string fixingindex_dcb, double atmrate, double& price);

void ParSwapRate(CDate today, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double& atmswaprate);

void FindATMVol(int starti, int endi, int num_fundingleg_cf, double atmswaprate, double price, CDate today, double strike, string callputflag, vector <double> fundingleg_notional, vector<CDate> fixing_date, vector<CDate> fundingleg_paydate, vector<double> fundingleg_df, vector<double> fwd1, vector<double> fwd2, vector<double> tau, vector<double> T, vector<int> sgn_1stcfincludeflag, double& flatvol);

double Bl(double K, double F, double v, double w);

double Caplet(double fwd, double k, double df, double tau, double t, double vol);

double Floorlet(double fwd, double k, double df, double tau, double t, double vol);

double DigitalCaplet(double fwd, double k, double df, double tau, double t, double vol);

double DigitalFloorlet(double fwd, double k, double df, double tau, double t, double vol);

double PlainOptionlet(string callputflag, double fwd, double k, double df, double tau, double t, double vol);

double DigitalSpreadCall(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k);

double DigitalSpreadPut(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k);

double SpreadCall(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k);

double SpreadPut(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k);

double Capletsum(double fwdvol0, double fwdvol1, double previouscap, double k, vector<double> df, vector<double> capletmaturity, vector<double> capletperiod, vector<double> fra);

double Findfwdvol(double fwdvol0, double fullcap, double previouscap, double k, vector<double> df, vector<double> capletmaturity, vector<double> capletperiod, vector<double> fra);

void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<CDate>> swapmaturity_date, vector<vector<double>> mktswaptionvol, vector<CDate> targetswapmaturity_date, CDate option_maturity, CDate swap_maturity, double& swaptionvol);

void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<double>> mktswaptionvol, vector<CDate> targetswapmaturity_date, CDate option_maturity, CDate swap_maturity, double& swaptionvol);

void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<double>> mktswaptionvol, vector<double> targetswapmaturity_t, CDate option_maturity, double swap_t, double& swaptionvol);

void PlainSwaptionPrice(CCurrency crcy, CDate today, CDate option_maturity, CDate settlement_date, double notional, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<double> fundingleg_notional, vector<double> fundingleg_mult, vector<double> fundingleg_spread, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, bool estflag, string payout_type, bool exercise_flag, string delivery_type, double swaptionvol, double& atmswaprate, double& optionprice);

void SwaptionPrice_Bl(CDate today, CDate option_maturity, CDate settlement_date, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double swaptionvol, double atmswaprate, double& optionprice);

void PlainSwapFundinglegFixedPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , vector<double> fixinghistory_rate
    , int num_fundingleg_cf
    , int fundingleg_cf_starti
    , vector<CDate> fixing_date
    , vector<CDate> fundingleg_calc_startdate
    , vector<CDate> fundingleg_calc_enddate
    , vector<CDate> fundingleg_paydate
    , string fundingleg_dcb
    , vector<double> fundingleg_notional
    , vector<double> fundingleg_mult
    , vector<double> fundingleg_spread
    , double& floatingleg_fixed_price
);

void PlainSwapFundinglegNonFixedPrice
(
    CDate today,
    double ondf,
    double settlement_date_df,
    CInterpolation zc,
    int num_remained_fundingleg_cf,
    vector<double> _fundingleg_notional,
    vector<double> _fundingleg_mult,
    vector<double> _fundingleg_spread,
    CDate fundingleg_calc_startdate0,
    double fundinglegcalcstartT0,
    vector<CDate> _fundingleg_paydate,
    vector<double> fundinglegtau,
    vector<double> fundinglegT,
    double& floatinglegnonfixedprice
);

void PlainSwapFundinglegNonFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_fundingleg_cf, vector<double> _fundingleg_notional, vector<double> _fundingleg_mult, vector<double> _fundingleg_spread, CDate fundingleg_calc_startdate0, double fundinglegcalcstartT0, vector<CDate> _fundingleg_calc_startdate, vector<CDate> _fundingleg_calc_enddate, vector<CDate> _fundingleg_paydate, vector<double> fundinglegtau, vector<double> fundinglegT, double& floatinglegnonfixedprice);

void PlainSwapCouponlegFixedPrice
(
    CDate today,
    double ondf,
    double settlement_date_df,
    CInterpolation zc,
    int num_remained_couponleg_cf,
    vector<bool> _ra_flag,
    vector<double> _couponleg_notional,
    vector<double> _couponleg_couponrate,
    vector<CDate> _couponleg_paydate,
    vector<double> couponlegtau,
    vector<double> couponlegT,
    double& fixedlegfixedprice
);

void PlainSwapZeroCouponlegFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void CMSSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void CMSSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void FixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void FixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void FixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void QuantoCMSSteepnerFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void CMSSlopeFixedAccrualPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> avg_fixeddenom
    , vector<vector<bool>> avg_fixedflag
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _caprates
    , vector<double> _floorrates
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrate1count
    , vector<double>& fixedrate2count

);

void CMSSlopeFixedAccrualPrice_XXX
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> avg_fixeddenom
    , vector<vector<bool>> avg_fixedflag
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _caprates
    , vector<double> _floorrates
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrate1count
    , vector<double>& fixedrate2count
    , double& fixedlegaccrualfixedprice
);

void CMSSlopeFixedAccrualPriceForFixedfloater(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count);

void CMSSlopeFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count);

void CMSSlopeFixedPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count);

void CDCMSSpreadDualFixedAccrualPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> fixedrate0
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _caprates
    , vector<double> _floorrates
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , double& fixedlegaccrualfixedprice
);

void CDCMSSpreadDualFixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate0, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void CDCMSSpreadDualLeverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate0, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _couponleg_ra_mult, vector<double> _couponleg_ra_sp, vector<double> _couponleg_ra_cap, vector<double> _couponleg_ra_floor, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void QuantoCMSSpreadDualLeverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed,
    vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponleg_notional,
    vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult,
    vector<double> _rahigh_bdry1, vector<bool> _rahigh_bdryin_flag1, vector<double> _ralow_bdry1, vector<bool> _ralow_bdryin_flag1,
    vector<double> _rahigh_bdry2, vector<bool> _rahigh_bdryin_flag2, vector<double> _ralow_bdry2, vector<bool> _ralow_bdryin_flag2,
    vector<double> _couponleg_ra_mult, vector<double> _couponleg_ra_sp, vector<double> _couponleg_ra_cap, vector<double> _couponleg_ra_floor, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void QuantoCMSSpreadDualLeverageFloatingAccrualPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> fixedrate3
    , vector<vector<double>> fixedrate4
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _couponlegindex3_mult
    , vector<double> _couponlegindex4_mult
    , vector<double> _rahigh_bdry1
    , vector<bool> _rahigh_bdryin_flag1
    , vector<double> _ralow_bdry1
    , vector<bool> _ralow_bdryin_flag1
    , vector<double> _rahigh_bdry2
    , vector<bool> _rahigh_bdryin_flag2
    , vector<double> _ralow_bdry2
    , vector<bool> _ralow_bdryin_flag2
    , vector<double> _couponleg_ra_mult
    , vector<double> _couponleg_ra_sp
    , vector<double> _couponleg_ra_cap
    , vector<double> _couponleg_ra_floor
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<vector<double>> f_rate
    , double& fixedlegaccrualfixedprice
);

void CMSSpreadNonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>>index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void CMSSpreadNonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag
    , vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag
    , vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void CMSSpreadNonFixedAccrualPrice_Repli(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void SwapRatewithVoladj(CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_estm, vector<int> num_index1_coup_cf,
    vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol, double input_corr, vector<double> index1_t,
    vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T,
    vector<vector<vector<double>>> index1_payT, vector<double>, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& swaprate);

void PlainCMSSpreadRangeAccrualSwap(CCurrency crcy, CDate today, CDate settlement_date, double notional, CInterpolation zc, int num_fundingleg_cf, int fundingleg_cf_starti, vector<CDate> fixing_date, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<double> fundingleg_notional, vector<double> fundingleg_mult, vector<double> fundingleg_spread, int couponleg_cf_starti, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> rahigh_bdryoh, vector<int> rahigh_bdryoh_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<double> ralow_bdryoh, vector<int> ralow_bdryoh_flag
    , vector<double> fixinghistory_rate, double& floatinglegprice, double& fixedlegprice, double& atmswaprate, double& swapprice);

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, int d, vector<int> num_fixing, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT
    , vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void PlainSwapFundinglegNonFixedPricexy(CDate callnotice_date, double callnotice_t, double ondf, CInterpolation zc, int couponleg_cf_starti, int num_remained_fundingleg_cf, vector<double> _fundingleg_notional, vector<double> _fundingleg_mult, vector<double> _fundingleg_spread, CDate fundingleg_calc_startdate0, double fundinglegcalcstartT0, vector<CDate> _fundingleg_paydate, vector<double> fundinglegtau, vector<double> fundinglegT, double& floatinglegnonfixedprice);

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT
    , vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void CMSSpreadNonFixedAccrualPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedZeroPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedZeroPricexy2(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice, double& fixedlegpartialnonfixedprice);

void NonFixedZeroPricexy3(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedAccrualPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk
    , vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice);

void NonFixedAccrualZeroPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry
    , vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double dP, double& fixedlegnonfixedprice);

void QuantoFloatingRateCMSFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<double> _couponfloatingfixinghistory_rate, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedlegaccrualfixedprice, vector<double>& fixedlegaccrualfixednumber);

void BasisFixedAccrualPrice
(
    CDate today,
    double ondf,
    double settlement_date_df,
    CInterpolation zc,
    int num_remained_couponleg_cf,
    vector<int> num_fixing,
    vector<int> num_fixed,
    vector<vector<double>> fixedrate1,
    vector<vector<double>> fixedrate2,
    vector<bool> _ra_flag,
    vector<double> _couponleg_notional,
    vector<double> _couponleg_couponrate,
    vector<double> _couponleg_spread,
    vector<double> _couponlegindex1_mult,
    vector<double> _couponlegindex2_mult,
    vector<double> _rahigh_bdry,
    vector<bool> _rahigh_bdryin_flag,
    vector<double> _ralow_bdry,
    vector<bool> _ralow_bdryin_flag,
    vector<CDate> _couponleg_paydate,
    vector<double> couponlegtau,
    vector<double> couponlegT,
    double& fixedlegaccrualfixedprice
);

void LIBORSpreadDualFixedAccrualPrice
(
    CDate today,
    double ondf,
    double settlement_date_df,
    CInterpolation zc,
    int num_remained_couponleg_cf,
    vector<int> num_fixing,
    vector<int> num_fixed,
    vector<vector<double>> fixedrate1,
    vector<vector<double>> fixedrate2,
    vector<bool> _ra_flag,
    vector<double> _couponleg_notional,
    vector<double> _couponleg_couponrate,
    vector<double> _couponleg_spread,
    vector<double> _couponlegindex1_mult,
    vector<double> _couponlegindex2_mult,
    vector<double> _caprates,
    vector<double> _floorrates,
    vector<double> _rahigh_bdry,
    vector<bool> _rahigh_bdryin_flag,
    vector<double> _ralow_bdry,
    vector<bool> _ralow_bdryin_flag,
    vector<CDate> _couponleg_paydate,
    vector<double> couponlegtau,
    vector<double> couponlegT,
    double& fixedlegaccrualfixedprice
);

void LIBORSpreadDualFixedAccrualPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void LIBORSpreadDualFixedAccrualPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void BasisTripleFixedAccrualPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> fixedrate3
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _couponlegindex3_mult
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<double> _rahigh_bdry2
    , vector<bool> _rahigh_bdry2in_flag
    , vector<double> _ralow_bdry2
    , vector<bool> _ralow_bdry2in_flag
    , vector<double> _rahigh_bdry3
    , vector<bool> _rahigh_bdry3in_flag
    , vector<double> _ralow_bdry3
    , vector<bool> _ralow_bdry3in_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , double& fixedlegaccrualfixedprice
);

void BasisDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void BasisDualSpreadFixedAccrualPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> fixedrate3
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _couponlegindex3_mult
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<double> _rahigh_bdry2
    , vector<bool> _rahigh_bdry2in_flag
    , vector<double> _ralow_bdry2
    , vector<bool> _ralow_bdry2in_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , double& fixedlegaccrualfixedprice
);

void BasisTripleFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _rahigh_bdry1, vector<bool> _rahigh_bdryin1_flag, vector<double> _ralow_bdry1, vector<bool> _ralow_bdryin1_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void BasisDualSpreadTripleFixedAccrualZCPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> fixedrate3
    , vector<vector<double>> fixedrate4
    , vector<vector<double>> fixedrate5
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate
    , vector<double> _couponleg_spread
    , double changeCondIdx
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _couponlegindex3_mult
    , vector<double> _couponlegindex4_mult
    , vector<double> _couponlegindex5_mult
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<double> _rahigh_bdry2
    , vector<bool> _rahigh_bdry2in_flag
    , vector<double> _ralow_bdry2
    , vector<bool> _ralow_bdry2in_flag
    , vector<double> _rahigh_bdry3
    , vector<bool> _rahigh_bdry3in_flag
    , vector<double> _ralow_bdry3
    , vector<bool> _ralow_bdry3in_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , double& fixedlegaccrualfixedprice
);

void BasisDualSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _rahigh_bdry2, vector<bool> _rahigh_bdry2in_flag, vector<double> _ralow_bdry2, vector<bool> _ralow_bdry2in_flag, vector<double> _caps, vector<double> _multis, vector<double> _sprds, vector<double> _floors, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void SpreadDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate21, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate31, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex21_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex31_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void SpreadDualFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void CDFXDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void StartEndVolFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void StartEndVolFixedPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<bool> arr_flag, double& fixedlegaccrualfixedprice);

void QuantoDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void QuantoDualFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> couponlegtau, double& fixedlegaccrualfixedprice);

void GeneralSwapCouponlegFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation indexzc, int num_remained_couponleg_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice);

void GeneralSwapCouponlegNonFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation indexzc, int num_remained_couponleg_cf, vector<int> num_index_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<vector<double>> index_tau, vector<CDate> index_startdate, vector<double> index_startT, vector<vector<CDate>> index_paydate, vector<vector<double>> index_payT, vector<double> couponindex_fixingT, vector<double> index_T, vector<double> index_vol, vector<vector<double>>  index_d, vector<double> index_t, double& fixedlegnonfixedprice);

void ComboLeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate1count, vector<double>& fixedrate2count, vector<double>& fixedrate3count, vector<double>& fixedrate4count);

void ComboLeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate1count, vector<double>& fixedrate2count, vector<double>& fixedrate3count, vector<double>& fixedrate4count, vector<double> _rahigh_bdry, vector<double> _ralow_bdry, vector<bool> _rahigh_bdryin_flag, vector<bool> _ralow_bdryin_flag, int& tmpfixedincount);

void LeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate2count, vector<double>& fixedrate4count);

void QuantoLeveragedAverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount);

void QuantoLeveragedAverageFixedAccrualPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice);

void QuantoLeveragedAverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg);

void QuantoLeveragedAverageFixedAccrualZeroPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_fixed_couponleg_cf
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<double> couponleg_couponrate
    , vector<bool> ra_flag
    , vector<double> couponleg_notional
    , vector<double> couponleg_spread
    , vector<double> couponlegindex1_mult
    , vector<double> couponlegindex2_mult
    , vector<double> caprates
    , vector<double> floorrates
    , vector<double> rahigh_bdry
    , vector<bool> rahigh_bdryin_flag
    , vector<double> ralow_bdry
    , vector<bool> ralow_bdryin_flag
    , double& fixedlegfixedprice
);

void QLevAveFixedAccZeroPrice_XXX
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_fixed_couponleg_cf
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<double> couponleg_couponrate
    , vector<bool> ra_flag
    , vector<double> couponleg_notional
    , vector<double> couponleg_spread
    , vector<double> couponlegindex1_mult
    , vector<double> couponlegindex2_mult
    , vector<double> caprates
    , vector<double> floorrates
    , vector<double> rahigh_bdry
    , vector<bool> rahigh_bdryin_flag
    , vector<double> ralow_bdry
    , vector<bool> ralow_bdryin_flag
    , double& fixedlegfixedprice
);

void QLevAveFixedAccZeroPrice_XXX
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_fixed_couponleg_cf
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<double> couponleg_couponrate
    , vector<bool> ra_flag
    , vector<double> couponleg_notional
    , vector<double> couponleg_spread
    , vector<double> couponlegindex1_mult
    , vector<double> couponlegindex2_mult
    , vector<double> caprates
    , vector<double> floorrates
    , vector<double> rahigh_bdry
    , vector<bool> rahigh_bdryin_flag
    , vector<double> ralow_bdry
    , vector<bool> ralow_bdryin_flag
    , double& fixedlegfixedprice
);

void QuantoLeveragedAverageFixedAccrualPrice2
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<vector<double>> rngchkfixedrate2
    , vector<double> _rngchkindex_mult
    , vector<double> _rngchkindex2_mult
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _caprates
    , vector<double> _floorrates
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrateincount
    , vector<double>& fixedavg

);

void QLevAveFixedAccrualPrice2_XXX
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_remained_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<vector<double>> rngchkfixedrate2
    , vector<double> _rngchkindex_mult
    , vector<double> _rngchkindex2_mult
    , vector<bool> _ra_flag
    , vector<double> _couponleg_notional
    , vector<double> _couponleg_spread
    , vector<double> _couponlegindex1_mult
    , vector<double> _couponlegindex2_mult
    , vector<double> _caprates
    , vector<double> _floorrates
    , vector<double> _rahigh_bdry
    , vector<bool> _rahigh_bdryin_flag
    , vector<double> _ralow_bdry
    , vector<bool> _ralow_bdryin_flag
    , vector<CDate> _couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrateincount
    , vector<double>& fixedavg
);

void QuantoAvgSpreadLiborRAFixedZeroCouponPrice
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_past_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<bool> ra_flag
    , vector<double> couponleg_notional
    , vector<double> couponleg_spread
    , vector<double> couponlegindex1_mult
    , vector<double> couponlegindex2_mult
    , vector<double> caprates
    , vector<double> floorrates
    , vector<double> rahigh_bdry
    , vector<bool> rahigh_bdryin_flag
    , vector<double> ralow_bdry
    , vector<bool> ralow_bdryin_flag
    , vector<CDate> couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrateincount
    , vector<double>& fixedavg
    , double& fixedlegaccrualfixedprice

);

void QAvgSprLibRAFixedZC_XXX
(
    CDate today
    , double ondf
    , double settlement_date_df
    , CInterpolation zc
    , int num_past_couponleg_cf
    , vector<int> num_fixing
    , vector<int> num_fixed
    , vector<vector<double>> fixedrate1
    , vector<vector<double>> fixedrate2
    , vector<vector<double>> rngchkfixedrate
    , vector<bool> ra_flag
    , vector<double> couponleg_notional
    , vector<double> couponleg_spread
    , vector<double> couponlegindex1_mult
    , vector<double> couponlegindex2_mult
    , vector<double> caprates
    , vector<double> floorrates
    , vector<double> rahigh_bdry
    , vector<bool> rahigh_bdryin_flag
    , vector<double> ralow_bdry
    , vector<bool> ralow_bdryin_flag
    , vector<CDate> couponleg_paydate
    , vector<double> couponlegtau
    , vector<double> couponlegT
    , vector<double>& fixedrateincount
    , vector<double>& fixedavg
    , double& fixedlegaccrualfixedprice
);

void QuantoAvgSpreadLiborRAFixedZeroCouponPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_past_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<vector<double>> rngchkfixedrate2,
    vector<double> rngchkindex_mult, vector<double> rngchkindex2_mult, vector<bool> ra_flag, vector<double> couponleg_notional, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> caprates, vector<double> floorrates, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg, double& fixedlegaccrualfixedprice);

void QuantoAvgSpreadLiborRAFixedPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_past_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<double> couponleg_couponrate, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<vector<double>> rngchkfixedrate2,
    vector<double> rngchkindex_mult, vector<double> rngchkindex2_mult, vector<bool> ra_flag, vector<double> couponleg_notional, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> caprates, vector<double> floorrates, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg, double& fixedlegaccrualfixedprice);

void PlainIROption_HW1F(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y, vector<double>& dydparam);

void mrqcof_HW1F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&));

void mrqmin_HW1F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&), double& alamda, bool& singularflag);

void rstar_HW1F(vector<double>& c, vector<double>& AtT, vector<double>& BtT, double& rstar);

void gaussj(vector<vector<double>>& a, vector<vector<double>>& b, bool& singularflag);

void mrqmin_HW1Fsig(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param1, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs_sig(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&), double& alamda, bool& singularflag);

void mrqcof_HW1Fsig(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param1, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs_sig(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&));

void PlainIROption_HW1Fsig(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, double param1, double& param0, double& y, double& dydparam0);

void mrqmin_HW1Fa(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param0, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs_a(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&), double& alamda, bool& singularflag);

void mrqcof_HW1Fa(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param0, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs_a(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&));

void PlainIROption_HW1Fa(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, double param0, double& param1, double& y, double& dydparam1);

void PlainIROption_HW2F
(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y, vector<double>& dydparam);

void PlainIROption_HW2F(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y);

void mrqcof_HW2F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&));

void mrqmin_HW2F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&), double& alamda, bool& singularflag);

void ybar_HW2F(vector<vector<double>>& lambij, vector<double>& BbtT, vector<double>& ybar);

double PlainIROption_HW2F_SE(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, vector<double>& param);

double amotsa_HW2F(int& idum, double& tt, vector<vector<double>>& p, vector<double>& yy, vector<double>& psum, vector<double>& pb, double& yb, vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, double funcs(vector<string>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&), const int ihi, double& yhi, const double fac);

void amebsa_HW2F(int& idum, double& tt, vector<vector<double>>& p, vector<double>& yy, vector<double>& pb, double& yb, const double ftol, vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, double funcs(vector<string>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&), int& iter, const double temptr);









