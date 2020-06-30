#pragma once

#include "Instrument.h"
#include "Interpolation.h"

using namespace std;

void zerocurve
(
	CCurrency crcy,
	CDate today,
	vector<string> tenor,
	vector<string> type,
	vector<string> fixedrateleg_freq,
	vector<string> floatingrateleg_freq,
	vector<double> mktrate,
	vector<CDate>& maturity,
	vector<double>& df,
	vector<double>& zero
);

void zerocurve
(
	CCurrency crcy,
	CDate today,
	vector<string> tenor,
	vector<string> type,
	vector<string> fixedrateleg_freq,
	vector<string> floatingrateleg_freq,
	vector<double> mktrate,

	vector<CDate>& startdate,
	vector<CDate>& maturity,

	vector<double>& df,
	vector<double>& zero
);

void interpolatedamt
(
	double notional,
	double startzero,
	double endzero,
	vector<double> coupon,
	vector<double> interval,
	CDate startdate,
	CDate enddate,
	vector<CDate> interdate,
	int startindx,
	int num_coup,
	double& sum
);

void findzerorate
(
	double notional,
	double amt,
	double startzero,
	vector<double> coupon,
	vector<double> interval,
	CDate startdate,
	CDate enddate,
	vector<CDate> interdate,
	int startindx,
	double& zero
);

void findzdr
(
	CCurrency crcy,
	CDate today,
	vector<string> tenor,
	vector<string> type,
	vector<string> fixedrateleg_freq,
	vector<string> floatingrateleg_freq,
	vector<double> mktrate,
	double dr,
	vector<vector<double>>& dzdr
);

void swappoint2marketrate
(
	bool swappt_mode,
	CCurrency crcy,
	CDate today,
	double fxspot,
	vector<string> swappointtenor,
	vector<double> swap_point,
	CInterpolation zc,
	vector<CDate>& fxswapstartdate,
	vector<CDate>& fxswapmaturity,
	vector<double>& fxrate,
	vector<string>& resulttenor,
	vector<double>& marketquote
);

void basiszerocurve
(
	CCurrency crcy,
	CDate today,
	vector<CDate> startdate,
	vector<CDate> mtrty,
	CInterpolation zc,
	double fxrate,

	CCurrency index1_crcy,
	string index1_tenor,
	string index1_type,

	CCurrency index2_crcy,
	string index2_tenor,
	string index2_type,

	vector<string> spread_tenor,
	vector<string> spread_type,

	vector<string> index1_freq,
	vector<string> index2_freq,

	vector<double> mkt_spread,

	vector<CDate>& basis_startdate,
	vector<CDate>& basis_mtrty,

	vector<double>& basis_df,
	vector<double>& basis_zero,

	vector<CDate>& spread_startdate,
	vector<CDate>& spread_mtrty,

	vector<double>& spread_df,
	vector<double>& spread_zero
);

void basiszero
(
	CDate today,
	CDate spotdate,
	double spotdf,
	vector<CDate> tenor,

	int num_mktspread,

	vector<vector<CDate>> start,
	vector<vector<CDate>> end,

	vector<vector<CDate>> fixingend,

	vector<vector<double>> df,

	CInterpolation zc,

	string indexdcb,

	vector<double> amt,

	vector<double>& x,

	vector<double>& fvec,

	vector<vector<double>>& fjac
);

void findbasiszeromnewt
(
	CDate today,
	CDate spotdate,
	double spotdf,
	vector<CDate> tenor,
	int num_mktspread,
	vector<vector<CDate>> start,
	vector<vector<CDate>> end,
	vector<vector<CDate>> fixingend,
	vector<vector<double>> df,

	CInterpolation zc,

	string dcb,

	vector<double> amt,

	const int ntrial,

	vector<double>& x,

	const double tolx,
	const double tolf
);

void findbasisdzdr
(
	CCurrency crcy,
	CDate today,
	vector<CDate> startdate,
	vector<CDate> mtrty,

	CInterpolation zc,

	double fxrate,

	CCurrency index1_crcy,
	string index1_tenor,
	string index1_type,

	CCurrency index2_crcy,
	string index2_tenor,
	string index2_type,

	vector<string> spread_tenor,
	vector<string> spread_type,

	vector<string> index1_freq,
	vector<string> index2_freq,

	vector<double> mkt_spread,
	double dr,
	vector<vector<double>>& dzdr
);

void cszerocurve
(
	CCurrency crcy,
	CDate today,
	vector<string> tenor,
	vector<string> type,
	vector<string> fixedrateleg_freq,
	vector<string> floatingrateleg_freq,
	vector<double> mktrate,
	CCurrency index_crcy,
	string index_tenor,
	string index_type,
	CInterpolation fzcdf,
	CInterpolation fzcest,
	vector<CDate>& startdate,
	vector<CDate>& maturity,
	vector<double>& df,
	vector<double>& zero
);

