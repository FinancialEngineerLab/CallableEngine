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

