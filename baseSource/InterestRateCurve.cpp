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
#include "ConstantNUmbers.h"
#include "Currency.h"

#include "Date.h"
#include "Calendar.h"

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
)
{
	CInstrument depo(crcy, "DEPO"), swap(crcy, "SWAP"), sfut(crcy, "SFUT"), fxswap(crcy, "FX"), bond(crcy, "BOND");
	
	CDate *depoholiday = depo.get_holiday(),
		*swapholiday = swap.get_holiday(),
		*sfutholiday = swap.get_holiday(),
		*fxswapholiday = fxswap.get_holiday(),
		*fxswapstartholiday = fxswap.get_startholiday(),
		*bondholiday = bond.get_holiday();

	int numdepoholiday = depo.get_numholiday(),
		numswapholiday = swap.get_numholiday(),
		numsfutholiday = sfut.get_numholiday(),
		numfxswapholiday = fxswap.get_numholiday(),
		numfxswapstartholiday = fxswap.get_numstartholiday(),
		numbondholiday = bond.get_numholiday();

	string depostub = depo.get_stub(),
		depodirection = depo.get_direction(),
		depoconv = depo.get_conv(),
		depoadjflag = depo.get_adjflag(),
		depopayin = depo.get_payin(),
		deposetin = depo.get_setin(),
		depobasis = depo.get_basis(),
		swapstub = swap.get_stub(),
		swapdirection = swap.get_direction(),
		swapconv = swap.get_conv(),
		swapadjflag = swap.get_adjflag(),
		swappayin = swap.get_payin(),
		swapsetin = swap.get_setin(),
		swapbasis = swap.get_basis(),
		sfutconv = sfut.get_conv(),
		sfutbasis = sfut.get_basis(),
		fxswapconv = fxswap.get_conv(),
		fxswapbasis = fxswap.get_basis(),
		bondstub = bond.get_stub(),
		bonddirection = bond.get_direction(),
		bondconv = bond.get_conv(),
		bondadjflag = bond.get_adjflag(),
		bondpayin = bond.get_payin(),
		bondsetin = bond.get_setin(),
		bondbasis = bond.get_basis();

	int depospotlag = depo.get_spotlag(),
		depofixlag = -depospotlag,
		swapspotflag = swap.get_spotlag(),
		swapfixlag = -swapspotflag,
		futmm = sfut.get_futmm(),
		futnth = sfut.get_futnth(),
		futday = sfut.get_futday(),
		sfutfixlag = -sfut.get_spotlag(),
		fxswaplag = fxswap.get_spotlag(),
		bondspotlag = bond.get_spotlag(),
		bondfixlag = -bondspotlag;

	int i, j, num_mktrate = int(tenor.size());

	vector<CDate> calstartdate(num_mktrate), calenddate(num_mktrate), paydate(num_mktrate);

	vector<vector<double>> amt(num_mktrate);

	CDate spotdate, fxswapspotdate;
	
	ShiftBusDate(today, depoholiday, numdepoholiday, depospotlag, spotdate);
	ShiftBusDate(today, fxswapstartholiday, numfxswapstartholiday, fxswaplag, fxswapspotdate);

	double notional = 1.0, spotdf = 1.0, fxswapspotdf = 1.0;

	for (i = 0; i < num_mktrate; i++)
	{
		if (type[i] == "DEPO")
		{
			amt[i] = vector<double>(1);

			depositeschedule
			(
				today
				, tenor[i]
				, depoholiday
				, numdepoholiday
				, depostub
				, depodirection
				, depoconv
				, depoadjflag
				, depopayin
				, depospotlag
				, calstartdate[i]
				, calenddate[i]
				, paydate[i]
			);
			
			depositeamount
			(
				notional
				, calstartdate[i]
				, paydate[i]
				, mktrate[i]
				, depobasis
				, amt[i][0]
			);
			
			df.push_back(spotdf*notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
			if (calstartdate[i] < spotdate)
			{
				spotdf = df[i];
			}

			maturity.push_back(paydate[i]);

		}
		else if (type[i] == "SFUT")
		{
			amt[i] = vector<double>(1);
			CDate tmplasttradedate;
			shortfuturesschedule
			(
				tenor[i],
				sfutholiday,
				numsfutholiday,
				futmm,
				futnth,
				futday,
				sfutconv,
				sfutfixlag,
				calstartdate[i],
				calenddate[i],
				tmplasttradedate
			);

			paydate[i] = calenddate[i];

			shortfuturesamount(notional, calstartdate[i], paydate[i], mktrate[i], sfutbasis, amt[i][0]);

			if (calstartdate[i]<=maturity[i-1])
			{
				CInterpolation zc(maturity, zero);
				double start_ds = zc(calstartdate[i]);
				double start_df = exp(-start_ds*cvg(today,calstartdate[i],"ACT/365"));
				df.push_back(spotdf*notional / amt[i][0]);
				zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
				maturity.push_back(paydate[i]);

			}
			else
			{
				maturity.push_back(paydate[i]);
				double t0 = cvg(today, maturity[i - 1], "ACT/365"),
					t1 = cvg(today, calstartdate[i], "ACT/365"),
					t2 = cvg(today, maturity[i], "ACT/365");
				zero.push_back(((t2-t0)*log(amt[i][0]/notional)+zero[i-1]*(t2-t1)+t1)/(t2*(t2-t0)-t1*(t1-t0)));
				df.push_back(exp(-zero[i] * cvg(today, paydate[i], "ACT/365")));
			}

		}
		else if (type[i] == "FX")
		{
			amt[i] = vector<double>(1);
			fxswapschedule(today, tenor[i], fxswapholiday, numfxswapholiday, fxswapstartholiday, numfxswapstartholiday, fxswapconv, fxswaplag, calstartdate[i], calenddate[i]);
			paydate[i] = calenddate[i];
			amt[i][0] = notional * (1.0 + mktrate[i] * cvg(calstartdate[i], paydate[i], fxswapbasis));
			df.push_back(fxswapspotdf*notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
			if (calstartdate[i]<fxswapspotdate)
			{
				fxswapspotdf = df[i];
				maturity.push_back(paydate[i]);
				CInterpolation zc(maturity, zero);
				spotdf = exp(-zc(spotdate)*cvg(today, spotdate, "ACT/365"));
			}
		}
		else if (type[i] == "SWAP")
		{
			CInterpolation zc(maturity, zero);
			
			int num_prevschdule = int(maturity.size()),
				num_fixedrateleg_schedule,
				num_fixedrateleg_freq = Nummonthayear / int(YMD2I(fixedrateleg_freq[i])[1]);
			
			double tmpaccamt = 0.0
				, tmprmdamt = 0.0;
			
			vector<vector<CDate>> fixedrateleg_schedule(3);
			
			calstartdate[i] = spotdate;

			num_fixedrateleg_schedule = findnumschedule(calstartdate[i], tenor[i], swapstub, swapdirection, fixedrateleg_freq[i]);
			
			for (j=0;j<3;j++)
			{
				fixedrateleg_schedule[j] = vector<CDate>(num_fixedrateleg_schedule);
			}

			fixedlegcashflowschedule
			(
				calstartdate[i]
				, tenor[i]
				, swapholiday
				, numswapholiday
				, swapstub
				, swapdirection
				, swapconv
				, fixedrateleg_freq[i]
				, swapadjflag
				, num_fixedrateleg_schedule
				, fixedrateleg_schedule
			);

			maturity.push_back(fixedrateleg_schedule[0][num_fixedrateleg_schedule - 1]);

			vector<double> coupon(num_fixedrateleg_schedule), interval(num_fixedrateleg_schedule);

			for (j=0;j<num_fixedrateleg_schedule;j++)
			{
				coupon[j] = mktrate[i] * cvg(fixedrateleg_schedule[1][j], fixedrateleg_schedule[2][j], swapbasis);
				interval[j] = cvg(today, fixedrateleg_schedule[0][j], "ACT/365");
			}
			for (j = 0; j < num_fixedrateleg_schedule; j++)
			{
				if (fixedrateleg_schedule[0][j]<=maturity[i-1])
				{
					tmpaccamt = tmpaccamt + exp(-zc(fixedrateleg_schedule[0][j])*interval[j])*notional*coupon[j];
				}
				else
				{
					break;
				}
			}
			tmprmdamt = notional * spotdf - tmpaccamt;

			zero.push_back(0.0);

			findzerorate(notional, tmprmdamt, zero[i-1], coupon, interval, maturity[i-1], maturity[i], fixedrateleg_schedule[0], j, zero[i]);

			df.push_back(exp(-zero[i] * interval[num_fixedrateleg_schedule - 1]));
		}
		else
		{
		}
	}		
}

void zerocurve
(
	CCurrency crcy,
	CDate today,
	vector<string> tenor,
	vector<string> type,
	vector<string> fixedrateleg_freq,
	vector<string> floatingrateleg_freq,
	vector<double> mktrate,

	vector<CDate>& startdate, // ***
	vector<CDate>& maturity,

	vector<double>& df,
	vector<double>& zero
)
{
	//ofstream foutdf("df.txt");

	CInstrument depo(crcy, "DEPO"), swap(crcy, "SWAP"), sfut(crcy, "SFUT"), fxswap(crcy, "FX"), bond(crcy, "BOND");

	CDate* depoholiday = depo.get_holiday(),
		*swapholiday = swap.get_holiday(),
		*sfutholiday = swap.get_holiday(),
		*fxswapholiday = fxswap.get_holiday(),
		*fxswapstartholiday = fxswap.get_startholiday(),
		*bondholiday = bond.get_holiday();

	int numdepoholiday = depo.get_numholiday(),
		numswapholiday = swap.get_numholiday(),
		numsfutholiday = sfut.get_numholiday(),
		numfxswapholiday = fxswap.get_numholiday(),
		numfxswapstartholiday = fxswap.get_numstartholiday(),
		numbondholiday = bond.get_numholiday();

	string depostub = depo.get_stub(),
		depodirection = depo.get_direction(),
		depoconv = depo.get_conv(),
		depoadjflag = depo.get_adjflag(),
		depopayin = depo.get_payin(),
		deposetin = depo.get_setin(),
		depobasis = depo.get_basis(),
		swapstub = swap.get_stub(),
		swapdirection = swap.get_direction(),
		swapconv = swap.get_conv(),
		swapadjflag = swap.get_adjflag(),
		swappayin = swap.get_payin(),
		swapsetin = swap.get_setin(),
		swapbasis = swap.get_basis(),
		sfutconv = sfut.get_conv(),
		sfutbasis = sfut.get_basis(),
		fxswapconv = fxswap.get_conv(),
		fxswapbasis = fxswap.get_basis(),
		bondstub = bond.get_stub(),
		bonddirection = bond.get_direction(),
		bondconv = bond.get_conv(),
		bondadjflag = bond.get_adjflag(),
		bondpayin = bond.get_payin(),
		bondsetin = bond.get_setin(),
		bondbasis = bond.get_basis();

	int depospotlag = depo.get_spotlag(),
		depofixlag = -depospotlag,
		swapspotflag = swap.get_spotlag(),
		swapfixlag = -swapspotflag,
		futmm = sfut.get_futmm(),
		futnth = sfut.get_futnth(),
		futday = sfut.get_futday(),
		sfutfixlag = -sfut.get_spotlag(),
		fxswaplag = fxswap.get_spotlag(),
		bondspotlag = bond.get_spotlag(),
		bondfixlag = -bondspotlag;

	int i, j, num_mktrate = int(tenor.size());

	vector<CDate> calstartdate(num_mktrate), calenddate(num_mktrate), paydate(num_mktrate);

	vector<vector<double>> amt(num_mktrate);
	CDate spotdate, fxswapspotdate, bondspotdate;

	ShiftBusDate(today, depoholiday, numdepoholiday, depospotlag, spotdate);
	ShiftBusDate(today, fxswapstartholiday, numfxswapstartholiday, fxswaplag, fxswapspotdate);
	ShiftBusDate(today, bondholiday, numbondholiday, bondspotlag, bondspotdate);

	double notional = 1.0, spotdf = 1.0, fxswapspotdf = 1.0;

	for (i = 0; i < num_mktrate; i++)
	{
		if (type[i] == "DEPO")
		{
			amt[i] = vector<double>(1);
			depositeschedule(today, tenor[i], depoholiday, numdepoholiday, depostub, depodirection, depoconv, depoadjflag, depopayin, depospotlag, calstartdate[i], calenddate[i], paydate[i]);
			depositeamount(notional, calstartdate[i], paydate[i], mktrate[i], depobasis, amt[i][0]);
			df.push_back(spotdf*notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
			
			if (calstartdate[i] < spotdate)
			{
				spotdf = df[i];
			}
			maturity.push_back(paydate[i]);
			startdate.push_back(calstartdate[i]);
		}
		else if (type[i] == "SFUT")
		{
			amt[i] = vector<double>(1);
			CDate tmplasttradedate;
			shortfuturesschedule(tenor[i], sfutholiday, numsfutholiday, futmm, futnth, futday, sfutconv, sfutfixlag, calstartdate[i], calenddate[i], tmplasttradedate);
			paydate[i] = calenddate[i];
			shortfuturesamount(notional, calstartdate[i], paydate[i], mktrate[i], sfutbasis, amt[i][0]);

			if (calstartdate[i] <= maturity[i - 1])
			{
				CInterpolation zc(maturity, zero);
				double start_ds = zc(calstartdate[i]);
				double start_df = exp(-start_ds * cvg(today, calstartdate[i], "ACT/365"));
				df.push_back(spotdf*notional / amt[i][0]);
				zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
				maturity.push_back(paydate[i]);
				startdate.push_back(calstartdate[i]);

			}
			else
			{
				maturity.push_back(paydate[i]);
				double t0 = cvg(today, maturity[i - 1], "ACT/365"),
					t1 = cvg(today, calstartdate[i], "ACT/365"),
					t2 = cvg(today, maturity[i], "ACT/365");
				zero.push_back(((t2 - t0)*log(amt[i][0] / notional) + zero[i - 1] * (t2 - t1) + t1) / (t2*(t2 - t0) - t1 * (t1 - t0)));
				df.push_back(exp(-zero[i] * cvg(today, paydate[i], "ACT/365")));
				startdate.push_back(calstartdate[i]);

			}

		}
		else if (type[i] == "FX")
		{
			amt[i] = vector<double>(1);
			fxswapschedule(today, tenor[i], fxswapholiday, numfxswapholiday, fxswapstartholiday, numfxswapstartholiday, fxswapconv, fxswaplag, calstartdate[i], calenddate[i]);
			paydate[i] = calenddate[i];
			amt[i][0] = notional * (1.0 + mktrate[i] * cvg(calstartdate[i], paydate[i], fxswapbasis));
			df.push_back(fxswapspotdf*notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
			if (calstartdate[i] < fxswapspotdate)
			{
				fxswapspotdf = df[i];
				maturity.push_back(paydate[i]);

				startdate.push_back(calstartdate[i]);
				
				CInterpolation zc(maturity, zero);
				spotdf = exp(-zc(spotdate)*cvg(today, spotdate, "ACT/365"));
			}
		}
		else if (type[i] == "SWAP")
		{
			CInterpolation zc(maturity, zero);

			int num_prevschdule = int(maturity.size()),
				num_fixedrateleg_schedule,
				num_fixedrateleg_freq = Nummonthayear / int(YMD2I(fixedrateleg_freq[i])[1]);

			double tmpaccamt = 0.0,
				tmprmdamt;

			vector<vector<CDate>> fixedrateleg_schedule(3);

			calstartdate[i] = spotdate;

			num_fixedrateleg_schedule = findnumschedule(calstartdate[i], tenor[i], swapstub, swapdirection, fixedrateleg_freq[i]);

			for (j = 0; j < 3; j++)
			{
				fixedrateleg_schedule[j] = vector<CDate>(num_fixedrateleg_schedule);
			}

			fixedlegcashflowschedule(calstartdate[i], tenor[i], swapholiday, numswapholiday, swapstub, swapdirection, swapconv, fixedrateleg_freq[i], swapadjflag, num_fixedrateleg_schedule, fixedrateleg_schedule);

			maturity.push_back(fixedrateleg_schedule[0][num_fixedrateleg_schedule - 1]);

			startdate.push_back(calstartdate[i]);

			vector<double> coupon(num_fixedrateleg_schedule), interval(num_fixedrateleg_schedule);

			for (j = 0; j < num_fixedrateleg_schedule; j++)
			{
				coupon[j] = mktrate[i] * cvg(fixedrateleg_schedule[1][j], fixedrateleg_schedule[2][j], swapbasis);
				interval[j] = cvg(today, fixedrateleg_schedule[0][j], "ACT/365");
			}
			for (j = 0; j < num_fixedrateleg_schedule; j++)
			{
				if (fixedrateleg_schedule[0][j] <= maturity[i - 1])
				{
					tmpaccamt = tmpaccamt + exp(-zc(fixedrateleg_schedule[0][j])*interval[j])*notional*coupon[j];
				}
				else
				{
					break;
				}
			}
			tmprmdamt = notional * spotdf - tmpaccamt;

			zero.push_back(0.0);

			findzerorate(notional, tmprmdamt, zero[i - 1], coupon, interval, maturity[i - 1], maturity[i], fixedrateleg_schedule[0], j, zero[i]);

			df.push_back(exp(-zero[i] * interval[num_fixedrateleg_schedule - 1]));
		}
		else if(type[i] == "KRWKTB")
		{
			amt[i] = vector<double>(1);
			depositeschedule
			(
				today,
				tenor[i],
				bondholiday,
				numbondholiday,
				bondstub,
				bonddirection,
				bondconv,
				bondadjflag,
				bondpayin,
				bondspotlag,
				calstartdate[i],
				calenddate[i],
				paydate[i]
			);

			depositeamount
			(
				notional,
				calstartdate[i],
				paydate[i],
				mktrate[i],
				bondbasis,
				amt[i][0]
			);

			if (calstartdate[i]<bondspotdate)
			{
				spotdf = df[i];
			}

			df.push_back(spotdf*notional/amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));

			if (calstartdate[i]<bondspotdate)
			{
				spotdf = df[i];
			}

			maturity.push_back(paydate[i]);
			startdate.push_back(calstartdate[i]);

		}
		else
		{
			CInterpolation zc(maturity, zero);

			char yy[2], cpn[5], yyyy[65];

			int ny, ncpn, bondnamei, k = -1, k0;

			double ytm, unitprice = 10000.0, targetbondprice = 0.0, tmpdf = 0.0, accrueddf = 0.0, tmpbondprice = 0.0, tmprmdamt;

			bondnamei = int(type[i].find_first_of("y")) - 3;

			ny = int(type[i]._Copy_s(yy, 2, 2, bondnamei + 4));

			ncpn = int(type[i]._Copy_s(cpn, 5, 5, bondnamei + 4));

			ny = atoi(yy);

			string temptenor;

			_itoa_s(ny, yyyy, 65, 10);

			temptenor.append(yyyy);
			temptenor.append("y");

			CDate mat(tenor[i]);
			CDate issuedate = DateAdd('y', -ny, mat);

			vector<vector<CDate>> bond_schedule(3);

			int num_bondcaashflowschedule = int(ny * Nummonthayear / int(YMD2I(fixedrateleg_freq[i])[1]));

			for (j=0; j<3; j++)
			{
				bond_schedule[j] = vector<CDate>(num_bondcaashflowschedule);
			}

			fixedlegcashflowschedule
			(
				issuedate,
				temptenor,
				bondholiday,
				numbondholiday,
				bondstub,
				bonddirection,
				"NOTHING",
				fixedrateleg_freq[i],
				"NOTHING",
				num_bondcaashflowschedule,
				bond_schedule
			);

			maturity.push_back(bond_schedule[0][num_bondcaashflowschedule - 1]);
			startdate.push_back(bondspotdate);

			vector<double> coupon(num_bondcaashflowschedule), interval(num_bondcaashflowschedule);

			for (j=0; j<num_bondcaashflowschedule; j++)
			{
				coupon[j] = double(atoi(cpn)) / unitprice / double(num_bondcaashflowschedule / ny) * 0.1;
				interval[j] = cvg(today, bond_schedule[0][j], "ACT/365");
			}

			ytm = mktrate[i] / double(num_bondcaashflowschedule / ny);

			for (j=0; j<num_bondcaashflowschedule; j++)
			{
				if (bond_schedule[0][j]>=today)
				{
					k = k + 1;
					tmpdf = pow(1.0 / (1.0 + ytm), k);
					targetbondprice = targetbondprice + tmpdf * coupon[j];
				}
			}

			k0 = num_bondcaashflowschedule - k - 1;
			targetbondprice = targetbondprice + tmpdf;

			if (k0 == 0)
			{
				accrueddf = 1.0 / (1.0 + ytm * (bond_schedule[2][k0] - bondspotdate) / (bond_schedule[2][k0] - bond_schedule[1][k0]));
			}
			else
			{
				accrueddf = 1.0 / (1.0 + ytm * (bond_schedule[2][k0] - bondspotdate) / (bond_schedule[2][k0] - bond_schedule[2][k0 - 1]));
			}

			targetbondprice *= accrueddf;
			targetbondprice = int(targetbondprice * unitprice) / 100.;
			targetbondprice = targetbondprice * exp(-zc(bondspotdate) * cvg(today, bondspotdate, "ACT/365"));

			for (j=k0; j<num_bondcaashflowschedule; j++)
			{
				if (bond_schedule[0][j]<=maturity[i-1])
				{
					tmpbondprice = tmpbondprice + exp(-zc(bond_schedule[0][j]) * interval[j]) * coupon[j];
				}
				else
				{
					break;
				}
			}

			tmprmdamt = (targetbondprice / 100. - tmpbondprice) * notional;

			zero.push_back(0.0);

			findzerorate(notional, tmprmdamt, zero[i - 1], coupon, interval, maturity[i - 1], maturity[i], bond_schedule[0], j, zero[i]);

			df.push_back(exp(-zero[i] * interval[num_bondcaashflowschedule - 1]));

		}
	}
}

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
)
{
	sum = 0.0;
	int i = 0;
	for (i = startindx; i<num_coup; i++)
	{
		sum = sum + notional * coupon[i] * exp(-(((enddate - interdate[i])*startzero + (interdate[i] - startdate)*endzero) / (enddate - startdate))*interval[i]);
	}
	sum = sum + notional * exp(-endzero * interval[num_coup - 1]);
}

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
	double &zero
)
{
	int num_coup = int(interdate.size());
	double left = -1.0, right = 1.0, mid = 0.5 * (left + right), ans;
	interpolatedamt(notional, startzero, mid, coupon, interval, startdate, enddate, interdate, startindx, num_coup, ans);
	ans = ans - amt;

	while (fabs(ans)>Error)
	{
		if (ans > 0.0)
		{
			left = mid;
			mid = 0.5 * (left + right);
			interpolatedamt(notional, startzero, mid, coupon, interval, startdate, enddate, interdate, startindx, num_coup, ans);
			ans = ans - amt;
		}
		else if (ans > 0.0)
		{
			right = mid;
			mid = 0.5 * (left + right);
			interpolatedamt(notional, startzero, mid, coupon, interval, startdate, enddate, interdate, startindx, num_coup, ans);
			ans = ans - amt;
		}
		else
		{
			break;
		}
	}
	zero = mid;
}

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
	vector<vector<double>> &dzdr
)
{
	int i, j, num_mktrate = int(tenor.size());
	vector<double> tmpmktrate(num_mktrate), df, zero;
	vector<CDate> mtrty, startdate;

	zerocurve(crcy, today, tenor, type, fixedrateleg_freq, floatingrateleg_freq, mktrate, startdate, mtrty, df, zero);

	for (i = 0; i < num_mktrate; i++)
	{
		tmpmktrate = mktrate;
		for (j=0;j<num_mktrate; j++)
		{
			if (i==j)
			{
				tmpmktrate[j] = tmpmktrate[j] + dr;
			}
		}

		vector<double> tmpdf, tmpzero;
		vector<CDate> tmpmtrty, tmpstartdate;

		zerocurve(crcy, today, tenor, type, fixedrateleg_freq, floatingrateleg_freq, tmpmktrate, tmpstartdate, tmpmtrty, tmpdf, tmpzero);

		for (j=i; j<num_mktrate; j++)
		{
			dzdr[i][j] = (tmpzero[j] - zero[j]) / dr;
		}
	}

}

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
)
{
	CInstrument fxswap(crcy, "FX");
	
	CDate* fxswapholiday = fxswap.get_holiday(), * fxswapstartholiday = fxswap.get_startholiday();

	int numfxswapholiday = fxswap.get_numholiday(), numfxswapstartholiday = fxswap.get_numstartholiday(), fxswaplag = fxswap.get_spotlag();

	string fxswapconv = fxswap.get_conv(), fxswapbasis = fxswap.get_basis();

	int i, num_mktquote = int(swappointtenor.size());

	CDate fxswapspotdate;

	ShiftBusDate
	(
		today,
		fxswapstartholiday,
		numfxswapstartholiday,
		fxswaplag,
		fxswapspotdate
	);

	vector<double> newswap_point;
	vector<CDate> calstartdate(num_mktquote), calenddate(num_mktquote);

	for (i=0; i<num_mktquote; i++)
	{
		fxswapschedule
		(
			today,
			swappointtenor[i],
			fxswapholiday,
			numfxswapholiday,
			fxswapstartholiday,
			numfxswapstartholiday,
			fxswapconv,
			fxswaplag,
			calstartdate[i],
			calenddate[i]
		);

		if (i>0)
		{
			if (calenddate[i]>calenddate[i-1])
			{
				resulttenor.push_back(swappointtenor[i]);
				fxswapstartdate.push_back(calstartdate[i]);
				fxswapmaturity.push_back(calenddate[i]);

				newswap_point.push_back(swap_point[i]);

				fxrate.push_back(0.0);

				marketquote.push_back(0.0);
			}
		}
		else
		{
			resulttenor.push_back(swappointtenor[i]);
			fxswapstartdate.push_back(calstartdate[i]);
			fxswapmaturity.push_back(calenddate[i]);

			newswap_point.push_back(swap_point[i]);

			fxrate.push_back(0.0);

			marketquote.push_back(0.0);

		}
	}

	num_mktquote = int(resulttenor.size());
	double swappointsum = 0.0, df1, df2, cap, t1, t2, tau1, tau2;

	for (i = num_mktquote-1; i>=0; i--)
	{
		if (fxswapstartdate[i]<fxswapspotdate)
		{
			swappointsum = swappointsum + newswap_point[i];
			fxrate[i] = fxspot - swappointsum;
		}
		else
		{
			fxrate[i] = fxspot + newswap_point[i];
		}
	}

	if (swappt_mode)
	{
		for (i=0; i<num_mktquote; i++)
		{
			if (fxswapstartdate[i]<fxswapspotdate)
			{
				if (i<num_mktquote-1)
				{
					t1 = cvg(today, fxswapstartdate[i], "ACT/365");
					t2 = cvg(today, fxswapmaturity[i], "ACT/365");

					tau1 = cvg(today, fxswapstartdate[i], fxswapbasis);
					tau2 = cvg(today, fxswapmaturity[i], fxswapbasis);

					df1 = exp(-zc(fxswapstartdate[i]) * t1);
					df2 = exp(-zc(fxswapmaturity[i]) * t2);

					if (fxswapstartdate[i+1]<fxswapspotdate)
					{
						cap = fxrate[i + 1] / fxrate[i] * df1 / df2;
					}
					else
					{
						cap = fxspot / fxrate[i] * df1 / df2;
					}

					marketquote[i] = (cap - 1.0) / (tau2 - tau1);
				}
			}
			else
			{
				t1 = cvg(today, fxswapstartdate[i], "ACT/365");
				t2 = cvg(today, fxswapmaturity[i], "ACT/365");

				tau1 = cvg(today, fxswapstartdate[i], fxswapbasis);
				tau2 = cvg(today, fxswapmaturity[i], fxswapbasis);

				df1 = exp(-zc(fxswapstartdate[i]) * t1);
				df2 = exp(-zc(fxswapmaturity[i]) * t2);

				cap = fxrate[i] / fxspot * df1 / df2;
				marketquote[i] = (cap - 1.0) / (tau2 - tau1);
			}
		}
	}
	else
	{
		for (i=0; i<num_mktquote; i++)
		{
			if (fxswapstartdate[i] < fxswapspotdate)
			{
				if (i<num_mktquote-1)
				{
					t1 = cvg(today, fxswapstartdate[i], "ACT/365");
					t2 = cvg(today, fxswapmaturity[i], "ACT/365");

					tau1 = cvg(today, fxswapstartdate[i], fxswapbasis);
					tau2 = cvg(today, fxswapmaturity[i], fxswapbasis);

					df1 = exp(-zc(fxswapstartdate[i]) * t1);
					df2 = exp(-zc(fxswapmaturity[i]) * t2);

					if (fxswapstartdate[i + 1] < fxswapspotdate)
					{
						cap = fxrate[i + 1] * df2 / (fxrate[i] * df1);
					}
					else
					{
						cap = fxspot*df2 / (fxrate[i] * df1);
					}

					marketquote[i] = (1.0 / cap - 1.0) / (tau2 - tau1);
				}
			}
			else
			{
				t1 = cvg(today, fxswapstartdate[i], "ACT/365");
				t2 = cvg(today, fxswapmaturity[i], "ACT/365");

				tau1 = cvg(today, fxswapstartdate[i], fxswapbasis);
				tau2 = cvg(today, fxswapmaturity[i], fxswapbasis);

				df1 = exp(-zc(fxswapstartdate[i]) * t1);
				df2 = exp(-zc(fxswapmaturity[i]) * t2);

				cap = fxspot * df2 / (fxspot * df1);

				marketquote[i] = (1.0 / cap - 1.0) / (tau2 - tau1);
			}
		}
	}
}

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
)
{
	CInstrument 
		basis(crcy, "BASIS"), 
		depo(crcy, "DEPO"), 
		fxswap(crcy, "FX"), 
		index1(index1_crcy, index1_type), 
		index2(index2_crcy, index2_type);

	CDate
		* basisholiday = basis.get_holiday(),
		* fxswapholiday = fxswap.get_holiday(),
		* fxswapstartholiday = fxswap.get_startholiday(),

		* index1holiday = index1.get_holiday(),
		* index2holiday = index2.get_holiday(),

		* index1fixingholiday = index1.get_fixingholiday(),
		* index2fixingholiday = index2.get_fixingholiday(),

		* depoholiday = depo.get_holiday();

	int numbasisholiday = basis.get_numholiday(),
		numfxswapholiday = fxswap.get_numholiday(),
		numfxswapstartholiday = fxswap.get_numstartholiday(),
		numindex1holiday = index1.get_numholiday(),
		numindex2holiday = index2.get_numholiday(),
		numindex1fixingholiday = index1.get_numfixingholiday(),
		numindex2fixingholiday = index2.get_numfixingholiday(),
		numdepoholiday = depo.get_numholiday();

	string
		basisstub = basis.get_stub(),
		index1stub = index1.get_stub(),
		index2stub = index2.get_stub(),

		basisdirection = basis.get_direction(),

		index1direction = index1.get_direction(),
		index2direction = index2.get_direction(),

		basisconv = basis.get_conv(),
		index1conv = index1.get_conv(),
		index2conv = index2.get_conv(),

		basisadjflag = basis.get_adjflag(),
		basispayin = basis.get_payin(),
		basissetin = basis.get_setin(),

		basisdcb = basis.get_basis(),
		index1dcb = index1.get_basis(),
		index2dcb = index2.get_basis(),

		fxswapconv = fxswap.get_conv(),
		fxswapbasis = fxswap.get_basis(),
		index1setin = index1.get_setin(),
		index2setin = index2.get_setin();

	int basisspotlag = basis.get_spotlag(),
		basisfixlag = -basisspotlag,

		fxswaplag = fxswap.get_spotlag(),
		index1spotlag = index1.get_spotlag(),
		index2spotlag = index2.get_spotlag();

	string dcb = "ACT/365";

	int i, j, num_mktspread = int(spread_tenor.size());

	vector<CDate>
		calstartdate(num_mktspread),
		calenddate(num_mktspread),
		paydate(num_mktspread);

	CDate
		spotdate,
		fxswapspotdate,
		index1_spotdate,
		index2_spotdate;

	ShiftBusDate
	(
		today,
		basisholiday,
		numbasisholiday,
		basisspotlag,
		spotdate
	);

	ShiftBusDate
	(
		today,
		fxswapstartholiday,
		numfxswapstartholiday,
		fxswaplag,
		fxswapspotdate
	);

	ShiftBusDate
	(
		today,
		index1holiday,
		numindex1holiday,
		index1spotlag,
		index1_spotdate
	);

	ShiftBusDate
	(
		today,
		index2holiday,
		numindex2holiday,
		index2spotlag,
		index2_spotdate
	);

	double
		notional = 1.0,
		spotdf = 1.0,
		fxswapspotdf = 1.0,
		index1tempsum = 0.0,
		index2tempsum = 0.0;

	vector<int> ymd;

	int
		indx1_tenor,
		indx2_tenor,
		prev_num_index1_cf = 0,
		prev_num_index2_cf = 0,

		num_index1_cf,
		num_index2_cf;

	ymd = YMD2I(index1_tenor);

	indx1_tenor = ymd[0] * Nummonthayear + ymd[1];

	ymd = YMD2I(index2_tenor);

	indx2_tenor = ymd[0] * Nummonthayear + ymd[1];

	vector<double>
		index1fra,
		index2fra,
		index1partamt,
		index2partamt,
		index1amt,
		index2amt,
		index1cf, 
		index2cf,
		index1df,
		index2df;

	vector<CDate> basis_maturity;

	for (i=0; i<num_mktspread; i++)
	{
		basis_startdate.push_back(index2_spotdate);

		CDate tempmatyrity, temprealmaturity;

		ymd = YMD2I(spread_tenor[i]);

		int m = ymd[0] * Nummonthayear + ymd[1];

		basis_maturity.push_back(DateAdd('m', m, basis_startdate[i]));

		ConvReflDate(basis_maturity[i], basisholiday, numbasisholiday, basisconv, temprealmaturity);

		basis_mtrty.push_back(temprealmaturity);

		basis_zero.push_back(0.0);
	}

	CInterpolation bszc(basis_mtrty, basis_zero), sprdzc = bszc + zc;

	spread_startdate = vector<CDate>(int(sprdzc.get_xds_size()));

	merge(basis_startdate.begin(), basis_startdate.end(), startdate.begin(), startdate.end(), spread_startdate.begin());

	spotdf = exp(-sprdzc(index2_spotdate) * cvg(today, index2_spotdate, dcb));

	vector<vector<CDate>> start, end, fixingend;
	vector<vector<double>> df;

	CDate valdate;

	for (i=0; i<num_mktspread; i++)
	{
		num_index1_cf = findnumschedule(index1_spotdate, spread_tenor[i], index1stub, index1direction, indx1_tenor);
		num_index2_cf = findnumschedule(index2_spotdate, spread_tenor[i], index2stub, index2direction, indx2_tenor);

		int indx1_freq, indx2_freq;

		ymd = YMD2I(index1_freq[i]);

		indx1_freq = ymd[0] * Nummonthayear + ymd[1];

		ymd = YMD2I(index2_freq[i]);

		indx2_freq = ymd[0] * Nummonthayear + ymd[1];

		int
			indx1_freqtopay = indx1_freq / indx1_tenor,
			indx2_freqtopay = indx2_freq / indx2_tenor;

		vector<CDate>
			index1_calc_startdate(num_index1_cf),
			index1_calc_enddate(num_index1_cf),

			index1_paydate(num_index1_cf),
			index1_fixingdate(num_index1_cf),

			fixingindex1_maturity(num_index1_cf),

			index2_calc_startdate(num_index2_cf),
			index2_calc_enddate(num_index2_cf),

			index2_paydate(num_index2_cf),
			index2_fixingdate(num_index2_cf),

			fixingindex2_maturity(num_index2_cf);

		calculationstartschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index1_tenor,
			basisadjflag,
			num_index1_cf,
			index1_calc_startdate
		);

		calculationendschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index1_tenor,
			basisadjflag,
			num_index1_cf,
			index1_calc_enddate
		);

		paymentschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index1_tenor,
			num_index1_cf,
			index1_paydate
		);

		fixingschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,

			index1fixingholiday,
			numindex1fixingholiday,

			basisstub,
			basisdirection,
			basisconv,

			index1_tenor,
			basisadjflag,

			index1setin,
			
			-index1spotlag,
			num_index1_cf,

			index1_fixingdate

		);

		index1tempsum = 0.0;

		for (j = prev_num_index1_cf; j < num_index1_cf; j++)
		{
			index1df.push_back(exp(-zc(index1_paydate[j]) * cvg(today, index1_paydate[j], dcb)));

			ShiftBusDate(index1_fixingdate[j], index1fixingholiday, numindex1fixingholiday, index1spotlag, fixingindex1_maturity[j]);

			findmaturity(fixingindex1_maturity[j], index1_tenor, depoholiday, numdepoholiday, index1conv, fixingindex1_maturity[j]);

			index1fra.push_back((exp(-zc(index1_calc_startdate[j]) * cvg(today, index1_calc_startdate[j], dcb) + zc(fixingindex1_maturity[j]) * cvg(today, fixingindex1_maturity[j], dcb)) - 1.0) / cvg(index1_calc_startdate[j], fixingindex1_maturity[j], index1dcb));
		}

		vector<double> index1cf;

		// int tmpprev_num_index1_cf;

		for (j=0; j<num_index1_cf; j++)
		{
			index1cf.push_back((index1fra[j] + mkt_spread[j]) * cvg(index1_calc_startdate[j], index1_calc_enddate[j], basisdcb));

			/*
			if (j % indx1_freqtopay==indx1_freqtopay-1)
			{
				index1tempamt
			}
			*/

			index1tempsum = index1tempsum + index1cf[j] * index1df[j];
		}

		index1partamt.push_back(index1tempsum);
		index1amt.push_back(index1tempsum);

		calculationstartschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index2_tenor,
			basisadjflag,
			num_index2_cf,
			index2_calc_startdate
		);

		calculationendschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index2_tenor,
			basisadjflag,
			num_index2_cf,
			index2_calc_enddate
		);

		paymentschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,
			basisstub,
			basisdirection,
			basisconv,
			index2_tenor,
			num_index2_cf,
			index2_paydate
		);

		fixingschedule
		(
			basis_startdate[i],
			basis_maturity[i],
			basisholiday,
			numbasisholiday,

			index2fixingholiday,
			numindex2fixingholiday,

			basisstub,
			basisdirection,
			basisconv,
			index2_tenor,
			basisadjflag,

			index2setin,
			-index2spotlag,
			num_index2_cf,
			index2_fixingdate
		);

		for (j=prev_num_index2_cf; j<num_index2_cf; j++)
		{
			index2df.push_back(exp(-zc(index2_paydate[j]) * cvg(today, index2_paydate[j], dcb)));
		}

		vector<CDate> tmpstart, tmpend, tmpfixingend;

		vector<double> tmpdf;

		for (j=0; j<num_index2_cf; j++)
		{
			ShiftBusDate(index2_fixingdate[j], index2fixingholiday, numindex2fixingholiday, index2spotlag, fixingindex2_maturity[j]);

			findmaturity(fixingindex2_maturity[j], index2_tenor, depoholiday, numdepoholiday, index2conv, fixingindex2_maturity[j]);

			tmpstart.push_back(index2_calc_startdate[j]);

			tmpend.push_back(index2_calc_enddate[j]);

			tmpfixingend.push_back(fixingindex2_maturity[j]);

			tmpdf.push_back(index2df[j]);
		}

		start.push_back(tmpstart);
		end.push_back(tmpend);

		fixingend.push_back(tmpfixingend);

		df.push_back(tmpdf);

		prev_num_index1_cf = num_index1_cf;
		prev_num_index2_cf = num_index2_cf;

	}

	vector<double> x(num_mktspread);

	for (i = 0; i < num_mktspread; i++)
	{
		x[i] = mkt_spread[i];
	}

	const long ntrial = 100000;
	const double tolx = Error;
	const double tolf = Error;

	findbasiszeromnewt
	(
		today,
		index2_spotdate,
		spotdf,
		basis_maturity,
		num_mktspread,
		start,
		end,
		fixingend,
		df,
		zc,
		index2dcb,
		index1partamt,
		ntrial,
		basis_zero,
		tolx,
		tolf
	);

	for (i=0; i<num_mktspread; i++)
	{
		basis_df.push_back(exp(-basis_zero[i] * cvg(today, basis_maturity[i], dcb)));
	}

	ofstream foutbszc("bszc.txt");

	foutbszc << "Basis Zero Curve" << endl;
	foutbszc << setprecision(15);

	for (i=0; i<num_mktspread; i++)
	{
		foutbszc << basis_df[i] << "\t" << basis_zero[i] * 100 << "\t" << basis_maturity[i].get_str_date() << endl;
	}

}


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

	vector<double> &x,

	vector<double> &fvec,

	vector<vector<double>> &fjac
)
{
	int i, j, k, l;

	string dcb = "ACT/365";

	vector<vector<double>>
		start_tenor(num_mktspread),
		fixingend_tenor(num_mktspread),
		start_today(num_mktspread),
		fixingend_today(num_mktspread),
		start_fixingend(num_mktspread),
		start_end(num_mktspread),
		zc_start(num_mktspread),
		zc_fixingend(num_mktspread),
		dfsp_start(num_mktspread),
		dfsp_fixingend(num_mktspread),
		fra(num_mktspread),
		y_start(num_mktspread),
		y_fixingend(num_mktspread);

	vector<double> btwtenor(num_mktspread);

	vector<vector<vector<double>>> dy_startdx(num_mktspread), dy_fixingenddx(num_mktspread);

	vector<int> add_num_cf(num_mktspread), num_cf(num_mktspread);

	for (i=0; i<num_mktspread; i++)
	{
		num_cf[i] = int(start[i].size());
		add_num_cf[i] = num_cf[i];

		if (i>0)
		{
			add_num_cf[i] = num_cf[i] - num_cf[i - 1];
		}

		start_tenor[i] = vector<double>(add_num_cf[i]);

		fixingend_tenor[i] = vector<double>(add_num_cf[i]);

		start_today[i] = vector<double>(add_num_cf[i]);

		fixingend_today[i] = vector<double>(add_num_cf[i]);

		start_fixingend[i] = vector<double>(add_num_cf[i]);

		start_end[i] = vector<double>(add_num_cf[i]);

		zc_start[i] = vector<double>(add_num_cf[i]);
		zc_fixingend[i] = vector<double>(add_num_cf[i]);

		dfsp_start[i] = vector<double>(add_num_cf[i]);
		dfsp_fixingend[i] = vector<double>(add_num_cf[i]);

		fra[i] = vector<double>(add_num_cf[i]);

		y_start[i] = vector<double>(add_num_cf[i]);
		y_fixingend[i] = vector<double>(add_num_cf[i]);

		dy_startdx[i] = vector<vector<double>>(add_num_cf[i]);
		dy_fixingenddx[i] = vector<vector<double>>(add_num_cf[i]);

		for (j=0; j<add_num_cf[i]; j++)
		{
			dy_startdx[i][j] = vector<double>(num_mktspread, 0.0);
			dy_fixingenddx[i][j] = vector<double>(num_mktspread, 0.0);
		}
	}

	btwtenor[0] = cvg(spotdate, tenor[0], dcb);

	for (i=1; i<num_mktspread; i++)
	{
		btwtenor[i] = cvg(tenor[i - 1], tenor[i], dcb);
	}

	for (j=0; j<add_num_cf[0]; j++)
	{
		start_tenor[0][j] = cvg(spotdate, start[0][j], dcb);
		fixingend_tenor[0][j] = cvg(spotdate, fixingend[0][j], dcb);
		start_today[0][j] = cvg(today, start[0][j], dcb);
		fixingend_today[0][j] = cvg(today, fixingend[0][j], dcb);
		
		y_start[0][j] = x[0];
		y_fixingend[0][j] = x[0];

		dy_startdx[0][j][0] = 1.0;
		dy_fixingenddx[0][j][0] = 1.0;

	}

	y_start[0][0] = 0.0;
	dy_startdx[0][0][0] = 0.0;

	if (fixingend[0][add_num_cf[0]-1]>tenor[0])
	{
		fixingend_tenor[0][add_num_cf[0] - 1] = cvg(tenor[0], fixingend[0][add_num_cf[0] - 1], dcb);

		if (num_mktspread>1)
		{
			y_fixingend[0][add_num_cf[0] - 1] = (x[1] - x[0]) / btwtenor[1] * fixingend_tenor[0][add_num_cf[0] - 1] + x[0];

			dy_fixingenddx[0][add_num_cf[0] - 1][0] = 1.0 - fixingend_tenor[0][add_num_cf[0] - 1] / btwtenor[1];

			dy_fixingenddx[0][add_num_cf[0] - 1][1] = fixingend_tenor[0][add_num_cf[0] - 1] / btwtenor[1];
		}
	}

	for (i=1; i<num_mktspread; i++)
	{
		for (j=0; j<add_num_cf[i]; j++)
		{
			start_tenor[i][j] = cvg(tenor[i-1], start[i][num_cf[i-1]+j], dcb);
			fixingend_tenor[i][j] = cvg(tenor[i - 1], fixingend[i][num_cf[i - 1] + j], dcb);

			start_today[i][j] = cvg(today, start[i][num_cf[i - 1] + j], dcb);
			fixingend_today[0][j] = cvg(today, fixingend[i][num_cf[i - 1] + j], dcb);

			y_start[i][j] = (x[i]- x[i-1])/btwtenor[i]*start_tenor[i][j]+x[i-1];
			y_fixingend[i][j] = (x[i] - x[i - 1]) / btwtenor[i] * fixingend_tenor[i][j] + x[i - 1];

			dy_startdx[i][j][i - 1] = 1.0 - start_tenor[i][j] / btwtenor[i];
			dy_startdx[i][j][i] = start_tenor[i][j] / btwtenor[i];

			dy_fixingenddx[i][j][i-1] = 1.0 - fixingend_tenor[i][j]/btwtenor[i];
			dy_fixingenddx[i][j][i] = fixingend_tenor[i][j] / btwtenor[i];
		}

		if (fixingend[i][add_num_cf[i]-1]>tenor[i])
		{
			fixingend_tenor[i][add_num_cf[i] - 1] = cvg(tenor[i], fixingend[i][add_num_cf[i] - 1], dcb);

			if (i<num_mktspread-1)
			{
				y_fixingend[i][add_num_cf[i] - 1] = (x[i + 1] - x[i]) / btwtenor[i + 1] * fixingend_tenor[i][add_num_cf[i] - 1] + x[i];

				dy_fixingenddx[i][add_num_cf[i] - 1][i] = 1.0 - fixingend_tenor[i][add_num_cf[i] - 1] / btwtenor[i + 1];

				dy_fixingenddx[i][add_num_cf[i] - 1][i+1] = fixingend_tenor[i][add_num_cf[i] - 1] / btwtenor[i + 1];

			}
			else
			{
				y_fixingend[i][add_num_cf[i] - 1] = x[i];
				dy_fixingenddx[i][add_num_cf[i] - 1][i] = 1.0;
			}
		}
	}

	for (j=0;j<add_num_cf[0]; j++)
	{
		start_fixingend[0][j] = cvg(start[0][j], fixingend[0][j], indexdcb);
		start_end[0][j] = cvg(start[0][j], end[0][j], indexdcb);

		zc_start[0][j] = zc(start[0][j]);
		zc_fixingend[0][j] = zc(fixingend[0][j]);

		dfsp_start[0][j] = exp(-(y_start[0][j] + zc_start[0][j]) * start_today[0][j]);
		dfsp_fixingend[0][j] = exp(-(y_fixingend[0][j] + zc_fixingend[0][j]) * fixingend_today[0][j]);

		fra[0][j] = (dfsp_start[0][j] / dfsp_fixingend[0][j] - 1.0) / start_fixingend[0][j];
	}

	for (i=1; i<num_mktspread; i++)
	{
		for (j = 0; j < add_num_cf[i]; j++)
		{
			start_fixingend[i][j] = cvg(start[i][num_cf[i - 1] + j], fixingend[i][num_cf[i - 1] + j], indexdcb);
			start_end[i][j] = cvg(start[i][num_cf[i - 1] + j], end[i][num_cf[i - 1] + j], indexdcb);

			zc_start[i][j] = zc(start[i][num_cf[i - 1] + j]);
			zc_fixingend[i][j] = zc(fixingend[i][num_cf[i - 1] + j]);

			dfsp_start[i][j] = exp(-(y_start[i][j] + zc_start[i][j]) * start_today[i][j]);
			dfsp_fixingend[i][j] = exp(-(y_fixingend[i][j] + zc_fixingend[i][j]) * fixingend_today[i][j]);

			fra[i][j] = (dfsp_start[i][j] / dfsp_fixingend[i][j] - 1.0) / start_fixingend[i][j];
		}
	}

	for (i=0; i<num_mktspread; i++)
	{
		fvec[i] = 0.0;

		for (l=0; l<add_num_cf[0];l++)
		{
			fvec[i] = fvec[i] + fra[0][l] * start_end[0][l] * df[0][l];
		}

		for (k=1; k<=i; k++)
		{
			for (l=0; l<add_num_cf[k]; l++)
			{
				fvec[i] = fvec[i] + fra[k][l] * start_end[k][l] * df[k][num_cf[k - 1] + l];
			}
		}
		fvec[i] = fvec[i] - amt[i];
	}

	for (i=0; i<num_mktspread; i++)
	{
		for (j = 0; j < num_mktspread; j++)
		{
			fjac[i][j] = 0.0;

			for (l=0; l<add_num_cf[0];l++)
			{
				fjac[i][j] = fjac[i][j] + (fra[k][l] * start_fixingend[0][l] + 1.0) * (-start_today[0][l] * dy_startdx[0][l][j] + fixingend_today[0][l] * dy_fixingenddx[0][l][j]) / start_fixingend[0][l] * start_end[0][l] * df[0][l];
			}


			for (k=1; k<=i; k++)
			{
				for (l=0; l<add_num_cf[k]; l++)
				{
					fjac[i][j] = fjac[i][j] + 
						(fra[k][l] * start_fixingend[k][l] + 1.0) *
						(-start_today[k][l] * dy_startdx[k][l][j] + fixingend_today[k][l] * dy_fixingenddx[k][l][j]) /
						start_fixingend[k][l] * start_end[k][l] * df[k][num_cf[k - 1] + l];
				}
			}
		}
	}
}

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
)
{
	int k, i;
	double errx, errf, d;
	int n = int(x.size());

	vector<int>  indx(n);

	vector<double> p(n), fvec(n);
	vector<vector<double>> fjac(n);

	for (i = 0; i < n; i++)
	{
		fjac[i] = vector<double>(n);
	}

	for (k = 0; k < ntrial; k++)
	{
		basiszero
		(
			today,
			spotdate,
			spotdf,
			tenor,
			num_mktspread,
			start,
			end,
			fixingend,
			df,
			zc,
			dcb,
			amt,
			x,
			fvec,
			fjac
		);

		errf = 0.0;

		for (i = 0; i < n; i++)
		{
			errf += fabs(fvec[i]);
		}

		if (errf <= tolf)
		{
			return;
		}

		for (i = 0; i < n; i++)
		{
			p[i] = -fvec[i];
		}

		ludcmp(fjac, indx, d);
		lubksb(fjac, indx, p);

		errx = 0.0;

		for (i = 0; i < n; i++)
		{
			errx += fabs(p[i]);
			x[i] += p[i];
		}
		if (errx < tolx)
		{
			return;
		}
	}

	return;
}

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
	vector<vector<double>> &dzdr
)
{
	int i, j, num_mktspread = int(spread_tenor.size());

	vector<double> tmpmkt_spread(num_mktspread), basis_df, basis_zero, spread_df, spread_zero;

	vector<CDate> basis_startdate, basis_mtrty, spread_mtrty, spread_startdate;

	basiszerocurve
	(
		crcy,
		today,
		startdate,
		mtrty,
		zc,
		fxrate,
		index1_crcy,
		index1_tenor,
		index1_type,
		index2_crcy,
		index2_tenor,
		index2_type,
		spread_tenor,
		spread_type,
		index1_freq,
		index2_freq,
		mkt_spread,
		basis_startdate,
		basis_mtrty,
		basis_df,
		basis_zero,
		spread_startdate,
		spread_mtrty,
		spread_df,
		spread_zero
	);

	for (i=0; i<num_mktspread; i++)
	{
		tmpmkt_spread = mkt_spread;

		for (j=0; j<num_mktspread; j++)
		{
			if (i==j)
			{
				tmpmkt_spread[j] = tmpmkt_spread[j] + dr;
			}
		}

		vector<double> tmpmktrate(num_mktspread), tmpbasis_df, tmpbasis_zero, tmpspread_df, tmpspread_zero;

		vector<CDate> tmpbasis_startdate, tmpbasis_mtrty, tmpspread_mtrty, tmpspread_startdate;

		basiszerocurve
		(
			crcy,
			today,
			startdate,
			mtrty,
			zc,
			fxrate,
			index1_crcy,
			index1_tenor,
			index1_type,
			index2_crcy,
			index2_tenor,
			index2_type,
			spread_tenor,
			spread_type,
			index1_freq,
			index2_freq,
			tmpmkt_spread,
			tmpbasis_startdate,
			tmpbasis_mtrty,
			tmpbasis_df,
			tmpbasis_zero,
			tmpspread_startdate,
			tmpspread_mtrty,
			tmpspread_df,
			tmpspread_zero
		);

		for (j=i; j<num_mktspread; j++)
		{
			dzdr[i][j] = (tmpbasis_zero[j] - basis_zero[j]) / dr;
		}

	}
}

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
	vector<CDate> &startdate,
	vector<CDate> &maturity,
	vector<double> &df,
	vector<double> &zero
)
{
	CInstrument
		depo(crcy, "DEPO"),
		swap(crcy, "SWAP"),
		sfut(crcy, "SFUT"),
		fxswap(crcy, "FX"),
		bond(crcy, "BOND"),
		index(index_crcy, index_type);

	CDate
		* depoholiday = depo.get_holiday(),
		* swapholiday = swap.get_holiday(),
		* sfutholiday = swap.get_holiday(),
		* fxswapholiday = fxswap.get_holiday(),
		* fxswapstartholiday = fxswap.get_startholiday(),
		* bondholiday = bond.get_holiday(),
		* indexholiday = index.get_holiday(),
		* indexfixingholiday = index.get_fixingholiday();

	int
		numdepoholiday = depo.get_numholiday(),
		numswapholiday = swap.get_numholiday(),
		numsfutholiday = sfut.get_numholiday(),
		numfxswapholiday = fxswap.get_numholiday(),
		numfxswapstartholiday = fxswap.get_numstartholiday(),
		numbondholiday = bond.get_numholiday(),
		numindexholiday = index.get_numholiday(),
		numindexfixingholiday = index.get_numfixingholiday();

	string
		depostub = depo.get_stub(),
		depodirection = depo.get_direction(),
		depoconv = depo.get_conv(),
		depoadjflag = depo.get_adjflag(),
		depopayin = depo.get_payin(),
		deposetin = depo.get_setin(),
		depobasis = depo.get_basis(),
		swapstub = swap.get_stub(),
		swapdirection = swap.get_direction(),
		swapconv = swap.get_conv(),
		swapadjflag = swap.get_adjflag(),
		swappayin = swap.get_payin(),
		swapsetin = swap.get_setin(),
		swapbasis = swap.get_basis(),

		sfutconv = sfut.get_conv(),
		sfutbasis = sfut.get_basis(),

		fxswapconv = fxswap.get_conv(),
		fxswapbasis = fxswap.get_basis(),

		bondstub = bond.get_stub(),
		bonddirection = bond.get_direction(),
		bondconv = bond.get_conv(),
		bondadjflag = bond.get_adjflag(),
		bondpayin = bond.get_payin(),
		bondsetin = bond.get_setin(),
		bondbasis = bond.get_basis(),

		indexstub = index.get_stub(),
		indexdirection = index.get_direction(),
		indexconv = index.get_conv(),
		indexdcb = index.get_basis(),
		indexsetin = index.get_setin();

	int
		depospotlag = depo.get_spotlag(),
		depofixlag = -depospotlag,
		swapspotlag = swap.get_spotlag(),
		swapfixlag = -swapspotlag,
		futmm = sfut.get_futmm(),
		futnth = sfut.get_futnth(),
		futday = sfut.get_futday(),
		sfutfixlag = -sfut.get_spotlag(),
		fxswaplag = fxswap.get_spotlag(),
		bondspotlag = bond.get_spotlag(),
		bondfixlag = -bondspotlag,
		bondquotedigit = bond.get_bondquotedigit(),
		indexspotlag = index.get_spotlag();

	int i, j, num_mktrate = int(tenor.size()), num_index_cf;

	vector<CDate> calstartdate(num_mktrate), calendate(num_mktrate), paydate(num_mktrate);

	vector<vector<double>> amt(num_mktrate);

	CDate spotdate, fxswapspotdate, bondspotdate;

	ShiftBusDate(today, depoholiday, numdepoholiday, depospotlag, spotdate);
	ShiftBusDate(today, fxswapstartholiday, numfxswapstartholiday, fxswaplag, fxswapspotdate);
	ShiftBusDate(today, bondholiday, numbondholiday, bondspotlag, bondspotdate);

	double notional = 1.0, spotdf = 1.0, fxswapspotdf = 1.0;

	vector<int> ymd;

	ymd = YMD2I(index_tenor);

	int indx_tenor = ymd[0] * Nummonthayear + ymd[1];

	for (i=0; i<num_mktrate; i++)
	{
		if (type[i]=="DEPO")
		{
			amt[i] = vector<double>(1);

			depositeschedule
			(
				today,
				tenor[i],
				depoholiday,
				numdepoholiday,
				depostub,
				depodirection,
				depoconv,
				depoadjflag,
				depopayin,
				depospotlag,
				calstartdate[i],
				calendate[i],
				paydate[i]
			);

			depositeamount
			(
				notional,
				calstartdate[i],
				paydate[i],
				mktrate[i],
				depobasis,
				amt[i][0]
			);

			df.push_back(spotdf * notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));

			if (calstartdate[i] < spotdate)
			{
				spotdf = df[i];
			}

			maturity.push_back(paydate[i]);
			startdate.push_back(calstartdate[i]);

		}
		else if (type[i] == "SFUT")
		{
			amt[i] = vector<double>(1);

			CDate tmplasttradedate;

			shortfuturesschedule
			(
				tenor[i],
				sfutholiday,
				numsfutholiday,
				futmm,
				futnth,
				futday,
				sfutconv,
				sfutfixlag,
				calstartdate[i],
				calendate[i],
				tmplasttradedate
			);

			paydate[i] = calendate[i];

			shortfuturesamount
			(
				notional,
				calstartdate[i],
				paydate[i],
				mktrate[i],
				sfutbasis,
				amt[i][0]
			);

			if (calstartdate[i]<=maturity[i-1])
			{
				CInterpolation zc(maturity, zero);

				double start_ds = zc(calstartdate[i]);
				double start_df = exp(-start_ds * cvg(today, calstartdate[i], "ACT/365"));

				df.push_back(start_df * notional / amt[i][0]);
				zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));
				maturity.push_back(paydate[i]);
				startdate.push_back(calstartdate[i]);
			}
			else
			{
				maturity.push_back(paydate[i]);
				double
					t0 = cvg(today, maturity[i - 1], "ACT/365"),
					t1 = cvg(today, calstartdate[i], "ACT/365"),
					t2 = cvg(today, maturity[i], "ACT/365");

				zero.push_back(((t2 - t0) * log(amt[i][0] / notional) + zero[i - 1] * (t2 - t1) * t1) / (t2 * (t2 - t0) - t1 * (t1 - t0)));
				df.push_back(exp(-zero[i] * cvg(today, paydate[i], "ACT/365")));
				startdate.push_back(calstartdate[i]);
			}
		}
		else if (type[i] == "FX")
		{
			amt[i] = vector<double>(1);

			fxswapschedule
			(
				today,
				tenor[i],
				fxswapholiday,
				numfxswapholiday,
				fxswapstartholiday,
				numfxswapstartholiday,
				fxswapconv,
				fxswaplag,
				calstartdate[i],
				calendate[i]
			);

			amt[i][0] = notional * (1.0 + mktrate[i] * cvg(calstartdate[i], paydate[i], fxswapbasis));

			df.push_back(fxswapspotdf* notional / amt[i][0]);
			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));

			if (calstartdate[i]<fxswapspotdate)
			{
				fxswapspotdf = df[i];
			}

			maturity.push_back(paydate[i]);
			startdate.push_back(calstartdate[i]);
			CInterpolation zc(maturity, zero);
			spotdf = exp(-zc(spotdate) * cvg(today, spotdate, "ACT/365"));

		}
		else if (type[i] == "SWAP")
		{
			CInterpolation zc(maturity, zero);
			int
				num_prevschedule = int(maturity.size()),
				num_fixedrateleg_schedule,
				num_fixedrateleg_freq = Nummonthayear / int(YMD2I(fixedrateleg_freq[i])[1]);

			double
				tmpaccamt = 0.0,
				tmprmdamt;

			vector<vector<CDate>> fixedrateleg_schedule(3);

			calstartdate[i] = spotdate;

			num_fixedrateleg_schedule = findnumschedule
										(
											calstartdate[i],
											tenor[i],
											swapstub,
											swapdirection,
											fixedrateleg_freq[i]
										);

			num_index_cf = findnumschedule
							(
								calstartdate[i],
								tenor[i],
								indexstub,
								indexdirection,
								indx_tenor
							);

			for (j=0; i<3; j++)
			{
				fixedrateleg_schedule[j] = vector<CDate>(num_fixedrateleg_schedule);
			}

			fixedlegcashflowschedule
			(
				calstartdate[i],
				tenor[i],
				swapholiday,
				numswapholiday,
				swapstub,
				swapdirection,
				swapconv,
				fixedrateleg_freq[i],
				swapadjflag,
				num_fixedrateleg_schedule,
				fixedrateleg_schedule
			);

			maturity.push_back(fixedrateleg_schedule[0][num_fixedrateleg_schedule - 1]);
			startdate.push_back(calstartdate[i]);

			vector<CDate>
				index_calc_startdate(num_index_cf),
				index_calc_enddate(num_index_cf),
				index_paydate(num_index_cf),
				index_fixingdate(num_index_cf),
				fixingindex_maturity(num_index_cf);

			calculationstartschedule
			(
				calstartdate[i],
				tenor[i],
				swapholiday,
				numswapholiday,
				swapstub,
				swapdirection,
				swapconv,
				index_tenor,
				swapadjflag,
				num_index_cf,
				index_calc_startdate
			);

			calculationendschedule
			(
				calstartdate[i],
				tenor[i],
				swapholiday,
				numswapholiday,
				swapstub,
				swapdirection,
				swapconv,
				index_tenor,
				swapadjflag,
				num_index_cf,
				index_calc_enddate
			);

			paymentschedule
			(
				calstartdate[i],
				tenor[i],
				swapholiday,
				numswapholiday,
				swapstub,
				swapdirection,
				swapconv,
				index_tenor,
				num_index_cf,
				index_paydate
			);

			fixingschedule
			(
				calstartdate[i],
				tenor[i],
				swapholiday,
				numswapholiday,
				indexfixingholiday,
				numindexfixingholiday,
				swapstub,
				swapdirection,
				swapconv,
				index_tenor,
				swapadjflag,
				indexsetin,
				-indexspotlag,
				num_index_cf,
				index_fixingdate
			);

			double famt = -notional * exp(-fzcdf(spotdate) * cvg(today, spotdate, "ACT/365"));

			vector<double> 
				indexfra(num_index_cf),
				indexdf(num_index_cf);

			for (j=0; j<num_index_cf; j++)
			{
				ShiftBusDate
				(
					index_fixingdate[j],
					indexfixingholiday,
					numindexfixingholiday,
					indexspotlag,
					fixingindex_maturity[j]
				);

				findmaturity
				(
					fixingindex_maturity[j],
					index_tenor,
					indexholiday,
					numindexfixingholiday,
					indexconv,
					fixingindex_maturity[j]
				);

				indexfra[j] 
					=  (
						exp(
							-fzcest(index_calc_startdate[j]) * cvg(today, index_calc_startdate[j], "ACT/365")
							+ fzcest(fixingindex_maturity[j]) * cvg(today, fixingindex_maturity[j], "ACT/365")
							) - 1.0
						) / cvg(index_calc_startdate[j], fixingindex_maturity[j], indexdcb);

				indexdf[j] = exp(-fzcdf(index_paydate[j]) * cvg(today, index_paydate[j], "ACT/365"));

				famt = famt + notional * indexfra[j] * cvg(index_calc_startdate[j], index_calc_enddate[j], indexdcb) * indexdf[j];
			}

			famt = famt + notional * indexdf[num_index_cf - 1];

			vector<double>
				coupon(num_fixedrateleg_schedule),
				interval(num_fixedrateleg_schedule);

			for (j=0; j<num_fixedrateleg_schedule; j++)
			{
				coupon[j] = mktrate[i] * cvg(fixedrateleg_schedule[1][j], fixedrateleg_schedule[2][j], swapbasis);
				interval[j] = cvg(today, fixedrateleg_schedule[0][j], "ACT/365");
			}

			for (j = 0; j < num_fixedrateleg_schedule; j++)
			{
				if (fixedrateleg_schedule[0][j] <= maturity[i-1])
				{
					tmpaccamt = tmpaccamt + exp(-zc(fixedrateleg_schedule[0][j]) * interval[j]) * notional * coupon[j];
				}
				else
				{
					break;
				}
			}

			tmprmdamt = notional * spotdf - tmpaccamt + famt * exp(fzcdf(spotdate) * cvg(today, spotdate, "ACT/365")) * spotdf;

			zero.push_back(0.0);

			findzerorate
			(
				notional,
				tmprmdamt,
				zero[i-1],
				coupon,
				interval,
				maturity[i-1],
				maturity[i],
				fixedrateleg_schedule[0],
				j,
				zero[i]
			);

			df.push_back(exp(-zero[i] * interval[num_fixedrateleg_schedule - 1]));
		}
		else if (type[i] == "KRWKTB")
		{
			amt[i] = vector<double>(1);

			depositeschedule
			(
				today,
				tenor[i],
				bondholiday,
				numbondholiday,
				bondstub,
				bonddirection,
				bondconv,
				bondadjflag,
				bondpayin,
				bondspotlag,
				calstartdate[i],
				calendate[i],
				paydate[i]
			);

			depositeamount
			(
				notional,
				calstartdate[i],
				paydate[i],
				mktrate[i],
				bondbasis,
				amt[i][0]
			);

			df.push_back(spotdf* notional / amt[i][0]);

			zero.push_back(-log(df[i]) / cvg(today, paydate[i], "ACT/365"));

			if (calstartdate[i] < bondspotdate)
			{
				spotdf = df[i];
			}

			maturity.push_back(paydate[i]);
			startdate.push_back(calstartdate[i]);
		}
		else
		{
			CInterpolation zc(maturity, zero);

			char
				yy[2],
				cpn[4],
				yyyy[65];

			int
				ny,
				ncpn,
				bondnamei,
				k = -1,
				k0;

			double
				ytm,
				unitprice = 10000.0,
				targetbondprice = 0.0,
				tmpdf = 0.0,
				accrueddf = 0.0,
				tmpbondprice = 0.0,
				tmprmdamt;

			bondnamei = int(type[i].find_first_of("y")) - 3;
			
			ny = int(type[i]._Copy_s(yy, 2, 2, bondnamei + 1));
			
			ncpn = int(type[i]._Copy_s(cpn, 4, 4, bondnamei + 4));
			
			ny = atoi(yy);

			string temptenor;

			_itoa_s(ny, yyyy, 65, 10);

			temptenor.append(yyyy);
			temptenor.append("y");

			CDate mat(tenor[i]);
			CDate issuedate = DateAdd('y', -ny, mat);

			vector<vector<CDate>> bond_schedule(3);

			int num_bondcashflowschedule = int(ny * Nummonthayear / int(YMD2I(fixedrateleg_freq[i])[1]));

			for (j=0; j<3; j++)
			{
				bond_schedule[j] = vector<CDate>(num_bondcashflowschedule);
			}

			fixedlegcashflowschedule
			(
				issuedate
				,temptenor
				,bondholiday
				,numbondholiday
				,bondstub
				,bonddirection
				,"NOTHING"
				,fixedrateleg_freq[i]
				,"NOTHING"
				,num_bondcashflowschedule
				,bond_schedule
			);

			maturity.push_back(bond_schedule[0][num_bondcashflowschedule - 1]);
			startdate.push_back(bondspotdate);

			vector<double>
				coupon(num_bondcashflowschedule),
				interval(num_bondcashflowschedule);

			for (j=0; j<num_bondcashflowschedule; j++)
			{
				coupon[j] = double(atoi(cpn)) / unitprice / double(num_bondcashflowschedule / ny);
				interval[j] = cvg(today, bond_schedule[0][j], "ACT/365");
			}

			ytm = mktrate[i] / double(num_bondcashflowschedule / ny);

			for (j=0; j<num_bondcashflowschedule; j++)
			{
				if (bond_schedule[0][j]>=today)
				{
					k = k + 1;
					tmpdf = pow(1.0 / (1.0 + ytm), k);
					targetbondprice = targetbondprice + tmpdf * coupon[j];
				}
			}

			k0 = num_bondcashflowschedule - k - 1;
			targetbondprice = targetbondprice + tmpdf;

			if (k0 == 0)
			{
				accrueddf = 1.0 / (1.0 + ytm + (bond_schedule[2][k0] - bondspotdate) / (bond_schedule[2][k0] - bond_schedule[1][k0]));
			}
			else
			{
				accrueddf = 1.0 / (1.0 + ytm + (bond_schedule[2][k0] - bondspotdate) / (bond_schedule[2][k0] - bond_schedule[2][k0 - 1]));
			}

			targetbondprice *= accrueddf;
			targetbondprice = int(targetbondprice * unitprice) / 100;
			targetbondprice = targetbondprice * exp(-zc(bondspotdate) * cvg(today, bondspotdate, "ACT/365")) * 100;

			for (j=k0; j<num_bondcashflowschedule; j++)
			{
				if (bond_schedule[0][j]<=maturity[i-1])
				{
					tmpbondprice = tmpbondprice + exp(-zc(bond_schedule[0][j]) * interval[j]) * coupon[j];
				}
				else
				{
					break;
				}
			}

			tmprmdamt = (targetbondprice / 100. - tmpbondprice) * notional;

			zero.push_back(0.0);

			findzerorate(notional, tmprmdamt, zero[i - 1], coupon, interval, maturity[i - 1], maturity[i], bond_schedule[0], j, zero[i]);

			df.push_back(exp(-zero[i] * interval[num_bondcashflowschedule - 1]));
		}
	}

}




