#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "PlainPricer.h"

// Added by hokim 2020-04-13
#include "Currency.h"

// Added by hokim 2020-04-02

#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <numeric>

// 


void PlainSwapPrice
(
    CCurrency crcy, CDate today
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
)
{
    //ofstream fout("D:/HOKim_BackUp/000_Calc/Callprob/callprob.txt");

    string dcb = "ACT/365";

    double ondf = 1.0;

    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;

    double
        settlement_date_df
        = exp
        (
            -zc(settlement_date)
            * cvg(today, settlement_date, dcb)
        )
        , level = 0.0
        , levelcvg = 0.0;

    int
        i
        , starti = 0
        , endi
        , num_couponleg_cf = int(couponleg_calc_startdate.size())
        , num_fundingleg_cf = int(fundingleg_calc_startdate.size());

    vector<double>
        couponleg_df(num_couponleg_cf, 1.0)
        , coupon(num_couponleg_cf)
        , fra(num_fundingleg_cf, 0.0)
        , fundingleg_df(num_fundingleg_cf, 1.0);

    while(today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_fundingleg_cf - 1)
        {
            break;
        }
    }

    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i]
            = exp
            (
                -zc(fundingleg_paydate[i])
                * cvg(today, fundingleg_paydate[i], dcb)
            );

        fra[i] = fixinghistory_rate[i];
        //if(fra[i]>0.0) floatinglegprice=floatinglegprice+fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];

        if(fra[i] > MinFixingRate)
        {
            floatinglegprice
                = floatinglegprice
                + fundingleg_notional[i]
                * (
                    fundingleg_mult[i]
                    * fra[i]
                    + fundingleg_spread[i]
                    )
                * cvg
                (
                    fundingleg_calc_startdate[i]
                    , fundingleg_calc_enddate[i]
                    , fundingleg_dcb
                )
                * fundingleg_df[i];
        }

        if (today <= fixing_date[i])
        {
            break;
        }
    }

    if (today == fundingleg_paydate[num_fundingleg_cf - 1])
    {
        floatinglegprice = 0.0;
    }

    endi = i;

    if(endi < num_fundingleg_cf)
    {
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");

            ShiftBusDate
            (
                today
                , temp.get_holiday()
                , temp.get_numholiday()
                , 1
                , ondate
            );

            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }

    settlement_date_df
        = settlement_date_df
        / ondf;

    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i]
            = exp
            (
                -zc(fundingleg_paydate[i])
                * cvg(today, fundingleg_paydate[i], dcb)
            );

        fra[i]
            = (
                exp
                (
                    -zc(fundingleg_calc_startdate[i])
                    * cvg(today, fundingleg_calc_startdate[i], dcb)
                    + zc(fixingindex_maturity[i])
                    * cvg(today, fixingindex_maturity[i], dcb)
                )
                - 1.0
                )
            / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        floatinglegprice
            = floatinglegprice
            + (estflag)
            *fundingleg_notional[i]
            * fundingleg_mult[i]
            * fra[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        floatinglegprice
            = floatinglegprice
            + fundingleg_notional[i]
            * fundingleg_spread[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        //fout << setprecision(18) << fixingindex_maturity[i].get_str_date() << "\t" << fra[i] << endl;

    }

    if (endi < num_fundingleg_cf)
    {
        floatinglegprice
            = floatinglegprice
            + (1 - estflag)
            * (
                fundingleg_notional[endi]
                * exp
                (
                    -zc(fundingleg_calc_startdate[endi])
                    * cvg(today, fundingleg_calc_startdate[endi], dcb)
                )
                - fundingleg_notional[num_fundingleg_cf - 1]
                * fundingleg_df[num_fundingleg_cf - 1]
                );
    }

    starti = 0;

    while(today >= couponleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_couponleg_cf - 1)
        {
            break;
        }
    }

    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i]
            = exp
            (
                -zc(couponleg_paydate[i])
                * cvg(today, couponleg_paydate[i], dcb)
            );

        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);

        level
            = level
            + levelcvg * couponleg_df[i];

        coupon[i]
            = couponleg_couponrate[i]
            * levelcvg;

        fixedlegprice
            = fixedlegprice
            + couponleg_notional[i]
            * coupon[i]
            * couponleg_df[i];
    }

    if (today == couponleg_paydate[num_couponleg_cf - 1])
    {
        fixedlegprice = 0.0;
    }

    if (starti < num_couponleg_cf)
    {
        atmswaprate = floatinglegprice / level / notional;
    }

    floatinglegprice
        = floatinglegprice
        / ondf
        / settlement_date_df;

    fixedlegprice
        = fixedlegprice
        / ondf
        / settlement_date_df;

    swapprice
        = (
            floatinglegprice
            - fixedlegprice
            );
}

// Created by Hokim 2020-03-27
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
)
{
    //ofstream fout("D:/HOKim_BackUp/000_Calc/Callprob/callprob.txt");

    string dcb = "ACT/365";

    double ondf = 1.0;

    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;

    // created by hokim 2020-04-09
    double floatinglegaccuralprice = 0.0;
    double flixedlegaccuralprice = 0.0;
    double accrualswapprice = 0.0;

    double
        settlement_date_df
        = exp
        (
            -zc(settlement_date)
            * cvg(today, settlement_date, dcb)
        )
        , level = 0.0
        , levelcvg = 0.0;

    int
        i
        , starti = 0
        , endi
        , num_couponleg_cf = int(couponleg_calc_startdate.size())
        , num_fundingleg_cf = int(fundingleg_calc_startdate.size());

    vector<double>
        couponleg_df(num_couponleg_cf, 1.0)
        , coupon(num_couponleg_cf)
        , fra(num_fundingleg_cf, 0.0)
        , fundingleg_df(num_fundingleg_cf, 1.0);

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_fundingleg_cf - 1)
        {
            break;
        }
    }

    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i]
            = exp
            (
                -zc(fundingleg_paydate[i])
                * cvg(today, fundingleg_paydate[i], dcb)
            );

        fra[i] = fixinghistory_rate[i];
        //if(fra[i]>0.0) floatinglegprice=floatinglegprice+fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];

        if (fra[i] > MinFixingRate)
        {
            floatinglegprice
                = floatinglegprice
                + fundingleg_notional[i]
                * (
                    fundingleg_mult[i]
                    * fra[i]
                    + fundingleg_spread[i]
                    )
                * cvg
                (
                    fundingleg_calc_startdate[i]
                    , fundingleg_calc_enddate[i]
                    , fundingleg_dcb
                )
                * fundingleg_df[i];

            // created by hokim 2020-04-09
            floatinglegaccuralprice
                = floatinglegaccuralprice
                + fundingleg_notional[i]
                * (
                    fundingleg_mult[i]
                    * fra[i]
                    + fundingleg_spread[i]
                    )
                * cvg
                (
                    fundingleg_calc_startdate[i]
                    , today
                    , fundingleg_dcb
                )
                * fundingleg_df[i];
        }

        if (today <= fixing_date[i])
        {
            break;
        }
    }

    if (today == fundingleg_paydate[num_fundingleg_cf - 1])
    {
        floatinglegprice = 0.0;
    }

    endi = i;

    if (endi < num_fundingleg_cf)
    {
        //if(today==fixing_date[i] && fra[i]>0.0)
        //if(today==fixing_date[i] && fra[i]>MinFixingRate)
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");

            ShiftBusDate
            (
                today
                , temp.get_holiday()
                , temp.get_numholiday()
                , 1
                , ondate
            );

            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }

    settlement_date_df
        = settlement_date_df
        / ondf;

    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i]
            = exp
            (
                -zc(fundingleg_paydate[i])
                * cvg(today, fundingleg_paydate[i], dcb)
            );

        fra[i]
            = (
                exp
                (
                    -zc(fundingleg_calc_startdate[i])
                    * cvg(today, fundingleg_calc_startdate[i], dcb)
                    + zc(fixingindex_maturity[i])
                    * cvg(today, fixingindex_maturity[i], dcb)
                )
                - 1.0
                )
            / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        floatinglegprice
            = floatinglegprice
            + (estflag)
            *fundingleg_notional[i]
            * fundingleg_mult[i]
            * fra[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        floatinglegprice
            = floatinglegprice
            + fundingleg_notional[i]
            * fundingleg_spread[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        //fout << setprecision(18) << fixingindex_maturity[i].get_str_date() << "\t" << fra[i] << endl;

    }
    if (endi < num_fundingleg_cf)
    {
        floatinglegprice
            = floatinglegprice
            + (1 - estflag)
            * (
                fundingleg_notional[endi]
                * exp
                (
                    -zc(fundingleg_calc_startdate[endi])
                    * cvg(today, fundingleg_calc_startdate[endi], dcb)
                )
                - fundingleg_notional[num_fundingleg_cf - 1]
                * fundingleg_df[num_fundingleg_cf - 1]
                );
    }

    starti = 0;

    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_couponleg_cf - 1)
        {
            break;
        }
    }

    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i]
            = exp
            (
                -zc(couponleg_paydate[i])
                * cvg(today, couponleg_paydate[i], dcb)
            );

        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);

        level
            = level
            + levelcvg
            * couponleg_df[i];

        coupon[i]
            = couponleg_couponrate[i]
            * levelcvg;

        fixedlegprice
            = fixedlegprice
            + couponleg_notional[i]
            * coupon[i]
            * couponleg_df[i];
    }

    flixedlegaccuralprice
        = couponleg_notional[starti]
        * couponleg_couponrate[starti]
        * cvg(couponleg_calc_startdate[starti], today, couponleg_dcb)
        * couponleg_df[starti];

    if (today == couponleg_paydate[num_couponleg_cf - 1])
    {
        fixedlegprice = 0.0;
    }

    if (starti < num_couponleg_cf)
    {
        atmswaprate = floatinglegprice / level / notional;
    }

    //floatinglegprice
    //  = (floatinglegprice + fundingleg_df[num_fundingleg_cf - 1] * fundingleg_notional[num_fundingleg_cf - 1])
    //  / ondf
    //  / settlement_date_df;

    //fixedlegprice
    //  = (fixedlegprice + couponleg_df[num_couponleg_cf - 1] * couponleg_notional[num_couponleg_cf - 1])
    //  / ondf
    //  / settlement_date_df;

    swapprice
        = (
            floatinglegprice
            - fixedlegprice
            );

    accrualswapprice = floatinglegaccuralprice - flixedlegaccuralprice;
}
/* 2018/09/17 CRS pricer 수정
 *  1.  기준 커브(zc1) : settle currency, 고정금리(cs krw), 상대 커브 설정(zc2) : 변동금리 (usd)
 */
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
    , double& floatinglegprice
    , double& fixedlegprice
    , double& atmswaprate
    , double& swapprice
)
{
    ofstream fout("test.txt");
    string dcb = "ACT/365";
    double ondf = 1.0;
    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;
    double settlement_date_df = exp(-zc1(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0;
    int i, starti = 0, endi, num_couponleg_cf = int(couponleg_calc_startdate.size()), num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> couponleg_df(num_couponleg_cf, 1.0), coupon(num_couponleg_cf), fra(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0);

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1) break;
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc2(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];
        if (fra[i] > MinFixingRate) floatinglegprice = floatinglegprice + fundingleg_notional[i] * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) floatinglegprice = 0.0;
    endi = i;
    if (endi < num_fundingleg_cf)
    {
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(coupcrcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc2(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc2(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = (exp(-zc2(fundingleg_calc_startdate[i]) * cvg(today, fundingleg_calc_startdate[i], dcb) + zc2(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        floatinglegprice = floatinglegprice + (estflag)*fundingleg_notional[i] * fundingleg_mult[i] * fra[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        floatinglegprice = floatinglegprice + fundingleg_notional[i] * fundingleg_spread[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        fout << setprecision(18) << fixingindex_maturity[i].get_str_date() << "\t" << fra[i] << endl;
    }
    if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + (1 - estflag) * (fundingleg_notional[endi] * exp(-zc2(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)) - fundingleg_notional[num_fundingleg_cf - 1] * fundingleg_df[num_fundingleg_cf - 1]);

    //if(starti<=0)
    if (today <= start_date)
    {
        double fundinglegstart_date_df = exp(-zc2(start_date) * cvg(today, start_date, dcb));
        floatinglegprice = floatinglegprice + fundingleg_notional[num_fundingleg_cf] * fundinglegstart_date_df;
    }
    if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + fundingleg_notional[num_fundingleg_cf + 1] * fundingleg_df[num_fundingleg_cf - 1];
    floatinglegprice = fxspot * floatinglegprice;
    starti = 0;
    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_couponleg_cf - 1) break;
    }
    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc1(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + couponleg_notional[i] * levelcvg * couponleg_df[i];
        coupon[i] = couponleg_couponrate[i] * levelcvg;
        fixedlegprice = fixedlegprice + couponleg_notional[i] * coupon[i] * couponleg_df[i];
    }
    double tmpstartend = 0.0;
    //if(starti<=0)
    if (today <= start_date)
    {
        double couponlegstart_date_df = exp(-zc1(start_date) * cvg(today, start_date, dcb));
        tmpstartend = tmpstartend + couponleg_notional[num_couponleg_cf] * couponlegstart_date_df;
    }
    if (starti < num_couponleg_cf) tmpstartend = tmpstartend + couponleg_notional[num_couponleg_cf + 1] * couponleg_df[num_couponleg_cf - 1];
    fixedlegprice = fixedlegprice + tmpstartend;

    if (today == couponleg_paydate[num_couponleg_cf - 1]) fixedlegprice = 0.0;
    if (starti < num_couponleg_cf) atmswaprate = (-floatinglegprice - tmpstartend) / level;
    floatinglegprice = floatinglegprice / ondf / settlement_date_df;
    fixedlegprice = fixedlegprice / ondf / settlement_date_df;
    swapprice = (floatinglegprice + fixedlegprice);
}


/* 2018/09/17 CRS pricer 고정-고정 추가
 *  1.  기준 커브(zc1) : settle currency, 고정금리(cs krw), 상대 커브 설정(zc2) : 고정금리 (usd)
 */
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
)
{
    ofstream fout("test_fixed2price.txt");
    string dcb = "ACT/365";
    double ondf = 1.0;
    couponleg1price = 0.0;
    couponleg2price = 0.0;
    atmswaprate = 0.0;
    double settlement_date_df = exp(-zc1(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0, level2 = 0.0, level2cvg = 0.0;;
    int i, starti_leg1 = 0, starti_leg2 = 0, num_couponleg1_cf = int(couponleg1_calc_startdate.size()), num_couponleg2_cf = int(couponleg2_calc_startdate.size());
    vector<double> coupon1(num_couponleg1_cf), couponleg1_df(num_couponleg1_cf, 1.0), coupon2(num_couponleg2_cf), couponleg2_df(num_couponleg2_cf, 1.0);

    while (today >= couponleg2_paydate[starti_leg2])
    {
        starti_leg2 = starti_leg2 + 1;
        if (starti_leg2 >= num_couponleg2_cf - 1) break;
    }
    for (i = starti_leg2; i < num_couponleg2_cf; i++)
    {
        couponleg2_df[i] = exp(-zc2(couponleg2_paydate[i]) * cvg(today, couponleg2_paydate[i], dcb));
        level2cvg = cvg(couponleg2_calc_startdate[i], couponleg2_calc_enddate[i], couponleg2_dcb);
        level2 = level2 + couponleg2_notional[i] * level2cvg * couponleg2_df[i];
        coupon2[i] = couponleg2_couponrate[i] * level2cvg;
        couponleg2price = couponleg2price + couponleg2_notional[i] * coupon2[i] * couponleg2_df[i];
    }
    if (today == couponleg2_paydate[num_couponleg2_cf - 1]) couponleg2price = 0.0;

    settlement_date_df = settlement_date_df / ondf;

    double tmpstartend_leg2 = 0.0;
    if (today <= start_date)
    {
        double couponleg2start_date_df = exp(-zc2(start_date) * cvg(today, start_date, dcb));
        tmpstartend_leg2 = tmpstartend_leg2 + couponleg2_notional[num_couponleg2_cf] * couponleg2start_date_df;
    }
    if (starti_leg2 < num_couponleg2_cf) tmpstartend_leg2 = tmpstartend_leg2 + couponleg2_notional[num_couponleg2_cf + 1] * couponleg2_df[num_couponleg2_cf - 1];
    couponleg2price = fxspot * (couponleg2price + tmpstartend_leg2);          //(원금 + 쿠폰)*spot

    while (today >= couponleg1_paydate[starti_leg1])
    {
        starti_leg1 = starti_leg1 + 1;
        if (starti_leg1 >= num_couponleg1_cf - 1) break;
    }
    for (i = starti_leg1; i < num_couponleg1_cf; i++)
    {
        couponleg1_df[i] = exp(-zc1(couponleg1_paydate[i]) * cvg(today, couponleg1_paydate[i], dcb));
        levelcvg = cvg(couponleg1_calc_startdate[i], couponleg1_calc_enddate[i], couponleg1_dcb);
        level = level + couponleg1_notional[i] * levelcvg * couponleg1_df[i];
        coupon1[i] = couponleg1_couponrate[i] * levelcvg;
        couponleg1price = couponleg1price + couponleg1_notional[i] * coupon1[i] * couponleg1_df[i];
    }

    double tmpstartend_leg1 = 0.0;
    //if(starti<=0)
    if (today <= start_date)
    {
        double couponlegstart_date_df = exp(-zc1(start_date) * cvg(today, start_date, dcb));
        tmpstartend_leg1 = tmpstartend_leg1 + couponleg1_notional[num_couponleg1_cf] * couponlegstart_date_df;
    }
    if (starti_leg1 < num_couponleg1_cf) tmpstartend_leg1 = tmpstartend_leg1 + couponleg1_notional[num_couponleg1_cf + 1] * couponleg1_df[num_couponleg1_cf - 1];
    couponleg1price = couponleg1price + tmpstartend_leg1;

    if (today == couponleg1_paydate[num_couponleg1_cf - 1]) couponleg1price = 0.0;
    //if(starti_leg1<num_couponleg1_cf) atmswaprate=(-couponleg2price-tmpstartend)/level;
    if (starti_leg2 < num_couponleg2_cf) atmswaprate = (-couponleg1price - tmpstartend_leg2 * fxspot) / (level2 * fxspot);
    couponleg2price = couponleg2price / ondf / settlement_date_df;
    couponleg1price = couponleg1price / ondf / settlement_date_df;
    swapprice = (couponleg2price + couponleg1price);
}

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
)
{
    ofstream fout("test.txt");
    string dcb = "ACT/365";
    double ondf = 1.0;
    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;
    double settlement_date_df = exp(-cszc(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0;
    int i, starti = 0, endi, num_couponleg_cf = int(couponleg_calc_startdate.size()), num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> couponleg_df(num_couponleg_cf, 1.0), coupon(num_couponleg_cf), fra(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0);

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1) break;
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];
        //if(fra[i]>0.0) floatinglegprice=floatinglegprice+fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];
        if (fra[i] > MinFixingRate) floatinglegprice = floatinglegprice + fundingleg_notional[i] * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        //floatinglegprice=floatinglegprice+fundingleg_notional[i]*(fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) floatinglegprice = 0.0;
    endi = i;
    if (endi < num_fundingleg_cf)
    {
        //if(today==fixing_date[i] && fra[i]>0.0)
        //if(today==fixing_date[i] && fra[i]>MinFixingRate)
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(coupcrcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = (exp(-zc(fundingleg_calc_startdate[i]) * cvg(today, fundingleg_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        //fra[i]=(exp(-basiszc(fundingleg_calc_startdate[i])*cvg(today,fundingleg_calc_startdate[i],dcb)+basiszc(fixingindex_maturity[i])*cvg(today,fixingindex_maturity[i],dcb))-1.0)/cvg(fundingleg_calc_startdate[i],fixingindex_maturity[i],fixingindex_dcb);
        //floatinglegprice=floatinglegprice+(estflag)*fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];
        floatinglegprice = floatinglegprice + (estflag)*fundingleg_notional[i] * fundingleg_mult[i] * fra[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        floatinglegprice = floatinglegprice + fundingleg_notional[i] * fundingleg_spread[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        fout << setprecision(18) << fixingindex_maturity[i].get_str_date() << "\t" << fra[i] << endl;
    }
    if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + (1 - estflag) * (fundingleg_notional[endi] * exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)) - fundingleg_notional[num_fundingleg_cf - 1] * fundingleg_df[num_fundingleg_cf - 1]);

    //if(starti<=0)
    if (today <= start_date)
    {
        double fundinglegstart_date_df = exp(-zc(start_date) * cvg(today, start_date, dcb));
        floatinglegprice = floatinglegprice + fundingleg_notional[num_fundingleg_cf] * fundinglegstart_date_df;
    }
    if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + fundingleg_notional[num_fundingleg_cf + 1] * fundingleg_df[num_fundingleg_cf - 1];
    floatinglegprice = fxspot * floatinglegprice;
    starti = 0;
    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_couponleg_cf - 1) break;
    }
    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-cszc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + couponleg_notional[i] * levelcvg * couponleg_df[i];
        coupon[i] = couponleg_couponrate[i] * levelcvg;
        fixedlegprice = fixedlegprice + couponleg_notional[i] * coupon[i] * couponleg_df[i];
    }
    double tmpstartend = 0.0;
    //if(starti<=0)
    if (today <= start_date)
    {
        double couponlegstart_date_df = exp(-cszc(start_date) * cvg(today, start_date, dcb));
        tmpstartend = tmpstartend + couponleg_notional[num_couponleg_cf] * couponlegstart_date_df;
    }
    if (starti < num_couponleg_cf) tmpstartend = tmpstartend + couponleg_notional[num_couponleg_cf + 1] * couponleg_df[num_couponleg_cf - 1];
    fixedlegprice = fixedlegprice + tmpstartend;

    if (today == couponleg_paydate[num_couponleg_cf - 1]) fixedlegprice = 0.0;
    if (starti < num_couponleg_cf) atmswaprate = (-floatinglegprice - tmpstartend) / level;
    floatinglegprice = floatinglegprice / ondf / settlement_date_df;
    fixedlegprice = fixedlegprice / ondf / settlement_date_df;
    swapprice = (floatinglegprice + fixedlegprice);
}


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
)
{
    string dcb = "ACT/365";
    double ondf = 1.0;
    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;

    double
        settlement_date_df =
        exp(-zc(settlement_date) * cvg(today, settlement_date, dcb))
        , level = 0.0
        , levelcvg = 0.0;

    int
        i
        , starti = 0
        , endi
        , num_couponleg_cf = int(couponleg_calc_startdate.size())
        , num_fundingleg_cf = int(fundingleg_calc_startdate.size());

    vector<double>
        couponleg_df(num_couponleg_cf, 1.0),
        coupon(num_couponleg_cf),
        fra(num_fundingleg_cf, 0.0),
        fundingleg_df(num_fundingleg_cf, 1.0);

    double imsi_var0 = 0.0;
    double imsi_var1 = 0.0;
    double imsi_var2 = 0.0;
    double imsi_var3 = 0.0;
    double imsi_var4 = 0.0;
    double imsi_var5 = 0.0;
    double imsi_var6 = 0.0;
    double imsi_var7 = 0.0;
    double imsi_var8 = 0.0;
    double imsi_var9 = 0.0;
    double imsi_var10 = 0.0;
    double imsi_var11 = 0.0;
    double imsi_var12 = 0.0;

    double capitalized_interests = 0.0;
    vector<double> levelcvg_vec(num_couponleg_cf);

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1)
        {
            break;
        }
    }

    for (i = starti; i < num_fundingleg_cf; i++)
    {

        imsi_var0 = zc(fundingleg_paydate[i]);
        imsi_var1 = cvg(today, fundingleg_paydate[i], dcb);


        fundingleg_df[i] =
            exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));

        fra[i] = fixinghistory_rate[i];

        imsi_var2 = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);

        if (fra[i] > MinFixingRate)
        {
            floatinglegprice =
                floatinglegprice
                + fundingleg_notional[i]
                * (
                    fundingleg_mult[i] * fra[i]
                    + fundingleg_spread[i]
                    )
                * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
                * fundingleg_df[i];
        }

        if (today <= fixing_date[i])
        {
            break;
        }
    }

    if (today == fundingleg_paydate[num_fundingleg_cf - 1])
    {
        floatinglegprice = 0.0;
    }

    endi = i;

    if (endi < num_fundingleg_cf)
    {
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));

            // ondate, ondf 계산 hokim
        }
    }

    settlement_date_df = settlement_date_df / ondf;

    for (i = endi; i < num_fundingleg_cf; i++)
    {
        // [
        imsi_var3 = zc(fundingleg_paydate[i]);
        imsi_var4 = cvg(today, fundingleg_paydate[i], dcb);
        // ]

        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));

        // [
        imsi_var5 = zc(fundingleg_calc_startdate[i]);
        imsi_var6 = cvg(today, fundingleg_calc_startdate[i], dcb);

        imsi_var7 = zc(fixingindex_maturity[i]);
        imsi_var8 = cvg(today, fixingindex_maturity[i], dcb);

        imsi_var9 = cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        // ]

        fra[i] = (exp(-zc(fundingleg_calc_startdate[i]) * cvg(today, fundingleg_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        // [
        imsi_var10 = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        // ]

        floatinglegprice = floatinglegprice + (estflag)*fundingleg_notional[i] * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];

        //floatinglegprice=floatinglegprice+(estflag)*fundingleg_notional[i]*fundingleg_mult[i]*fra[i]*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];

        floatinglegprice = floatinglegprice + fundingleg_notional[i] * fundingleg_spread[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];

        //fout<<setprecision(18)<<fixingindex_maturity[i].get_str_date()<<"\t"<<fra[i]<<endl;
    }
    if (endi < num_fundingleg_cf)
    {
        // [
        imsi_var11 = zc(fundingleg_calc_startdate[endi]);
        imsi_var12 = cvg(today, fundingleg_calc_startdate[endi], dcb);

        // ]

        floatinglegprice = floatinglegprice + (1 - estflag) * (fundingleg_notional[endi] * exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)) - fundingleg_notional[num_fundingleg_cf - 1] * fundingleg_df[num_fundingleg_cf - 1]);
        //floatinglegprice = floatinglegprice + (1 - estflag)*(fundingleg_notional[endi] * exp(-zc(fundingleg_calc_startdate[endi])*cvg(today, fundingleg_calc_startdate[endi], dcb)));
        //floatinglegprice = floatinglegprice + (1 - estflag)*(fundingleg_notional[num_fundingleg_cf - 1] * fundingleg_df[num_fundingleg_cf - 1]);
    }
    // Start by hokim
    //floatinglegprice = floatinglegprice + fundingleg_notional[num_fundingleg_cf - 1] * fundingleg_df[num_fundingleg_cf - 1];
    // End by hokim
    starti = 0;

    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_couponleg_cf - 1) break;
    }

    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + levelcvg * couponleg_df[i];
        //coupon[i]=couponleg_couponrate[i]*levelcvg;
        //fixedlegprice=fixedlegprice+couponleg_notional[i]*coupon[i]*couponleg_df[i];
    }

    for (i = 0; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);

        //level=level+levelcvg*couponleg_df[i];
        coupon[i] = couponleg_couponrate[i] * levelcvg;

        //fixedlegprice=fixedlegprice+couponleg_notional[i]*coupon[i]*couponleg_df[i];
        fixedlegprice =
            fixedlegprice
            + couponleg_notional[i]
            * coupon[i];

        // MUREX ZERO COUPON INCLUDING COMPOUNDING hokim
        // [[[
        //if (i != num_couponleg_cf-1)
        //{
        //levelcvg_vec[i + 1] = cvg(couponleg_calc_startdate[i+1], couponleg_calc_enddate[i+1], couponleg_dcb);
        //capitalized_interests = fixedlegprice;
        //fixedlegprice = fixedlegprice*(1 + coupon[i] * levelcvg_vec[i + 1]);
        //}
        // ]]]
    }
    fixedlegprice = fixedlegprice * couponleg_df[num_couponleg_cf - 1];
    // Start by hokim
    //fixedlegprice = fixedlegprice + couponleg_notional[num_couponleg_cf - 1] * couponleg_df[num_couponleg_cf - 1];
    // End by hokim
    if (today == couponleg_paydate[num_couponleg_cf - 1])
    {
        fixedlegprice = 0.0;
    }

    if (starti < num_couponleg_cf)
    {
        atmswaprate = floatinglegprice / level / notional;
    }

    floatinglegprice =
        floatinglegprice
        / ondf
        / settlement_date_df;

    fixedlegprice =
        fixedlegprice
        / ondf
        / settlement_date_df;

    swapprice = (floatinglegprice - fixedlegprice);
}


// 2018-08-08 created by YJ
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
)
{
    ofstream foutcallprob("D:/HOKim_BackUp/000_Calc/Callprob/callprob.txt");

    string dcb = "ACT/365";

    double ondf = 1.0;

    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;

    double
        settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb))
        , level = 0.0
        , levelcvg = 0.0
        , couponleg_df = 0.0;

    int
        i
        , couponleg_idx
        , starti = 0
        , endi
        , num_couponleg_cf = int(couponleg_calc_startdate.size())
        , num_fundingleg_cf = int(fundingleg_calc_startdate.size());

    vector<double>
        fra(num_fundingleg_cf, 0.0)
        , fundingleg_df(num_fundingleg_cf, 1.0);     //couponleg_df(num_couponleg_cf,1.0), coupon(num_couponleg_cf),

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_fundingleg_cf - 1)
        {
            break;
        }
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];

        if (fra[i] > MinFixingRate)
        {
            floatinglegprice = floatinglegprice + fundingleg_notional[i] * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        }

        if (today <= fixing_date[i])
        {
            break;
        }
    }

    if (today == fundingleg_paydate[num_fundingleg_cf - 1])
    {
        floatinglegprice = 0.0;
    }

    endi = i;

    if (endi < num_fundingleg_cf)
    {
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));

        fra[i] =
            (
                exp
                (
                    -zc(fundingleg_calc_startdate[i])
                    * cvg(today, fundingleg_calc_startdate[i], dcb)
                    + zc(fixingindex_maturity[i])
                    * cvg(today, fixingindex_maturity[i], dcb)
                ) - 1.0
                )
            / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        floatinglegprice =
            floatinglegprice
            + (estflag)
            *fundingleg_notional[i]
            * (
                fundingleg_mult[i]
                * fra[i]
                + fundingleg_spread[i]
                )
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        floatinglegprice =
            floatinglegprice
            + fundingleg_notional[i]
            * fundingleg_spread[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];
    }

    if (endi < num_fundingleg_cf)
    {
        floatinglegprice =
            floatinglegprice
            + (1 - estflag)
            * (
                fundingleg_notional[endi]
                * exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)
                )
                - fundingleg_notional[num_fundingleg_cf - 1]
                * fundingleg_df[num_fundingleg_cf - 1]);
    }

    starti = 0;

    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_couponleg_cf - 1)
        {
            break;
        }
    }

    couponleg_idx = num_couponleg_cf - 1;

    couponleg_df = exp(-zc(couponleg_paydate[couponleg_idx]) * cvg(today, couponleg_paydate[couponleg_idx], dcb));

    levelcvg = cvg(couponleg_calc_startdate[0], couponleg_calc_enddate[couponleg_idx], couponleg_dcb);

    fixedlegprice =
        (
            pow
            (
            (
                1 + couponleg_couponrate[couponleg_idx]
                )
                , levelcvg
            )
            - 1
            )
        * couponleg_notional[couponleg_idx] * couponleg_df;

    if (today == couponleg_paydate[num_couponleg_cf - 1])
    {
        fixedlegprice = 0.0;
    }

    //if(starti<num_couponleg_cf) atmswaprate=floatinglegprice/level/notional;

    floatinglegprice =
        floatinglegprice
        / ondf
        / settlement_date_df;

    fixedlegprice =
        fixedlegprice
        / ondf
        / settlement_date_df;

    swapprice =
        (
            floatinglegprice
            - fixedlegprice
            );
}

// created by HoKim 2020-03-27
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
)
{
    ofstream foutcallprob("D:/HOKim_BackUp/000_Calc/Callprob/callprob.txt");

    string dcb = "ACT/365";

    double ondf = 1.0;

    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;

    // created by hokim 2020-04-13
    double floatinglegaccuralprice = 0.0;
    double flixedlegaccuralprice = 0.0;
    double accrualswapprice = 0.0;

    double
        settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb))
        , level = 0.0
        , levelcvg = 0.0
        , couponleg_df = 0.0;

    int
        i
        , couponleg_idx
        , starti = 0
        , endi
        , num_couponleg_cf = int(couponleg_calc_startdate.size())
        , num_fundingleg_cf = int(fundingleg_calc_startdate.size());

    vector<double>
        fra(num_fundingleg_cf, 0.0)
        , fundingleg_df(num_fundingleg_cf, 1.0);        //couponleg_df(num_couponleg_cf,1.0), coupon(num_couponleg_cf),

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_fundingleg_cf - 1)
        {
            break;
        }
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];

        if (fra[i] > MinFixingRate)
        {
            floatinglegprice
                = floatinglegprice
                * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i])
                * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
                + fundingleg_notional[i]
                * fundingleg_df[i];

            // created by hokim 2020-04-13
            floatinglegaccuralprice
                = floatinglegaccuralprice
                + fundingleg_notional[i]
                * (fundingleg_mult[i] * fra[i] + fundingleg_spread[i])
                * cvg(fundingleg_calc_startdate[i], today, fundingleg_dcb)
                * fundingleg_df[i];

        }

        if (today <= fixing_date[i])
        {
            break;
        }
    }

    if (today == fundingleg_paydate[num_fundingleg_cf - 1])
    {
        floatinglegprice = 0.0;
    }

    endi = i;

    if (endi < num_fundingleg_cf)
    {
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;

    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));

        fra[i] =
            (
                exp
                (
                    -zc(fundingleg_calc_startdate[i])
                    * cvg(today, fundingleg_calc_startdate[i], dcb)
                    + zc(fixingindex_maturity[i])
                    * cvg(today, fixingindex_maturity[i], dcb)
                ) - 1.0
                )
            / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);

        floatinglegprice =
            floatinglegprice
            + (estflag)
            *fundingleg_notional[i]
            * (
                fundingleg_mult[i]
                * fra[i]
                + fundingleg_spread[i]
                )
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];

        floatinglegprice =
            floatinglegprice
            + fundingleg_notional[i]
            * fundingleg_spread[i]
            * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
            * fundingleg_df[i];
    }

    if (endi < num_fundingleg_cf)
    {
        floatinglegprice =
            floatinglegprice
            + (1 - estflag)
            * (
                fundingleg_notional[endi]
                * exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)
                )
                - fundingleg_notional[num_fundingleg_cf - 1]
                * fundingleg_df[num_fundingleg_cf - 1]);
    }

    starti = 0;

    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;

        if (starti >= num_couponleg_cf - 1)
        {
            break;
        }
    }

    couponleg_idx = num_couponleg_cf - 1;

    couponleg_df = exp(-zc(couponleg_paydate[couponleg_idx]) * cvg(today, couponleg_paydate[couponleg_idx], dcb));

    levelcvg = cvg(couponleg_calc_startdate[0], couponleg_calc_enddate[couponleg_idx], couponleg_dcb);

    fixedlegprice =
        (
            pow
            (
            (
                1 + couponleg_couponrate[couponleg_idx]
                )
                , levelcvg
            )
            - 1
            )
        * couponleg_notional[couponleg_idx] * couponleg_df;

    flixedlegaccuralprice
        = couponleg_notional[starti]
        * couponleg_couponrate[starti]
        * cvg(couponleg_calc_startdate[starti], today, couponleg_dcb)
        * couponleg_df;

    if (today == couponleg_paydate[num_couponleg_cf - 1])
    {
        fixedlegprice = 0.0;
    }

    //if(starti<num_couponleg_cf) atmswaprate=floatinglegprice/level/notional;

    // Including principal at maturity By Hokim 2020-03-27
    //floatinglegprice =
    //  (floatinglegprice + fundingleg_df[num_fundingleg_cf - 1]* fundingleg_notional[num_fundingleg_cf - 1])
    //  / ondf
    //  / settlement_date_df;

    //fixedlegprice =
    //  (fixedlegprice + couponleg_notional[couponleg_idx] * couponleg_df)
    //  / ondf
    //  / settlement_date_df;

    swapprice =
        (
            floatinglegprice
            - fixedlegprice
            );
}

void PlainSwapPrice(
    CCurrency crcy, CDate today, CDate settlement_date, double notional, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double couponleg_couponrate, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, double fixinghistory_rate, bool estflag, double& floatinglegprice, double& fixedlegprice, double& atmswaprate, double& swapprice)
{
    //  ofstream fout("test.txt");
    string dcb = "ACT/365";
    double ondf = 1.0;
    fixedlegprice = 0.0;
    floatinglegprice = 0.0;
    atmswaprate = 0.0;
    double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0;
    int i, starti = 0, endi, num_couponleg_cf = int(couponleg_calc_startdate.size()), num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> couponleg_df(num_couponleg_cf, 1.0), coupon(num_couponleg_cf), fra(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0);

    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1) break;
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate;
        //if(fra[i]>0.0) floatinglegprice=floatinglegprice+fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];
        if (fra[i] > MinFixingRate) floatinglegprice = floatinglegprice + notional * (fra[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) floatinglegprice = 0.0;
    endi = i;
    if (endi < num_fundingleg_cf)
    {
        //if(today==fixing_date[i] && fra[i]>0.0)
        //if(today==fixing_date[i] && fra[i]>MinFixingRate)
        //if(today==fixing_date[endi] && fra[endi]>0)
        if (today == fixing_date[endi] && fra[endi] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    double tmp = 0;
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        tmp = zc(fixingindex_maturity[i]);
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = (exp(-zc(fundingleg_calc_startdate[i]) * cvg(today, fundingleg_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        //floatinglegprice=floatinglegprice+(estflag)*fundingleg_notional[i]*(fundingleg_mult[i]*fra[i]+fundingleg_spread[i])*cvg(fundingleg_calc_startdate[i],fundingleg_calc_enddate[i],fundingleg_dcb)*fundingleg_df[i];
        floatinglegprice = floatinglegprice + (estflag)*notional * fra[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        //fout<<setprecision(18)<<fixingindex_maturity[i].get_str_date()<<"\t"<<fra[i]<<endl;
    }
    if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + (1 - estflag) * (notional * exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)) - notional * fundingleg_df[num_fundingleg_cf - 1]);
    starti = 0;
    while (today >= couponleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_couponleg_cf - 1) break;
    }
    for (i = starti; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + levelcvg * couponleg_df[i];
        coupon[i] = couponleg_couponrate * levelcvg;
        fixedlegprice = fixedlegprice + notional * coupon[i] * couponleg_df[i];
    }
    if (today == couponleg_paydate[num_couponleg_cf - 1]) fixedlegprice = 0.0;
    if (starti < num_couponleg_cf) atmswaprate = floatinglegprice / level / notional;
    floatinglegprice = floatinglegprice / ondf / settlement_date_df;
    fixedlegprice = fixedlegprice / ondf / settlement_date_df;
    swapprice = (floatinglegprice - fixedlegprice);
}


void bondpricezc(CDate today, CDate spotdate, CInterpolation zc, double coupon, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double& bondd, double& bondc)
{
    int j, k = -1, k0;
    double tmpdf, accrued, spotdf = exp(-zc(spotdate) * cvg(today, spotdate, "ACT/365"));
    bondd = 0.0;
    for (j = 0; j < num_bondcashflowschedule; j++)
    {
        if (bond_schedule[0][j] >= today)
        {
            k = k + 1;
            //tmpdf=pow(1.0/(1.0+ytm/ncpy),k);
            //bondvalue=bondvalue+tmpdf*coupon;
            tmpdf = exp(-zc(bond_schedule[0][j]) * cvg(today, bond_schedule[0][j], "ACT/365"));
            bondd = bondd + tmpdf * coupon / ncpy;
        }
    }
    k0 = num_bondcashflowschedule - k - 1;
    bondd = bondd + tmpdf;
    bondd /= spotdf;

    //if(k0==0) accrued=coupon/ncpy*(bond_schedule[2][k0]-spotdate)/(bond_schedule[2][k0]-bond_schedule[1][k0]);
    //else accrued=coupon/ncpy*(bond_schedule[2][k0]-spotdate)/(bond_schedule[2][k0]-bond_schedule[2][k0-1]);
    if (k0 == 0) accrued = coupon / ncpy * (spotdate - bond_schedule[1][k0]) / (bond_schedule[2][k0] - bond_schedule[1][k0]);
    else accrued = coupon / ncpy * (spotdate - bond_schedule[2][k0 - 1]) / (bond_schedule[2][k0] - bond_schedule[2][k0 - 1]);
    bondc = bondd - accrued;//<-------------여기

}

double findbondytmzc(CDate today, CDate spotdate, CInterpolation zc, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double bondp)
{
    double left = 0.0, right = 100.0, ytm = 0.5 * (left + right);
    double ans, bondd, bondc;
    bondpricezc(today, spotdate, zc, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
    ans = bondc - bondp;
    while (fabs(ans) > Error&& fabs(ytm - left) > Error&& fabs(ytm - right) > Error)
    {
        if (ans < 0.0)
        {
            left = ytm;
            ytm = 0.5 * (left + right);
            bondpricezc(today, spotdate, zc, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
            ans = bondc - bondp;
        }
        else if (ans > 0.0)
        {
            right = ytm;
            ytm = 0.5 * (left + right);
            bondpricezc(today, spotdate, zc, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
            ans = bondc - bondp;
        }
        else break;
    }
    return ytm;
}

void bondpriceytm(CDate today, CDate settledate, double coupon, double ytm, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double& bondd, double& bondc)
{
    int j, k = -1, k0;
    double tmpdf, accrued, accrueddf;//, spotdf=exp(-zc(spotdate)*cvg(today,spotdate,"ACT/365"))
    bondd = 0.0;
    for (j = 0; j < num_bondcashflowschedule; j++)
    {
        if (bond_schedule[0][j] >= settledate)
        {
            k = k + 1;
            tmpdf = pow(1.0 / (1.0 + ytm / ncpy), k);
            bondd = bondd + tmpdf * coupon / ncpy;
        }
    }
    k0 = num_bondcashflowschedule - k - 1;
    bondd = bondd + tmpdf;

    if (k0 == 0)
    {
        accrueddf = 1.0 / (1.0 + ytm / ncpy * (bond_schedule[2][k0] - settledate) / (bond_schedule[2][k0] - bond_schedule[1][k0]));
        accrued = coupon / ncpy * (settledate - bond_schedule[1][k0]) / (bond_schedule[2][k0] - bond_schedule[1][k0]);
    }
    else
    {
        accrueddf = 1.0 / (1.0 + ytm / ncpy * (bond_schedule[2][k0] - settledate) / (bond_schedule[2][k0] - bond_schedule[2][k0 - 1]));
        accrued = coupon / ncpy * (settledate - bond_schedule[2][k0 - 1]) / (bond_schedule[2][k0] - bond_schedule[2][k0 - 1]);
    }
    bondd *= accrueddf;//<-------------여기
    bondc = bondd - accrued;//<-------------여기
}

double findbondytm(CDate today, CDate spotdate, CDate settledate, double coupon, double ncpy, int num_bondcashflowschedule, vector<vector<CDate>> bond_schedule, double bondp)
{
    double left = 0.0, right = 100.0, ytm = 0.5 * (left + right);
    double ans, bondd, bondc;
    bondpriceytm(today, settledate, coupon, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
    ans = bondc - bondp;
    while (fabs(ans) > Error&& fabs(ytm - left) > Error&& fabs(ytm - right) > Error)
    {
        if (ans > 0.0)
        {
            left = ytm;
            ytm = 0.5 * (left + right);
            bondpriceytm(today, settledate, coupon, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
            ans = bondc - bondp;
        }
        else if (ans < 0.0)
        {
            right = ytm;
            ytm = 0.5 * (left + right);
            bondpriceytm(today, settledate, coupon, ytm, ncpy, num_bondcashflowschedule, bond_schedule, bondd, bondc);
            ans = bondc - bondp;
        }
        else break;
    }
    return ytm;
}

void CapFloorForwardVolMatrix(CCurrency crcy, CDate today, CDate* indexholiday, int num_indexholiday, CInterpolation zc, int num_strike, double centervol, double volinterval, vector<string> tenor, vector<string> freq, vector<double> capatmvol, vector<double>& capstrike, vector<CDate>& cpltmtrty, vector<vector<double>>& forwardvol)
{
    CInstrument cap(crcy, "CAP");
    CDate* holiday = cap.get_holiday(), * fixingholiday = cap.get_fixingholiday();
    int numholiday = cap.get_numholiday(), numfixingholiday = cap.get_numfixingholiday();
    string capstub = cap.get_stub(), capdirection = cap.get_direction(), capconv = cap.get_conv(), capadjflag = cap.get_adjflag(), cappayin = cap.get_payin(), capsetin = cap.get_setin(), capdcb = cap.get_basis(), dcb = "ACT/365";
    int i, j, k, num_capkeytenor = int(tenor.size()), spotlag = cap.get_spotlag(), fixlag = -spotlag;
    CDate spotdate, tmpdate;
    vector<CDate> maturity(num_capkeytenor);
    vector<int> num_caplet(num_capkeytenor, 0), num_prevcaplet(num_capkeytenor, 0);
    vector<double> fra, df, tau, T;
    ShiftBusDate(today, holiday, numholiday, spotlag, spotdate);
    capstrike[num_strike] = centervol;

    for (i = 0; i < num_strike; i++)
    {
        capstrike[num_strike + i + 1] = centervol + (i + 1) * volinterval;
        capstrike[num_strike - (i + 1)] = centervol - (i + 1) * volinterval;
    }
    for (i = 0; i < num_capkeytenor; i++)
    {
        num_prevcaplet[i] = int(cpltmtrty.size());
        num_caplet[i] = findnumschedule(spotdate, tenor[i], capstub, capdirection, freq[i]);
        vector<CDate> calcstartdate(num_caplet[i]), calcenddate(num_caplet[i]), paydate(num_caplet[i]), fixingdate(num_caplet[i]), fraenddate(num_caplet[i]), fracalcstartdate(num_caplet[i]);//, fracalcenddate(num_caplet[i])
        calculationstartschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], calcstartdate);
        calculationendschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], calcenddate);

        calculationstartschedule(spotdate, tenor[i], indexholiday, num_indexholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], fracalcstartdate);
        //calculationendschedule(spotdate,tenor[i],indexholiday,num_indexholiday,capstub,capdirection,capconv,freq[i],capadjflag,num_caplet[i],fracalcenddate);

        paymentschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], num_caplet[i], paydate);
        fixingschedule(spotdate, tenor[i], holiday, numholiday, fixingholiday, numfixingholiday, capstub, capdirection, capconv, freq[i], capadjflag, capsetin, fixlag, num_caplet[i], fixingdate);
        for (j = num_prevcaplet[i]; j < num_caplet[i]; j++)
        {
            ShiftBusDate(fixingdate[j], fixingholiday, numfixingholiday, spotlag, tmpdate);
            findmaturity(tmpdate, freq[i], holiday, numholiday, capconv, tmpdate);
            df.push_back(exp(-zc(calcenddate[j]) * cvg(today, calcenddate[j], dcb)));
            //fra.push_back((exp(-zc(calcstartdate[j])*cvg(today,calcstartdate[j],dcb)+zc(tmpdate)*cvg(today,tmpdate,dcb))-1.0)/cvg(calcstartdate[j],tmpdate,capdcb));
            fra.push_back((exp(-zc(fracalcstartdate[j]) * cvg(today, fracalcstartdate[j], dcb) + zc(tmpdate) * cvg(today, tmpdate, dcb)) - 1.0) / cvg(fracalcstartdate[j], tmpdate, capdcb));
            tau.push_back(cvg(calcstartdate[j], calcenddate[j], capdcb));
            T.push_back(cvg(today, fixingdate[j], dcb));
            cpltmtrty.push_back(fixingdate[j]);
        }
    }
    forwardvol = vector<vector<double>>(num_caplet[num_capkeytenor - 1]);
    for (i = 0; i < num_caplet[num_capkeytenor - 1]; i++) forwardvol[i] = vector<double>(2 * num_strike + 1);
    for (i = 0; i < num_caplet[0]; i++) for (j = -num_strike; j <= num_strike; j++) forwardvol[i][num_strike + j] = capatmvol[0];
    double prevcap, fullcap, fwdvol0, fwdvol1;
    for (j = -num_strike; j <= num_strike; j++)
    {
        for (i = 1; i < num_capkeytenor; i++)
        {
            prevcap = 0.0;
            fullcap = 0.0;
            for (k = 1; k < num_prevcaplet[i]; k++) prevcap = prevcap + Caplet(fra[k], capstrike[num_strike + j], df[k], tau[k], T[k], forwardvol[k][num_strike + j]);
            for (k = 1; k < num_caplet[i]; k++) fullcap = fullcap + Caplet(fra[k], capstrike[num_strike + j], df[k], tau[k], T[k], capatmvol[i]);
            fwdvol0 = forwardvol[num_prevcaplet[i] - 1][num_strike + j];
            int tmpnumcaplet = num_caplet[i] - num_prevcaplet[i] + 1;
            vector<double> tmpfra(tmpnumcaplet), tmpdf(tmpnumcaplet), tmpcapletperiod(tmpnumcaplet), tmpcapletmaturity(tmpnumcaplet);
            for (k = 0; k < tmpnumcaplet; k++)
            {
                tmpfra[k] = fra[num_prevcaplet[i] - 1 + k];
                tmpdf[k] = df[num_prevcaplet[i] - 1 + k];
                tmpcapletperiod[k] = tau[num_prevcaplet[i] - 1 + k];
                tmpcapletmaturity[k] = T[num_prevcaplet[i] - 1 + k];
            }
            fwdvol1 = Findfwdvol(fwdvol0, fullcap, prevcap, capstrike[num_strike + j], tmpdf, tmpcapletmaturity, tmpcapletperiod, tmpfra);
            forwardvol[num_caplet[i] - 1][num_strike + j] = fwdvol1;
            for (k = num_prevcaplet[i]; k < num_caplet[i]; k++) forwardvol[k][num_strike + j] = (fwdvol1 - fwdvol0) / (T[num_caplet[i] - 1] - T[num_prevcaplet[i] - 1]) * (T[k] - T[num_prevcaplet[i] - 1]) + fwdvol0;
            fwdvol0 = fwdvol1;
        }
    }
}

void CapFloorForwardVolMatrix(CCurrency crcy, CDate today, CDate* indexholiday, int num_indexholiday, CInterpolation zc, vector<string> tenor, vector<string> freq, vector<vector<double>> capvolskew, vector<double> capstrike, vector<CDate>& cpltmtrty, vector<vector<double>>& forwardvol)
{
    CInstrument cap(crcy, "CAP");
    CDate* holiday = cap.get_holiday(), * fixingholiday = cap.get_fixingholiday();
    int numholiday = cap.get_numholiday(), numfixingholiday = cap.get_numfixingholiday();
    string capstub = cap.get_stub(), capdirection = cap.get_direction(), capconv = cap.get_conv(), capadjflag = cap.get_adjflag(), cappayin = cap.get_payin(), capsetin = cap.get_setin(), capdcb = cap.get_basis(), dcb = "ACT/365";
    int i, j, k, num_capkeytenor = int(tenor.size()), spotlag = cap.get_spotlag(), fixlag = -spotlag, num_strike = int(capstrike.size());
    CDate spotdate, tmpdate;
    vector<CDate> maturity(num_capkeytenor);
    vector<int> num_caplet(num_capkeytenor, 0), num_prevcaplet(num_capkeytenor, 0);
    vector<double> fra, df, tau, T;
    ShiftBusDate(today, holiday, numholiday, spotlag, spotdate);

    //ofstream foutschedule("schedule.txt");
    //foutschedule<<"cal_start\t"<<"cal_end\t"<<"fixing\t"<<"indexmat\t"<<"fra\t"<<"pay\t"<<"df\t"<<"tau"<<endl;

    for (i = 0; i < num_capkeytenor; i++)
    {
        num_prevcaplet[i] = int(cpltmtrty.size());
        num_caplet[i] = findnumschedule(spotdate, tenor[i], capstub, capdirection, freq[i]);
        vector<CDate> calcstartdate(num_caplet[i]), calcenddate(num_caplet[i]), paydate(num_caplet[i]), fixingdate(num_caplet[i]), fraenddate(num_caplet[i]), fracalcstartdate(num_caplet[i]);//, fracalcenddate(num_caplet[i])
        calculationstartschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], calcstartdate);
        calculationendschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], calcenddate);

        calculationstartschedule(spotdate, tenor[i], indexholiday, num_indexholiday, capstub, capdirection, capconv, freq[i], capadjflag, num_caplet[i], fracalcstartdate);
        //calculationendschedule(spotdate,tenor[i],indexholiday,num_indexholiday,capstub,capdirection,capconv,freq[i],capadjflag,num_caplet[i],fracalcenddate);

        paymentschedule(spotdate, tenor[i], holiday, numholiday, capstub, capdirection, capconv, freq[i], num_caplet[i], paydate);
        fixingschedule(spotdate, tenor[i], holiday, numholiday, fixingholiday, numfixingholiday, capstub, capdirection, capconv, freq[i], capadjflag, capsetin, fixlag, num_caplet[i], fixingdate);
        for (j = num_prevcaplet[i]; j < num_caplet[i]; j++)
        {
            ShiftBusDate(fixingdate[j], fixingholiday, numfixingholiday, spotlag, tmpdate);
            findmaturity(tmpdate, freq[i], holiday, numholiday, capconv, tmpdate);
            df.push_back(exp(-zc(calcenddate[j]) * cvg(today, calcenddate[j], dcb)));
            //fra.push_back((exp(-zc(calcstartdate[j])*cvg(today,calcstartdate[j],dcb)+zc(tmpdate)*cvg(today,tmpdate,dcb))-1.0)/cvg(calcstartdate[j],tmpdate,capdcb));
            fra.push_back((exp(-zc(fracalcstartdate[j]) * cvg(today, fracalcstartdate[j], dcb) + zc(tmpdate) * cvg(today, tmpdate, dcb)) - 1.0) / cvg(fracalcstartdate[j], tmpdate, capdcb));
            tau.push_back(cvg(calcstartdate[j], calcenddate[j], capdcb));
            T.push_back(cvg(today, fixingdate[j], dcb));
            cpltmtrty.push_back(fixingdate[j]);
            //foutschedule<<calcstartdate[j].get_str_date()<<"\t"<<calcenddate[j].get_str_date()<<"\t"<<fixingdate[j].get_str_date()<<"\t"<<tmpdate.get_str_date()<<"\t"<<setprecision(12)<<fra[j]<<"\t"<<paydate[j].get_str_date()<<"\t"<<setprecision(12)<<df[j]<<"\t"<<setprecision(12)<<tau[j]<<endl;
        }
    }
    forwardvol = vector<vector<double>>(num_caplet[num_capkeytenor - 1]);
    for (i = 0; i < num_caplet[num_capkeytenor - 1]; i++) forwardvol[i] = vector<double>(num_strike);
    for (i = 0; i < num_caplet[0]; i++) for (j = 0; j < num_strike; j++) forwardvol[i][j] = capvolskew[0][j];
    double prevcap, fullcap, fwdvol0, fwdvol1;
    for (j = 0; j < num_strike; j++)
    {
        for (i = 1; i < num_capkeytenor; i++)
        {
            prevcap = 0.0;
            fullcap = 0.0;
            for (k = 1; k < num_prevcaplet[i]; k++) prevcap = prevcap + Caplet(fra[k], capstrike[j], df[k], tau[k], T[k], forwardvol[k][j]);
            for (k = 1; k < num_caplet[i]; k++) fullcap = fullcap + Caplet(fra[k], capstrike[j], df[k], tau[k], T[k], capvolskew[i][j]);
            fwdvol0 = forwardvol[num_prevcaplet[i] - 1][j];
            int tmpnumcaplet = num_caplet[i] - num_prevcaplet[i] + 1;
            vector<double> tmpfra(tmpnumcaplet), tmpdf(tmpnumcaplet), tmpcapletperiod(tmpnumcaplet), tmpcapletmaturity(tmpnumcaplet);
            for (k = 0; k < tmpnumcaplet; k++)
            {
                tmpfra[k] = fra[num_prevcaplet[i] - 1 + k];
                tmpdf[k] = df[num_prevcaplet[i] - 1 + k];
                tmpcapletperiod[k] = tau[num_prevcaplet[i] - 1 + k];
                tmpcapletmaturity[k] = T[num_prevcaplet[i] - 1 + k];
            }
            fwdvol1 = Findfwdvol(fwdvol0, fullcap, prevcap, capstrike[j], tmpdf, tmpcapletmaturity, tmpcapletperiod, tmpfra);
            forwardvol[num_caplet[i] - 1][j] = fwdvol1;
            for (k = num_prevcaplet[i]; k < num_caplet[i]; k++) forwardvol[k][j] = (fwdvol1 - fwdvol0) / (T[num_caplet[i] - 1] - T[num_prevcaplet[i] - 1]) * (T[k] - T[num_prevcaplet[i] - 1]) + fwdvol0;
            fwdvol0 = fwdvol1;
        }
    }
}

void PlainCapFloorDigitalPrice(CCurrency crcy, CDate today, CDate settlement_date, double notional, double strike, double inputvol, string callputflag, bool firstcashflowincludeflag, double firstfixingrate, CInterpolation zc, CInterpolation vol, vector<double> fundingleg_notional, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, vector<CDate> fundinglegindex_calc_startdate, string fundingleg_dcb, vector<double> fundingleg_mult, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, double& atmrate, double& flatvol, double& price)
{
    string dcb = "ACT/365";
    price = 0.0;
    double ondf = 1.0;
    double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb));
    int i, starti = 0, endi, num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> fra(num_fundingleg_cf, 0.0), fwd1(num_fundingleg_cf, 0.0), fwd2(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0), tau(num_fundingleg_cf), T(num_fundingleg_cf);
    vector<int> sgn_1stcfincludeflag(num_fundingleg_cf);
    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1) break;
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];
        fwd1[i] = fundingleg_mult[i] * (sign(i) * fra[i] + (1 - sign(i)) * ((1.0 - sign(firstfixingrate)) * fra[i] + sign(firstfixingrate) * firstfixingrate));
        tau[i] = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        T[i] = cvg(today, fixing_date[i], dcb);
        sgn_1stcfincludeflag[i] = sign(i + firstcashflowincludeflag);
        price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd1[i], strike, fundingleg_df[i], tau[i], T[i], sign(inputvol) * inputvol + (1.0 - sign(inputvol)) * vol(fixing_date[i], strike));
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) price = 0.0;
    endi = i;
    if (endi < num_fundingleg_cf)
    {
        //if(today==fixing_date[i] && fra[i]>0.0)
        if (today == fixing_date[i] && fra[i] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        //fra[i]=(exp(-zc(fundingleg_calc_startdate[i])*cvg(today,fundingleg_calc_startdate[i],dcb)+zc(fixingindex_maturity[i])*cvg(today,fixingindex_maturity[i],dcb))-1.0)/cvg(fundingleg_calc_startdate[i],fixingindex_maturity[i],fixingindex_dcb);
        fra[i] = (exp(-zc(fundinglegindex_calc_startdate[i]) * cvg(today, fundinglegindex_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundinglegindex_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        fwd2[i] = fundingleg_mult[i] * (sign(i) * fra[i] + (1 - sign(i)) * ((1.0 - sign(firstfixingrate)) * fra[i] + sign(firstfixingrate) * firstfixingrate));
        tau[i] = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        T[i] = cvg(today, fixing_date[i], dcb);
        sgn_1stcfincludeflag[i] = sign(i + firstcashflowincludeflag);
        price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], sign(inputvol) * inputvol + (1.0 - sign(inputvol)) * vol(fixing_date[i], strike));
    }

    vector<double> couponrate(num_fundingleg_cf, 0.0), mult(num_fundingleg_cf, 1.0), spread(num_fundingleg_cf, 0.0);
    double floatinglegprice, fixedlegprice, swapprice;
    PlainSwapPrice(crcy, today, settlement_date, notional, zc, fundingleg_calc_startdate, fundingleg_calc_enddate, fundingleg_paydate, fundingleg_dcb, fundingleg_notional, couponrate, fundingleg_calc_startdate, fundingleg_calc_enddate, fundingleg_paydate, fundingleg_dcb, fundingleg_notional, mult, spread, fixing_date, fixingindex_maturity, fixingindex_dcb, fixinghistory_date, fixinghistory_rate, true, floatinglegprice, fixedlegprice, atmrate, swapprice);
    FindATMVol(starti, endi, num_fundingleg_cf, atmrate, price, today, strike, callputflag, fundingleg_notional, fixing_date, fundingleg_paydate, fundingleg_df, fwd1, fwd2, tau, T, sgn_1stcfincludeflag, flatvol);

    price = price / ondf / settlement_date_df;

}

void PlainCapFloorDigitalPrice2(CCurrency crcy, CDate today, CDate settlement_date, double notional, double strike, double inputvol, string callputflag, bool firstcashflowincludeflag, double firstfixingrate, CInterpolation zc, CInterpolation vol, vector<double> fundingleg_notional, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, vector<CDate> fundinglegindex_calc_startdate, string fundingleg_dcb, vector<double> fundingleg_mult, vector<CDate> fixing_date, vector<CDate> fixingindex_calc_startdate, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, double& atmrate, double& flatvol, double& price)
{
    string dcb = "ACT/365";
    price = 0.0;
    double ondf = 1.0;
    double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb));
    int i, starti = 0, endi, num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> fra(num_fundingleg_cf, 0.0), fwd1(num_fundingleg_cf, 0.0), fwd2(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0), tau(num_fundingleg_cf), T(num_fundingleg_cf);
    vector<int> sgn_1stcfincludeflag(num_fundingleg_cf);
    while (today >= fundingleg_paydate[starti])
    {
        starti = starti + 1;
        if (starti >= num_fundingleg_cf - 1) break;
    }
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = fixinghistory_rate[i];
        fwd1[i] = fundingleg_mult[i] * (sign(i) * fra[i] + (1 - sign(i)) * ((1.0 - sign(firstfixingrate)) * fra[i] + sign(firstfixingrate) * firstfixingrate));
        tau[i] = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        T[i] = cvg(today, fixing_date[i], dcb);
        sgn_1stcfincludeflag[i] = sign(i + firstcashflowincludeflag);
        if (T[i] <= 0.0)
        {
            if (callputflag == "CALL") price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * max(fwd1[i] - strike, 0.0) * tau[i] * fundingleg_df[i];
            else if (callputflag == "PUT") price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * max(strike - fwd1[i], 0.0) * tau[i] * fundingleg_df[i];
            else if (callputflag == "DCALL")
            {
                if (fwd1[i] >= strike) price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * tau[i] * fundingleg_df[i];
            }
            else
            {
                if (fwd1[i] <= strike) price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * tau[i] * fundingleg_df[i];
            }

        }
        else { if (fra[i] > MinFixingRate) price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd1[i], strike, fundingleg_df[i], tau[i], T[i], sign(inputvol) * inputvol + (1.0 - sign(inputvol)) * vol(fixing_date[i], strike)); }
        if (today <= fixing_date[i])
        {
            if (fra[i] <= MinFixingRate) price = 0.0;
            break;
        }
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) price = 0.0;
    endi = i;
    if (endi < num_fundingleg_cf)
    {
        //if(today==fixing_date[i] && fra[i]>0.0)
        if (today == fixing_date[i] && fra[i] > MinFixingRate)
        {
            endi = i + 1;
            CDate ondate;
            CInstrument temp(crcy, "DEPO");
            ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
            ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
        }
    }
    settlement_date_df = settlement_date_df / ondf;
    for (i = endi; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        //fra[i]=(exp(-zc(fundingleg_calc_startdate[i])*cvg(today,fundingleg_calc_startdate[i],dcb)+zc(fixingindex_maturity[i])*cvg(today,fixingindex_maturity[i],dcb))-1.0)/cvg(fundingleg_calc_startdate[i],fixingindex_maturity[i],fixingindex_dcb);
        fra[i] = (exp(-zc(fixingindex_calc_startdate[i]) * cvg(today, fixingindex_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fixingindex_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        fwd2[i] = fundingleg_mult[i] * (sign(i) * fra[i] + (1 - sign(i)) * ((1.0 - sign(firstfixingrate)) * fra[i] + sign(firstfixingrate) * firstfixingrate));
        tau[i] = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        T[i] = cvg(today, fixing_date[i], dcb);
        sgn_1stcfincludeflag[i] = sign(i + firstcashflowincludeflag);
        price = price + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], sign(inputvol) * inputvol + (1.0 - sign(inputvol)) * vol(fixing_date[i], strike));
    }

    vector<double> couponrate(num_fundingleg_cf, 0.0), mult(num_fundingleg_cf, 1.0), spread(num_fundingleg_cf, 0.0);
    double floatinglegprice, fixedlegprice, swapprice;
    PlainSwapPrice(crcy, today, settlement_date, notional, zc, fundingleg_calc_startdate, fundingleg_calc_enddate, fundingleg_paydate, fundingleg_dcb, fundingleg_notional, couponrate, fundingleg_calc_startdate, fundingleg_calc_enddate, fundingleg_paydate, fundingleg_dcb, fundingleg_notional, mult, spread, fixing_date, fixingindex_maturity, fixingindex_dcb, fixinghistory_date, fixinghistory_rate, true, floatinglegprice, fixedlegprice, atmrate, swapprice);
    FindATMVol(starti, endi, num_fundingleg_cf, atmrate, price, today, strike, callputflag, fundingleg_notional, fixing_date, fundingleg_paydate, fundingleg_df, fwd1, fwd2, tau, T, sgn_1stcfincludeflag, flatvol);

    price = price / ondf / settlement_date_df;

}

void CapPrice_Bl(CCurrency crcy, CDate today, CDate settlement_date, CInterpolation zc, CInterpolation vol, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<CDate> fixing_date, vector<CDate> fixingindex_calc_startdate, vector<CDate> fixingindex_maturity, string fixingindex_dcb, double atmrate, double& price)
{
    string dcb = "ACT/365";
    price = 0.0;
    double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb));
    int i, num_fundingleg_cf = int(fundingleg_calc_startdate.size());
    vector<double> fra(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0), tau(num_fundingleg_cf), T(num_fundingleg_cf);
    for (i = 0; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        fra[i] = (exp(-zc(fixingindex_calc_startdate[i]) * cvg(today, fixingindex_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fixingindex_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
        tau[i] = cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb);
        T[i] = cvg(today, fixing_date[i], dcb);
        price = price + PlainOptionlet("CALL", fra[i], atmrate, fundingleg_df[i], tau[i], T[i], vol(fixing_date[i], atmrate));
    }
    price = price / settlement_date_df;
}


void ParSwapRate(CDate today, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double& atmswaprate)
{
    string dcb = "ACT/365";
    atmswaprate = 0.0;
    double floatinglegprice = 0.0, level = 0.0, levelcvg = 0.0;
    int i, num_couponleg_cf = int(couponleg_calc_startdate.size());
    vector<double> couponleg_df(num_couponleg_cf, 1.0);

    for (i = 0; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + levelcvg * couponleg_df[i];
    }
    floatinglegprice = (exp(-zc(couponleg_calc_startdate[0]) * cvg(today, couponleg_calc_startdate[0], dcb)) - couponleg_df[num_couponleg_cf - 1]);
    atmswaprate = floatinglegprice / level;
}


void FindATMVol(int starti, int endi, int num_fundingleg_cf, double atmswaprate, double price, CDate today, double strike, string callputflag, vector <double> fundingleg_notional, vector<CDate> fixing_date, vector<CDate> fundingleg_paydate, vector<double> fundingleg_df, vector<double> fwd1, vector<double> fwd2, vector<double> tau, vector<double> T, vector<int> sgn_1stcfincludeflag, double& flatvol)
{
    int i = 0, maxnum = 5000, count = 0;
    double ans, tmpprice0 = 0.0, tmpprice = 0.0, left = 0.0, right = 10.0;
    flatvol = (left + right) * 0.5;
    for (i = starti; i < num_fundingleg_cf; i++)
    {
        tmpprice0 = tmpprice0 + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd1[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) tmpprice0 = 0.0;
    tmpprice = tmpprice0;
    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
    ans = tmpprice - price;
    if (callputflag == "DCALL")
    {
        if (atmswaprate >= strike)
        {
            while (fabs(ans) > Error&& fabs(flatvol - left) > Error&& fabs(flatvol - right) > Error&& count < maxnum)
            {
                if (ans > 0.0)
                {
                    left = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
                    ans = tmpprice - price;
                }
                else if (ans < 0.0)
                {
                    right = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
                    ans = tmpprice - price;
                }
                else break;
                count = count + 1;
            }
        }
        else
        {
            while (fabs(ans) > Error&& fabs(flatvol - left) > Error&& fabs(flatvol - right) > Error&& count < maxnum)
            {
                if (ans < 0.0)
                {
                    left = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                    ans = tmpprice - price;
                }
                else if (ans > 0.0)
                {
                    right = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                    ans = tmpprice - price;
                }
                else break;
                count = count + 1;
            }
        }
    }
    else if (callputflag == "DPUT")
    {
        if (atmswaprate >= strike)
        {
            while (fabs(ans) > Error&& fabs(flatvol - left) > Error&& fabs(flatvol - right) > Error&& count < maxnum)
            {
                if (ans < 0.0)
                {
                    left = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
                    ans = tmpprice - price;
                }
                else if (ans > 0.0)
                {
                    right = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag, fwd2[i], strike, fundingleg_df[i], tau[i], T[i], flatvol);
                    ans = tmpprice - price;
                }
                else break;

                count = count + 1;
            }
        }
        else
        {
            while (fabs(ans) > Error&& fabs(flatvol - left) > Error&& fabs(flatvol - right) > Error&& count < maxnum)
            {
                if (ans > 0.0)
                {
                    left = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                    ans = tmpprice - price;
                }
                else if (ans < 0.0)
                {
                    right = flatvol;
                    flatvol = 0.5 * (left + right);
                    tmpprice = tmpprice0;
                    for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                    ans = tmpprice - price;
                }
                else break;
                count = count + 1;
            }
        }
    }
    else
    {
        while (fabs(ans) > Error&& fabs(flatvol - left) > Error&& fabs(flatvol - right) > Error)
        {
            if (ans < 0.0)
            {
                left = flatvol;
                flatvol = 0.5 * (left + right);
                tmpprice = tmpprice0;
                for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                ans = tmpprice - price;
            }
            else if (ans > 0.0)
            {
                right = flatvol;
                flatvol = 0.5 * (left + right);
                tmpprice = tmpprice0;
                for (i = endi; i < num_fundingleg_cf; i++) tmpprice = tmpprice + fundingleg_notional[i] * sgn_1stcfincludeflag[i] * PlainOptionlet(callputflag,fwd2[i],strike,fundingleg_df[i],tau[i],T[i],flatvol);
                ans = tmpprice - price;
            }
            else break;
        }
    }
}


double Bl(double K, double F, double v, double w)
{
    return F * w * cumnormal_Pol_App(w * (log(F / K) + v * v * 0.5) / v) - K * w * cumnormal_Pol_App(w * (log(F / K) - v * v * 0.5) / v);
}

double Caplet(double fwd, double k, double df, double tau, double t, double vol)
{
    double ans = 0.0;
    if (fwd > 0.0)
    {
        if (t > 0.0) ans = df * tau * Bl(k, fwd, vol * sqrt(t), 1);
        else if (t >= -tau) ans = max(fwd - k, 0.0) * tau * df;
    }
    return ans;
}

double Floorlet(double fwd, double k, double df, double tau, double t, double vol)
{
    double ans = 0.0;
    if (fwd > 0.0)
    {
        if (t > 0.0) ans = df * tau * Bl(k, fwd, vol * sqrt(t), -1);
        else if (t >= -tau) ans = max(k - fwd, 0.0) * tau * df;
    }
    else ans = max(k - fwd, 0.0) * tau * df;
    return ans;
}

double DigitalCaplet(double fwd, double k, double df, double tau, double t, double vol)
{
    double ans = 0.0, sig = vol * sqrt(t), notional = 0.01;
    if (fwd > 0.0)
    {
        if (t > 0.0) ans = notional * df * tau * cumnormal((log(fwd / k) - sig * sig * 0.5) / sig);
        else if (t >= -tau && fwd >= k) ans = notional * tau * df;
    }
    return ans;
}

double DigitalFloorlet(double fwd, double k, double df, double tau, double t, double vol)
{
    double ans = 0.0, sig = vol * sqrt(t), notional = 0.01;
    if (fwd > 0.0)
    {
        if (t > 0.0) ans = notional * df * tau * cumnormal((-log(fwd / k) + sig * sig * 0.5) / sig);
        else if (t >= -tau && fwd <= k) ans = notional * tau * df;
    }
    else ans = notional * tau * df;
    return ans;
}

double PlainOptionlet(string callputflag, double fwd, double k, double df, double tau, double t, double vol)
{
    double ans = 0.0;
    if (callputflag == "CALL") ans = Caplet(fwd, k, df, tau, t, vol);
    else if (callputflag == "PUT") ans = Floorlet(fwd, k, df, tau, t, vol);
    else if (callputflag == "DCALL") ans = DigitalCaplet(fwd, k, df, tau, t, vol);
    else ans = DigitalFloorlet(fwd, k, df, tau, t, vol);
    return ans;
}

double DigitalSpreadCall(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k)
{
    int i;
    double sum = 0.0, notional = 1.0, M1 = log(f1) - 0.5 * sig1 * sig1 * t, M2 = log(f2) - 0.5 * sig2 * sig2 * t, S1 = sig1 * sqrt(t), S2 = sig2 * sqrt(t);
    if (k >= 0.0) for (i = 0; i < NumGLQ; i++) sum = sum + GLQw[i] * cumnormal_Pol_App(((y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i]) - log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k)) / (S2 * sqrt(1.0 - corr * corr)));
    else for (i = 0; i < NumGLQ; i++) sum
        = sum + GLQw[i] * cumnormal_Pol_App((-y1 - M1 - corr * S1 * M_SQRT2 * GLQx[i] + log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k)) / (S1 * sqrt(1.0 - corr * corr)));
    sum = sum / M_SQRTPI;
    return df * notional * sum;
}

double DigitalSpreadPut(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k)
{
    int i;
    double sum = 0.0, notional = 1.0, M1 = log(f1) - 0.5 * sig1 * sig1 * t, M2 = log(f2) - 0.5 * sig2 * sig2 * t, S1 = sig1 * sqrt(t), S2 = sig2 * sqrt(t);
    if (k >= 0.0) for (i = 0; i < NumGLQ; i++) sum = sum + GLQw[i] * cumnormal_Pol_App((log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k) - (y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i])) / (S2 * sqrt(1.0 - corr * corr)));
    else for (i = 0; i < NumGLQ; i++) sum = sum + GLQw[i] * cumnormal_Pol_App((-log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k) + (y1 + M1 + corr * S1 * M_SQRT2 * GLQx[i])) / (S1 * sqrt(1.0 - corr * corr)));
    sum = sum / M_SQRTPI;
    return df * notional * sum;
}

double SpreadCall(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k)
{
    int i;
    double ans = 0.0, sum = 0.0, notional = 1.0, M1, M2, S1 = sig1 * sqrt(t), S2 = sig2 * sqrt(t);
    if (f1 > 0.0 && f2 > 0.0)
    {
        M1 = log(f1) - 0.5 * sig1 * sig1 * t;
        M2 = log(f2) - 0.5 * sig2 * sig2 * t;
        if (k >= 0.0)
        {
            for (i = 0; i < NumGLQ; i++)
                sum = sum + GLQw[i] *
                (exp(M2 + corr * S2 * M_SQRT2 * GLQx[i] + 0.5 * S2 * S2 * (1.0 - corr * corr)) * cumnormal_Pol_App(((y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i] + S2 * S2 * (1.0 - corr * corr)) - log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k)) / (S2 * sqrt(1.0 - corr * corr)))
                    - (exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k) * cumnormal_Pol_App(((y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i]) - log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k)) / (S2 * sqrt(1.0 - corr * corr))));
        }
        else
        {
            for (i = 0; i < NumGLQ; i++)
                sum = sum + GLQw[i] *
               ((exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k) * cumnormal_Pol_App((-y1 - M1 - corr * S1 * M_SQRT2 * GLQx[i] + log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k)) / (S1 * sqrt(1.0 - corr * corr)))
               - exp(M1 + corr * S1 * M_SQRT2 * GLQx[i] + 0.5 * S1 * S1 * (1.0 - corr * corr)) * cumnormal_Pol_App((-y1 - M1 - corr * S1 * M_SQRT2 * GLQx[i] - S1 * S1 * (1.0 - corr * corr) + log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k)) / (S1 * sqrt(1.0 - corr * corr))));
        }
        sum = sum / M_SQRTPI;
        ans = df * notional * sum;
    }
    else if (f1 > 0.0 && f2 <= 0.0)
    {
        if (f2 - k > 0) ans = notional * df * Bl(f2 - k, f1, S1, -1);
    }
    else if (f1 <= 0.0 && f2 > 0.0)
    {
        if (f1 + k > 0) ans = notional * df * Bl(f1 + k, f2, S2, 1);
        else ans = notional * df * (f2 - f1 - k);
    }
    else
    {
        ans = notional * df * max(f2 - f1 - k,0.0);
    }
    return ans;
}

double SpreadPut(double f1, double f2, double sig1, double sig2, double y1, double y2, double corr, double df, double t, double k)
{
    int i;
    double sum = 0.0, notional = 1.0, M1 = log(f1) - 0.5 * sig1 * sig1 * t, M2 = log(f2) - 0.5 * sig2 * sig2 * t, S1 = sig1 * sqrt(t), S2 = sig2 * sqrt(t);
    if (k >= 0.0)
    {
        for (i = 0; i < NumGLQ; i++)
            sum = sum + GLQw[i] *
            (-exp(M2 + corr * S2 * M_SQRT2 * GLQx[i] + 0.5 * S2 * S2 * (1.0 - corr * corr)) * cumnormal_Pol_App((log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k) - (y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i] + S2 * S2 * (1.0 - corr * corr))) / (S2 * sqrt(1.0 - corr * corr)))
                + (exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k) * cumnormal_Pol_App((log(exp(y1 + M1 + S1 * M_SQRT2 * GLQx[i]) + k) - (y2 + M2 + corr * S2 * M_SQRT2 * GLQx[i])) / (S2 * sqrt(1.0 - corr * corr))));
    }
    else
    {
        for (i = 0; i < NumGLQ; i++)
            sum = sum + GLQw[i] *
           (-(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k) * cumnormal_Pol_App((y1 + M1 + corr * S1 * M_SQRT2 * GLQx[i] - log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k)) / (S1 * sqrt(1.0 - corr * corr)))
           + exp(M1 + corr * S1 * M_SQRT2 * GLQx[i] + 0.5 * S1 * S1 * (1.0 - corr * corr)) * cumnormal_Pol_App((y1 + M1 + corr * S1 * M_SQRT2 * GLQx[i] + S1 * S1 * (1.0 - corr * corr) - log(exp(y2 + M2 + S2 * M_SQRT2 * GLQx[i]) - k)) / (S1 * sqrt(1.0 - corr * corr))));
    }
    sum = sum / M_SQRTPI;
    return df * notional * sum;
}

double Capletsum(double fwdvol0, double fwdvol1, double previouscap, double k, vector<double> df, vector<double> capletmaturity, vector<double> capletperiod, vector<double> fra)
{
    int i, numcaplet = int(df.size());
    vector<double> fwdvol(numcaplet);
    double sum = previouscap;
    for (i = 1; i < numcaplet; i++)
    {
        fwdvol[i] = (fwdvol1 - fwdvol0) / (capletmaturity[numcaplet - 1] - capletmaturity[0]) * (capletmaturity[i] - capletmaturity[0]) + fwdvol0;
        sum = sum + Caplet((fra[i]), k, (df[i]), (capletperiod[i]), (capletmaturity[i]), (fwdvol[i]));
    }
    return sum;
}

double Findfwdvol(double fwdvol0, double fullcap, double previouscap, double k, vector<double> df, vector<double> capletmaturity, vector<double> capletperiod, vector<double> fra)
{
    double left = 0.0, right = 10.0, fwdvol1 = 0.5 * (left + right);
    double ans = Capletsum(fwdvol0, fwdvol1, previouscap, k, df, capletmaturity, capletperiod, fra) - fullcap;
    while (fabs(ans) > Error&& fabs(fwdvol1 - left) > Error&& fabs(fwdvol1 - right) > Error)
    {
        if (ans < 0.0)
        {
            left = fwdvol1;
            fwdvol1 = 0.5 * (left + right);
            ans = Capletsum(fwdvol0, fwdvol1, previouscap, k, df, capletmaturity, capletperiod, fra) - fullcap;
        }
        else if (ans > 0.0)
        {
            right = fwdvol1;
            fwdvol1 = 0.5 * (left + right);
            ans = Capletsum(fwdvol0, fwdvol1, previouscap, k, df, capletmaturity, capletperiod, fra) - fullcap;
        }
        else break;
    }
    return fwdvol1;
}

void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<CDate>> swapmaturity_date, vector<vector<double>> mktswaptionvol, vector<CDate> targetswapmaturity_date, CDate option_maturity, CDate swap_maturity, double& swaptionvol)
{
    int i, num_option = int(optionmaturity_date.size()), num_swap = int(targetswapmaturity_date.size());
    vector<CInterpolation> volsrf(num_option);
    vector<double> volx(num_swap);
    for (i = 0; i < num_swap; i++)
    {
        volsrf[i] = CInterpolation(optionmaturity_date, mktswaptionvol[i]);
        volx[i] = volsrf[i](option_maturity);
    }
    CInterpolation vol(targetswapmaturity_date, volx);
    swaptionvol = vol(swap_maturity);
}

void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<double>> mktswaptionvol, vector<CDate> targetswapmaturity_date, CDate option_maturity, CDate swap_maturity, double& swaptionvol)
{
    int i, num_option = int(optionmaturity_date.size()), num_swap = int(targetswapmaturity_date.size());
    vector<CInterpolation> volsrf(num_option);
    vector<double> volx(num_swap);
    for (i = 0; i < num_swap; i++)
    {
        volsrf[i] = CInterpolation(optionmaturity_date, mktswaptionvol[i]);
        volx[i] = volsrf[i](option_maturity);
    }
    CInterpolation vol(targetswapmaturity_date, volx);
    swaptionvol = vol(swap_maturity);
}

//murex case : find swaption volatility - Option : 날짜로 찾음, Swap : 옵션일 다음날 부터 스케쥴 작성하여 기간으로 찾음
void FindSwaptionVol(vector<CDate> optionmaturity_date, vector<vector<double>> mktswaptionvol, vector<double> targetswapmaturity_t, CDate option_maturity, double swap_t, double& swaptionvol)
{
    int i, num_option = int(optionmaturity_date.size()), num_swap = int(targetswapmaturity_t.size());
    vector<CInterpolation> volsrf(num_option);
    vector<double> volx(num_swap);

    for (i = 0; i < num_swap; i++)
    {
        volsrf[i] = CInterpolation(optionmaturity_date, mktswaptionvol[i]);
        volx[i] = volsrf[i](option_maturity);
    }
    CInterpolation vol(targetswapmaturity_t, volx);
    swaptionvol = vol(swap_t);
}

void PlainSwaptionPrice(CCurrency crcy, CDate today, CDate option_maturity, CDate settlement_date, double notional, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<double> fundingleg_notional, vector<double> fundingleg_mult, vector<double> fundingleg_spread, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, bool estflag, string payout_type, bool exercise_flag, string delivery_type, double swaptionvol, double& atmswaprate, double& optionprice)
{
    atmswaprate = 0.0;
    optionprice = 0.0;
    double floatinglegprice = 0.0, fixedlegprice = 0.0, sumspread = 0.0;
    if (today < option_maturity)
    {
        string dcb = "ACT/365";
        double ondf = 1.0;
        atmswaprate = 0.0;
        double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0;
        int i, starti = 0, endi, num_couponleg_cf = int(couponleg_calc_startdate.size()), num_fundingleg_cf = int(fundingleg_calc_startdate.size());
        vector<double> couponleg_df(num_couponleg_cf, 1.0), coupon(num_couponleg_cf), fra(num_fundingleg_cf, 0.0), fundingleg_df(num_fundingleg_cf, 1.0);

        while (today >= fundingleg_paydate[starti])
        {
            starti = starti + 1;
            if (starti >= num_fundingleg_cf - 1) break;
        }
        for (i = starti; i < num_fundingleg_cf; i++)
        {
            fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
            fra[i] = fixinghistory_rate[i];
            //if(fra[i]>0.0)
            if (fra[i] > MinFixingRate)
            {
                floatinglegprice = floatinglegprice + fundingleg_mult[i] * fra[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
                sumspread = sumspread + fundingleg_spread[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
            }
            if (today <= fixing_date[i]) break;
        }
        if (today == fundingleg_paydate[num_fundingleg_cf - 1]) floatinglegprice = 0.0;
        endi = i;
        if (endi < num_fundingleg_cf)
        {
            //if(today==fixing_date[i] && fra[i]>0.0)
            if (today == fixing_date[i] && fra[i] > MinFixingRate)
            {
                endi = i + 1;
                CDate ondate;
                CInstrument temp(crcy, "DEPO");
                ShiftBusDate(today, temp.get_holiday(), temp.get_numholiday(), 1, ondate);
                ondf = exp(-zc(ondate) * cvg(today, ondate, dcb));
            }
        }
        settlement_date_df = settlement_date_df / ondf;
        for (i = endi; i < num_fundingleg_cf; i++)
        {
            fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
            fra[i] = (exp(-zc(fundingleg_calc_startdate[i]) * cvg(today, fundingleg_calc_startdate[i], dcb) + zc(fixingindex_maturity[i]) * cvg(today, fixingindex_maturity[i], dcb)) - 1.0) / cvg(fundingleg_calc_startdate[i], fixingindex_maturity[i], fixingindex_dcb);
            floatinglegprice = floatinglegprice + (estflag)*fundingleg_mult[i] * fra[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
            sumspread = sumspread + fundingleg_spread[i] * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        }
        if (endi < num_fundingleg_cf) floatinglegprice = floatinglegprice + (1 - estflag) * (exp(-zc(fundingleg_calc_startdate[endi]) * cvg(today, fundingleg_calc_startdate[endi], dcb)) - fundingleg_df[num_fundingleg_cf - 1]);
        starti = 0;
        while (today >= couponleg_paydate[starti])
        {
            starti = starti + 1;
            if (starti >= num_couponleg_cf - 1) break;
        }
        for (i = starti; i < num_couponleg_cf; i++)
        {
            couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
            levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
            level = level + levelcvg * couponleg_df[i];
            coupon[i] = couponleg_couponrate[i] * levelcvg;
            fixedlegprice = fixedlegprice + coupon[i] * couponleg_df[i];
        }
        if (today == couponleg_paydate[num_couponleg_cf - 1]) fixedlegprice = 0.0;
        if (starti < num_couponleg_cf)
        {
            atmswaprate = floatinglegprice / level;
            sumspread = sumspread / level;
            //fixedlegprice=fixedlegprice/level;
        }
        fixedlegprice = fixedlegprice / ondf / settlement_date_df;

        if (payout_type == "PAYER")
        {
            if (couponleg_couponrate[0] - sumspread > 0) optionprice = Bl(couponleg_couponrate[0] - sumspread, atmswaprate, swaptionvol * sqrt(cvg(today, option_maturity, dcb)), 1) * level / settlement_date_df;
            //else optionprice=(atmswaprate+sumspread-fixedlegprice)*level/settlement_date_df;
            else optionprice = (atmswaprate + sumspread) * level / ondf / settlement_date_df - fixedlegprice;
        }
        else if (payout_type == "RECEIVER")
        {
            if (couponleg_couponrate[0] - sumspread > 0) optionprice = Bl(couponleg_couponrate[0] - sumspread, atmswaprate, swaptionvol * sqrt(cvg(today, option_maturity, dcb)), -1) * level / settlement_date_df;
        }
        else
        {
            if (couponleg_couponrate[0] - sumspread > 0) optionprice = (Bl(couponleg_couponrate[0] - sumspread,atmswaprate,swaptionvol * sqrt(cvg(today,option_maturity,dcb)),1) + Bl(couponleg_couponrate[0] - sumspread,atmswaprate,swaptionvol * sqrt(cvg(today,option_maturity,dcb)),-1)) * level / settlement_date_df;
            else optionprice = (atmswaprate + sumspread) * level / ondf / settlement_date_df - fixedlegprice;
        }
        optionprice = optionprice * notional;
    }
}


void SwaptionPrice_Bl(CDate today, CDate option_maturity, CDate settlement_date, CInterpolation zc, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, double swaptionvol, double atmswaprate, double& optionprice)
{
    optionprice = 0.0;
    string dcb = "ACT/365";
    double settlement_date_df = exp(-zc(settlement_date) * cvg(today, settlement_date, dcb)), level = 0.0, levelcvg = 0.0;
    int i, num_couponleg_cf = int(couponleg_calc_startdate.size());
    vector<double> couponleg_df(num_couponleg_cf, 1.0), coupon(num_couponleg_cf);

    for (i = 0; i < num_couponleg_cf; i++)
    {
        couponleg_df[i] = exp(-zc(couponleg_paydate[i]) * cvg(today, couponleg_paydate[i], dcb));
        levelcvg = cvg(couponleg_calc_startdate[i], couponleg_calc_enddate[i], couponleg_dcb);
        level = level + levelcvg * couponleg_df[i];
    }
    optionprice = Bl(atmswaprate, atmswaprate, swaptionvol * sqrt(cvg(today, option_maturity, dcb)), 1) * level / settlement_date_df;
}

/*
void PlainCMSSpreadRangeAccrualSwap(CCurrency crcy, CDate today, CDate settlement_date, double notional, CInterpolation zc, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<double> fundingleg_notional, vector<double> fundingleg_mult, vector<double> fundingleg_spread, vector<CDate> fixing_date, vector<CDate> fixingindex_maturity, string fixingindex_dcb, vector<CDate> fixinghistory_date, vector<double> fixinghistory_rate, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, CDate *couponleg_holidays, int num_couponlegholidays, string couponleg_dcb, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> rahigh_bdryoh, vector<int> rahigh_bdryoh_flag, vector<double>
ralow_bdry, vector<bool> ralow_bdryin_flag, vector<double> ralow_bdryoh, vector<int> ralow_bdryoh_flag, vector<CDate> couponindexfixinghistory_date, vector<double> couponindex1fixinghistory_rate, vector<double> couponindex2fixinghistory_rate, vector<bool> ra_flag, vector<CCurrency> couponlegindex1_crcys, vector<string> couponlegindex1_tenors, vector<string> couponlegindex1_types, vector<string> couponlegindex1fixed_freqs, vector<CDate *> couponlegindex1_holidayss, vector<int> num_couponlegindex1holidayss, vector<CCurrency> couponlegindex2_crcys, vector<string> couponlegindex2_tenors, vector<string> couponlegindex2_types, vector<string> couponlegindex2fixed_freqs, vector<CDate *> couponlegindex2_holidayss, vector<int> num_couponlegindex2holidayss, int nshift, int nlockout, bool rafirstdatein_flag, bool ralastdatein_flag, int raoptm_prd1, int raoptm_frq1, int raoptm_prd2, int raoptm_frq2, double corr, vector<CDate> swaptionoptionmaturity_date, vector<string> swaptionswaptenor, vector<vector<double>> mktswapti
nvol, double &floatinglegprice, double &fixedlegprice, double &atmswaprate, double &swapprice)
{
    int num_couponleg_cf=int(couponleg_paydate.size());
    vector<double> temp_couponrate(num_couponleg_cf,0.0);
    double tempfloatinglegprice, tempfixedlegprice, tempatmswaprate, tempswapprice;
    PlainSwapPrice(crcy,today,settlement_date,notional,zc,couponleg_calc_startdate,couponleg_calc_enddate,couponleg_paydate,couponleg_dcb,couponleg_notional,temp_couponrate,fundingleg_calc_startdate,fundingleg_calc_enddate,fundingleg_paydate,fundingleg_dcb,fundingleg_notional,fundingleg_mult,fundingleg_spread,fixing_date,fixingindex_maturity,fixingindex_dcb,fixinghistory_date,fixinghistory_rate,true,floatinglegprice,tempfixedlegprice,tempatmswaprate,tempswapprice);

    string dcb="ACT/365";
    CDate raprd1=DateAdd(m,raoptm_prd1,today), raprd2=DateAdd(m,raoptm_prd2,today);

    int couponleg_cf_starti=0;
    while(today>=couponleg_paydate[couponleg_cf_starti])
    {
        couponleg_cf_starti=couponleg_cf_starti+1;
        if(couponleg_cf_starti>=num_couponleg_cf-1) break;
    }

    int i, j, k, couponindexfixinghistory_starti;
    vector<CDate> couponleg_calc_shiftedstartdate, couponleg_calc_shiftedenddate, lockout_startdate;
    couponindexfixinghistory_starti=0;
    vector<double> fixedlegprices(num_couponleg_cf-couponleg_cf_starti,0.0);
    vector<double> tempr1, tempr2, swapfra1, swapfra2;
    vector<int> num_fixing, num_lockout, num_realfixing;
    vector<vector<CDate>> index1_cf(3), index2_cf(3);
    string index1_dcb, index1_stub, index1_direction, index1_conv, index1_adjflag, index2_dcb, index2_stub, index2_direction, index2_conv, index2_adjflag;
    int num_fixedrate_in=0, tempi, tempfixingstarti=0, tempfixingendi=0, num_index1_coup_cf, num_index2_coup_cf, index1_spotlag, index2_spotlag;
    double level1, level2;
    for(i=couponleg_cf_starti;i<num_couponleg_cf;i++)
    {
        if(!ra_flag[i])
        {
            vector<CDate> temp_startdate(1), temp_enddate(1), temp_paydate(1);
            vector<double> temp_notional(1), temp_rate(1);
            temp_startdate[0]=couponleg_calc_startdate[i];
            temp_enddate[0]=couponleg_calc_enddate[i];
            temp_paydate[0]=couponleg_paydate[i];
            temp_notional[0]=couponleg_notional[i];
            temp_rate[0]=couponleg_couponrate[i];
            PlainSwapPrice(crcy,today,settlement_date,notional,zc,temp_startdate,temp_enddate,temp_paydate,couponleg_dcb,temp_notional,temp_rate,fundingleg_calc_startdate,fundingleg_calc_enddate,fundingleg_paydate,fundingleg_dcb,fundingleg_notional,fundingleg_mult,fundingleg_spread,fixing_date,fixingindex_maturity,fixingindex_dcb,fixinghistory_date,fixinghistory_rate,true,tempfloatinglegprice,fixedlegprices[i-couponleg_cf_starti],tempatmswaprate,tempswapprice);
        }
        Else
        {
            CDate nowdate, startdate1, startdate2;
            CInstrument index1(couponlegindex1_crcys[i],couponlegindex1_types[i]), index2(couponlegindex2_crcys[i],couponlegindex2_types[i]);
            index1_dcb=index1.get_basis();
            index1_stub=index1.get_stub();
            index1_direction=index1.get_direction();
            index1_conv=index1.get_conv();
            index1_adjflag=index1.get_adjflag();
            index1_spotlag=index1.get_spotlag();
            index2_dcb=index2.get_basis();
            index2_stub=index2.get_stub();
            index2_direction=index2.get_direction();
            index2_conv=index2.get_conv();
            index2_adjflag=index2.get_adjflag();
            index2_spotlag=index2.get_spotlag();

            int fixedrateincount=0, tempcount=0;
            couponleg_calc_shiftedstartdate.push_back(couponleg_calc_startdate[i]);
            couponleg_calc_shiftedenddate.push_back(couponleg_calc_enddate[i]);
            tempi=int(couponleg_calc_shiftedstartdate.size())-1;
            ShiftBusDate(couponleg_calc_shiftedstartdate[tempi],couponleg_holidays,num_couponlegholidays,nshift,couponleg_calc_shiftedstartdate[tempi]);
            ShiftBusDate(couponleg_calc_shiftedenddate[tempi],couponleg_holidays,num_couponlegholidays,nshift,couponleg_calc_shiftedenddate[tempi]);
            lockout_startdate.push_back(couponleg_calc_shiftedenddate[tempi]);
            ShiftBusDate(lockout_startdate[tempi],couponleg_holidays,num_couponlegholidays,nlockout,lockout_startdate[tempi]);
            num_fixing.push_back(int(couponleg_calc_shiftedenddate[tempi]-couponleg_calc_shiftedstartdate[tempi])+(ralastdatein_flag-1)+rafirstdatein_flag);
            num_lockout.push_back(int(couponleg_calc_shiftedenddate[tempi]-lockout_startdate[tempi])+(ralastdatein_flag-1));
            num_realfixing.push_back(num_fixing[tempi]-num_lockout[tempi]);
            if(int(couponindexfixinghistory_date.size())>0 && couponindexfixinghistory_starti<int(couponindexfixinghistory_date.size()))
            {
                while(couponindexfixinghistory_date[couponindexfixinghistory_starti]<couponleg_calc_shiftedstartdate[tempi])
                {
                    couponindexfixinghistory_starti=couponindexfixinghistory_starti+1;
                    if(couponindexfixinghistory_starti>=int(couponindexfixinghistory_date.size())) break;
                }
                for(j=couponindexfixinghistory_starti;j<int(couponindexfixinghistory_date.size()) && (couponindexfixinghistory_date[j]<couponleg_calc_shiftedenddate[tempi] || (ralastdatein_flag&&couponindexfixinghistory_date[j]<=couponleg_calc_shiftedenddate[tempi]));j++)
                {
                    if(is_holiday(couponindexfixinghistory_date[j],couponlegindex1_holidayss[i],num_couponlegindex1holidayss[i])){ if(int(tempr1.size())>0) tempr1.push_back(tempr1[int(tempr1.size())-1]);}
                    else tempr1.push_back(couponindex1fixinghistory_rate[j]);
                    if(is_holiday(couponindexfixinghistory_date[j],couponlegindex2_holidayss[i],num_couponlegindex2holidayss[i])){ if(int(tempr2.size())>0) tempr2.push_back(tempr2[int(tempr2.size())-1]);}
                    else tempr2.push_back(couponindex2fixinghistory_rate[j]);
                }
                tempfixingendi=int(tempr1.size());
            }
            for(j=tempfixingstarti;j<tempfixingendi;j++)
            {
                if(j<=tempfixingstarti+num_realfixing[tempi]-1 && rangeincheck(couponlegindex1_mult[i]*tempr1[j]+couponlegindex2_mult[i]*tempr2[j],rahigh_bdry[i],ralow_bdry[i],rahigh_bdryin_flag[i],ralow_bdryin_flag[i]) && tempr1[j]>MinFixingRate)
                {
                    fixedrateincount=fixedrateincount+1;
                    tempcount=tempcount+1;
                }
                if(j>=tempfixingstarti+num_realfixing[tempi] && j<=tempfixingstarti+num_fixing[tempi]-1 && rangeincheck(couponlegindex1_mult[i]*tempr1[tempfixingstarti+num_realfixing[tempi]-1]+couponlegindex2_mult[i]*tempr2[tempfixingstarti+num_realfixing[tempi]-1],rahigh_bdry[i],ralow_bdry[i],rahigh_bdryin_flag[i],ralow_bdryin_flag[i]))
                {
                    fixedrateincount=fixedrateincount+1;
                    tempcount=tempcount+1;
                }
            }
            tempfixingstarti=tempfixingendi;
            if(tempcount<num_realfixing[tempi])
            {
                for(j=tempcount;j<num_realfixing[tempi];j++)
                {
                    nowdate=DateAdd(d,j,couponleg_calc_shiftedstartdate[tempi]);
                    ShiftBusDate(nowdate,couponlegindex1_holidayss[i],num_couponlegindex1holidayss[i],index1_spotlag,startdate1);
                    num_index1_coup_cf=findnumschedule(startdate1,couponlegindex1_tenors[i],index1_stub,index1_direction,couponlegindex1fixed_freqs[i]);
                    ShiftBusDate(nowdate,couponlegindex2_holidayss[i],num_couponlegindex2holidayss[i],index2_spotlag,startdate2);
                    num_index2_coup_cf=findnumschedule(startdate2,couponlegindex2_tenors[i],index2_stub,index2_direction,couponlegindex2fixed_freqs[i]);
                    for(k=0;k<3;k++)
                    {
                        index1_cf[k]=vector<CDate>(num_index1_coup_cf);
                        index2_cf[k]=vector<CDate>(num_index2_coup_cf);
                    }
                    fixedlegcashflowschedule(nowdate,couponlegindex1_tenors[i],couponlegindex1_holidayss[i],num_couponlegindex1holidayss[i],index1_stub,index1_direction,index1_conv,couponlegindex1fixed_freqs[i],index1_adjflag,index1_cf);
                    level1=0.0;
                    for(k=0;k<num_index1_coup_cf;k++) level1=level1+cvg(index1_cf[1][i],index1_cf[2][i],index1_dcb)*exp(-zc(index1_cf[0][i])*cvg(today,index1_cf[0][i],dcb));
                    swapfra1.push_back((exp(-zc(index1_cf[1][0])*cvg(today,index1_cf[1][0],dcb))-exp(-zc(index1_cf[0][num_index1_coup_cf-1])*cvg(today,index1_cf[0][num_index1_coup_cf-1],dcb)))/level1);
                    fixedlegcashflowschedule(nowdate,couponlegindex2_tenors[i],couponlegindex2_holidayss[i],num_couponlegindex2holidayss[i],index2_stub,index2_direction,index2_conv,couponlegindex2fixed_freqs[i],index2_adjflag,index2_cf);
                    level2=0.0;
                    for(k=0;k<num_index2_coup_cf;k++) level2=level2+cvg(index2_cf[1][i],index2_cf[2][i],index2_dcb)*exp(-zc(index2_cf[0][i])*cvg(today,index2_cf[0][i],dcb));
                    swapfra2.push_back((exp(-zc(index2_cf[1][0])*cvg(today,index2_cf[1][0],dcb))-exp(-zc(index2_cf[0][num_index2_coup_cf-1])*cvg(today,index2_cf[0][num_index2_coup_cf-1],dcb)))/level2);
                }
            }
            else if(tempcount<num_fixing[tempi])
            {
                for(j=tempcount;j<num_fixing[tempi];j++)
                {

                }
            }
            fixedlegprices[i-couponleg_cf_starti]=couponleg_notional[i]*double(fixedrateincount)/double(num_fixing[tempi])*couponleg_couponrate[i]*cvg(couponleg_calc_startdate[i],couponleg_calc_enddate[i],couponleg_dcb)*exp(-zc(couponleg_paydate[i])*cvg(today,couponleg_paydate[i],dcb));
        }
    }
}
*/

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
)
{
    int i;
    string dcb = "ACT/365";
    vector<double> fundingleg_df(num_fundingleg_cf);
    floatingleg_fixed_price = 0.0;

    if (today < fundingleg_paydate[num_fundingleg_cf - 1])
    {
        for (i = fundingleg_cf_starti; i < num_fundingleg_cf; i++)
        {
            fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));

            if (fixinghistory_rate[i] > MinFixingRate)
            {
                floatingleg_fixed_price
                    = floatingleg_fixed_price
                    + fundingleg_notional[i]
                    * (fundingleg_mult[i] * fixinghistory_rate[i] + fundingleg_spread[i])
                    * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb)
                    * fundingleg_df[i];
            }

            if (today <= fixing_date[i])
            {
                break;
            }
        }

        floatingleg_fixed_price = floatingleg_fixed_price / ondf / settlement_date_df;
    }
}

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
)
{
    int i;
    string dcb = "ACT/365";
    vector<double> fundingleg_df(num_remained_fundingleg_cf);
    floatinglegnonfixedprice = 0.0;

    if (num_remained_fundingleg_cf > 0)
    {
        if (today < _fundingleg_paydate[num_remained_fundingleg_cf - 1])
        {
            // calculate funding leg spread
            for (i = 0; i < num_remained_fundingleg_cf; i++)
            {
                fundingleg_df[i] = exp(-zc(_fundingleg_paydate[i]) * fundinglegT[i]);
                floatinglegnonfixedprice
                    = floatinglegnonfixedprice
                    + _fundingleg_notional[i]
                    * _fundingleg_spread[i]
                    * fundinglegtau[i]
                    * fundingleg_df[i];
            }

            // calculate funding leg price by short end minus long end
            if (num_remained_fundingleg_cf > 0)
            {
                floatinglegnonfixedprice
                    = floatinglegnonfixedprice
                    + (
                        _fundingleg_notional[0]
                        * exp(-zc(fundingleg_calc_startdate0) * fundinglegcalcstartT0)
                        - _fundingleg_notional[num_remained_fundingleg_cf - 1]
                        * fundingleg_df[num_remained_fundingleg_cf - 1]
                        );
            }

            floatinglegnonfixedprice = floatinglegnonfixedprice / ondf / settlement_date_df;
        }
    }
}

void PlainSwapFundinglegNonFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_fundingleg_cf, vector<double> _fundingleg_notional, vector<double> _fundingleg_mult, vector<double> _fundingleg_spread, CDate fundingleg_calc_startdate0, double fundinglegcalcstartT0, vector<CDate> _fundingleg_calc_startdate, vector<CDate> _fundingleg_calc_enddate, vector<CDate> _fundingleg_paydate, vector<double> fundinglegtau, vector<double> fundinglegT, double& floatinglegnonfixedprice)
{
    int i;
    string dcb = "ACT/365";
    vector<double> fundingleg_df(num_remained_fundingleg_cf);
    floatinglegnonfixedprice = 0.0;
    if (num_remained_fundingleg_cf > 0)
    {
        if (today < _fundingleg_paydate[num_remained_fundingleg_cf - 1])
        {
            for (i = 0; i < num_remained_fundingleg_cf; i++)
            {
                fundingleg_df[i] = exp(-zc(_fundingleg_paydate[i]) * fundinglegT[i]);
                floatinglegnonfixedprice = floatinglegnonfixedprice + _fundingleg_notional[i] * _fundingleg_spread[i] * fundinglegtau[i] * fundingleg_df[i];
            }

            //if(num_remained_fundingleg_cf>0) floatinglegnonfixedprice=floatinglegnonfixedprice+(_fundingleg_notional[0]*exp(-zc(fundingleg_calc_startdate0)*fundinglegcalcstartT0)-_fundingleg_notional[num_remained_fundingleg_cf-1]*fundingleg_df[num_remained_fundingleg_cf-1]);

            for (i = 0; i < num_remained_fundingleg_cf; i++)
            {
                floatinglegnonfixedprice = floatinglegnonfixedprice + _fundingleg_notional[i] * _fundingleg_mult[i] * ((exp(-zc(_fundingleg_calc_startdate[i]) * cvg(today, _fundingleg_calc_startdate[i], dcb) + zc(_fundingleg_calc_enddate[i]) * cvg(today, _fundingleg_calc_enddate[i], dcb)) - 1.0) / fundinglegtau[i]) * fundinglegtau[i] * fundingleg_df[i];
            }

            floatinglegnonfixedprice = floatinglegnonfixedprice / ondf / settlement_date_df;
        }
    }
}

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
)
{
    int i;
    fixedlegfixedprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (!_ra_flag[i])
        {
            fixedlegfixedprice
                = fixedlegfixedprice
                + exp(-zc(_couponleg_paydate[i]) * couponlegT[i])
                * _couponleg_notional[i]
                * couponlegtau[i]
                * _couponleg_couponrate[i];
        }
    }

    fixedlegfixedprice = fixedlegfixedprice / ondf / settlement_date_df;
}

void PlainSwapZeroCouponlegFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i;
    fixedlegfixedprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (!_ra_flag[i])
        {
            fixedlegfixedprice = fixedlegfixedprice + _couponleg_notional[i] * couponlegtau[i] * _couponleg_couponrate[i];
        }
    }
    fixedlegfixedprice = exp(-zc(_couponleg_paydate[num_remained_couponleg_cf - 1]) * couponlegT[num_remained_couponleg_cf - 1]) * fixedlegfixedprice / ondf / settlement_date_df;
}

////////////////////////////////////////////////////////////////////////////////////
void CMSSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i, j;//, n=min(num_remained_couponleg_cf,2);
    fixedlegfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    //vector<double> fixedlegfixedcf(n,0.0);
    //for(i=0;i<n;i++)
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            //fixedlegfixedcf[i]=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(fixedrateincount/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
            //fixedlegfixedprice=fixedlegfixedprice+fixedlegfixedcf[i];
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedrateincount / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegfixedprice = fixedlegfixedprice + tempprice;
        }
    }
    fixedlegfixedprice = fixedlegfixedprice / ondf / settlement_date_df;
}
////////////////////////////////////////////////////////////////////////////////////
void CMSSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i, j;
    fixedlegfixedprice = 0.0;
    double fixedrateincount = 0.0, ra_coupon_fixedprice = 0.0, fixed_coupon_fixedprice = 0.0, tempprice = 0.0;
    bool temprangeincheckflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            fixedrateincount = 0.0;
            temprangeincheckflag = false;
            if (num_fixed[i] > 0) { if (fixedrate1[i][num_fixed[i] - 1] <= MinFixingRate) num_fixed[i] = num_fixed[i] - 1; }
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (couponindex_fixeddate[i][j] <= lockout_startdate[i])
                {
                    temprangeincheckflag = rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]);
                    fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
                }
                else fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedrateincount / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            ra_coupon_fixedprice = ra_coupon_fixedprice + tempprice;
        }
        else
        {
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (_couponleg_couponrate[i] + _couponleg_spread[i]);
            fixed_coupon_fixedprice = fixed_coupon_fixedprice + tempprice;
        }
    }
    fixedlegfixedprice = (ra_coupon_fixedprice + fixed_coupon_fixedprice) / ondf / settlement_date_df;
}
/////////////////////////////////////////////////////////////////////////////////////
void FixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i, j;//, n=min(num_remained_couponleg_cf,2);
    fixedlegfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    //vector<double> fixedlegfixedcf(n,0.0);
    //for(i=0;i<n;i++)
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            //fixedlegfixedcf[i]=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(fixedrateincount/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
            //fixedlegfixedprice=fixedlegfixedprice+fixedlegfixedcf[i];
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedrateincount / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegfixedprice = fixedlegfixedprice + tempprice;
        }
    }
    fixedlegfixedprice = fixedlegfixedprice / ondf / settlement_date_df;
}
////////////////////////////////////////////////////////////////////////////////////
void FixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i, j;
    fixedlegfixedprice = 0.0;
    double fixedrateincount = 0.0, ra_coupon_fixedprice = 0.0, fixed_coupon_fixedprice = 0.0, tempprice = 0.0;
    bool temprangeincheckflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            fixedrateincount = 0.0;
            temprangeincheckflag = false;
            if (num_fixed[i] > 0) { if (fixedrate1[i][num_fixed[i] - 1] <= MinFixingRate) num_fixed[i] = num_fixed[i] - 1; }
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (couponindex_fixeddate[i][j] <= lockout_startdate[i])
                {
                    temprangeincheckflag = rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]);
                    fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
                }
                else fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedrateincount / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            ra_coupon_fixedprice = ra_coupon_fixedprice + tempprice;
        }
        else
        {
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (_couponleg_couponrate[i] + _couponleg_spread[i]);
            fixed_coupon_fixedprice = fixed_coupon_fixedprice + tempprice;
        }
    }
    fixedlegfixedprice = (ra_coupon_fixedprice + fixed_coupon_fixedprice) / ondf / settlement_date_df;
}
/////////////////////////////////////////////////////////////////////////////////////
void FixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<CDate> couponleg_calc_shiftedstartdate, vector<CDate> couponleg_calc_shiftedenddate, vector<int> num_fixing, vector<int> num_fixed, vector<CDate> lockout_startdate, vector<vector<CDate>> couponindex_fixeddate, vector<vector<double>> fixedrate1, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i, j;
    fixedlegfixedprice = 0.0;
    double fixedrateincount = 0.0, ra_coupon_fixedprice = 0.0, fixed_coupon_fixedprice = 0.0, tempprice = 0.0;
    bool temprangeincheckflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            fixedrateincount = 0.0;
            temprangeincheckflag = false;
            if (num_fixed[i] > 0) { if (fixedrate1[i][num_fixed[i] - 1] <= MinFixingRate) num_fixed[i] = num_fixed[i] - 1; }
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (couponindex_fixeddate[i][j] <= lockout_startdate[i])
                {
                    temprangeincheckflag = rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]);
                    fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
                }
                else fixedrateincount = fixedrateincount + temprangeincheckflag * 1.0;
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(fixedrateincount/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
            tempprice = _couponleg_notional[i] * couponlegtau[i] * (fixedrateincount / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            ra_coupon_fixedprice = ra_coupon_fixedprice + tempprice;
        }
        else
        {
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(_couponleg_couponrate[i]+_couponleg_spread[i]);
            tempprice = _couponleg_notional[i] * couponlegtau[i] * (_couponleg_couponrate[i] + _couponleg_spread[i]);
            fixed_coupon_fixedprice = fixed_coupon_fixedprice + tempprice;
        }
    }
    fixedlegfixedprice = exp(-zc(_couponleg_paydate[num_remained_couponleg_cf - 1]) * couponlegT[num_remained_couponleg_cf - 1]) * (ra_coupon_fixedprice + fixed_coupon_fixedprice) / ondf / settlement_date_df;
}
/////////////////////////////////////////////////////////////////////////////////////
void QuantoCMSSteepnerFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + min(max(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponlegindex1_mult[i] * _couponleg_spread[i], _floorrates[i]), _caprates[i]);
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

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

)
{
    int i, j;
    //fixedlegaccrualfixedprice=0.0;
    //fixedrate1count=0.0, fixedrate2count=0.0;//, tempprice=0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        //fixedrate1count=0.0;
        //fixedrate2count=0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (avg_fixedflag[i][j] && j <= 0)
                {
                    fixedrate1count[i] = fixedrate1count[i] + _couponlegindex1_mult[i] * fixedrate1[i][j];
                    fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i] * fixedrate2[i][j];
                }
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*min(max(fixedrate1count/fixedrate2count-1+_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice+tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
}

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
)
{
    int i, j;

    int NumAvgCount = 0;

    fixedlegaccrualfixedprice = 0.0;

    double tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        //fixedrate1count=0.0;
        //fixedrate2count=0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (avg_fixedflag[i][j])
                {
                    fixedrate1count[i] = fixedrate1count[i] + _couponlegindex1_mult[i] * fixedrate1[i][j];
                    fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i] * fixedrate2[i][j];

                    NumAvgCount = NumAvgCount + 1;
                }
            }

            fixedrate1count[i] = fixedrate1count[i] / NumAvgCount;
            fixedrate2count[i] = fixedrate2count[i] / NumAvgCount;

            NumAvgCount = 0;

            tempprice
                = exp
                (
                    -zc(_couponleg_paydate[num_remained_couponleg_cf - 1])
                    * couponlegT[num_remained_couponleg_cf - 1]
                )
                * _couponleg_notional[i]
                * couponlegtau[i]
                * min
                (
                    max
                    (
                        fixedrate1count[i] + fixedrate2count[i]
                        + _couponleg_spread[i]
                        , _floorrates[i]
                    )
                    , _caprates[i]
                );

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}
/*2016-12-12 생성 by seoyeojin
 * KRW CMS Fixed Floater 일 경우 => 예를 들어 7.1*(Average 10Y CMS -3.5%)
 *index2가 상수로 고정일 경우 index2_mult= -7.1*3.5%를 입력
 */

void CMSSlopeFixedAccrualPriceForFixedfloater(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count)
{
    int i, j;
    //fixedlegaccrualfixedprice=0.0;
    //fixedrate1count=0.0, fixedrate2count=0.0;//, tempprice=0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        //fixedrate1count=0.0;
        //fixedrate2count=0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (avg_fixedflag[i][j] && j <= 0)
                {
                    fixedrate1count[i] = fixedrate1count[i] + _couponlegindex1_mult[i] * fixedrate1[i][j];
                    fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i];//*fixedrate2[i][j];
                }
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*min(max(fixedrate1count/fixedrate2count-1+_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice+tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
}


void CMSSlopeFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count)
{
    int i, j;
    //fixedlegaccrualfixedprice=0.0;
    //fixedrate1count=0.0, fixedrate2count=0.0;//, tempprice=0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        //fixedrate1count=0.0;
        //fixedrate2count=0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (avg_fixedflag[i][j] && j == 0)
                {
                    fixedrate1count[i] = _couponlegindex1_mult[i] * fixedrate1[i][j];
                    fixedrate2count[i] = _couponlegindex2_mult[i] * fixedrate2[i][j];
                }
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*min(max(fixedrate1count/fixedrate2count-1+_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice+tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
}

void CMSSlopeFixedPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> avg_fixeddenom, vector<vector<bool>> avg_fixedflag, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrate1count, vector<double>& fixedrate2count)
{
    int i, j;
    //fixedlegaccrualfixedprice=0.0;
    //fixedrate1count=0.0, fixedrate2count=0.0;//, tempprice=0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        //fixedrate1count=0.0;
        //fixedrate2count=0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (avg_fixedflag[i][j] && j == 0)
                {
                    fixedrate1count[i] = _couponlegindex1_mult[i] * fixedrate1[i][j];
                    fixedrate2count[i] = _couponlegindex2_mult[i] * fixedrate2[i][j];
                }
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*min(max(fixedrate1count/fixedrate2count-1+_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice+tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
}


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
)
{
    int i, j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0
        , tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            fixedrate0[i][j]
                            , _rahigh_bdry[i]
                            , _ralow_bdry[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                        )
                        &&
                        rangeincheck
                        (
                            _couponlegindex1_mult[i] * fixedrate1[i][j]
                            + _couponlegindex2_mult[i] * fixedrate2[i][j]
                            , _caprates[i]
                            , _floorrates[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
            }

            tempprice
                = exp(-zc(_couponleg_paydate[i]) * couponlegT[i])
                * _couponleg_notional[i]
                * couponlegtau[i]
                * (_couponleg_spread[i] + _couponleg_couponrate[i] * fixedrateincount / num_fixing[i]);

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }

    fixedlegaccrualfixedprice
        = fixedlegaccrualfixedprice
        / ondf
        / settlement_date_df;
}


void CDCMSSpreadDualFixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate0, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(fixedrate0[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _caprates[i], _floorrates[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;//min(max(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponlegindex1_mult[i]*_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(_couponleg_spread[i]+_couponleg_couponrate[i]*fixedrateincount/num_fixing[i]);
            tempprice = _couponleg_notional[i] * couponlegtau[i] * (_couponleg_spread[i] + _couponleg_couponrate[i] * fixedrateincount / (num_fixing[i] + num_fixed[i]));
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
}



//2017-07-13 추가 by yj
void CDCMSSpreadDualLeverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate0, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _couponleg_ra_mult, vector<double> _couponleg_ra_sp, vector<double> _couponleg_ra_cap, vector<double> _couponleg_ra_floor, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(fixedrate0[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _caprates[i], _floorrates[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;//min(max(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponlegindex1_mult[i]*_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (_couponleg_spread[i] + _couponleg_couponrate[i] * min(max(_couponleg_ra_mult[i] * fixedrateincount / num_fixing[i] + _couponleg_ra_sp[i], _couponleg_ra_floor[i]), _couponleg_ra_cap[i]));
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

//2017-08-23 추가 by yj
void QuantoCMSSpreadDualLeverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed,
    vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponleg_notional,
    vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult,
    vector<double> _rahigh_bdry1, vector<bool> _rahigh_bdryin_flag1, vector<double> _ralow_bdry1, vector<bool> _ralow_bdryin_flag1,
    vector<double> _rahigh_bdry2, vector<bool> _rahigh_bdryin_flag2, vector<double> _ralow_bdry2, vector<bool> _ralow_bdryin_flag2,
    vector<double> _couponleg_ra_mult, vector<double> _couponleg_ra_sp, vector<double> _couponleg_ra_cap, vector<double> _couponleg_ra_floor, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry1[i], _ralow_bdry1[i], _rahigh_bdryin_flag1[i], _ralow_bdryin_flag1[i]) &&
                    rangeincheck(_couponlegindex3_mult[i] * fixedrate3[i][j] + _couponlegindex4_mult[i] * fixedrate4[i][j], _rahigh_bdry2[i], _ralow_bdry2[i], _rahigh_bdryin_flag2[i], _ralow_bdryin_flag2[i]))
                    fixedrateincount = fixedrateincount + 1.0;//min(max(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponlegindex1_mult[i]*_couponleg_spread[i],_floorrates[i]),_caprates[i]);
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (_couponleg_spread[i] + _couponleg_couponrate[i] * min(max(_couponleg_ra_mult[i] * fixedrateincount / num_fixing[i] + _couponleg_ra_sp[i], _couponleg_ra_floor[i]), _couponleg_ra_cap[i]));
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}


//2018-02-22 추가 by yj
//_couponleg_coupon rate 사용안함, f_rate : Max[0%, (KRW IRS 10Y + 0.85%)*N/M 중 floating rate에 해당하는 부분
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
)
{
    int i, j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0,
        tempprice = 0.0,
        struct_rate = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (
                    rangeincheck(
                        _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j],
                        _rahigh_bdry1[i],
                        _ralow_bdry1[i],
                        _rahigh_bdryin_flag1[i],
                        _ralow_bdryin_flag1[i]
                    )
                    &&
                    rangeincheck(
                        _couponlegindex3_mult[i] * fixedrate3[i][j] + _couponlegindex4_mult[i] * fixedrate4[i][j],
                        _rahigh_bdry2[i],
                        _ralow_bdry2[i],
                        _rahigh_bdryin_flag2[i],
                        _ralow_bdryin_flag2[i]
                    )
                    )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }

                struct_rate = f_rate[i][0];
            }

            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i])
                * _couponleg_notional[i]
                * couponlegtau[i]
                * min(
                    max(
                        _couponleg_couponrate[i] + _couponleg_ra_mult[i] * (struct_rate + _couponleg_spread[i]) + _couponleg_ra_sp[i],
                        _couponleg_ra_floor[i]
                    ),
                    _couponleg_ra_cap[i]
                )
                * fixedrateincount / num_fixing[i];

            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

///////////////////////////////////////////////////////////////
void CMSSpreadNonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>>index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        nonfixedrateincount = 0.0;
        index1_fra[i] = vector<double>(num_estm[i]);
        index2_fra[i] = vector<double>(num_estm[i]);
        index1_adjrate[i] = vector<double>(num_estm[i]);
        index2_adjrate[i] = vector<double>(num_estm[i]);
        for (j = 0; j < num_estm[i]; j++)
        {
            level1 = 0.0;
            level2 = 0.0;
            for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
            for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
            index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
            index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
            index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
            index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], couponindex_fixingT[i][j], num_index2_coup_cf[i], index2_d[i], index2_t[i]);
            tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
            fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice;
        }
        for (j = num_realfixing[i]; j < num_fixing[i]; j++) fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice;
        fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
        fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
}

void CMSSpreadNonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index2_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate[i] = vector<double>(num_estm[i]);
            index2_adjrate[i] = vector<double>(num_estm[i]);
            for (j = 0; j < num_estm[i]; j++)
            {
                level1 = 0.0;
                level2 = 0.0;
                for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
                index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
                index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], couponindex_fixingT[i][j], num_index2_coup_cf[i], index2_d[i], index2_t[i]);
                tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
                fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
            }
            fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
}
/*

void NonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double>
lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double &fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice=0.0;
    double nonfixedrateincount=0.0, level1, tmpprice=0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf,0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf);
    bool tmpflag=false;
    for(i=0;i<num_remained_couponleg_cf;i++)
    {
        nonfixedrateincount=0.0;
        index1_fra[i]=vector<double>(num_estm[i]);
        index1_adjrate[i]=vector<double>(num_estm[i]);
        for(j=0;j<num_estm[i];j++)
        {
            level1=0.0;
            for(k=0;k<num_index1_coup_cf[i];k++) level1=level1+index1_tau[i][j][k]*exp(-zc(index1_paydate[i][j][k])*index1_payT[i][j][k]);
            index1_fra[i][j]=(exp(-zc(index1_startdate[i][j])*index1_startT[i][j])-exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i]-1])*index1_T[i][j]))/level1;
            index1_adjrate[i][j]=index1_fra[i][j]+convadj(index1_fra[i][j],index1_vol[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
            //tmpprice=ls[i]*(SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],lk[i])-SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],_ralow_bdry[i]))+hs[i]*(SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],hk[i])-SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],_rahigh_bdry[i]));
            tmpprice=ls[i]*(SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],lk[i])-SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],_ralow_bdry[i]))+hs[i]*(SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],hk[i])-SpreadCall(-_couponlegindex2_mult[i]*index2_adjrate[i][j],_couponlegindex1_mult[i]*index1_adjrate[i][j],index2_vol[i][j],index1_vol[i][j],0.0,0.0,input_corr,1.0,couponindex_fixingT[i][j],_rahigh_bdry[i]));
            fixedlegnonfixedcf[i]=fixedlegnonfixedcf[i]+tmpprice;
        }
        for(j=num_realfixing[i];j<num_fixing[i];j++) fixedlegnonfixedcf[i]=fixedlegnonfixedcf[i]+tmpprice;
        fixedlegnonfixedcf[i]=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(fixedlegnonfixedcf[i]/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
        fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];
    }
    fixedlegnonfixedprice=fixedlegnonfixedprice/ondf/settlement_date_df;
}
*/
void NonFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag
    , vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate_ld(num_remained_couponleg_cf), index1_adjrate_lu(num_remained_couponleg_cf), index1_adjrate_rd(num_remained_couponleg_cf), index1_adjrate_ru(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate_ld[i] = vector<double>(num_estm[i]);
            index1_adjrate_lu[i] = vector<double>(num_estm[i]);
            index1_adjrate_rd[i] = vector<double>(num_estm[i]);
            index1_adjrate_ru[i] = vector<double>(num_estm[i]);
            if (_ralow_bdry[i] <= 0.0)
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ld[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_lu[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_rd[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ru[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);

                    tmpprice = ls[i] * (_ralow_bdry[i] - lk[i]) + hs[i] * (PlainOptionlet("CALL", index1_adjrate_rd[i][j], hk[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]) - PlainOptionlet("CALL", index1_adjrate_rd[i][j], _rahigh_bdry[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            else
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ld[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_lu[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_rd[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ru[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);

                    tmpprice = ls[i] * (PlainOptionlet("CALL",index1_adjrate_ld[i][j],lk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ld[i][j]) - PlainOptionlet("CALL",index1_adjrate_lu[i][j],_ralow_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_lu[i][j])) + hs[i] * (PlainOptionlet("CALL",index1_adjrate_rd[i][j],hk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]) - PlainOptionlet("CALL",index1_adjrate_rd[i][j],_rahigh_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
}

void NonFixedAccrualZeroPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag
    , vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate_ld(num_remained_couponleg_cf), index1_adjrate_lu(num_remained_couponleg_cf), index1_adjrate_rd(num_remained_couponleg_cf), index1_adjrate_ru(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate_ld[i] = vector<double>(num_estm[i]);
            index1_adjrate_lu[i] = vector<double>(num_estm[i]);
            index1_adjrate_rd[i] = vector<double>(num_estm[i]);
            index1_adjrate_ru[i] = vector<double>(num_estm[i]);
            if (_ralow_bdry[i] <= 0.0)
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ld[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_lu[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_rd[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ru[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);

                    tmpprice = ls[i] * (_ralow_bdry[i] - lk[i]) + hs[i] * (PlainOptionlet("CALL", index1_adjrate_rd[i][j], hk[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]) - PlainOptionlet("CALL", index1_adjrate_rd[i][j], _rahigh_bdry[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            else
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ld[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_lu[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_rd[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ru[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);

                    tmpprice = ls[i] * (PlainOptionlet("CALL",index1_adjrate_ld[i][j],lk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ld[i][j]) - PlainOptionlet("CALL",index1_adjrate_lu[i][j],_ralow_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_lu[i][j])) + hs[i] * (PlainOptionlet("CALL",index1_adjrate_rd[i][j],hk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]) - PlainOptionlet("CALL",index1_adjrate_rd[i][j],_rahigh_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
                //fixedlegnonfixedcf[i]=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*couponlegtau[i]*(fixedlegnonfixedcf[i]/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
            fixedlegnonfixedcf[i] = _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = exp(-zc(_couponleg_paydate[num_remained_couponleg_cf - 1]) * couponlegT[num_remained_couponleg_cf - 1]) * fixedlegnonfixedprice / ondf / settlement_date_df;
}



void CMSSpreadNonFixedAccrualPrice_Repli(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> _rahigh_bdryoh, vector<bool> _rahigh_bdryin_flag, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> _ralow_bdryoh, vector<bool> _ralow_bdryin_flag, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0, Dtp;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index2_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate[i] = vector<double>(num_estm[i]);
            index2_adjrate[i] = vector<double>(num_estm[i]);
            Dtp = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]);
            for (j = 0; j < num_estm[i]; j++)
            {
                level1 = 0.0;
                level2 = 0.0;
                for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
                index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
                index1_adjrate[i][j] = index1_fra[i][j] + convadj_strd(index1_fra[i][j], index1_vol[i][j], (couponlegT[i] - index1_startT[i][j]) / index1_tau[i][j][0], couponindex_fixingT[i][j], index1_t[i], num_index1_coup_cf[i], level1, Dtp);
                index2_adjrate[i][j] = index2_fra[i][j] + convadj_strd(index2_fra[i][j], index2_vol[i][j], (couponlegT[i] - index2_startT[i][j]) / index2_tau[i][j][0], couponindex_fixingT[i][j], index2_t[i], num_index2_coup_cf[i], level2, Dtp);
                tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
                fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
            }
            fixedlegnonfixedcf[i] = Dtp * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
}
///////////////////////////////////////////////////////////////

void SwapRatewithVoladj(CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_estm, vector<int> num_index1_coup_cf,
    vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol, double input_corr, vector<double> index1_t,
    vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T,
    vector<vector<vector<double>>> index1_payT, vector<double>, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& swaprate)
{
    int i, j, k;

    double level1, Dtp;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf);

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        index1_fra[i] = vector<double>(num_estm[i]);
        index1_adjrate[i] = vector<double>(num_estm[i]);
        Dtp = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]);
        for (j = 0; j < num_estm[i]; j++)
        {
            level1 = 0.0;
            for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
            index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
            index1_adjrate[i][j] = index1_fra[i][j] + convadj_strd(index1_fra[i][j], index1_vol[i][j], (couponlegT[i] - index1_startT[i][j]) / index1_tau[i][j][0], couponindex_fixingT[i][j], index1_t[i], num_index1_coup_cf[i], level1, Dtp);
        }
    }

}

void PlainCMSSpreadRangeAccrualSwap(CCurrency crcy, CDate today, CDate settlement_date, double notional, CInterpolation zc, int num_fundingleg_cf, int fundingleg_cf_starti, vector<CDate> fixing_date, vector<CDate> fundingleg_calc_startdate, vector<CDate> fundingleg_calc_enddate, vector<CDate> fundingleg_paydate, string fundingleg_dcb, vector<double> fundingleg_notional, vector<double> fundingleg_mult, vector<double> fundingleg_spread, int couponleg_cf_starti, vector<CDate> couponleg_calc_startdate, vector<CDate> couponleg_calc_enddate, vector<CDate> couponleg_paydate, string couponleg_dcb, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> rahigh_bdryoh, vector<int> rahigh_bdryoh_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<double> ralow_bdryoh, vector<int> ralow_bdryoh_flag
    , vector<double> fixinghistory_rate, double& floatinglegprice, double& fixedlegprice, double& atmswaprate, double& swapprice)
{
    int i;
    string dcb = "ACT/365";
    vector<double> fundingleg_df(num_fundingleg_cf);
    floatinglegprice = 0.0;
    for (i = fundingleg_cf_starti; i < num_fundingleg_cf; i++)
    {
        fundingleg_df[i] = exp(-zc(fundingleg_paydate[i]) * cvg(today, fundingleg_paydate[i], dcb));
        if (fixinghistory_rate[i] > MinFixingRate) floatinglegprice = floatinglegprice + fundingleg_notional[i] * (fundingleg_mult[i] * fixinghistory_rate[i] + fundingleg_spread[i]) * cvg(fundingleg_calc_startdate[i], fundingleg_calc_enddate[i], fundingleg_dcb) * fundingleg_df[i];
        if (today <= fixing_date[i]) break;
    }
    if (today == fundingleg_paydate[num_fundingleg_cf - 1]) floatinglegprice = 0.0;
    //  if(i<num_fundingleg_cf)
}

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional
    , vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        nonfixedrateincount = 0.0;
        index1_fra[i] = vector<double>(num_estm[i]);
        index2_fra[i] = vector<double>(num_estm[i]);
        index1_adjrate[i] = vector<double>(num_estm[i]);
        index2_adjrate[i] = vector<double>(num_estm[i]);
        for (j = 0; j < num_estm[i]; j++)
        {
            level1 = 0.0;
            level2 = 0.0;
            for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
            for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
            index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
            index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
            index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
            index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], couponindex_fixingT[i][j], num_index2_coup_cf[i], index2_d[i], index2_t[i]);
            tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
            fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice;
        }
        for (j = num_realfixing[i]; j < num_fixing[i]; j++) fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice;
        fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
        fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, int d, vector<int> num_fixing, vector<int> num_realfixing, vector<int> num_estm, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT
    , vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        nonfixedrateincount = 0.0;
        index1_fra[i] = vector<double>(num_estm[i]);
        index2_fra[i] = vector<double>(num_estm[i]);
        index1_adjrate[i] = vector<double>(num_estm[i]);
        index2_adjrate[i] = vector<double>(num_estm[i]);
        for (j = 0; j < num_estm[i]; j = j + d)
        {
            level1 = 0.0;
            level2 = 0.0;
            for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
            for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
            index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
            index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
            index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
            index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], couponindex_fixingT[i][j], num_index2_coup_cf[i], index2_d[i], index2_t[i]);
            tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
            fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + d * tmpprice;
        }
        fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + (num_estm[i] - j) * tmpprice;
        for (j = num_realfixing[i]; j < num_fixing[i]; j = j + d) fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + d * tmpprice;
        fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + (num_fixing[i] - j) * tmpprice;
        fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
        fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

//void PlainSwapFundinglegNonFixedPricexy(CDate today, double callnotice_t, double ondf, double settlement_date_df, CInterpolation zc, int couponleg_cf_starti, int num_remained_fundingleg_cf, vector<double> _fundingleg_notional, vector<double> _fundingleg_mult, vector<double> _fundingleg_spread, CDate fundingleg_calc_startdate0, double fundinglegcalcstartT0, vector<CDate> _fundingleg_paydate, vector<double> fundinglegtau, vector<double> fundinglegT, double &floatinglegnonfixedprice)
void PlainSwapFundinglegNonFixedPricexy(CDate callnotice_date, double callnotice_t, double ondf, CInterpolation zc, int couponleg_cf_starti, int num_remained_fundingleg_cf, vector<double> _fundingleg_notional, vector<double> _fundingleg_mult, vector<double> _fundingleg_spread, CDate fundingleg_calc_startdate0, double fundinglegcalcstartT0, vector<CDate> _fundingleg_paydate, vector<double> fundinglegtau, vector<double> fundinglegT, double& floatinglegnonfixedprice)
{
    int i;
    double imsi_0 = 0.0;
    double imsi_1 = 0.0;

    double imsi_2 = 0.0;
    double imsi_3 = 0.0;

    double imsi_4 = 0.0;
    double imsi_5 = 0.0;

    string dcb = "ACT/365";
    vector<double> fundingleg_df(num_remained_fundingleg_cf);
    floatinglegnonfixedprice = 0.0;
    if (num_remained_fundingleg_cf > 0)
    {
        if (callnotice_date < _fundingleg_paydate[num_remained_fundingleg_cf - 1])
        {
            for (i = couponleg_cf_starti; i < num_remained_fundingleg_cf; i++)
            {
                imsi_0 = zc(_fundingleg_paydate[i]);
                imsi_1 = (fundinglegT[i] - callnotice_t);

                fundingleg_df[i] = exp(-zc(_fundingleg_paydate[i]) * (fundinglegT[i] - callnotice_t));

                floatinglegnonfixedprice = floatinglegnonfixedprice + _fundingleg_notional[i] * _fundingleg_spread[i] * fundinglegtau[i] * fundingleg_df[i];
            }
            if (num_remained_fundingleg_cf > 0)
            {
                imsi_2 = zc(fundingleg_calc_startdate0);
                imsi_3 = (fundinglegcalcstartT0 - callnotice_t);

                // org
                imsi_4 = _fundingleg_notional[0] * exp(-zc(fundingleg_calc_startdate0) * (fundinglegcalcstartT0 - callnotice_t));
                // renew
                imsi_4 = _fundingleg_notional[0] * exp(-zc(fundingleg_calc_startdate0) * (fundinglegcalcstartT0));
                imsi_4 = _fundingleg_notional[num_remained_fundingleg_cf - 1] * fundingleg_df[num_remained_fundingleg_cf - 1];
                // org
                floatinglegnonfixedprice = floatinglegnonfixedprice + (_fundingleg_notional[0] * exp(-zc(fundingleg_calc_startdate0) * (fundinglegcalcstartT0 - callnotice_t)) - _fundingleg_notional[num_remained_fundingleg_cf - 1] * fundingleg_df[num_remained_fundingleg_cf - 1]);
            }
            //floatinglegnonfixedprice=floatinglegnonfixedprice/ondf/settlement_date_df;
            floatinglegnonfixedprice = floatinglegnonfixedprice / ondf;
        }
    }
}

void CMSSpreadNonFixedAccrualPricexy(int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT, vector<vector<vector<double>>> index2_payT
    , vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index2_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate[i] = vector<double>(num_estm[i]);
            index2_adjrate[i] = vector<double>(num_estm[i]);
            for (j = 0; j < num_estm[i]; j++)
            {
                level1 = 0.0;
                level2 = 0.0;
                for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * index2_payT[i][j][k]);
                index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * index2_startT[i][j]) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * index2_T[i][j])) / level2;
                index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], couponindex_fixingT[i][j], num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], couponindex_fixingT[i][j], num_index2_coup_cf[i], index2_d[i], index2_t[i]);
                tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j], _rahigh_bdry[i]));
                fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
            }
            fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void CMSSpreadNonFixedAccrualPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<int> num_index2_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<CDate>>index2_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<vector<CDate>>> index2_paydate, vector<vector<double>> index1_vol, vector<vector<double>> index2_vol, double input_corr, vector<double> index1_t, vector<double> index2_t, vector<vector<double>> index1_d, vector<vector<double>> index2_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<vector<double>>> index2_tau, vector<vector<double>> index1_startT, vector<vector<double>> index2_startT, vector<vector<double>> index1_T, vector<vector<double>> index2_T, vector<vector<vector<double>>> index1_payT
    , vector<vector<vector<double>>> index2_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, level2, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index2_fra(num_remained_couponleg_cf), index1_adjrate(num_remained_couponleg_cf), index2_adjrate(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index2_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate[i] = vector<double>(num_estm[i]);
            index2_adjrate[i] = vector<double>(num_estm[i]);
            for (j = 0; j < num_estm[i]; j++)
            {
                level1 = 0.0;
                level2 = 0.0;
                for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * (index1_payT[i][j][k] - callnotice_t));
                for (k = 0; k < num_index2_coup_cf[i]; k++) level2 = level2 + index2_tau[i][j][k] * exp(-zc(index2_paydate[i][j][k]) * (index2_payT[i][j][k] - callnotice_t));
                index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * (index1_startT[i][j] - callnotice_t)) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * (index1_T[i][j] - callnotice_t))) / level1;
                index2_fra[i][j] = (exp(-zc(index2_startdate[i][j]) * (index2_startT[i][j] - callnotice_t)) - exp(-zc(index2_paydate[i][j][num_index2_coup_cf[i] - 1]) * (index2_T[i][j] - callnotice_t))) / level2;
                index1_adjrate[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                index2_adjrate[i][j] = index2_fra[i][j] + convadj(index2_fra[i][j], index2_vol[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index2_coup_cf[i], index2_d[i], index2_t[i]);
                tmpprice = ls[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j] - callnotice_t, lk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j] - callnotice_t, _ralow_bdry[i])) + hs[i] * (SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j] - callnotice_t, hk[i]) - SpreadCall(-_couponlegindex2_mult[i] * index2_adjrate[i][j], _couponlegindex1_mult[i] * index1_adjrate[i][j], index2_vol[i][j], index1_vol[i][j], 0.0, 0.0, input_corr, 1.0, couponindex_fixingT[i][j] - callnotice_t, _rahigh_bdry[i]));
                fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
            }
            fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * (couponlegT[i] - callnotice_t)) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void NonFixedPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i;
    fixedlegnonfixedprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);

    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * (couponlegT[i] - callnotice_t)) * _couponleg_notional[i] * couponlegtau[i] * _couponleg_couponrate[i];
        fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void NonFixedZeroPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i;

    double imsi_0 = 0.0;
    double imsi_1 = 0.0;

    fixedlegnonfixedprice = 0.0;
    //vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf,0.0);

    vector<double> fixedlegnonfixedcf(num_couponleg_cf, 0.0);

    bool tmpflag = false;

    double fullzerocoupon = 0.0, partialzerocoupon = 0.0;

    for (i = 0; i < num_couponleg_cf; i++)
    {
        imsi_0 = couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];

        fullzerocoupon = fullzerocoupon + couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];

        //fixedlegnonfixedcf[i]=exp(-zc(_couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*_couponleg_notional[i]*couponlegtau[i]*_couponleg_couponrate[i];
        //fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];

        //fixedlegnonfixedcf[i]=exp(-zc(couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*couponleg_notional[i]*couponlegtau[i]*couponleg_couponrate[i];
        //fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];

    }
    fullzerocoupon = fullzerocoupon * exp(-zc(couponleg_paydate[num_couponleg_cf - 1]) * (couponlegT[num_couponleg_cf - 1] - callnotice_t));
    //fullzerocoupon = fullzerocoupon*exp(-zc(couponleg_paydate[num_couponleg_cf - 1])*(couponlegT[num_couponleg_cf - 1]));
    //fullzerocoupon = fullzerocoupon;

    for (i = 0; i < startk; i++)
    {
        imsi_1 = couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];

        partialzerocoupon = partialzerocoupon + couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];

        //fixedlegnonfixedcf[i]=exp(-zc(_couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*_couponleg_notional[i]*couponlegtau[i]*_couponleg_couponrate[i];
        //fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];

        //fixedlegnonfixedcf[i]=exp(-zc(couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*couponleg_notional[i]*couponlegtau[i]*couponleg_couponrate[i];
        //fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];

    }

    //partialzerocoupon = partialzerocoupon*exp(zc(couponleg_paydate[startk - 1])*(couponlegT[startk - 1] - callnotice_t));
    // [[[ hokim 20191104 zc -> -zc
    //partialzerocoupon = partialzerocoupon*exp(-zc(couponleg_paydate[startk - 1])*(couponlegT[startk - 1] - callnotice_t));
    // ]]]

    fixedlegnonfixedprice = fullzerocoupon - partialzerocoupon;
    //fixedlegnonfixedprice = fullzerocoupon;


    /*for(i=startk;i<num_couponleg_cf;i++)
    {
        fixedlegnonfixedcf[i]=exp(-zc(couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*couponleg_notional[i]*couponlegtau[i]*couponleg_couponrate[i];
        fixedlegnonfixedprice=fixedlegnonfixedprice+fixedlegnonfixedcf[i];
    }
    */

    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}


void NonFixedZeroPricexy2(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice, double& fixedlegpartialnonfixedprice)
{
    int i;
    fixedlegnonfixedprice = 0.0;
    //vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf,0.0);

    bool tmpflag = false;

    double fullzerocoupon = 0.0, partialzerocoupon = 0.0;

    for (i = 0; i < num_couponleg_cf; i++)
    {
        fullzerocoupon = fullzerocoupon + couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];
    }
    fullzerocoupon = fullzerocoupon * exp(-zc(couponleg_paydate[num_couponleg_cf - 1]) * (couponlegT[num_couponleg_cf - 1] - callnotice_t));

    for (i = 0; i < startk; i++)
    {
        partialzerocoupon = partialzerocoupon + couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];
    }

    fixedlegpartialnonfixedprice = partialzerocoupon / ondf;

    fixedlegnonfixedprice = fullzerocoupon / ondf;
}

// 2018-08-08 created by YJSEO : compound coupon rate
void NonFixedZeroPricexy3(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_couponleg_cf, vector<double> couponleg_notional, vector<double> couponleg_couponrate, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    //int i;
    fixedlegnonfixedprice = 0.0;

    bool tmpflag = false;

    double fullzerocoupon = 0.0, partialzerocoupon = 0.0;

    /*for(i=0;i<num_couponleg_cf;i++)
    {
        fullzerocoupon=fullzerocoupon+couponleg_notional[i]*couponlegtau[i]*couponleg_couponrate[i];
    }*/
    fullzerocoupon = couponleg_notional[num_couponleg_cf - 1] * (pow((1 + couponleg_couponrate[num_couponleg_cf - 1]), couponlegtau[num_couponleg_cf - 1]) - 1);
    fullzerocoupon = fullzerocoupon * exp(-zc(couponleg_paydate[num_couponleg_cf - 1]) * (couponlegT[num_couponleg_cf - 1] - callnotice_t));

    /*for(i=0;i<startk;i++)
    {
        partialzerocoupon=partialzerocoupon+couponleg_notional[i]*couponlegtau[i]*couponleg_couponrate[i];
    }*/

    partialzerocoupon = couponleg_notional[startk - 1] * (pow((1 + couponleg_couponrate[startk - 1]), couponlegtau[startk - 1]) - 1);
    fixedlegnonfixedprice = fullzerocoupon - partialzerocoupon;

    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void NonFixedAccrualPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry, vector<double> lk
    , vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate_ld(num_remained_couponleg_cf), index1_adjrate_lu(num_remained_couponleg_cf), index1_adjrate_rd(num_remained_couponleg_cf), index1_adjrate_ru(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate_ld[i] = vector<double>(num_estm[i]);
            index1_adjrate_lu[i] = vector<double>(num_estm[i]);
            index1_adjrate_rd[i] = vector<double>(num_estm[i]);
            index1_adjrate_ru[i] = vector<double>(num_estm[i]);
            if (_ralow_bdry[i] <= 0.0)
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * (index1_payT[i][j][k] - callnotice_t));
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * (index1_startT[i][j] - callnotice_t)) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * (index1_T[i][j] - callnotice_t))) / level1;

                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ld[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_lu[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_rd[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ru[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);

                    tmpprice = ls[i] * (_ralow_bdry[i] - lk[i]) + hs[i] * (PlainOptionlet("CALL", index1_adjrate_rd[i][j], hk[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]) - PlainOptionlet("CALL", index1_adjrate_rd[i][j], _rahigh_bdry[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]));
                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            else
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ld[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_lu[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_rd[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ru[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);

                    tmpprice = ls[i] * (PlainOptionlet("CALL",index1_adjrate_ld[i][j],lk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ld[i][j]) - PlainOptionlet("CALL",index1_adjrate_lu[i][j],_ralow_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_lu[i][j])) + hs[i] * (PlainOptionlet("CALL",index1_adjrate_rd[i][j],hk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]) - PlainOptionlet("CALL",index1_adjrate_rd[i][j],_rahigh_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            fixedlegnonfixedcf[i] = exp(-zc(_couponleg_paydate[i]) * (couponlegT[i] - callnotice_t)) * _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf;
}

void NonFixedAccrualZeroPricexy(double callnotice_t, int startk, double ondf, CInterpolation zc, int num_remained_couponleg_cf, vector<bool> _ra_flag, vector<int> num_fixing, vector<int> num_estm, vector<vector<int>> num_estm_period, vector<int> num_index1_coup_cf, vector<vector<CDate>> index1_startdate, vector<vector<vector<CDate>>> index1_paydate, vector<vector<double>> index1_vol_ld, vector<vector<double>> index1_vol_lu, vector<vector<double>> index1_vol_rd, vector<vector<double>> index1_vol_ru, vector<double> index1_t, vector<vector<double>> index1_d, vector<vector<double>> couponindex_fixingT, vector<vector<vector<double>>> index1_tau, vector<vector<double>> index1_startT, vector<vector<double>> index1_T, vector<vector<vector<double>>> index1_payT, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<double> hk, vector<double> hs, vector<double> _ralow_bdry
    , vector<double> lk, vector<double> ls, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double dP, double& fixedlegnonfixedprice)
{
    int i, j, k;
    fixedlegnonfixedprice = 0.0;
    double nonfixedrateincount = 0.0, level1, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<vector<double>> index1_fra(num_remained_couponleg_cf), index1_adjrate_ld(num_remained_couponleg_cf), index1_adjrate_lu(num_remained_couponleg_cf), index1_adjrate_rd(num_remained_couponleg_cf), index1_adjrate_ru(num_remained_couponleg_cf);
    bool tmpflag = false;
    for (i = startk; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            nonfixedrateincount = 0.0;
            index1_fra[i] = vector<double>(num_estm[i]);
            index1_adjrate_ld[i] = vector<double>(num_estm[i]);
            index1_adjrate_lu[i] = vector<double>(num_estm[i]);
            index1_adjrate_rd[i] = vector<double>(num_estm[i]);
            index1_adjrate_ru[i] = vector<double>(num_estm[i]);
            if (_ralow_bdry[i] <= 0.0)
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * (index1_payT[i][j][k] - callnotice_t));
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * (index1_startT[i][j] - callnotice_t)) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * (index1_T[i][j] - callnotice_t))) / level1;

                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ld[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_lu[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_rd[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j], index1_vol_ru[i][j], (couponindex_fixingT[i][j] - callnotice_t), num_index1_coup_cf[i], index1_d[i], index1_t[i]);

                    tmpprice = ls[i] * (_ralow_bdry[i] - lk[i]) + hs[i] * (PlainOptionlet("CALL", index1_adjrate_rd[i][j], hk[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]) - PlainOptionlet("CALL", index1_adjrate_rd[i][j], _rahigh_bdry[i], 1.0, 1.0, couponindex_fixingT[i][j], index1_vol_ru[i][j]));
                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
            else
            {
                for (j = 0; j < num_estm[i]; j++)
                {
                    level1 = 0.0;
                    for (k = 0; k < num_index1_coup_cf[i]; k++) level1 = level1 + index1_tau[i][j][k] * exp(-zc(index1_paydate[i][j][k]) * index1_payT[i][j][k]);
                    index1_fra[i][j] = (exp(-zc(index1_startdate[i][j]) * index1_startT[i][j]) - exp(-zc(index1_paydate[i][j][num_index1_coup_cf[i] - 1]) * index1_T[i][j])) / level1;
                    index1_adjrate_ld[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ld[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_lu[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_lu[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_rd[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_rd[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);
                    index1_adjrate_ru[i][j] = index1_fra[i][j] + convadj(index1_fra[i][j],index1_vol_ru[i][j],couponindex_fixingT[i][j],num_index1_coup_cf[i],index1_d[i],index1_t[i]);

                    tmpprice = ls[i] * (PlainOptionlet("CALL",index1_adjrate_ld[i][j],lk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ld[i][j]) - PlainOptionlet("CALL",index1_adjrate_lu[i][j],_ralow_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_lu[i][j])) + hs[i] * (PlainOptionlet("CALL",index1_adjrate_rd[i][j],hk[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]) - PlainOptionlet("CALL",index1_adjrate_rd[i][j],_rahigh_bdry[i],1.0,1.0,couponindex_fixingT[i][j],index1_vol_ru[i][j]));

                    fixedlegnonfixedcf[i] = fixedlegnonfixedcf[i] + tmpprice * num_estm_period[i][j];
                }
            }
                //fixedlegnonfixedcf[i]=exp(-zc(_couponleg_paydate[i])*(couponlegT[i]-callnotice_t))*_couponleg_notional[i]*couponlegtau[i]*(fixedlegnonfixedcf[i]/num_fixing[i]*_couponleg_couponrate[i]+_couponleg_spread[i]);
            fixedlegnonfixedcf[i] = _couponleg_notional[i] * couponlegtau[i] * (fixedlegnonfixedcf[i] / num_fixing[i] * _couponleg_couponrate[i] + _couponleg_spread[i]);
            fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];
        }
    }
    fixedlegnonfixedprice = exp(-zc(_couponleg_paydate[num_remained_couponleg_cf - 1]) * (couponlegT[num_remained_couponleg_cf - 1] - callnotice_t)) * fixedlegnonfixedprice / ondf;
    fixedlegnonfixedprice = min(fixedlegnonfixedprice, INT_MAX * dP);
}



void QuantoFloatingRateCMSFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<double> _couponfloatingfixinghistory_rate, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedlegaccrualfixedprice, vector<double>& fixedlegaccrualfixednumber)
{
    int i, j;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            fixedlegaccrualfixednumber[i] = fixedrateincount;
            if (num_fixed[i] == num_fixing[i]) fixedlegaccrualfixedprice[i] = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * min(max((_couponlegindex1_mult[i] * _couponfloatingfixinghistory_rate[i] + _couponleg_couponrate[i] + _couponleg_spread[i]) * fixedlegaccrualfixednumber[i] / num_fixing[i], _floorrates[i]), _caprates[i]) / ondf / settlement_date_df;
        }
    }
}

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
)
{
    int
        i,
        j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0,
        tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i],
                            _rahigh_bdry[i],
                            _ralow_bdry[i],
                            _rahigh_bdryin_flag[i],
                            _ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
            }

            tempprice
                = exp(-zc(_couponleg_paydate[i]) * couponlegT[i])
                * _couponleg_notional[i]
                * _couponleg_couponrate[i]
                * couponlegtau[i]
                * fixedrateincount
                / num_fixing[i];

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }

    fixedlegaccrualfixedprice
        = fixedlegaccrualfixedprice
        / ondf
        / settlement_date_df;
}

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
)
{
    int
        i,
        j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0,
        tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i],
                            _rahigh_bdry[i],
                            _ralow_bdry[i],
                            _rahigh_bdryin_flag[i],
                            _ralow_bdryin_flag[i]
                        )
                        &&
                        rangeincheck
                        (
                            fixedrate1[i][j],
                            _caprates[i],
                            _floorrates[i],
                            _rahigh_bdryin_flag[i],
                            _ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
            }

            tempprice
                = exp(-zc(_couponleg_paydate[num_remained_couponleg_cf - 1]) * couponlegT[num_remained_couponleg_cf - 1])
                * _couponleg_notional[i]
                * _couponleg_couponrate[i]
                * couponlegtau[i]
                * fixedrateincount
                / num_fixing[i];

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }

    fixedlegaccrualfixedprice
        = fixedlegaccrualfixedprice
        / ondf
        / settlement_date_df;
}

void LIBORSpreadDualFixedAccrualPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(fixedrate3[i][j], _caprates[i], _floorrates[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}


void LIBORSpreadDualFixedAccrualPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex3_mult[i] * fixedrate3[i][j], _caprates[i], _floorrates[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}


// 0<= CMT10Y - CMY5Y, USDCMS10Y <=5.5%, CMT10Y <=5.5% - 2017/02/02 by yeojin
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
)
{
    int
        i
        , j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0
        , tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                //if(rangeincheck(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponleg_spread[i],_rahigh_bdry[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i]*fixedrate3[i][j],_rahigh_bdry[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i])) fixedrateincount=fixedrateincount+1.0;

                if
                    (
                        rangeincheck
                        (
                            _couponlegindex1_mult[i] * fixedrate1[i][j]
                            + _couponlegindex2_mult[i] * fixedrate2[i][j]
                            + _couponleg_spread[i], _rahigh_bdry[i]
                            , _ralow_bdry[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                        )
                        &&
                        rangeincheck
                        (
                            _couponlegindex2_mult[i] * fixedrate2[i][j]
                            , _rahigh_bdry2[i]
                            , _ralow_bdry2[i]
                            , _rahigh_bdry2in_flag[i]
                            , _ralow_bdry2in_flag[i]
                        )
                        &&
                        rangeincheck
                        (
                            _couponlegindex3_mult[i] * fixedrate3[i][j]
                            , _rahigh_bdry3[i]
                            , _ralow_bdry3[i]
                            , _rahigh_bdry3in_flag[i]
                            , _ralow_bdry3in_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
            }

            tempprice =
                exp
                (
                    -zc(_couponleg_paydate[i])
                    * couponlegT[i]
                )
                * _couponleg_notional[i]
                * _couponleg_couponrate[i]
                * couponlegtau[i]
                * fixedrateincount
                / num_fixing[i];

            fixedlegaccrualfixedprice =
                fixedlegaccrualfixedprice
                + tempprice;
        }
    }
    fixedlegaccrualfixedprice =
        fixedlegaccrualfixedprice
        / ondf
        / settlement_date_df;
}

void BasisDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i] * fixedrate3[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

// 2019-03-28 CMS1-CMS2 > 0, CMS1-CMS3 >-0.1
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
)
{
    int i
        , j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0
        , tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (
                    rangeincheck
                    (
                        _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i]
                        , _rahigh_bdry[i]
                        , _ralow_bdry[i]
                        , _rahigh_bdryin_flag[i]
                        , _ralow_bdryin_flag[i]
                    )
                    &&
                    rangeincheck
                    (
                        _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex3_mult[i] * fixedrate3[i][j]
                        , _rahigh_bdry2[i]
                        , _ralow_bdry2[i]
                        , _rahigh_bdry2in_flag[i]
                        , _ralow_bdry2in_flag[i]
                    )
                    )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
            }

            tempprice =
                exp(-zc(_couponleg_paydate[i]) * couponlegT[i])
                * _couponleg_notional[i]
                * _couponleg_couponrate[i]
                * couponlegtau[i]
                * fixedrateincount
                / num_fixing[i];

            fixedlegaccrualfixedprice =
                fixedlegaccrualfixedprice
                + tempprice;
        }
    }

    fixedlegaccrualfixedprice =
        fixedlegaccrualfixedprice
        / ondf
        / settlement_date_df;
}
//generated on 2018-03-08
void BasisTripleFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _rahigh_bdry1, vector<bool> _rahigh_bdryin1_flag, vector<double> _ralow_bdry1, vector<bool> _ralow_bdryin1_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i] * fixedrate3[i][j], _rahigh_bdry1[i], _ralow_bdry1[i], _rahigh_bdryin1_flag[i], _ralow_bdryin1_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / (num_fixing[i] + num_fixed[i]);
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

// 2019-05-10 modified yj : 경과이자 제대로 계산
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
)
{
    int
        i
        , j;

    fixedlegaccrualfixedprice = 0.0;

    double
        fixedrateincount = 0.0
        , tempprice = 0.0;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;

        if (_ra_flag[i])
        {
            if (i < changeCondIdx)
            {
                for (j = 0; j < num_fixed[i]; j++)
                {
                    if
                        (
                            rangeincheck
                            (
                                _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j]
                                , _rahigh_bdry[i]
                                , _ralow_bdry[i]
                                , _rahigh_bdryin_flag[i]
                                , _ralow_bdryin_flag[i]
                            )
                            && rangeincheck
                            (
                                _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex3_mult[i] * fixedrate3[i][j]
                                , _rahigh_bdry2[i]
                                , _ralow_bdry2[i]
                                , _rahigh_bdry2in_flag[i]
                                , _ralow_bdry2in_flag[i]
                            )
                            && rangeincheck
                            (
                                _couponlegindex1_mult[i] * fixedrate1[i][j]
                                , _rahigh_bdry3[i]
                                , _ralow_bdry3[i]
                                , _rahigh_bdry3in_flag[i]
                                , _ralow_bdry3in_flag[i]
                            )
                            )
                    {
                        fixedrateincount = fixedrateincount + 1.0;
                    }
                }

                tempprice
                    = _couponleg_couponrate[i]
                    * fixedrateincount
                    / num_fixing[i]
                    * _couponleg_notional[i]
                    * couponlegtau[i];

            }
            else
            {
                for (j = 0; j < num_fixed[i]; j++)
                {
                    if
                    (
                        rangeincheck(
                            _couponlegindex4_mult[i] * fixedrate4[i][j] + _couponlegindex5_mult[i] * fixedrate5[i][j]
                            , _rahigh_bdry[i]
                            , _ralow_bdry[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                            )
                        && rangeincheck(
                            _couponlegindex4_mult[i] * fixedrate4[i][j] + _couponlegindex3_mult[i] * fixedrate3[i][j]
                            , _rahigh_bdry2[i]
                            , _ralow_bdry2[i], _rahigh_bdry2in_flag[i], _ralow_bdry2in_flag[i]
                            )
                        && rangeincheck(
                            _couponlegindex4_mult[i] * fixedrate4[i][j]
                            , _rahigh_bdry3[i]
                            , _ralow_bdry3[i]
                            , _rahigh_bdry3in_flag[i]
                            , _ralow_bdry3in_flag[i]
                        )
                            )
                    {
                        fixedrateincount = fixedrateincount + 1.0;
                    }
                }
            //num_fixed[i-1]-1 : 매 이자계산 기간 시작일로부터 1 영업일전
            tempprice
                = max(
                    fixedrate4[i - 1][num_fixed[i - 1] - 1] + _couponleg_couponrate[i]
                    ,0.0
                    )
                * fixedrateincount / num_fixing[i]
                * _couponleg_notional[i]
                * couponlegtau[i];
            }

                //tempprice = _couponleg_notional[i] * max(_couponleg_couponrate[i] + _couponleg_spread[i],0.0) * couponlegtau[i] * fixedrateincount / num_fixing[i];           //평가일이 거래후 1년이상 일때 체크

                fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice + tempprice;

        }
    }
    fixedlegaccrualfixedprice
        = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

//void BasisTripleFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double &fixedlegaccrualfixedprice)
//{
//  int i, j;
//  fixedlegaccrualfixedprice=0.0;
//  double fixedrateincount=0.0, tempprice=0.0;
//  for(i=0;i<num_remained_couponleg_cf;i++)
//  {
//      fixedrateincount=0.0;
//      if(_ra_flag[i])
//      {
//          for(j=0;j<num_fixed[i];j++)
//          {
//              if(rangeincheck(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponleg_spread[i],_rahigh_bdry[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i]*fixedrate3[i][j],_caprates[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i])) fixedrateincount=fixedrateincount+1.0;
//          }
//          tempprice=_couponleg_notional[i]*_couponleg_couponrate[i]*couponlegtau[i]*fixedrateincount/(num_fixing[i]+num_fixed[i]);
//          fixedlegaccrualfixedprice=fixedlegaccrualfixedprice+tempprice;
//      }
//  }
//  fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;
//}

//2018-01-23 추가 ktb - usd >-0.4 and ktb - eur spread >04 if then Max[ 1.9%, Min(7%, (2 X KTB 10Y) ? 2.18%) ] X P/M
void BasisDualSpreadFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> _rahigh_bdry2, vector<bool> _rahigh_bdry2in_flag, vector<double> _ralow_bdry2, vector<bool> _ralow_bdry2in_flag, vector<double> _caps, vector<double> _multis, vector<double> _sprds, vector<double> _floors, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    double index1;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        index1 = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex3_mult[i] * fixedrate3[i][j], _rahigh_bdry2[i], _ralow_bdry2[i], _rahigh_bdry2in_flag[i], _ralow_bdry2in_flag[i])) fixedrateincount = fixedrateincount + 1.0;
                index1 = fixedrate1[i][j];
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * (_couponleg_couponrate[i] + min(_caps[i], max(_multis[i] * index1 + _sprds[i], _floors[i]))) * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void SpreadDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate21, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate31, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex21_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex31_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponlegindex21_mult[i] * fixedrate21[i][j] + _couponleg_spread[i], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i] * fixedrate3[i][j] + _couponlegindex31_mult[i] * fixedrate31[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

//generated on 2017-09-27
void SpreadDualFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponleg_spread[i], _caprates[i], _floorrates[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex3_mult[i] * fixedrate3[i][j] + _couponlegindex4_mult[i] * fixedrate4[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / (num_fixing[i]);
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void CDFXDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void StartEndVolFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

/* made on 2019-04-19
 * abs(start-end) 평균 레잇 넣음
*/
void StartEndVolFixedPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<bool> arr_flag, double& fixedlegaccrualfixedprice)
{
    int i, j, k = 0;
    fixedlegaccrualfixedprice = 0.0;

    double fixedrate = 0.0, tempprice = 0.0, arr_value = 0.0, adv_value = 0.0, n = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrate = 0.0;
        n = 0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (arr_flag[k]) {
                    n = n + 1;
                    adv_value = arr_value;
                    arr_value = fixedrate1[i][j];
                    if (n > 1) {
                        fixedrate = fixedrate + _couponlegindex1_mult[i] * abs(arr_value - adv_value);
                    }
                }
                k++;
                //rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
                if (j == num_fixed[i] - 1) fixedrate = max(min(_rahigh_bdry[i], fixedrate / (n - 1)), _ralow_bdry[i]);
            }

            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * fixedrate;

            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void QuantoDualFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            tempprice = exp(-zc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * _couponleg_couponrate[i] * couponlegtau[i] * fixedrateincount / num_fixing[i];
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void QuantoDualFixedAccrualZeroCouponPrice(CDate today, double ondf, double settlement_date_df, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<double> couponlegtau, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double fixedrateincount = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        fixedrateincount = 0.0;
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex1_mult[i] * fixedrate1[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i]) && rangeincheck(_couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount = fixedrateincount + 1.0;
            }
            //tempprice=exp(-zc(_couponleg_paydate[i])*couponlegT[i])*_couponleg_notional[i]*_couponleg_couponrate[i]*couponlegtau[i]*fixedrateincount/num_fixing[i];
            tempprice = _couponleg_notional[i] * couponlegtau[i] * (_couponleg_couponrate[i] * fixedrateincount / num_fixing[i] + _couponleg_spread[i]);
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;
}

void GeneralSwapCouponlegFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation indexzc, int num_remained_couponleg_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegfixedprice)
{
    int i;
    fixedlegfixedprice = 0.0;
    double fixed_coupon_fixedprice = 0.0, tempprice = 0.0;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        tempprice = exp(-indexzc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (_couponleg_couponrate[i] + _couponleg_spread[i]);
        fixed_coupon_fixedprice = fixed_coupon_fixedprice + tempprice;
    }
    fixedlegfixedprice = fixed_coupon_fixedprice / ondf / settlement_date_df;
}

void GeneralSwapCouponlegNonFixedPrice(CDate today, double ondf, double settlement_date_df, CInterpolation indexzc, int num_remained_couponleg_cf, vector<int> num_index_cf, vector<double> _couponleg_notional, vector<double> _couponleg_couponrate, vector<double> _couponleg_spread, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<vector<double>> index_tau, vector<CDate> index_startdate, vector<double> index_startT, vector<vector<CDate>> index_paydate, vector<vector<double>> index_payT, vector<double> couponindex_fixingT, vector<double> index_T, vector<double> index_vol, vector<vector<double>>  index_d, vector<double> index_t, double& fixedlegnonfixedprice)
{
    int i, j;
    fixedlegnonfixedprice = 0.0;
    double level, tmpprice = 0.0;
    vector<double> fixedlegnonfixedcf(num_remained_couponleg_cf, 0.0);
    vector<double> index_fra(num_remained_couponleg_cf), index_adjrate(num_remained_couponleg_cf);

    //ofstream foutr("r.txt");

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        level = 0.0;
        for (j = 0; j < num_index_cf[i]; j++) level = level + index_tau[i][j] * exp(-indexzc(index_paydate[i][j]) * index_payT[i][j]);
        index_fra[i] = (exp(-indexzc(index_startdate[i]) * index_startT[i]) - exp(-indexzc(index_paydate[i][num_index_cf[i] - 1]) * index_T[i])) / level;
        index_adjrate[i] = index_fra[i] + convadj(index_fra[i], index_vol[i], couponindex_fixingT[i], num_index_cf[i], index_d[i], index_t[i]);

        fixedlegnonfixedcf[i] = exp(-indexzc(_couponleg_paydate[i]) * couponlegT[i]) * _couponleg_notional[i] * couponlegtau[i] * (index_adjrate[i] + _couponleg_spread[i]);

        fixedlegnonfixedprice = fixedlegnonfixedprice + fixedlegnonfixedcf[i];

        //foutr<<i<<"\t"<<index_fra[i]*100.0<<"\t"<<index_adjrate[i]*100.0<<endl;
    }
    fixedlegnonfixedprice = fixedlegnonfixedprice / ondf / settlement_date_df;
}


void ComboLeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate1count, vector<double>& fixedrate2count, vector<double>& fixedrate3count, vector<double>& fixedrate4count)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                fixedrate1count[i] = fixedrate1count[i] + _couponlegindex1_mult[i] * fixedrate1[i][j];
                fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i] * fixedrate2[i][j];
                fixedrate3count[i] = fixedrate3count[i] + _couponlegindex3_mult[i] * fixedrate3[i][j];
                fixedrate4count[i] = fixedrate4count[i] + _couponlegindex4_mult[i] * fixedrate4[i][j];
            }
        }
    }
}

void ComboLeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate3, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex3_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate1count, vector<double>& fixedrate2count, vector<double>& fixedrate3count, vector<double>& fixedrate4count, vector<double> _rahigh_bdry, vector<double> _ralow_bdry, vector<bool> _rahigh_bdryin_flag, vector<bool> _ralow_bdryin_flag, int& tmpfixedincount)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(_couponlegindex4_mult[i] * fixedrate4[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) tmpfixedincount = tmpfixedincount + 1;
                fixedrate1count[i] = fixedrate1count[i] + _couponlegindex1_mult[i] * fixedrate1[i][j];
                fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i] * fixedrate2[i][j];
                fixedrate3count[i] = fixedrate3count[i] + _couponlegindex3_mult[i] * fixedrate3[i][j];
                fixedrate4count[i] = fixedrate4count[i] + _couponlegindex4_mult[i] * fixedrate4[i][j];
            }
        }
    }
}


void LeveragedAverageFixedPrice(int num_remained_couponleg_cf, vector<int> num_fixed, vector<vector<double>> fixedrate2, vector<vector<double>> fixedrate4, vector<bool> _ra_flag, vector<double> _couponlegindex2_mult, vector<double> _couponlegindex4_mult, vector<double>& fixedrate2count, vector<double>& fixedrate4count)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                fixedrate2count[i] = fixedrate2count[i] + _couponlegindex2_mult[i] * fixedrate2[i][j];
                fixedrate4count[i] = fixedrate4count[i] + _couponlegindex4_mult[i] * fixedrate4[i][j];
            }
        }
    }
}

void QuantoLeveragedAverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                //if(rangeincheck(fixedrate1[i][j],_rahigh_bdry[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i])) fixedrateincount=fixedrateincount+min(max(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponlegindex1_mult[i]*_couponleg_spread[i],_floorrates[i]),_caprates[i]);
                fixedrateincount[i] = fixedrateincount[i] + _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponlegindex1_mult[i] * _couponleg_spread[i];   //

            }
        }
    }
}

void QuantoLeveragedAverageFixedAccrualPrice1(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, double& fixedlegaccrualfixedprice)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                //if(rangeincheck(fixedrate1[i][j],_rahigh_bdry[i],_ralow_bdry[i],_rahigh_bdryin_flag[i],_ralow_bdryin_flag[i])) fixedrateincount=fixedrateincount+min(max(_couponlegindex1_mult[i]*fixedrate1[i][j]+_couponlegindex2_mult[i]*fixedrate2[i][j]+_couponlegindex1_mult[i]*_couponleg_spread[i],_floorrates[i]),_caprates[i]);
                fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + min(max((_couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j]) / num_fixing[i] + _couponleg_spread[i], _floorrates[i]), _caprates[i]) * _couponleg_notional[i];

            }
        }
    }
}

void QuantoLeveragedAverageFixedAccrualPrice(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_remained_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<bool> _ra_flag, vector<double> _couponleg_notional, vector<double> _couponleg_spread, vector<double> _couponlegindex1_mult, vector<double> _couponlegindex2_mult, vector<double> _caprates, vector<double> _floorrates, vector<double> _rahigh_bdry, vector<bool> _rahigh_bdryin_flag, vector<double> _ralow_bdry, vector<bool> _ralow_bdryin_flag, vector<CDate> _couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg)
{
    int i, j;
    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(rngchkfixedrate[i][j], _rahigh_bdry[i], _ralow_bdry[i], _rahigh_bdryin_flag[i], _ralow_bdryin_flag[i])) fixedrateincount[i] = fixedrateincount[i] + 1;
                fixedavg[i] = fixedavg[i] + _couponlegindex1_mult[i] * fixedrate1[i][j] + _couponlegindex2_mult[i] * fixedrate2[i][j] + _couponlegindex1_mult[i] * _couponleg_spread[i];

            }
        }
    }
}

//2018-07-16 fixed zero price 생성
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
)
{
    int i, j;
    double fixedrateincount = 0.0, rng_coup = 0.0;

    for (i = 0; i < num_fixed_couponleg_cf; i++)
    {
        rng_coup = 0.0;
        fixedrateincount = 0.0;

        if (ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            rngchkfixedrate[i][j]
                            , rahigh_bdry[i]
                            , ralow_bdry[i]
                            , rahigh_bdryin_flag[i]
                            , ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
                rng_coup
                    = rng_coup
                    + (
                        couponlegindex1_mult[i]
                        * fixedrate1[i][j]
                        + couponlegindex2_mult[i]
                        * fixedrate2[i][j]
                        );

            }
        }

        fixedlegfixedprice
            = fixedlegfixedprice
            + (
                couponleg_couponrate[i]
                + min
                (
                    max
                    (
                        rng_coup / num_fixed[i] + couponleg_spread[i]
                        , floorrates[i]
                    )
                    , caprates[i]
                )
                * fixedrateincount
                / num_fixed[i]
                )
            * couponleg_notional[i];
    }
}

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
)
{
    int
        i
        , j;

    double
        fixedrateincount = 0.0
        , rng_coup = 0.0;

    for (i = 0; i < num_fixed_couponleg_cf; i++)
    {
        rng_coup = 0.0;

        fixedrateincount = 0.0;

        if (ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            rngchkfixedrate[i][j]
                            , rahigh_bdry[i]
                            , ralow_bdry[i]
                            , rahigh_bdryin_flag[i]
                            , ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount = fixedrateincount + 1.0;
                }
                rng_coup
                    = rng_coup
                    + (
                        couponlegindex1_mult[i]
                        * fixedrate1[i][j]
                        + couponlegindex2_mult[i]
                        * fixedrate2[i][j]
                        );

            }
        }

        fixedlegfixedprice
            = fixedlegfixedprice
            + (
                couponleg_couponrate[i]
                + min
                (
                    max
                    (
                        rng_coup / num_fixed[i] + couponleg_spread[i]
                        , floorrates[i]
                    )
                    , caprates[i]
                )
                * fixedrateincount
                / num_fixed[i]
                )
            * couponleg_notional[i];

        cout << fixedlegfixedprice << endl;
    }
}

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

)
{
    int i, j;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            _rngchkindex_mult[i] * rngchkfixedrate[i][j]
                            + _rngchkindex2_mult[i] * rngchkfixedrate2[i][j]
                            , _rahigh_bdry[i]
                            , _ralow_bdry[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount[i] = fixedrateincount[i] + 1;
                }

                fixedavg[i]
                    = fixedavg[i]
                    + _couponlegindex1_mult[i] * fixedrate1[i][j]
                    + _couponlegindex2_mult[i] * fixedrate2[i][j]
                    + _couponlegindex1_mult[i] * _couponleg_spread[i];

            }
        }
    }
}

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

)
{
    int i, j;

    for (i = 0; i < num_remained_couponleg_cf; i++)
    {
        if (_ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            _rngchkindex_mult[i] * rngchkfixedrate[i][j]
                            + _rngchkindex2_mult[i] * rngchkfixedrate2[i][j]
                            , _rahigh_bdry[i]
                            , _ralow_bdry[i]
                            , _rahigh_bdryin_flag[i]
                            , _ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount[i] = fixedrateincount[i] + 1;
                }

                fixedavg[i]
                    = fixedavg[i]
                    + _couponlegindex1_mult[i] * fixedrate1[i][j]
                    + _couponlegindex2_mult[i] * fixedrate2[i][j]
                    + _couponlegindex1_mult[i] * _couponleg_spread[i];

            }
        }
    }
}


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

)
{
    int
        i
        , j;

    fixedlegaccrualfixedprice = 0.0;

    double tempprice = 0.0;

    for (i = 0; i < num_past_couponleg_cf; i++)
    {
        fixedrateincount[i] = 0.0;

        if (ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            rngchkfixedrate[i][j]
                            , rahigh_bdry[i]
                            , ralow_bdry[i]
                            , rahigh_bdryin_flag[i]
                            , ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount[i] = fixedrateincount[i] + 1;
                }

                fixedavg[i]
                    = fixedavg[i]
                    + couponlegindex1_mult[i] * fixedrate1[i][j]
                    + couponlegindex2_mult[i] * fixedrate2[i][j]
                    + couponleg_spread[i];
            }

            fixedavg[i]
                = min
                (
                    max
                    (
                        fixedavg[i] + couponleg_spread[i]
                        , floorrates[i]
                    )
                    , caprates[i]
                )
                / (num_fixing[i]);

            tempprice
                = couponleg_notional[i]
                * fixedavg[i]
                * couponlegtau[i]
                * fixedrateincount[i]
                / (num_fixing[i]);

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }

}

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

)
{
    int
        i
        , j;

    fixedlegaccrualfixedprice = 0.0;

    double tempprice = 0.0;

    for (i = 0; i < num_past_couponleg_cf; i++)
    {
        fixedrateincount[i] = 0.0;

        if (ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if
                    (
                        rangeincheck
                        (
                            rngchkfixedrate[i][j]
                            , rahigh_bdry[i]
                            , ralow_bdry[i]
                            , rahigh_bdryin_flag[i]
                            , ralow_bdryin_flag[i]
                        )
                        )
                {
                    fixedrateincount[i] = fixedrateincount[i] + 1;
                }

                fixedavg[i]
                    = fixedavg[i]
                    + couponlegindex1_mult[i] * fixedrate1[i][j]
                    + couponlegindex2_mult[i] * fixedrate2[i][j];
            }

            fixedavg[i] = fixedavg[i] / num_fixing[i];

            fixedavg[i]
                = min
                (
                    max
                    (
                        fixedavg[i]
                        + couponleg_spread[i]
                        , floorrates[i]
                    )
                    , caprates[i]
                )
                * fixedrateincount[i]
                / num_fixing[i];

            tempprice
                = couponleg_notional[i]
                * fixedavg[i]
                * couponlegtau[i];

            fixedlegaccrualfixedprice
                = fixedlegaccrualfixedprice
                + tempprice;
        }
    }

}


/* 20170216 */
void QuantoAvgSpreadLiborRAFixedZeroCouponPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_past_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<vector<double>> rngchkfixedrate2,
    vector<double> rngchkindex_mult, vector<double> rngchkindex2_mult, vector<bool> ra_flag, vector<double> couponleg_notional, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> caprates, vector<double> floorrates, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg, double& fixedlegaccrualfixedprice)
{
    int i, j;
    fixedlegaccrualfixedprice = 0.0;
    double tempprice = 0.0;
    for (i = 0; i < num_past_couponleg_cf; i++)
    {
        fixedrateincount[i] = 0.0;
        if (ra_flag[i])
        {
            for (j = 0; j < num_fixed[i]; j++)
            {
                if (rangeincheck(rngchkindex_mult[i] * rngchkfixedrate[i][j] + rngchkindex2_mult[i] * rngchkfixedrate2[i][j], rahigh_bdry[i], ralow_bdry[i], rahigh_bdryin_flag[i], ralow_bdryin_flag[i])) fixedrateincount[i] = fixedrateincount[i] + 1;
                fixedavg[i] = fixedavg[i] + couponlegindex1_mult[i] * fixedrate1[i][j] + couponlegindex2_mult[i] * fixedrate2[i][j] + couponleg_spread[i];
            }
            fixedavg[i] = min(max(fixedavg[i] + couponleg_spread[i], floorrates[i]), caprates[i]) / (num_fixing[i]);
            tempprice = couponleg_notional[i] * fixedavg[i] * couponlegtau[i] * fixedrateincount[i] / (num_fixing[i]);
            fixedlegaccrualfixedprice = fixedlegaccrualfixedprice + tempprice;
        }
    }
    //fixedlegaccrualfixedprice=fixedlegaccrualfixedprice/ondf/settlement_date_df;

}

/* 20190529 생성 */
void QuantoAvgSpreadLiborRAFixedPrice2(CDate today, double ondf, double settlement_date_df, CInterpolation zc, int num_past_couponleg_cf, vector<int> num_fixing, vector<int> num_fixed, vector<double> couponleg_couponrate, vector<vector<double>> fixedrate1, vector<vector<double>> fixedrate2, vector<vector<double>> rngchkfixedrate, vector<vector<double>> rngchkfixedrate2,
    vector<double> rngchkindex_mult, vector<double> rngchkindex2_mult, vector<bool> ra_flag, vector<double> couponleg_notional, vector<double> couponleg_spread, vector<double> couponlegindex1_mult, vector<double> couponlegindex2_mult, vector<double> caprates, vector<double> floorrates, vector<double> rahigh_bdry, vector<bool> rahigh_bdryin_flag, vector<double> ralow_bdry, vector<bool> ralow_bdryin_flag, vector<CDate> couponleg_paydate, vector<double> couponlegtau, vector<double> couponlegT, vector<double>& fixedrateincount, vector<double>& fixedavg, double& fixedlegaccrualfixedprice)
{
    int i = num_past_couponleg_cf, j;
    fixedlegaccrualfixedprice = 0.0;
    double tempprice = 0.0;
    //  for (i = num_past_couponleg_cf-1; i<num_past_couponleg_cf; i++)
    //  {
    fixedrateincount[i] = 0.0;
    if (ra_flag[i])
    {
        for (j = 0; j < num_fixed[i]; j++)
        {
            if (rangeincheck(rngchkindex_mult[i] * rngchkfixedrate[i][j] + rngchkindex2_mult[i] * rngchkfixedrate2[i][j], rahigh_bdry[i], ralow_bdry[i], rahigh_bdryin_flag[i], ralow_bdryin_flag[i])) fixedrateincount[i] = fixedrateincount[i] + 1;
            fixedavg[i] = fixedavg[i] + couponlegindex1_mult[i] * fixedrate1[i][j] + couponlegindex2_mult[i] * fixedrate2[i][j] + couponleg_spread[i];
        }
        fixedavg[i] = min(max(fixedavg[i] / num_fixed[i] + couponleg_spread[i], floorrates[i]), caprates[i]);
        tempprice = couponleg_notional[i] * couponlegtau[i] * fixedavg[i] * fixedrateincount[i] / num_fixed[i];
    }
    else {
        tempprice = couponleg_notional[i] * couponlegtau[i] * couponleg_couponrate[i];
    }
    //  }
    fixedlegaccrualfixedprice = tempprice * exp(-zc(couponleg_paydate[i]) * couponlegT[i]);
    fixedlegaccrualfixedprice = fixedlegaccrualfixedprice / ondf / settlement_date_df;

}

void PlainIROption_HW1F(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y, vector<double>& dydparam)
{
    if (flag == "CAP")
    {
        int i, ncaplet = int(tau.size());
        double BtT, sigp, h, dsigpdsig, dsigpda, dhdsig, dhda, exp2at, expatT, lnPtTX, sqrtexp2at;
        y = 0.0;
        dydparam[0] = 0.0;
        dydparam[1] = 0.0;
        double a = min(max(pow(param[0], 2.0), Error), 1.0 / Error), sig = min(max(pow(param[1], 2.0), Error), 1.0 / Error);
        for (i = 1; i < ncaplet; i++)
        {
            expatT = exp(-a * (t[i + 1] - t[i]));
            exp2at = exp(-2.0 * a * t[i]);
            lnPtTX = log(df[i + 1] * (1.0 + x * tau[i]) / df[i]);
            sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
            BtT = (1.0 - expatT) / a;
            dsigpdsig = 2.0 * param[1] * sqrtexp2at * BtT;
            sigp = sig * sqrtexp2at * BtT;
            h = lnPtTX / sigp + 0.5 * sigp;
            dhdsig = -lnPtTX * dsigpdsig / (sigp * sigp) + 0.5 * dsigpdsig;
            dsigpda = 0.5 * param[0] * ((2.0 * a * t[i] + 1.0) * exp2at - 1.0) * BtT / (sqrtexp2at * a * a) + param[0] * sig * sqrtexp2at * ((a * (t[i + 1] - t[i]) + 1.0) * expatT - 1.0) / (a * a);
            dhda = -lnPtTX * dsigpda / (sigp * sigp) + 0.5 * dsigpda;
            y = y + df[i] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * (1.0 + x * tau[i]) * cumnormal_Pol_App(-h);
            dydparam[0] = dydparam[0] + df[i] * max(pdfnormal(-h + sigp), Error) * (-dhda + dsigpda) - df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(-h), Error) * (-dhda);
            dydparam[1] = dydparam[1] + df[i] * max(pdfnormal(-h + sigp), Error) * (-dhdsig + dsigpdsig) - df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(-h), Error) * (-dhdsig);
        }
    }
    else
    {
        int i, nswaplet = int(tau.size());
        double sigp, h, dsigpdsig, dsigpda, dhdsig, dhda, exp2at, expatT, lnPtTX, sqrtexp2at, rstar, X, denomrstar = 0.0, nomrstara = 0.0, nomrstarsig = 0.0, dXda, dXdsig, drstarda, drstardsig;
        y = 0.0;
        dydparam[0] = 0.0;
        dydparam[1] = 0.0;
        double a = min(max(pow(param[0],2.0),Error),1.0 / Error), sig = min(max(pow(param[1],2.0),Error),1.0 / Error);
        vector<double> AtT(nswaplet), BtT(nswaplet), c(nswaplet), dAda(nswaplet), dAdsig(nswaplet), dBda(nswaplet), expBrstar(nswaplet);
        exp2at = exp(-2.0 * a * t[0]);
        sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
        for (i = 0; i < nswaplet; i++)
        {
            c[i] = x * tau[i];
            expatT = exp(-a * (t[i + 1] - t[0]));
            BtT[i] = (1.0 - expatT) / a;
            AtT[i] = df[i + 1] / df[0] * max(exp(BtT[i] * fM - 0.25 * pow(sig * BtT[i],2.0) * (1.0 - exp2at) / a),Error);
            dBda[i] = 2.0 * param[0] * (-BtT[i] + (t[i + 1] - t[0]) * expatT) / a;
            dAda[i] = AtT[i] * (dBda[i] * fM + 2.0 * param[0] * (pow(0.5 * sig * BtT[i] / a,2.0) * (1.0 - exp2at) - 0.5 * t[0] * exp2at * pow(sig * BtT[i],2.0) / a) - 0.5 * pow(sig,2.0) * (1.0 - exp2at) * BtT[i] / a * dBda[i]);
            dAdsig[i] = -param[1] * AtT[i] * sig * (1.0 - exp2at) * pow(BtT[i],2.0) / a;
        }
        c[nswaplet - 1] = c[nswaplet - 1] + 1.0;
        rstar_HW1F(c,AtT,BtT,rstar);
        for (i = 0; i < nswaplet; i++)
        {
            expBrstar[i] = exp(-BtT[i] * rstar);
            nomrstara = nomrstara + c[i] * (dAda[i] - AtT[i] * rstar * dBda[i]) * expBrstar[i];
            nomrstarsig = nomrstarsig + c[i] * dAdsig[i] * expBrstar[i];
            denomrstar = denomrstar + c[i] * AtT[i] * BtT[i] * expBrstar[i];
        }
        drstarda = nomrstara / denomrstar;
        drstardsig = nomrstarsig / denomrstar;
        for (i = 0; i < nswaplet; i++)
        {
            X = AtT[i] * expBrstar[i];
            dsigpdsig = 2.0 * param[1] * sqrtexp2at * BtT[i];
            sigp = sig * sqrtexp2at * BtT[i];
            lnPtTX = log(df[i + 1] / (df[0] * X));
            h = lnPtTX / sigp + 0.5 * sigp;
            dXda = (dAda[i] - AtT[i] * dBda[i] * rstar - AtT[i] * BtT[i] * drstarda) * expBrstar[i];
            dXdsig = (dAdsig[i] - AtT[i] * BtT[i] * drstardsig) * expBrstar[i];
            dhdsig = -lnPtTX * dsigpdsig / (sigp * sigp) - dXdsig / X / sigp + 0.5 * dsigpdsig;
            dsigpda = 0.5 * param[0] * ((2.0 * a * t[0] + 1.0) * exp2at - 1.0) * BtT[i] / (sqrtexp2at * a * a) + sig * sqrtexp2at * dBda[i];
            dhda = -lnPtTX * dsigpda / (sigp * sigp) - dXda / (sigp * X) + 0.5 * dsigpda;
            y = y + c[i] * (X * df[0] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * cumnormal_Pol_App(-h));
            dydparam[0] = dydparam[0] + c[i] * (df[0] * (dXda * cumnormal_Pol_App(-h + sigp) + X * max(pdfnormal(-h + sigp),Error) * (-dhda + dsigpda)) + df[i + 1] * max(pdfnormal(-h),Error) * dhda);
            dydparam[1] = dydparam[1] + c[i] * (df[0] * (dXdsig * cumnormal_Pol_App(-h + sigp) + X * max(pdfnormal(-h + sigp),Error) * (-dhdsig + dsigpdsig)) - df[i + 1] * max(pdfnormal(-h),Error) * (-dhdsig));
        }
    }
}

void mrqcof_HW1F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&))
{
    int i, j, k, l, m, mfit = 0, ndata = int(x.size()), ma = int(a.size());
    double ymod, wt, sig2i, dy;
    vector<double> dyda(ma);
    for (j = 0; j < ma; j++) if (ia[j]) mfit++;
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
        beta[j] = 0.0;
    }
    chisq = 0.0;
    for (i = 0; i < ndata; i++)
    {
        funcs(flag[i], df[i], t[i], tau[i], fM[i], x[i], a, ymod, dyda);
        sig2i = 1.0 / (sig[i] * sig[i]);
        dy = y[i] - ymod;
        for (j = 0, l = 0; l < ma; l++)
        {
            if (ia[l])
            {
                wt = dyda[l] * sig2i;
                for (k = 0, m = 0; m < l + 1; m++) if (ia[m]) alpha[j][k++] += wt * dyda[m];
                beta[j++] += dy * wt;
            }
        }
        chisq += dy * dy * sig2i;
    }
    for (j = 1; j < mfit; j++) for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];

}

void mrqmin_HW1F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&), double& alamda, bool& singularflag)
{
    static int mfit;
    static double ochisq;
    int j, k, l, ma = int(a.size());
    static vector<vector<double>> oneda(ma);
    for (j = 0; j < ma; j++) oneda[j] = vector<double>(1);
    static vector<double> atry(ma), beta(ma), da(ma);
    if (alamda < 0.0)
    {
        mfit = 0;
        for (j = 0; j < ma; j++) if (ia[j]) mfit++;
        alamda = 0.001;
        mrqcof_HW1F(flag, df, t, tau, fM, x, y, sig, a, ia, alpha, beta, chisq, funcs);
        ochisq = chisq;
        for (j = 0; j < ma; j++) atry[j] = a[j];
    }
    vector<vector<double>> temp(mfit);
    for (j = 0; j < ma; j++) temp[j] = vector<double>(mfit);
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
        covar[j][j] = alpha[j][j] * (1.0 + alamda);
        for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
        oneda[j][0] = beta[j];
    }
    gaussj(temp, oneda, singularflag);
    if (singularflag != true)
    {
        for (j = 0; j < mfit; j++)
        {
            for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k];
            da[j] = oneda[j][0];
        }
        if (alamda == 0.0)
        {
            covsrt(covar, ia, mfit);
            covsrt(alpha, ia, mfit);
            return;
        }
        for (j = 0, l = 0; l < ma; l++) if (ia[l]) atry[l] = a[l] + da[j++];
        mrqcof_HW1F(flag, df, t, tau, fM, x, y, sig, atry, ia, alpha, beta, chisq, funcs);
        if (chisq < ochisq)
        {
            alamda *= 0.1;
            ochisq = chisq;
            for (j = 0; j < mfit; j++)
            {
                for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
                beta[j] = da[j];
            }
            for (l = 0; l < ma; l++) a[l] = atry[l];
        }
        else
        {
            alamda *= 10.0;
            chisq = ochisq;
        }
    }
}

void rstar_HW1F(vector<double>& c, vector<double>& AtT, vector<double>& BtT, double& rstar)
{
    double temperror = 1.0e-12;
    int n = int(c.size()), i;
    double left0 = -100.0, left, right0 = 100.0, right, ans, sum = 0.0, ans0, shift = 2.0;
    double sumr = 0.0, suml = 0.0;
    for (i = 0; i < n; i++) sumr = sumr + c[i] * AtT[i] * exp(-BtT[i] * right0);
    for (i = 0; i < n; i++) suml = suml + c[i] * AtT[i] * exp(-BtT[i] * left0);
    
    while (sumr > 1.0)
    {
        right0 = right0 * shift;
        sumr = 0.0;
        for (i = 0; i < n; i++) sumr = sumr + c[i] * AtT[i] * exp(-BtT[i] * right0);
    }
    while(suml < 1.0)
    {
        left0 = left0 * shift;
        suml = 0.0;
        for (i = 0; i < n; i++) suml = suml + c[i] * AtT[i] * exp(-BtT[i] * left0);
    }

    left = left0;
    right = right0;
    rstar = 0.5 * (left0 + right0);
    for (i = 0; i < n; i++) sum = sum + c[i] * AtT[i] * exp(-BtT[i] * rstar);
    ans = sum - 1.0;
    ans0 = ans - 1.0;
    while (fabs(ans) > Error&& rstar > left0 + Error && rstar < right0 - Error && ans0 != ans)
    {
        ans0 = ans;
        if (ans > 0.0)
        {
            left = rstar;
            rstar = 0.5 * (left + right);
            sum = 0.0;
            for (i = 0; i < n; i++) sum = sum + c[i] * AtT[i] * exp(-BtT[i] * rstar);
            ans = sum - 1.0;
        }
        else if (ans < 0.0)
        {
            right = rstar;
            rstar = 0.5 * (left + right);
            sum = 0.0;
            for (i = 0; i < n; i++) sum = sum + c[i] * AtT[i] * exp(-BtT[i] * rstar);
            ans = sum - 1.0;
        }
        else break;
    }
}

void gaussj(vector<vector<double>>& a, vector<vector<double>>& b, bool& singularflag)
{
    int i, icol, irow, j, k, l, ll, n = int(a.size()), m = 0;
    if (int(b.size()) > 0) m = int(b[0].size());
    double big, dum, pivinv;
    vector<int> indxc(n), indxr(n), ipiv(n);
    for (j = 0; j < n; j++) ipiv[j] = 0;
    for (i = 0; i < n; i++)
    {
        big = 0.0;
        for (j = 0; j < n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (k = 0; k < n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j][k]) >= big)
                        {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l = 0; l < n; l++) swap(a[irow][l], a[icol][l]);
            for (l = 0; l < m; l++) swap(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0)
        {
            cout << "gaussj:Singular Matrix" << endl;
            singularflag = true;
        }
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 0; l < n; l++) a[icol][l] *= pivinv;
        for (l = 0; l < m; l++) b[icol][l] *= pivinv;
        for (ll = 0; ll < n; ll++)
        {
            if (ll != icol)
            {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
                for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
            }
        }
    }
    for (l = n - 1; l >= 0; l--) if (indxr[l] != indxc[l]) for (k = 0; k < n; k++) swap(a[k][indxr[l]], a[k][indxc[l]]);
}


void mrqmin_HW1Fsig(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param1, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs_sig(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&), double& alamda, bool& singularflag)
{
    static int mfit;
    static double ochisq;
    int j, k, l, ma = 1;
    static vector<vector<double>> oneda(ma);
    for (j = 0; j < ma; j++) oneda[j] = vector<double>(1);
    static vector<double>  beta(ma), da(ma);
    static double atry;
    if (alamda < 0.0)
    {
        mfit = 0;
        for (j = 0; j < ma; j++) if (ia[j]) mfit++;
        alamda = 0.001;
        mrqcof_HW1Fsig(flag, df, t, tau, fM, x, y, sig, a, param1, ia, alpha, beta, chisq, funcs_sig);
        ochisq = chisq;
        atry = a;
    }
    vector<vector<double>> temp(mfit);
    for (j = 0; j < ma; j++) temp[j] = vector<double>(mfit);
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
        covar[j][j] = alpha[j][j] * (1.0 + alamda);
        for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
        oneda[j][0] = beta[j];
    }
    gaussj(temp, oneda, singularflag);
    if (singularflag != true)
    {
        for (j = 0; j < mfit; j++)
        {
            for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k];
            da[j] = oneda[j][0];
        }
        if (alamda == 0.0)
        {
            covsrt(covar, ia, mfit);
            covsrt(alpha, ia, mfit);
            return;
        }
        for (j = 0, l = 0; l < ma; l++) if (ia[l]) atry = a + da[j++];
        mrqcof_HW1Fsig(flag, df, t, tau, fM, x, y, sig, atry, param1, ia, alpha, beta, chisq, funcs_sig);
        if (chisq < ochisq)
        {
            alamda *= 0.1;
            ochisq = chisq;
            for (j = 0; j < mfit; j++)
            {
                for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
                beta[j] = da[j];
            }
            a = atry;
        }
        else
        {
            alamda *= 10.0;
            chisq = ochisq;
        }
    }
}

void mrqcof_HW1Fsig(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param1, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs_sig(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&))
{
    int i, j, k, l, m, mfit = 0, ndata = int(x.size()), ma = 1;
    double ymod, wt, sig2i, dy;
    double dyda;
    for (j = 0; j < ma; j++) if (ia[j]) mfit++;
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
        beta[j] = 0.0;
    }
    chisq = 0.0;
    for (i = 0; i < ndata; i++)
    {
        funcs_sig(flag[i], df[i], t[i], tau[i], fM[i], x[i], param1, a, ymod, dyda);
        sig2i = 1.0 / (sig[i] * sig[i]);
        dy = y[i] - ymod;
        for (j = 0, l = 0; l < ma; l++)
        {
            if (ia[l])
            {
                wt = dyda * sig2i;
                for (k = 0, m = 0; m < l + 1; m++) if (ia[m]) alpha[j][k++] += wt * dyda;
                beta[j++] += dy * wt;
            }
        }
        chisq += dy * dy * sig2i;
    }
    for (j = 1; j < mfit; j++) for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];

}


void PlainIROption_HW1Fsig(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, double param1, double& param0, double& y, double& dydparam0)
{
    if (flag == "CAP")
    {
        int i, ncaplet = int(tau.size());
        double BtT, sigp, h, dsigpda, dhda, exp2at, expatT, lnPtTX, sqrtexp2at;
        y = 0.0;
        dydparam0 = 0.0;
        double a = min(max(pow(param0, 2.0), Error), 1.0 / Error), sig = min(max(pow(param1, 2.0), Error), 1.0 / Error);
        for (i = 1; i < ncaplet; i++)
        {
            expatT = exp(-a * (t[i + 1] - t[i]));
            exp2at = exp(-2.0 * a * t[i]);
            lnPtTX = log(df[i + 1] * (1.0 + x * tau[i]) / df[i]);
            sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
            BtT = (1.0 - expatT) / a;
            sigp = sig * sqrtexp2at * BtT;
            h = lnPtTX / sigp + 0.5 * sigp;
            dsigpda = 0.5 * param0 * ((2.0 * a * t[i] + 1.0) * exp2at - 1.0) * BtT / (sqrtexp2at * a * a) + param0 * sig * sqrtexp2at * ((a * (t[i + 1] - t[i]) + 1.0) * expatT - 1.0) / (a * a);
            dhda = -lnPtTX * dsigpda / (sigp * sigp) + 0.5 * dsigpda;
            y = y + df[i] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * (1.0 + x * tau[i]) * cumnormal_Pol_App(-h);
            dydparam0 = dydparam0 + df[i] * max(pdfnormal(-h + sigp), Error) * (-dhda + dsigpda) - df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(-h), Error) * (-dhda);
        }
    }
    else
    {
        int i, nswaplet = int(tau.size());
        double sigp, h, dsigpda, dhda, exp2at, expatT, lnPtTX, sqrtexp2at, rstar, X, denomrstar = 0.0, nomrstara = 0.0, dXda, drstarda;
        y = 0.0;
        dydparam0 = 0.0;
        double a = min(max(pow(param0,2.0),Error),1.0 / Error), sig = min(max(pow(param1,2.0),Error),1.0 / Error);
        vector<double> AtT(nswaplet), BtT(nswaplet), c(nswaplet), dAda(nswaplet), dBda(nswaplet), expBrstar(nswaplet);
        exp2at = exp(-2.0 * a * t[0]);
        sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
        for (i = 0; i < nswaplet; i++)
        {
            c[i] = x * tau[i];
            expatT = exp(-a * (t[i + 1] - t[0]));
            BtT[i] = (1.0 - expatT) / a;
            AtT[i] = df[i + 1] / df[0] * max(exp(BtT[i] * fM - 0.25 * pow(sig * BtT[i],2.0) * (1.0 - exp2at) / a),Error);
            dBda[i] = 2.0 * param0 * (-BtT[i] + (t[i + 1] - t[0]) * expatT) / a;
            dAda[i] = AtT[i] * (dBda[i] * fM + 2.0 * param0 * (pow(0.5 * sig * BtT[i] / a,2.0) * (1.0 - exp2at) - 0.5 * t[0] * exp2at * pow(sig * BtT[i],2.0) / a) - 0.5 * pow(sig,2.0) * (1.0 - exp2at) * BtT[i] / a * dBda[i]);
        }
        c[nswaplet - 1] = c[nswaplet - 1] + 1.0;
        rstar_HW1F(c,AtT,BtT,rstar);
        for (i = 0; i < nswaplet; i++)
        {
            expBrstar[i] = exp(-BtT[i] * rstar);
            nomrstara = nomrstara + c[i] * (dAda[i] - AtT[i] * rstar * dBda[i]) * expBrstar[i];
            denomrstar = denomrstar + c[i] * AtT[i] * BtT[i] * expBrstar[i];
        }
        drstarda = nomrstara / denomrstar;
        for (i = 0; i < nswaplet; i++)
        {
            X = AtT[i] * expBrstar[i];
            sigp = sig * sqrtexp2at * BtT[i];
            lnPtTX = log(df[i + 1] / (df[0] * X));
            h = lnPtTX / sigp + 0.5 * sigp;
            dXda = (dAda[i] - AtT[i] * dBda[i] * rstar - AtT[i] * BtT[i] * drstarda) * expBrstar[i];
            dsigpda = 0.5 * param0 * ((2.0 * a * t[0] + 1.0) * exp2at - 1.0) * BtT[i] / (sqrtexp2at * a * a) + sig * sqrtexp2at * dBda[i];
            dhda = -lnPtTX * dsigpda / (sigp * sigp) - dXda / (sigp * X) + 0.5 * dsigpda;
            y = y + c[i] * (X * df[0] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * cumnormal_Pol_App(-h));
            dydparam0 = dydparam0 + c[i] * (df[0] * (dXda * cumnormal_Pol_App(-h + sigp) + X * max(pdfnormal(-h + sigp),Error) * (-dhda + dsigpda)) + df[i + 1] * max(pdfnormal(-h),Error) * dhda);
        }
    }
}

void mrqmin_HW1Fa(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param0, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs_a(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&), double& alamda, bool& singularflag)
{
    static int mfit;
    static double ochisq;
    int j, k, l, ma = 1;
    static vector<vector<double>> oneda(ma);
    for (j = 0; j < ma; j++) oneda[j] = vector<double>(1);
    static vector<double>  beta(ma), da(ma);
    static double atry;
    if (alamda < 0.0)
    {
        mfit = 0;
        for (j = 0; j < ma; j++) if (ia[j]) mfit++;
        alamda = 0.001;
        mrqcof_HW1Fa(flag, df, t, tau, fM, x, y, sig, a, param0, ia, alpha, beta, chisq, funcs_a);
        ochisq = chisq;
        atry = a;
    }
    vector<vector<double>> temp(mfit);
    for (j = 0; j < ma; j++) temp[j] = vector<double>(mfit);
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
        covar[j][j] = alpha[j][j] * (1.0 + alamda);
        for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
        oneda[j][0] = beta[j];
    }
    gaussj(temp, oneda, singularflag);
    if (singularflag != true)
    {
        for (j = 0; j < mfit; j++)
        {
            for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k];
            da[j] = oneda[j][0];
        }
        if (alamda == 0.0)
        {
            covsrt(covar, ia, mfit);
            covsrt(alpha, ia, mfit);
            return;
        }
        for (j = 0, l = 0; l < ma; l++) if (ia[l]) atry = a + da[j++];
        mrqcof_HW1Fa(flag, df, t, tau, fM, x, y, sig, atry, param0, ia, alpha, beta, chisq, funcs_a);
        if (chisq < ochisq)
        {
            alamda *= 0.1;
            ochisq = chisq;
            for (j = 0; j < mfit; j++)
            {
                for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
                beta[j] = da[j];
            }
            a = atry;
        }
        else
        {
            alamda *= 10.0;
            chisq = ochisq;
        }
    }
}


void mrqcof_HW1Fa(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, double& a, double& param0, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs_a(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, double, double&, double&, double&))
{
    int i, j, k, l, m, mfit = 0, ndata = int(x.size()), ma = 1;
    double ymod, wt, sig2i, dy;
    double dyda;
    for (j = 0; j < ma; j++) if (ia[j]) mfit++;
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
        beta[j] = 0.0;
    }
    chisq = 0.0;
    for (i = 0; i < ndata; i++)
    {
        funcs_a(flag[i], df[i], t[i], tau[i], fM[i], x[i], param0, a, ymod, dyda);
        sig2i = 1.0 / (sig[i] * sig[i]);
        dy = y[i] - ymod;
        for (j = 0, l = 0; l < ma; l++)
        {
            if (ia[l])
            {
                wt = dyda * sig2i;
                for (k = 0, m = 0; m < l + 1; m++) if (ia[m]) alpha[j][k++] += wt * dyda;
                beta[j++] += dy * wt;
            }
        }
        chisq += dy * dy * sig2i;
    }
    for (j = 1; j < mfit; j++) for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];

}


void PlainIROption_HW1Fa(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, double param0, double& param1, double& y, double& dydparam1)
{
    if (flag == "CAP")
    {
        int i, ncaplet = int(tau.size());
        double BtT, sigp, h, dsigpdsig, dhdsig, exp2at, expatT, lnPtTX, sqrtexp2at;
        y = 0.0;
        dydparam1 = 0.0;
        double a = min(max(pow(param0, 2.0), Error), 1.0 / Error), sig = min(max(pow(param1, 2.0), Error), 1.0 / Error);
        for (i = 1; i < ncaplet; i++)
        {
            expatT = exp(-a * (t[i + 1] - t[i]));
            exp2at = exp(-2.0 * a * t[i]);
            lnPtTX = log(df[i + 1] * (1.0 + x * tau[i]) / df[i]);
            sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
            BtT = (1.0 - expatT) / a;
            dsigpdsig = 2.0 * param1 * sqrtexp2at * BtT;
            sigp = sig * sqrtexp2at * BtT;
            h = lnPtTX / sigp + 0.5 * sigp;
            dhdsig = -lnPtTX * dsigpdsig / (sigp * sigp) + 0.5 * dsigpdsig;
            y = y + df[i] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * (1.0 + x * tau[i]) * cumnormal_Pol_App(-h);
            dydparam1 = dydparam1 + df[i] * max(pdfnormal(-h + sigp), Error) * (-dhdsig + dsigpdsig) - df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(-h), Error) * (-dhdsig);
        }
    }
    else
    {
        int i, nswaplet = int(tau.size());
        double sigp, h, dsigpdsig, dhdsig, exp2at, expatT, lnPtTX, sqrtexp2at, rstar, X, denomrstar = 0.0, nomrstara = 0.0, nomrstarsig = 0.0, dXdsig, drstardsig;
        y = 0.0;
        dydparam1 = 0.0;
        double a = min(max(pow(param0,2.0),Error),1.0 / Error), sig = min(max(pow(param1,2.0),Error),1.0 / Error);
        vector<double> AtT(nswaplet), BtT(nswaplet), c(nswaplet), dAdsig(nswaplet), expBrstar(nswaplet);
        exp2at = exp(-2.0 * a * t[0]);
        sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
        for (i = 0; i < nswaplet; i++)
        {
            c[i] = x * tau[i];
            expatT = exp(-a * (t[i + 1] - t[0]));
            BtT[i] = (1.0 - expatT) / a;
            AtT[i] = df[i + 1] / df[0] * max(exp(BtT[i] * fM - 0.25 * pow(sig * BtT[i],2.0) * (1.0 - exp2at) / a),Error);
            dAdsig[i] = -param1 * AtT[i] * sig * (1.0 - exp2at) * pow(BtT[i],2.0) / a;
        }
        c[nswaplet - 1] = c[nswaplet - 1] + 1.0;
        rstar_HW1F(c,AtT,BtT,rstar);
        for (i = 0; i < nswaplet; i++)
        {
            expBrstar[i] = exp(-BtT[i] * rstar);
            nomrstarsig = nomrstarsig + c[i] * dAdsig[i] * expBrstar[i];
            denomrstar = denomrstar + c[i] * AtT[i] * BtT[i] * expBrstar[i];
        }
        drstardsig = nomrstarsig / denomrstar;
        for (i = 0; i < nswaplet; i++)
        {
            X = AtT[i] * expBrstar[i];
            dsigpdsig = 2.0 * param1 * sqrtexp2at * BtT[i];
            sigp = sig * sqrtexp2at * BtT[i];
            lnPtTX = log(df[i + 1] / (df[0] * X));
            h = lnPtTX / sigp + 0.5 * sigp;
            dXdsig = (dAdsig[i] - AtT[i] * BtT[i] * drstardsig) * expBrstar[i];
            dhdsig = -lnPtTX * dsigpdsig / (sigp * sigp) - dXdsig / X / sigp + 0.5 * dsigpdsig;
            y = y + c[i] * (X * df[0] * cumnormal_Pol_App(-h + sigp) - df[i + 1] * cumnormal_Pol_App(-h));
            dydparam1 = dydparam1 + c[i] * (df[0] * (dXdsig * cumnormal_Pol_App(-h + sigp) + X * max(pdfnormal(-h + sigp),Error) * (-dhdsig + dsigpdsig)) - df[i + 1] * max(pdfnormal(-h),Error) * (-dhdsig));
        }
    }
}


void PlainIROption_HW2F
(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y, vector<double>& dydparam)
{
    double err = Error * (1.0e+10);//1.0e-5;
    if (flag == "CAP")
    {
        int i, ncaplet = int(tau.size());
        double expatT, exp2at, expbtT, exp2bt, expapbt, voltT, dvoltTda, dvoltTdb, dvoltTdsig, dvoltTdeta, dvoltTdrho, lnPtTX, pid1, pid2;
        y = 0.0;
        dydparam[0] = 0.0;
        dydparam[1] = 0.0;
        dydparam[2] = 0.0;
        dydparam[3] = 0.0;
        dydparam[4] = 0.0;
        //double a=min(max(pow(param[0],2.0),Error),1.0/Error), sig=min(max(pow(param[1],2.0),Error),1.0/Error), b=min(max(pow(param[2],2.0),Error),1.0/Error), eta=min(max(pow(param[3],2.0),Error),1.0/Error), rho=min(max(1.0-2.0/(1.0+exp(param[4])),-1.0+Error*100000000),1.0-Error*100000000);
        double a = min(max(pow(param[0], 2.0), err), 1.0 / err), sig = min(max(pow(param[1], 2.0), err), 1.0 / err), b = min(max(pow(param[2], 2.0), err), 1.0 / err), eta = min(max(pow(param[3], 2.0), err), 1.0 / err), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])), -1.0 + err), 1.0 - err);
        for (i = 1; i < ncaplet; i++)
        {

            expatT = min(max(exp(-a * (t[i + 1] - t[i])), err), 1.0 / err);
            exp2at = min(max(exp(-2.0 * a * t[i]), err), 1.0 / err);
            expbtT = min(max(exp(-b * (t[i + 1] - t[i])), err), 1.0 / err);
            exp2bt = min(max(exp(-2.0 * b * t[i]), err), 1.0 / err);
            expapbt = min(max(exp(-(a + b) * t[i]), err), 1.0 / err);

            /*
                        expatT=min(max(exp(-a*(t[i+1]-t[i])),Error),1.0/Error);
                        exp2at=min(max(exp(-2.0*a*t[i]),Error),1.0/Error);
                        expbtT=min(max(exp(-b*(t[i+1]-t[i])),Error),1.0/Error);
                        exp2bt=min(max(exp(-2.0*b*t[i]),Error),1.0/Error);
                        expapbt=min(max(exp(-(a+b)*t[i]),Error),1.0/Error);

                        expatT=exp(-a*(t[i+1]-t[i]));
                        exp2at=exp(-2.0*a*t[i]);
                        expbtT=exp(-b*(t[i+1]-t[i]));
                        exp2bt=exp(-2.0*b*t[i]);
                        expapbt=exp(-(a+b)*t[i]);
            */

            voltT = sqrt(sig * sig / (2.0 * a * a * a) * (1.0 - expatT) * (1.0 - expatT) * (1.0 - exp2at) + eta * eta / (2.0 * b * b * b) * (1.0 - expbtT) * (1.0 - expbtT) * (1.0 - exp2bt) + 2.0 * rho * sig * eta / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt));

            lnPtTX = -log(df[i + 1] * (1.0 + x * tau[i]) / df[i]);

            dvoltTda = (sig * sig / (a * a * a) * (1.0 - expatT) * (-1.5 / a * (1.0 - expatT) * (1.0 - exp2at) + expatT * (t[i + 1] - t[i]) * (1.0 - exp2at) + (1.0 - expatT) * exp2at * t[i]) + 2.0 * rho * sig * eta / (a * b * (a + b)) * (1.0 - expbtT) * (-(1.0 / a + 1.0 / (a + b)) * (1.0 - expatT) * (1.0 - expapbt) + expatT * (t[i + 1] - t[i]) * (1.0 - expapbt) + (1.0 - expatT) * expapbt * t[i])) * param[0] / voltT;

            dvoltTdsig = (sig / (a * a * a) * (1.0 - expatT) * (1.0 - expatT) * (1.0 - exp2at) + 2.0 * rho * eta / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt)) * param[1] / voltT;

            dvoltTdb = (eta * eta / (b * b * b) * (1.0 - expbtT) * (-1.5 / b * (1.0 - expbtT) * (1.0 - exp2bt) + expbtT * (t[i + 1] - t[i]) * (1.0 - exp2bt) + (1.0 - expbtT) * exp2bt * t[i]) + 2.0 * rho * sig * eta / (a * b * (a + b)) * (1.0 - expatT) * (-(1.0 / b + 1.0 / (a + b)) * (1.0 - expbtT) * (1.0 - expapbt) + expbtT * (t[i + 1] - t[i]) * (1.0 - expapbt) + (1.0 - expbtT) * expapbt * t[i])) * param[2] / voltT;

            dvoltTdeta = (eta / (b * b * b) * (1.0 - expbtT) * (1.0 - expbtT) * (1.0 - exp2bt) + 2.0 * rho * sig / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt)) * param[3] / voltT;

            dvoltTdrho = sig * eta / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt) * (1.0 - rho * rho) * 0.5 / voltT;

            y = y + df[i] * cumnormal_Pol_App(lnPtTX / voltT + 0.5 * voltT) - df[i + 1] * (1.0 + x * tau[i]) * cumnormal_Pol_App(lnPtTX / voltT - 0.5 * voltT);

            pid1 = df[i] * max(pdfnormal(lnPtTX / voltT + 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) + 0.5);
            pid2 = df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(lnPtTX / voltT - 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) - 0.5);

            dydparam[0] = dydparam[0] + (pid1 - pid2) * dvoltTda;
            dydparam[1] = dydparam[1] + (pid1 - pid2) * dvoltTdsig;
            dydparam[2] = dydparam[2] + (pid1 - pid2) * dvoltTdb;
            dydparam[3] = dydparam[3] + (pid1 - pid2) * dvoltTdeta;
            dydparam[4] = dydparam[4] + (pid1 - pid2) * dvoltTdrho;
        }
    }
    else
    {
        double maxaaa = 10.0,minaaa = 0.0001;
        double maxd = 10000.0, mind = -10000.0;
        int i, j, nswaplet = int(tau.size());
        double expat, exp2at, expbt, exp2bt, sqrtexp2at, sqrtexp2bt, expapbt, siga, etab, rhoab, rhoab1, mua, mub, Vt, dexpatda, dexp2atda, dexpbtdb, dexp2btdb, dsqrtexp2atda, dsqrtexp2btdb, dexpapbtda, dexpapbtdb, dsigada, dsigadsig, detabdb, detabdeta, drhoabda, drhoabdb, drhoabdrho, drhoab1da, drhoab1db, drhoab1drho, dmuada, dmuadb, dmuadsig, dmuadeta, dmuadrho, dmubda, dmubdb, dmubdsig, dmubdeta, dmubdrho, dVtda, dVtdb, dVtdsig, dVtdeta, dVtdrho, yy, denom, nomdybarda, nomdybardb, nomdybardsig, nomdybardeta, nomdybardrho, tmpexpBbtTybar, dyyda, dyydb, dyydsig, dyydeta, dyydrho;
        y = 0.0;
        dydparam[0] = 0.0;
        dydparam[1] = 0.0;
        dydparam[2] = 0.0;
        dydparam[3] = 0.0;
        dydparam[4] = 0.0;
        double a = min(max(pow(param[0],2.0),minaaa),maxaaa), sig = min(max(pow(param[1],2.0),minaaa),maxaaa), b = min(max(pow(param[2],2.0),minaaa),maxaaa), eta = min(max(pow(param[3],2.0),minaaa),maxaaa), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])),-0.99),0.99);
        vector<double> expaT(nswaplet,0.0), exp2aT(nswaplet,0.0), expatT(nswaplet,0.0), exp2atT(nswaplet,0.0), expbT(nswaplet,0.0), exp2bT(nswaplet,0.0), expbtT(nswaplet,0.0), exp2btT(nswaplet,0.0), expapbT(nswaplet,0.0), expapbtT(nswaplet,0.0), tT(nswaplet,0.0), PMtT(nswaplet,0.0), VT(nswaplet,0.0), VtT(nswaplet,0.0), AtT(nswaplet,0.0), BatT(nswaplet,0.0), c(nswaplet,0.0), BbtT(nswaplet,0.0), ybar(NumGLQ,0.0), xx(NumGLQ,0.0), h1(NumGLQ,0.0), dexpaTda(nswaplet,0.0), dexp2aTda(nswaplet,0.0), dexpatTda(nswaplet,0.0), dexp2atTda(nswaplet,0.0), dexpbTdb(nswaplet,0.0), dexp2bTdb(nswaplet,0.0), dexpbtTdb(nswaplet,0.0), dexp2btTdb(nswaplet,0.0), dexpapbTda(nswaplet,0.0), dexpapbTdb(nswaplet,0.0), dexpapbtTda(nswaplet,0.0), dexpapbtTdb(nswaplet,0.0), dVTda(nswaplet,0.0), dVTdb(nswaplet,0.0), dVTdsig(nswaplet,0.0), dVTdeta(nswaplet,0.0), dVTdrho(nswaplet,0.0), dVtTda(nswaplet,0.0), dVtTdb(nswaplet,0.0), dVtTdsig(nswaplet,0.0), dVtTdeta(nswaplet,0.0), dVtTdrho(nswaplet,0.0), dAtTda(nswaplet,0.0)
            , dAtTdb(nswaplet,0.0), dAtTdsig(nswaplet,0.0), dAtTdeta(nswaplet,0.0), dAtTdrho(nswaplet,0.0), dBatTda(nswaplet,0.0), dBbtTdb(nswaplet,0.0), dybarda(NumGLQ,0.0), dybardb(NumGLQ,0.0), dybardsig(NumGLQ,0.0), dybardeta(NumGLQ,0.0), dybardrho(NumGLQ,0.0), dxxda(NumGLQ,0.0), dxxdsig(NumGLQ,0.0), dh1da(NumGLQ,0.0), dh1db(NumGLQ,0.0), dh1dsig(NumGLQ,0.0), dh1deta(NumGLQ,0.0), dh1drho(NumGLQ,0.0);//dxxdb(NumGLQ), dxxdeta(NumGLQ), dxxdrho(NumGLQ),
        
        vector<vector<double>> lambij(nswaplet), expkappij(nswaplet), h2(nswaplet), dlambijda(nswaplet), dlambijdb(nswaplet), dlambijdsig(nswaplet), dlambijdeta(nswaplet), dlambijdrho(nswaplet), dexpkappijda(nswaplet), dexpkappijdb(nswaplet), dexpkappijdsig(nswaplet), dexpkappijdeta(nswaplet), dexpkappijdrho(nswaplet), dh2da(nswaplet), dh2db(nswaplet), dh2dsig(nswaplet), dh2deta(nswaplet), dh2drho(nswaplet);

        expat = exp(-a * t[0]);
        dexpatda = expat * (-t[0]);

        exp2at = expat * expat;
        dexp2atda = 2.0 * expat * dexpatda;

        sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
        dsqrtexp2atda = 0.25 * (a * dexp2atda - (1.0 - exp2at)) / (a * a * sqrtexp2at * sqrtexp2at * sqrtexp2at);

        expbt = exp(-b * t[0]);
        dexpbtdb = expbt * (-t[0]);

        exp2bt = expbt * expbt;
        dexp2btdb = 2.0 * expbt * dexpbtdb;

        sqrtexp2bt = sqrt((1.0 - exp2bt) / (2.0 * b));
        dsqrtexp2btdb = 0.25 * (b * dexp2btdb - (1.0 - exp2bt)) / (b * b * sqrtexp2bt * sqrtexp2bt * sqrtexp2bt);

        expapbt = exp(-(a + b) * t[0]);
        dexpapbtda = expapbt * (-t[0]);
        dexpapbtdb = dexpapbtda;

        siga = sig * sqrtexp2at;
        dsigada = sig * dsqrtexp2atda;
        dsigadsig = sqrtexp2at;

        etab = eta * sqrtexp2bt;
        detabdb = eta * dsqrtexp2btdb;
        detabdeta = sqrtexp2bt;

        rhoab = rho * (1.0 - expapbt) / (sqrtexp2at * sqrtexp2bt) / (a + b);
        drhoabda = rho / (sqrtexp2at * sqrtexp2bt) / (a + b) * (-(1.0 - expapbt) * (1.0 / (a + b) + 1.0 / sqrtexp2at * dsqrtexp2atda) - dexpapbtda);
        drhoabdb = rho / (sqrtexp2at * sqrtexp2bt) / (a + b) * (-(1.0 - expapbt) * (1.0 / (a + b) + 1.0 / sqrtexp2bt * dsqrtexp2btdb) - dexpapbtdb);
        drhoabdrho = (1.0 - expapbt) / (sqrtexp2at * sqrtexp2bt) / (a + b);

        rhoab1 = sqrt(1.0 - rhoab * rhoab);
        drhoab1da = -rhoab / rhoab1 * drhoabda;
        drhoab1db = -rhoab / rhoab1 * drhoabdb;
        drhoab1drho = -rhoab / rhoab1 * drhoabdrho;

        mua = -(sig * sig / (a * a) + rho * sig * eta / (a * b)) * (1.0 - expat) + 0.5 * sig * sig / (a * a) * (1.0 - exp2at) + rho * sig * eta / (b * (a + b)) * (1.0 - expapbt);
        dmuada = (2.0 * sig * sig / (a * a * a) + rho * sig * eta / (a * a * b)) * (1.0 - expat) + (sig * sig / (a * a) + rho * sig * eta / (a * b)) * dexpatda - sig * sig / (a * a * a) * (1.0 - exp2at) - 0.5 * sig * sig / (a * a) * dexp2atda - rho * sig * eta / (b * (a + b) * (a + b)) * (1.0 - expapbt) + rho * sig * eta / (b * (a + b)) * (-dexpapbtda);
        dmuadb = rho * sig * eta / (a * b * b) * (1.0 - expat) - rho * sig * eta / (b * (a + b)) * (1.0 / b + 1.0 / (a + b)) * (1.0 - expapbt) + rho * sig * eta / (b * (a + b)) * (-dexpapbtdb);
        dmuadsig = -(2.0 * sig / (a * a) + rho * eta / (a * b)) * (1.0 - expat) + sig / (a * a) * (1.0 - exp2at) + rho * eta / (b * (a + b)) * (1.0 - expapbt);
        dmuadeta = -rho * sig / (a * b) * (1.0 - expat) + rho * sig / (b * (a + b)) * (1.0 - expapbt);
        dmuadrho = -sig * eta / (a * b) * (1.0 - expat) + sig * eta / (b * (a + b)) * (1.0 - expapbt);

        mub = -(eta * eta / (b * b) + rho * sig * eta / (a * b)) * (1.0 - expbt) + 0.5 * eta * eta / (b * b) * (1.0 - exp2bt) + rho * sig * eta / (a * (a + b)) * (1.0 - expapbt);
        dmubda = rho * sig * eta / (a * a * b) * (1.0 - expbt) - rho * sig * eta / (a * (a + b)) * (1.0 / a + 1.0 / (a + b)) * (1.0 - expapbt) + rho * sig * eta / (a * (a + b)) * (-dexpapbtda);
        dmubdb = (2.0 * eta * eta / (b * b * b) + rho * sig * eta / (a * b * b)) * (1.0 - expbt) + (eta * eta / (b * b) + rho * sig * eta / (a * b)) * dexpbtdb - eta * eta / (b * b * b) * (1.0 - exp2bt) - 0.5 * eta * eta / (b * b) * dexp2btdb - rho * sig * eta / (a * (a + b) * (a + b)) * (1.0 - expapbt) + rho * sig * eta / (a * (a + b)) * (-dexpapbtdb);
        dmubdsig = -rho * eta / (a * b) * (1.0 - expbt) + rho * eta / (a * (a + b)) * (1.0 - expapbt);
        dmubdeta = -(2.0 * eta / (b * b) + rho * sig / (a * b)) * (1.0 - expbt) + eta / (b * b) * (1.0 - exp2bt) + rho * sig / (a * (a + b)) * (1.0 - expapbt);
        dmubdrho = -sig * eta / (a * b) * (1.0 - expbt) + sig * eta / (a * (a + b)) * (1.0 - expapbt);

        Vt = sig * sig / (a * a) * (t[0] + 2.0 / a * expat - 0.5 / a * exp2at - 1.5 / a) + eta * eta / (b * b) * (t[0] + 2.0 / b * expbt - 0.5 / b * exp2bt - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));
        dVtda = -2.0 * sig * sig / (a * a * a) * (t[0] + 2.0 / a * expat - 0.5 / a * exp2at - 1.5 / a) + sig * sig / (a * a) * (-2.0 / (a * a) * expat + 2.0 / a * dexpatda + 0.5 / (a * a) * exp2at - 0.5 / a * dexp2atda + 1.5 / (a * a)) - 2.0 * rho * sig * eta / (a * a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (-(expat - 1.0) / (a * a) + dexpatda / a + (expapbt - 1.0) / ((a + b) * (a + b)) - dexpapbtda / (a + b));
        dVtdb = -2.0 * eta * eta / (b * b * b) * (t[0] + 2.0 / b * expbt - 0.5 / b * exp2bt - 1.5 / b) + eta * eta / (b * b) * (-2.0 / (b * b) * expbt + 2.0 / b * dexpbtdb + 0.5 / (b * b) * exp2bt - 0.5 / b * dexp2btdb + 1.5 / (b * b)) - 2.0 * rho * sig * eta / (a * b * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (-(expbt - 1.0) / (b * b) + dexpbtdb / b + (expapbt - 1.0) / ((a + b) * (a + b)) - dexpapbtdb / (a + b));
        dVtdsig = 2.0 * sig / (a * a) * (t[0] + 2.0 / a * expat - 0.5 / a * exp2at - 1.5 / a) + 2.0 * rho * eta / (a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));
        dVtdeta = 2.0 * eta / (b * b) * (t[0] + 2.0 / b * expbt - 0.5 / b * exp2bt - 1.5 / b) + 2.0 * rho * sig / (a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));
        dVtdrho = 2.0 * sig * eta / (a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));

        for (j = 0; j < NumGLQ; j++)
        {
            xx[j] = M_SQRT2 * GLQx[j] * siga + mua;
            dxxda[j] = M_SQRT2 * GLQx[j] * dsigada + dmuada;
            dxxdsig[j] = M_SQRT2 * GLQx[j] * dsigadsig + dmuadsig;
        }
        //dxxdb=dmuadb;
        //dxxdeta=dmuadeta;
        //dxxdrho=dmuadrho;

        for (i = 0; i < nswaplet; i++)
        {
            expaT[i] = exp(-a * t[i + 1]);
            dexpaTda[i] = expaT[i] * (-t[i + 1]);

            expbT[i] = exp(-b * t[i + 1]);
            dexpbTdb[i] = expbT[i] * (-t[i + 1]);

            exp2aT[i] = expaT[i] * expaT[i];
            dexp2aTda[i] = 2.0 * expaT[i] * dexpaTda[i];

            exp2bT[i] = expbT[i] * expbT[i];
            dexp2bTdb[i] = 2.0 * expbT[i] * dexpbTdb[i];

            expapbT[i] = exp(-(a + b) * t[i + 1]);
            dexpapbTda[i] = expapbT[i] * (-t[i + 1]);
            dexpapbTdb[i] = dexpapbTda[i];

            tT[i] = t[i + 1] - t[0];

            expatT[i] = exp(-a * tT[i]);
            dexpatTda[i] = expatT[i] * (-tT[i]);

            expbtT[i] = exp(-b * tT[i]);
            dexpbtTdb[i] = expbtT[i] * (-tT[i]);

            exp2atT[i] = expatT[i] * expatT[i];
            dexp2atTda[i] = 2.0 * expatT[i] * dexpatTda[i];

            exp2btT[i] = expbtT[i] * expbtT[i];
            dexp2btTdb[i] = 2.0 * expbtT[i] * dexpbtTdb[i];

            expapbtT[i] = exp(-(a + b) * tT[i]);
            dexpapbtTda[i] = expapbtT[i] * (-tT[i]);
            dexpapbtTdb[i] = dexpapbtTda[i];

            PMtT[i] = df[i + 1] / df[0];

            VT[i] = sig * sig / (a * a) * (t[i + 1] + 2.0 / a * expaT[i] - 0.5 / a * exp2aT[i] - 1.5 / a) + eta * eta / (b * b) * (t[i + 1] + 2.0 / b * expbT[i] - 0.5 / b * exp2bT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));
            dVTda[i] = -2.0 * sig * sig / (a * a * a) * (t[i + 1] + 2.0 / a * expaT[i] - 0.5 / a * exp2aT[i] - 1.5 / a) + sig * sig / (a * a) * (-2.0 / (a * a) * expaT[i] + 2.0 / a * dexpaTda[i] + 0.5 / (a * a) * exp2aT[i] - 0.5 / a * dexp2aTda[i] + 1.5 / (a * a)) - 2.0 * rho * sig * eta / (a * a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (dexpaTda[i] / a - (expaT[i] - 1.0) / (a * a) - dexpapbTda[i] / (a + b) + (expapbT[i] - 1.0) / ((a + b) * (a + b)));
            dVTdb[i] = -2.0 * eta * eta / (b * b * b) * (t[i + 1] + 2.0 / b * expbT[i] - 0.5 / b * exp2bT[i] - 1.5 / b) + eta * eta / (b * b) * (-2.0 / (b * b) * expbT[i] + 2.0 / b * dexpbTdb[i] + 0.5 / (b * b) * exp2bT[i] - 0.5 / b * dexp2bTdb[i] + 1.5 / (b * b)) - 2.0 * rho * sig * eta / (a * b * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (dexpbTdb[i] / b - (expbT[i] - 1.0) / (b * b) - dexpapbTdb[i] / (a + b) + (expapbT[i] - 1.0) / ((a + b) * (a + b)));
            dVTdsig[i] = 2.0 * sig / (a * a) * (t[i + 1] + 2.0 / a * expaT[i] - 0.5 / a * exp2aT[i] - 1.5 / a) + 2.0 * rho * eta / (a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));
            dVTdeta[i] = 2.0 * eta / (b * b) * (t[i + 1] + 2.0 / b * expbT[i] - 0.5 / b * exp2bT[i] - 1.5 / b) + 2.0 * rho * sig / (a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));
            dVTdrho[i] = 2.0 * sig * eta / (a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));

            VtT[i] = sig * sig / (a * a) * (tT[i] + 2.0 / a * expatT[i] - 0.5 / a * exp2atT[i] - 1.5 / a) + eta * eta / (b * b) * (tT[i] + 2.0 / b * expbtT[i] - 0.5 / b * exp2btT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));
            dVtTda[i] = -2.0 * sig * sig / (a * a * a) * (tT[i] + 2.0 / a * expatT[i] - 0.5 / a * exp2atT[i] - 1.5 / a) + sig * sig / (a * a) * (-2.0 / (a * a) * expatT[i] + 2.0 / a * dexpatTda[i] + 0.5 / (a * a) * exp2atT[i] - 0.5 / a * dexp2atTda[i] + 1.5 / (a * a)) - 2.0 * rho * sig * eta / (a * a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (dexpatTda[i] / a - (expatT[i] - 1.0) / (a * a) - dexpapbtTda[i] / (a + b) + (expapbtT[i] - 1.0) / ((a + b) * (a + b)));
            dVtTdb[i] = -2.0 * eta * eta / (b * b * b) * (tT[i] + 2.0 / b * expbtT[i] - 0.5 / b * exp2btT[i] - 1.5 / b) + eta * eta / (b * b) * (-2.0 / (b * b) * expbtT[i] + 2.0 / b * dexpbtTdb[i] + 0.5 / (b * b) * exp2btT[i] - 0.5 / b * dexp2btTdb[i] + 1.5 / (b * b)) - 2.0 * rho * sig * eta / (a * b * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b)) + 2.0 * rho * sig * eta / (a * b) * (dexpbtTdb[i] / b - (expbtT[i] - 1.0) / (b * b) - dexpapbtTdb[i] / (a + b) + (expapbtT[i] - 1.0) / ((a + b) * (a + b)));
            dVtTdsig[i] = 2.0 * sig / (a * a) * (tT[i] + 2.0 / a * expatT[i] - 0.5 / a * exp2atT[i] - 1.5 / a) + 2.0 * rho * eta / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));
            dVtTdeta[i] = 2.0 * eta / (b * b) * (tT[i] + 2.0 / b * expbtT[i] - 0.5 / b * exp2btT[i] - 1.5 / b) + 2.0 * rho * sig / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));
            dVtTdrho[i] = 2.0 * sig * eta / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));

            c[i] = x * tau[i];

            BatT[i] = (1.0 - expatT[i]) / a;
            dBatTda[i] = -(dexpatTda[i] + BatT[i]) / a;

            BbtT[i] = (1.0 - expbtT[i]) / b;
            dBbtTdb[i] = -(dexpbtTdb[i] + BbtT[i]) / b;

            AtT[i] = min(max(PMtT[i] * exp(0.5 * (VtT[i] - VT[i] + Vt)),minaaa),maxaaa);
            dAtTda[i] = 0.5 * AtT[i] * (dVtTda[i] - dVTda[i] + dVtda);
            dAtTdb[i] = 0.5 * AtT[i] * (dVtTdb[i] - dVTdb[i] + dVtdb);
            dAtTdsig[i] = 0.5 * AtT[i] * (dVtTdsig[i] - dVTdsig[i] + dVtdsig);
            dAtTdeta[i] = 0.5 * AtT[i] * (dVtTdeta[i] - dVTdeta[i] + dVtdeta);
            dAtTdrho[i] = 0.5 * AtT[i] * (dVtTdrho[i] - dVTdrho[i] + dVtdrho);

            lambij[i] = vector<double>(NumGLQ);
            expkappij[i] = vector<double>(NumGLQ);
            h2[i] = vector<double>(NumGLQ);
            dlambijda[i] = vector<double>(NumGLQ);
            dlambijdb[i] = vector<double>(NumGLQ);
            dlambijdsig[i] = vector<double>(NumGLQ);
            dlambijdeta[i] = vector<double>(NumGLQ);
            dlambijdrho[i] = vector<double>(NumGLQ);
            dexpkappijda[i] = vector<double>(NumGLQ);
            dexpkappijdb[i] = vector<double>(NumGLQ);
            dexpkappijdsig[i] = vector<double>(NumGLQ);
            dexpkappijdeta[i] = vector<double>(NumGLQ);
            dexpkappijdrho[i] = vector<double>(NumGLQ);
            dh2da[i] = vector<double>(NumGLQ);
            dh2db[i] = vector<double>(NumGLQ);
            dh2dsig[i] = vector<double>(NumGLQ);
            dh2deta[i] = vector<double>(NumGLQ);
            dh2drho[i] = vector<double>(NumGLQ);


            for (j = 0; j < NumGLQ; j++)
            {
                lambij[i][j] = c[i] * AtT[i] * exp(-BatT[i] * xx[j]);
                dlambijda[i][j] = c[i] * (dAtTda[i] * exp(-BatT[i] * xx[j]) + AtT[i] * exp(-BatT[i] * xx[j]) * (-dBatTda[i] * xx[j] - BatT[i] * dxxda[j]));
                dlambijdb[i][j] = c[i] * (dAtTdb[i] * exp(-BatT[i] * xx[j]) + AtT[i] * exp(-BatT[i] * xx[j]) * (-BatT[i] * dmuadb));
                dlambijdsig[i][j] = c[i] * (dAtTdsig[i] * exp(-BatT[i] * xx[j]) + AtT[i] * exp(-BatT[i] * xx[j]) * (-BatT[i] * dxxdsig[j]));
                dlambijdeta[i][j] = c[i] * (dAtTdeta[i] * exp(-BatT[i] * xx[j]) + AtT[i] * exp(-BatT[i] * xx[j]) * (-BatT[i] * dmuadeta));
                dlambijdrho[i][j] = c[i] * (dAtTdrho[i] * exp(-BatT[i] * xx[j]) + AtT[i] * exp(-BatT[i] * xx[j]) * (-BatT[i] * dmuadrho));

                expkappij[i][j] = exp(-BbtT[i] * (mub - 0.5 * rhoab1 * rhoab1 * etab * etab * BbtT[i] + rhoab * etab * GLQx[j] * M_SQRT2));
                dexpkappijda[i][j] = expkappij[i][j] * (-BbtT[i] * (dmubda - rhoab1 * drhoab1da * etab * etab * BbtT[i] + drhoabda * etab * GLQx[j] * M_SQRT2));
                dexpkappijdb[i][j] = expkappij[i][j] * (-dBbtTdb[i] * (mub - 0.5 * rhoab1 * rhoab1 * etab * etab * BbtT[i] + rhoab * etab * GLQx[j] * M_SQRT2) - BbtT[i] * (dmubdb - rhoab1 * drhoab1db * etab * etab * BbtT[i] - rhoab1 * rhoab1 * etab * detabdb * BbtT[i] - 0.5 * rhoab1 * rhoab1 * etab * etab * dBbtTdb[i] + drhoabdb * etab * GLQx[j] * M_SQRT2 + rhoab * detabdb * GLQx[j] * M_SQRT2));
                dexpkappijdsig[i][j] = expkappij[i][j] * (-BbtT[i] * dmubdsig);
                dexpkappijdeta[i][j] = expkappij[i][j] * (-BbtT[i] * (dmubdeta - rhoab1 * rhoab1 * etab * detabdeta * BbtT[i] + rhoab * detabdeta * GLQx[j] * M_SQRT2));
                dexpkappijdrho[i][j] = expkappij[i][j] * (-BbtT[i] * (dmubdrho - rhoab1 * drhoab1drho * etab * etab * BbtT[i] + drhoabdrho * etab * GLQx[j] * M_SQRT2));
            }
        }
        for (j = 0; j < NumGLQ; j++)
        {
            lambij[nswaplet - 1][j] = lambij[nswaplet - 1][j] + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]);
            dlambijda[nswaplet - 1][j] = dlambijda[nswaplet - 1][j] + dAtTda[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) * (-dBatTda[nswaplet - 1] * xx[j] - BatT[nswaplet - 1] * dxxda[j]);
            dlambijdb[nswaplet - 1][j] = dlambijdb[nswaplet - 1][j] + dAtTdb[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) * (-BatT[nswaplet - 1] * dmuadb);
            dlambijdsig[nswaplet - 1][j] = dlambijdsig[nswaplet - 1][j] + dAtTdsig[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) * (-BatT[nswaplet - 1] * dxxdsig[j]);
            dlambijdeta[nswaplet - 1][j] = dlambijdeta[nswaplet - 1][j] + dAtTdeta[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) * (-BatT[nswaplet - 1] * dmuadeta);
            dlambijdrho[nswaplet - 1][j] = dlambijdrho[nswaplet - 1][j] + dAtTdrho[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]) * (-BatT[nswaplet - 1] * dmuadrho);
        }
        ybar_HW2F(lambij,BbtT,ybar);
        for (j = 0; j < NumGLQ; j++)
        {
            h1[j] = (ybar[j] - mub) / (etab * rhoab1) - rhoab * GLQx[j] * M_SQRT2 / rhoab1;
            dh1da[j] = (dybarda[j] - dmubda) / (etab * rhoab1) - (ybar[j] - mub) / (etab * rhoab1 * rhoab1) * drhoab1da - drhoabda * GLQx[j] * M_SQRT2 / rhoab1 + rhoab * GLQx[j] * M_SQRT2 / (rhoab1 * rhoab1) * drhoab1da;
            dh1db[j] = (dybardb[j] - dmubdb) / (etab * rhoab1) - (ybar[j] - mub) / (etab * etab * rhoab1) * detabdb - (ybar[j] - mub) / (etab * rhoab1 * rhoab1) * drhoab1db - drhoabdb * GLQx[j] * M_SQRT2 / rhoab1 + rhoab * GLQx[j] * M_SQRT2 / (rhoab1 * rhoab1) * drhoab1db;
            dh1dsig[j] = (dybardsig[j] - dmubdsig) / (etab * rhoab1);
            dh1deta[j] = (dybardeta[j] - dmubdeta) / (etab * rhoab1) - (ybar[j] - mub) / (etab * etab * rhoab1) * detabdeta;
            dh1drho[j] = (dybardrho[j] - dmubdrho) / (etab * rhoab1) - (ybar[j] - mub) / (etab * rhoab1 * rhoab1) * drhoab1drho - drhoabdrho * GLQx[j] * M_SQRT2 / rhoab1 + rhoab * GLQx[j] * M_SQRT2 / (rhoab1 * rhoab1) * drhoab1drho;

            denom = 0.0;
            nomdybarda = 0.0;
            nomdybardb = 0.0;
            nomdybardsig = 0.0;
            nomdybardeta = 0.0;
            nomdybardrho = 0.0;

            for (i = 0; i < nswaplet; i++)
            {
                tmpexpBbtTybar = exp(-BbtT[i] * ybar[j]);
                denom = denom + lambij[i][j] * tmpexpBbtTybar * BbtT[i];
                nomdybarda = nomdybarda + tmpexpBbtTybar * dlambijda[i][j];
                nomdybardb = nomdybardb + tmpexpBbtTybar * (dlambijdb[i][j] - lambij[i][j] * dBbtTdb[i] * ybar[j]);
                nomdybardsig = nomdybardsig + tmpexpBbtTybar * dlambijdsig[i][j];
                nomdybardeta = nomdybardeta + tmpexpBbtTybar * dlambijdeta[i][j];
                nomdybardrho = nomdybardrho + tmpexpBbtTybar * dlambijdrho[i][j];
            }

            dybarda[j] = nomdybarda / denom;
            dybardb[j] = nomdybardb / denom;
            dybardsig[j] = nomdybardsig / denom;
            dybardeta[j] = nomdybardeta / denom;
            dybardrho[j] = nomdybardrho / denom;

            yy = 0.0;
            dyyda = 0.0;
            dyydb = 0.0;
            dyydsig = 0.0;
            dyydeta = 0.0;
            dyydrho = 0.0;

            for (i = 0; i < nswaplet; i++)
            {
                h2[i][j] = h1[j] + BbtT[i] * etab * rhoab1;
                dh2da[i][j] = dh1da[j] + BbtT[i] * etab * drhoab1da;
                dh2db[i][j] = dh1db[j] + dBbtTdb[i] * etab * rhoab1 + BbtT[i] * detabdb * rhoab1 + BbtT[i] * etab * drhoab1db;
                dh2dsig[i][j] = dh1dsig[j];
                dh2deta[i][j] = dh1deta[j] + BbtT[i] * detabdeta * rhoab1;
                dh2drho[i][j] = dh1drho[j] + BbtT[i] * etab * drhoab1drho;

                yy = yy - lambij[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]);
                dyyda = dyyda - dlambijda[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * dexpkappijda[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * expkappij[i][j] * max(pdfnormal(-h2[i][j]),Error) * (-dh2da[i][j]);
                dyydb = dyydb - dlambijdb[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * dexpkappijdb[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * expkappij[i][j] * max(pdfnormal(-h2[i][j]),Error) * (-dh2db[i][j]);
                dyydsig = dyydsig - dlambijdsig[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * dexpkappijdsig[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * expkappij[i][j] * max(pdfnormal(-h2[i][j]),Error) * (-dh2dsig[i][j]);
                dyydeta = dyydeta - dlambijdeta[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * dexpkappijdeta[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * expkappij[i][j] * max(pdfnormal(-h2[i][j]),Error) * (-dh2deta[i][j]);
                dyydrho = dyydrho - dlambijdrho[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * dexpkappijdrho[i][j] * cumnormal_Pol_App(-h2[i][j]) - lambij[i][j] * expkappij[i][j] * max(pdfnormal(-h2[i][j]),Error) * (-dh2drho[i][j]);
            }

            yy = yy + cumnormal_Pol_App(-h1[j]);
            dyyda = dyyda + max(pdfnormal(-h1[j]),Error) * (-dh1da[j]);
            dyydb = dyydb + max(pdfnormal(-h1[j]),Error) * (-dh1db[j]);
            dyydsig = dyydsig + max(pdfnormal(-h1[j]),Error) * (-dh1dsig[j]);
            dyydeta = dyydeta + max(pdfnormal(-h1[j]),Error) * (-dh1deta[j]);
            dyydrho = dyydrho + max(pdfnormal(-h1[j]),Error) * (-dh1drho[j]);

            y = y + GLQw[j] * yy;
            dydparam[0] = dydparam[0] + GLQw[j] * dyyda;
            dydparam[1] = dydparam[1] + GLQw[j] * dyydb;
            dydparam[2] = dydparam[2] + GLQw[j] * dyydsig;
            dydparam[3] = dydparam[3] + GLQw[j] * dyydeta;
            dydparam[4] = dydparam[4] + GLQw[j] * dyydrho;
        }
        y = y * df[0] / M_SQRTPI;
        dydparam[0] = min(max(dydparam[0] * df[0] / M_SQRTPI * param[0] * 2.0,-10000.0),10000.0);
        dydparam[1] = min(max(dydparam[1] * df[0] / M_SQRTPI * param[1] * 2.0,-10000.0),10000.0);
        dydparam[2] = min(max(dydparam[2] * df[0] / M_SQRTPI * param[2] * 2.0,-10000.0),10000.0);
        dydparam[3] = min(max(dydparam[3] * df[0] / M_SQRTPI * param[3] * 2.0,-10000.0),10000.0);
        dydparam[4] = min(max(dydparam[4] * df[0] / M_SQRTPI * (1.0 - rho * rho) * 0.5,-10000.0),10000.0);
    }
}

///////////////////////////////////////////////////////////////////////////////
void PlainIROption_HW2F(const string flag, const vector<double> df, const vector<double> t, const vector<double> tau, const double fM, const double x, vector<double>& param, double& y)
{
    double err = Error * (1.0e+10);//1.0e-5;
    if (flag == "CAP")
    {
        int i, ncaplet = int(tau.size());
        double expatT, exp2at, expbtT, exp2bt, expapbt, voltT, lnPtTX, pid1, pid2;
        y = 0.0;
        double a = min(max(pow(param[0], 2.0), err), 1.0 / err), sig = min(max(pow(param[1], 2.0), err), 1.0 / err), b = min(max(pow(param[2], 2.0), err), 1.0 / err), eta = min(max(pow(param[3], 2.0), err), 1.0 / err), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])), -1.0 + err), 1.0 - err);
        for (i = 1; i < ncaplet; i++)
        {

            expatT = min(max(exp(-a * (t[i + 1] - t[i])), err), 1.0 / err);
            exp2at = min(max(exp(-2.0 * a * t[i]), err), 1.0 / err);
            expbtT = min(max(exp(-b * (t[i + 1] - t[i])), err), 1.0 / err);
            exp2bt = min(max(exp(-2.0 * b * t[i]), err), 1.0 / err);
            expapbt = min(max(exp(-(a + b) * t[i]), err), 1.0 / err);

            voltT = sqrt(sig * sig / (2.0 * a * a * a) * (1.0 - expatT) * (1.0 - expatT) * (1.0 - exp2at) + eta * eta / (2.0 * b * b * b) * (1.0 - expbtT) * (1.0 - expbtT) * (1.0 - exp2bt) + 2.0 * rho * sig * eta / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt));

            lnPtTX = -log(df[i + 1] * (1.0 + x * tau[i]) / df[i]);

            y = y + df[i] * cumnormal_Pol_App(lnPtTX / voltT + 0.5 * voltT) - df[i + 1] * (1.0 + x * tau[i]) * cumnormal_Pol_App(lnPtTX / voltT - 0.5 * voltT);

            pid1 = df[i] * max(pdfnormal(lnPtTX / voltT + 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) + 0.5);
            pid2 = df[i + 1] * (1.0 + x * tau[i]) * max(pdfnormal(lnPtTX / voltT - 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) - 0.5);

        }
    }
    else
    {
        double maxaaa = 10.0,minaaa = 0.0001;
        double maxd = 10000.0, mind = -10000.0;
        int i, j, nswaplet = int(tau.size());
        double expat, exp2at, expbt, exp2bt, sqrtexp2at, sqrtexp2bt, expapbt, siga, etab, rhoab, rhoab1, mua, mub, Vt, yy;
        y = 0.0;
        double a = min(max(pow(param[0],2.0),minaaa),maxaaa), sig = min(max(pow(param[1],2.0),minaaa),maxaaa), b = min(max(pow(param[2],2.0),minaaa),maxaaa), eta = min(max(pow(param[3],2.0),minaaa),maxaaa), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])),-0.99),0.99);
        vector<double> expaT(nswaplet,0.0), exp2aT(nswaplet,0.0), expatT(nswaplet,0.0), exp2atT(nswaplet,0.0), expbT(nswaplet,0.0), exp2bT(nswaplet,0.0), expbtT(nswaplet,0.0), exp2btT(nswaplet,0.0), expapbT(nswaplet,0.0), expapbtT(nswaplet,0.0), tT(nswaplet,0.0), PMtT(nswaplet,0.0), VT(nswaplet,0.0), VtT(nswaplet,0.0), AtT(nswaplet,0.0), BatT(nswaplet,0.0), c(nswaplet,0.0), BbtT(nswaplet,0.0), ybar(NumGLQ,0.0), xx(NumGLQ,0.0), h1(NumGLQ,0.0);
        vector<vector<double>> lambij(nswaplet), expkappij(nswaplet), h2(nswaplet);

        expat = exp(-a * t[0]);
        exp2at = expat * expat;
        sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
        expbt = exp(-b * t[0]);
        exp2bt = expbt * expbt;
        sqrtexp2bt = sqrt((1.0 - exp2bt) / (2.0 * b));
        expapbt = exp(-(a + b) * t[0]);
        siga = sig * sqrtexp2at;
        etab = eta * sqrtexp2bt;
        rhoab = rho * (1.0 - expapbt) / (sqrtexp2at * sqrtexp2bt) / (a + b);
        rhoab1 = sqrt(1.0 - rhoab * rhoab);

        mua = -(sig * sig / (a * a) + rho * sig * eta / (a * b)) * (1.0 - expat) + 0.5 * sig * sig / (a * a) * (1.0 - exp2at) + rho * sig * eta / (b * (a + b)) * (1.0 - expapbt);
        mub = -(eta * eta / (b * b) + rho * sig * eta / (a * b)) * (1.0 - expbt) + 0.5 * eta * eta / (b * b) * (1.0 - exp2bt) + rho * sig * eta / (a * (a + b)) * (1.0 - expapbt);

        Vt = sig * sig / (a * a) * (t[0] + 2.0 / a * expat - 0.5 / a * exp2at - 1.5 / a) + eta * eta / (b * b) * (t[0] + 2.0 / b * expbt - 0.5 / b * exp2bt - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));

        for (j = 0; j < NumGLQ; j++)
        {
            xx[j] = M_SQRT2 * GLQx[j] * siga + mua;
        }

        for (i = 0; i < nswaplet; i++)
        {
            expaT[i] = exp(-a * t[i + 1]);
            expbT[i] = exp(-b * t[i + 1]);
            exp2aT[i] = expaT[i] * expaT[i];
            exp2bT[i] = expbT[i] * expbT[i];
            expapbT[i] = exp(-(a + b) * t[i + 1]);

            tT[i] = t[i + 1] - t[0];

            expatT[i] = exp(-a * tT[i]);
            expbtT[i] = exp(-b * tT[i]);
            exp2atT[i] = expatT[i] * expatT[i];
            exp2btT[i] = expbtT[i] * expbtT[i];
            expapbtT[i] = exp(-(a + b) * tT[i]);

            PMtT[i] = df[i + 1] / df[0];

            VT[i] = sig * sig / (a * a) * (t[i + 1] + 2.0 / a * expaT[i] - 0.5 / a * exp2aT[i] - 1.5 / a) + eta * eta / (b * b) * (t[i + 1] + 2.0 / b * expbT[i] - 0.5 / b * exp2bT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));
            VtT[i] = sig * sig / (a * a) * (tT[i] + 2.0 / a * expatT[i] - 0.5 / a * exp2atT[i] - 1.5 / a) + eta * eta / (b * b) * (tT[i] + 2.0 / b * expbtT[i] - 0.5 / b * exp2btT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));

            c[i] = x * tau[i];

            BatT[i] = (1.0 - expatT[i]) / a;
            BbtT[i] = (1.0 - expbtT[i]) / b;
            AtT[i] = min(max(PMtT[i] * exp(0.5 * (VtT[i] - VT[i] + Vt)),minaaa),maxaaa);

            lambij[i] = vector<double>(NumGLQ);
            expkappij[i] = vector<double>(NumGLQ);
            h2[i] = vector<double>(NumGLQ);

            for (j = 0; j < NumGLQ; j++)
            {
                lambij[i][j] = c[i] * AtT[i] * exp(-BatT[i] * xx[j]);
                expkappij[i][j] = exp(-BbtT[i] * (mub - 0.5 * rhoab1 * rhoab1 * etab * etab * BbtT[i] + rhoab * etab * GLQx[j] * M_SQRT2));
            }
        }
        for (j = 0; j < NumGLQ; j++)
        {
            lambij[nswaplet - 1][j] = lambij[nswaplet - 1][j] + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]);
        }
        ybar_HW2F(lambij,BbtT,ybar);
        for (j = 0; j < NumGLQ; j++)
        {
            h1[j] = (ybar[j] - mub) / (etab * rhoab1) - rhoab * GLQx[j] * M_SQRT2 / rhoab1;

            yy = 0.0;
            for (i = 0; i < nswaplet; i++)
            {
                h2[i][j] = h1[j] + BbtT[i] * etab * rhoab1;
                yy = yy - lambij[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]);
            }
            yy = yy + cumnormal_Pol_App(-h1[j]);
            y = y + GLQw[j] * yy;
        }
        y = y * df[0] / M_SQRTPI;
    }
}



///////////////////////////////////////////////////////////////////////////////


void mrqcof_HW2F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& alpha, vector<double>& beta, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&))
{
    int i, j, k, l, m, mfit = 0, ndata = int(x.size()), ma = int(a.size());
    double ymod, wt, sig2i, dy;
    vector<double> dyda(ma);
    for (j = 0; j < ma; j++) if (ia[j]) mfit++;
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
        beta[j] = 0.0;
    }
    chisq = 0.0;
    for (i = 0; i < ndata; i++)
    {
        funcs(flag[i], df[i], t[i], tau[i], fM[i], x[i], a, ymod, dyda);
        sig2i = 1.0 / (sig[i] * sig[i]);
        dy = y[i] - ymod;
        for (j = 0, l = 0; l < ma; l++)
        {
            if (ia[l])
            {
                wt = dyda[l] * sig2i;
                for (k = 0, m = 0; m < l + 1; m++) if (ia[m]) alpha[j][k++] += wt * dyda[m];
                beta[j++] += dy * wt;
            }
        }
        chisq += dy * dy * sig2i;
    }
    for (j = 1; j < mfit; j++) for (k = 0; k < j; k++) alpha[k][j] = alpha[j][k];

}

void mrqmin_HW2F(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& fM, vector<double>& x, vector<double>& y, vector<double>& sig, vector<double>& a, vector<bool>& ia, vector<vector<double>>& covar, vector<vector<double>>& alpha, double& chisq, void funcs(const string, const vector<double>, const vector<double>, const vector<double>, const double, const double, vector<double>&, double&, vector<double>&), double& alamda, bool& singularflag)
{
    static int mfit;
    static double ochisq;
    int j, k, l, ma = int(a.size());
    static vector<vector<double>> oneda(ma);
    for (j = 0; j < ma; j++) oneda[j] = vector<double>(1);
    static vector<double> atry(ma), beta(ma), da(ma);
    if (alamda < 0.0)
    {
        mfit = 0;
        for (j = 0; j < ma; j++) if (ia[j]) mfit++;
        alamda = 0.001;
        mrqcof_HW2F(flag, df, t, tau, fM, x, y, sig, a, ia, alpha, beta, chisq, funcs);
        ochisq = chisq;
        for (j = 0; j < ma; j++) atry[j] = a[j];
    }
    vector<vector<double>> temp(mfit);
    for (j = 0; j < ma; j++) temp[j] = vector<double>(mfit);
    for (j = 0; j < mfit; j++)
    {
        for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
        covar[j][j] = alpha[j][j] * (1.0 + alamda);
        for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k];
        oneda[j][0] = beta[j];
    }
    gaussj(temp, oneda, singularflag);
    if (singularflag != true)
    {
        for (j = 0; j < mfit; j++)
        {
            for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k];
            da[j] = oneda[j][0];
        }
        if (alamda == 0.0)
        {
            covsrt(covar, ia, mfit);
            covsrt(alpha, ia, mfit);
            return;
        }
        for (j = 0, l = 0; l < ma; l++) if (ia[l]) atry[l] = a[l] + da[j++];
        mrqcof_HW2F(flag, df, t, tau, fM, x, y, sig, atry, ia, alpha, beta, chisq, funcs);
        if (chisq < ochisq)
        {
            alamda *= 0.1;
            ochisq = chisq;
            for (j = 0; j < mfit; j++)
            {
                for (k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
                beta[j] = da[j];
            }
            for (l = 0; l < ma; l++) a[l] = atry[l];
        }
        else
        {
            alamda *= 10.0;
            chisq = ochisq;
        }
    }
}

/*
void ybar_HW2F(vector<double> &c, vector<double> &AtT, vector<double> &BatT, vector<double> &BbtT, vector<double> &x, vector<double> &ybar)
{
    double temperror=1.0e-12;
    int n=int(c.size()), m=int(x.size()), i, j;
    double left0=-100.0, left, right0=100.0, right, ans, sum=0.0, ans0, shift=2.0;
    double sumr=0.0, suml=0.0;
    for(j=0;j<m;j++)
    {
        left0=-100.0, right0=100.0, sum=0.0, sumr=0.0, suml=0.0;
        for(i=0;i<n;i++) sumr=sumr+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*right0);
        for(i=0;i<n;i++) suml=suml+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*left0);
        While (sumr > 1#)
        {
            right0=right0*shift;
            sumr=0.0;
            for(i=0;i<n;i++) sumr=sumr+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*right0);
        }
        While (suml < 1#)
        {
            left0=left0*shift;
            suml=0.0;
            for(i=0;i<n;i++) suml=suml+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*left0);
        }

        left=left0;
        right=right0;
        ybar[j]=0.5*(left0+right0);
        for(i=0;i<n;i++) sum=sum+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*ybar[j]);
        ans=sum-1.0;
        ans0=ans-1.0;
        while(fabs(ans)>Error && ybar[j]>left0+Error && ybar[j]<right0-Error && ans0!=ans)
        {
            ans0=ans;
            if(ans>0.0)
            {
                left=ybar[j];
                ybar[j]=0.5*(left+right);
                sum=0.0;
                for(i=0;i<n;i++) sum=sum+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*ybar[j]);
                ans=sum-1.0;
            }
            else if(ans<0.0)
            {
                right=ybar[j];
                ybar[j]=0.5*(left+right);
                sum=0.0;
                for(i=0;i<n;i++) sum=sum+c[i]*AtT[i]*exp(-BatT[i]*x[j]-BbtT[i]*ybar[j]);
                ans=sum-1.0;
            }
            else break;
        }
    }
}
*/

void ybar_HW2F(vector<vector<double>>& lambij, vector<double>& BbtT, vector<double>& ybar)
{
    double temperror = 1.0e-12;
    int n = int(BbtT.size()), i, j;
    double left0 = -100.0, left, right0 = 100.0, right, ans, sum = 0.0, ans0, shift = 2.0;
    double sumr = 0.0, suml = 0.0;
    for (j = 0; j < NumGLQ; j++)
    {
        left0 = -100.0, right0 = 100.0, sum = 0.0, sumr = 0.0, suml = 0.0;
        for (i = 0; i < n; i++) sumr = sumr + lambij[i][j] * exp(-BbtT[i] * right0);
        for (i = 0; i < n; i++) suml = suml + lambij[i][j] * exp(-BbtT[i] * left0);
        
        while(sumr > 1.0)
        {
            right0 = right0 * shift;
            sumr = 0.0;
            for (i = 0; i < n; i++) sumr = sumr + lambij[i][j] * exp(-BbtT[i] * right0);
        }
        
        while(suml < 1.0)
        {
            left0 = left0 * shift;
            suml = 0.0;
            for (i = 0; i < n; i++) suml = suml + lambij[i][j] * exp(-BbtT[i] * left0);
        }

        left = left0;
        right = right0;
        ybar[j] = 0.5 * (left0 + right0);
        for (i = 0; i < n; i++) sum = sum + lambij[i][j] * exp(-BbtT[i] * ybar[j]);
        ans = sum - 1.0;
        ans0 = ans - 1.0;
        while (fabs(ans) > Error&& ybar[j] > left0 + Error && ybar[j] < right0 - Error && ans0 != ans)
        {
            ans0 = ans;
            if (ans > 0.0)
            {
                left = ybar[j];
                ybar[j] = 0.5 * (left + right);
                sum = 0.0;
                for (i = 0; i < n; i++) sum = sum + lambij[i][j] * exp(-BbtT[i] * ybar[j]);
                ans = sum - 1.0;
            }
            else if (ans < 0.0)
            {
                right = ybar[j];
                ybar[j] = 0.5 * (left + right);
                sum = 0.0;
                for (i = 0; i < n; i++) sum = sum + lambij[i][j] * exp(-BbtT[i] * ybar[j]);
                ans = sum - 1.0;
            }
            else break;
        }
    }
}

double PlainIROption_HW2F_SE(vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, vector<double>& param)
{
    double se = 0.0;
    int k, ndata = int(flag.size());
    for (k = 0; k < ndata; k++)
    {
        double err = Error * (1.0e+10);//1.0e-5;
        if (flag[k] == "CAP")
        {
            int i, ncaplet = int(tau[k].size());
            double expatT, exp2at, expbtT, exp2bt, expapbt, voltT, lnPtTX, pid1, pid2;
            y[k] = 0.0;
            double a = min(max(pow(param[0], 2.0), err), 1.0 / err), sig = min(max(pow(param[1], 2.0), err), 1.0 / err), b = min(max(pow(param[2], 2.0), err), 1.0 / err), eta = min(max(pow(param[3], 2.0), err), 1.0 / err), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])), -1.0 + err), 1.0 - err);
            for (i = 1; i < ncaplet; i++)
            {

                expatT = min(max(exp(-a * (t[k][i + 1] - t[k][i])), err), 1.0 / err);
                exp2at = min(max(exp(-2.0 * a * t[k][i]), err), 1.0 / err);
                expbtT = min(max(exp(-b * (t[k][i + 1] - t[k][i])), err), 1.0 / err);
                exp2bt = min(max(exp(-2.0 * b * t[k][i]), err), 1.0 / err);
                expapbt = min(max(exp(-(a + b) * t[k][i]), err), 1.0 / err);

                voltT = sqrt(sig * sig / (2.0 * a * a * a) * (1.0 - expatT) * (1.0 - expatT) * (1.0 - exp2at) + eta * eta / (2.0 * b * b * b) * (1.0 - expbtT) * (1.0 - expbtT) * (1.0 - exp2bt) + 2.0 * rho * sig * eta / (a * b * (a + b)) * (1.0 - expatT) * (1.0 - expbtT) * (1.0 - expapbt));

                lnPtTX = -log(df[k][i + 1] * (1.0 + x[k] * tau[k][i]) / df[k][i]);

                y[k] = y[k] + df[k][i] * cumnormal_Pol_App(lnPtTX / voltT + 0.5 * voltT) - df[k][i + 1] * (1.0 + x[k] * tau[k][i]) * cumnormal_Pol_App(lnPtTX / voltT - 0.5 * voltT);

                pid1 = df[k][i] * max(pdfnormal(lnPtTX / voltT + 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) + 0.5);
                pid2 = df[k][i + 1] * (1.0 + x[k] * tau[k][i]) * max(pdfnormal(lnPtTX / voltT - 0.5 * voltT), Error) * (-lnPtTX / (voltT * voltT) - 0.5);
            }
            se = se + (blp[k] - y[k]) * (blp[k] - y[k]);
        }
        else
        {
            double maxaaa = 10.0,minaaa = 0.0001;
            double maxd = 10000.0, mind = -10000.0;
            int i, j, nswaplet = int(tau[k].size());
            double expat, exp2at, expbt, exp2bt, sqrtexp2at, sqrtexp2bt, expapbt, siga, etab, rhoab, rhoab1, mua, mub, Vt, yy, denom, tmpexpBbtTybar;
            y[k] = 0.0;
            double a = min(max(pow(param[0],2.0),minaaa),maxaaa), sig = min(max(pow(param[1],2.0),minaaa),maxaaa), b = min(max(pow(param[2],2.0),minaaa),maxaaa), eta = min(max(pow(param[3],2.0),minaaa),maxaaa), rho = min(max(1.0 - 2.0 / (1.0 + exp(param[4])),-0.99),0.99);
            vector<double> expaT(nswaplet,0.0), exp2aT(nswaplet,0.0), expatT(nswaplet,0.0), exp2atT(nswaplet,0.0), expbT(nswaplet,0.0), exp2bT(nswaplet,0.0), expbtT(nswaplet,0.0), exp2btT(nswaplet,0.0), expapbT(nswaplet,0.0), expapbtT(nswaplet,0.0), tT(nswaplet,0.0), PMtT(nswaplet,0.0), VT(nswaplet,0.0), VtT(nswaplet,0.0), AtT(nswaplet,0.0), BatT(nswaplet,0.0), c(nswaplet,0.0), BbtT(nswaplet,0.0), ybar(NumGLQ,0.0), xx(NumGLQ,0.0), h1(NumGLQ,0.0);
            vector<vector<double>> lambij(nswaplet), expkappij(nswaplet), h2(nswaplet);

            expat = exp(-a * t[k][0]);
            exp2at = expat * expat;
            sqrtexp2at = sqrt((1.0 - exp2at) / (2.0 * a));
            expbt = exp(-b * t[k][0]);
            exp2bt = expbt * expbt;
            sqrtexp2bt = sqrt((1.0 - exp2bt) / (2.0 * b));
            expapbt = exp(-(a + b) * t[k][0]);
            siga = sig * sqrtexp2at;
            etab = eta * sqrtexp2bt;
            rhoab = rho * (1.0 - expapbt) / (sqrtexp2at * sqrtexp2bt) / (a + b);
            rhoab1 = sqrt(1.0 - rhoab * rhoab);

            mua = -(sig * sig / (a * a) + rho * sig * eta / (a * b)) * (1.0 - expat) + 0.5 * sig * sig / (a * a) * (1.0 - exp2at) + rho * sig * eta / (b * (a + b)) * (1.0 - expapbt);
            mub = -(eta * eta / (b * b) + rho * sig * eta / (a * b)) * (1.0 - expbt) + 0.5 * eta * eta / (b * b) * (1.0 - exp2bt) + rho * sig * eta / (a * (a + b)) * (1.0 - expapbt);

            Vt = sig * sig / (a * a) * (t[k][0] + 2.0 / a * expat - 0.5 / a * exp2at - 1.5 / a) + eta * eta / (b * b) * (t[k][0] + 2.0 / b * expbt - 0.5 / b * exp2bt - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[k][0] + (expat - 1.0) / a + (expbt - 1.0) / b - (expapbt - 1.0) / (a + b));

            for (j = 0; j < NumGLQ; j++)
            {
                xx[j] = M_SQRT2 * GLQx[j] * siga + mua;
            }

            for (i = 0; i < nswaplet; i++)
            {
                expaT[i] = exp(-a * t[k][i + 1]);
                expbT[i] = exp(-b * t[k][i + 1]);
                exp2aT[i] = expaT[i] * expaT[i];
                exp2bT[i] = expbT[i] * expbT[i];
                expapbT[i] = exp(-(a + b) * t[k][i + 1]);

                tT[i] = t[k][i + 1] - t[k][0];

                expatT[i] = exp(-a * tT[i]);
                expbtT[i] = exp(-b * tT[i]);
                exp2atT[i] = expatT[i] * expatT[i];
                exp2btT[i] = expbtT[i] * expbtT[i];
                expapbtT[i] = exp(-(a + b) * tT[i]);

                PMtT[i] = df[k][i + 1] / df[k][0];

                VT[i] = sig * sig / (a * a) * (t[k][i + 1] + 2.0 / a * expaT[i] - 0.5 / a * exp2aT[i] - 1.5 / a) + eta * eta / (b * b) * (t[k][i + 1] + 2.0 / b * expbT[i] - 0.5 / b * exp2bT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (t[k][i + 1] + (expaT[i] - 1.0) / a + (expbT[i] - 1.0) / b - (expapbT[i] - 1.0) / (a + b));
                VtT[i] = sig * sig / (a * a) * (tT[i] + 2.0 / a * expatT[i] - 0.5 / a * exp2atT[i] - 1.5 / a) + eta * eta / (b * b) * (tT[i] + 2.0 / b * expbtT[i] - 0.5 / b * exp2btT[i] - 1.5 / b) + 2.0 * rho * sig * eta / (a * b) * (tT[i] + (expatT[i] - 1.0) / a + (expbtT[i] - 1.0) / b - (expapbtT[i] - 1.0) / (a + b));

                c[i] = x[k] * tau[k][i];

                BatT[i] = (1.0 - expatT[i]) / a;
                BbtT[i] = (1.0 - expbtT[i]) / b;
                AtT[i] = min(max(PMtT[i] * exp(0.5 * (VtT[i] - VT[i] + Vt)),minaaa),maxaaa);

                lambij[i] = vector<double>(NumGLQ);
                expkappij[i] = vector<double>(NumGLQ);
                h2[i] = vector<double>(NumGLQ);

                for (j = 0; j < NumGLQ; j++)
                {
                    lambij[i][j] = c[i] * AtT[i] * exp(-BatT[i] * xx[j]);
                    expkappij[i][j] = exp(-BbtT[i] * (mub - 0.5 * rhoab1 * rhoab1 * etab * etab * BbtT[i] + rhoab * etab * GLQx[j] * M_SQRT2));
                }
            }
            for (j = 0; j < NumGLQ; j++)
            {
                lambij[nswaplet - 1][j] = lambij[nswaplet - 1][j] + AtT[nswaplet - 1] * exp(-BatT[nswaplet - 1] * xx[j]);
            }
            ybar_HW2F(lambij,BbtT,ybar);
            for (j = 0; j < NumGLQ; j++)
            {
                h1[j] = (ybar[j] - mub) / (etab * rhoab1) - rhoab * GLQx[j] * M_SQRT2 / rhoab1;

                denom = 0.0;

                for (i = 0; i < nswaplet; i++)
                {
                    tmpexpBbtTybar = exp(-BbtT[i] * ybar[j]);
                    denom = denom + lambij[i][j] * tmpexpBbtTybar * BbtT[i];
                }


                yy = 0.0;
                for (i = 0; i < nswaplet; i++)
                {
                    h2[i][j] = h1[j] + BbtT[i] * etab * rhoab1;
                    yy = yy - lambij[i][j] * expkappij[i][j] * cumnormal_Pol_App(-h2[i][j]);
                }
                yy = yy + cumnormal_Pol_App(-h1[j]);
                y[k] = y[k] + GLQw[j] * yy;
            }
            y[k] = y[k] * df[k][0] / M_SQRTPI;
            se = se + (blp[k] - y[k]) * (blp[k] - y[k]);
        }
    }
    return se;
}

double amotsa_HW2F(int& idum, double& tt, vector<vector<double>>& p, vector<double>& yy, vector<double>& psum, vector<double>& pb, double& yb, vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, double funcs(vector<string>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&), const int ihi, double& yhi, const double fac)
{
    int j;
    double fac1, fac2, yflu, ytry;

    int ndim = int(p[0].size());
    vector<double> ptry(ndim);
    fac1 = (1.0 - fac) / double(ndim);
    fac2 = fac1 - fac;
    for (j = 0; j < ndim; j++) ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    ytry = funcs(flag, df, t, tau, x, y, blp, ptry);
    ////////
    //cout<<ytry<<"\t"<<yb<<endl;
    ////////
    if (ytry <= yb)
    {
        for (j = 0; j < ndim; j++) pb[j] = ptry[j];
        yb = ytry;
    }
    yflu = ytry - tt * log(ran1(idum));
    if (yflu < yhi)
    {
        yy[ihi] = ytry;
        yhi = yflu;
        for (j = 0; j < ndim; j++)
        {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    return yflu;
}

void amebsa_HW2F(int& idum, double& tt, vector<vector<double>>& p, vector<double>& yy, vector<double>& pb, double& yb, const double ftol, vector<string>& flag, vector<vector<double>>& df, vector<vector<double>>& t, vector<vector<double>>& tau, vector<double>& x, vector<double>& y, vector<double>& blp, double funcs(vector<string>&, vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&), int& iter, const double temptr)
{
    /*///////////////////////////////////////////////////
        ofstream foutpara("inputparamaa.dat");
        foutpara<<"idum="<<"\t"<<idum<<endl;
        foutpara<<"tt="<<"\t"<<tt<<endl;

        int nprow=int(p.size()), npcol=int(p[0].size());
        foutpara<<"p:"<<"\t"<<nprow<<"x"<<npcol<<endl;
        for(int ii=0;ii<nprow;ii++)
        {
            for(int jj=0;jj<npcol;jj++)
            {
                foutpara<<"p["<<ii<<"]["<<jj<<"]="<<p[ii][jj]<<"\t";
            }
            foutpara<<endl;
        }

        int nyy=int(yy.size());
        foutpara<<"yy:"<<"\t"<<nyy<<endl;
        for(int ii=0;ii<nyy;ii++) foutpara<<"yy["<<ii<<"]="<<yy[ii]<<"\t";
        foutpara<<endl;

        int npb=int(pb.size());
        foutpara<<"pb:"<<"\t"<<npb<<endl;
        for(int ii=0;ii<npb;ii++) foutpara<<"pb["<<ii<<"]="<<pb[ii]<<"\t";
        foutpara<<endl;

        foutpara<<"yb="<<"\t"<<yb<<endl;
        foutpara<<"ftol="<<"\t"<<ftol<<endl;

        int nflag=int(flag.size());
        foutpara<<"flag="<<"\t"<<nflag<<endl;
        for(int ii=0;ii<nflag;ii++) foutpara<<"flag["<<ii<<"]="<<flag[ii]<<"\t";
        foutpara<<endl;

        int ndfrow=int(df.size());
        foutpara<<"df:\t"<<"ndfrow="<<ndfrow<<endl;
        for(int ii=0;ii<ndfrow;ii++)
        {
            int ndfcol=int(df[ii].size());
            foutpara<<"df:\t"<<"ndfcol="<<ndfcol<<endl;
            for(int jj=0;jj<ndfcol;jj++)
            {
                foutpara<<"df["<<ii<<"]["<<jj<<"]="<<df[ii][jj]<<"\t";
            }
            foutpara<<endl;
        }

        int ntrow=int(t.size());
        foutpara<<"t:\t"<<"ntrow="<<ntrow<<endl;
        for(int ii=0;ii<ntrow;ii++)
        {
            int ntcol=int(t[ii].size());
            foutpara<<"t:\t"<<"ntcol="<<ntcol<<endl;
            for(int jj=0;jj<ntcol;jj++)
            {
                foutpara<<"t["<<ii<<"]["<<jj<<"]="<<t[ii][jj]<<"\t";
            }
            foutpara<<endl;
        }

        int ntaurow=int(tau.size());
        foutpara<<"tau:\t"<<"ntaurow="<<ntaurow<<endl;
        for(int ii=0;ii<ntaurow;ii++)
        {
            int ntaucol=int(tau[ii].size());
            foutpara<<"tau:\t"<<"ntaucol="<<ntaucol<<endl;
            for(int jj=0;jj<ntaucol;jj++)
            {
                foutpara<<"tau["<<ii<<"]["<<jj<<"]="<<tau[ii][jj]<<"\t";
            }
            foutpara<<endl;
        }

        int nx=int(x.size());
        foutpara<<"x:"<<"\t"<<nx<<endl;
        for(int ii=0;ii<nx;ii++) foutpara<<"x["<<ii<<"]="<<x[ii]<<"\t";
        foutpara<<endl;

        int ny=int(y.size());
        foutpara<<"y:"<<"\t"<<ny<<endl;
        for(int ii=0;ii<ny;ii++) foutpara<<"y["<<ii<<"]="<<y[ii]<<"\t";
        foutpara<<endl;

        int nblp=int(blp.size());
        foutpara<<"blp:"<<"\t"<<nblp<<endl;
        for(int ii=0;ii<nblp;ii++) foutpara<<"blp["<<ii<<"]="<<blp[ii]<<"\t";
        foutpara<<endl;

        foutpara<<"iter="<<"\t"<<iter<<endl;
        foutpara<<"temptr="<<"\t"<<temptr<<endl;

    *////////////////////////////////////////////////////
    int i, ihi, ilo, j, n;
    double rtol, yhi, ylo, ynhi, ysave, yt, ytry;

    int mpts = int(p.size());
    int ndim = int(p[0].size());
    vector<double> psum(ndim);
    /*/////////////////////////////////////////////////////
        int npsum=int(psum.size());
        foutpara<<"psum:"<<"\t"<<npsum<<endl;
        for(int ii=0;ii<npsum;ii++) foutpara<<"psum["<<ii<<"]="<<psum[ii]<<"\t";
        foutpara<<endl;
    */////////////////////////////////////////////////////
    tt = -temptr;
    get_psum(p, psum);
    for (;;)
    {
        ilo = 0;
        ihi = 1;
        ynhi = ylo = yy[0] + tt * log(ran1(idum));
        yhi = yy[1] + tt * log(ran1(idum));
        if (ylo > yhi)
        {
            ihi = 0;
            ilo = 1;
            ynhi = yhi;
            yhi = ylo;
            ylo = ynhi;
        }
        for (i = 3; i <= mpts; i++)
        {
            yt = yy[i - 1] + tt * log(ran1(idum));
            if (yt <= ylo)
            {
                ilo = i - 1;
                ylo = yt;
            }
            if (yt > yhi)
            {
                ynhi = yhi;
                ihi = i - 1;
                yhi = yt;
            }
            else if (yt > ynhi) ynhi = yt;
        }
        rtol = 2.0 * fabs(yhi - ylo) / (fabs(yhi) + fabs(ylo));
        if (rtol < ftol || iter < 0)
        {
            SWAP(yy[0], yy[ilo]);
            for (n = 0; n < ndim; n++) SWAP(p[0][n], p[ilo][n]);
            break;
        }
        iter -= 2;
        ytry = amotsa_HW2F(idum, tt, p, yy, psum, pb, yb, flag, df, t, tau, x, y, blp, funcs, ihi, yhi, -1.0);
        if (ytry <= ylo) ytry = amotsa_HW2F(idum, tt, p, yy, psum, pb, yb, flag, df, t, tau, x, y, blp, funcs, ihi, yhi, 2.0);
        else if (ytry <= ynhi)
        {
            ysave = yhi;
            ytry = amotsa_HW2F(idum, tt, p, yy, psum, pb, yb, flag, df, t, tau, x, y, blp, funcs, ihi, yhi, 0.5);
            if (ytry >= ysave)
            {
                for (i = 0; i < mpts; i++)
                {
                    if (i != ilo)
                    {
                        for (j = 0; j < ndim; j++)
                        {
                            psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                            p[i][j] = psum[j];
                        }
                        yy[i] = funcs(flag, df, t, tau, x, y, blp, psum);
                    }
                }
                iter -= ndim;
                get_psum(p, psum);
            }
        }
        else ++iter;
    }
}



