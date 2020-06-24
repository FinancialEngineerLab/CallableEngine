#pragma once

#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <vector>
#include "Date.h"
#include <algorithm>
#include <typeinfo>

using namespace std;

class CInterpolation
{
private:

	vector<double> xs;
	vector<double> ys;
	vector<CDate> xds;

	vector<double> x1s;
	vector<double> x2s;
	vector<CDate> xd1s;
	vector<CDate> xd2s;
	vector<vector<double>> zs;

public:

	CInterpolation();
	~CInterpolation();

	CInterpolation(vector<double> Xs, vector<double> Ys);
	CInterpolation(vector<CDate> Xs, vector<double> Ys);
	double operator()(double x);
	double operator()(CDate xd);

	CInterpolation(vector<double> X1s, vector<double> X2s, vector<vector<double>> Zs);
	CInterpolation(vector<CDate> X1s, vector<double> X2s, vector<vector<double>> Zs);
	CInterpolation(vector<CDate> X1s, vector<CDate> X2s, vector<vector<double>> Zs);

	double operator()(double x1, double x2);
	double operator()(CDate xd1, CDate xd2);
	double operator()(CDate xd1, double x2);

	int get_xs_size();
	int get_xds_size();
	int get_x1s_size();
	int get_xd1s_size();
	int get_x2s_size();
	int get_xd2s_size();

	vector<double> get_xs();
	vector<double> get_ys();
	vector<CDate> get_xds();
	vector<double> get_x1s();
	vector<CDate> get_xd1s();
	vector<double> get_x2s();
	vector<CDate> get_xd2s();
	vector<vector<double>> get_zs();

	CInterpolation operator+(CInterpolation& zc);

};

#endif