#pragma once

#ifndef _DATE_H_
#define _DATE_H_

#include <string>
#include <vector>
#include <stdlib.h>
#include "MyUtility.h"

using namespace std;;

class CDate
{
private:
	int y;
	int m;
	int d;
	int totaldays;

	double t;

	string yyyymmdd;

	vector<int> ymd;
	
public:
	CDate();
	~CDate();
	CDate(int Y, int M, int D);
	CDate(string yyyymmdd);
	CDate(vector<int> YMD);
	CDate(int totaldays);

	string get_str_date();
	vector<int> get_vec_date();

	int get_int_date();
	int get_year();
	int get_month();
	int get_day();

	void set_year(int year);
	void set_month(int month);
	void set_day(int day);

	double get_dbl_date();

	int operator[](int i);

	double operator-(const CDate & date) const;

	bool operator==(const CDate & date);
	bool operator!=(const CDate & date);
	bool operator<(const CDate & date);
	bool operator>(const CDate & date);
	bool operator>=(const CDate & date);
	bool operator<=(const CDate & date);

};

double cvg(CDate t_st, CDate t_end, string dcb);

int get_monthdays(int year, int month);

int get_totaldays(int year, int month, int day);

int is_leapyear(int year);

vector<int> totaldays2YYYYMMDD(int totaldays);

vector<int> get_yyyymmdd(string yyyymmdd);

string set_yyyymmdd(vector<int>YYYYMMDD);

int get_weekday(CDate date);

CDate DateAdd(char ymd, int numymd, CDate date);

int get_month_from_MMMYY(string str);

int get_year_from_MMMYY(string str);

#endif // !_DATE_H_




