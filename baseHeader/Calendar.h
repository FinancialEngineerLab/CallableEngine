//#pragma once

#ifndef  _CALENDAR_H
#define _CALENDAR_H

#include <string>
#include <vector>
#include "Date.h"
#include "ConstantNumbers.h"

using namespace std;

class CCalendar
{
private:
	string holidaysname;
	int numholidays;
	vector<CDate> holidays;

public:
	CCalendar();
	~CCalendar();
	CCalendar(string HOLIDAYSNAME);

	string get_holidaysname();
	int get_numholidays();
	CDate operator[](int i);

};

bool is_holiday(CDate date, CCalendar holidays);

void ShiftBusDate(CDate date, CCalendar holidays, int numshift, CDate& newdate);

void ConvReflDate(CDate date, CCalendar holidays, string conv, CDate& newdate);

#endif 


