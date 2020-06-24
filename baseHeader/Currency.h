#pragma once

#ifndef _CURRENCY_H_
#define _CURRENCY_H_

#include <string>
#include <stdlib.h>

using namespace std;

class CCurrency
{
private:
	string fullname;
	string symbol;

public:
	CCurrency();
	~CCurrency();
	CCurrency(string FULLNAME, string SYMBOL);
	CCurrency(string SYMBOL);
	string get_fullname();
	string get_symbol();

};


#endif // !_CURRENCY_H_
