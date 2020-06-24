#include "Currency.h"

CCurrency::CCurrency()
{

}

CCurrency::CCurrency(string FULLNAME, string SYMBOL)
{
	fullname = FULLNAME;
	symbol = SYMBOL;
}

CCurrency::CCurrency(string SYMBOL)
{
	symbol = SYMBOL;
	fullname = SYMBOL;
}

string CCurrency::get_fullname()
{
	return fullname;
}

string CCurrency::get_symbol()
{
	return symbol;
}

CCurrency::~CCurrency()
{

}