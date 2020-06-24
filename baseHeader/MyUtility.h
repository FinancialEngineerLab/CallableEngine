#pragma once
#include <vector>
#include <string>
#include <stdlib.h>
#include <iomanip>

using namespace std;

double rnd();

double	snrnd();

vector<double> msnrnd(vector<vector<double>> corr);

double Forwardzero(double r0, double r1, double t0, double t1);

double Forwardvol(double v0, double v1, double t0, double t1);

vector<vector<double>> choldc(int n, vector<vector<double>> corr);

void ludcmp(vector<vector<double>>& a, vector<int>& indx, double& d);

void lubksb(vector<vector<double>>& a, vector<int>& indx, vector<double>& b);

void LSRegression
(
	vector<vector<vector<double>>>& tsimulpath,
	vector<int>& callableflag,
	vector<double>& sum,
	int deg,
	int nk,
	int redempk
);

double timeintervalunituze(double left, double right, double mid);

double pdfnormal(double x);

double cumnormal(double x);

double cumnormal_Pol_App(double x);

double integral(double (*g)(double x), double a, double b);

vector<vector<double>> CubicSplineCoeff(vector<double> xs, vector<double> ys);

void SWAP(double& x, double& y);

void SWAP(int& x, int& y);

inline void get_psum(vector<vector<double>>& p, vector<double>& psum);

void amoeba(vector<vector<double>>& p, vector<double>& y, const double ftol, double funk(vector<double>&), int& nfunk);

double amotry
(
	vector<vector<double>>& p,
	vector<double>& y,
	vector<double>& psum,

	double funk(vector<double>&),

	const int ihi,
	const double fac
);

void gaussj(vector<vector<double>>& a, vector<vector<double>>& b);

void covsrt(vector<vector<double>>& covar, vector<bool>& ia, const int mfit);

void mrqcof
(
	vector<double>& x,
	vector<double>& y,
	vector<double>& sig,
	vector<double>& a,

	vector<bool>& ia,

	vector<vector<double>>& alpha,

	vector<double>& beta,

	double& chisq,

	void funcs(const double, vector<double>&, double&, vector<double>&)
);

void mrqmin
(
	vector<double>& x,
	vector<double>& y,
	vector<double>& sig,
	vector<double>& a,

	vector<bool>& ia,

	vector<vector<double>>& covar,

	vector<vector<double>>& alpha,

	double& chisq,

	void funcs(const double, vector<double>&, double&, vector<double>&),

	double& alamda
);

int sign(int n);

double sign(double n);

vector<int> YMD2I(string ymd);

double convadj
(
	double y,
	double v,
	double T,
	int n,
	vector<double> d,
	vector<double> t
);

double convadj
(
	double y,
	double v,
	double T,
	int n,
	vector<double> d,
	double t
);

double convadj_strd
(
	double R0s,
	double sig,
	double delta,
	double tau,
	double oneoverq,
	int n,
	double L0,
	double Dtp
);

bool rangeincheck
(
	double r,
	double h,
	double l,

	bool h_in_flag,
	bool l_in_flag
);

void corvec
(
	int n,
	vector<double> std,
	vector<double> ex,
	vector<vector<double>> mtrx,
	vector<double>& xv
);

void level_HW2F
(
	int n,
	vector<double> tau,
	vector<double> coef,
	vector<double> BatT,
	vector<double> BbtT,

	double x,
	double y,
	double& level
);

void level_HW1F
(
	int n,
	vector<double> tau,
	vector<double> coef,
	vector<double> BatT,

	double x,
	double& level
);

void subcomp_next
(
	int n,
	int k,
	vector<int>& a,
	bool* more,
	int* h,
	int* t
);

void comp_next
(
	int n,
	int k,
	vector<int>& a,
	bool* more,
	int* h,
	int* t
);

void find_all_subcomp
(
	int n,
	int k,
	vector<vector<int>>& comp
);

double poly
(
	int k,
	vector<double> x,
	vector<int> a
);

void LSRegression
(
	vector<vector<vector<double>>>& X,
	vector<vector<double>>& hX,
	vector<bool>& callableflag,
	int ncall,
	int Npath,
	int deg,
	int nx,
	double& sum
);

double B_G2PP(double a, double t, double T);

double V_G2PP
(
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
);

double A_G2PP
(
	double PMt,
	double PMT,
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
);

double mu_G2PP
(
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
);


double sig_G2PP(double sig, double a, double t, double T);

double rho_G2PP(double sig, double eta, double a, double b, double rho, double t, double T);

void svdcmp(vector<vector<double>>& a, vector<double>& w, vector<vector<double>>& v);

double pythag(const double a, const double b);

void indexx(vector<double>& arr, vector<int>& indx);

void jacobi
(
	vector<vector<double>>& a,
	vector<double>& d,
	vector<vector<double>>& v,
	int& nrot
);

void eigsrt(vector<double>& d, vector<vector<double>>& v);

double ran1(int& idum);
