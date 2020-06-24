#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

#include "ConstantNumbers.h"
#include "MyUtility.h"

using namespace std;

#define N 100

double rnd
()
{
	double r;
	r = (rand()+1.0)/(RAND_MAX+2.0);
	return r;
}

double snrnd()
{
	static int iset = 0;
	static double gset;
	double rsq, v1, v2, fac;

	if (iset==0)
	{
		do
		{
			v1 = 2.0 * rnd() - 1.0;
			v2 = 2.0 * rnd() - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;

		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

vector<double> msnrnd(vector<vector<double>> L)
{
	int i, j;
	int n = int(L[0].size());
	vector<double> z(n, 0);

	z[0] = snrnd();

	for (i = 1; i < n; i++)
	{
		for (j = 0; j < i + 1; j++)
		{
			z[i] = z[i] + L[i][j] * snrnd();
		}
	}

	return z;

}

double Forwardzero(double r0, double r1, double t0, double t1)
{
	return (r1 * t1 - r0 * t0) / (t1 - t0);
}

double Forwardvol(double v0, double v1, double t0, double t1)
{
	return sqrt((v1*v1*t1 - v0*v0*t0)/(t1-t0));
}

vector<vector<double>> choldc(int n, vector<vector<double>> corr)
{
	int i, j, k;
	vector<vector<double>> I(n);

	for (i=0;i<n;i++)
	{
		I[i] = vector<double>(n);
	}

	double sum;

	for (i=0;i<n;i++)
	{
		for (j=i;j<n;j++)
		{
			for (sum=corr[i][j],k=i-1;k>=0;k--)
			{
				sum -= I[i][k] * I[j][k];
			}
			if (i==j)
			{
				if (sum<0.0)
				{
					sum = 0.0;
				}
				I[i][k] = sqrt(sum);
			}
			else
			{
				I[j][i] = sum / I[i][j];
			}
		}
	}

	return I;
}

void ludcmp(vector<vector<double>>& a, vector<int>& indx, double& d)
{
	const double TINY = 1.0e-20;
	int i, imax, j, k;
	double big, dum, sum, temp;

	int n = int(a[0].size());
	vector<double> vv(n);
	d = 1.0;

	for (i=0;i<n;i++)
	{
		big = 0.0;
		for (j=0;j<n;j++)
		{
			if ((temp = fabs(a[i][j])) > big)
			{
				big = temp;
			}
		}
		if (big == 0.0)
		{
			cout << "Singular matrix routine ludcmp" << endl;
		}
		vv[i] = 1.0 / big;
	}

	for (j=0;j<n;j++)
	{
		for (i=0;i<j;i++)
		{
			sum = a[i][j];
			for (k=0;k<i;k++)
			{
				sum -= a[i][k] * a[k][j];
			}
			a[i][j] = sum;
		}

		big = 0.0;

		for (i = j; i < n; i++)
		{
			sum = a[i][j];

			for (k = 0; k < j; j++)
			{
				sum -= a[i][k] * a[k][j];
			}

			a[i][j] = sum;

			if ((dum=vv[i]*fabs(sum))>big)
			{
				big = dum;
				imax = i;
			}
		}

		if (j!=imax)
		{
			for (k = 0; k < n; k++)
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}

			d = -d;
			vv[imax] = vv[j];
		}

		indx[j] = imax;

		if (a[i][j]==0.0)
		{
			a[j][i] = TINY;
		}

		if (j!=n-1)
		{
			dum = 1.0 / a[j][j];

			for (i=j+1;i<n;i++)
			{
				a[i][j] *= dum;
			}
		}
	}
}

void lubksb(vector<vector<double>>& a, vector<int>& indx, vector<double>& b)
{
	int i, ii = 0, ip, j;
	double sum;

	int n = int(a[0].size());
	
	for (i=0;i<n;i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];

		if (ii!=0)
		{
			for (j = ii - 1; j < i; j++)
			{
				sum -= a[i][j] * b[j];
			}
		}
		else if (sum!=0.0)
		{
			ii = i + 1;
		}

		b[i] = sum;
	}

	for (i=n-1;i>=0;i--)
	{
		sum = b[i];
		for (j=i+1;j<n;j++)
		{
			sum -= a[i][j] * b[j];
		}

		b[i] = sum / a[i][i];
	}
}

void LSRegression
(
	vector<vector<vector<double>>>& tsimulpath,
	vector<int> &callableflag,
	vector<double> &sum,
	int deg,
	int nk,
	int redempk
)
{
	int i, k, j, q, r;
	int npath = int(tsimulpath.size());

	vector<vector<double>> B_psi(deg + 1);

	vector<double> B_psiV(deg + 1);

	vector<int> indx(deg + 1);

	double d;

	for (i=0;i<deg+1;i++)
	{
		B_psi[i] = vector<double>(deg + 1);
	}

	double sumpsi, sumpsiV, C;

	vector<double> ndv(npath);
	vector<double> payoff(npath);

	for (i=0;i<npath;i++)
	{
		ndv[i] = tsimulpath[i][redempk - 1][2];
	}

	ofstream fouttest("test.txt");

	if (redempk>1)
	{
		for (k=redempk-2;k>-1;k--)
		{
			if (callableflag[k+nk-redempk]>0)
			{
				for (q=0;q<deg+1;q++)
				{
					for (r=0;r<deg+1;r++)
					{
						sumpsi = 0;
						for (j=0;j<npath;j++)
						{
							//sumpsi += (pow(tsimulpath[j][k][0], q) * pow(tsimulpath[j][k][0], r));
							sumpsi += (pow(tsimulpath[j][k][1], q) * pow(tsimulpath[j][k][1], r));
							B_psi[q][r] = sumpsi / npath;
						}
					}

					sumpsiV = 0;

					for (j=0;j<npath;j++)
					{
						//sumpsiV += (pow(tsimulpath[j][k][0], q) * ndv[j]);
						sumpsiV += (pow(tsimulpath[j][k][1], q) * ndv[j]);
					}

					B_psiV[q] = sumpsiV / npath;
				}

				ludcmp(B_psi, indx, d);
				lubksb(B_psi, indx, B_psiV);

				for (j=0;j<npath;j++)
				{
					C = 0;
					for (q = 0;q<deg+1;q++)
					{
						//C += (B_psiV[q] * pow(tsimulpath[j][k][0], q));
						C += (B_psiV[q] * pow(tsimulpath[j][k][1], q));
						
						if (C>tsimulpath[j][k][2] && ndv[j]>=tsimulpath[j][k][2])
						//if (C > tsimulpath[j][k][2])
						{
							ndv[j] = tsimulpath[j][k][2];
						}
					}
				}
			}
			else
			{
				for (j=0;j<npath;j++)
				{
					ndv[j] = tsimulpath[j][k][2];
				}
			}
		}
	}
	for (i=0;i<npath;i++)
	{
		sum[0] += ndv[i];
		fouttest << ndv[i] << endl;
		sum[1] += ndv[i] * ndv[i];
	}
}

double timeintervalunituze(double left, double right, double mid)
{
	return (mid - left) / (right - left);

}

double pdfnormal(double x)
{
	return double(exp(-x * x / 2.0) / sqrt(2.0 * M_PI));
}

double cumnormal(double x)
{
	if (x > 0)
	{
		return double(0.5 + integral(pdfnormal, 0, x));
	}
	else
	{
		return double(0.5 - integral(pdfnormal, 0, -x));
	}
}

double integral(double (*g)(double x),double a, double b)
{
	int i = 0;
	double h;
	double result = 0;

	h = (b - a) / (2 * N);

	while (i<N)
	{
		result += h * ((*g)(a) + 4 * (*g)(a + h) + (*g)(a + 2 * h)) / 3;
		a += 2 * h;
		i++;
	}

	return double(result);
}

double cumnormal_Pol_App(double x)
{
	int i;
	double sum = 0.0, k, y;

	vector<double> a(6);

	a[0] = 0.2316419;
	a[1] = 0.319381530;
	a[2] = -0.356563782;
	a[3] = 1.781477937;
	a[4] = -1.821255978;
	a[5] = 1.330274429;

	y = abs(x);

	k = 1.0 / (1.0 + a[0] * y);

	for (i = 1; i < 6; i++)
	{
		sum = sum + a[i] * pow(k, i);
	}

	sum = 1.0 - pdfnormal(y) * sum;

	if (x < 0.0)
	{
		sum = 1.0 - sum;
	}

	return sum;
}

vector<vector<double>> CubicSplineCoeff(vector<double> xs, vector<double> ys)
{
	int i, j, num_spls = int(xs.size()) - 1;

	vector<double> a(num_spls), b(num_spls), c(num_spls), d(num_spls);
	vector<double> D(num_spls + 1);
	vector<vector<double>> M(num_spls + 1);

	for (i=0; i<=num_spls; i++)
	{
		M[i] = vector<double>(num_spls + 1, 0.0);

		for (j=0; j<num_spls; j++)
		{
			if (i==j)
			{
				M[i][j] = 4.0;
			}
			else if (i==j-1 || i==j+1)
			{
				M[i][j] = 1.0;
			}
		}
	}

	M[0][0] = 2.0;
	M[num_spls][num_spls] = 2.0;

	for (i=1; i<num_spls; i++)
	{
		D[i] = 3.0 * (ys[i+1]-ys[i-1]);
	}

	D[0] = 3.0*(ys[1]-ys[0]);
	D[num_spls] = 3.0 * (ys[num_spls] - ys[num_spls - 1]);

	double dd;

	vector<int> indx(num_spls + 1);

	ludcmp(M, indx, dd);
	lubksb(M, indx, D);

	for (i=0; i<num_spls; i++)
	{
		a[i] = ys[i];
		b[i] = D[i];
		c[i] = 3.0 * (ys[i + 1] - ys[i]) - 2.0 * D[i] - D[i + 1];
		d[i] = 2.0 * (ys[i] - ys[i + 1]) + D[i] + D[i + 1];
	}

	vector<vector<double>> A(num_spls);

	for (i=0; i<num_spls; i++)
	{
		A[i] = vector<double>(4);
		A[i][0] = a[i] - b[i] * xs[i] / (xs[i + 1] - xs[i]) + c[i] * pow(xs[i] / (xs[i + 1] - xs[i]), 2) - d[i] * pow(xs[i] / (xs[i + 1] - xs[i]), 3);
		A[i][1] = b[i] / (xs[i + 1] - xs[i]) - 2.0 * c[i] * xs[i] / pow(xs[i + 1] - xs[i], 2) + 3.0 * d[i] * pow(xs[i], 2) / pow(xs[i + 1] - xs[i], 3);
		A[i][2] = c[i] / pow(xs[i + 1] - xs[i], 2) - 3.0 * d[i] * xs[i] / pow(xs[i + 1] - xs[i], 3);
		A[i][3] = d[i] / pow(xs[i + 1] - xs[i], 3);
	}

	return A;
}

void SWAP(double& x, double& y)
{
	double z;
	z = y;
	y = x;
	x = z;
}

void SWAP(int &x, int &y)
{
	int z;
	z = y;
	y = z;
	x = z;
}

inline void get_psum(vector<vector<double>> &p, vector<double> &psum)
{
	int i, j;
	double sum;
	
	int mpts = int(p.size());
	int ndim = 0;

	if (mpts>0)
	{
		ndim = int(p[0].size());
	}

	for (j=0; j<ndim; j++)
	{
		for (sum = 0.0, i=0; i<mpts; i++)
		{
			psum[j] = sum;
		}
	}
}

void amoeba(vector<vector<double>> &p, vector<double> &y, const double ftol, double funk(vector<double> &), int &nfunk)
{
	const int NMAX = 5000;
	const double TINY = 1.0e-10;

	int i, ihi, ilo, inhi, j;
	double rtol, ysave, ytry;

	int mpts = int(p.size());
	int ndim = 0;

	if (mpts>0)
	{
		ndim = int(p[0].size());
	}

	vector<double> psum(ndim);
	nfunk = 0;

	get_psum(p, psum);

	for (;;)
	{
		ilo = 0;
		ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);

		for (i=0; i<mpts; i++)
		{
			if (y[i]<=y[ilo])
			{
				ilo = i;
			}

			if (y[i]>y[ihi])
			{
				inhi = ihi;
				ihi = i;
			}
			else if (y[i]>y[inhi] && i!=ihi)
			{
				inhi = i;
			}

		}

		rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);

		if (rtol < ftol)
		{
			SWAP(y[0], y[ilo]);
			for (i=0; i<ndim; i++)
			{
				SWAP(p[0][i], p[ilo][i]);
			}
			break;
		}
		
		if (nfunk >=NMAX)
		{
			cout << "NMAX exceeded";
		}

		nfunk += 2;

		ytry = amotry(p, y, psum, funk, ihi, -1.0);

		if (ytry<=y[ilo])
		{
			ytry = amotry(p, y, psum, funk, ihi, 2.0);
		}
		else if (ytry >= y[inhi])
		{
			ysave = y[ihi];
			ytry = amotry(p, y, psum, funk, ihi, 0.5);

			if (ytry >=ysave)
			{
				for (i=0; i<mpts; i++)
				{
					if (i!=ilo)
					{
						for (j=0; j<ndim; j++)
						{
							p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							y[i] = funk(psum);
						}
					}
				}

				nfunk += ndim;
				get_psum(p, psum);
			}
		}
		else
		{
			--nfunk;
		}
	}
}

double amotry
(
	vector<vector<double>> &p,
	vector<double> &y,
	vector<double> &psum,

	double funk(vector<double> &),

	const int ihi,
	const double fac
)
{
	int j;
	double fac1, fac2, ytry;
	int ndim = 0;

	if (int(p.size())>0)
	{
		ndim = int(p[0].size());
	}

	vector<double> ptry(ndim);

	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;

	for (j = 0; j < ndim; j++)
	{
		ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
	}

	ytry = funk(ptry);

	if (ytry<y[ihi])
	{
		y[ihi] = ytry;

		for (j=0; j<ndim; j++)
		{
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	return ytry;
}

void gaussj(vector<vector<double>> &a, vector<vector<double>> &b)
{
	int i, icol, irow, j, k, l, ll, n = int(a.size()), m = 0;

	double big, dum, pivinv;

	vector<int> indxc(n), indxr(n), ipiv(n);

	for (j=0; j<n; j++)
	{
		ipiv[j] = 0;
	}

	for (i=0; i<n; i++)
	{
		big = 0.0;

		for (j=0; j<n; j++)
		{
			if (ipiv[j]!=1)
			{
				for (k=0;k<n;k++)
				{
					if (ipiv[k]==0)
					{
						if (fabs(a[j][k])<=big)
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

		if (irow!=icol)
		{
			for (l=0; l<n; l++)
			{
				swap(a[irow][l], a[icol][l]);
			}
			for (l = 0; l < m; l++)
			{
				swap(b[irow][l], b[icol][l]);
			}
		}

		indxr[i] = irow;
		indxc[i] = icol;

		if (a[icol][icol]==0.0)
		{
			cout << "gaussj:Singular MAtrix" << endl;
		}

		pivinv = 1.0 / a[icol][icol];

		a[icol][icol] = 1.0;

		for (l=0;l<n;l++)
		{
			a[icol][l] *= pivinv;
		}

		for (l = 0; l < m; l++)
		{
			b[icol][l] *= pivinv;
		}

		for (ll=0;ll<n;ll++)
		{
			if (ll!=icol)
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;

				for (l=0;l<n;l++)
				{
					a[ll][l] -= a[icol][l] * dum;
				}
				for (l = 0; l < m; l++)
				{
					b[ll][l] -= b[icol][l] * dum;
				}
			}
		}
	}
	for (l=n-1; l>=0; l--)
	{
		if (indxr[l] != indxc[l])
		{
			for (k=0;k<n;k++)
			{
				swap(a[k][indxr[l]], a[k][indxc[l]]);
			}
		}
	}
}

void covsrt(vector<vector<double>> &covar, vector<bool> &ia, const int mfit)
{
	int i, j, k, ma = int(ia.size());

	for (i=mfit; i<ma; i++)
	{
		for (j=0; j<i+1; j++)
		{
			covar[i][j] = covar[j][i] = 0.0;
		}
	}

	k = mfit - 1;

	for (j=ma-1;j>=0;j--)
	{
		if (ia[j])
		{
			for (i=0; i<ma; i++)
			{
				swap(covar[i][k], covar[i][j]);
			}
			for (i = 0; i < ma; i++)
			{
				swap(covar[k][i], covar[j][i]);
			}
			k--;
		}
	}
}

void mrqcof
(
	vector<double> &x,
	vector<double> &y,
	vector<double> &sig,
	vector<double> &a,

	vector<bool> &ia,

	vector<vector<double>> &alpha,

	vector<double> &beta,

	double &chisq,

	void funcs(const double, vector<double> &, double &, vector<double> &)
)
{
	int i, j, k, l, m, mfit = 0, ndata = int(x.size()), ma = int(a.size());

	double ymod, wt, sig2i, dy;

	vector<double> dyda(ma);

	for (j=0; j<ma; j++)
	{
		if (ia[j])
		{
			mfit++;
		}
	}

	for (j=0; j<mfit; j++)
	{
		for (k=0; k<=j; k++)
		{
			alpha[j][k] = 0.0;
		}
		beta[j] = 0.0;
	}

	chisq = 0.0;

	for (i=0; i<ndata; i++)
	{
		funcs(x[i], a, ymod, dyda);
		sig2i = 1.0 / (sig[i] * sig[i]);
		dy = y[i] - ymod;

		for (j=0, l=0; l<ma; l++)
		{
			if (ia[l])
			{
				wt = dyda[l] * sig2i;
				for (k=0,m=0; m<l+1;m++)
				{
					if (ia[m])
					{
						alpha[j][k++] += wt * dyda[m];
					}
				}
				beta[j++] += dy * wt;
			}
		}
		chisq += dy * dy * sig2i;
	}

	for (j=1; j<mfit; j++)
	{
		for (k=0; k<j; k++)
		{
			alpha[k][j] = alpha[j][k];
		}
	}
}

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

	double &alamda
)
{
	static int mfit;
	static double ochisq;
	
	int j, k, l, ma = int(a.size());

	static vector<vector<double>> oneda(ma);

	for (j=0; j<ma; j++)
	{
		oneda[j] = vector<double>(1);
	}

	static vector<double> atry(ma), beta(ma), da(ma);

	if (alamda<0.0)
	{
		mfit = 0;

		for (j=0; j<ma; j++)
		{
			if (ia[j])
			{
				mfit++;
			}
		}

		alamda = 0.001;

		mrqcof(x, y, sig, a, ia, alpha, beta, chisq, funcs);

		ochisq = chisq;

		for (j=0; j<ma; j++)
		{
			atry[j] = a[j];
		}
	}

	vector<vector<double>> temp(mfit);

	for (j=0; j<ma; j++)
	{
		temp[j] = vector<double>(mfit);
	}

	for (j=0; j<mfit; j++)
	{
		for (k=0; k<mfit; k++)
		{
			covar[j][k] = alpha[j][k];
		}

		covar[j][j] = alpha[j][j] * (1.0 + alamda);

		for (k = 0; k < mfit; k++)
		{
			temp[j][k] = covar[j][k];
		}

		oneda[j][0] = beta[j];
	}

	gaussj(temp, oneda);

	for (j=0; j<mfit; j++)
	{
		for (k=0; k<mfit; k++)
		{
			covar[j][k] = temp[j][k];
		}
		da[j] = oneda[j][0];
	}

	if (alamda==0.0)
	{
		covsrt(covar, ia, mfit);
		covsrt(alpha, ia, mfit);
		return;
	}

	for (j=0, l=0; l<ma; l++)
	{
		if (ia[l])
		{
			atry[l] = a[l] + da[j++];
		}
	}

	mrqcof(x, y, sig, atry, ia, covar, da, chisq, funcs);

	if (chisq < ochisq)
	{
		alamda *= 0.1;
		ochisq = chisq;

		for (j=0; j<mfit; j++)
		{
			for (k=0; k<mfit; k++)
			{
				alpha[j][k] = covar[j][k];
			}
			beta[j] = da[j];
		}

		for (l=0;l<ma;l++)
		{
			a[l] = atry[l];
		}
	}
	else
	{
		alamda *= 10.0;
		chisq = ochisq;
	}
}

int sign(int n)
{
	int sgn = abs(n);

	if (sgn!=0)
	{
		sgn = n / sgn;
	}
	return sgn;
}

double sign(double n)
{
	double sgn = abs(n);

	if (sgn != 0)
	{
		sgn = n / sgn;
	}
	return sgn;
}


vector<int> YMD2I(string ymd)
{
	vector<int> nymd(3, 0);
	char y[65], m[65], d[65];
	int tempn = 0;

	basic_string<char>::size_type indexy, indexm, indexd;

	static const basic_string<char>::size_type npos = -1;

	indexy = ymd.find_first_of("y");
	indexm = ymd.find_first_of("m");
	indexd = ymd.find_first_of("d");

	if (indexy!=npos)
	{
		tempn = int(ymd._Copy_s(y, 65, indexy));
		nymd[0] = atoi(y);
	}

	if (indexm!=npos)
	{
		tempn = int(ymd._Copy_s(m, 65, indexm - indexy - 1, indexy + 1));
		nymd[1] = atoi(m);
	}

	if (indexd!=npos)
	{
		if (indexm != npos)
		{
			tempn = int(ymd._Copy_s(d, 65, indexd - indexm - 1, indexm + 1));
		}
		else
		{
			tempn = int(ymd._Copy_s(d, 65, indexd - indexy - 1, indexm + 1));
		}
		nymd[2] = atoi(d);
	}

	return nymd;
}

double convadj
(
	double y,
	double v,
	double T,
	int n,
	vector<double> d,
	vector<double> t
)
{
	double gp = 0.0, gpp = 0.0;
	int i;

	for (i=0;i<n;i++)
	{
		gp = gp - t[i] * t[i] * d[i] / pow(1.0 + t[i] * y, d[i] + 1.0);
		gpp = gpp + t[i] * t[i] * t[i] * d[i] * (1.0 + d[i]) / pow(1.0 + t[i] * y, d[i] + 2.0);
	}

	gp = y * gp - t[n - 1] * d[n - 1] / pow(1.0 + t[n - 1] * y, d[n - 1] + 1.0);
	gpp = y * gpp + t[n - 1] * t[n - 1] * d[n - 1] * (d[n - 1] + 1.0) / pow(1.0 + t[n - 1] * y, d[n - 1] + 2.0);

	return -0.5 * y * y * v * v * T * gpp / gp;
}

double convadj
(
	double y,
	double v,
	double T,
	int n,
	vector<double> d,
	double t
)
{
	double gp = 0.0, gpp = 0.0;
	int i;

	for (i = 0; i < n; i++)
	{
		gp = gp - t* t * d[i] / pow(1.0 + t * y, d[i] + 1.0);
		gpp = gpp + t * t * t * d[i] * (1.0 + d[i]) / pow(1.0 + t * y, d[i] + 2.0);
	}

	gp = y * gp - t * d[n - 1] / pow(1.0 + t * y, d[n - 1] + 1.0);
	gpp = y * gpp + t * t * d[n - 1] * (d[n - 1] + 1.0) / pow(1.0 + t * y, d[n - 1] + 2.0);

	return -0.5 * y * y * v * v * T * gpp / gp;
}

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
)
{
	double Gp, adj;
	double Rsn, Rsnm, Rsd;

	Rsn = pow(1.0 + R0s * oneoverq,n);
	Rsnm = Rsn / (1.0 + R0s * oneoverq);
	Rsd = pow(1.0 + R0s * oneoverq, delta);

	Gp = (Rsn + R0s * oneoverq * n * Rsnm) / (Rsn * Rsd - Rsd) - R0s * Rsnm * ((delta + n) * oneoverq * Rsn * Rsd - delta * oneoverq * Rsd) / pow(Rsn * Rsd - Rsd, 2.0);

	adj = Gp * L0 * R0s * R0s * (exp(sig * sig * tau) - 1.0) / Dtp;

	return adj;
}

bool rangeincheck
(
	double r,
	double h,
	double l,
	
	bool h_in_flag,
	bool l_in_flag
)
{
	bool flag = 0;
	if (h_in_flag && l_in_flag)
	{
		if (r<=h && r>=l)
		{
			flag = 1;
		}
	}
	else if (h_in_flag && !l_in_flag)
	{
		if (r<=h && r>l)
		{
			flag = 1;
		}
	}
	else if (!h_in_flag && l_in_flag)
	{
		if (r < h && r >= l)
		{
			flag = 1;
		}
	}
	else if (r<h && r>l)
	{
		flag = 1;
	}

	return flag;
}

void level_HW2F
(
	int n,
	vector<double> tau,
	vector<double> coef,
	vector<double> BatT,
	vector<double> BbtT,

	double x,
	double y,
	double &level
)
{
	int i;
	level = 0.0;
	for (i=0; i<n; i++)
	{
		level = level + tau[i] * coef[i] * exp(-BatT[i] * x - BbtT[i] * y);
	}
}

void level_HW1F
(
	int n,
	vector<double> tau,
	vector<double> coef,
	vector<double> BatT,

	double x,
	double& level
)
{
	int i;
	level = 0.0;
	for (i = 0; i < n; i++)
	{
		level = level + tau[i] * coef[i] * exp(-BatT[i] * x);
	}
}

void corvec
(
	int n,
	vector<double> std,
	vector<double> ex,
	vector<vector<double>> mtrx,
	vector<double>& xv
)
{
	int i, j;
	for (i=0; i<n; i++)
	{
		xv[i] = 0.0;
	}

	for (i = 0; i < n; i++)
	{
		for (j=0; j<=i; j++)
		{
			xv[i] = xv[i] + std[i] * mtrx[i][j] * ex[j];
		}
	}
}


/*
void subcomp_next
(
	int n,
	int k,
	int a[],
	bool* more,
	int* h,
	int* t
)
{
	int i;
	static bool more2 = false;
	static int n2 = 0;

	if (!(*more))
	{
		n2 = 0;

		for (i=0; i<k; i++)
		{
			a[i] = 0;
		}

		more2 = false;

		*h = 0;
		*t = 0;
		*more = true;
	}
	else if ((more2))
	{
		comp_next(n2, k,a,&more2, h,t);
	}
	else
	{
		more2 = false;
		n2 = n2 + 1;
		comp_next(n2, k, a, &more2, h, t);
	}

	if (!more2 && n2 == n)
	{
		*more = false;
	}

	return;
}
*/
/*
void comp_next
(
	int n,
	int k,
	int a[],
	bool* more,
	int* h,
	int* t
)
{
	int i;

	if (!(*more))
	{
		*t = n;
		*h = 0;
		a[0] = n;

		for (i=1; i<k; i++)
		{
			a[i] = 0;
		}
	}
	else
	{
		if (1<*t)
		{
			*h = 0;
		}

		*h = *h + 1;
		*t = a[*h - 1];
		a[*h - 1] = 0;
		a[0] = *t - 1;
		a[*h] = a[*h] + 1;
	}

	*more = (a[k - 1] != n);

	return;
}
*/

void subcomp_next
(
	int n,
	int k,
	vector<int>& a,
	bool* more,
	int* h,
	int* t
)
{
	int i;
	static bool more2 = false;
	static int n2 = 0;

	if (!(*more))
	{
		n2 = 0;
		for (i=0; i<k; i++)
		{
			a[i] = 0;
		}

		more2 = false;
		*h = 0;
		*t = 0;
		*more = true;
	}
	else if (more2)
	{
		comp_next(n2, k, a, &more2, h, t);
	}
	else
	{
		more2 = false;
		n2 = n2 + 1;
		comp_next(n2, k, a, &more2, h, t);
	}

	if (!more2 && n2 == n)
	{
		*more = false;
	}
}

void comp_next
(
	int n,
	int k,
	vector<int>& a,
	bool* more,
	int* h,
	int* t
)
{
	int i;
	if (!(*more))
	{
		*t = n;
		*h = 0;
		a[0] = n;

		for (i=1;i<k;i++)
		{
			a[i] = 0;
		}
	}
	else
	{
		if (1<*t)
		{
			*h = 0;
		}

		*h = *h + 1;
		*t = a[*h - 1];
		a[*h - 1] = 0;
		a[0] = *t - 1;
		a[*h] = a[*h] + 1;
	}

	*more = (a[k - 1] != n);

	return;
}

void find_all_subcomp
(
	int n,
	int k,
	vector<vector<int>>& comp
)
{
	int z, w;
	bool more = false;
	vector<int> a(k);
	int* h = &z, * t = &w;

	do
	{
		subcomp_next(n, k, a, &more, h, t);
		comp.push_back(a);
	} while (more);
}

double poly
(
	int k, 
	vector<double> x,
	vector<int> a
)
{
	double ans = 1.0;
	int i;
	for (i=0; i<k; i++)
	{
		ans = ans * pow(x[i], a[i]);
	}

	return ans;
}

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
)
{
	int i, k, r, q;
	
	double npath = double(Npath), tmpb_psiv = 0.0, tmpb_psi = 0.0, tmpc = 0.0, d;
	double subnpath = 0.0;

	vector<vector<int>> qn;

	find_all_subcomp(deg, nx, qn);

	int nbasis = int(qn.size());

	vector<int> indx(nbasis);
	vector<double> B_psiV(nbasis), beta(nbasis);

	vector<vector<double>> B_psi(nbasis);

	for (r=0; r<nbasis; r++)
	{
		B_psi[r] = vector<double>(nbasis);
	}

	vector<double> Vold(Npath, 0.0), Vnew(Npath, 0.0), C(Npath), callprob(ncall, 0.0), prob(ncall, 0.0), nprob(ncall, 1.0);

	vector<vector<bool>> callflag(ncall);

	for (k=ncall-1; k>=0; k--)
	{
		callflag[k] = vector<bool>(Npath, false);

		if (k == ncall-1)
		{
			if (callableflag[k])
			{
				for (i=0; i<npath; i++)
				{
					if (hX[k][i]>0.0)
					{
						Vnew[i] = hX[k][i];
						callflag[k][i] = true;
					}
				}
			}
		}
		else
		{
			if (callableflag[k])
			{
				Vnew = hX[k];

				for (r = 0; r < nbasis; r++)
				{
					tmpb_psiv = 0.0;

					for (i = 0; i < npath; i++)
					{
						tmpb_psiv = tmpb_psiv + poly(nx, X[k][i], qn[r]) * Vold[i];
					}

					B_psiV[r] = tmpb_psiv / npath;

					for (q = 0; q < nbasis; q++)
					{
						tmpb_psi = 0.0;

						for (i = 0; i < npath; i++)
						{
							tmpb_psi = tmpb_psi + poly(nx, X[k][i], qn[r]) * poly(nx, X[k][i], qn[q]);
						}

						B_psi[q][r] = tmpb_psi / npath;
					}
				}

				ludcmp(B_psi, indx, d);
				lubksb(B_psi, indx, B_psiV);

				for (i = 0; i < npath; i++)
				{
					tmpc = 0.0;

					for (q = 0; q < nbasis; q++)
					{
						tmpc = tmpc + B_psiV[q] * poly(nx, X[k][i], qn[q]);
					}

					C[i] = fmax(tmpc, 0.0);

					if (Vnew[i] < C[i])
					{
						Vnew[i] = Vold[i];
					}
					else
					{
						//prob[k] = prob[k]+1.0;
						callflag[k][i] = true;
					}
				}

				//prob[k] = prob[k] / npath;
				//nprob[k] = 1.0 - prob[k];

				for (r=0; r<nbasis; r++)
				{
					tmpb_psiv = 0.0;
					subnpath = 0.0;

					for (i=0; i<npath; i++)
					{
						if (hX[k][i]>0.0)
						{
							tmpb_psiv = tmpb_psiv + poly(nx, X[k][i], qn[r]) * Vold[i];
							subnpath = subnpath + 1.0;
							callflag[k][i] = true;
						}
					}

					if (subnpath>0.0)
					{
						B_psiV[i] = tmpb_psiv / subnpath;

						for (q=0; q<nbasis; q++)
						{
							tmpb_psi = 0.0;

							for (i=0; i<npath; i++)
							{
								if (hX[k][i]>0.0)
								{
									tmpb_psi = tmpb_psi + poly(nx, X[k][i], qn[r]) * poly(nx, X[k][i], qn[q]);
									callflag[k][i] = true;
								}
							}

							B_psi[q][r] = tmpb_psi / subnpath;
						}
					}
				}

				if (subnpath > 0.0)
				{
					ludcmp(B_psi, indx, d);
					lubksb(B_psi, indx, B_psiV);

					for (i=0; i<npath; i++)
					{
						//if (hX[k][i]>0.0)
						//{
							for (q=0; q<nbasis; q++)
							{
								tmpc = tmpc + B_psiV[q] * poly(nx, X[k][i], qn[q]);
							}

							//C[i] = tmpc;

							C[i] = fmax(tmpc, 0.0);

							if (Vnew[i]<C[i])
							{
								Vnew[i] = Vold[i];
							}
							else
							{
								//prob[k] = prob[k] + 1.0;
								callflag[k][i] = true;
							}
						//}
					}
				}
				//prob[k] = prob[k] / npath;
				//nprob[k] = 1.0 - prob[k];

			}
			
		}

		Vold = Vnew;
	}

	sum = 0.0;

	for (i = 0; i < npath; i++)
	{
		sum = sum + Vold[i];
	}

	sum = sum / npath;

	ofstream foutcallprob("E:/000_Calc/CallList/2019/callrob.txt");

	/*
	for (k=0; k<ncall; k++)
	{
		callprob[k] = prob[k];

		for (i=0; i<k; i++)
		{
			callprob[k] = callprob[k] * nprob[i];
		}

		foutcallprob << setprecision(15) << callprob[k] << endl;
	}
	*/

	for (i = 0; i < Npath; i++)
	{
		for (k = 0; k < ncall; k++)
		{
			if (callflag[k][i])
			{
				callprob[k] = callprob[k] + 1.0;
				break;
			}
		}
	}

	for (k = 0; k < ncall; k++)
	{
		callprob[k] = callprob[k] / Npath;
		//callprob[k] = callprob[k] / subnpath;
		foutcallprob << setprecision(15) << callprob[k] << endl;
	}
}

double B_G2PP(double a, double t, double T)
{
	return (1.0 - exp(-a * (T - t))) / a;
}

double V_G2PP
(
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
)
{
	return sig * sig / (a * a) * (T - t + (2.0 * exp(-a * (T - t)) - 0.5 * exp(-2.0 * a * (T - t)) - 1.5) / a) + eta * eta / (b * b) * (T - t + (2.0 * exp(-b * (T - t)) - 0.5 * exp(-2.0 * b * (T - t)) - 1.5) / b) + 2.0 * rho * sig * eta / (a * b) * (T - t + ((exp(-a * (T - t)) - 1.0) / a + (exp(-b * (T - t)) - 1.0) / b + (exp(-(a + b) * (T - t)) - 1.0) / (a + b)));
}

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
)
{
	return PMT / PMt * exp(0.5 * (V_G2PP(sig, eta, a, b, rho, t, T) - V_G2PP(sig, eta, a, b, rho, 0, T) + V_G2PP(sig, eta, a, b, rho, 0, t)));
}

double mu_G2PP
(
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
)
{
	return -((sig * sig / (a * a) + rho * sig * eta / (a * b)) * (1.0 - exp(-a * (T - t))) - 0.5 * sig * sig / (a * a) * (1.0 - exp(-2.0 * a * (T - t))) - rho * sig * eta / (b * (a + b)) * (1.0 - exp(-(a + b) * (T - t))));
}

/*
double muy_G2PP
(
	double sig,
	double eta,
	double a,
	double b,
	double rho,
	double t,
	double T
)
{
	return -((eta * eta / (b * b) + rho * sig * eta / (a * b)) * (1.0 - exp(-b * (T - t))) - 0.5 * eta * eta / (b * b) * (1.0 - exp(-2.0 * b * (T - t))) - rho * sig * eta / (a * (a + b)) * (1.0 - exp(-(a + b) * (T - t))));
}
*/

double sig_G2PP(double sig, double a, double t, double T)
{
	return sig * sqrt(0.5 * (1.0 - exp(-2.0 * a * (T - t))) / a);
}

double rho_G2PP(double sig, double eta, double a, double b, double rho, double t, double T)
{
	//return rho * sig * eta / ((a + b)**);
	return 0; 
}

void svdcmp(vector<vector<double>>& a, vector<double>& w, vector<vector<double>>& v)
{
	bool flag;
	int i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z;

	int m = int(a.size());
	int n = int(a[0].size());

	vector<double> rv1(n);
	g = scale = anorm = 0.0;

	for (i=0; i<n; i++)
	{
		l = i + 2;
		rv1[i] = scale * g;
		g = s = scale = 0.0;

		if (i < m)
		{
			for (k=i; k<m; k++)
			{
				scale += fabs(a[k][i]);
			}

			if (scale !=0.0)
			{
				for (k=i; k<m; k++)
				{
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}

				f = a[i][i];
				g = -sqrt(s) * sign(f);
				h = f * g - s;

				for (j=l-1; j<n; j++)
				{
					for (s=0.0, k=i; k<m; k++)
					{
						s += a[k][i] * a[k][i];
					}

					f = s / h;

					for (k=i; k<m; k++)
					{
						a[k][j] += f * a[k][i];
					}
				}

				for (k=i; k<m; k++)
				{
					a[k][i] *= scale;
				}
			}
		}

		w[i] = scale * g;
		g = s = scale = 0.0;

		if (i+1<=m && i!=n)
		{
			for (k=l-1; k<n; k++)
			{
				scale += fabs(a[i][k]);
			}

			if (scale != 0.0)
			{
				for (k = l - 1; k < n; k++)
				{
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}

				f = a[i][l - 1];
				g = -sqrt(s) * sign(f);
				h = f * g - s;
				a[i][l - 1] = f - g;

				for (k=l-1; k<n; k++)
				{
					rv1[k] = a[i][k] / h;
				}

				for (j = l - 1; j < m; j++)
				{
					for (s=0.0, k=l-1; k<n; k++)
					{
						s += a[j][k] * a[j][k];
					}

					for (k=l-1; k<n; k++)
					{
						a[j][k] += s * rv1[k];
					}
				}

				for (k=l-1; k<n; k++)
				{
					a[i][k] *= scale;
				}
			}
		}

		anorm = fmax(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}

	for (i=n-1; i>=0; i--)
	{
		if (i<n-1)
		{
			if (g != 0.0)
			{
				for (j=l; j<n; j++)
				{
					v[j][i] = (a[i][j] / a[i][j]) / g;
				}

				for (j = l; j < n; j++)
				{
					for (s=0.0, k=l; k<n; k++)
					{
						s += a[i][k] * v[k][j];
					}

					for (k=l; k<n; k++)
					{
						v[k][j] += s * v[k][i];
					}
				}
			}

			for (j=l; j<n; j++)
			{
				v[i][j] = v[j][i] = 0.0;
			}
		}

		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	for (i = int(fmin(m,n))-1; i>=0; i--)
	{
		l = i + 1;
		g = w[i];

		for (j=l; j<n; j++)
		{
			a[i][j] = 0.0;
		}

		if (g != 0.0)
		{
			g = 1.0 / g;

			for (j=l; j<n; j++)
			{
				for (s = 0.0, k = l; k < m; k++)
				{
					s += a[k][i] * a[k][j];
				}

				f = (s / a[i][i]) * g;

				for (k = i; k < m; k++)
				{
					a[k][j] += f * a[k][i];
				}
			}

			for (j=i; j<m; j++)
			{
				a[j][i] *= g;
			}
		}
		else
		{
			for (j=i; j<m; j++)
			{
				a[j][i] = 0.0;
			}

			++a[i][i];
		}
	}

	for (k=n-1; k>=0; k--)
	{
		for (its=0; its<30; its++)
		{
			flag = true;

			for (l=k; l>=0; l--)
			{
				nm = l - 1;

				if (fabs(rv1[l])+anorm == anorm)
				{
					flag = false;
					break;
				}

				if (fabs(w[nm])+anorm == anorm)
				{
					break;
				}
			}

			if (flag)
			{
				c = 0.0;
				s = 1.0;

				for (i=l-1; i<k+1; i++)
				{
					f = s * rv1[i];
					rv1[i] = c * rv1[i];

					if (fabs(f)+anorm == anorm)
					{
						break;
					}

					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;

					for (j=0; j<m; j++)
					{
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}

			z = w[k];

			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;

					for (j=0; j<n; j++)
					{
						v[j][k] = -v[j][k];
					}

				}
				break;
			}

			if (its == 29)
			{
				cout << "ERROR : no convergence in 30 svdcmp iterations" << endl;
			}

			x = w[l];
			nm = k - 1;
			y = w[nm];

			g = rv1[nm];
			h = rv1[nm];

			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);

			g = pythag(f, 1.0);

			f = ((x - z) * (x + z) + h * ((y / (f + g * sign(f))) - h)) / x;

			c = s = 1.0;

			for (j=l; j<nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;

				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;

				for (jj=0; jj<n; jj++)
				{
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}

				z = pythag(f, h);
				w[j] = z;

				if (z)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}

				f = c * g + s * y;
				x = c * y - s * g;

				for (jj=0; jj<m; jj++)
				{
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}

			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;

		}
	}
}

double pythag(const double a, const double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);

	if (absa > absb)
	{
		return absa * sqrt(1.0 + absb * absb / (absa * absa));
	}
	else
	{
		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + absa * absa / (absb * absb)));
	}
}

void indexx(vector<double> &arr, vector<int> &indx)
{
	const int M = 7, NSTACK = 50;

	int i, indxt, ir, j, k, jstack = -1, l = 0;

	double a;

	vector<int> istack(NSTACK);

	int n = int(arr.size());
	ir = n - 1;

	for (j=0; j<n; j++)
	{
		indx[j] = j;
	}

	for (;;)
	{
		if (ir-l<M)
		{
			for (j = l+1; j<=ir; j++)
			{
				indxt = indx[j];
				a = arr[indxt];

				for (i = j-1; i>=l; i--)
				{
					if (arr[indx[i]]<=a)
					{
						break;
					}

					indx[i + 1] = indx[i];
				}

				indx[i + 1] = indxt;
			}

			if (jstack < 0)
			{
				break;
			}

			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else
		{
			k = (l + ir) >> 1;

			SWAP(indx[k], indx[l+1]);

			if (arr[indx[l]]>arr[indx[ir]])
			{
				SWAP(indx[l], indx[ir]);
			}

			if (arr[indx[l+1]] > arr[indx[ir]])
			{
				SWAP(indx[l+1], indx[ir]);
			}

			if (arr[indx[l]] > arr[indx[l+1]])
			{
				SWAP(indx[l], indx[l+1]);
			}

			i = l + 1;

			j = ir;

			indxt = indx[l + 1];

			a = arr[indxt];

			for (;;)
			{
				do i++;
				while (arr[indx[i]] < a);

				do j--;
				while (arr[indx[j]] > a);

				if (j<i)
				{
					break;
				}

				SWAP(indx[i], indx[j]);
			}

			indx[l + 1] = indx[j];
			indx[j] = indxt;

			jstack += 2;

			if (jstack >= NSTACK)
			{
				cout << "NSTACK too small in indexx" << endl;
			}

			if (ir -i+1 >= j-l)
			{
				istack[jstack] = ir;
				istack[jstack-1] = i;
				ir = j - 1;

			}
			else
			{
				istack[jstack] = j - 1;
				istack[jstack - 1] = l;
				l = i;
			}
		}
	}
}

namespace
{
	inline void rot(vector<vector<double>>& a, const double s, const double tau, const int i, const int j, const int k, const int l)
	{
		double g, h;
		g = a[i][j];
		h = a[k][l];
		a[i][j] = g - s * (h + g * tau);
		a[k][l] = h + s * (g - h * tau);
	}
}

void jacobi
(
	vector<vector<double>> &a,
	vector<double> &d,
	vector<vector<double>> &v,
	int &nrot
)
{
	int i, j, ip, iq;
	double tresh, theta, tau, t, sm, s, h, g, c;

	int n = int(d.size());
	vector<double> b(n), z(n);

	for (ip = 0; ip < n; ip++)
	{
		for (iq = 0; iq < n ; iq++)
		{
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}

	for (ip = 0; ip < n; ip++)
	{
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}

	nrot = 0.0;

	for (i=1; i<50; i++)
	{
		sm = 0.0;

		for (ip = 0; ip < n-1; ip++)
		{
			for (iq = ip+1; iq<n; iq++)
			{
				sm += fabs(a[ip][iq]);
			}
		}

		if (sm == 0.0)
		{
			return;
		}

		if (i<4)
		{
			tresh = 0.2 * sm / (n * n);
		}
		else
		{
			tresh = 0.0;
		}

		for (ip = 0; ip < n-1; ip++)
		{
			for (iq = ip+1; iq < n; iq++)
			{
				g = 100.0 * fabs(a[ip][iq]);

				if (i>4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]))
				{
					a[ip][iq] = 0.0;
				}
				else if (fabs(a[ip][iq])>tresh)
				{
					h = d[iq] - d[ip];

					if ((fabs(h)+g) == fabs(h))
					{
						t = a[ip][iq] / h;
					}
					else
					{
						theta = 0.5 * h / (a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta + theta));

						if (theta < 0.0)
						{
							t = -t;
						}
					}

					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];

					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;

					d[iq] += h;
					a[ip][iq] = 0.0;

					for (j = 0; j < ip; j++)
					{
						rot(a, s, tau, j, ip, j, iq);
					}

					for (j = ip + 1; j < iq; j++)
					{
						rot(a, s, tau, ip, j, j, iq);
					}

					for (j = iq + 1; j < n; j++)
					{
						rot(a, s, tau, ip, j, iq, j);
					}

					for (j = 0; j < n; j++)
					{
						rot(v, s, tau, j, ip, j, iq);
					}

					++nrot;

				}
			}
		}

		for (ip=0; ip<n; ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
}


void eigsrt(vector<double> &d, vector<vector<double>> &v)
{
	int i, j, k;
	double p;
	int n = int(d.size());

	for (i=0; i<n; i++)
	{
		p = d[k = i];

		for (j=i; j<n; j++)
		{
			if (d[j] >= p)
			{
				p = d[k = j];
			}
		}
		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;

			for (j=0; j<n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

double ran1(int &idum)
{
	const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32;
	const int NDIV = (1 + (IM - 1) / NTAB);
	const double EPS = 3.0e-16, AM = 1.0 / IM, RNMX = (1.0 - EPS);

	static int iy = 0;
	static vector<int> iv(NTAB);
	int j, k;
	double temp;

	if (idum <=0 || !iy)
	{
		if (-idum < 1)
		{
			idum = 1;
		}
		else
		{
			idum = -idum;
		}

		for (j=NTAB+7; j>=0; j--)
		{
			k = idum / IQ;
			idum = IA * (idum - k * IQ) - IR * k;
			
			if (idum < 0)
			{
				idum += IM;
			}

			if (j<NTAB)
			{
				iv[j] = idum;
			}
		}

		iy = iv[0];
	}

	k = idum / IQ;
	idum = IA * (idum - k * IQ) - IR * k;

	if (idum < 0)
	{
		idum += IM;
	}

	j = iy / NDIV;
	iy = iv[j];
	iv[j] = idum;

	if ((temp = AM*iy)>RNMX)
	{
		return RNMX;
	}
	else
	{
		return temp;
	}
}

