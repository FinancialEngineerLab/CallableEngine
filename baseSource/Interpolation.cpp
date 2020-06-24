# include "Interpolation.h"

CInterpolation::CInterpolation()
{
}

CInterpolation::CInterpolation(vector<double> Xs, vector<double> Ys)
{
	xs = Xs;
	ys = Ys;
}

CInterpolation::CInterpolation(vector<CDate> Xs, vector<double> Ys)
{
	xds = Xs;
	ys = Ys;
}

CInterpolation::CInterpolation(vector<double> X1s, vector<double> X2s, vector<vector<double>> Zs)
{
	x1s = X1s;
	x2s = X2s;
	zs = Zs;
}

CInterpolation::CInterpolation(vector<CDate> Xd1s, vector<CDate> Xd2s, vector<vector<double>> Zs)
{
	xd1s = Xd1s;
	xd2s = Xd2s;
	zs = Zs;
}

CInterpolation::CInterpolation(vector<CDate> Xd1s, vector<double> X2s, vector<vector<double>> Zs)
{
	xd1s = Xd1s;
	x2s = X2s;
	zs = Zs;
}

double CInterpolation::operator()(double x)
{
	double y = ys[0];
	int i, n = int(xs.size());

	if (x<=xs[0])
	{
		y = ys[0];
	}
	else if( x>= xs[n-1])
	{
		y = ys[n - 1];
	}
	else
	{
		for (i=1; i<n; i++)
		{
			if (x==xs[i-1])
			{
				y = ys[i - 1];
			}
			else if (x > xs[i-1] && x < xs[i])
			{
				y = ((ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1])) * (x - xs[i]) + ys[i];
			}
		}
	}

	return y;
}

double CInterpolation::operator()(CDate xd)
{
	double y = ys[0];

	int i, n = int(xds.size());

	if (xd<=xds[0])
	{
		y = ys[0];
	}
	else if (xd >= xds[n-1])
	{
		y = ys[n - 1];
	}
	else
	{
		for (i=1; i<n; i++)
		{
			if (xd==xds[i-1])
			{
				y = ys[i - 1];
			}
			else if(xd > xds[i-1] && xd < xds[i])
			{
				y = ((ys[i] - ys[i - 1]) / (xds[i] - xds[i - 1])) * (xd - xds[i]) + ys[i];
			}
		}
	}

	return y;
}

double CInterpolation::operator()(double x1, double x2)
{
	double z = zs[0][0];

	int i, j, x1n = int(x1s.size()), x2n = int(x2s.size());

	if (x1 <= x1s[0])
	{
		if (x2 <= x2s[0])
		{
			z = zs[0][0];
		}
		else if (x2 >= x2s[x2n-1])
		{
			z = zs[0][x2n - 1];
		}
		else
		{
			for (j=1; j<x2n; j++)
			{
				if (x2==x2s[j-1])
				{
					z = zs[0][j - 1];
					break;
				}
				else if (x2 > x2s[j-1] && x2 < x2s[j])
				{
					z = ((zs[0][j] - zs[0][j - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[0][j];
					break;
				}
			}
		}
	}
	else if (x1 >=x1s[x1n-1])
	{
		if (x2 <= x2s[0])
		{
			z = zs[x1n - 1][0];
		}
		else if (x2 >= x2s[x2n-1])
		{
			z = zs[x1n - 1][x2n - 1];
		}
		else
		{
			for (j=1; j<x2n; j++)
			{
				if (x2==x2s[j-1])
				{
					z = zs[x1n - 1][j - 1];
					break;
				}
				else if (x2>x2s[j-1] && x2<x2s[j])
				{
					z = ((zs[x1n - 1][j] - zs[x1n - 1][j - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[x1n - 1][j];
					break;
				}
			}
		}
	}
	else
	{
		if (x2 <= x2s[0])
		{
			for (i=1; i<x1n; i++)
			{
				if (x1==x1s[i-1])
				{
					z = zs[i - 1][0];
					break;
				}
				else if (x1>x1s[i-1] && x1<x1s[i])
				{
					z = ((zs[i][0] - zs[i - 1][0]) / (x1s[i] - x1s[i - 1])) * (x1 - x1s[i]) + zs[i][0];
					break;
				}
			}
		}
		else if (x2>=x2s[x2n-1])
		{
			for (i=1; i<x1n; i++)
			{
				if (x1==x1s[i-1])
				{
					z = zs[i - 1][x2n - 1];
					break;
				}
				else if (x1>x1s[i-1] && x1<x1s[i])
				{
					z = ((zs[i][x2n - 1] - zs[i - 1][x2n - 1]) / (x2s[i] - x1s[i - 1])) * (x1 - x1s[i]) + zs[i][x2n - 1];
					break;
				}
			}
		}
		else
		{
			for (i=1; i<x1n; i++)
			{
				if (x1==x1s[i-1])
				{
					for (j=1; j<x2n; j++)
					{
						if (x2 == x2s[j - 1])
						{
							z = zs[i - 1][j - 1];
							break;
						}
						else if (x2 > x2s[j - 1] && x2 < x2s[j])
						{
							z = ((zs[i - 1][j] - zs[i - 1][j - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[i - 1][j];
							break;
						}
					}
				}
				else if (x1>x1s[i-1] && x1<x1s[i])
				{
					double lambda1;
					lambda1 = (x1s[i] - x1) / (x1s[i] - x1s[i - 1]);

					for (j=1; j<x2n; j++)
					{
						if (x2 == x2s[j-1])
						{
							z = ((zs[i][j - 1] - zs[i - 1][j - 1]) / (x1s[i] - x1s[i - 1])) * (x1 - x1s[i]) + zs[i][j - 1];
							break;
						}
						else if (x2 > x2s[j-1] && x2 < x2s[j])
						{
							double lambda2;
							lambda2 = (x2s[j] - x2) / (x2s[j] - x2s[j - 1]);

							z = (1.0 - lambda1) * (1.0 - lambda2) * zs[i][j] + (1.0 - lambda1) * lambda2 * zs[i][j - 1] + lambda1 * (1.0 - lambda2) * zs[i - 1][j] + lambda1 * lambda2 * zs[i - 1][j - 1];
							break;
						}
					}
				}
			}
		}
	}

	return z;
}

double CInterpolation::operator()(CDate xd1, CDate xd2)
{
	double z = zs[0][0];
	int i, j, x1n = int(xd1s.size()), x2n = int(xd2s.size());

	if (xd1 <= xd1s[0])
	{
		if (xd2 <=xd2s[0])
		{
			z = zs[0][0];
		}
		else if (xd2 >= xd2s[x2n-1])
		{
			z = zs[0][x2n - 1];
		}
		else
		{
			for (j=1; j<x2n; j++)
			{
				if (xd2 == xd2s[j-1])
				{
					z = zs[0][j - 1];
					break;
				}
				else if (xd2 > xd2s[j-1] && xd2 < xd2s[j])
				{
					z = ((zs[0][j] - zs[0][j - 1]) / (xd2s[j] - xd2s[j - 1])) * (xd2 - xd2s[j]) + zs[0][j];
					break;
				}
			}
		}
	}
	else if (xd1 >=xd1s[x1n-1])
	{
		if (xd2 <= xd2s[0])
		{
			z = zs[x1n-1][0];
		}
		else if (xd2 >= xd2s[x2n - 1])
		{
			z = zs[x1n-1][x2n - 1];
		}
		else
		{
			for (j = 1; j < x2n; j++)
			{
				if (xd2 == xd2s[j - 1])
				{
					z = zs[x1n-1][j - 1];
					break;
				}
				else if (xd2 > xd2s[j - 1] && xd2 < xd2s[j])
				{
					z = ((zs[x1n-1][j] - zs[x1n-1][j - 1]) / (xd2s[j] - xd2s[j - 1])) * (xd2 - xd2s[j]) + zs[x1n-1][j];
					break;
				}
			}
		}
	}
	else
	{
		if (xd2 <= xd2s[0])
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1==xd1s[i-1])
				{
					z = zs[i - 1][0];
					break;
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					z = ((zs[i][0] - zs[i - 1][0]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][0];
					break;
				}
			}
		}
		else if (xd2 >=xd2s[x2n-1])
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1==xd1s[i-1])
				{
					z = zs[i - 1][x2n - 1];
					break;
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					z = ((zs[i][x2n - 1] - zs[i - 1][x2n - 1]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][x2n - 1];
					break;
				}
			}
		}
		else
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1 == xd1s[i-1])
				{
					for (j=1; j<x2n; j++)
					{
						if (xd2==xd2s[j-1])
						{
							z = zs[i - 1][j - 1];
							break;
						}
						else if (xd2 > xd2s[j-1] && xd2 < xd2s[j])
						{
							z = ((zs[i - 1][j] - zs[i - 1][j - 1]) / (xd2s[j] - xd2s[j - 1])) * (xd2 - xd2s[j]) + zs[i - 1][j];
							break;
						}
					}
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					double lambda1;
					lambda1 = (xd1s[i] - xd1) / (xd1s[i] - xd1s[i - 1]);

					for (j=1; j<x2n; j++)
					{
						if (xd2 == xd2s[j-1])
						{
							z = ((zs[i][j - 1] - zs[i - 1][j - 1]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][j - 1];
							break;
						} 
						else if (xd2 > xd2s[j-1] && xd2 < xd2s[j])
						{
							double lambda2;
							lambda2 = (xd2s[j] - xd2) / (xd2s[j] - xd2s[j - 1]);

							z = (1.0 - lambda1) * (1.0 - lambda2) * zs[i][j] + (1.0 - lambda1) * lambda2 * zs[i][j - 1] + lambda1 * (1.0 - lambda2) * zs[i - 1][j] + lambda1 * lambda2 * zs[i - 1][j - 1];
							break;
						}
					}
				}
			}
		}
	}

	return z;
}

double CInterpolation::operator()(CDate xd1, double x2)
{
	double z = zs[0][0];
	int i, j, x1n = int(xd1s.size()), x2n = int(x2s.size());

	if (xd1 <= xd1s[0])
	{
		if (x2 <= x2s[0])
		{
			z = zs[0][0];
		}
		else if (x2>=x2s[x2n-1])
		{
			z = zs[0][x2n-1];
		}
		else
		{
			for (j=1; j<x2n; j++)
			{
				if (x2 == x2s[j-1])
				{
					z = zs[0][j - 1];
					break;
				}
				else if (x2 > x2s[j-1] && x2 < x2s[j])
				{
					z = ((zs[0][j] - zs[0][z - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[0][j];
					break;
				}
			}
		}
	}
	else if (xd1>=xd1s[x1n-1])
	{
		if (x2 <= x2s[0])
		{
			z = zs[x1n - 1][0];
		}
		else if (x2 >= x2s[x2n-1])
		{
			z = zs[x1n - 1][x2n - 1];
		}
		else
		{
			for (j=1; j<x2n; j++)
			{
				if (x2 == x2s[j-1])
				{
					z = zs[x1n - 1][j - 1];
					break;
				}
				else if (x2 > x2s[j-1] && x2 < x2s[j])
				{
					z = ((zs[x1n - 1][j] - zs[x1n - 1][j - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[x1n - 1][j];
					break;
				}
			}
		}
	}
	else
	{
		if (x2<=x2s[0])
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1 == xd1s[i-1])
				{
					z = zs[i - 1][0];
					break;
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					z = ((zs[i][0] - zs[i - 1][0]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][0];
					break;
				}
			}
		}
		else if (x2 >= x2s[x2n-1])
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1 == xd1s[i-1])
				{
					z = zs[i - 1][x2n - 1];
					break;
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					z = ((zs[i][x2n - 1] - zs[i - 1][x2n - 1]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][x2n - 1];
					break;
				}
			}
		}
		else
		{
			for (i=1; i<x1n; i++)
			{
				if (xd1 == xd1s[i-1])
				{
					for (j=1; j<x2n; j++)
					{
						if (x2 == x2s[j-1])
						{
							z = zs[i - 1][j - 1];
							break;
						}
						else if (x2 > x2s[j-1] && x2 < x2s[j])
						{
							z = ((zs[i - 1][j] - zs[i - 1][j - 1]) / (x2s[j] - x2s[j - 1])) * (x2 - x2s[j]) + zs[i - 1][j];
							break;
						}
					}
				}
				else if (xd1 > xd1s[i-1] && xd1 < xd1s[i])
				{
					double lambda1;
					lambda1 = (xd1s[i] - xd1) / (xd1s[i] - xd1s[i - 1]);

					for (j=1; j<x2n; j++)
					{
						if (x2 == x2s[j-1])
						{
							z = ((zs[i][j - 1] - zs[i - 1][j - 1]) / (xd1s[i] - xd1s[i - 1])) * (xd1 - xd1s[i]) + zs[i][j - 1];
							break;
						}
						else if (x2 > x2s[j-1] && x2 < x2s[j])
						{
							double lambda2;
							lambda2 = (x2s[j] - x2) / (x2s[j] - x2s[j - 1]);

							z = (1.0 - lambda1) * (1.0 - lambda2) * zs[i][j] + (1.0 - lambda1) * lambda2 - zs[i][j - 1] + lambda1 * (1.0 - lambda2) * zs[i - 1][j] + lambda1 * lambda2 * zs[i - 1][j - 1];
							break;
						}
					}
				}
			}
		}
	}

	return z;
}


CInterpolation CInterpolation::operator+(CInterpolation& zc)
{
	int i;

	if (zc.zs.size() > 0)
	{
		int j;
		int newintp1_size, newintp2_size;

		if (xd1s.size() > 0)
		{
			if (xd2s.size() > 0)
			{
				newintp1_size = int(xd1s.size() + zc.xd1s.size());
				newintp2_size = int(xd2s.size() + zc.xd2s.size());

				vector<CDate> newxd1s(newintp1_size), newxd2s(newintp2_size);

				merge(xd1s.begin(), xd1s.end(), zc.xd1s.begin(), zc.xd1s.end(), newxd1s.begin());
				merge(xd2s.begin(), xd2s.end(), zc.xd2s.begin(), zc.xd2s.end(), newxd2s.begin());

				vector<vector<double>> newzs(newintp1_size);

				for (i=0; i<newintp1_size; i++)
				{
					newzs[i] = vector<double>(newintp2_size);

					for (j=0; i<newintp2_size; j++)
					{
						newzs[i][j] = (*this)(newxd1s[i], newxd2s[j]) + zc(newxd1s[i], newxd2s[j]);
					}
				}

				CInterpolation newzc(newxd1s, newxd2s, newzs);
				return newzc;

			}
			else
			{
				newintp1_size = int(xd1s.size() + zc.xd1s.size());
				newintp2_size = int(x2s.size() + zc.x2s.size());

				vector<CDate> newxd1s(newintp1_size);
				vector<double> newx2s(newintp2_size);

				merge(xd1s.begin(), xd1s.end(), zc.xd1s.begin(), zc.xd1s.end(), newxd1s.begin());
				merge(x2s.begin(), x2s.end(), zc.x2s.begin(), zc.x2s.end(), newx2s.begin());

				vector<vector<double>> newzs(newintp1_size);

				for (i = 0; i < newintp1_size; i++)
				{
					newzs[i] = vector<double>(newintp2_size);

					for (j = 0; i < newintp2_size; j++)
					{
						newzs[i][j] = (*this)(newxd1s[i], newx2s[j]) + zc(newxd1s[i], newx2s[j]);
					}
				}

				CInterpolation newzc(newxd1s, newx2s, newzs);
				return newzc;

			}
		}
		else
		{
			newintp1_size = int(x1s.size() + zc.x1s.size());
			newintp2_size = int(x2s.size() + zc.x2s.size());

			vector<double> newx1s(newintp1_size), newx2s(newintp2_size);

			merge(x1s.begin(), x1s.end(), zc.x1s.begin(), zc.x1s.end(), newx1s.begin());
			merge(x2s.begin(), x2s.end(), zc.x2s.begin(), zc.x2s.end(), newx2s.begin());

			vector<vector<double>> newzs(newintp1_size);

			for (i = 0; i < newintp1_size; i++)
			{
				newzs[i] = vector<double>(newintp2_size);

				for (j = 0; i < newintp2_size; j++)
				{
					newzs[i][j] = (*this)(newx1s[i], newx2s[j]) + zc(newx1s[i], newx2s[j]);
				}
			}

			CInterpolation newzc(newx1s, newx2s, newzs);
			return newzc;

		}
	}
	else
	{
		int newintp_size = int(ys.size() + zc.ys.size());

		vector<double> newys(newintp_size);

		if (xs.size()>0)
		{
			vector<double> newxs(newintp_size);

			merge(xs.begin(), xs.end(), zc.xs.begin(), zc.xs.end(), newxs.begin());

			for (i=0; i<newintp_size; i++)
			{
				newys[i] = (*this)(newxs[i]) + zc(newxs[i]);
			}

			CInterpolation newzc(newxs, newys);
			return newzc;
		}
		else
		{
			vector<CDate> newxds(newintp_size);

			merge(xds.begin(), xds.end(), zc.xds.begin(), zc.xds.end(), newxds.begin());

			if (zc.xds.size()>0)
			{
				for (i=0; i<newintp_size; i++)
				{
					newys[i] = (*this)(newxds[i]) + zc(newxds[i]);
				}

				CInterpolation newzc(newxds, newys);
				return newzc;
			}
			else
			{
				CInterpolation newzc = (*this);
				return newzc;
			}
		}
	}
}

int CInterpolation::get_xs_size()
{
	return int(xs.size());
}

int CInterpolation::get_xds_size()
{
	return int(xds.size());
}

int CInterpolation::get_x1s_size()
{
	return int(x1s.size());
}

int CInterpolation::get_xd1s_size()
{
	return int(xd1s.size());
}

int CInterpolation::get_x2s_size()
{
	return int(x2s.size());
}

int CInterpolation::get_xd2s_size()
{
	return int(xd2s.size());
}

vector<double> CInterpolation::get_xs()
{
	return xs;
}

vector<double> CInterpolation::get_ys()
{
	return ys;
}

vector<CDate> CInterpolation::get_xds()
{
	return xds;
}

vector<double> CInterpolation::get_x1s()
{
	return x1s;
}

vector<CDate> CInterpolation::get_xd1s()
{
	return xd1s;
}

vector<double> CInterpolation::get_x2s()
{
	return x2s;
}

vector<CDate> CInterpolation::get_xd2s()
{
	return xd2s;
}

vector<vector<double>> CInterpolation::get_zs()
{
	return zs;
}

CInterpolation::~CInterpolation()
{

}