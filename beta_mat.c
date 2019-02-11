//
//	Bマトリックスの計算
//

#include "typedef.h"
#include "datadef.h"

/* ヤコビアンの逆行列を計算する。ともに行列値を返す。 */
double invJ_2D (double *iJ, const double *J)
{
	double detJ;
	detJ = J [0] * J [3] - J [2] * J [1];
	iJ [0] =  J [3] / detJ;
	iJ [1] = -J [1] / detJ;
	iJ [2] = -J [2] / detJ;
	iJ [3] =  J [0] / detJ;
	return detJ;
}

double invJ_3D (double *iJ, const double *J)
{
	double detJ;
	detJ =    J [0] * (J [4] * J [8] - J [5] * J [7])
		+ J [1] * (J [5] * J [6] - J [3] * J [8])
		+ J [2] * (J [3] * J [7] - J [4] * J [6]);
	iJ [0] = (J [4] * J [8] - J [5] * J [7]) /detJ;
	iJ [1] = (J [7] * J [2] - J [8] * J [1]) /detJ;
	iJ [2] = (J [1] * J [5] - J [2] * J [4]) /detJ;
	iJ [3] = (J [5] * J [6] - J [3] * J [8]) /detJ;
	iJ [4] = (J [8] * J [0] - J [6] * J [2]) /detJ;
	iJ [5] = (J [2] * J [3] - J [0] * J [5]) /detJ;
	iJ [6] = (J [3] * J [7] - J [4] * J [6]) /detJ;
	iJ [7] = (J [6] * J [1] - J [7] * J [0]) /detJ;
	iJ [8] = (J [0] * J [4] - J [1] * J [3]) /detJ;
	return detJ;
}

/* ヤコビアンの計算 形状関数の微分形Hに要素の節点座標値をかけることで得られる。 */
void calc_Jmat (double *J, const double *H, const Element *e)
{
	int i, j, k;
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	memset (J, 0, sizeof (double) * ndim * ndim);
	for (i = 0; i < ndim; i++)
		for (j = 0; j < ndim; j++)
			for (k = 0; k < nnode; k++)
				J [i * ndim + j] += H [i * nnode + k] * node [e->conn [k]].coord [j];
}

/*	二次元配列の掛け算	*/
void mat_mul2 (double *A, const double *B, const double *C, int sz1, int sz2, int sz3)
{
	int i, j, k;
	memset (A, 0, sizeof (double) * sz1 * sz3);
	for (i = 0; i< sz1; i++)
		for (j = 0; j < sz2; j++)
			for (k = 0; k < sz3; k++)
				A [i * sz2 + j] += B [i * sz3 + k] * C [k * sz2 + j];
}

//	2次元要素のAマトリックスを作る。
//	副産物として得られるヤコビアンは要素計算情報に格納される。
void calc_PL_iJH (const Element *e, const double *N, const double *H, int nint)
{
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	ElementCalcInfo *cp = e->cinfo + nint;
	double J [ndim * ndim], iJ [ndim * ndim];

	calc_Jmat (J, H, e);
	cp->detJ = invJ_2D (iJ, J);
	mat_mul2 (cp->iJH, iJ, H, ndim, nnode, ndim);
}

//	2次元平面応力・ひずみ要素のBマトリックスを求める。
void calc_PL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
{
	int i;
	int nnode = e->info->info1->nnode;
	int size = ndim * nnode;
	ElementCalcInfo *cp = e->cinfo + nint;
	calc_PL_B_mat (e, N, H, nint);
	// Bマトリックスの作成。一つ飛ばし(つまり自由度分)に値を代入している。
	for (i = 0; i < nnode; i++) {
		B [0 * size + i * 2 + 0] = B [2 * size + i * 2 + 1] = cp->iJH [i];
		B [2 * size + i * 2 + 0] = B [1 * size + i * 2 + 1] = cp->iJH [i + nnode];
		B [1 * size + i * 2 + 0] = B [0 * size + i * 2 + 1] = 0;
	}
}

//軸対称要素のBマトリックスを求める。
void calc_AXSOL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
{
	int i;
	int nnode = e->info->info1->nnode;
	int size = ndim * nnode;
	ElementCalcInfo *cp = e->cinfo + nint;
	double cx;

	calc_PL_B_mat (e, N, H, nint);
	// 積分点での半径を求める。
	cx = 0;
	for (i = 0; i < nnode; i++)
		cx += N [i] * node [e->conn [i]].coord [0];
	cp->iJH [nnode * 3] = cx;
	for (i = 0; i < nnode; i++) {
		B [0 * size + i * 2 + 0] = B [3 * size + i * 2 + 1] = cp->iJH [i];
		B [3 * size + i * 2 + 0] = B [1 * size + i * 2 + 1] = cp->iJH [i + nnode];
		B [2 * size + i * 2 + 0]                            = cp->iJH [i + nnode * 2] = N [i] / cx;
		B [1 * size + i * 2 + 0] = B [0 * size + i * 2 + 1] = B [2 * size + i * 2 + 1] = 0;
	}
}


