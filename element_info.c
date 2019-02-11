//
//要素情報構造体関連のデーターとルーチンの定義
//

#include "typedef.h"
#include "datadef.h"
#include "shapefunc.h"
#include "integral.h"

//MODULE INIT_ELEMENT_INFO
/*
 * 要素情報テーブルの定義
*/
ElementTypeInfo etinfo [] = {
	/* number   name                prob    dim dof stat scl add */
	{ 0,        "None",             PT_None,  0,  0,   0,  0, 0, }, 
	{ PT_PSS,   "Plane stress",     PT_PSS,   2,  2,   3,  2, 0, }, 
	{ PT_PSN,   "Plane strain",     PT_PSN,   2,  2,   3,  2, 0, }, 
	{ PT_AXSOL, "Aximentric solid", PT_AXSOL, 2,  2,   4,  3, 1, },
};

ElemIntegralInfo eiinfo [] = {
/*        num nod ord int dim sf1   sf2    N1    dN1   ipt1    wt1     */
	{  0,  0, 0,  0,  0,  NULL, NULL,  NULL, NULL, NULL,   NULL,   },
	{  1,  2, 1,  2,  1,  LL1,  dLL1,  NULL, NULL, LL1_it, LL1_wt, },
	{  2,  3, 2,  3,  1,  LL2,  dLL2,  NULL, NULL, LL2_it, LL2_wt, },
	{  3,  4, 3,  4,  1,  LL3,  dLL3,  NULL, NULL, LL3_it, LL3_wt, },
	{  4,  3, 1,  1,  2,  Tr1,  dTr1,  NULL, NULL, Tr1_it, Tr1_wt, },
	{  5,  6, 2,  3,  2,  Tr2,  dTr2,  NULL, NULL, Tr2_it, Tr2_wt, },
	{  6, 10, 3,  7,  2,  Tr3L, dTr3L, NULL, NULL, Tr3_it, Tr3_wt, },
	{  7,  9, 3,  7,  2,  Tr3S, dTr3S, NULL, NULL, Tr3_it, Tr3_wt, },
	{  8,  4, 1,  4,  2,  Q1,   dQ1,   NULL, NULL, Q1_it,  Q1_wt,  },
	{  9,  9, 2,  9,  2,  Q2L,  dQ2L,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 10,  8, 2,  9,  2,  Q2S,  dQ2S,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 11, 16, 3, 16,  2,  Q3L,  dQ3L,  NULL, NULL, Q3_it,  Q3_wt,  },
	{ 12, 12, 3, 16,  2,  Q3S,  dQ3S,  NULL, NULL, Q3_it,  Q3_wt,  },
	{ 13,  6, 2,  1,  2,  Tr2,  dTr2,  NULL, NULL, Tr1_it, Tr1_wt, },
	{ 14, 10, 3,  3,  2,  Tr3L, dTr3L, NULL, NULL, Tr2_it, Tr2_wt, },
	{ 15,  9, 3,  3,  2,  Tr3S, dTr3S, NULL, NULL, Tr2_it, Tr2_wt, },
	{ 16,  4, 1,  1,  2,  Q1,   dQ1,   NULL, NULL, Q0_it,  Q0_wt,  },
	{ 17,  9, 2,  4,  2,  Q2L,  dQ2L,  NULL, NULL, Q1_it,  Q1_wt,  },
	{ 18,  8, 2,  4,  2,  Q2S,  dQ2S,  NULL, NULL, Q1_it,  Q1_wt,  },
	{ 19, 16, 3,  9,  2,  Q3L,  dQ3L,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 20, 12, 3,  9,  2,  Q3S,  dQ3S,  NULL, NULL, Q2_it,  Q2_wt,  },
};

ElementInfo einfo [] = {
/*        num name    vert face info1         info2         info3 */
	{  0, NULL,      0,  0,  NULL,         NULL,         NULL,        },
	{  1, "Tria1",   3,  3,  &eiinfo [ 4], NULL,         &eiinfo [1], },
	{  2, "Tria2",   3,  3,  &eiinfo [ 5], NULL,         &eiinfo [2], },
	{  3, "Tria3L",  3,  3,  &eiinfo [ 6], NULL,         &eiinfo [3], },
	{  4, "Tria3S",  3,  3,  &eiinfo [ 7], NULL,         &eiinfo [3], },
	{  5, "Quad1",   4,  4,  &eiinfo [ 8], NULL,         &eiinfo [1], },
	{  6, "Quad2L",  4,  4,  &eiinfo [ 9], NULL,         &eiinfo [2], },
	{  7, "Quad2S",  4,  4,  &eiinfo [10], NULL,         &eiinfo [2], },
	{  8, "Quad3L",  4,  4,  &eiinfo [11], NULL,         &eiinfo [3], },
	{  9, "Quad3S",  4,  4,  &eiinfo [12], NULL,         &eiinfo [3], },
	{ 10, "Tria2R",  3,  3,  &eiinfo [13], NULL,         &eiinfo [2], },
	{ 11, "Tria3LR", 3,  3,  &eiinfo [14], NULL,         &eiinfo [3], },
	{ 12, "Tria3SR", 3,  3,  &eiinfo [15], NULL,         &eiinfo [3], },
	{ 13, "Quad1R",  4,  4,  &eiinfo [16], NULL,         &eiinfo [1], },
	{ 14, "Quad2LR", 4,  4,  &eiinfo [17], NULL,         &eiinfo [2], },
	{ 15, "Quad2SR", 4,  4,  &eiinfo [18], NULL,         &eiinfo [2], },
	{ 16, "Quad3LR", 4,  4,  &eiinfo [19], NULL,         &eiinfo [3], },
	{ 17, "Quad3SR", 4,  4,  &eiinfo [20], NULL,         &eiinfo [3], },
};

int netinfo = sizeof (etinfo) / sizeof (etinfo [0]);
int neiinfo = sizeof (eiinfo) / sizeof (eiinfo [0]);
int neinfo  = sizeof (einfo)  / sizeof (einfo [0]);

/*
 * 関数の定義
 */

/* 要素積分計算情報構造体の初期化 */
void init_einfo ()
{
		// 0番目のエントリーは使用しないので、1を足して1番目からアクセスするようにしている。
	ElemIntegralInfo *eip = eiinfo + 1;
	int i, j;
	for (i = 1; i < neiinfo; i++, eip++) {
		int nnode = eip->nnode;
		int ndim = eip->ndim;
		int ntint = eip->ntint;
		double *ipt = eip->ipt;
		// マトリックス領域の割り当て
		eip->N  = malloc (sizeof (double) * ntint * nnode);
		eip->dN = malloc (sizeof (double) * ntint * nnode * ndim);
		for (j = 0; j < eip->ntint; j++) {
			(*eip->sf)  (eip->N  + j * nnode,        ipt + j * ndim);
			(*eip->dsf) (eip->dN + j * nnode * ndim, ipt + j * ndim);
		}
	}
}

