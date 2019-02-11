/* 境界条件の作成
   荷重は計算用配列に入れる
 */

#include "typedef.h"
#include "datadef.h"

/* 荷重データーの計算 */
static void calc_load_value (const LoadData *p, double *load)
{
	double v;
	int i, ndim;
	ndim = etinfo [prob].ndim;
	v = 1;
	switch (p->val_type % 10) {
	case BC_VAL:
	case BC_NORMAL:
		for (i = 0; i < ndim; i++)
			load [i] = p->value [i];
		break;
	case BC_VEC:
		for (i = 0, v = 0; i < ndim; i++)
			v += p->value [i] * p->value [i];
		v = sqrt (v);
		// 0ベクトルなら値の設定をせずルーチンから抜ける
		if (fabs (v) <= tolerance) break;
	case BC_VAL_VEC: // fall thrurh
		for (i = 0; i < ndim; i++)
			load [i] = p->value [ndim] * p->value [i] / v;
		break;
	}
}

// 節点荷重の設定
static void set_pload (const LoadData *p)
{
	int i, j, pos, rpos;
	double f;
	for (i = 0; i < mdof; i++) {
		if (fabs (p->load [i]) < tolerance) continue;
		for (j = 0; j < p->ndata; j++) {
			pos = p->data [j] * mdof + i;
			rpos = rev_line_info [pos];
			if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
			f = p->load [i];
			// 軸対称問題の場合、荷重は必ず半径を掛けるため、荷重点のX座標を取り出すようにしている
			if (prob == PT_AXSOL) f *= coord [p->data [j] * mdof];
			lhs_value [rpos] = rhs_value [rpos] += f;
		}
	}
}

/* 辺荷重の設定*/
static void calc_edge_load (LoadData *p, const Element *e, int *face, int nint)
{
	int i, j, k;
	ElemIntegralInfo *info = e->info->info3;
	int ndim = etinfo [prob].ndim, ndof = etinfo [prob].ndof;
	int nnode = info->nnode;
	int pos, rpos;
	double f, J, x [ndim], n [ndim];
	double *H = info->dN + nint * nnode;
	double *N = info->N + nint * nnode;
	double wt = info->wt [nint];
	// 接線ベクトルとヤコビアンを求める
	memset (x, 0, sizeof (double) * ndim);
	for (J = i = 0; i < ndim; i++) {
		for (j = 0; j < nnode; j++) {
			x [i] += H [j] * node [e->conn [face [j]]].coord [i];
		}
		J += x [i] * x [i];
	}
	J = sqrt (J);
	if (fabs(J) < tolerance) return;
	for (i = 0; i < ndim; i++) x [i] /= J;	// 単位ベクトルにする
	if (p->val_type % 10 == BC_NORMAL) {		// 面に垂直な荷重の向きを求める
		if (ndim == 2)
			n [0] = x [1], n [1] = -x [0];
		if (ndim == 3)
			;	// 未実装
		for (i = 0; i < nnode; i++) {
			for (j = 0; j < ndof; j++) {
				pos = e->conn [face [i]] * ndof + j;
				rpos = rev_line_info [pos];
				if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
				f = J * wt * N [i] * (p->load [0] * n [j] + p->load [1] * x [j]);
				if (prob == PT_AXSOL) f *= x [0] * J;
				lhs_value [rpos] = rhs_value [rpos] += f;
			}
		}
	} else {
		for (i = 0; i < nnode; i++) {
			for (j = 0; j < ndof; j++) {
				pos = e->conn [face [i]] * ndof + j;
				rpos = rev_line_info [pos];
				if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
				f = J * wt * N [i] * p->load [j];
				if (prob == PT_AXSOL) f *= x [0] * J;
				lhs_value [rpos] = rhs_value [rpos] += f;
			}
		}
	}
}

static void set_fload (LoadData *p)
{
	int i, j, nnode, *elist;
	Element *e;
	for (i = 0; i < p->ndata; i++) {
		e = elem + p->data [i];
		nnode = e->info->info3->nnode;
		elist = e->info->edge_list + nnode * p->data [i * 2 + 1];
		for (j = 0; j < e->info->info3->ntint; j++)
			calc_edge_load (p, e, elist, j);
	}
}

// 積分点ごとの体積力を求める
static void calc_bload (const LoadData *p, const Element *e, int nint)
{
	int i, j;
	ElemIntegralInfo *info = e->info->info4;
	int ndim = etinfo [prob].ndim, nnode = info->nnode;
	int ndof = etinfo [prob].ndof;
	int pos, rpos;
	double f, detJ, cx, J [ndim * ndim], iJ [ndim * ndim];
	double *H = info->dN + nint * nnode * ndim;
	double *N = info->N + nint * nnode;
	double wt = info->wt [nint];
	// ヤコビアンを求める。
	calc_Jmat (J, H, e);
	detJ = invJ_2D (iJ, J);
	f = detJ * wt;
	// 自重を考慮。
	if (p->val_type / 10 == BC_DENSITY) f *= e->part->mat->density;
	// 軸対称問題では積分点の半径をかける必要がある。
	if (prob == PT_AXSOL) {
		cx = 0;
		for (j = 0; j < nnode; j++)	// 積分点での半径を求める
			cx += N [j] * coord [e->conn [j] * mdof];
		f *= cx;
	}
	for (j = 0; j < nnode; j++)
		for (i = 0; i < ndof; i++) {
			pos = e->conn [j] * ndof + i;
			rpos = rev_line_info [pos];
			if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
			if (fabs (p->load [i]) < tolerance) continue;
			lhs_value [rpos] = rhs_value [rpos] += p->load [i] * N [j] * f;
		}
}

static void set_bload (const LoadData *p)
{
	int i, j;
	for (i = 0; i < p->ndata; i++) {
		Element *e = &elem [p->data [i]];
		for (j = 0; j < e->info->info4->nnode; j++)
			calc_bload (p, e, j);
	}
}

void calc_load ()
{
	int i;
	LoadData *p;
	for (i = 0, p = pload; i < npload; i++, p++) {
		calc_load_value (p, p->load);
		set_pload (p);
	}
	for (i = 0, p = fload; i < nfload; i++, p++) {
		calc_load_value (p, p->load);
		set_fload (p);
	}
	for (i = 0, p = bload; i < nbload; i++, p++) {
		calc_load_value (p, p->load);
		set_bload (p);
	}
}

