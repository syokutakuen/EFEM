//
//	全体剛性マトリックスを組み立てる。
//

#include "typedef.h"
#include "datadef.h"

//	ひずみ変位行列と要素剛性の宣言。結構でかいのでここで宣言する。
static double __DB__ [6 * 64 * sizeof (double)];
static double __B__ [6 * 64 * sizeof (double)];

/* 要素に接続していない節点を探す */
static void set_connected_dof ()
{
	int i, j, k, nid;
	int ndof = etinfo [prob].ndof;
	Element *p;

	// 初めに水子節点であると仮定する。
	for (i = 0; i < ntnode * mdof; i++)
		nodal_dof_list [i] = KLDOF_ORPHAN;
	// コネクティビティのチェック。要素に接続する節点にKLDOF_NONEをセット
	for (i = 0, p = elem; i < ntelem; i++, p++)
		for (j = 0; j < p->info->info1->nnode; j++) {
			nid = p->conn [j];
			for (k = 0; k < ndof; k++)
				nodal_dof_list [nid * mdof + k] = KLDOF_NONE;
		}
}

/* 境界条件を適用する自由度に印をつける */
void set_bc_pos ()
{
	int i, j, k, pos;
	FixDisp *fd;

	for (i = 0, fd = fdisp; i < nfdisp; i++, fd++) {
		for (j = 0; j < mdof; j++) {
			if (!fd->flags [j]) continue;
			for (k = 0; k < fd->ndata; k++) {
				int np = fd->data [k];
				pos = np * mdof + j;
				//
				nodal_dof_list [pos] = (fd->value [j] > 0)? KLDOF_NONZERO: KLDOF_ZERO;
				if (fd->value [j] > 0) disp [pos] += fd->value [j];
			}
		}
	}
}

/* 非拘束自由度・拘束自由度・水子節点の順に節点番号と自由度を並べる */
void create_rev_line_info ()
{
	int i, j, kltype, pos, count;
	// 最初と最後の位置は黙って決まるのでこの時点で設定
	rank_line_info [KLDOF_NONE][0] = 0;
	rank_line_info [KLDOF_ORPHAN][1] = ntnode * mdof;
	for (kltype = KLDOF_NONE, count = 0; kltype < KLDOF_ORPHAN; kltype++) {
		rank_line_info [kltype][1] = rank_line_info [kltype][0];
		// 拘束タイプの検索と配置。配列の[1]-[0]が個数になる。
		// 従って各データーの終点位置は実際の終点位置の一つ前にある。
		// よって順次アクセスする場合はCの配列に倣うようにしている。
		for (i = 0; i < ntnode; i++) {
			for (j = 0; j < mdof; j++) {
				pos = i * mdof + j;
				if (nodal_dof_list [pos] == kltype) {
					decode_rev_info [count] = pos;
					rev_line_info [pos] = count++;
				}
			}
		}
		// 次のデーターの起点位置は前のデーターの「終点位置」に等しい。
		rank_line_info [kltype + 1][0] = rank_line_info [kltype][1] = count;
	}
	rank_line_info [KLDOF_ORPHAN][0] = rank_line_info [KLDOF_ZERO][1];
}

//	ライン情報の初期化
void init_K_line_info ()
{
	int i, j;
	K_line_info *p;
	for (i = j = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		while (rev_line_info [j] >= rank_line_info [KLDOF_NONE][1])
			j++;
		p->start = p->end = i;
		j++;
	}
}

//	フルマトリックス用の初期化を行う。
void set_range_all ()
{
	int i, size = rank_line_info [KLDOF_NONE][1] - rank_line_info [KLDOF_NONE][0];
	K_line_info *p;
	K_line_info tmp = {
		0,
		rank_line_info [KLDOF_NONE][0],
		rank_line_info [KLDOF_NONE][1],
		size, NULL,
	};

	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		*p = tmp;
		p->pos = size * i;
	}
}

/* 各要素のコネクティビティ―から各ラインにおけるバンド幅を求める。 */
void set_range_by_elem ()
{
	int i, j, k, pos, minpos, maxpos;
	int ndof = etinfo [prob].ndof;
	K_line_info *p, *pp;
	Element *e;

	for (i = 0, e = elem; i < ntelem; i++, e++) {
		ElemIntegralInfo *info = e->info->info1;
		int nnode = info->nnode;
		minpos = ntnode * mdof;
		maxpos = 0;
		// 要素におけるコネクティビティの範囲を探す
		for (j = 0; j < nnode; j++) {
			for (k = 0; k < ndof; k++) {
				pos = rev_line_info [e->conn [j] * mdof + k];
				if (pos >= rank_line_info [KLDOF_NONE][1]) continue;
				if (minpos > pos) minpos = pos;
				if (maxpos < pos) maxpos = pos;
			}
		}
		for (j = 0; j < nnode; j++) {
			for (k = 0; k < ndof; k++) {
				pos = rev_line_info [e->conn [j] * mdof + k];
				if (pos >= rank_line_info [KLDOF_NONE][1]) continue;
				p = line_info + pos;
				if (sys_range_mode == RANGE_BEGIN_END && p->start > minpos)
					p->start = minpos;
				if (p->end < maxpos + 1) p->end = maxpos + 1;
			}
		}
	}
	line_info [0].size = line_info [0].end - line_info [0].start;
	for (i = 1, pp = line_info, p = pp + 1; i < rank_line_info [KLDOF_NONE][1]; i++, p++, pp++) {
		if (pp->end > p->end) p->end = pp->end;
		p->size = p->end - p->start;
	}
}

/* 全体剛性マトリックスのサイズを計算する */
int set_K_size ()
{
	int i, size;
	K_line_info *p;

	for (i = size = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		p->pos = size;
		size += p->size;
	}
	return size;
}

/*
 * 全体剛性のセットアップ
 *   以下に使用メモリが最小限になるようにするか考えた
 * 算法
 *   全自由度を求める
 *   ライン情報のメモリーを割り当てる
 *   ライン情報の設定を行う。初期位置は対角線上にあり、拘束はない状態とする。
 */
void create_K_info (FILE *fout)
{
	int i, size;
	K_line_info *p;

	size = ntnode * mdof;
	nodal_dof_list = calloc (sizeof (int), size);
	disp = calloc (sizeof (double), size);
	force = calloc (sizeof (double), size);
	rev_line_info = calloc (sizeof (int), size);
	decode_rev_info = calloc (sizeof (int), size);

	set_connected_dof ();
	set_bc_pos ();
	create_rev_line_info ();

	line_info = calloc (sizeof (K_line_info), rank_line_info [KLDOF_NONE][1]);
	rhs_value = calloc (sizeof (double), rank_line_info [KLDOF_NONE][1]);
	lhs_value = calloc (sizeof (double), rank_line_info [KLDOF_NONE][1]);
	init_K_line_info ();

	if (sys_range_mode == RANGE_ALL)
		set_range_all ();
	else
		set_range_by_elem ();

	size = set_K_size ();
	fprintf (stderr, "Total matrix size: " LLU0 "bytes\n", size * sizeof (double));
	fprintf (fout, "\nTotal matrix size: " LLU0 "bytes\n", size * sizeof (double));

	sysk = calloc (sizeof (double), size);
	if (sysk == NULL) {
		fputs ("Error: Cannot allocate system stiff matricx\n", stderr);
		free_data ();
		exit (-11);
	}

	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++)
		p->p = sysk + p->pos - p->start;
}

//	各積分点ごとの要素剛性マトリックスを全体剛性に繰り込む
void attach_B_mat (const Element *e, const double *D, int nint)
{
	int i, j, k, kx, ky, ypos, ex, ey, dofx, dofy;
	double Ke, thick;
	ElemIntegralInfo *info = e->info->info1;
	int nvstat = etinfo [prob].nvstat;
	int nnode = info->nnode;
	int ndof = etinfo [prob].ndof;
	int size = nnode * ndof;
	ElementCalcInfo *cp = e->cinfo + nint;
	// 要素剛性に必要な値の計算。積分点での重みとヤコビアンの積を求める。
	double factor = info->wt [nint] * cp->detJ;
	// B^T*Dを計算する。
	memset (__DB__, 0, size * nvstat * sizeof (double));
	mat_mul2 (__DB__, D, __B__, nvstat, size, nvstat);
	// 要素剛性の計算。容量が大きいため一度に全部を計算せず一つ一つ計算する。
	for (i = 0; i < size; i++) {
		ex = e->conn [i / ndof]; dofx = i % ndof;
		kx = rev_line_info [ex * mdof + dofx];
		if (kx >= rank_line_info [KLDOF_NONZERO][0]) continue;
		for (j = 0; j < size; j++) {
			ey = e->conn [j / ndof]; dofy = j % ndof;
			ypos = ey * mdof + dofy;
				ky = rev_line_info [ypos];
			// ここでKe[i][j]のみ計算する
			Ke = 0;
			for (k = 0; k < nvstat; k++)
				Ke += __B__ [k * size + i] * __DB__ [k * size + j];
			// 最終的な要素剛性の計算。平面応力・歪問題と軸対称問題ではかける係数が異なる。
			// 係数は積分点ごとの重みとヤコビアンの積だが、軸対称問題では半径もかける。
			switch (prob) {
			case PT_PSS:
				thick = ((PSS_GeomData *) (e->part->geom))->thick;
				Ke *= thick;
			case PT_PSN:
				Ke *= factor;
				break;
			case PT_AXSOL:
				Ke *= cp->iJH [nnode * 3] * factor;
				break;
			}
			// 要素剛性を全体剛性に組み込む
			// バンドマトリックス計算までは「そのまま」位置を探して組み込むが、
			// 対称マトリックスは対角線の上だけに値を入れるため、位置を調べて該当する値だけをセットする。
			if (ky < rank_line_info [KLDOF_NONE][1]) {
				if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END) {
					line_info [kx].p [ky] += Ke;
				} else if (sys_range_mode == RANGE_SYMM_UPPER) {
					if (ky >= kx)
						line_info [kx].p [ky] += Ke;
				}
			// 変位拘束条件のかかる自由度は、要素剛性と変位をかけたものを荷重ベクトルへ足し込む。
			} else if (ky < rank_line_info [KLDOF_NONZERO][1])
				rhs_value [kx] -= Ke * disp [ypos];
		}
	}
}

//	Bマトリックスを作成するルーチンのテーブル。
typedef void (*FP_calc_B) (double *, const Element *, const double *, const double *, int);
FP_calc_B calc_B [] = {
	NULL,
	calc_PL_B_mat,
	calc_PL_B_mat,
	calc_AXSOL_B_mat,
	NULL,
};

//	最終的な全体剛性の組み立て。
void calc_K_mat (FILE *fout)
{
	int i, j;
	int nvstat = etinfo [prob].nvstat;
	create_K_info (fout);
	calc_load ();
	for (i = 0; i < ntelem; i++) {
		Element *e = elem + i;
		ElemIntegralInfo *info = e->info->info1;
		double *N = info->N;
		double *H = info->dN;
		double *D = e->part->init_mat;
		int ndim = etinfo [prob].ndim, nnode = info->nnode;
		for (j = 0; j < info->ntint; j++) {
			(*calc_B [prob]) (__B__, e, N, H, j);
			attach_B_mat (e, D, j);
			N += nnode;
			H += nnode * ndim;
		}
	}
}

