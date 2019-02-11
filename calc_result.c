//
//
//

#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "datadef.h"

static void calc_PL_stress_strain ()
{
	int i, j, k, l, m, nnode, ntint, ndof, nvstat, ndim;
	double f, *D, *wt, *H1, *H2, *H3, *stress, *strain;
	Element *e;
	ElementCalcInfo *cp;

	nvstat = etinfo [prob].nvstat;
	ndof = etinfo [prob].ndof;
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		ElemIntegralInfo *info = e->info->info1;
		nnode = info->nnode; ntint = info->ntint;
		ndim = etinfo [prob].ndim;
		wt = info->wt;
		IsoMatData *mat = (IsoMatData *)(e->part->mat);
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			H1 = cp->iJH; H2 = H1 + nnode; H3 = H2 + nnode;
			strain = cp->strain; stress = cp->stress;
			memset (strain, 0, sizeof (double) * 6);
			memset (stress, 0, sizeof (double) * 6);
			//fprintf (stderr, "%5d\n", prob);
			// 応力とひずみの計算。
			for (k = 0; k < nnode; k++) {
				int pos1 = e->conn [k] * mdof;
				int pos2 = pos1 + 1;
				strain [0] += H1 [k] * disp [pos1];
				strain [1] += H2 [k] * disp [pos2];
				strain [3] += H2 [k] * disp [pos1] + H1 [k] * disp [pos2];
				if (prob == PT_AXSOL)
					strain [2] += H3 [k] * disp [pos1];
			}
			if (prob == PT_PSS || prob == PT_PSN)
				strain [2] = strain [3];
			D = e->part->init_mat;
			for (k = 0; k < nvstat; k++) {
				for (l = 0; l < nvstat; l++, D++)
					stress [k] += *D * strain [l];
			}
			if (prob == PT_PSS || prob == PT_PSN) {
				stress [3] = stress [2];
				switch (prob) {
					case PT_PSS:
						stress [2] = 0;
						strain [2] = mat->poisson / mat->young * (stress [0] + stress [1]);
						break;
					case PT_PSN:
						stress [2] = mat->poisson * (stress [0] + stress [1]);
						strain [2] = 0;
						break;
				}
			}
			// 節点力の計算。外力と反力が得られる。
			f = cp->detJ * wt [j];
			for (k = 0; k < nnode; k++) {
				double f2 = f;
				int pos1 = e->conn [k] * ndof;
				int pos2 = pos1 + 1;
				if (prob == PT_AXSOL) f2 *= cp->iJH [nnode * 3];
				force [pos1] += f2 * (H1 [k] * stress [0] + H2 [k] * stress [3]);
				force [pos2] += f2 * (H2 [k] * stress [1] + H1 [k] * stress [3]);
				if (prob == PT_AXSOL) force [pos1] += f2 * (H3 [k] * stress [2]);
			}
			// 変形勾配テンソルの計算
			// F=(J^-1H)^Tで求められる。
			memset (cp->FT, 0, sizeof (double) * 9);
			for (k = 0; k < ndim; k++)
				for (l = 0; l < ndim; l++)
					for (m = 0; m < nnode; m++) {
						int node_id = e->conn [m];
						double xi = node [node_id].coord [l] + disp [node_id * ndim + l];
						cp->FT [l * 3 + k] += cp->iJH [k * nnode + m] * xi;
					}
			if (prob <= PT_AXSOL) cp->FT [8] = 1;
		}
	}
}

void calc_result ()
{
	int i, j, pos2, pos;
	// 節点変位と荷重の整列
	for (i = 0; i < ntnode; i++) {
		for (j = 0; j < mdof; j++) {
			pos = i * mdof + j;
			if (pos < rank_line_info [KLDOF_NONE][1]) {
				pos2 = decode_rev_info [pos];
				disp [pos2] = rhs_value [pos];
			}
		}
	}
	calc_PL_stress_strain ();
}

/* 結果出力 */
void write_result (FILE *fout)
{
	int i, j, k, l;
	double *p;
	Element *e;
	ElementCalcInfo *cp;

	fprintf (fout, "Nodal displacement\n%10s%15s%15s%15s\n", "Number", "X-dir", "Y-dir", "Z-dir");
	for (i = 0, p = disp; i < ntnode; i++, p += mdof) {
		fprintf (fout, "%10d", node [i].number);
		for (j = 0; j < mdof; j++)
			fprintf (fout, "%15g", p [j]);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
	fprintf (fout, "Reaction/External force\n%10s%15s%15s%15s\n", "Number", "X-dir", "Y-dir", "Z-dir");
	for (i = 0, p = force; i < ntnode; i++, p += mdof) {
		fprintf (fout, "%10d", node [i].number);
		for (j = 0; j < mdof; j++)
			fprintf (fout, "%15g", p [j]);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
	fprintf (fout, "Small strain\n%10s%10s%15s%15s%15s%15s%15s%15s\n",
			"Number", "IPT", "X-dir", "Y-dir", "Z-dir", "XY-shar", "YZ-shar", "ZX-shar");
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			fprintf (fout, "%10d%10d", e->number, j + 1);
			for (k = 0; k < 6; k++)
				fprintf (fout, "%15g", cp->strain [k]);
			fputc ('\n', fout);
		}
	}
	fputc ('\n', fout);
	fprintf (fout, "Cauchy stress\n%10s%10s%15s%15s%15s%15s%15s%15s\n",
			"Number", "IPT", "X-dir", "Y-dir", "Z-dir", "XY-shar", "YZ-shar", "ZX-shar");
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			fprintf (fout, "%10d%10d", e->number, j + 1);
			for (k = 0; k < 6; k++)
				fprintf (fout, "%15g", cp->stress [k]);
			fputc ('\n', fout);
		}
	}
	fputc ('\n', fout);
	fprintf (fout, "Displacement gradient tensor\n%10s%10s%15d%15d%15d\n", "Number", "IPT", 1, 2, 3);
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			for (k = 0; k < 3; k++) {
				if (k == 0)
					fprintf (fout, "%10d%10d", e->number, j + 1);
				else
					fprintf (fout, "%20c", ' ');
				for (l = 0; l < 3; l++)
					fprintf (fout, "%15g", cp->FT [k * 3 + l]);
				fputc ('\n', fout);
			}
		}
	}
	fputc ('\n', fout);
}

