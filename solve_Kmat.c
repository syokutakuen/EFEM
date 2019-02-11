//
//方程式の求解ルーチン
//

#include <stdio.h>
#include <math.h>

#include "typedef.h"
#include "datadef.h"
#include "memfree.h"

//  ンドマトリックス形式の方程式を解く
void solve_band_mat (double tol)
{
	int i, j, k, size;
	double pivot, Kji;

	// 前進消去
	size = rank_line_info [KLDOF_NONE][1];
	for (i = 0; i < size; i++) {
		pivot = line_info [i].p [i];
		if (fabs (pivot) < tol) {
			fprintf (stderr, "Error: pivot tolerance under %g at %d\n", pivot, i);
			free_data ();
			exit (-12);
		}
		for (j = i + 1; j < line_info [i].end; j++) {
			//if (line_info [j].start > i) continue;
			if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END)
				Kji = line_info [j].p [i];
			else if (sys_range_mode == RANGE_SYMM_UPPER)
				Kji = line_info [i].p [j];
			if (fabs (Kji) < tol) continue;
			// フルマトリックスとバンドマトリックスはピボットの直下から計算するが
			// 上半分は対角要素から計算するようにする
			if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END)
				for (k = i + 1; k < j; k++)
					line_info [j].p [k] -= line_info [i].p [k] * Kji / pivot;
			for (k = j; k < line_info [i].end; k++)
				line_info [j].p [k] -= line_info [i].p [k] * Kji / pivot;
			rhs_value [j] -= rhs_value [i] * Kji / pivot;
		}
	}
	//fputs ("end of forword deletion\n", stderr);
	// 後退代入
	for (i = size - 1; i >= 0; i--) {
		pivot = line_info [i].p [i];
		for (j = i + 1; j < line_info [i].end; j++) {
			Kji = line_info [i].p [j];
			if (fabs (Kji) < tol) continue;
			rhs_value [i] -= rhs_value [j] * Kji;
		}
		rhs_value [i] /= pivot;
		line_info [i].p [i] /= pivot;
		//fprintf (stderr, "%5d%15g%15g\n", i, rhs_value [i], pivot);
	}
}
