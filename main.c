//
//  EFEM
//    unEasy FEM program
//
//    written by syokutakuen
//

#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "datadef.h"
#include "data_io.h"
#include "material_mat.h"
#include "function.h "
#include "memfree.h"

int main (int argc, char **argv)
{
	FILE *fin, *fout;
	if (argc != 3) {
		fprintf (stderr, "usage: %s input output\n", argv [0]);
		exit(-1);
	}
	fprintf (stderr, "Start of analysis\n");
	init_einfo ();
	fprintf (stderr, "End of initializing\n");

	fin = fopen (argv [1], "r");
	if (fin == NULL) {
		fprintf (stderr, "cannot open input file: %s", argv [1]);
		return -1;
	}
	fout = fopen (argv [2], "w");
	if (fout == NULL) {
		fprintf (stderr, "cannot open output file: %s", argv [2]);
		return -1;
	}
	io_data (fin, fout);
	io_bc (fin, fout);
	fprintf (stderr, "End of reading data\n");

	calc_D_mat();
	sys_range_mode = RANGE_ALL;
	sys_range_mode = RANGE_BEGIN_END;
	sys_range_mode = RANGE_SYMM_UPPER;

	calc_K_mat (fout);
	fprintf (stderr, "End of construct stiffness matrix\n");
	solve_band_mat (tolerance);
	calc_result ();
	write_result (fout);
	fprintf (stderr, "End of solve probrem\n");

	free_data ();
	fprintf (stderr, "End of analysis\n");
	return 0;
}

