/*
 * メモリの解放
 */

#include <stdlib.h>
#include <memory.h>
#include "typedef.h"
#include "datadef.h"

/* indexまでの要素内データーを全削除し、要素リストを削除する。 */
void free_elem_data (int index)
{
	int i, j;
	Element *p;
	ElementCalcInfo *cp;
	for (i = 0, p = elem; i <= index; i++, p++) {
		int ntint = p->info->info1->ntint;
		free (p->conn);
		for (j = 0, cp = p->cinfo; j < ntint; j++, cp++)
			free (cp->iJH);
		free (p->cinfo);
	}
}

void free_set_data (SetData *p)
{
	if (p == NULL) return;
	free (p->name);
	free (p->data);
}

void free_mesh_data ()
{
	int i;
	free (title); title = NULL;
	free (part); part = NULL;
	for (i = 0; i < nmat; i++) {
		if (mat [i] != NULL)
			free (mat [i]->name);
		free (mat [i]);
	}
	free (mat);
	nmat = 0; mat = NULL;
	for (i = 0; i < ngeom; i++) {
		if (geom [i] != NULL)
			free (geom [i]->name);
		free (geom [i]);
	}
	free (geom);
	ngeom = 0; geom = NULL;
	free (coord);
	free (node);
	ntnode = 0; node = NULL;
	free_elem_data (ntelem - 1);
	free (elem);
	ntelem = 0; elem = NULL;
	for (i = 0; i < nset; i++)
		free_set_data (set + i);
	free (set);
	nset = 0; set = NULL;
}

void free_load (int ndata, LoadData *p)
{
	int i;
	for (i = 0; i < ndata; i++)
		free_set_data ((SetData *)(p + i));
}

void free_bc_data ()
{
	int i;
	for (i = 0; i < nfdisp; i++)
		free_set_data ((SetData *)(fdisp + i));
	free (fdisp);
	nfdisp = 0; fdisp = NULL;
	free_load (npload, pload);
	free (pload);
	npload = 0; pload = NULL;
	free_load (nfload, fload);
	free (fload);
	nfload = 0; fload = NULL;
	free_load (nbload, bload);
	free (bload);
	nbload = 0; bload = NULL;
}

void free_sys_info ()
{
	free (line_info);
	free (nodal_dof_list);
	free (rev_line_info);
	free (decode_rev_info);
	free (sysk);
	free (rhs_value);
	free (lhs_value);
	free (disp);
	free (force);
}

void free_data ()
{
	free_mesh_data ();
	free_bc_data ();
	free_sys_info ();
}
