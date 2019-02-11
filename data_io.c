//  IOルーチン群

#include "memory.h"

static int nline = 0;
static char buffer [201];

//  文字列入力関数群
static void read_line (FILE *fin)
{
	char *p;
	fgets (buffer, sizeof (buffer), fin), nline++;
	while (*buffer == '#')
		fgets (buffer, sizeof (buffer), fin), nline++;
	for (p = buffer; !*p && *p != '#'; p++)
		;
	if (*p == '#')
		for (; !*p; p++) *p = 0;
	for (p = buffer; *p; p++)
		if (*p < 0x20) *p = 0;
}

static void read_str (FILE *fin, char **s)
{
	read_line (fin);
	*s = malloc (strlen(buffer)+1);
	strcpy (*s, buffer);
}

//
int iJH_size (const Element *p)
{
	return etinfo [prob].ijh_scale * p->info->info1->nnode + etinfo [prob].ijh_add;
}

void alloc_cinfo (Element *p, int ijh_size)
{
	int i;
	ElementCalcInfo *cp;
	int ntint = p->info->info1->ntint;
	p->cinfo = calloc (sizeof (ElementCalcInfo), ntint);
	if (p->cinfo == NULL) free_elem_data (p - elem);
	for (i = 0, cp = p->cinfo; i < ntint; i++, cp++) {
		cp->iJH = calloc (sizeof (double), ijh_size + 22);
		if (cp->iJH == NULL) free_elem_data (i);
		cp->stress = cp->iJH    + ijh_size;
		cp->strain = cp->stress + 6;
		cp->FT     = cp->strain + 6;
		cp->stat   = cp->FT     + 9;
	}
}

void attach_entities (int type, int ndata, int *data)
{
	int i, n, *ip;
	switch (type) {
	case SET_NODE:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_node (ip, *ip);
		break;
	case SET_SEG:
		break;
	case SET_ELEM:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_elem (ip, *ip);
		break;
	case SET_FACE:
		for (i = 0, ip = data; i < ndata; i++, ip++) {
			attach_elem (ip, *ip);
			n = *ip++;
			if ((*ip <= 0) || (*ip > elem [n].info->nface))
				error (-6, "Illegal edge/face number");
			*ip -= 1;
		}
		break;
	case SET_PART:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_part (ip, *ip);
		break;
	default:
		error (-6, "Illegal set type");
	}
}

//  データー入出力関数群
static void read_header (FILE *fin)
{
	read_str (fin, &title);
	read_line (fin);
	sscanf (buffer, "%d%d%d%d%d%d%d", &prob, &npart, &nmat, &ngeom, &ntnode, &ntelem, &nset);
	if (prob <= 0 && prob > PT_AXSOL) error (-2, "Illegal progrem setting.");
	if (npart < 0) error (-2, "Illegal part count.");
	if (nmat  < 0) error (-2, "Illegal material table count.");
	if (ngeom < 0) error (-2, "Illegal geometric table count.");
	if (ntnode < 0) error (-2, "Illegal node count");
	if (ntelem < 0) error (-2, "Illegal element count.");
	if (nset  < 0) error (-2, "Illegal set list count.");
}

static void write_header (FILE *fout)
{
	fprintf (fout, "\n\nEFEM/C ver1.06\n");
	fprintf (fout, "\tby H.Hoshino\n\n");
	fprintf (fout, "Analysis title\n");
	fprintf (fout, "\t%s\n\n", title);
	fprintf (fout, "Analysis type\n");
	fprintf (fout, "\t%s\n\n", etinfo [prob].name);
}

static void alloc_fe_data ()
{
	part = calloc (sizeof (PartData), npart);
	if (part == NULL) goto Error;
	mat  = calloc (sizeof (MatData *), nmat);
	if (mat == NULL)  goto Error;
	geom = calloc (sizeof (GeomData *), ngeom);
	if (geom == NULL) goto Error;
	node = calloc (sizeof (Node), ntnode);
	coord = node [0].coord = calloc (sizeof (double), ntnode * mdof);
	if (node == NULL) goto Error;
	elem = calloc (sizeof (Element), ntelem);
	if (elem == NULL) goto Error;
	set  = calloc (sizeof (SetData), nset);
	if (set == NULL) goto Error;
	return;
Error:
	free_mesh_data ();
	error  (-3, "No enough memory.");
}

static void read_mat (FILE *fin)
{
	int i;
	for (i = 0; i < nmat; i++) {
		MatData m;
		IsoMatData *mp;
		read_line (fin);
		sscanf (buffer, "%d%d%lg", &m.number, &m.type, &m.density);
		read_str (fin, &m.name);
		mp = malloc (sizeof (IsoMatData));
		*(MatData*)mp = m;
		read_line (fin);
		sscanf (buffer, "%lg%lg", &mp->young, &mp->poisson);
		mp->shar_mod = mp->young / 2 / (1 + mp->poisson);
		mat [i] = (MatData*)mp;
	}
}

static void write_mat (FILE *fout)
{
	int i;
	fprintf (fout, "Material table\n");
	fprintf (fout, "\tnumber of table: %d\n", nmat);
	for (i = 0; i < nmat; i++) {
		IsoMatData *p = (IsoMatData*)mat [i];
		fprintf (fout, "\n\t%-20s%d\n", "Table number: ", p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",           p->name);
		fprintf (fout, "\t\t%-20s%g\n", "Mass density:",   p->density);
		fprintf (fout, "\t\t%-20s%g\n", "Young modulas:",  p->young);
		fprintf (fout, "\t\t%-20s%g\n", "Poisson ratio:",  p->poisson);
	}
	fputc ('\n', fout);
}

static void read_geom (FILE *fin)
{
	int i;
	GeomData g, *gp;
	for (i = 0; i < nmat; i++) {
		read_line (fin);
		sscanf (buffer, "%d", &g.number);
		g.type = prob;
		read_str (fin, &g.name);
		if (prob == PT_PSS)
			gp = malloc (sizeof (PSS_GeomData));
		else
			gp = malloc (sizeof (GeomData));
		*gp = g;
		geom [i] = gp;
		if (prob == PT_PSS) {
			PSS_GeomData *pg = (PSS_GeomData *) gp;
			read_line (fin);
			sscanf (buffer, "%lg", &pg->thick);
		}
	}
}

static void write_geom (FILE *fout)
{
	int i;
	fprintf (fout, "Geometrical table\n");
	fprintf (fout, "\tnumber of table: %d\n", ngeom);
	for (i = 0; i < ngeom; i++) {
		GeomData *p = geom [i];
		fprintf (fout, "\n\t%-20s%d\n", "Table number:", p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",         p->name);
		fprintf (fout, "\t\t%-20s%s\n", "Type:",         etinfo [p->type].name);
		if (prob == PT_PSS) {
			PSS_GeomData *pg = (PSS_GeomData *) geom [i];
			fprintf (fout, "\t\t%-20s%g\n", "Thick:", pg->thick);
		}
	}
	fputc ('\n', fout);
}

static void read_part (FILE *fin)
{
	int i, j, m, g;
	PartData *p;
	for (i = 0; i < npart; i++) {
		p = part + i;
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &m, &g);
		read_str (fin, &p->name);
		for (j = 0; j < nmat; j++) {
			if (mat [j]->number == m) {
				p->mat = mat [i];
				break;
			}
		}
		if (j == nmat) error (-4, "Illegal material table number");
		for (j = 0; j < ngeom; j++) {
			if (geom [j]->number == m) {
				p->geom = geom [i];
				break;
			}
		}
		if (j == ngeom) error (-4, "Illegal geometrical table number");
	}
}

static void write_part (FILE *fout)
{
	int i;
	fprintf (fout, "Part information\n");
	fprintf (fout, "\tnumber of part: %d\n", npart);
	for (i = 0; i < npart; i++) {
		PartData *p = part + i;
		fprintf (fout, "\n\t%-20s%d\n", "Part Number ",       p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",              p->name);
		fprintf (fout, "\t\t%-20s%s\n", "Material table:",    p->mat->name);
		fprintf (fout, "\t\t%-20s%s\n", "Geometrical table:", p->geom->name);
	}
	fputc ('\n', fout);
}

static void read_node (FILE *fin)
{
	int i;
	double *dp;
	Node *p;
	for (i = 0, p = node, dp = coord; i < ntnode; i++, p++, dp += mdof) {
		p->coord = dp;
		read_line (fin);
		sscanf (buffer, "%d%lg%lg", &p->number, &p->coord [0], &p->coord[1]);
	}
}

static void write_node (FILE *fout)
{
	int i, j;
	double *dp;
	Node *p;
	char *sp [] = { "X-coord", "Y-coord", "Z-coord", };
	fprintf (fout, "%s\n" "\t%s : %d\n", "Node list", "number of node", ntnode);
	fprintf (fout, "%10s", "number");
	for (i = 0; i < mdim; i++)
		fprintf (fout, "%15s", sp [i]);
	for (i = 0, p = node; i < ntnode; i++, p++) {
		fprintf (fout, "%10d", p->number);
		for (j = 0, dp = p->coord; j < mdim; j++, dp++)
			fprintf (fout, "%15.6g", *dp);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
}

static void read_elem (FILE *fin)
{
	int i, j, ip, ie, ptr, pos;
	Element *p;
	char *cp;

	for (i = 0, p = elem; i < ntelem; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%d%n", &p->number, &ip, &ie, &pos);
		ptr = search_part (ip);
		if (ptr == npart) error (-5, "Illegal part table number");
		p->part = part + ptr;
		ptr = search_einfo (ie);
		if (ptr == neinfo) error (-5, "Illegal element infomation table number");
		p->info = einfo + ptr;
		p->conn = calloc (sizeof (int), p->info->info1->nnode);
		if (p->conn == NULL) free_elem_data (i);
		alloc_cinfo (p, iJH_size (p));
		cp = buffer + pos;
		for (j = 0; j < p->info->info1->nnode; j++) {
			int pos1, n;
			sscanf (cp, "%d%n", &n, &pos1);
			if (pos == pos1) {
				read_line (fin);
				cp = buffer;
				continue;
			}
			pos += pos1;
			cp += pos1;
			attach_node (p->conn + j, n);
		}
	}
}

static void write_elem (FILE *fout)
{
	int i, j, k, nnode, index;
	Element *p;
	fprintf (fout, "%s\n" "\t%s : %d\n\n", "Element list", "number of element", ntelem);
	fprintf (fout, "%10s%10s%10s %s\n", "number", "part", "type", "connectivity");
	for (i = 0, p = elem; i < ntelem; i++, p++) {
		fprintf (fout, "%10d%10d%10s", p->number, p->part->number, p->info->name);
		nnode = p->info->info1->nnode;
		if (nnode <= 8) {
			for (j = 0; j < nnode; j++)
				fprintf (fout, "%10d", node [p->conn [j]].number);
			fputc ('\n', fout);
		} else {
			for (j = index = 0; j < 8; j++, index++)
				fprintf (fout, "%10d", node [p->conn [index]].number);
			fputc ('\n', fout);
			nnode -= 8;
			for (j = 0; j < nnode / 10; j++) {
				fprintf (fout, "%10c", " ");
				for (k = 0; k < 10; k++, index++)
					fprintf (fout, "%10d", node [p->conn [index]].number);
				fputc ('\n', fout);
			}
			nnode %= 10;
			if (nnode > 0) {
				for (j = 0; j < nnode; j++, index++)
					fprintf (fout, "%10d", node [p->conn [index]].number);
				putc ('\n', fout);
			}
		}
	}
	fputc ('\n', fout);
}

static void read_int_list (FILE *fin, int ndata, int *ip)
{
	int pos, idata;
	char *cp;
	idata = 0;
	while (idata < ndata) {
		read_line (fin);
		cp = buffer;
		while (cp < buffer + strlen (buffer) && idata < ndata) {
			sscanf (cp, "%d%n", ip, &pos), cp += pos, idata++;
			//fprintf (stdout, "%5d%5d\n", *ip, pos);
			ip++;
		}
	}
}


static void read_set (FILE *fin, int nset, SetData *p)
{
	int i, scale, *ip;
	for (i = 0; i < nset; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &p->type, &p->ndata);
		if ((p->type <= SET_NONE) || (p->type > SET_PART))
			error (-6, "Illegal set type");
		scale = 1;
		if (p->type == SET_FACE) scale = 2;
		read_str (fin, &p->name);
		p->data = calloc (sizeof (int), p->ndata * scale);
		read_int_list (fin, p->ndata * scale, p->data);
		attach_entities (p->type, p->ndata, p->data);
	}
}

static void write_part_list (FILE *fout, int ndata, const int *ip)
{
	int i;
	for (i = 0; i < ndata; i++) {
		fprintf (fout, "%10d", part [*ip++].number);
		if (i % 8 == 0) fputc ('\n', fout);
	}
	if (i % 8) fputs ("\n", fout);
}

static void write_node_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++)
			fprintf (fout, "%10d", node [*ip++].number);
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++)
			fprintf (fout, "%10d", node [*ip++].number);
		fputs ("\n", fout);
	}
}

static void write_elem_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++)
			fprintf (fout, "%10d", elem [*ip++].number);
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++)
			fprintf (fout, "%10d", elem [*ip++].number);
		fputs ("\n", fout);
	}
}

static void write_face_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++) {
			fprintf (fout, "%10d:%d", elem [ip [0]].number, ip [1] + 1);
			ip += 2;
		}
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++) {
			fprintf (fout, "%10d:%d", elem [ip [0]].number, ip [1] + 1);
			ip += 2;
		}
		fputc ('\n', fout);
	}
}

static void write_set (FILE *fout, int nset, SetData *p)
{
	int i;
	fprintf (fout, "%s\n",          "Set data");
	fprintf (fout, "\t%s : %d\n\n", "number of set", nset);
	for (i = 0; i < nset; i++, p++) {
		fprintf (fout, "\t%20s%d\n", "Table number:", p->number);
		fprintf (fout, "\t%20s%s\n", "Name:",         p->name);
		fprintf (fout, "\t%20s\n",   "Data list");
		switch (p->type) {
		case SET_NODE:
			fputs ("Node list", fout);
			write_node_list (fout, p->ndata, p->data);
			break;
		case SET_SEG:
			break;
		case SET_ELEM:
			fputs ("Element list", fout);
			write_elem_list (fout, p->ndata, p->data);
			break;
		case SET_FACE:
			fputs ("Boundary edge/face list", fout);
			write_face_list (fout, p->ndata, p->data);
			break;
		case SET_PART:
			fputs ("Part list", fout);
			write_part_list (fout, p->ndata, p->data);
			break;
		default:
			error (-6, "Illegal set type");
		}
	}
	fputc ('\n', fout);
}

void io_data (FILE *fin, FILE *fout)
{
	read_header   (fin);
	write_header  (fout);
	alloc_fe_data ();

	read_mat   (fin);
	read_geom  (fin);
	read_part  (fin);
	write_part (fout);
	write_mat  (fout);
	write_geom (fout);

	read_node  (fin);
	write_node (fout);
	read_elem  (fin);
	write_elem (fout);

	read_set  (fin, nset, set);
	write_set (fout, nset, set);
	fflush (fout);
}

//  境界条件の入出力
static void read_bc_header (FILE *fin)
{
	read_line (fin);
	sscanf (buffer, "%d%d%d%d", &nfdisp, &npload, &nfload, &nbload);
	if (nfdisp < 0) error  (-7, "Illegal fixed displacement data count.");
	if (npload < 0) error  (-7, "Illegal point load data count.");
	if (nfload < 0) error  (-7, "Illegal face load data count.");
	if (nbload < 0) error  (-7, "Illegal body load data count.");
}

static void write_bc_header (FILE *fout)
{
	fprintf (fout, "%s\n\t%s: %d\n\t%s: %d\n\t%s: %d\n\t%s: %d\n\n",
		"Boundary condition tables",
		"Count of fixed disp table", nfdisp,
		"Count of point load table", npload,
		"Count of distributed load table", nfload,
		"Count of volume load table", nbload);
}

static void alloc_bc_data ()
{
	fdisp = calloc (sizeof (FixDisp), nfdisp);
	if (fdisp == NULL) goto Error;
	pload = calloc (sizeof (LoadData), npload);
	if (pload == NULL) goto Error;
	fload = calloc (sizeof (LoadData), nfload);
	if (fload == NULL) goto Error;
	bload = calloc (sizeof (LoadData), nbload);
	if (bload == NULL) goto Error;
	return;
Error:
	free_bc_data ();
	error (-7, "no enough memory");
}

static void read_fdisp (FILE *fin)
{
	int i, j, pos, *ip;
	char *cp;
	double *dp;
	FixDisp *p;

	for (i = 0, p = fdisp; i < nfdisp; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%n", &p->number, &p->ndata, &pos);
		for (j = 0, ip = p->flags, cp = buffer + pos; j < etinfo [prob].ndof; j++, ip++) {
			sscanf (cp, "%d%n", ip, &pos);
			cp += pos;
		}
		read_str (fin, &p->name);
		read_line (fin);
		for (j = 0, dp = p->value, cp = buffer; j < etinfo [prob].ndof; j++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
		sscanf (buffer, "%lg%lg", &p->value [0], &p->value [1]);
		p->data = calloc (sizeof (int), p->ndata);
		read_int_list (fin, p->ndata, p->data);
		p->set_type = SET_NODE;
		for (j = 0; j < p->ndata; j++)
			attach_node (p->data + j, p->data [j]);
	}
}

static void set_load_value (FILE *fin, int nload, LoadData *p)
{
	char *cp;
	double *dp;
	int i, pos;

	read_line (fin);
	if ((p->val_type % 10== BC_VAL) || (p->val_type % 10 == BC_NORMAL && p->load_type % 10 == BC_FACE)) {
		for (i = 0, dp = p->value, cp = buffer; i < etinfo [prob].ndim; i++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
	} else if (p->val_type % 10 <= BC_VAL_VEC) {
		for (i = 0, dp = p->value, cp = buffer; i < etinfo [prob].ndim + 1; i++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
	} else {
		free_bc_data ();
		error (-8, "illegal point load type");
	}
}

static void read_load (FILE *fin, int nload, LoadData *p, int ltype)
{
	int i, j, ndata, scale, *ip;

	for (i = 0; i < nload; i++, p++) {
		memset (p, 0, sizeof (LoadData));
		p->load_type = ltype;
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &ndata, &p->val_type);
		read_str (fin, &p->name);
		set_load_value (fin, i, p);
		p->ndata = ndata;
		p->data = calloc (sizeof (int), ndata);
		scale = 1;
		if (ltype == BC_FACE) scale = 2;
		read_int_list (fin, ndata * scale, p->data);
		attach_entities (ltype, p->ndata, p->data);
	}
}

static void read_bc (FILE *fin)
{
	read_fdisp (fin);
	read_load (fin, npload, pload, BC_POINT);
	read_load (fin, nfload, fload, BC_FACE);
	read_load (fin, nbload, bload, BC_BODY);
}

static void write_fdisp (FILE *fout)
{
	int i, j;
	FixDisp *p;
	char *sp[] = { "X-dir", "Y-dir", "Z-dir", "RX-dir", "RY-dir", "RZ-dir", };
	for (i = 0, p = fdisp; i < nfdisp; i++, p++) {
		fprintf (fout, "\n%20s%d\n%20s%s\n", "Table number:", p->number, "Name:", p->name);
		fprintf (fout, "%20s\n", "Constraint dof");
		for (j = 0; j <  etinfo [prob].ndof; j++)
			if (p->flags [j]) fprintf (fout, "%28s:%g\n", sp [j], p->value[j]);
		fprintf (fout, "Apply nodes\n");
		write_node_list (fout, p->ndata, p->data);
		fputc ('\n', fout);
	}
}

typedef void (*fp_write_int_list) (FILE *fout, int, const int *);

static void write_load (FILE *fout, int ndata, const LoadData *p, fp_write_int_list fp)
{
	int i, j;
	char *type1_sp  [] = { "None", "Direction", "Director and value", "Director mult value", };
	char *type2a_sp [] = { "None", "par length", "par area", };
	char *type2b_sp [] = { "None", "gravity", "uniform load", };
	char **type2_sp    = (p->set_type == BC_BODY)? type2b_sp: type2a_sp;
	char *type3_sp  [] = { "None", "at node", "at segment", "at face", "at body", "at part", };

	for (i = 0; i < ndata; i++, p++) {
		fprintf (fout, "%20s%d\n" "%20s%s\n", "Table number:", p->number, "Name:", p->name);
		fprintf (fout, "%20s:%s %s %s\n",
				"Load types", type1_sp [p->val_type % 10],
				type2_sp [p->val_type / 10], type3_sp [p->set_type]);
		fprintf (fout, "%20s", "Load value and vector:");
		for (j = 0; j <  etinfo [prob].ndof; j++)
			fprintf (fout, "%15g", p->value[j]);
		fputc ('\n', fout);
		fprintf (fout, "Apply entitis\n");
		(*fp) (fout, p->ndata, p->data);
		fputc ('\n', fout);
	}
}

static void write_bc (FILE *fout)
{
	fprintf (fout, "%s\n\t%s: %d\n\n", "Fixed displacement data", "Number of tables", nfdisp);
	write_fdisp (fout);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Point load data", "Number of tables", npload);
	write_load (fout, npload, pload, write_node_list);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Distributed load data", "Number of tables", nfload);
	write_load (fout, nfload, fload, write_face_list);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Volume load data", "Number of tables", nbload);
	write_load (fout, nbload, bload, write_elem_list);
}

void io_bc (FILE *fin, FILE *fout)
{
	read_bc_header (fin);
	write_bc_header (fout);
	alloc_bc_data ();
	fflush (fout);
	read_bc (fin);
	write_bc (fout);
	fflush (fout);
}



