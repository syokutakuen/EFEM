//入出力ルーチンの宣言

int iJH_size (const Element *p);
void alloc_cinfo (Element *p, int ijh_size);
void attach_entities (int type, int ndata, int *data);

void io_data (FILE *fin, FILE *fout);
void io_bc (FILE *fin, FILE *fout);
