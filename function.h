//
//  関数の宣言
//  数が少ないものをピックアップ

// element_info.h
void init_einfo ();

// setup_total_K_mat.c
void create_K_info (FILE *fout);
void calc_K_mat (FILE *fout);

//
void calc_load ();

//
void solve_band_mat (double tol);

//
void calc_result ();
void write_result (FILE *fout)
