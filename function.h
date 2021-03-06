//
//  関数の宣言
//  数が少ないものをピックアップ

// element_info.h
void init_einfo ();

// setup_total_K_mat.c
void create_K_info (FILE *fout);
void calc_K_mat (FILE *fout);

// setup_mech_bc.c
void calc_load ();

// solve_Kmat.c
void solve_band_mat (double tol);

// calc_result.c
void calc_result ();
void write_result (FILE *fout)
