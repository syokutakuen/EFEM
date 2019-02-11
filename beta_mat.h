//
//
//

double invJ_2D (double *iJ, const double *J);
double invJ_3D (double *iJ, const double *J);

//  ヤコビアンの計算 形状関数の微分形Hに要素の節点座標値をかけることで得られる。
//  修正予定
void calc_Jmat (double *J, const double *H, const Element *e);

//  二次元配列の掛け算
void mat_mul2 (double *A, const double *B, const double *C, int sz1, int sz2, int sz3);


//	2次元要素のAマトリックスを作る。
void calc_PL_iJH (const Element *e, const double *N, const double *H, int nint);

//  2次元力学解析用要素のBマトリックスの作成
void calc_PL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
void calc_AXSOL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
