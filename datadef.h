//
//	大域変数の定義
//

#if !defined (__datadef_H__)
#define __datadef_H__

//	大域空間での座標系の数と自由度の数
EXTERN int mdof;
EXTERN int mdim;

//	解析タイトル
EXTERN char *title;

//	汎用許容値
EXTERN double tolerance = 1e-10;

//	問題番号。 PT_*で始まる列挙子の何れかがが入る。
EXTERN int prob;

//	パーツや物性値、幾何特性値の配列
EXTERN int npart, nmat, ngeom;
EXTERN PartData *part;
EXTERN MatData **mat;
EXTERN GeomData **geom;

//	節点と要素の配列の宣言
EXTERN int ntnode, ntelem;
EXTERN Node *node;
EXTERN Element *elem;
EXTERN double *coord;

//	セットデーターの宣言
EXTERN int nset;
EXTERN SetData *set;

//	境界条件データーの宣言
EXTERN int nfdisp, npload, nfload, nbload;
EXTERN FixDisp *fdisp;
EXTERN LoadData *pload, *fload, *bload;

//	全体剛性マトリックスへの数値の格納方法
//	PROB_*で始まる列挙子のいずれかが入る。
EXTERN int sys_range_mode;

//	ライン情報構造体の変数
EXTERN K_line_info *line_info;

//	全体剛性の実データーと変位と荷重
EXTERN double *sysk, *lhs_value, *rhs_value;

//	計算結果を格納
EXTERN double *disp, *force;

//	節点番号と自由度から全体剛性のインデックスを逆に引くリスト
EXTERN int *rev_line_info, *decode_rev_info;
//	各境界条件の開始位置と終了位置
EXTERN int rank_line_info [4][2];

//	節点が参照する要素数と各自由度のタイプ
EXTERN int *nodal_dof_list;

//	要素情報テーブルの宣言
EXTERN ElementTypeInfo etinfo [];
EXTERN ElemIntegralInfo eiinfo [];
EXTERN ElementInfo einfo [];
EXTERN int netinfo;
EXTERN int neiinfo;
EXTERN int neinfo;

#else
#endif
