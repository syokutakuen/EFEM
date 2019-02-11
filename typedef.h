//
//	有限要素解析に必要なデーター型の定義
//

//多重定義の禁止
#if !defined (__TYPE_DEF_H__)

//	物性値の型
enum {
	MatNone = 0. 
	MatMechIso,
	MatMechOrtho,
	MatMechAniso,
	MatHeathIso,
	MatHeathOrtho,
	MathHeatAniso,
};

//	物性値データー。基本型のほか、等方体構造解析のためのデーターを定義する。
typedef struct MatData *pMatData;
typedef struct {
	int number;			//	一連番号
	char *name;			//	名前
	int type;			//	この物性値の型。今のところ未使用
	pMatData ref_mat;	//	他のデーターを参照するためのフィールド
	double ref_temp;	//	参照温度。解析はこの温度を基準とする。
	double density;		//	質量密度
} MatData;

//	等方弾性体の定義
//	MatDataにデーターを追加する形で定義する。
typedef struct {
	int number;
	char *name;
	int type;
	pMatData ref_mat;
	double ref_temp;
	double density;
	double young, poisson, shar_mod;
} MechIsoMatData;

//	等方性熱物性
typedef struct {
	int number;
	char *name;
	int type;
	pMatData ref_mat;
	double ref_temp;
	double density;
	double heat_ratio;
	double heat_trans_ratio;
} HeatIsoMatData;

// 幾何特性値
typedef struct {
	int number;
	char *name;
	int type;
} GeomData;

// パーツデーター。物性値と幾何特性値を持つ。要素から参照される。
typedef struct {
	int number;
	char *name;
	MatData *mat;
	GeomData *geom;
	double *init_mat; // 初期物性マトリックス
} PartData;

//	節点データー
//	座標値がポインタなのは、ヤコビアンの計算ルーチンが他でも転用できることが判明したため。
typedef struct {
	int number;
	double *coord;
} Node;

/* 要素タイプ構造体 */
typedef struct {
	int number;
	char *name;
	int prob, ndim, ndof, nvstat, ijh_scale, ijh_add;
} ElementTypeInfo;

/* 形状関数のポインタ */
typedef void (*SF)(double *, const double *);

/* 要素積分計算情報構造体の定義 */
struct ElementInfo;
typedef struct {
	int number;
	int nnode, norder, ntint, ndim; // 要素あたりの節点数、次数、全積分点数、次元数
	SF sf, dsf;   // 形状関数へのポインター
	double *N, *dN, *ipt, *wt; // 形状関数と形状関数の微分形、積分点情報と重みの配列
} ElemIntegralInfo;

/* 要素情報構造体の定義。各要素から参照される。
 * 計算情報が二つあるが、要素全体と各表面の形状関数を含むためである。
 * 要素表面の計算は分布荷重の計算に用いる。
 */
typedef struct {
	int number;
	char *name;
	int nvert, nface;  // 要素あたりの節点数、頂点(隅節点)の数、表面の数
	ElemIntegralInfo *info1, *info2, *info3; // 要素全体と表面の計算情報
} ElementInfo;

/* 計算情報データー。構造解析に特化したデーターの組み方をしているが、それ以外にも使用可能である。 */
typedef struct {
	double detJ, sigy, *iJH;  // ヤコビアンの行列値とJ-1*Hの結果。
	double *stress, *strain, *FT; // 応力とひずみ、変形勾配
	double *stat;   // 状態変数(温度など)
} ElementCalcInfo;

/* 要素データー */
typedef struct {
	int number;
	PartData *part;
	ElementInfo *info;
	int *conn;
	ElementCalcInfo *cinfo;
} Element;

/* セットデーター */
/* 最初にセットに含まれるデーター型を定義する。
 * このようなものは数値そのものを使用せず、シンボル定数を使用する。
 */
enum { SET_NONE=0, SET_NODE, SET_SEG, SET_ELEM, SET_FACE, SET_PART, };
typedef struct {
	int number;
	char *name;
	int type;
	int ndata, *data;
} SetData;

/* 
 * 境界条件の定義
 */
/* 拘束条件の定義 */
typedef struct {
	int number;
	char *name;
	int flags [6];		// 拘束条件をかける場合、1をセットする。
	double value [6];	// 拘束条件の値。
	SetData *set;		// セットを使用する。この場合名前を使わない。
} FixDisp;

/* 荷重条件の定義 */
enum { BC_NONE, BC_VAL, BC_VEC, BC_VAL_VEC, };	// 無効(0)、各成分の値、方向と大きさ（方向は正規化される)、 方向と大きさの積
enum { BC_LENGTH, BC_AREA = 10, };		// 単位長さ(線分のみ)と単位面積ごとの入力。
						// 線分に荷重をかける要素で有効。ソリッドシェル・固体要素では無視される。
enum { BC_POINT=1, BC_FACE, BC_BODY, };		// 荷重をかける対象。節点・表面・体積(要素全体)

typedef struct {
	int number;
	char *name;
	int load_type;
	int val_type;	// 上の3つの列挙体要素の格納先
	double value [4], load [3];	// 入力値と実際の荷重
	SetData *set;			// セット。名前は使わない。
} LoadData;

#else
#define __TYPE_DEF_H__
#endif
