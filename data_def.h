//
//	有限要素解析に必要なデーター型の定義
//

//多重定義の禁止
#if !defined (__DATA_DEF_H__)

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

#else
#define __DATA_DEF_H__
#endif
