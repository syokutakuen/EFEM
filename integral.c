/* 
 * 線要素の積分点座標
 *  各形状ごと次数ごとに個別の配列としているが、
 *  引用するときは次数が決まっているので、
 *  この方が効率がいいだろうと判断したことによる。
 */
/* 順に1～3次の積分点位置の定義。 */
#define __LL1_it  0.577350269189626
#define __LL2_it  0.774596669241483
#define __LL3_it1 0.861136311594053
#define __LL3_it2 0.339981043584856

/* 積分点の位置の配列。各次数用に個別の配列として定義する。 */
double LineL0_intpt [] = { 0, };
double LineL1_intpt [] = { -__LL1_it, __LL1_it, };
double LineL2_intpt [] = { -__LL2_it, 0, __LL2_it, };
double LineL3_intpt [] = { -__LL3_it1, -__LL3_it2, __LL3_it2, __LL3_it1, };

/* 3次要素のための重みの値。無理数なのでマクロ定義する。 */
#define __LL3_wt1 0.347854845137454
#define __LL3_wt2 0.652145154862456

/* 0～3次の重みの配列。積分点位置同様、個別の配列とする。 */
double LineL0_weight [] = { 2 };
double LineL1_weight [] = { 1, 1, };
double LineL2_weight [] = { 5./9, 8./9, 5./9, };
double LineL3_weight [] = { __LL3_wt1, __LL3_wt2, __LL3_wt2, __LL3_wt1, };

/* 三角形要素の積分点と次数。三角形要素には0次のものはないので省略する。 */
double tria1_intpt [] = { 1./3., 1./3., };
double tria2_intpt [] = { 0.5, 0.5, 0, 0.5, 0.5, 0, };
double tria3_intpt [] = { 1./3., 1./3., 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 1, 0, 0, 1, 0, 0, };
double tria1_weight [] = { .5, };
double tria2_weight [] = { 1./6., 1./6., 1./6., };
double tria3_weight [] = { 9./40, 1./15, 1./15, 1./15, 1./40., 1./40., 1./40., };

/* 四角形要素の積分点位置と重み。線要素のデーターから作成した。 */
double quad0_intpt [] = { 0, 0, };
double quad1_intpt [] = {
	-__LL1_it, -__LL1_it, __LL1_it, -__LL1_it,
	-__LL1_it,  __LL1_it, __LL1_it,  __LL1_it,
};
double quad2_intpt [] = {
	-__LL2_it, -__LL2_it, 0, -__LL2_it, __LL2_it, -__LL2_it,
	-__LL2_it, 0,         0, 0,         __LL2_it, 0,
	-__LL2_it,  __LL2_it, 0,  __LL2_it, __LL2_it,  __LL2_it,
};
double quad3_intpt [] = {
	-__LL3_it1, -__LL3_it1, -__LL3_it2, -__LL3_it1, __LL3_it2, -__LL3_it1,  __LL3_it1, -__LL3_it1,
	-__LL3_it1, -__LL3_it2, -__LL3_it2, -__LL3_it2, __LL3_it2, -__LL3_it2,  __LL3_it1, -__LL3_it2,
	-__LL3_it1,  __LL3_it2, -__LL3_it2,  __LL3_it2, __LL3_it2,  __LL3_it2,  __LL3_it1,  __LL3_it2,
	-__LL3_it1,  __LL3_it1, -__LL3_it2,  __LL3_it1, __LL3_it2,  __LL3_it1,  __LL3_it1,  __LL3_it1,
};

double quad0_weight [] = { 4, };
double quad1_weight [] = { 1, 1, 1, 1, };
double quad2_weight [] = { 25./81, 40./81, 25./81, 40./81, 64./81, 40./81, 25./81, 40./81, 25./81, };

#define __LL3_wt11 (__LL3_wt1*__LL3_wt1)
#define __LL3_wt12 (__LL3_wt1*__LL3_wt2)
#define __LL3_wt21 __LL3_wt12
#define __LL3_wt22 (__LL3_wt2*__LL3_wt2)
double quad3_weight [] = {
	__LL3_wt11, __LL3_wt12, __LL3_wt12, __LL3_wt11,
	__LL3_wt21, __LL3_wt22, __LL3_wt22, __LL3_wt21,
	__LL3_wt21, __LL3_wt22, __LL3_wt22, __LL3_wt21,
	__LL3_wt11, __LL3_wt12, __LL3_wt12, __LL3_wt11,
};

//
// 春日屋の方法による積分点の定義
// 1~2次の積分点のみ定義する。
// 四角形1次の積分点と重み。4点則
// この公式にはもう一つ解ガウスの数値積分と同じである。
#define __QK1_it 0.8164965809277261
double Line1K_intpt [] = {
	0, -__QK1_it, __QK1_it, 0, 0, __QK1_it, -__QK1_it, 0,
};

// 積分点での重みはガウスのそれと同じである。
double *ine1K_weight = quad1_weight;

// 四角形二次の積分点と重み
// ガウスの数値積分点と異なり8点則である。
#derine __QK2_it1 0.6831300510639732
#derine __QK2_it2 0.8819171036881969
double Line2K_intpt [] = {
	 0,         -__QK2_it1, __QK2_it1,  0,          0,         __QK2_it1, -__QK2_it1, 0,
	-__QK2_it2, -__QK2_it2, __QK2_it2, -__QK2_it2, -__QK2_it2, __QK2_it2,  __QK2_it2, __QK2_it2,
};

#derine __QK2_wt1 0.8163265306122449
#derine __QK2_wt2 0.1836734693877551
double Line2K_intpt [] = {
	__QK2_wt1, __QK2_wt1, __QK2_wt1, __QK2_wt1,
	__QK2_wt2, __QK2_wt2, __QK2_wt2, __QK2_wt2,
};
