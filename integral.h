//
//
//

#if !defined (__INTEGRAL_H__)
#define __INTEGRAL_H__

//  積分点の位置の配列。各次数用に個別の配列として定義する。
EXTERN double LineL0_intpt [];
EXTERN double LineL1_intpt [];
EXTERN double LineL2_intpt [];
EXTERN double LineL3_intpt [];

//  0～3次の重みの配列。積分点位置同様、個別の配列とする。
EXTERN double LineL0_weight [];
EXTERN double LineL1_weight [];
EXTERN double LineL2_weight [];
EXTERN double LineL3_weight [];

//  三角形要素の積分点と次数。三角形要素には0次のものはないので省略する。
EXTERN double tria1_intpt [];
EXTERN double tria2_intpt [];
EXTERN double tria3_intpt [];
EXTERN double tria1_weight [];
EXTERN double tria2_weight [];
EXTERN double tria3_weight [];

//  四角形要素の積分点位置と重み。線要素のデーターから作成した。
EXTERN double quad0_intpt [];
EXTERN double quad1_intpt [];
EXTERN double quad2_intpt [];
EXTERN double quad3_intpt [];
EXTERN double quad0_weight [];
EXTERN double quad1_weight [];
EXTERN double quad2_weight [];
EXTERN double quad3_weight [];

// 春日屋の方法による積分点の定義
EXTERN double Line1K_intpt [];
EXTERN double *ine1K_weight;
EXTERN double Line2K_intpt [];
EXTERN double Line2K_intpt [];

#endif
