//
// 1～2次元形状関数の定義。
//

//
// ラグランジュ補完による線要素の形状関数
// アンダースコアがついているのは内部で使用するため
//

static void __line1L (double *N, double t)
{
	N [0] = (1-t)/2;
	N [1] = (1+t)/2;
}

static void __line2L (double *N, double t)
{
	N [0] = t*(t-1)/2;
	N [1] = 1-t*t;
	N [2] = t*(t+1)/2;
}

static void __line3L (double *N, double t)
{
	N [0] = -(1-t)*(1-9*t*t)/16;
	N [1] = 9*(1-3*t)*(1-t*t)/16;
	N [2] = 9*(1+3*t)*(1-t*t)/16;
	N [3] = -(1+t)*(1-9*t*t)/16;
}

static void __dline1L (double *N, double t)
{
	N [0] = -0.5;
	N [1] =  0.5;
}

static void __dline2L (double *N, double t)
{
	N [0] = (2*t-1)/2;
	N [1] = -2*t;
	N [2] = (2*t+1)/2;
}

static void __dline3L (double *N, double t)
{
	N [0] = (1+18*t-27*t*t)/16;
	N [1] = 9*(-3-2*t+9*t*t)/16;
	N [2] = 9*(3-2*t-9*t*t)/16;
	N [3] = (-1+18*t+27*t*t)/16;
}

//	外部から呼び出すためのラッパー
void line1L (double *N, const double *t)  { __line1L (N, *t); }
void line2L (double *N, const double *t)  { __line2L (N, *t); }
void line3L (double *N, const double *t)  { __line3L (N, *t); }
void dline1L (double *N, const double *t) { __dline1L (N, *t); }
void dline2L (double *N, const double *t) { __dline2L (N, *t); }
void dline3L (double *N, const double *t) { __dline3L (N, *t); }

//二次元要素の形状関数
//三角形要素
void tri1 (double *N, const double *t)
{
	N [0] = t [0];
	N [1] = t [1];
	N [2] = 1 - t [0] - t [1];
}

void tri2 (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [0] = t1 * (2 * t1 - 1);
	N [1] = t2 * (2 * t2 - 1);
	N [2] = t3 * (2 * t3 - 1);
	N [3] = 4 * t1 * t2;
	N [4] = 4 * t2 * t3;
	N [5] = 4 * t3 * t1;
}

void tri3L (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [0] = t1 * (3 * t1 - 1) * (3 * t1 - 2) / 2;
	N [1] = t2 * (3 * t2 - 1) * (3 * t2 - 2) / 2;
	N [2] = t3 * (3 * t3 - 1) * (3 * t3 - 2) / 2;
	N [3] = 9 * t1 * t2 * (3 *t1 -1) / 2;
	N [4] = 9 * t1 * t2 * (3 *t2 -1) / 2;
	N [5] = 9 * t2 * t3 * (3 *t2 -1) / 2;
	N [6] = 9 * t2 * t3 * (3 *t3 -1) / 2;
	N [7] = 9 * t3 * t1 * (3 *t3 -1) / 2;
	N [8] = 9 * t3 * t1 * (3 *t1 -1) / 2;
	N [9] = 27 * t1 * t2 * t3;
}

// 中心節点が無い三角形三次要素
//この種の要素は例えば「有限要素法全解」（1990、パーソナルメディア）にあるが、精度が悪いため開発した。
//完全三次要素の多項式からt1t2の項を覗くことで得られる。
void tri3S (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [0] = t1 * (3 * t1 - 1) * (3 * t1 - 2) / 2;
	N [1] = t2 * (3 * t2 - 1) * (3 * t2 - 2) / 2;
	N [2] = t3 * (3 * t3 - 1) * (3 * t3 - 2) / 2;
	N [3] = 9 * t1 * t2 * (2 * t1 - t2) / 2;
	N [4] = 9 * t1 * t2 * (2 * t2 - t1) / 2;
	N [5] = 9 * t2 * t3 * (2 * t2 - t3) / 2;
	N [6] = 9 * t2 * t3 * (2 * t3 - t2) / 2;
	N [7] = 9 * t3 * t1 * (2 * t3 - t1) / 2;
	N [8] = 9 * t3 * t1 * (2 * t1 - t3) / 2;
}

void dtri1 (double *N, const double *t)
{
	N [0] = 1;
	N [1] = 0;
	N [2] = -1;
	N [3] = 0;
	N [4] = 1;
	N [5] = -1;

}

void dtri2 (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [ 0] = 4 * t1 - 1;
	N [ 1] = 0;
	N [ 2] = -4 * t3 + 1;
	N [ 3] = 4 * t2;
	N [ 4] = -4 * t2;
	N [ 5] = 4 * (t3 - t1);

	N [ 6] = 0;
	N [ 7] = 4 * t2 - 1;
	N [ 8] = -4 * t3 + 1;
	N [ 9] = 4 * t1;
	N [10] = 4 * (t3 - t2);
	N [11] = -4 * t1;
}

void dtri3L (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [ 0] = (27 * t1 * t1 - 18 * t1 + 2) / 2;
	N [ 1] = 0;
	N [ 2] = -(27 * t3 * t3 - 18 * t3 + 2) / 2;
	N [ 3] = 9 * t2 * (6 *t1 -1) / 2;
	N [ 4] = 9 * t2 * (3 *t2 -1) / 2;
	N [ 5] = -9 * t2 * (3 *t2 -1) / 2;
	N [ 6] = -9 * t2 * (6 *t3 -1) / 2;
	N [ 7] = 9 * (3 * t3 * t3 - 6 * t1 * t3 + t1 - t3) / 2;
	N [ 8] = 9 * (6 * t3 * t3 - 3 * t1 * t1 + t1 - t3) / 2;
	N [0] = L1 [0] * L2 [0];
	N [ 9] = 27 * t2 * (t3 - t1);

	N [10] = 0;
	N [11] = (27 * t2 * t2 - 18 * t2 + 2) / 2;
	N [12] = -(27 * t3 * t3 - 18 * t3 + 2) / 2;
	N [13] = 9 * t1 * (3 *t1 -1) / 2;
	N [14] = 9 * t1 * (6 *t2 -1) / 2;
	N [15] = 9 * (6 * t2 * t3 - 3 * t2 * t2 + t2 - t3) / 2;
	N [16] = 9 * (3 * t3 * t3 - 6 * t2 * t3 + t2 - t3) / 2;
	N [17] = -9 * t1 * (6 *t3 -1) / 2;
	N [18] = -9 * t1 * (3 *t1 -1) / 2;
	N [19] = 27 * t1 * (t3 - t2);
}

void tri3S (double *N, const double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [ 0] =  (27 * t1 * t1 - 18 * t1 + 2) / 2;
	N [ 1] = 0;
	N [ 2] = -(27 * t3 * t3 - 18 * t3 + 2) / 2;
	N [ 3] =   9 * t2 * (4 * t1 - t2) / 2;
	N [ 4] =  18 * t2 * (t2 - t1) / 2;
	N [ 5] = -18 * t2 * (t2 - t3) / 2;
	N [ 6] =  -9 * t2 * (4 * t3 - t2) / 2;
	N [ 7] =   9 * (2 * t3 * t3 + t1 * t1 - 6 * t3 * t1) / 2;
	N [ 8] =  -9 * (t3 * t3 + 2 * t1 * t1 - 6 * t3 * t1) / 2;

	N [ 9] = 0;
	N [10] =  (27 * t2 * t2 - 18 * t2 + 2) / 2;
	N [11] = -(27 * t3 * t3 - 18 * t3 + 2) / 2;
	N [12] =  18 * t1 * (t1 - t2) / 2;
	N [13] =   9 * t1 * (4 * t2 - t1) / 2;
	N [14] =  -9 * (t3 * t3 + 2 * t2 * t2 - 6 * t2 * t3) / 2;
	N [15] =   9 * (2 * t3 * t3 + t2 * t2 - 6 * t2 * t3) / 2;
	N [16] =  -9 * t1 * (4 * t3 - t1) / 2;
	N [17] = -18 * t1 * (t1 - t3) / 2;
}

//	四角形ラグランジュ要素の実装
//1次要素
//コネクティビティ
// 4---3
// |   |
// 1---2
void __quad1 (double *N, const double *t, const double *L1, const double *L2)
{
	N [0] = L1 [0] * L2 [0];
	N [1] = L1 [1] * L2 [0];
	N [2] = L1 [1] * L2 [1];
	N [3] = L1 [0] * L2 [1];
}

void quad1 (double *N, const double *t)
{
	double L1 [2], L2 [2];
	__line1L (L1, t [0]);
	__line1L (L2, t [1]);
	__quad1L (N, t, L1, L2);
}

//四角形要素の微分係数を求める
//	一々式を書くのは鬱陶しいので、既存の関数を使いまわすことにした。
//	こうして同じようなことを延々と行うことがなくなる。
void dquad1 (double *N, const double *t)
{
	double L1 [2], L2 [2];
	__dline1L (L1, t [0]);
	__line1L  (L2, t [1]);
	__quad1L (N, t, L1, L2);

	__line1L  (L1, t [0]);
	__dline1L (L2, t [1]);
	__quad1L (N + 4, t, L1, L2);
}

//2次要素
//コネクティビティ
// 4--7--3
// |     |
// 8  9  6
// |     |
// 1--5--2
void __quad2L (double *N, const double *t, const double *L1, const double *L2)
{
	N [0] = L1 [0] * L2 [0];
	N [1] = L1 [2] * L2 [0];
	N [2] = L1 [2] * L2 [2];
	N [3] = L1 [0] * L2 [2];

	N [4] = L1 [1] * L2 [0];
	N [5] = L1 [2] * L2 [1];
	N [6] = L1 [1] * L2 [2];
	N [7] = L1 [0] * L2 [1];
	N [8] = L1 [1] * L2 [1];
}

void quad2L (double *N, const double *t)
{
	double L1 [3], L2 [3];
	__line2L (L1, t [0]);
	__line2L (L2, t [1]);
	__quad2L (N, t, L1, L2);
}

void dquad2L (double *N, const double *t)
{
	double L1 [3], L2 [3];
	__dline2L (L1, t [0]);
	__line2L  (L2, t [1]);
	__quad2L (N, t, L1, L2);

	__line2L  (L1, t [0]);
	__dline2L (L2, t [1]);
	__quad2L (N + 9, t, L1, L2);
}

//3次要素
//コネクティビティ
//  4--10---9--3
//  |          |
// 11  16  15  8
//  |          |
// 12  13  14  7
//  |          |
//  1---5---6--2
void __quad3L (double *N, const double *t, const double *L1, const double *L2)
{
	N [ 0] = L1 [0] * L2 [0];
	N [ 1] = L1 [3] * L2 [0];
	N [ 2] = L1 [3] * L2 [3];
	N [ 3] = L1 [0] * L2 [3];

	N [ 4] = L1 [1] * L2 [0];
	N [ 5] = L1 [2] * L2 [0];
	N [ 6] = L1 [3] * L2 [2];
	N [ 7] = L1 [3] * L2 [1];

	N [ 8] = L1 [0] * L2 [1];
	N [ 9] = L1 [0] * L2 [2];
	N [10] = L1 [3] * L2 [2];
	N [11] = L1 [3] * L2 [1];

	N [12] = L1 [1] * L2 [1];
	N [13] = L1 [2] * L2 [1];
	N [14] = L1 [2] * L2 [2];
	N [15] = L1 [1] * L2 [2];
}

void quad3L (double *N, const double *t)
{
	double L1 [4], L2 [4];
	__line3L (L1, t [0]);
	__line3L (L2, t [1]);
	__quad3L (N, t, L1, L2);
}

void dquad3 (double *N, const double *t)
{
	double L1 [2], L2 [2];
	__dline3L (L1, t [0]);
	__line3L  (L2, t [1]);
	__quad3L (N, t, L1, L2);

	__line3L  (L1, t [0]);
	__dline3L (L2, t [1]);
	__quad3L (N + 16, t, L1, L2);
}

//	セレンディピティ要素の実装
//	具体的な内容はO.C.Zienkiewiczの定本による
//	この方法は合理的だが、使用者がほとんど見当たらないのが謎。
//	一次要素はラグランジュ要素と同じなので割愛。

void __quad2S (double *N, const double *t, const double *L1, const double *L2, double L12, double L22)
{
	// 中間節点の形状関数。
	N [4] = L12    * L2 [0];
	N [5] = L1 [1] * L22;
	N [6] = L12    * L2 [1];
	N [7] = L1 [0] * L22;

	// 隅節点の形状関数を定義
	// 初めに一次要素の隅節点の形状関数を求める
	double NL [4];
	__quad1L (NL, t, L1, L2);
	N [0] = NL [0] - (N [4] + N [7]) / 2;
	N [1] = NL [1] - (N [5] + N [4]) / 2;
	N [2] = NL [2] - (N [6] + N [5]) / 2;
	N [3] = NL [3] - (N [7] + N [6]) / 2;
}

void quad2S (double *N, const double *t)
{
	double t1 = t [0], t2 = t [1];
	double L1 [2], L2 [2];
	double L21 [3], L22 [3];
	__line1L (L1, t1);
	__line1L (L2, t2);
	__line2L (L21, t1);
	__line2L (L22, t2);
	__quad2S (N, t, L1, L2, L21 [1], L22 [1]);
}

void quad2S (double *N, const double *t)
{
	double t1 = t [0], t2 = t [1];
	double L1 [2], L2 [2], NL [4];
	double L21 [3], L22 [3];

	__dline1L (L1, t1);
	__line1L  (L2, t2);
	__dline2L (L21, t1);
	__line2L  (L22, t2);
	__quad2S  (N, t, L1, L2, L21 [1], L22 [1]);

	__line1L  (L1, t1);
	__dline1L (L2, t2);
	__line2L  (L21, t1);
	__dline2L (L22, t2);
	__quad2S  (N + 8, t, L1, L2, L21 [1], L22 [1]);
}

void __quad3S (double *N, const double *t, const double *L1, const double *L2, const double *L13, const double *L23)
{
	// 中間節点の形状関数。
	N [ 4] = L13 [1] * L2  [0];
	N [ 5] = L13 [2] * L2  [0];
	N [ 6] = L1 [1]  * L23 [2];
	N [ 7] = L1 [1]  * L23 [1];
	N [ 8] = L13 [2] * L2  [1];
	N [ 9] = L13 [1] * L2  [1];
	N [10] = L1 [0]  * L23 [1];
	N [11] = L1 [0]  * L23 [2];

	// 隅節点の形状関数を定義
	// 初めに一次要素の隅節点の形状関数を求める
	double NL [4];
	__quad1L (NL, t, L11, L21);
	N [0] = NL [0] - (2 * (N [ 4] + N [11]) + (N [ 5] + N [10])) / 3;
	N [1] = NL [1] - (2 * (N [ 6] + N [ 5]) + (N [ 7] + N [ 4])) / 3;
	N [2] = NL [2] - (2 * (N [ 8] + N [ 7]) + (N [ 9] + N [ 6])) / 3;
	N [3] = NL [3] - (2 * (N [10] + N [ 9]) + (N [11] + N [ 8])) / 3;
}

void quad3S (double *N, const double *t)
{
	double L1 [2], L2 [2], L13 [4], L23 [4];
	__line1L (L1, t [0]);
	__line1L (L2, t [1]);
	__line3L (L13, t);
	__line3L (L23, t);
	__quad3S (N, t, L1, L2, L13, L23);
}

void dquad3S (double *N, const double *t)
{
	double L1 [2], L2 [2], L13 [4], L23 [4];
	__dline1L (L1, t [0]);
	__line1L  (L2, t [1]);
	__dline3L (L13, t);
	__line3L  (L23, t);
	__quad3S  (N, t, L1, L2, L13, L23);

	__line1L  (L1, t [0]);
	__dline1L (L2, t [1]);
	__line3L  (L13, t);
	__dline3L (L23, t);
	__quad3S  (N + 12, t, L1, L2, L13, L23);
}


