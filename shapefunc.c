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
void line1L (double *N, double *t)  { __line1L (N, *t); }
void line2L (double *N, double *t)  { __line2L (N, *t); }
void line3L (double *N, double *t)  { __line3L (N, *t); }
void dline1L (double *N, double *t) { __dline1L (N, *t); }
void dline2L (double *N, double *t) { __dline2L (N, *t); }
void dline3L (double *N, double *t) { __dline3L (N, *t); }

//二次元要素の形状関数
//三角形要素
void tri1 (double *N, double *t)
{
	N [0] = t [0];
	N [1] = t [1];
	N [2] = 1 - t [0] - t [1];
}

void tri2 (double *N, double *t)
{
	double t1 = t[0], t2 = t [1], t3 = 1 - t1 - t2;
	N [0] = t1 * (2 * t1 - 1);
	N [1] = t2 * (2 * t2 - 1);
	N [2] = t3 * (2 * t3 - 1);
	N [3] = 4 * t1 * t2;
	N [4] = 4 * t2 * t3;
	N [5] = 4 * t3 * t1;
}

void tri3L (double *N, double *t)
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
void tri3S (double *N, double *t)
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

void dtri1 (double *N, double *t)
{
	N [0] = 1;
	N [1] = 0;
	N [2] = -1;
	N [3] = 0;
	N [4] = 1;
	N [5] = -1;

}

void dtri2 (double *N, double *t)
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

void dtri3L (double *N, double *t)
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

void tri3S (double *N, double *t)
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

