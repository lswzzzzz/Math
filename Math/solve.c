#include "solve.h"
#include <assert.h>
#include <math.h>
#include "stdio.h"

const double PI = 3.1415926535898;

double getPow(double x, double k){
	if (x <= 0){
		return -pow(-x, k);
	}
	else{
		return pow(x, k);
	}
};

//   1/4g^2*t^4 + (gy-r^2)*t^2 - (2r^3*t)/v + y^2 + x^2 - (r^4) / v^2 = 0
int getSita(double x, double y, double g, double v, double r, double *sita)
{
	double a, b, c, d, e;
	a = 0.25 * g * g;
	b = 0;
	double r2 = r * r;
	c = g * v - r2;
	d = -(2 * r2 * r) / v;
	e = y * y + x * x - (r2 * r2) / (v * v);
	double result[4];
	int length = biquadratic(a, b, c, d, e, result);
	double realResult[4];
	int len = 0;
	for (int i = 0; i < length; i++){
		if (result[i] > 0){
			realResult[len++] = result[i];
		}
	}
	for (int i = 0; i < len; i++){
		double tansita = (y + g * pow(realResult[i], 2) * 0.5) / x;
		sita[i] = atan(tansita);
	}
	return len;
}

//http://www.xieguofang.cn/Maths/Quartic/Quartic_Xie_formulae_real_form.htm
int biquadratic(double a, double b, double c, double d, double e, double* result){
	b = b / 4;
	c = c / 6;
	d = d / 4;
	double H = b * b - a * c;
	double I = a * e - 4 * b * d + 3 * c * c;
	double G = a * a * d - 3 * a * b * c + 2 * getPow(b, 3);
	double J = (4 * getPow(H, 3) - a * a * H * I - G * G) / getPow(a, 3);
	double delta = getPow(I, 3) - 27 * J * J;
	int length = 0;
	if (G != 0.0 && I == 0.0 && J == 0.0){
		assert(H > 0);
		int sgn = G > 0 ? 1 : -1;
		length = 2;
		result[0] = (-b + sgn * sqrt(H)) / a;
		result[1] = (-b - 3 * sgn * sqrt(H)) / a;
	}
	else if (G == 0){
		double square = 12 * H * H - a * a * I;
		assert(square > 0);
		square = sqrt(square);
		double bigSquare;
		bigSquare = 3 * H + square;
		if (bigSquare > 0){
			result[length++] = (-b + sqrt(bigSquare)) / a;
			result[length++] = (-b - sqrt(bigSquare)) / a;
		}
		bigSquare = 3 * H - square;
		if (bigSquare > 0){
			result[length++] = (-b + sqrt(bigSquare)) / a;
			result[length++] = (-b - sqrt(bigSquare)) / a;
		}
	}
	else if (delta < 0){
		double t = a / 2 * (getPow(-J+getPow(-delta/27, 2), 1.0/3.0)) + getPow(-J-getPow(-delta/27, 2), 1.0/3.0) + H;
		double square = fabs(G) / sqrt(t) - t + 3 * H;
		int sgn = G > 0 ? 1 : -1;
		length = 2;
		result[0] = (-b - sgn * sqrt(t) + sqrt(square)) / a;
		result[1] = (-b - sgn * sqrt(t) - sqrt(square)) / a;
	}
	else {
		I = fabs(I);
		int sgn = G > 0 ? 1 : -1;
		double sita = acos(-J / sqrt(pow(I, 3) / 27));
		double y1 = a * sqrt(I / 3) * cos(sita / 3) + H;
		double y2 = a * sqrt(I / 3) * cos(sita / 3 + 2 * PI / 3) + H;
		double y3 = a * sqrt(I / 3) * cos(sita / 3 - 2 * PI / 3) + H;
		assert(y1 > 0);
		assert(y2 > 0);
		assert(y3 > 0);
		y1 = sqrt(y1);
		y2 = sqrt(y2);
		y3 = sqrt(y3);
		int s = -sgn;
		result[length++] = (-b + s * y1 + y2 + y3) / a;
		result[length++] = (-b + s * y1 - y2 - y3) / a;
		result[length++] = (-b - s * y1 + y2 - y3) / a;
		result[length++] = (-b - s * y1 - y2 + y3) / a;
	}
	return length;
};

///http://baike.baidu.com/picture/606391/606391/0/bd315c6034a85edf250e387e4b540923dc5475e9#aid=0&pic=fcfaaf51f3deb48f11246783f21f3a292df5781a
int thrice(double a, double b, double c, double d, double* result){
	double A = b * b - 3 * a * c;
	double B = b * c - 9 * a * d;
	double C = c * c - 3 * b * d;
	double delta = B * B - 4 * A * C;
	int length;
	if (A == B && B == 0){
		length = 1;
		result[0] = -b / (3 * a);
	}
	else if (delta > 0){
		double y1, y2;
		y1 = A*b + 3 * a * (-B + sqrt(delta)) * 0.5;
		y2 = A*b + 3 * a * (-B - sqrt(delta)) * 0.5;
		double powy1, powy2;
		powy1 = getPow(y1, 1.0 / 3);
		powy2 = getPow(y2, 1.0 / 3);
		length = 1;
		result[0] = (-b - powy1 - powy2) / (3 * a);
	}
	else if (delta == 0){
		length = 2;
		assert(A != 0.0);
		double K = B / A;
		result[0] = -b / a + K;
		result[1] = -K * 0.5;
	}
	else{
		length = 3;
		double T = (2 * A* b - 3 * a * B) * 0.5 / (getPow(A, 3 / 2));
		assert(T > -1.0);
		assert(T < 1.0);
		assert(A > 0);
		double arccos = acos(T);
		double sitac = cos(arccos / 3);
		double sitas = sin(arccos / 3);
		double d3a = 1 / (3 * a);
		result[0] = (-b - 2 * sqrt(A) * sitac) * d3a;
		result[1] = (-b + sqrt(A)* (sitac + sqrt(3) * sitas)) * d3a;
		result[2] = (-b + sqrt(A)* (sitac - sqrt(3) * sitas)) * d3a;
	}
	return length;
}

int quadratic(double a, double b, double c, double* result){
	double delta = b * b - 4 * a * c;
	assert(delta >= 0.0);
	int length;
	if (delta == 0){
		length = 1;
		result[0] = -b / a * 0.5;
	}
	else{
		length = 2;
		result[0] = (-b + sqrt(delta)) / a * 0.5;
		result[1] = (-b - sqrt(delta)) / a * 0.5;
	}
	return length;
}