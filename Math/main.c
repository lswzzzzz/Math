#pragma once
#pragma execution_character_set("utf-8")
#include "solve.h"

int main(){
	int length;
	double result[10];
	//(x-1)*(x-2)*(x-3)*(x-4) = 0
	//==> x^4 -10 * x^3 + 35 * x^2 - 50 * x + 24 = 0
	length = biquadratic(1, -10, 35, -50, 24, &result);
	printf("One element and four equations\n");
	for (int i = 0; i < length; i++){
		printf("%f\n", result[i]);
	}
	printf("\n");
	printf("One element and three equations\n");
	//(y - 3.3)*(y-2.2)*(y-1.1)=0
	//==> y^3 - 6.6*y^2 + 13.31 * y - 7.986 = 0
	length = thrice(1, -6.6, 13.31, -7.986, result);
	for (int i = 0; i < length; i++){
		printf("%f\n", result[i]);
	}
	printf("\n");
	printf("One element and two equations\n");
	//(x-2)^2 = 0;
	//==> x^2 - 4 * x + 4 = 0
	length = quadratic(1, -4, 4, result);
	for (int i = 0; i < length; i++){
		printf("%f\n", result[i]);
	}
	printf("\n");
	getchar();
	return 0;
}