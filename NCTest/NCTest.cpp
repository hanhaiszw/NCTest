// NCTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include "gf.hpp"
#include "NCUtils.h"

using std::cout;
using std::endl;

void test01() {
	gf_init(8, 0x00000187);
	int a = gf_add(20, 10);
	gf_uninit();
	cout << a << endl;
}

void test02() {
	byte** bArray = new byte*[2];
	for (int i = 0; i < 2; i++) {
		bArray[i] = new byte[2];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			if (i == j) {
				bArray[i][j] = 1;
			}
			else {
				bArray[i][j] = 0;
			}
		}
	}

	int rank = NCUtils::getRank(bArray, 2, 2);
	cout << "rank = " << rank << endl;

	for (int i = 0; i < 2; i++) {
		delete[] bArray[i];
	}
	delete bArray;
}

void test03() {
	byte** bArray = new byte*[2];
	/*for (int i = 0; i < 2; i++) {
		bArray[i] = new byte[2];
	}
	bArray[0][0] = 123;
	bArray[0][1] = 17;
	bArray[1][0] = 250;
	bArray[1][1] = 47;*/

	bArray[0] = new byte[2]{ 123,17 };
	bArray[1] = new byte[2]{ 250,47 };
	

	byte** ret = NCUtils::inverse(bArray, 2);
	/*
	 185 193
	 25   77
	 */
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << (int)ret[i][j] << "  ";
		}
		cout << endl;
	}

	for (int i = 0; i < 2; i++) {
		delete[] bArray[i];
		delete[] ret[i];
	}
	delete bArray;
	delete ret;
}


void test04() {
	byte** a = new byte*[2];
	byte** b = new byte*[2];
	/*for (int i = 0; i < 2; i++) {
	a[i] = new byte[2];
	b[i] = new byte[2];
	}*/
	a[0] = new byte[2]{ 112,232 };
	a[1] = new byte[2]{ 23,110 };
	b[0] = new byte[2]{ 11,31 };
	b[1] = new byte[2]{ 9,252 };

	//byte a[2][2] = { 112,232,23,110 };
	//byte b[2][2] = { 11,31,9,252 };

	byte** result = NCUtils::multiply(a, 2, 2, b, 2, 2);
	/*
	 45  47
	 145 109
	*/
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++)
		{
			cout << (int)result[i][j] << "  ";
		}
		cout << endl;
	}
	//cout << endl;
	for (int i = 0; i < 2; i++) {
		delete[] a[i];
		delete[] b[i];
		delete[] result[i];
	}
	delete a;
	delete b;
	delete result;
}

int main()
{
	test04();
	system("pause");
	return 0;
}

