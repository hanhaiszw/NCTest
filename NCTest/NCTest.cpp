// NCTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include "gf.hpp"
#include "NCUtils.h"
#include <string>
#include <cmath>
#include <engine.h>
#include<string.h>
#include <time.h>
#include <stdlib.h>
#include <random>

#include <sstream>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using namespace std;
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
	delete[] bArray;
	delete[] ret;
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
	delete[] a;
	delete[] b;
	delete[] result;
}

void test05() {
	int row = 4;
	int col = 4;
	auto ret = NCUtils::generateRandMatrix(row, col);
	cout << "---------------------------------------------------" << endl;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			cout << (int)ret[i][j] << "   ";
		}
		cout << endl;
	}
	//cout << endl;
	cout << "rank = " << NCUtils::getRank(ret) << endl;
}

void test06() {
	vector<byte> v1{ 139, 151, 236, 112, 146, 23, 6};
	vector<byte> v2{ 246, 33,  185, 89,  47,  190,78 };
	vector<byte> v3{ 139, 151, 236, 112, 146, 23, 6 };
	vbArray varray{ v1,v2,v3 };

	/*for (int i = 0; i < 5; i++) {
		vector<byte> v(7, 0);
		v[i] = 1;
		varray.push_back(v);
	}*/
	int rank = NCUtils::getRank(varray);
	cout << "rank = " << rank << endl;
}

void test07() {
	unsigned int a = 0xFFFFFFFF;
	cout << a << endl;

}

void test08() {
	byte b1[][6] = {{ 1,0,0,0,0,0 },
					{ 0,1,0,0,0,0 },
					{ 0,0,1,0,0,0 },
					{ 0,0,0,1,0,0 },
					{ 0,0,0,0,1,0 },
					{ 0,0,0,0,0,1 } };

	vector<vector<byte>> vMatrix(6);
	/*for (auto& v : vMatrix) {
		v.resize(6);
	}*/
	/*for (int i = 0; i < 6; i++) {
			vector<byte> v(6,0);
			v[i] = 1;
			vMatrix.push_back(v);
		}*/

	vMatrix[0] = { 0,0,0,0,1,0 };
	vMatrix[1] = { 0,0,1,0,0,0 };
	vMatrix[2] = { 0,1,0,0,0,0 };
	vMatrix[3] = { 0,0,0,0,0,1 };
	vMatrix[4] = { 1,0,0,0,0,0 };
	vMatrix[5] = { 0,0,0,1,0,0 };


	auto ret = NCUtils::inverse(vMatrix);


	for (auto& v : ret) {
		for (auto& i : v) {
			cout << (int)i << "   ";
		}
		cout << endl;
	}
	
	cout << NCUtils::getRank(vMatrix) << endl;

	vbArray varr;
	varr.push_back({ 23,34,234 });
	varr.push_back({ 117,91,5 });
	varr.push_back({ 234,221,110 });
	auto ret1 = NCUtils::inverse(varr);
	for (auto& v : ret1) {
		for (auto& i : v) {
			cout << (int)i << "   ";
		}
		cout << endl;
	}

}



void test09() {
	Engine* m_pEngine;
	m_pEngine = engOpen(NULL);
	if (m_pEngine == NULL) {
		cout << "error!" << endl;
		exit(-1);
	}
	engSetVisible(m_pEngine,false);

	/*engEvalString(m_pEngine, "x=0:0.01:2*pi;");
	engEvalString(m_pEngine, "y=sin(x);");
	engEvalString(m_pEngine, "figure; plot(x,y,'g');");*/

	engEvalString(m_pEngine, "global m,poly");
	engEvalString(m_pEngine, "m=8;poly=391;gftable(m,poly);");
	engEvalString(m_pEngine, "gf_a=gf(23,m,poly);gf_b=gf(123,m,poly);");
	engEvalString(m_pEngine, "gf_result=gf_a*gf_b;");

	mxArray * result = engGetVariable(m_pEngine, "gf_result");
	
	//cout << string(mxArrayToString(result)) << endl;
	// 无法获取到result中的值   囧

	system("pause");
	engClose(m_pEngine);
}


void test10() {
	int num = ceil(10.5);
	cout << num << endl;

	// 保留两位小数  转化为string
	double pi = 3.14659265359;
	stringstream stream;
	stream << fixed << setprecision(2) << pi;
	string s = stream.str();
	cout << s << endl;
}


void test11() {
	// 会有叮咚提示音
	cout << "hello world\a" << endl;


}

int main()
{
	test11();
	
	system("pause");
	return 0;
}

