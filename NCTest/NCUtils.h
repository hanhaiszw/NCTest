/*
 * NCUtils.h
 * 编码工具类
 *  Created on: 2017年11月29日
 *      Author: Administrator
 */
#pragma once
#include "GF.h"
#include <vector>
using std::vector;

//C++中 没有byte关键字
typedef unsigned char byte;

class NCUtils {
private:
	static GF gf;
	NCUtils() = delete;
public:	
	//两矩阵相乘
    static byte** multiply(byte** matrix1, int row1, int col1, byte** matrix2,
            int row2, int col2);
    //只对方阵求逆
    static byte** inverse(byte** matrix, int nK);
	//求秩
    static int getRank(byte** matrix, int nRow, int nCol);
	// 生成随机数据 0-256 
	// 数据量比较小  使用vector较方便  不用再手动释放空间
	static vector<vector<byte>> generateRandMatrix(int row, int col);
	static int getRank(vector<vector<byte>> matrix);
};



