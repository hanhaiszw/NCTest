/*
 * NCUtils.h
 * ���빤����
 *  Created on: 2017��11��29��
 *      Author: Administrator
 */
#pragma once
#include "GF.h"
//C++�� û��byte�ؼ���
typedef unsigned char byte;

class NCUtils {
private:
	static GF gf;
	NCUtils() = delete;
public:	
	//���������
    static byte** multiply(byte** matrix1, int row1, int col1, byte** matrix2,
            int row2, int col2);
    //ֻ�Է�������
    static byte** inverse(byte** matrix, int nK);
	//����
    static int getRank(byte** matrix, int nRow, int nCol);
};


