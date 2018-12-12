/*
 * NCUtils.cpp
 *
 *  Created on: 2017年11月29日
 *      Author: Administrator
 */
#include "stdafx.h"
#include "NCUtils.h"


// 有限域运算库
GF NCUtils::gf;

byte** NCUtils::multiply(byte** matrix1, int row1, int col1, byte** matrix2,
        int row2, int col2) {
    if (col1 != row2) {
        return NULL;
    }

    byte **result = new byte*[row1];
    for (int i = 0; i < row1; ++i) {
        result[i] = new byte[col2];
    }
    byte temp = 0;
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < col2; ++j) {
            temp = 0;
            for (int k = 0; k < col1; ++k) {
                temp = gf.add(temp, gf.mul(matrix1[i][k], matrix2[k][j]));
            }
            result[i][j] = temp;
        }
    }
    return result;
}

byte** NCUtils::inverse(byte** matrix, int nK) {
    //先计算矩阵的秩  不是满秩不能求秩
    int rank = getRank(matrix, nK, nK);
    if (rank != nK) {
        return NULL;
    }
    int k = nK;
    int nCol = nK;

    unsigned int **M = new unsigned int *[k];
    for (int i = 0; i < k; ++i) {
        M[i] = new unsigned int[k];
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            M[i][j] = matrix[i][j];  // Copy the coefficient to M.
        }
    }

	/*
    //unsigned int IM[k][k];
    unsigned int **IM = new unsigned int *[k];
    for (int i = 0; i < k; ++i) {
        IM[i] = new unsigned int[k];
    }*/
	byte **IM = new byte*[k];
	for (int i = 0; i < k; ++i) {
		IM[i] = new byte[k];
	}

    // Init
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            if (i == j) {
                IM[i][j] = 1;
            } else {
                IM[i][j] = 0;
            }
        }
    }
    /************************************************************************/
    /* Step 1. Change to a lower triangle matrix.                           */
    /************************************************************************/
    for (int i = 0; i < nCol; i++) {
        for (int j = i + 1; j < nCol; j++) {
            // Now, the main element must be nonsingular.
            GFType temp = gf.div(M[j][i], M[i][i]);

            for (int z = 0; z < nCol; z++) {
                M[j][z] = gf.add(M[j][z], gf.mul(temp, M[i][z]));
                IM[j][z] = gf.add(IM[j][z], gf.mul(temp, IM[i][z]));
            }
        }
    }
    /************************************************************************/
    /* Step 2. Only the elements on the diagonal are non-zero.                  */
    /************************************************************************/
    for (int i = 1; i < nCol; i++) {
        for (int j = 0; j < i; j++) {
            GFType temp = gf.div(M[j][i], M[i][i]);
            for (int z = 0; z < nCol; z++) {
                M[j][z] = gf.add(M[j][z], gf.mul(temp, M[i][z]));
                IM[j][z] = gf.add(IM[j][z], gf.mul(temp, IM[i][z]));
            }
        }
    }
    /************************************************************************/
    /* Step 3. The elements on the diagonal are 1.                  */
    /************************************************************************/
    for (int i = 0; i < nCol; i++) {
        if (M[i][i] != 1) {
            GFType temp = M[i][i];
            for (int z = 0; z < nCol; z++) {
                M[i][z] = gf.div(M[i][z], temp);
                IM[i][z] = gf.div(IM[i][z], temp);
            }
        }
    }
    
	/*
	// 可以直接返回IM
	// 为了保证通用性，这里转化为byte
	// 转化结果为byte
    byte** result = new byte*[nK];
    for (int i = 0; i < nK; ++i) {
        result[i] = new byte[nK];
    }
    for (int i = 0; i < nK; ++i) {
        for (int j = 0; j < nK; ++j) {
            result[i][j] = (byte) IM[i][j];
        }
    }
	*/

    //释放数组
    for (int i = 0; i < k; ++i) {
        delete[] M[i];
        //delete[] IM[i];
    }
    delete[] M;
    //delete[] IM;
    //return result;
	return IM;
}

int NCUtils::getRank(byte** matrix, int nRow, int nCol) {
    unsigned int **M = new unsigned int *[nRow];
    for (int i = 0; i < nRow; ++i) {
        M[i] = new unsigned int[nCol];
    }

    //unsigned int test = 0;
    for (int i = 0; i < nRow; i++) {
        for (int j = 0; j < nCol; j++) {
            //test = pData[i * nCol + j];
            M[i][j] = matrix[i][j];
        }
    }

    // Define a variable to record the position of the main element.
    int yPos = 0;

    for (int i = 0; i < nRow; i++) {
        // Find the main element which must be non-zero.
        bool bFind = false;
        for (int x = yPos; x < nCol; x++) {
            for (int k = i; k < nRow; k++) {
                if (M[k][x] != 0) {
                    // Exchange the two vectors.
                    for (int x = 0; x < nCol; x++) {
                        byte nVal = M[i][x];
                        M[i][x] = M[k][x];
                        M[k][x] = nVal;
                    }                      // We have exchanged the two vectors.
                    bFind = true;
                    break;
                }
            }
            if (bFind == true) {
                yPos = x;
                break;
            }
        }

        for (int j = i + 1; j < nRow; j++) {
            // Now, the main element must be nonsingular.
            unsigned int temp = gf.div(M[j][yPos], M[i][yPos]);
            for (int z = 0; z < nCol; z++) {
                M[j][z] = (byte) (gf.add(M[j][z],
                        gf.mul(temp, M[i][z])));
            }
        }
        //
        yPos++;

    }

    // The matrix becomes a scalar matrix. we need to make more elements become 0 with elementary transformations.
    yPos = 0;
    for (int i = 1; i < nRow; i++) {
        for (int j = 0; j < nCol; j++) {
            if (M[i][j] != 0) {
                // the main element is found.
                yPos = j;
                break;
            }
        }
        for (int k = 0; k < i; k++) {
            unsigned int temp = gf.div(M[k][yPos], M[i][yPos]);
            for (int z = 0; z < nCol; z++) {
                M[k][z] = (byte) (gf.add(M[k][z],
                        gf.mul(temp, M[i][z])));
            }
        }
    }

    int nRank = 0;
    // Get the rank.
    for (int i = 0; i < nRow; i++) {
        int nNonzero = 0;
        for (int j = 0; j < nCol; j++) {
            if (M[i][j] != 0) {
                nNonzero++;
            }
        }
        // If there is only one nonzero element in the new matrix, it is concluded an original packet is leaked.
        if (nNonzero > 0) {
            // Leaked.
            nRank++;
        }
    }

    //释放内存
    for (int i = 0; i < nRow; ++i) {
        delete[] M[i];
    }
    delete[] M;

    return nRank;
}


