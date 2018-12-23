/*
 * NCUtils.cpp
 *
 *  Created on: 2017年11月29日
 *      Author: Administrator
 */
#include "stdafx.h"
#include "NCUtils.h"
#include <time.h>
#include <stdlib.h>

#define random(x) (rand()%x)
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
   
	/************************************************************************/
	/* Start to get the inverse matrix!                                     */
	/************************************************************************/
	int **N = new int*[nK];
	for (int i = 0; i < nK; i++)
	{
		N[i] = new int[2 * nK];
	}

	// 一半写入原矩阵   一半写入单位矩阵
	for (int i = 0; i < nK; i++)
	{
		for (int j = 0; j < nK; j++)
		{
			N[i][j] = matrix[i][j];
		}


		for (int j = nK; j < 2 * nK; j++) {
			if (i == j - nK) {
				N[i][j] = 1;
			}
			else {
				N[i][j] = 0;
			}
		}
	}

	/************************************************************************/
	/* Step 1. Change to a lower triangle matrix.                           */
	/************************************************************************/
	for (int i = 0; i < nK; i++)
	{
		if (N[i][i] == 0) {
			// recode this line 
			// exchange
			for (int z = i + 1; z < nK; z++) {
				if (N[z][i] != 0) {
					for (int x = 0; x < 2 * nK; x++) {
						int temp = N[i][x];
						N[i][x] = N[z][x];
						N[z][x] = temp;
					}
					break;
				}
			}
		}

		for (int j = i + 1; j < nK; j++) {
			// Now, the main element must be nonsingular.
			GFType temp = gf.div(N[j][i], N[i][i]);
			for (int z = 0; z < 2 * nK; z++) {
				N[j][z] = gf.add(N[j][z], gf.mul(temp, N[i][z]));
			}
		}

	}



	/************************************************************************/
	/* Step 2. Only the elements on the diagonal are non-zero.                  */
	/************************************************************************/
	for (int i = 1; i < nK; i++)
	{
		for (int k = 0; k < i; k++) {
			GFType temp = gf.div(N[k][i], N[i][i]);
			for (int z = 0; z < 2 * nK; z++) {
				N[k][z] = gf.add(N[k][z], gf.mul(temp, N[i][z]));
			}
		}
	}

	/************************************************************************/
	/* Step 3. The elements on the diagonal are 1.                  */
	/************************************************************************/
	for (int i = 0; i < nK; i++) {
		if (N[i][i] != 1) {
			GFType temp = N[i][i];
			for (int z = 0; z < 2 * nK; z++) {
				N[i][z] = gf.div(N[i][z], temp);
			}
		}
	}


	/************************************************************************/
	/* Get the new matrix.                                                  */
	/************************************************************************/
	byte** result = new byte*[nK];
	for (int i = 0; i < nK; i++) {
		result[i] = new byte[nK];
	}
	for (int i = 0; i < nK; i++) {
		for (int j = 0; j < nK; j++) {
			result[i][j] = (byte)N[i][j + nK];
		}
	}

	for (int i = 0; i < nK; i++) {
		delete[] N[i];
	}
	delete[] N;

	return result;
}

// 不能使用于 nRow > nCol的情况
int NCUtils::getRank(byte** matrix, int nRow, int nCol) {
    GFType **M = new GFType *[nRow];
    for (int i = 0; i < nRow; ++i) {
        M[i] = new GFType[nCol];
    }

    //GFType test = 0;
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
            GFType temp = gf.div(M[j][yPos], M[i][yPos]);
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
            GFType temp = gf.div(M[k][yPos], M[i][yPos]);
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

vbArray NCUtils::generateRandMatrix(int row, int col)
{
	
	vbArray result(row);
	for (auto& v : result) {
		v.resize(col);
	}
	
	// 为了防止获取随机数过快，出现种子一样的情况
	static int old_seed = 0;
	int seed = time(nullptr);
	if (seed <= old_seed) {
		seed = old_seed + 1;
	}
	old_seed = seed;

	srand(seed);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i][j] = random(256);
		}
	}

	return result;
}



int NCUtils::getRank(vbArray matrix)
{
	int row = matrix.size();
	if (row == 0) {
		return 0;
	}
	int col = matrix[0].size();
	byte** mat = new byte*[row];
	for (int i = 0; i < row; i++) {
		mat[i] = new byte[col];
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			mat[i][j] = matrix[i][j];
		}
	}
	int rank = getRank(mat, row, col);

	// 释放内存
	for (int i = 0; i < row; i++) {
		delete[] mat[i];
	}
	delete[] mat;

	return rank;
}


vbArray NCUtils::inverse(vbArray& matrix) {
	int row = matrix.size();
	int col = matrix[0].size();
	byte** mat = new byte*[row];
	for (int i = 0; i < row; i++) {
		mat[i] = new byte[col];
	}

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			mat[i][j] = matrix[i][j];
		}
	}


	auto ret = inverse(mat, row);
	if (ret == NULL) {
		return vbArray();
	}

	vbArray vb(row);
	for (auto& v : vb) {
		v.resize(row);
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < row; j++) {
			vb[i][j] = ret[i][j];
		}
	}



	// 释放内存
	for (int i = 0; i < row; i++) {
		delete[] mat[i];
		delete[] ret[i];
	}
	delete[] mat;
	delete[] ret;

	return vb;

}