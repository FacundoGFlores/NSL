/*****************************************************************************
Name		: 	linearMatrix.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the implementation of the linearMatrix class.
License     :
    			Copyright (C) 2012 Flores, Facundo Gabriel

    			This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.

    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************** */

#include <iostream>

#include <fstream>

#include <cstdlib>

#include <ctime>

#include <vector>

#include <linearMatrix.hpp>

using namespace std;

#define FALSE 0
#define TRUE 1

/**
 * We test the sizes entered.
 * @param W the width set
 * @param H the height set
 * @return 0: non errors 1: errors
 */
bool linearMatrix::Test_Sizes(const int W, const int H)
{
	if(W < 1 || W > 32500)
		return FALSE;
	if(H < 1 || H> 32500)
		return FALSE;
	return TRUE;
}

/**
 * A constructor. In this case we have to
 * initialize the dimension of the vector.
 * Then we can generate random matrices
 * using the method of the class.
 * Don't worry if you want to use a 1D vector
 * just set Height to 1.
 * @param Width the width of the matrix.
 * @param Height the height of the matrix
 */
linearMatrix::linearMatrix(const int Width, const int Height)
{
	if(Test_Sizes(Width, Height))
	{
		myWidth = Width;
		myHeight = Height;
		myVector.resize(myWidth * myHeight);
	}
	else
	{
		cout<<"Bad sizes!"<<endl;
		exit(1);
	}

}

/**
 * A copy constructor. In this case we send
 * a constant reference to another linearMatrix
 * object, and the constructor copy all of the
 * attributes in "that".
 * @param that a constant linearMatrix reference
 */
linearMatrix::linearMatrix(const linearMatrix &that)
{
	myWidth = that.myWidth;
	myHeight = that.myHeight;
	myVector = that.myVector;
}

/**
 * An overload of the "=" operator. This is part
 * of the rule of three in C++. We have to give
 * to the user the assignment possibility. Use:
 * linearMatrix = linearMatrix.
 * @param that a constant linearMatrix reference
 */
linearMatrix &linearMatrix :: operator =(const linearMatrix& that)
{
	if(this != &that)
	{
		myHeight = that.myHeight;
		myWidth = that.myWidth;

		myVector = that.myVector;
	}
	return *this;
}

/**
 * We generate a random matrix using the
 * Width and Height set.
 */
void linearMatrix::Generate_Random_Matrix(void)
{
	srand(time(NULL));
	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			// We are putting a number between 0 and 99
			myVector[i * myHeight + j] = rand() % 100;
		}
}

/**
 * We generate a diagonal dominant matrix
 * using the Width and Height set.
 * This is a symmetric matrix for a general
 * purpose.
 */
void linearMatrix::Generate_Diagonal_Dominant(void)
{
	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			if(i == j) //Diagonal item
				myVector[i * myHeight + j] = 500;
			else
				myVector[i * myHeight + j] = 0.025;
		}
}

/**
 * If you have got a little matrix(less than 100
 * elements) the function will write it into
 * the console. If you have got a big matrix,
 * the function will show you a message.
 */
void linearMatrix::Show_Little_Matrix(void) const
{
	cout.precision(5);
	if(myHeight < 11 && myWidth < 11)
	{
		for(int i = 0; i < myWidth; i++)
		{
			for(int j = 0; j < myHeight; j++)
			{
				cout<<myVector[i * myHeight+ j]<<" ";
			}
			cout << endl;
		}
	}
	else
	{
		cout << "I can show only little matrices in the console " << endl;
	}
}

/**
 * We write a binary file called "Path". If it does
 * not exist then it creates a new one. If it
 * already exists, it replace another one file.
 * @param Path the name of the binary file
 */
void linearMatrix::Write_File_Matrix(const char *Path)
{
    ofstream outFILE(Path, ios::out | ofstream::binary);
    outFILE.write(reinterpret_cast<const char *>(&myWidth), sizeof(float));
    outFILE.write(reinterpret_cast<const char *>(&myHeight), sizeof(float));
    outFILE.write(reinterpret_cast<const char *>(&myVector[0]), myVector.size() * sizeof(float));
}

/**
 * We read a binary file, if it does not exist
 * then the program exits.
 * @param Path
 */
void linearMatrix::Read_File_Matrix(const char *Path)
{
    float myfloat;

    ifstream inFILE(Path, ios::in | ifstream::binary);

    if(!inFILE)
    {
        cerr << "Cannot open the file" << endl;
        exit(1);
    }

    inFILE.read(reinterpret_cast<char *>(&myWidth), sizeof(float));

    inFILE.read(reinterpret_cast<char *>(&myHeight), sizeof(float));
    int i = 0;
    while(inFILE.read(reinterpret_cast<char *>(&myfloat), sizeof(float)))
    {
    	myVector[i] = myfloat;
    	i++;
    }
}

/**
 * The function returns the stored matrix
 * @return the matrix
 */
vector<float> linearMatrix::Get_Vector() const
{
	return myVector;
}

/**
 * The function returns the width of the
 * stored matrix
 * @return Width
 */
int linearMatrix::Get_Width() const
{
	return myWidth;
}

/**
 * The function returns the height of the
 * stored matrix
 * @return Height
 */
int linearMatrix::Get_Height() const
{
	return myHeight;
}

/**
 * We put a c++ vector(std::vector) inside
 * our linearMatrix object, it has to be
 * created previously.
 * @param v a std::vector, that is the matrix
 * @param Width the Width of the vector passed
 * @param Height the "Height of the vector passed"
 */
void linearMatrix::Set_Matrix(const vector<float> &v, const int Width, const int Height)
{
	if(Test_Sizes(Width, Height))
	{
		myWidth = Width;
		myHeight = Height;
		myVector = v;
	}
	else
	{
		cout<<"Bad sizes!"<<endl;
		exit(1);
	}
}

/**
 * We generate a null matrix using the
 * Width and Height set.
 */
void linearMatrix::Generate_Null_Matrix(void)
{
	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			myVector[i * myHeight + j] = 0;
		}
}

/**
 * The function returns the diagonal
 * vector of the square matrix
 * @return
 */
vector<float> linearMatrix::Get_Diagonal() const
{
	if(myHeight != myWidth)
	{
		cout << "Cannot obtain diagonal matrix" << endl;
		exit(1);
	}
	vector<float> D;
	D.resize(myHeight);
	for(int i = 0; i < myHeight; i++)
	{
		D[i] = myVector[i * myHeight + i];
	}

	return D;
}

/**
 * The function returns the transpose of the stored square matrix
 * @return the transpose matrix
 */
vector<float> linearMatrix::Get_Transpose_Vector() const
{
	if(myHeight != myWidth)
	{
		cout << "It is not a square matrix" << endl;
		exit(1);
	}

	vector<float> transpose;
	transpose.resize(myHeight * myWidth);

	for(int i = 0; i < myHeight - 2; i++)
		for(int j = i + 1; j < myWidth - 1; j++)
		{
			transpose[i * myHeight + j] = myVector[i * myHeight + j];
			transpose[j * myWidth + i] = myVector[i * myHeight + j];
		}

	return transpose;
}
/**
 * We generate an identity matrix.
 */
void linearMatrix::Generate_Idendity_Matrix(void)
{
	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			if(i != j)
				myVector[i * myHeight + j] = 0;
			else
				myVector[i * myHeight + j] = 1;
		}
}

/**
 * The function returns a LU matrix
 * from the pre-existent
 * @return the LU
 */
linearMatrix linearMatrix::Get_LU() const
{
	if(myHeight != myWidth)
	{
		cout << "Cannot obtain LU matrix" << endl;
		exit(1);
	}

	vector<float> LU;
	LU.resize(myHeight * myWidth);

	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			if(i == j) //Diagonal item
				LU[i * myHeight + j] = 0;
			else
				LU[i * myHeight + j] = myVector[i * myHeight + j];
		}

	linearMatrix LU_MATRIX(myWidth, myHeight);
	LU_MATRIX.Set_Matrix(LU, myWidth, myHeight);

	return LU_MATRIX;
}

/**
 * The function returns a L matrix
 * from the pre-existent
 * @return the L
 */
linearMatrix linearMatrix::Get_L() const
{
	if(myHeight != myWidth)
	{
		cout << "Cannot obtain L matrix" << endl;
		exit(1);
	}

	vector<float> L;
	L.resize(myHeight * myWidth);

	for(int i = 0; i < myWidth; i++)
		for(int j = 0; j < myHeight; j++)
		{
			if(i > j)
				L[i * myHeight + j] = 0;
			else
				L[i * myHeight + j] = myVector[i * myHeight + j];
		}

	linearMatrix L_MATRIX(myWidth, myHeight);
	L_MATRIX.Set_Matrix(L, myWidth, myHeight);

	return L_MATRIX;
}


/**
 * We put an element into the matrix given "x" and "y" indexes.
 * @param value The data to insert into the matrix
 * @param x The x index
 * @param y The y index
 */
void linearMatrix::Set_Element(const float value, const int x, const int y)
{
	if(x >= 0 && x < myWidth)
		if(y >=0 && y < myHeight)
			myVector[y * myHeight + x] = value;
		else
		{
			cout << "Index out of range!" << endl;
		}
	else
	{
		cout << "Index out of range" << endl;
	}
}

/**
 * We take al element given "x" and "y" indexes
 * @param x The x index
 * @param y The y index
 * @return the element requested
 */
float linearMatrix::Get_Element(const int x, const int y)
{
	if(x >= 0 && x < myWidth)
		if(y >=0 && y < myHeight)
			return myVector[y * myHeight + x];
		else
		{
			cout << "Index out of range" << endl;
			exit(1);
		}
	else
	{
		cout << "Index out of range" << endl;
		exit(1);
	}
}
