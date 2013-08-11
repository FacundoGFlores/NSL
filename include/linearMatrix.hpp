/*****************************************************************************
Name		: 	linearMatrix.hpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	linearMatrix class. It helps us leading with
				1D linear vectors. It's useful when we work
				with the cuBlas library. Remember there are a
				linearAlgebra.hpp which handle some kind of
				basic linear algebra and that library uses it.
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


#ifndef LINEARMATRIX_HPP_
#define LINEARMATRIX_HPP_

#include <vector>
using std::vector;

class linearMatrix
{
public:

	/**
	 *  A default constructor.
	 */
	linearMatrix()
	{
		myWidth = 0;
		myHeight = 0;
		myVector.resize(0);
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
	linearMatrix(const int Width, const int Heigth);

	/**
	 * A copy constructor. In this case we send
	 * a constant reference to another linearMatrix
	 * object, and the constructor copy all of the
	 * attributes in "that".
	 * @param that a constant linearMatrix reference
	 */
	linearMatrix(const linearMatrix &that);

	/**
	 * An overload of the "=" operator. This is part
	 * of the rule of three in C++. We have to give
	 * to the user the assignment possibility. Use:
	 * linearMatrix = linearMatrix.
	 * @param that a constant linearMatrix reference
	 */
	linearMatrix& operator=(const linearMatrix& that);

	/**
	 * We read a binary file, if it does not exist
	 * then the program exits.
	 * @param Path
	 */
	void Read_File_Matrix(const char *Path);

	/**
	 * We write a binary file called "Path". If it does
	 * not exist then it creates a new one. If it
	 * already exists, it replace another one file.
	 * @param Path the name of the binary file
	 */
	void Write_File_Matrix(const char *Path);

	/**
	 * We generate a diagonal dominant matrix
	 * using the Width and Height set.
	 * This is a symmetric matrix for a general
	 * purpose.
	 */
	void Generate_Diagonal_Dominant(void);

	/**
	 * We generate an identity matrix.
	 */
	void Generate_Idendity_Matrix(void);

	/**
	 * We generate a random matrix using the
	 * Width and Height set.
	 */
	void Generate_Random_Matrix(void);

	/**
	 * We generate a null matrix using the
	 * Width and Height set.
	 */
	void Generate_Null_Matrix(void);

	/**
	 * We put a c++ vector(std::vector) inside
	 * our linearMatrix object, it has to be
	 * created previously.
	 * @param v a std::vector, that is the matrix
	 * @param Width the Width of the vector passed
	 * @param Height the "Height of the vector passed"
	 */
	void Set_Matrix(const vector<float> &v, const int Width, const int Height);

	/**
	 * We put an element into the matrix given "x" and "y" indexes.
	 * @param value The data to insert into the matrix
	 * @param x The x index
	 * @param y The y index
	 */
	void Set_Element(const float value, const int x, const int y);

	/**
	 * We take al element given "x" and "y" indexes
	 * @param x The x index
	 * @param y The y index
	 * @return the element requested
	 */
	float Get_Element(const int x, const int y);


	/**
	 * If you have got a little matrix(less than 100
	 * elements) the function will write it into
	 * the console. If you have got a big matrix,
	 * the function will show you a message.
	 */
	void Show_Little_Matrix(void) const;

	/**
	 * The function returns the stored matrix
	 * @return the matrix
	 */
	vector<float> Get_Vector() const;

	/**
	 * The function returns the transpose of the stored matrix
	 * @return the transpose matrix
	 */
	vector<float> Get_Transpose_Vector() const;

	/**
	 * The function returns the width of the
	 * stored matrix
	 * @return Width
	 */
	int Get_Width() const;

	/**
	 * The function returns the height of the
	 * stored matrix
	 * @return Height
	 */
	int Get_Height() const;

	/**
	 * The function returns the diagonal
	 * vector of the square matrix
	 * @return Diagonal vector
	 */
	vector<float> Get_Diagonal() const;

	/**
	 * The function returns a LU matrix
	 * from the pre-existent
	 * @return the LU
	 */
	linearMatrix Get_LU() const;

	/**
	 * The function returns a L matrix
	 * from the pre-existent
	 * @return the L
	 */
	linearMatrix Get_L() const;



private:
	int myWidth, myHeight; // Width and Height
	vector<float> myVector; // The matrix

	/**
	 * We test the sizes entered.
	 * @param Width the width set
	 * @param Height the height set
	 * @return 0: non errors 1: errors
	 */
	bool Test_Sizes(const int Width, const int Height);
};



#endif /* LINEARMATRIX_HPP_ */
