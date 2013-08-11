/*****************************************************************************
Name		: 	cg.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is implementation of the conjugate gradient.
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

#include <cmath>

#include <vector>

#include <sys/time.h>

#include <linearAlgebra.hpp>

#include <cg.hpp>

using namespace std;

/**
 * Solve a linear system of equations using the
 * conjugate gradient method.
 * @param Matrix_A Matrix A
 * @param Vector_B Vector B
 * @param Vector_X Vector X
 * @param e allowed error
 * @param MaxIter allowed iterations
 * @return the time spent
 */
long npConjugate_Gradient(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{
    timeval start;
    timeval end;

    long seconds, useconds, final;

	int k = 0;

	int N = Vector_B.Get_Width();

	linearMatrix tmp(N, 1);
	linearMatrix r(N, 1);
	linearMatrix p(N, 1);
	linearMatrix w(N, 1);


    gettimeofday(&start, 0);


	tmp = lAmatrixvecmult(Matrix_A, Vector_X);

	r = lAvectorvectorsub(Vector_B, tmp);

	float deltha = lAdotproduct(r, r);
	float myerror = sqrt(deltha);
	float deltha_old;
	float beta = 0;
	float alpha = 0;

	float norm_b = sqrt(lAdotproduct(Vector_B, Vector_B));
	float new_e = e * norm_b;


	while(myerror > new_e  && k < MaxIter)
	{
		k++;

		if(k == 1)
			p = r;
		else
		{
			beta = deltha / deltha_old;
			tmp = lAscalarvecmult(p, beta);
			p = lAvectorvectorsum(r, tmp);
		}

		w = lAmatrixvecmult(Matrix_A, p);

		alpha = deltha / lAdotproduct(p, w);

		tmp = lAscalarvecmult(p, alpha);

		Vector_X = lAvectorvectorsum(Vector_X, tmp);

		tmp = lAscalarvecmult(w, alpha);
		r = lAvectorvectorsub(r, tmp);

		deltha_old = deltha;
		deltha = lAdotproduct(r, r);

		myerror = sqrt(deltha);
	}

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;

	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	return final;
}

