/*****************************************************************************
Name		: 	bicg.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is implementation of the bi-conjugate gradient.
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

#include <bicg.hpp>

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
long npBicg(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{
	long seconds, useconds, final;
	timeval start;
	timeval end;

	int k = 0;

	int N = Vector_B.Get_Width();

	linearMatrix Matrix_A_T(N, N);
	vector<float> t = Matrix_A.Get_Transpose_Vector();
	Matrix_A_T.Set_Matrix(t, N, N);

	linearMatrix v(N, 1);
	linearMatrix r(N, 1);
	linearMatrix tmp(N, 1);

	gettimeofday(&start, 0);

	//r = b-Ax
	tmp = lAmatrixvecmult(Matrix_A, Vector_X);
	r = lAvectorvectorsub(Vector_B, tmp);

	linearMatrix rhat(N, 1);
	rhat = r;

	float rho_old = 1;
	float rho = lAdotproduct(r,r);

	linearMatrix p(N, 1);
	linearMatrix phat(N, 1);

	p.Generate_Null_Matrix();
	phat.Generate_Null_Matrix();

	float alpha;
	float betha;
	float dotprod;

	float lasterror = sqrt(lAdotproduct(r, r));

	float norm_b = sqrt(lAdotproduct(Vector_B, Vector_B));
	float new_e = e * norm_b;


	while(lasterror > new_e  && k < MaxIter)
	{
		k++;

		betha = rho / rho_old;

		//p = r + betha * p
		tmp = lAscalarvecmult(p, betha);
		p = lAvectorvectorsum(r, tmp);

		//phat = rhat + betha * phat
		tmp = lAscalarvecmult(phat, betha);
		phat = lAvectorvectorsum(rhat, tmp);

		//v = Ap
		v = lAmatrixvecmult(Matrix_A, p);

		//alpha = rhok / <phat,v>
		dotprod = lAdotproduct(phat, v);
		alpha = rho / dotprod;

		//x = x + alpha * p
		tmp = lAscalarvecmult(p, alpha);
		Vector_X = lAvectorvectorsum(Vector_X, tmp);

		//r = r - alpha * v
		tmp = lAscalarvecmult(v, alpha);
		r = lAvectorvectorsub(r, tmp);

		//rhat = rhat - alpha * A^T * phat
		tmp = lAmatrixvecmult(Matrix_A_T, phat);
		tmp = lAscalarvecmult(tmp, alpha);
		rhat = lAvectorvectorsub(rhat, tmp);


		rho_old = rho;
		rho = lAdotproduct(rhat, r);

		lasterror = sqrt(lAdotproduct(r,r));
	}

	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;

	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	return final;
}

