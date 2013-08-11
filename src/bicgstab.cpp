/*****************************************************************************
Name		: 	bicgstab.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is implementation of the bi-conjugate gradient stable-
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

#include <bicgstab.hpp>

using namespace std;

/**
 * Solve a linear system of equations using the
 * bi-conjugate gradient stable method.
 * @param Matrix_A Matrix A
 * @param Vector_B Vector B
 * @param Vector_X Vector X
 * @param e allowed error
 * @param MaxIter allowed iterations
 * @return the time spent
 */
long npBicgstab(const linearMatrix &Matrix_A, const linearMatrix &Vector_B, linearMatrix &Vector_X,
						const float e, const int MaxIter)
{

    timeval start;
    timeval end;

    long seconds, useconds, final;

	int N = Vector_B.Get_Width();

	linearMatrix tmp(N, 1); linearMatrix tmp1(N, 1); linearMatrix tmp2(N, 1);
	linearMatrix r(N, 1); linearMatrix rhat(N, 1); linearMatrix rhat0(N, 1);
	linearMatrix v(N, 1); linearMatrix p(N, 1);
	linearMatrix s(N, 1); linearMatrix t(N, 1);

	v.Generate_Null_Matrix();
	p.Generate_Null_Matrix();


    gettimeofday(&start, 0);


	tmp = lAmatrixvecmult(Matrix_A, Vector_X);
	r = lAvectorvectorsub(Vector_B, tmp);

	rhat0 = r;

	float betha;
	float alpha = 1;
	float w = 1;
	float rho_old = 1;
	float aux_norm;

	int k = 0;

	float rho = lAdotproduct(rhat0, r);
	float norm_b = sqrt(lAdotproduct(Vector_B, Vector_B));
	float new_e = e * norm_b;



	while(sqrt(lAdotproduct(r, r)) > new_e  && k < MaxIter)
	{
		k++;

		betha = (rho / rho_old) * (alpha / w);

		tmp1 = lAscalarvecmult(v, w);
		tmp2 = lAvectorvectorsub(p, tmp1);
		tmp = lAscalarvecmult(tmp2, betha);
		p = lAvectorvectorsum(r, tmp);

		v = lAmatrixvecmult(Matrix_A, p);

		alpha = rho / lAdotproduct(rhat0, v);

		tmp = lAscalarvecmult(v, alpha);
		s = lAvectorvectorsub(r, tmp);
		t = lAmatrixvecmult(Matrix_A, s);

		aux_norm = sqrt(lAdotproduct(t, t));
		w = lAdotproduct(t, s) / aux_norm;

		rho_old = rho;
		rho = lAdotproduct(rhat0, t);
		rho = rho * (-1) * w;

		tmp1 = lAscalarvecmult(s, w);
		tmp2 = lAscalarvecmult(p, alpha);
		tmp = lAvectorvectorsum(tmp1, tmp2);
		Vector_X = lAvectorvectorsum(Vector_X, tmp);

		tmp1 = lAscalarvecmult(t, w);
		r = lAvectorvectorsub(s, tmp1);
	}
	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;

	useconds = end.tv_usec - start.tv_usec;

	final = ((seconds) * 1000 + useconds/1000.0) + 0.5;

	if(k == MaxIter)
		cout << "Can't Solve that system of equations!" << endl;

	return final;
}

