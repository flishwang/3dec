/*
This file is part of 3Dec, an accurate base-calling software for sequences.

Copyright (c) 2015, Bo Wang, Academy of Mathematics and Systems Science,
Chinese Academy of Sciences, Beijing 100190, China


This Source Code Form is subject to the terms of the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International Public License.
If a copy of the licence was not distributed with this file, you can obtain
one at http://creativecommons.org/licenses/by-nc-sa/4.0/
*/

#include "classtype.h"
#include "purity.h"

#include <iostream>


#ifndef _cacc_funs
#define _cacc_funs

using namespace Eigen;
using namespace std;

typedef Matrix<bool,Dynamic,Dynamic> MatrixXb;
typedef Matrix<mylong,Dynamic,Dynamic> MatrixXl;

//typedef float myfloat;
//Purity
template <typename myfloat>
int purity(Matrix<myfloat,Dynamic,Dynamic> &, MatrixXf&, Matrixidx &,char);
template <typename myfloat>
int purity(Matrix<myfloat,Dynamic,Dynamic> &,MatrixXf &);
template <typename myfloat>
int purity(Matrix<myfloat,Dynamic,Dynamic> &,Matrixidx &);


void pick(Ref<MatrixXf>, Ref<MatrixXb>, MatrixXf &);
void pick(MatrixXf &, Ref<MatrixXb>);

void picki(Ref<MatrixXf>, Ref<MatrixXl>, MatrixXf &);
void picki(MatrixXf &, Ref<MatrixXl>);

// Correct ACC from the first CIF using the second CIF.
int correct_acc(Cif &,Cif &,MatrixXf &loc,float ks=0.7,float blocksize=250,float buffersize=15);



/////////////////////////////////////////////
//min,max
template <typename mytype>
    mytype mymin(mytype a,mytype b)
{
    if(a<=b) return a;
    else return b;
}

template <typename mytype>
mytype mymax(mytype a,mytype b)
{
    if(a>=b) return a;
    else return b;
}




#endif