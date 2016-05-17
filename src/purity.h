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
#include <iostream>


using namespace Eigen;
using namespace std;

#ifndef cacc_purity
#define cacc_purity

const int chnpur=4;

template <typename myfloat>
int purity(Matrix<myfloat,Dynamic,Dynamic> &a, MatrixXf &pur, Matrixidx &tid,char purflag=1)
{
    typedef Matrix<myfloat,Dynamic,Dynamic> mymat1;
    mylong rows=(mylong)a.rows();
    mylong cols=(mylong)a.cols();

    if(cols==chnpur)
    {
        mymat1 tma=mymat1::Constant(rows,1,0);
        tid=Matrixidx::Constant(rows,1,0);
        for(int8_t k=0;k<chnpur;k++)
            for(mylong i=0;i<rows;i++)
                if(a.coeff(i,k)>tma.coeff(i))
                    tma(i)=a.coeff(i,k),tid(i)=k;
        if(purflag)
            pur=(tma.template cast<float>().array()/a.cwiseAbs().rowwise().sum().array().max(0.000001f)).max(0);
    }
    else if(rows==chnpur)
    {
        mymat1 tma=mymat1::Constant(1,cols,0);
        tid=Matrixidx::Constant(1,cols,0);
        for(int8_t k=0;k<chnpur;k++)
            for(mylong i=0;i<cols;i++)
                if(a.coeff(k,i)>tma.coeff(i))
                    tma(i)=a.coeff(k,i),tid(i)=k;
        if(purflag)
            pur=(tma.template cast<float>().array()/a.cwiseAbs().colwise().sum().array().max(0.000001f)).max(0);
    }
    else return -1;
    return 0;
}

template <typename myfloat>
int purity(Matrix<myfloat,Dynamic,Dynamic> &a,MatrixXf &pur)
{
    Matrixidx tid;
    return purity(a,pur,tid);
}

template <typename myfloat>
int purity( Matrix<myfloat,Dynamic,Dynamic> &a,Matrixidx &tid)
{
    MatrixXf pur;
    return purity(a,pur,tid,0);
}




#endif
