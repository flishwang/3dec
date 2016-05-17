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
#include "cacc_funs.h"


#include <iostream>


void pick(Ref<MatrixXf>z,Ref<MatrixXb>idx,MatrixXf &b)
{
    mylong j,n;
    MatrixXf a;
    n=(mylong)(idx.array()>0).count();
    if(idx.cols()==1)
    {
        a.resize(n,z.cols());
        j=0;
        for(mylong i=0;i<n;i++)
        {
            while(!idx.coeff(j,0)) j++;
            a.row(i)=z.row(j),j++;
        }
    }
    else
    {
        a.resize(z.rows(),n);
        j=0;
        for(mylong i=0;i<n;i++)
        {
            while(!idx.coeff(0,j)) j++;
            a.col(i)=z.col(j),j++;
        }
    }
    b=a;
}

void picki(Ref<MatrixXf>z,Ref<MatrixXl>idx,MatrixXf &b)
{
    mylong n;
    MatrixXf a;

    if(idx.cols()==1)
    {
        n=(mylong)idx.rows();
        a.resize(n,z.cols());
        
        for(mylong i=0;i<n;i++)
        {
            a.row(i)=z.row(idx(i,0));
        }
    }
    else
    {
        n=(mylong)idx.cols();
        a.resize(z.rows(),n);
        for(mylong i=0;i<n;i++)
        {
            a.col(i)=z.col(idx(0,i));
        }
    }
    b=a;
}




void pick(MatrixXf &z,Ref<MatrixXb>idx)
{
    pick(z,idx,z);
}
void picki(MatrixXf &z,Ref<MatrixXl>idx)
{
    picki(z,idx,z);
}