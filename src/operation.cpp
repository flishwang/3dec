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

using namespace Eigen;
using namespace std;

dym Cif::topdim()
{
    return (dim+1)%3;
}

int Cif::copycif(Cif &zz,int copydata)
{
    if(copydata==1)
        data=zz.data;
    
    dim=zz.dim;
    for(int i=0;i<3;i++)
        size[i]=zz.size[i];
    return 1;
}


void Cif::assignslice(dym tdim, mylong index, MatrixXf &sli)
{
    if (tdim==topdim())
    {
        Map<MatrixXf> mapdata(data.data(),size[tdim],size[(tdim+1)%3]*size[(tdim+2)%3]);
        mapdata.row(index)=Map<VectorXf>(sli.data(),size[(tdim+1)%3]*size[(tdim+2)%3]);
        return;
    }
    rotate(tdim);
    data.col(index)=Map<VectorXf>(sli.data(),size[(tdim+1)%3]*size[(tdim+2)%3]);
}

void Cif::getslice(dym tdim, mylong index,MatrixXf &sli)
{
    //dym dim=Cif::dim;
    if (tdim==topdim())
    {
        Map<MatrixXf> mapdata(data.data(),size[tdim],size[(tdim+1)%3]*size[(tdim+2)%3]);
        sli=mapdata.row(index);
        sli.resize(size[(tdim+1)%3],size[(tdim+2)%3]);
    }
	else if (tdim == dim)
	{
		sli = data.col(index);
		sli.resize(size[(dim + 1) % 3], size[(dim + 2) % 3]);
	}
	else
		cerr << "Error: cannot get slice of middle dimension\n",
		sli.resize(0, 0);
}


void Cif::rotate(dym tdim)
{
   // dym dim=Cif::dim;
	if (tdim == dim) {
		data.resize(size[(tdim + 1) % 3] * size[(tdim + 2) % 3], size[tdim]);
		return;
	}
	cerr << "Rotating to End:" << tdim << "...";
    if ((tdim-dim+3)%3==1)
        data.resize(size[tdim],size[(tdim+1)%3]*size[(tdim+2)%3]),
        data.transposeInPlace(),
        cerr<<"Succeed\n";
    else
		data.resize(size[(dim + 1) % 3] * size[(dim + 2) % 3], size[dim]),
        data.transposeInPlace(),
        data.resize(size[(tdim+1)%3]*size[(tdim+2)%3],size[tdim]),
        cerr<<"Succeed\n";
   	dim=tdim;
}
void Cif::rotate_top(dym tdim)
{
    rotate((tdim+2)%3);
}

void Cif::multi(dym tdim,MatrixXf &c)
{
  //  dym dim=Cif::dim;
    if(dim==tdim)
        data=data*c.transpose();
    else
        rotate_top(tdim),
        data.resize(size[tdim],size[(tdim+1)%3]*size[(tdim+2)%3]),
        data=c*data,
        data.resize(size[tdim]*size[(tdim+1)%3],size[(tdim+2)%3]);
}




