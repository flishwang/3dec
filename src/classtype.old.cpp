#include "classtype.h"

/*
class cif
{
typedef int dym;
private:
const dym dcluster=2, dchn=1, dcyc=0;
public:
// Three storage of data: dim=dcyc: cyc, chn, cluster;
//                        dim=dchn: chn, cluster, cyc;
//                        dim=dcluster: cluster,chn,cyc;
dym dim;
MatrixXf data;
unsigned long size[3];
// rotate:  rotate so that the para dim is at last place
void rotate(dym);
void rotate_top(dym);
// I/O operations  read/write data of para1 clusters, para2 cycles 
//                 from para3th clusters and para4th cycles.
void read(char *, long, long, long, long);
void cycread(char *, long, long ,long ,long);
void write(char *, long, long, long, long);
void cycwrite(char *, long, long, long, long);
// multi    multiply matrix at the given dimension
void multi(dym,MatrixXf &);
// get the para slice in the dim
MatrixXf slice(dym,long);
};
*/
// Todo list: read, cycread, write, cycwrite

MatrixXf cif::slice(dym tdim, long index)
{
    if ((tdim-dim)%3==1)
        return data.block<size[(dim-1)%3],size[dim]>(size[(dim-1)%3]*index,0).resize((dim-1)%3],size[dim]);
	rotate(tdim);
    return data.col(index).resize(size[(dim+1)%3],size[(dim+2)%3]);
}
void cif::rotate(dym tdim)
{
    if(tdim==dim) return;
    else if ((tdim-dim)%3==1)
        data.resize(size[tdim],size[dim]*size[(dim-1)%3]).transposeInPlace();
    else
    	data.transposeInPlace().resize(size[(tdim-1)%3]*size[(tdim-2)%3],size[tdim]);
   	dim=tdim;
}
void cif::rotate_top(dym tdim)
{
 	 rotate((tdim-1)%3);
}
void multi(dym tdim,MatrixXf &c)
{
 	 if(dim==tdim)
 	     data=data*c.transpose();
     else
         rotate_top(tdim),data=c*data;
}
void read(char *file, long ncluster, long ncyc, long scluster, long scyc)
{
 }
 
void read(char *file)
{
 	 
}



