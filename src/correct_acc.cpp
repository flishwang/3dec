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

#include <stdio.h>
#include <math.h>
#include <vector>

using namespace::std;
using namespace::Eigen;

class caccenv:public Cif
{
public:
    float *startx;
    MatrixXi loc;
    VectorXf width;
    float thr;
    VectorXf qd;
    int blocksize;
    int buffersize;
    int block;
    mylong minx;
    mylong miny;
    mylong maxx;
    mylong maxy;
    Map<MatrixXf> zcluster;
    Map<MatrixXf> zzcluster;
    mylong cchn;
    mylong ccyc;
    MatrixXb *inblock;
    MatrixXl *inbuffer;
    int pipe();
    caccenv(Cif &,Cif &,MatrixXf & ,float,float,float);
    ~caccenv();
    
    float findss(mylong , mylong ,float);
    float f(VectorXf &);
};

caccenv::caccenv(Cif &z,Cif &zz,MatrixXf &ll,float ks,float blocks,float buffers):
zcluster(z.data.data(),z.size[dcluster],z.size[dchn]*z.size[dcyc]),
zzcluster(zz.data.data(),zz.size[dcluster],zz.size[dchn]*zz.size[dcyc])
{
    VectorXf ww;
    float scale=0;
    copycif(z,0);
    if(size[dcluster]!=zz.size[dcluster])
        cerr<<"Data for correcting ACC not suitable, critical error may happen!\n";
    
    cerr<<"Preparing for ACC correction...\n";
    blocksize=(int)blocks;
    buffersize=(int)buffers;
    cchn=zz.size[dchn];
    ccyc=zz.size[dcyc];

    minx=(mylong)ll.col(0).minCoeff();
    miny=(mylong)ll.col(1).minCoeff();
    maxx=(mylong)ll.col(0).maxCoeff();
    maxy=(mylong)ll.col(1).maxCoeff();

    loc.resizeLike(ll);
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(mylong i=0;i<ll.rows();i++)
        loc(i,0)=(int)floor(ll.coeff(i,0)-minx+0.5f),
        loc(i,1)=(int)floor(ll.coeff(i,1)-miny+0.5f);
    maxx-=minx-1;
    maxy-=miny-2;
    minx=0;
    miny=0;
    
    if(maxx<maxy)
    {
        mylong z;
        z=maxx,maxx=maxy,maxy=z;
        loc.col(0).swap(loc.col(1));
    }
    
    width=VectorXf::Zero(ccyc);
    
    block=(int)ceil((float)(maxx-minx)/blocksize);
    startx=new float[block];
    inblock=new MatrixXb[block];
    inbuffer=new MatrixXl[block];
    
#ifdef _OPENMP
  
#pragma omp parallel for
#endif
    for(int i=0;i<block;i++)
    {
        startx[i]=(float)(minx+i*blocksize-buffersize);
        mylong n=0;
        float xloc;
        MatrixXl tmp(size[dcluster],1);
        for(mylong j=0;j<size[dcluster];j++)
            if(loc.coeff(j,0)>=startx[i]-buffersize&&loc.coeff(j,0)<startx[i]+blocksize+buffersize)
                tmp(n,0)=j,n++;
        inbuffer[i]=tmp.topRows(n);
        
        inblock[i].resize(n,1);
        for(mylong j=0;j<n;j++)
            xloc=loc(inbuffer[i].coeff(j,0),0),inblock[i](j,0)=xloc>=startx[i]&&xloc<startx[i]+blocksize;
    }
    
    qd.resize(size[dcluster]);
    mylong mid=(cchn-2)/2;
    
#ifdef _OPENMP
    int mt=omp_get_max_threads();
    MatrixXf twidth=MatrixXf::Zero(ccyc,mt);
    MatrixXf tscale=VectorXf::Zero(mt);
#pragma omp parallel for
#endif
    for(mylong i=0;i<size[dcluster];i++)
    {
        float dqi=1000000;
        float dmi=0;
        for(mylong j=0;j<ccyc;j++)
        {
            VectorXf tmp=zzcluster.array().row(i).segment(cchn*j,cchn).abs();
            float ttt;
            for(int t1=0;t1<cchn-1;t1++)
                for(int t2=t1+1;t2<cchn;t2++)
                    if(tmp.coeff(t1)>tmp.coeff(t2))
                        ttt=tmp.coeff(t1),tmp(t1)=tmp.coeff(t2),tmp(t2)=ttt;
#ifndef _OPENMP
            {
                width(j)+=tmp.coeff(mid);
                scale+=tmp.coeff(cchn-1);
            }
#else
            {
                twidth(j,omp_get_thread_num())+=tmp.coeff(mid);
                tscale(omp_get_thread_num())+=tmp.coeff(cchn-1);
            }
#endif
            if(tmp.coeff(cchn-1)>dmi)
                dmi=tmp.coeff(cchn-1);
            if(tmp.coeff(cchn-1)-tmp.coeff(cchn-2)<dqi)
                dqi=tmp.coeff(cchn-1)-tmp.coeff(cchn-2);
        }
        qd(i)=dqi/dmi;
    }
#ifdef _OPENMP
    width=twidth.rowwise().sum();
    scale=tscale.sum();
#endif
    
    width=width.array().inverse()*size[dcluster];
    thr=(-.75f+1.5f*ks)*width.sum()*(scale/size[dcluster]/ccyc);
}



int caccenv::pipe()
{
    data.resizeLike(zcluster);
    cerr<<"Clusters are divided into "<<block<<" blocks:\n";
    for(int i=0;i<block;i++)
        cerr<<'=';
    cerr<<endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int k=0;k<block;k++)
    {
        cerr<<'.';
        typedef Triplet<float> myt;
        vector<myt> list;
        list.reserve(inbuffer[k].rows()*3);
        MatrixXl mymap=MatrixXl::Zero(blocksize+2*buffersize+2,maxy+2);
        for(mylong i=0;i<inbuffer[k].rows();i++)
        {
            mylong t=inbuffer[k].coeff(i,0);
            mylong bx=loc.coeff(t,0)-startx[k]+buffersize;
            mylong by=loc.coeff(t,1);
            while(mymap(bx,by))
                by++;
            mymap(bx,by)=i;
        }
        for(mylong i=0;i<inbuffer[k].rows();i++)
        {
            list.push_back(myt(i,i,1));
            mylong i1=inbuffer[k].coeff(i,0);
            mylong j;
            if(qd(i1)<0.7)
                for(mylong j1=loc.coeff(i1,0)-startx[k]+buffersize-4;j1<=loc.coeff(i1,0)-startx[k]+buffersize+4&&j1<blocksize+2*buffersize;j1++)
                    for(mylong j2=loc.coeff(i1,1)-4;j2<=loc.coeff(i1,1)+4&&j2<maxy;j2++)
                        if(j1>=0&&j2>=0&&(j=mymap(j1,j2))&&j!=i)
                        {
                            float xx=findss(i1,inbuffer[k].coeff(j,0),thr);
			//				float xx = findss(i1, inbuffer[k].coeff(j, 0), 0);
			//				if (xx > 0)
			//					xx = findss(i1, inbuffer[k].coeff(j, 0), thr / 0.3 / ccyc/xx);

                            if(xx>0)
                                list.push_back(myt(i,j,xx));
                        }
            
        }
        SparseMatrix<float> acc((mylong)inbuffer[k].rows(),(mylong)inbuffer[k].rows());
        acc.setFromTriplets(list.begin(),list.end());
        SparseLU<SparseMatrix<float> > solver;
        solver.analyzePattern(acc);
        solver.factorize(acc);
        MatrixXf tmp;
        tmp.resize(inbuffer[k].rows(),size[dchn]*size[dcyc]);
        for(mylong i=0;i<inbuffer[k].rows();i++)
            tmp.row(i)=zcluster.row(inbuffer[k].coeff(i,0));
        tmp=solver.solve(tmp);
        for(mylong i=0;i<inbuffer[k].rows();i++)
            if(inblock[k].coeff(i,0))
                data.row(inbuffer[k].coeff(i,0))=tmp.row(i);
    }
    cerr<<"\nDone.\n";
    return 0;
}





caccenv::~caccenv()
{
    for(int i=0;i<block;i++)
        inblock[i].resize(0,0),
        inbuffer[i].resize(0,0);
    delete []startx;
	delete[]inblock;
	delete[]inbuffer;
}




float caccenv::f(VectorXf & a)
{
    float y=0;
    float t;
    for(int j=0;j<ccyc;j++)
    {
        t=0;
        float u=-1000000;
        for(int k=0;k<cchn;k++)
        {
            float s=a.coeff(k+j*cchn);
            t+=abs(s);
            if(s>u)
                u=s;
        }
        y+=width.coeff(j)*(t-abs(u));
    }
    return y;
}
float caccenv::findss(mylong i1,mylong i2,float thr)
{
    float left=0;
    float right=0.98f;
    float alpha=0.001f;
    float middle;
    float ftmp;
    VectorXf x=zzcluster.row(i1);
    VectorXf y=zzcluster.row(i2);
    VectorXf tmp;
    while(right-left>0.001f)
    {
        middle=left*0.618f+right*0.382f;
        tmp=x-middle*y;
        ftmp=f(tmp);
        tmp=x-(alpha+middle)*y;
        if(ftmp-f(tmp)>thr*alpha)
            left=middle;
        else
            right=middle;
    }
    return left;
}


int correct_acc(Cif &z,Cif &zz,MatrixXf &loc,float ks,float blocksize,float buffersize)
{
    z.rotate_top(dcluster);
    zz.rotate_top(dcluster);
    MatrixXf tmp1;
 
    caccenv env(z,zz,loc,ks,blocksize,buffersize);
    env.pipe();
   
    z.data=env.data;
    
    return 0;
}