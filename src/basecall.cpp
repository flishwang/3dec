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
#include <math.h>
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <stdio.h>

//#define _OPENMP
//       dim=dchn: [chn, cluster, cyc];
//   Eigen:        [cyc,cluster,chn]

//#define _test
//#define _test2


#ifdef _test
int checkblock=7;
int tt=45786;
#endif

class basecall: public Cif
{
public:
    
    MatrixXf neighbour;
    MatrixXf width;
    MatrixXf lambda;
    MatrixXf initcolor;
    Matrixidx idx;
    Map<MatrixXf> zig;
    Map<MatrixXf> zcyc;
    MatrixXf zchn;
    

    MatrixXf zresult;

    int currentiter;
    int itertime;
	int fastmode;

    mylong nblock;
    mylong cnblock;
    mylong cnstart;
    mylong cntotal;
    
    int bufferwidth;
    int largeblock;
    int block;
    int *blockstart;
    int *blocklength;
    int *bufferlength;
    int *bufferstart;
    int *cyclebelong;
    int *cyclepos;
    float stein;
    float wl;
    float colorscale;
    
    MatrixXf *color;
    MatrixXf *phasing;
    MatrixXf previous;
#ifdef _OPENMP
    MatrixXf *color_o;
    MatrixXf *phasing_o;
#endif
    
    basecall(Cif &,int iters, int blocksize, int buffer,mylong maxcluster);
    ~basecall();
    
    int pipe(); //
    int color_kmeans(MatrixXf &a);//
    int fix(int iter,int cb);//
    int inferphasing(int cb);//
    void solveleft(int cb);//
    void fixcolor(int cb);//
    void fixphasing(int cb);//
    void solve(int cb,int calright=1,int writedata=1);//
    void correctneighbour();
    void init();
    
    void minusconstant(int ccyc,MatrixXf &a,int minuspc=1,int minusback=1);
    void buildh(MatrixXf &,Ref<Matrixidx> ,int);
    void multiphasing(MatrixXf &,Ref<MatrixXf> );
    void multicolor(MatrixXf & ,Ref<MatrixXf> );
    
    float myscale(float);
    float mywidth(float);
};


int basecall::color_kmeans(MatrixXf & z)   //z: chn*(m)
{
    const int iters=5;
    const float step=1.25;
    MatrixXf tmp,tmp2;
    z.transposeInPlace();
    int chn=(int)z.rows();
    if(chn!=size[dchn]) return -1;
    MatrixXf z2,z3,tmp3;
    Matrix<bool,Dynamic,Dynamic> tor;
    float zmean;
    MatrixXf cc(4,4),cc0;
    Matrixidx tid;
    if(chn==4)
        cc<<3300,3300,0,0,
        600,2800,0,0,
        0,0,2800,2600,
        0,0,0,3300;
    else
        cc=MatrixXf::Identity(chn,chn);
    cc.transposeInPlace();
    tmp=z.colwise().norm();
    zmean=tmp.mean();
    tor=tmp.array()>zmean/2;
    pick(z,tor);
    pick(tmp,tor);
    
    for(int i=0;i<iters;i++)
    {
        z2.resize(chn,z.cols());
        for(int j=0;j<chn;j++)
        {
            tmp2=(z.array().colwise()*cc.array().col(j)).array().colwise().sum();
            z2.row(j)=tmp2.array()/tmp.array()/cc.col(j).norm();
        }
        purity(z2,tid);
        cc0=cc;
        for(int j=0;j<chn;j++)
        {
            tor=tid.array()==j;
            pick(z,tor,z2);
            cc.col(j)=z2.rowwise().mean();
        }
        if((cc-cc0).norm()<cc0.norm()/10000)
            break;
    }
    for(int j=0;j<chn;j++)
    {
        tor=tid.array()==j;
        pick(z,tor,z2);
        pick(tmp,tor,tmp2);
        
        float thr=0.2f;
        for(int i=0;i<iters/2;i++)
        {
            z3=(z2.array().colwise()*cc.row(j).transpose().array()).array().colwise().sum();
            z3=1-z3.array()/tmp2.array()/cc.row(j).norm();
            tor=z3.array()<thr;
            pick(z2,tor,z3);
            if(z3.cols()>2)
                tmp3=z3.rowwise().mean(),
                cc.col(j)+=step*(tmp3-cc.col(j));
            thr=thr*thr*3;
        }
    }
    initcolor.resize(chn,chn+1);
    initcolor.block(0,0,chn,chn)=cc;
    
    if(chn!=4) initcolor.col(chn).fill(0);
    else
    {
        initcolor.col(chn).segment(0,2)=cc.block(0,2,2,2).rowwise().mean();
        initcolor.col(chn).segment(2,2)=cc.block(2,0,2,2).rowwise().mean();
        for(int k=0;k<4;k++)
            initcolor.col(k)-=initcolor.col(chn);
    }
    
    return 0;
}

int basecall::inferphasing(int cb)
{
    VectorXf vtmp;
    float var_pre;
    float var_pha;
    float var_stop;
    float var_const;
    float var_dim;
    
    phasing[cb]=MatrixXf::Zero(blocklength[cb],bufferlength[cb]+2);
    
    if(cb==0)
        //phasing[cb].diagonal(bufferstart[cb]-blockstart[cb]).fill(1);
        var_pre=var_pha=0.002f,var_stop=0.001f,var_const=0,var_dim=0.997f;
    else
    {
        var_pre=phasing[0](blocklength[0]-1,blocklength[0])/phasing[0](blocklength[0]-1,blocklength[0]-1)/blocklength[0];
        var_pha=phasing[0](blocklength[0]-1,blocklength[0]-2)/phasing[0](blocklength[0]-1,blocklength[0]-1)/blocklength[0];
        var_stop=phasing[cb-1].col(bufferlength[cb-1]).mean();
        var_const=phasing[cb-1].col(bufferlength[cb-1]+1).mean();
        var_dim=phasing[cb-1].row(blocklength[cb-1]-1).head(bufferlength[cb-1]).array().sum()/phasing[cb-1].row(0).head(bufferlength[cb-1]).array().sum();
        var_dim=(var_dim+blocklength[cb-1]-1)/blocklength[cb-1];
        
    }
    
    for(int i=0;i<blocklength[cb];i++)
    {
        int left,right;
        if(i==0&&cb>0)
        {
            left=bufferstart[cb]-bufferstart[cb-1];
            vtmp=phasing[cb-1].row(blocklength[cb-1]-1).segment(left,bufferlength[cb-1]-left);
        }
        else if(i>0)
            vtmp=phasing[cb].row(i-1);
        
        if(i>0||cb>0)
        {
            phasing[cb].row(i).fill(0);
            phasing[cb](i,bufferlength[cb])=var_stop;
            phasing[cb](i,bufferlength[cb]+1)=var_const;
            
            right=mymin((int)vtmp.size(),bufferlength[cb]);
            phasing[cb].row(i).head(right)+=var_dim*var_pha*vtmp.head(right);
            
            right=mymin((int)vtmp.size(),bufferlength[cb]-1);
            phasing[cb].row(i).segment(1,right)+=var_dim*(1-var_pha-var_pre)*vtmp.head(right);
            
            right=mymin((int)vtmp.size(),bufferlength[cb]-2);
            phasing[cb].row(i).segment(2,right)+=var_dim*var_pre*vtmp.head(right);
        }
        else
            phasing[cb](0,0)=1;
        
    }
    return 0;


#ifdef _test
    if(cb==checkblock){
    cout<<"Inferphasing:"<<cb<<endl;
    cout<<phasing[cb].col(bufferlength[cb]).head(9).transpose()<<endl;
    cout<<phasing[cb].col(bufferlength[cb]+1).head(9).transpose()<<endl;
    cout<<phasing[cb].diagonal(blockstart[cb]-bufferstart[cb]).head(15).transpose()<<endl;
    cout<<"Inferphasing end\n";
    }
#endif

}

void basecall::fixcolor(int cb)
{
    MatrixXf tmp;
    MatrixXf zoo=zcyc.block(blockstart[cb],0,blocklength[cb],size[dcluster]*size[dchn]);
    MatrixXf h;
    Matrixidx tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
    minusconstant(blockstart[cb],zoo,1,0);
    buildh(h,tid,cb);
    tmp=phasing[cb].block(0,0,blocklength[cb],bufferlength[cb]);
    multiphasing(h,tmp);
    MatrixXf sh(size[dcluster],size[dchn]+1),szoo(size[dcluster],size[dchn]);
    sh.col(size[dchn]).fill(1);
    
    h.resize(blocklength[cb]*size[dcluster],size[dchn]);
    zoo.resize(blocklength[cb]*size[dcluster],size[dchn]);
    
    for(mylong i=0;i<size[dcluster];i++)
        szoo.row(i)=zoo.row(blocklength[cb]*i+(i%blocklength[cb])),
        sh.row(i).head(size[dchn])=h.row(blocklength[cb]*i+(i%blocklength[cb]));
    sh.array().colwise()*=width.col(cb).array();
    szoo.array().colwise()*=width.col(cb).array();
    float ttw=width.col(cb).array().square().sum();
    float tw=ttw*stein;
    MatrixXf left=sh.transpose()*sh;
    MatrixXf right=sh.transpose()*szoo+color[cb].transpose()*tw;
    left.diagonal().array()+=tw;
    left(size[dchn],size[dchn])+=itertime>1?tw*(itertime-currentiter+1-2*(currentiter>0)):0;
    
    LLT<MatrixXf> llth(left);
#ifndef _OPENMP
    color[cb]=llth.solve(right).transpose();
    color[cb].array().leftCols(size[dchn])*=colorscale/color[cb].diagonal().norm()*sqrt((float)size[dchn]);
#else
    color_o[cb]=llth.solve(right).transpose();
    color_o[cb].array().leftCols(size[dchn])*=colorscale/color_o[cb].diagonal().norm()*sqrt((float)size[dchn]);
#endif
    
#ifdef _test
    if(cb==checkblock)
    cout<<"Fix color:"<<cb<<'\n'<<color[cb]<<endl;
#endif
	if (ttw<1e30f && ttw>-1e30f)
	{
		char tch[100];
		sprintf(tch, "%8.3g", ttw / size[dcluster]);
		cerr << tch;
	}
	else
		cerr << "  NAN:" << cb << ' ';
}
void basecall::fixphasing(int cb)
{
    MatrixXf tmp;
    MatrixXf zoo=zcyc.block(blockstart[cb],0,blocklength[cb],size[dcluster]*size[dchn]);
    MatrixXf h;
    Matrixidx tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
    minusconstant(blockstart[cb],zoo,0,1);
    buildh(h,tid,cb);
    multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
    MatrixXf sh(size[dcluster],bufferlength[cb]+2),szoo(size[dcluster],blocklength[cb]);
    
    h.resize(bufferlength[cb]*size[dcluster],size[dchn]);
    zoo.resize(blocklength[cb]*size[dcluster],size[dchn]);
    
    for(mylong i=0;i<size[dcluster];i++)
        szoo.row(i)=zoo.col(i%size[dchn]).segment(blocklength[cb]*i,blocklength[cb]),
        sh.row(i).head(bufferlength[cb])=h.col(i%size[dchn]).segment(bufferlength[cb]*i,bufferlength[cb]),
        sh(i,bufferlength[cb])=previous.coeff(cb,size[dcluster]*(1%size[dchn])+i);
    sh.col(bufferlength[cb]+1).fill(colorscale);
    
    sh.array().colwise()*=width.col(cb).array();
    szoo.array().colwise()*=width.col(cb).array();
    float tw=width.col(cb).array().square().sum()*stein*colorscale*colorscale;
    MatrixXf left=sh.transpose()*sh;
    MatrixXf right=sh.transpose()*szoo+phasing[cb].transpose()*tw;
    left.diagonal().array()+=tw;
    left(bufferlength[cb]+1,bufferlength[cb]+1)+=itertime>1?tw*(itertime-currentiter+1-2*(currentiter>0)):0;
    
    LLT<MatrixXf> llth(left);
#ifndef _OPENMP
    phasing[cb]=llth.solve(right).transpose();
#else
    phasing_o[cb]=llth.solve(right).transpose();
#endif
    
#ifdef _test
 //   cout<<phasing[cb].col(bufferlength[cb]).head(9).transpose()<<endl;
 //   cout<<phasing[cb].col(bufferlength[cb]+1).head(9).transpose()<<endl;
 //   cout<<phasing[cb].diagonal(blockstart[cb]-bufferstart[cb]).head(9).transpose()<<endl;
    if(cb==checkblock){
    cout<<"Fix phasing:"<<cb<<'\n';
    cout<<phasing[cb].topLeftCorner(9,9)<<endl<<'\n';
    cout<<phasing[cb].bottomRightCorner(9,9)<<endl;
        cout<<endl;}
#endif
}


void basecall::solveleft(int cb)
{
    solve(cb,0,0);
}
void basecall::solve(int cb,int calright,int writedata)
{
    MatrixXf zoo;
    MatrixXf h,tmp;
    Matrixidx tid;
    MatrixXf zco;
    int zoostart=blockstart[cb];
    int zoowidth=blocklength[cb];
    

    if(calright>=0&&blockstart[cb]>bufferstart[cb])
        zoostart-=blockstart[cb]-bufferstart[cb],
        zoowidth+=blockstart[cb]-bufferstart[cb];
    if((calright==1/*||(calright==-1&&cnblock>0)*/)&&blockstart[cb]+blocklength[cb]<bufferstart[cb]+bufferlength[cb])
        zoowidth+=bufferstart[cb]+bufferlength[cb]-(blockstart[cb]+blocklength[cb]);
    zoo=zcyc.middleRows(zoostart,zoowidth);
    
        #ifdef _test2
            int tt=25786;
            int checkblock=1;
        #endif
        #ifdef _test2
            if(cnblock==0&&cb==checkblock){
            tmp.resize(zoowidth,4);
            cout <<"zoowidth,zoo.size "<<zoowidth <<' ' << zoo.size()<<endl;
            tmp<<zoo.col(tt),zoo.col(size[dcluster]+tt),zoo.col(size[dcluster]*2+tt),zoo.col(size[dcluster]*3+tt);
                cout<<"Iter "<<currentiter<<" cluster 45786 block:(Primary)"<<cb<<"\n"<<tmp.topRows(8).transpose()<<endl;}
        #endif
    
    minusconstant(zoostart,zoo);
    
        #ifdef _test2
            if(cnblock==0&&cb==checkblock){
                tmp.resize(zoowidth,4);
                tmp<<zoo.col(tt),zoo.col(size[dcluster]+tt),zoo.col(size[dcluster]*2+tt),zoo.col(size[dcluster]*3+tt);
                cout<<"Iter "<<currentiter<<" cluster 45786 block:(MinusConstant)"<<cb<<"\n"<<tmp.topRows(8).transpose()<<endl;}
        #endif
    
    zco=zoo.middleRows(blockstart[cb]-zoostart,blocklength[cb]);
    
    tmp=color[cyclebelong[zoostart]].leftCols(size[dchn]).inverse();
    for(int i=0;i<zoowidth;i++)
    {
        if(i>0&&cyclebelong[zoostart+i]!=cyclebelong[zoostart+i-1])
            tmp=color[cyclebelong[zoostart+i]].leftCols(size[dchn]).inverse();
        h=zoo.row(i);
        multicolor(h,tmp);
        zoo.row(i)=h;
    }
    
        #ifdef _test2
            if(cnblock==0&&cb==checkblock){
                tmp.resize(zoowidth,4);
                tmp<<zoo.col(tt),zoo.col(size[dcluster]+tt),zoo.col(size[dcluster]*2+tt),zoo.col(size[dcluster]*3+tt);
                cout<<"Iter "<<currentiter<<" cluster 45786 block:(Corret Color)"<<cb<<"\n"<<tmp.topRows(8).transpose()<<endl;}
        #endif
            
    tmp=MatrixXf::Zero(zoowidth,zoowidth);
    if(calright>=0)
    {
        if(blockstart[cb]>bufferstart[cb])
        {
            int cw=bufferstart[cb-1]+bufferlength[cb-1]-bufferstart[cb];
            int rw=blockstart[cb]-bufferstart[cb];
            tmp.topLeftCorner(rw,cw)=phasing[cb-1].block(blocklength[cb-1]-rw,bufferlength[cb-1]-cw,rw,cw);
            tmp.block(rw,0,blocklength[cb],zoowidth)=phasing[cb].leftCols(zoowidth);
        }
        else
            tmp.topLeftCorner(blocklength[cb],zoowidth)=phasing[cb].leftCols(zoowidth);
        
        if(calright==1&&blockstart[cb]+blocklength[cb]<bufferstart[cb]+bufferlength[cb])
        {
            int rw=bufferstart[cb]+bufferlength[cb]-blockstart[cb]-blocklength[cb];
            int cw=rw+blockstart[cb+1]-bufferstart[cb+1];
            tmp.bottomRightCorner(rw,cw)=phasing[cb+1].topLeftCorner(rw,cw);
        }
    }
    else
        tmp=phasing[cb].middleCols(blockstart[cb]-bufferstart[cb], blocklength[cb]);
    
        #ifdef _test
            if(cnblock==0&&cb==checkblock)
            cout<<"TMP MATRIX:"<<cb<<" width:"<<zoowidth<<endl<<tmp.topLeftCorner(8,8)<<"\n\n"<<tmp.bottomRightCorner(8,8).leftCols(8)<<endl;
        #endif
    
    tmp=tmp.inverse();
    
        #ifdef _test
            if(cnblock==0&&cb==checkblock)
            cout<<"TMP INVERSE MATRIX:"<<cb<<" width:"<<zoowidth<<endl<<tmp.topLeftCorner(8,8)<<"\n\n"<<tmp.bottomRightCorner(8,8).leftCols(8)<<endl;
        #endif
    
    zoo=tmp.middleRows(blockstart[cb]-zoostart,blocklength[cb])*zoo;
    
    
        #ifdef _test2
            if(cnblock==0&&cb==checkblock){
                tmp.resize(blocklength[cb],4);
                tmp<<zoo.col(tt),zoo.col(size[dcluster]+tt),zoo.col(size[dcluster]*2+tt),zoo.col(size[dcluster]*3+tt);
                cout<<"Iter "<<currentiter<<" cluster 45786 block:(Correct phasing)"<<cb<<" delay "<<blockstart[cb]-zoostart<<"\n"<<tmp.topRows(8).transpose()<<endl;}
        #endif
    
    if(currentiter==itertime-1&&writedata)
        zresult.block(blockstart[cb],0,blocklength[cb],size[dchn]*size[dcluster])=zoo;

    
    zoo.resize(size[dcluster]*blocklength[cb],size[dchn]);
    
    purity(zoo,tid);
    tid.resize(blocklength[cb],size[dcluster]);
    
#ifndef _OPENMP
    idx.block(blockstart[cb],0,blocklength[cb],size[dcluster])=tid;
#else
#pragma omp critical(idx)
    {
        idx.block(blockstart[cb],0,blocklength[cb],size[dcluster])=tid;
    }
#endif
    
        #ifdef _test
            if(cnblock==0&&cb==checkblock){
                cout<<"Idx(tt)"<<idx.cast<int>().col(tt).segment(blockstart[cb],blocklength[cb]).transpose()<<endl;}
        #endif
    
    //width,lambda
    if((!cb&&!fastmode)||!writedata||currentiter<itertime-1)
    {
        float meanscale=color[cb].maxCoeff();
        meanscale*=meanscale;
        int ww=blockstart[cb]-bufferstart[cb]+blocklength[cb];

#ifdef _OPENMP
#pragma omp critical(idx)
#endif
        {
            if(calright==1)
                tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
            else if(calright==0)
                tid=idx.block(bufferstart[cb],0,ww,size[dcluster]);
            else
                tid=idx.block(blockstart[cb],0,blocklength[cb],size[dcluster]);
        }
/*#else
        if(calright==1)
            tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
        else if(calright==0)
            tid=idx.block(bufferstart[cb],0,ww,size[dcluster]);
        else
            tid=idx.block(blockstart[cb],0,blocklength[cb],size[dcluster]);
#endif*/
        
            #ifdef _test2
                    if(cnblock==0&&cb==checkblock){
                        tmp.resize(blocklength[cb],1);
                        tmp=tid.col(tt).cast<float>();
                        cout<<"Iter "<<currentiter<<" cluster 45786 block:(IDX)"<<cb<<"\n"<<tmp.transpose()<<endl;}
            #endif
        
        if(calright==1)
        {
            buildh(h,tid,cb);
            multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
            multiphasing(h,phasing[cb].block(0,0,blocklength[cb],bufferlength[cb]));
        }
        else if(calright==0)
        {
            buildh(h,tid,cb);
            multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
            multiphasing(h,phasing[cb].block(0,0,blocklength[cb],ww));
        }
        else
        {
            buildh(h,tid,cb);
            multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
            multiphasing(h,phasing[cb].middleCols(blockstart[cb]-bufferstart[cb],blocklength[cb]));
        }
        
        int cycs=blocklength[cb];
        h.resize(cycs*size[dcluster],size[dchn]);
        zco.resize(cycs*size[dcluster],size[dchn]);
        
        for(mylong i=0;i<size[dcluster];i++)
        {
            Ref<MatrixXf> h0=h.block(i*cycs,0,cycs,size[dchn]);
            Ref<MatrixXf> h1=zco.block(cycs*i,0,cycs,size[dchn]);
            
                #ifdef _test2
                    if(i==tt&&cnblock==0&&cb==checkblock)
                        cout<<"H0 AND H1\n"<<h0.topRows(9).transpose()<<endl<<h1.topRows(9).transpose()<<endl;
                #endif
            
            float scale=(h0.array()*h1.array()).sum()/h0.array().square().sum();
            scale=scale<0.2f?0.2f:scale;
            scale=scale>5?5:scale;
            width(i,cb)=mywidth((h0*scale-h1).array().square().sum()/blocklength[cb]/size[dchn]/meanscale);
            lambda(i,cb)*=scale; 
        }
        
        lambda.array().col(cb)/=lambda.col(cb).mean();
    }
    
        #ifdef _test
            cout<<"Lambda(tt):"<<lambda(tt,cb)<<endl;
        #endif

}

void basecall::correctneighbour()
{
    MatrixXf nei(size[dchn],size[dchn]);
    MatrixXb tor;
    
    MatrixXf m1=zcyc.row(1);
    VectorXf mb;
    
    m1.resize(size[dcluster],size[dchn]);

    
    mb=m1.rowwise().maxCoeff();
    
    mb=mb.array();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int kk=0;kk<size[dchn]*size[dchn];kk++)
    {
        int i=kk/size[dchn],j=kk%size[dchn];
        float mymean=mb.mean();
        float mysum=mymean,mywidth=1;
        for(int k=0;k<size[dcluster];k++)
            if(idx.coeff(0,k)==i&&idx.coeff(1,k)==j&&mb.coeff(k)>mymean/2&&mb.coeff(k)<mymean*2)
                mysum+=mb.coeff(k),mywidth+=lambda(k,0);
        nei(i,j)=mysum/mywidth;
    }
    nei=0.1f+0.9f*nei.array();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=1;i<size[dcyc];i++)
        for(mylong j=0;j<size[dcluster];j++)
            zchn.array().row(size[dcyc]*j+i)/=nei.array().row(idx(i-1,j));
    
}

void basecall::buildh(MatrixXf &h,Ref<Matrixidx> pid,int cb)
{
    mylong cycs=(mylong)pid.rows();
    h=MatrixXf::Zero(cycs*size[dcluster],size[dchn]);
    for(int i=0;i<cycs;i++)
        for(mylong j=0;j<size[dcluster];j++)
            h(i+j*cycs,pid(i,j))=lambda.coeff(j,cb);
    h.resize(cycs,size[dcluster]*size[dchn]);
}


void basecall::multiphasing(MatrixXf &zoo,Ref<MatrixXf> p)
{
    zoo=p*zoo;
}

void basecall::multicolor(MatrixXf &zoo,Ref<MatrixXf> c)
{
    mylong cycs=(mylong)zoo.rows();
    zoo.resize(cycs*size[dcluster],size[dchn]);
    zoo=zoo*c.transpose();
    zoo.resize(cycs,size[dchn]*size[dcluster]);
}

void basecall::minusconstant(int ccyc,MatrixXf &a,int minuspc,int minusback)
{
    mylong cycles=(mylong)a.rows();
    
    for(int i=0;i<cycles;i++)
    {
        int ccb=cyclebelong[ccyc+i];
        int ccc=cyclepos[ccyc+i];
        if(minuspc)
            a.array().row(i)-=phasing[ccb].coeff(ccc,bufferlength[ccb])*previous.array().row(ccb)+phasing[ccb].coeff(ccc,bufferlength[ccb]+1)*colorscale;
        if(minusback)
            for(int j=0;j<size[dchn];j++)
                a.array().row(i).segment(size[dcluster]*j,size[dcluster])-=color[ccb].coeff(j,size[dchn]);
    }
    
}

float basecall::mywidth(float n)
{
    return exp(-n*wl);
}


#ifndef _OPENMP

int basecall::fix(int iter,int cb)
{

    currentiter=iter;
    
    if(iter==0)
    {
        if(cnblock==0)
            inferphasing(cb);
        if(cb==0)
            solveleft(cb);
        if(cb<block-1)
        {
            if(cnblock==0)
                inferphasing(cb+1);
            solveleft(cb+1);
        }
    }
    fixphasing(cb);
    fixcolor(cb);
    
    solve(cb);
    return 0;
}

int basecall::pipe()
{

    for(int k=0;k<nblock;k++)
    {
        cnblock=k;
        init();
        if(k) itertime=1;
        
        for(int i=0;i<itertime;i++,stein=stein>0.02?stein-0.0162:stein/1.5,wl=10*((i+3.0)))
        {
            cerr<<"\nBasecalling nblock "<<k+1<<", iteration "<<i+1<<endl;
            for(int j=0;j<block;j++)
                fix(i,j);
        }
        
        zcyc=zresult;
    //    correctneighbour();
        data.block(cnstart*size[dcyc],0,size[dcluster]*size[dcyc],size[dchn])=zchn;
    }

    size[dcluster]=cntotal;
    
    return 0;
}
#else

int basecall::pipe()
{
	if (fastmode)
		cerr << "Some steps are omitted to save calculation time.\n";
    for(int k=0;k<nblock;k++)
    {
        cnblock=k;
        init();
        
        if(cnblock) itertime=1;
        currentiter=0;
        if((cnblock<4 &&!fastmode) || cnblock<2)
            itertime++;
  /*
        if(cnblock==0)
        {
            inferphasing(0);
            solveleft(0);
            if(block>1)
                inferphasing(1),solveleft(1);
            fixphasing(0);
            fixcolor(0);
            phasing[0]=phasing_o[0];
            color[0]=color_o[0];
            for(int i=1;i<block;i++)
                inferphasing(i);
        }
   */
        if(cnblock==0)
            for(int i=0;i<block;i++)
                inferphasing(i);
        
        
        for(int i=0;i<itertime;i++)
        {
            currentiter=i;
            
            cerr<<"\nBasecalling nblock "<<cnblock+1<<", iteration "<<i+1<<endl;
            
		if (!fastmode || itertime>1 || !(cnblock%4))
		{

			if (currentiter == 0)
#pragma omp parallel for
				for (int j = 0; j < block; j++)
					solve(j, -1, 0);
#pragma omp parallel for
			for (int j = 0; j < block; j++)
				fixphasing(j);
#pragma omp parallel for
			for (int j = 0; j < block; j++)
				phasing[j] = phasing_o[j];
#pragma omp parallel for
			for (int j = 0; j < block; j++)
				fixcolor(j);
#pragma omp parallel for
			for (int j = 0; j < block; j++)
				color[j] = color_o[j];
		}
		else
			cerr << "Skip re-estimation.";

#pragma omp parallel for
            for(int j=0;j<block;j++)
                solve(j);
            stein=stein>0.02f?stein-0.0162f:stein/1.5f,wl=10*((i+4.0f));
        }

        zcyc=zresult;
		if (!fastmode)
			correctneighbour();
        data.block(cnstart*size[dcyc],0,size[dcluster]*size[dcyc],size[dchn])=zchn;
    }
    
    size[dcluster]=cntotal;
    return 0;
}
#endif


void basecall::init()
{
    
    cnstart=(mylong)floor((float)cntotal/nblock*cnblock+0.5f);
	size[dcluster] = (mylong)floor(((float)cntotal / nblock*(1 + cnblock)) + 0.5f) - cnstart;
    zchn=zig.block(size[dcyc]*cnstart,0,size[dcyc]*size[dcluster],size[dchn]);
    new (& zcyc) Map<MatrixXf>(zchn.data(),size[dcyc],size[dcluster]*size[dchn]);
    zresult.resizeLike(zcyc);

    if(cnblock==0)
    {
        MatrixXf tmp;
        tmp=zcyc.row(0);
        tmp.resize(size[dcluster],size[dchn]);
        color_kmeans(tmp);
        
        for(int i=0;i<block;i++)
            color[i]=initcolor;
        
        stein=.05f,wl=30;
    }
    else
        stein=.01f,wl=40;
    
    
    previous.resize(block,size[dcluster]*size[dchn]);
    colorscale=color[0].diagonal().norm()/sqrt((float)size[dchn]);
    
    width=MatrixXf::Ones(size[dcluster],block);
    lambda=MatrixXf::Ones(size[dcluster],block);
    idx.resize(size[dcyc],size[dcluster]);
    
    for(int i=0;i<block;i++)
    {
        if(i==0)
            previous.row(i).fill(0);
        else
        {
            previous.row(i)=previous.row(i-1);
            for(int j=bufferstart[i-1];j<bufferstart[i];j++)
            {
                previous.row(i)+=zcyc.row(j);
            }
        }
    }
    for(int i=1;i<block;i++)
        previous.array().row(i)/=(float)(bufferstart[i]>5?bufferstart[i]:5);
}



basecall::basecall(Cif &z,int iters, int blocksize, int buffer,mylong maxcluster):
    zig(z.data.data(),z.data.rows(),z.data.cols()),
    zcyc(z.data.data(),z.size[dcyc],z.size[dcluster]*z.size[dchn])
{
    copycif(z,0);
    data.resize(size[dcyc]*size[dcluster],size[dchn]);
    cntotal=size[dcluster];
	if (iters < 0)
		itertime = -iters, fastmode = 1;
	else
		itertime = iters, fastmode = 0;

    bufferwidth=buffer;
    
    block=(int)ceil((float)size[dcyc]/blocksize);
    if (block<1) block=1;
    
	if (maxcluster == 0) //default block numbers
		maxcluster = 625 * (buffer * 2 + blocksize + 2)*block;
	if (maxcluster<0)
        maxcluster=-size[dcluster]/maxcluster+2;
    
    nblock=(mylong)ceil((float)size[dcluster]/maxcluster);
    cerr<<"Clusters are divided into "<<nblock<<" blocks.\n";
    
    blocksize=size[dcyc]/block;
    
    blockstart=new int[block];
    blocklength=new int[block];
    bufferlength=new int[block];
    bufferstart=new int[block];
    cyclebelong=new int[size[dcyc]];
    cyclepos=new int[size[dcyc]];
    largeblock=size[dcyc]%block;
    
    color=new MatrixXf [block];
    phasing=new MatrixXf [block];
#ifdef _OPENMP
    color_o=new MatrixXf [block];
    phasing_o=new MatrixXf [block];
#endif
    
    blockstart[0]=0;
    bufferstart[0]=0;
    for(int i=0;i<block;i++)
    {
        if(i>0)
            blockstart[i]=blockstart[i-1]+blocklength[i-1],
            bufferstart[i]=mymax<int>(0,blockstart[i]-buffer);
        
        blocklength[i]= i<largeblock?blocksize+1:blocksize;
        bufferlength[i]= mymin<int>(blockstart[i]+blocklength[i]+buffer,size[dcyc])-bufferstart[i];
        
        for(int j=0;j<blocklength[i];j++)
            cyclebelong[blockstart[i]+j]=i,
            cyclepos[blockstart[i]+j]=j;
    }
    
}
basecall::~basecall()
{
    delete []blockstart;
    delete []blocklength;
    delete []bufferlength;
    delete []bufferstart;
    delete []cyclebelong;
    for(int i=0;i<block;i++)
        color[i].resize(0,0),phasing[i].resize(0,0);
    
#ifdef _OPENMP
    for(int i=0;i<block;i++)
        color_o[i].resize(0,0),phasing_o[i].resize(0,0);
#endif
    
}


int Cif::basecalling(Cif &z,int iters, int blocksize, int buffer,mylong maxcluster)
{
 //   chnpur=size[dchn];
    z.rotate(dchn);
    basecall ba(z,iters,blocksize,buffer,maxcluster);
    ba.pipe();
    
    copycif(z,0);
    data=ba.data;
    cerr<<"\nBasecalling done.\n";
    return 1;
}




