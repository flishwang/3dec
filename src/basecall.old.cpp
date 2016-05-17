#include "classtype.h"
#include "cacc_funs.h"
#include <math.h>
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>

//       dim=dchn: [chn, cluster, cyc];
//   Eigen:        [cyc,cluster,chn]


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
    
    int currentiter;
    int itertime;
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
    float stein;
    float wl;
    float colorscale;
    
    MatrixXf *color;
    MatrixXf *phasing;
    MatrixXf previous;
    
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
    
    void minusconstant(int cb,MatrixXf &a,int minusprevious=1,int minusback=1);
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
        
        float thr=0.2;
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
    if(cb>=block)
        return 0;
    
    phasing[cb]=MatrixXf::Zero(blocklength[cb],bufferlength[cb]+2);
    
    if(cb==0)
        phasing[cb].diagonal(bufferstart[cb]-blockstart[cb]).fill(1);
    else
    {
        float var_pre=phasing[0](blocklength[0]-1,blocklength[0])/phasing[0](blocklength[0]-1,blocklength[0]-1)/blocklength[0];
        float var_pha=phasing[0](blocklength[0]-1,blocklength[0]-2)/phasing[0](blocklength[0]-1,blocklength[0]-1)/blocklength[0];
        float var_stop=phasing[cb-1].col(bufferlength[cb-1]).mean();
        float var_const=phasing[cb-1].col(bufferlength[cb-1]+1).mean();
        float var_dim=sqrt(phasing[cb-1].row(2).head(bufferlength[cb-1]).array().sum()/phasing[cb-1].row(0).head(bufferlength[cb-1]).array().sum());
        

        
        for(int i=0;i<blocklength[cb];i++)
        {
            int left,right;
            if(i==0)
            {
                left=blocklength[cb-1]+blockstart[cb-1]-bufferstart[cb-1]-bufferwidth;
                vtmp=phasing[cb-1].row(blocklength[cb-1]-1).segment(left,bufferlength[cb-1]-left);
            }
            else
                vtmp=phasing[cb].row(i-1);
            
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
    }
    
#ifdef _test
    cout<<"Inferphasing:"<<cb<<endl;
    cout<<phasing[cb].col(bufferlength[cb]).head(9).transpose()<<endl;
    cout<<phasing[cb].col(bufferlength[cb]+1).head(9).transpose()<<endl;
    cout<<phasing[cb].diagonal(blockstart[cb]-bufferstart[cb]).head(9).transpose()<<endl;
#endif
    
    return 0;
}


void basecall::fixcolor(int cb)
{
    MatrixXf tmp;
    MatrixXf zoo=zcyc.block(blockstart[cb],0,blocklength[cb],size[dcluster]*size[dchn]);
    MatrixXf h;
    Matrixidx tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
    minusconstant(cb,zoo,1,0);
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
    left(size[dchn],size[dchn])+=tw*(itertime-currentiter+1-2*(currentiter>0));
    
    LLT<MatrixXf> llth(left);
    color[cb]=llth.solve(right).transpose();
    
#ifdef _test
    cout<<color[cb]<<endl<<endl;
#endif
    
    cerr<<"Basecalling nblock "<<cnblock<<", iteration "<<currentiter<<", block "<<cb;
    cerr<<", WidthRatio:"<<ttw/size[dcluster]<<"\n";
}
void basecall::fixphasing(int cb)
{
    MatrixXf tmp;
    MatrixXf zoo=zcyc.block(blockstart[cb],0,blocklength[cb],size[dcluster]*size[dchn]);
    MatrixXf h;
    Matrixidx tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
    minusconstant(cb,zoo,0,1);
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
    left(bufferlength[cb]+1,bufferlength[cb]+1)+=tw*(itertime-currentiter+1-2*(currentiter>0));
    
    
    LLT<MatrixXf> llth(left);
    phasing[cb]=llth.solve(right).transpose();
    
#ifdef _test
    cout<<phasing[cb].col(bufferlength[cb]).head(9).transpose()<<endl;
    cout<<phasing[cb].col(bufferlength[cb]+1).head(9).transpose()<<endl;
    cout<<phasing[cb].diagonal(blockstart[cb]-bufferstart[cb]).head(9).transpose()<<endl;
#endif
}


void basecall::solveleft(int cb)
{
    solve(cb,0,0);
}
void basecall::solve(int cb,int calright,int writedata)
{
    MatrixXf zoo=zcyc.block(blockstart[cb],0,blocklength[cb],size[dcluster]*size[dchn]);
    MatrixXf h,tmp;
    Matrixidx tid;
    minusconstant(cb,zoo);
    MatrixXf zco=zoo;
    
    float meanscale=color[cb].maxCoeff();
    meanscale*=meanscale;
    
    if(blockstart[cb]>bufferstart[cb])
    {
        int ww=blockstart[cb]-bufferstart[cb];
        Matrixidx leftidx=idx.block(bufferstart[cb],0,ww,size[dcluster]);
        buildh(h,leftidx,cb);
        multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
        multiphasing(h,phasing[cb].block(0,0,blocklength[cb],ww));
        zoo-=h;
    }
    if(calright && blockstart[cb]+blocklength[cb]<bufferstart[cb]+bufferlength[cb])
    {
        int ww=bufferstart[cb]+bufferlength[cb]-(blockstart[cb]+blocklength[cb]);
        Matrixidx rightidx=idx.block(blockstart[cb]+blocklength[cb],0,ww,size[dcluster]);
        buildh(h,rightidx,cb);
        multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
        multiphasing(h,phasing[cb].block(0,blockstart[cb]-bufferstart[cb]+blocklength[cb],blocklength[cb],ww));
        zoo-=h;
    }
    zoo=phasing[cb].block(0,blockstart[cb]-bufferstart[cb],blocklength[cb],blocklength[cb]).inverse()*zoo;
    tmp=color[cb].block(0,0,size[dchn],size[dchn]).inverse();
    multicolor(zoo,tmp);
    
    if(currentiter==itertime-1&&writedata)
        zcyc.block(blockstart[cb],0,blocklength[cb],size[dchn]*size[dcluster])=zoo;
    
    zoo.resize(size[dcluster]*blocklength[cb],size[dchn]);
    
    purity(zoo,tid);
    tid.resize(blocklength[cb],size[dcluster]);
    
    idx.block(blockstart[cb],0,blocklength[cb],size[dcluster])=tid;
    
#ifdef _test
    cout<<"Idx(100)"<<idx.cast<int>().col(100).segment(blockstart[cb],blocklength[cb]).transpose()<<endl;
#endif
    
    //width,lambda
    if(!cb||!writedata||currentiter<itertime-1)
    {
        if(calright)
        {
            tid=idx.block(bufferstart[cb],0,bufferlength[cb],size[dcluster]);
            buildh(h,tid,cb);
            multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
            multiphasing(h,phasing[cb].block(0,0,blocklength[cb],bufferlength[cb]));
        }
        else
        {
            int ww=blockstart[cb]-bufferstart[cb]+blocklength[cb];
            tid=idx.block(bufferstart[cb],0,ww,size[dcluster]);
            buildh(h,tid,cb);
            multicolor(h,color[cb].block(0,0,size[dchn],size[dchn]));
            multiphasing(h,phasing[cb].block(0,0,blocklength[cb],ww));
        }
        
        int cycs=blocklength[cb];
        h.resize(cycs*size[dcluster],size[dchn]);
        zco.resize(cycs*size[dcluster],size[dchn]);
        
        for(mylong i=0;i<size[dcluster];i++)
        {
            Ref<MatrixXf> h0=h.block(i*cycs,0,cycs,size[dchn]);
            Ref<MatrixXf> h1=zco.block(cycs*i,0,cycs,size[dchn]);
            
#ifdef _test
            if(i==100)
                cout<<"H0 AND H1\n"<<h0<<endl<<h1<<endl;
#endif
            
            float scale=(h0.array()*h1.array()).sum()/h0.array().square().sum();
            scale=scale<0.2?0.2:scale;
            scale=scale>5?5:scale;
            width(i,cb)=mywidth((h0*scale-h1).array().square().sum()/blocklength[cb]/size[dchn]/meanscale);
            lambda(i,cb)*=scale;
        }
        
        lambda.array().col(cb)/=lambda.col(cb).mean();
    }
    
#ifdef _test
    cout<<"Lambda(100):"<<lambda(100,cb)<<endl;
    
#endif

}

void basecall::correctneighbour()
{
    MatrixXf nei(size[dchn],size[dchn]);
    MatrixXb tor;
    

    MatrixXf m0=zcyc.row(0);
    MatrixXf m1=zcyc.row(1);
    VectorXf ma,mb;
    m0.resize(size[dcluster],size[dchn]);
    m1.resize(size[dcluster],size[dchn]);

    ma=m0.rowwise().maxCoeff();
    mb=m1.rowwise().maxCoeff();
    
    m0=mb.array()/ma.array();
    for(int i=0;i<size[dchn];i++)
        for(int j=0;j<size[dchn];j++)
        {
            tor=((idx.row(0).array()==i) && (idx.row(1).array()==j) && (ma.transpose().array()>(ma.mean()/2))).cast<bool>().transpose();
            pick(m0,tor,m1);
            nei(i,j)=m1.mean();
        }
    nei=0.1+0.9*nei.array();
    
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

void basecall::minusconstant(int cb,MatrixXf &a,int minusprevious,int minusback)
{
    mylong cycles=(mylong)a.rows();
    if(minusprevious)
    {
        MatrixXf coef=phasing[cb].col(bufferlength[cb]);
        a-=coef*previous.row(cb);
        a.colwise()-=phasing[cb].col(bufferlength[cb]+1)*colorscale;
    }
    if(minusback)
    {
        a.resize(size[dcluster]*cycles,size[dchn]);
        a.rowwise()-=color[cb].col(size[dchn]).transpose();
        a.resize(cycles,size[dcluster]*size[dchn]);
    }
}

float basecall::mywidth(float n)
{
    return exp(-n*wl);
}




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
            for(int j=0;j<block;j++)
                fix(i,j);
        
        correctneighbour();
        data.block(cnstart*size[dcyc],0,size[dcluster]*size[dcyc],size[dchn])=zchn;
    }

    size[dcluster]=cntotal;
    
    
    return 0;
}


void basecall::init()
{
    
    cnstart=floor((float)cntotal/nblock*cnblock+0.5);
    size[dcluster]=floor(((float)cntotal/nblock*(1+cnblock))+0.5)-cnstart;
    zchn=zig.block(size[dcyc]*cnstart,0,size[dcyc]*size[dcluster],size[dchn]);
    new (& zcyc) Map<MatrixXf>(zchn.data(),size[dcyc],size[dcluster]*size[dchn]);
    
    if(cnblock==0)
    {
        MatrixXf tmp;
        tmp=zcyc.row(0);
        tmp.resize(size[dcluster],size[dchn]);
        color_kmeans(tmp);
        
        for(int i=0;i<block;i++)
            color[i]=initcolor;
        
        stein=.05,wl=30;
    }
    else
        stein=.01,wl=40;
    
    
    previous.resize(block,size[dcluster]*size[dchn]);
    colorscale=color[0].leftCols(size[dchn]).mean();
    
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
        previous.array().row(i)/=bufferstart[i];
}



basecall::basecall(Cif &z,int iters, int blocksize, int buffer,mylong maxcluster):
    zig(z.data.data(),z.data.rows(),z.data.cols()),
    zcyc(z.data.data(),z.size[dcyc],z.size[dcluster]*z.size[dchn])
{
    copycif(z,0);
    data.resize(size[dcyc]*size[dcluster],size[dchn]);
    cntotal=size[dcluster];
    itertime=iters;
    bufferwidth=buffer;
    
    block=ceil((float)size[dcyc]/blocksize);
    
    if(maxcluster==0)
        maxcluster=size[dcluster]/(625*(buffer*2+blocksize+2)*block),
        maxcluster=size[dcluster]/(maxcluster>=1?maxcluster:1)+2;
    
    nblock=ceil((float)size[dcluster]/maxcluster);
    cerr<<"Clusters are divided into "<<nblock<<" blocks.\n";
    
    if (block<1) block=1;
    blocksize=size[dcyc]/block;
    
    blockstart=new int[block];
    blocklength=new int[block];
    bufferlength=new int[block];
    bufferstart=new int[block];
    
    largeblock=size[dcyc]%blocksize;
    
    color=new MatrixXf [block];
    phasing=new MatrixXf [block];
    
    blockstart[0]=0;
    bufferstart[0]=0;
    for(int i=0;i<block;i++)
    {
        if(i>0)
            blockstart[i]=blockstart[i-1]+blocklength[i-1],
            bufferstart[i]=mymax<int>(0,blockstart[i]-buffer);
        
        blocklength[i]= i<largeblock?blocksize+1:blocksize;
        bufferlength[i]= mymin<int>(blockstart[i]+blocklength[i]+buffer,size[dcyc])-bufferstart[i];
        
    }
    
}
basecall::~basecall()
{
    delete blockstart;
    delete blocklength;
    delete bufferlength;
    delete bufferstart;
    for(int i=0;i<block;i++)
        color[i].resize(0,0),phasing[i].resize(0,0);
}


int Cif::basecalling(Cif &z,int iters, int blocksize, int buffer,mylong maxcluster)
{
 //   chnpur=size[dchn];
    z.rotate(dchn);
    basecall ba(z,iters,blocksize,buffer,maxcluster);
    ba.pipe();
    
    copycif(z,0);
    data=ba.data;
    return 1;
}




