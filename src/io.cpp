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
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;



//todo list: cycle-based read; write;

class locfile
{
public:
    ifstream in;
    char  compress;
    float x;
    float y;
    int get();
    uint8_t cb;
    mylong wblock;
    mylong ib;
    mylong jb;
    mylong bufferi,buffern;
    int buffereof;
    char *buffer;
    locfile(Paras &,string &,int);
    void getcomp(uint8_t &);
    ~locfile()
    {
        in.close();
        delete [] buffer;
    }
};

void locfile::getcomp(uint8_t &a)
{
    if(bufferi==buffern)
    {
        in.read((char *)buffer,sizeof(uint8_t)*2000000);
        bufferi=0;
        buffern=(mylong)in.gcount();
        if(buffern==0)
        {
            buffereof=1;
            a=0;
            return;
        }
    }
    a=buffer[bufferi],bufferi++;
}


int locfile::get()
{
    uint8_t a;
    if(!compress)
    {
        in>>x>>y;
        return !in.eof();
    }
    else
    {
        while(!cb)
        {
            if(buffereof) return 0;
            getcomp(cb);
            ib++;
            if(ib==wblock)
                ib=0,jb++;
        }
        getcomp(a),x=ib*25+(float)a/10;
        getcomp(a),y=jb*25+(float)a/10;
        cb--;
        return 1;
    }
}

int locread(Paras &cac,string &name, MatrixXf &ll,int wide)
{
    
    locfile loc(cac,name,wide);
    mylong n=0;
    mylong crow=200000;
    if(loc.compress==-1) return -1;
    ll.resize(crow,2);
    while(loc.get())
    {
        ll(n,0)=loc.x;
        ll(n,1)=loc.y;
        n++;
        if(n>=crow)
            crow*=2,ll.conservativeResize(crow,2);
    }
    ll.conservativeResize(n,2);
    cerr<<"Succeed\n";
    return 0;
}

locfile::locfile(Paras &cac,string &name,int wide)
{
    wblock=wide;
    ib=-1;
    jb=0;
    cb=0;
    bufferi=0;
    buffern=0;
    buffereof=0;
    buffer=new char[2000010];
    string truename=cac.inpath+fs+name+cac.loctype;
	cerr << "Reading location file " << truename << "...";
    if(cac.loctype=="_pos.txt")
        compress=0,in.open(truename.data());
    else if(cac.loctype==".clocs")
        compress=1,in.open(truename.data(),ios::binary),
        in.seekg(5,ios_base::beg);
    else
    {
        cerr<<"Location type not recognized\n";
        compress=-1;
    }
    if(!in.is_open()) compress=-1;
}


template<typename myint> void Cif::read_single_cif(MatrixXf &data,string name,mylong ncluster,mylong ncyc,mylong bcluster,mylong bcyc)
{
    ifstream in(name.data(),ios::binary);
    myint *tdata=(myint *)malloc(sizeof(myint)*(ncluster*ncyc*4));
    mylong cdata=0;
    uint32_t clusterpercyc;
    uint16_t cycperfile;
	if (ncyc == 1)
		cerr << '.';
	else
		cerr<<"Reading CIF file "<<name<<"...";
    if(!in.is_open())
    {
        cerr<<"failed\n";
        data=MatrixXf::Zero(ncluster*4,ncyc);
        return;
    }
    in.ignore(7);
    in.read((char *)&cycperfile,2);
    in.read((char *)&clusterpercyc,4);
    if((mylong)cycperfile<ncyc+bcyc)
        cerr<<name<<"Warning: required Cycles exceeds existence in file\n";
    if((mylong)clusterpercyc<bcluster+ncluster)
        cerr<<"Warning: required Clusters exceeds existence in file\n";
    in.ignore(sizeof(myint)*4*bcyc*clusterpercyc);
    for(mylong i=0;i<ncyc;i++)
        for(int k=0;k<4;k++)
        {
            in.ignore(sizeof(myint)*bcluster);
            in.read((char *)&tdata[cdata],sizeof(myint)*ncluster);
            cdata+=ncluster;
            in.ignore(sizeof(myint)*(clusterpercyc-ncluster-bcluster));
        }
    Map< Matrix<myint,Dynamic,Dynamic> > tmpdata(tdata,ncluster*4,ncyc);
    data=tmpdata.template cast<float>();
    free(tdata);
	if (ncyc>1)
		cerr<<"Succeed\n";
}


int Cif::cifinfo(string name, char &intlength, mylong &cycle, mylong &clusters)
{
    ifstream in(name.data(),ios::binary);
    char buf[10]="";
    uint16_t tcyc[1];
    uint32_t tclu[1];
    in.get(buf,4);
    if (buf[0]!='C'|| buf[1]!='I'|| buf[2]!='F')
    {
        buf[3]=0;
        cerr<<name<<":Not CIF file\n";//<<(int)buf[0]<<(int)buf[1]<<(int)buf[2];
        return -1;
    }
    in.ignore(1);
    intlength=in.get();
    in.ignore(2);
    in.read((char *)tcyc,2);
    in.read((char *)tclu,4);

    cycle=tcyc[0];
    clusters=tclu[0];
    in.close();
    return 1;
}



void Cif::read(Paras &cac,string &name, mylong ncluster, mylong ncyc, mylong bcluster, mylong bcyc)
{
    char intlength;
    mylong cycleincif;
    mylong clusterpercyc;
    string inprefix=cac.inprefix;
    string insubfix=cac.insubfix;
    if(!cac.insep)
    {
        if(cifinfo(cac.inpath+fs+inprefix+name+insubfix+".cif",intlength,cycleincif,clusterpercyc)==-1)
            return;
        if(bcyc+ncyc>cycleincif||bcluster+ncluster>clusterpercyc)
            cerr<<"Warning: cluster or Cycle index exceeds limits\n",
            ncyc=0,ncluster=0;
        if(ncluster==0) ncluster=clusterpercyc-bcluster;
        if(ncyc==0) ncyc=cycleincif-bcyc;
        switch(intlength){
            case 1:
                read_single_cif<int8_t>(data,cac.inpath+fs+inprefix+name+insubfix+".cif",ncluster,ncyc,bcluster,bcyc);
                break;
            case 2:
                read_single_cif<int16_t>(data,cac.inpath+fs+inprefix+name+insubfix+".cif",ncluster,ncyc,bcluster,bcyc);
                break;
            case 3:
                read_single_cif<int32_t>(data,cac.inpath+fs+inprefix+name+insubfix+".cif",ncluster,ncyc,bcluster,bcyc);
                break;
            default:
                cerr<<"CIF Data precision not supported\n";
        }
    
    }
    else
    {

		cerr << "Reading CIF data of tile " << inprefix + name + insubfix<<" from " << ncyc << " files:\n";
		for (int i = 0; i<ncyc; i++)
			cerr << '=';
		cerr << '\n';

		char cycnum[10];
        sprintf(cycnum,"C%d.1",(int)bcyc+1);
        if(cifinfo(cac.inpath+fs+cycnum+fs+inprefix+name+insubfix+".cif",intlength,cycleincif,clusterpercyc)==-1)
            return;
        if(bcluster+ncluster>clusterpercyc)
            cerr<<"Warning: cluster requirement exceeds limits\n",ncluster=0;
        if(ncluster==0) ncluster=clusterpercyc-bcluster;
        data.resize(4*ncluster,ncyc);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i=0;i<ncyc;i++)
        {
			MatrixXf tdata;
			char cycnum[10];
            sprintf(cycnum,"C%d.1",(int)bcyc+1+i);
            switch(intlength){
                case 1:
                    read_single_cif<int8_t>(tdata,cac.inpath+fs+cycnum+fs+inprefix+name+insubfix+".cif",ncluster,1,bcluster,0);
                    break;
                case 2:
                    read_single_cif<int16_t>(tdata,cac.inpath+fs+cycnum+fs+inprefix+name+insubfix+".cif",ncluster,1,bcluster,0);
                    break;
                case 3:
                    read_single_cif<int32_t>(tdata,cac.inpath+fs+cycnum+fs+inprefix+name+insubfix+".cif",ncluster,1,bcluster,0);
                    break;
                default:
                    cerr<<"CIF Data precision not supported\n";
            }
            data.col(i)=tdata;
        }
		cerr << "\nSucceed\n";
        
    }
    size[dcluster]=ncluster;
    size[dcyc]=ncyc;
    size[dchn]=4;
    dim=dcyc;
}

void Cif::write_single_cif(string name,MatrixXf &data,mylong ncluster, mylong ncyc,int multiscale)
{
    const uint16_t scale=4000;
    ofstream out(name.data(),ios::binary);
    const int8_t head[8]={67,73,70,1,2,1,0};
    mylong tcyc=ncyc;
    mylong tcluster=ncluster;
    MatrixXf tmpdata=data;
    out.write((char *)head,7);
    out.write((char *)&tcyc,2);
    out.write((char *)&tcluster,4);
	if (multiscale)
        tmpdata*=scale;
	tmpdata = tmpdata.array().cwiseMin(32767) + 0.5;
    Matrix<int16_t,Dynamic,Dynamic> tmpdata2=tmpdata.cast<int16_t>();
    out.write((char *)tmpdata2.data(),2*ncluster*ncyc*4);
    out.close();

}

void Cif::write(Paras &cac,string &name, mylong bcyc)
{
    rotate(dcyc);
	int multiscale = 0;
	if ((data.array().abs() > 10).count() < (float)size[dcluster]*size[dcyc]*(0.5f * 4))
		multiscale = 1;
    if(!cac.outsep)
    {
		cerr << " Writing CIF file: " << cac.outpath + fs + cac.oprefix + name + cac.osubfix + ".cif" << "...";
        write_single_cif(cac.outpath+fs+cac.oprefix+name+cac.osubfix+".cif",data,size[dcluster],size[dcyc],multiscale);
		cerr << "succeed\n";
    }
    else
    {
		cerr << "Writing CIF data of tile " << cac.oprefix + name + cac.osubfix << " into " << size[dcyc] << " files:\n";
		for (int i = 0; i < size[dcyc]; i++)
			cerr << '=';
		cerr << '\n';
		
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i=0;i<size[dcyc];i++)
        {
			MatrixXf tmp;
			char cycnum[10];
            getslice(dcyc,i,tmp);

            sprintf(cycnum,"C%d.1",(int)bcyc+i+1);
            write_single_cif(cac.outpath+fs+cycnum+fs+cac.oprefix+name+cac.osubfix+".cif",tmp,size[dcluster],1,multiscale);
			cerr << '.';
        }

		cerr << "\nSucceed.\n";
    }
}


void Cif::writefq(Paras &cac, string &name)
{
	mylong mychn = size[dchn];
	string truename = cac.outpath + fs + cac.oprefix + name + cac.osubfix + ".fastq";
	ofstream out(truename.data());
	rotate(dchn);
	cerr << "Writing " << truename << "...";
	const float loc[9] = { 0, 1, 2, 4, 6, 8, 11, 14, 17 };
	const float inter[8] = { 1, 1, 2, 2, 2, 3, 3,  3 };
	char *bl = new char[(size[dcyc] + 2)*size[dcluster]], *phr = new char[(size[dcyc] + 2)*size[dcluster]];
	
#pragma omp parallel
	{
		int8_t id,id2,*fid=new int8_t[size[dcyc]];
		int8_t *sid=new int8_t[size[dcyc]];
		mylong tj;
		float *sr = new float[size[dcyc]];
		float *d = new float[size[dchn]];
#pragma omp for
		for (mylong i = 0; i < size[dcluster]; i++)
		{
			mylong c0 = i*(size[dcyc] + 2);
			double e,e2,tsr=0,tsr2=0,sig=0;
			bl[c0 + size[dcyc]] = phr[c0+size[dcyc]] = '\n';
			bl[c0 + size[dcyc]+1] = phr[c0 + size[dcyc]+1] = 0;
			for (mylong j = 0; j < size[dcyc]; j++)
			{
				tj = j + size[dcyc] * i;

				for (int k = 0; k < mychn; k++)
					d[k] = data.coeff(tj, k);

				//Find the first and second largest value and letter
				id = 0;
				e = d[0];
				for (int8_t k = 1; k < mychn; k++)
					if (e < d[k])
						e = d[k], id = k;
				id2 = 0;
				e2 = -100000;
				for (int8_t k = 0; k < mychn; k++)
					if (k != id&&e2 < d[k])
						id2 = k, e2 = d[k];
				
				sr[j] = e-e2;
				tsr+=e-e2;
				tsr2+=(e-e2)*(e-e2);
				sig += e;
				fid[j]=id;
				sid[j]=id2;
				switch (id)
				{
				case 0: bl[c0+j] = 'A'; break;
				case 1: bl[c0+j] = 'C'; break;
				case 2: bl[c0+j] = 'G'; break;
				case 3: bl[c0+j] = 'T'; break;
				default:bl[c0+j] = 'N';
				}
			
			}
			
			tsr/=size[dcyc];
			tsr2=tsr2/size[dcyc]-tsr*tsr;
			sig /= size[dcyc];
			for(mylong j=0;j<size[dcyc];j++)
			{
				int a[12];
				float b[12];
				int cloc;
				int mm=0;
				sr[j] *= 5;
				if (sr[j]>19.5) sr[j] = 19.5;
				if(sr[j]>=8)
					cloc=floor((sr[j]+1)/3+2);
				else if (sr[j]>=2)
					cloc=floor((sr[j]/2)+1);
				else cloc=floor(sr[j]);
				if (cloc > 7)
					cloc = 7;
				float beta = 1 - (sr[j] - loc[cloc]) / inter[cloc];
			
				a[mm] = cloc+1;
				b[mm] = beta;
				mm++;
				a[mm] = cloc+2;
				b[mm] = 1-beta;
				mm++;
				float term = sig > 3 ? sig : 3;
				term = term > 30 ? 30 : term;
				if (tsr2 > 1)
					tsr2 = 1;
				if (sig > 3)
					sig = 3;
				if (sig < 0)
					sig = 0;
				a[mm] = 10;      
				b[mm] = 1/(0.02+sqrt(tsr2));
				mm++;
				a[mm] = 11;    
				b[mm] = log(0.02 + tsr2);
				mm++;
				a[mm] = 12; 
				b[mm] = 1 / (0.01 + sig);
				mm++;
				a[mm] = 13; 
				b[mm] =  1 / (float)(size[dcyc]-j);
				mm++;
				a[mm] = 14;
				b[mm] = term;
				mm++;
				if (j <= 6)
				{
					a[mm] = 15 + j;
					b[mm] = 1;
					mm++;
				}

				int8_t cc = fid[j] * 3 + sid[j] - (sid[j] > fid[j]);
				if(j>0)
				{
					cc+=fid[j - 1] * 12;
					a[mm] = 22 + cc;
					b[mm] = 1;
					mm++;
				}
			
				a[mm] = tablelength;
				b[mm] = 1;
				mm++;
				float etmp=0;
				for(int kk=0;kk<mm;kk++)
					etmp+=cac.table[a[kk]-1]*b[kk];
				etmp=etmp>12?12:etmp;
				etmp = 33 + 10 * log10(1 + exp(etmp));
		//		if (etmp < 54)
		//			etmp -= (54 - etmp)*0.15;
				phr[c0+j] =(char)( etmp>=33?etmp:33 );
			}
		
		}
		delete[]d;
		delete[]sid;
		delete[]fid;
	}
	for (mylong i = 0; i < size[dcluster]; i++)
		out << "@Cluster_" << name << '_' << i + 1 << '\n' << bl+i*(size[dcyc] + 2) << "+\n" << phr+i*(size[dcyc] + 2);

	delete[]bl;
	delete[]phr;
	
	cerr << " succeed\n";
	out.close();
}


/*
void Cif::writefq(Paras &cac, string &name)
{
	mylong mychn = size[dchn];
	Matrixidx idx, idx2;
	string truename = cac.outpath + fs + cac.oprefix + name + cac.osubfix + ".fastq";
	ofstream out(truename.data());
	rotate(dchn);
	cerr << "Writing " << truename << "...";

	char *bl = new char[(size[dcyc] + 2)*size[dcluster]], *phr = new char[(size[dcyc] + 2)*size[dcluster]];

#pragma omp parallel
	{
		int8_t id,id2;
		mylong tj;
		
		float *d = new float[mychn], e,e2,etmp;
		float *dif = new float[size[dcyc]];
#pragma omp for
		for (mylong i = 0; i < size[dcluster]; i++)
		{
			mylong c0 = i*(size[dcyc] + 2);
			bl[c0 + size[dcyc]] = phr[c0+size[dcyc]] = '\n';
			bl[c0 + size[dcyc]+1] = phr[c0 + size[dcyc]+1] = 0;
			float sum = 0;
			for (mylong j = 0; j < size[dcyc]; j++)
			{
				tj = j + size[dcyc] * i;
				id = idx.coeff(tj);
				

				for (int k = 0; k < mychn; k++)
					d[k] = data.coeff(tj, k);

				//Find the first and second largest value and letter
				id = 0;
				e = d[0];
				for (int8_t k = 1; k < mychn; k++)
					if (e < d[k])
						e = d[k], id = k;
				id2 = 0;
				e2 = -100000;
				for (int8_t k = 0; k < mychn; k++)
					if (k != id&&e2 < d[k])
						id2 = k, e2 = d[k];
				dif[j] = e-e2;
				for (int8_t k = 0; k < mychn; k++)
					if (k != id&&k != id2)
						sum += d[k]>0?d[k]:-d[k];

				switch (id)
				{
				case 0: bl[c0+j] = 'A'; break;
				case 1: bl[c0+j] = 'C'; break;
				case 2: bl[c0+j] = 'G'; break;
				case 3: bl[c0+j] = 'T'; break;
				default:bl[c0+j] = 'N';
				}
				etmp =  sqrt(dif[j] / (sum / j / (mychn - 2)));
				etmp = etmp < 3.2 ? etmp * 9 - 2 : sqrt(etmp) * 9 + 10.7;
				if (etmp >= 34)
					etmp = 34 - 3 * etmp;
				etmp = etmp < 0 ? 0 : etmp;
				
				phr[c0 + j] =  (char)(33.0f+(etmp > 45.0f ? 45.0f : etmp));
			}
			sum /= mychn * (mychn-2);
		
		}
		delete[]d;
		delete[]dif;
	}
	for (mylong i = 0; i < size[dcluster]; i++)
		out << "@Cluster_" << name << '_' << i + 1 << '\n' << bl+i*(size[dcyc] + 2) << "+\n" << phr+i*(size[dcyc] + 2);

	delete[]bl;
	delete[]phr;
	
	cerr << " succeed\n";
	out.close();
}
*/
