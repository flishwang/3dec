/*
This file is part of 3Dec, an accurate base-calling software for sequences.

Copyright (c) 2015, Bo Wang, Academy of Mathematics and Systems Science,
Chinese Academy of Sciences, Beijing 100190, China


This Source Code Form is subject to the terms of the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International Public License.
If a copy of the licence was not distributed with this file, you can obtain
one at http://creativecommons.org/licenses/by-nc-sa/4.0/
*/

#ifndef _classtype
#define _classtype

#include <iostream>
#include <string>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#ifdef _OPENMP
#include <omp.h>

#else

#define omp_get_thread_num() 0
#endif

#include <inttypes.h>

#ifdef _WIN32
    const char fs='\\';
#else
    const char fs='/';
#endif

using namespace Eigen;
using namespace std;

const int tablelength = 70;
typedef int dym;
typedef int32_t mylong;
typedef Matrix<int8_t,Dynamic,Dynamic> Matrixidx;

class Paras
{
public:
    Paras()
    {
        insep=0;
        outsep=0;
        cif=1;
        currentend=1;
        totalend=1;
        loctype=".clocs";
        inprefix="";
        insubfix="";
        inpath=".";
        outpath=".";
        oprefix="";
        osubfix="";
		iters = 4;
        ocif=0;
        ofastq=0;
        oiter=0;
	donothing = 0;
        table << -5.26,-1.883,-0.2153, 1.726, 2.253, 2.066, 1.913,     0,     0,0.2074,-0.8086,0.1871,-1.242,0.5032,0.2255,-0.2874,     0,0.1121, 0.364,0.3755, 0.375,0.9716,0.4276,0.5448, 1.405, 1.454, 1.326, 1.702, 1.853, 2.039, 1.467, 1.018, 0.732, 1.601,0.8862, 1.239,0.2015,0.9255,0.2033, 1.204, 1.686,0.8465, 1.178,0.4473,0.1083,0.7897, 1.183, 0.129,0.8025,0.8445,0.2447, 0.986, 1.383,0.3707, 2.002, 1.489, 1.517,0.7776,0.6369, 1.221, 1.371,0.3425,0.6357, 1.021, 1.412, 1.489,0.7531,0.8672,-0.116,0.05525;
    }
    bool insep,outsep;
    bool cif;
    int currentend;
    int totalend;
	int iters;
    string loctype;
    string inpath,outpath;
    string inprefix;
    string insubfix;
    string oprefix;
    string osubfix;
    Matrix<float,tablelength,1> table;
    bool ocif,ofastq,oiter,donothing;
};

const dym dcluster=0, dchn=1, dcyc=2;

class Cif
{
private:
    template<typename myint> void read_single_cif(MatrixXf &,string ,mylong,mylong,mylong,mylong);
    void write_single_cif(string ,MatrixXf &,mylong,mylong,int);
    int cifinfo(string ,char &,mylong &,mylong &);
public:
    // Three storage of data:       ||     eigen expression
    //       dim=dcyc: [cyc, chn, cluster];(Default cif when read) || cluster chn cyc
    //       dim=dchn: [chn, cluster, cyc];                        || cyc cluster chn
    //       dim=dcluster: [cluster,chn,cyc];                      || cyc chn cluster
    // Expressions in eigen matrixes are reversed.
    dym dim;
    dym topdim();
    MatrixXf data;
    mylong size[3];
    // rotate:  rotate so that the para dim is at last place (eigen)
    void rotate(dym);
    void rotate_top(dym);
    // I/O operations  read/write data of para1 clusters, para2 cycles
    //                 from para3th clusters and para4th cycles.
    void read(Paras &, string &, mylong ncluster=0, mylong ncyc=0, mylong bcluster=0, mylong bcyc=0);
    void write(Paras &, string &,mylong bcyc=0);
    // multi    multiply matrix at the given dimension
    void multi(dym,MatrixXf &);
    // get the para slice in the dim
    void getslice(dym,mylong, MatrixXf &);
    void assignslice(dym,mylong, MatrixXf &);

    //
    int copycif(Cif &,int copydata=1);
    int basecalling(Cif &,int iters=4,int blocksize=14,int buffer=6,mylong maxcluster=0);
    
    void writefq(Paras &, string &);
    
};

int locread(Paras &,string &, MatrixXf &,int wide=82);


#endif
