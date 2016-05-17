
#include "classtype.h"
#include "cacc_funs.h"

#include "Eigen/Eigen"
#include <iostream>

using namespace Eigen;
using namespace std;

//getopt_long    #include <getopt.h>




int main(int args, char *argv[])
{
    time_t start,end;
    start=time(NULL);
    
    
    Cif z1;
    Paras par;
    string str;
    
    /*
     str=argv[1];
     z1.read(par,str);
     z1.basecalling(z1,4,14,6);
     str=argv[2];
     z1.writefq(par,str);
     z1.write(par,str);
     
     */
    
    MatrixXf ll;
    str=argv[1];
    z1.read(par,str);
    z1.basecalling(z1);
    
    
    
    str=argv[2];
    z1.write(par,str);
    
    end=time(NULL);
    cerr<<"Total time:"<<difftime(end,start)<<endl;
    
    str=argv[1];
    if(locread(par,str,ll)==-1)
    {
        par.loctype="_pos.txt";
        locread(par,str,ll);
    }
    
    correct_acc(z1,z1,ll);
    par.osubfix="_new";
    z1.write(par,str);
    
    end=time(NULL);
    cerr<<"Total time:"<<difftime(end,start)<<endl;
    return 0;
    
}
