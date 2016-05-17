#include "classtype.h"
#include "cacc_funs.h"

#include "Eigen/Eigen"
#include <iostream>
#include <fstream>

#include <time.h>

using namespace std;



int main(int args, char *argv[])
{
    time_t start,end;
    start=time(NULL);
    int cb,nb;
    
    Cif z1,z2;
    Paras par;
    string str;
    sscanf(argv[3],"%d",&cb);
    sscanf(argv[4],"%d",&nb);    

  
    MatrixXf ll;
    str=argv[1];
    z1.read(par,str);
    z1.basecalling(z1,4,cb,6,nb);
    str=argv[2];
    z1.writefq(par,str);
    
    end=time(NULL);
    cerr<<"Total time:"<<difftime(end,start)<<endl;
    return 0;
    
}
