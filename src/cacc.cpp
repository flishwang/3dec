
#include "classtype.h"
#include "cacc_funs.h"


#include <iostream>

using namespace Eigen;
using namespace std;

//getopt_long    #include <getopt.h>




int main(int args, char *argv[])
{
    time_t start,end;
    start=time(NULL);
    
    
    Cif z1,z2;
    Paras par;
    string str;
    

    
    MatrixXf ll;
    str=argv[1];
    z1.read(par,str);
    z2.basecalling(z1);
    
    end=time(NULL);
    cerr<<"Basecalling time:"<<difftime(end,start)<<endl;
    
    if(locread(par,str,ll)==-1)
    {
        par.loctype="_pos.txt";
        locread(par,str,ll);
    }
    correct_acc(z1,z2,ll);
    str=argv[2];
    z1.write(par,str);

    end=time(NULL);
    cerr<<"Total time:"<<difftime(end,start)<<endl;
    return 0;
    
}
