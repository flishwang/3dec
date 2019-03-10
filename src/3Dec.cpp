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

#include "Eigen/Eigen"
#include <iostream>
#include <fstream>

#include <time.h>
#include <string.h>

#ifndef _WIN32
#include <dirent.h>
#include <sys/stat.h>
#else
#include <direct.h>
#include "../include/dirent.windows.h"
#endif

void version()
{
	cout << "3Dec Version 0.03 by Bo Wang (flish_wang@sina.com)\n";
#ifdef _OPENMP
	cout << "Multi-thread version.\n";
#else
	cout << "Single-thread version.\n";
#endif
	cout << "\nPermission is granted under Creative Commons Attribution-NonCommercial-\nShareAlike 4.0 International Public License. ";
	cout << "If the licence is not attached \nwith the software, you can obtain it at: \nhttp://creativecommons.org/licenses/by-nc-sa/4.0/\n\n";
	
}

using namespace std;

string cifpattern[100],loctype="";
int loclength;
int patterns=0;
int cycles[100][2]={{0}};
int pthread;

void putusage(int mode=0)
{
	if (mode)
	{
		cout << "\n3Dec is an accurate base-caller for Illumina Sequencing platforms including";
		cout << "\nGA, GAII, Hiseq and Miseq. It includes the module of correcting spatial";
		cout << "\ncrosstalk. It reads Cluster Intensity Files (CIF) and cluster location files";
		cout << "\n(.clocs or _pos.txt); and outputs CIF files of corrected spatial crosstalk, or";
		cout << "\ncalled sequences in FASTQ format.\n\n";
		version();
	}
	cout << "\nType \"3Dec --help\" to show the help.\n";
    cout<<"Usage:\n";
    cout<<"3Dec [options]* {-t -q | -r} <name|pattern> [name|pattern]* ...\n";
    cout<<"    -t    outputs spatial-crosstalk-corrected CIF files.\n";
    cout<<"    -q    outputs called sequences in Fastq format.\n";
    cout<<"    -r    outputs called sequences in Fastq format, re-estimating matrices\n";
    cout<<"          after correcting spatial crosstalk (slower but more accurate).\n";
    cout<<"   name   Specifies the tile name to be processed.\n";
    cout<<" pattern  A pattern XYZ+ specifies all tile names beginning with XYZ.\n\n";
    cout<<"Options:\n";
	if (mode)
	{
		cout << "    -l    Specifies the subfix (or the expand name) of \n";
		cout << "          input location files in the next input arguments(Default .clocs)\n";
		cout << "    -L    short for [-l _pos.txt].\n";
		cout << "    -s    input CIFs are seperated (eg. when input Illumina Runfolder):\n";
		cout << "          each cycle in a subfolder.(Default: input intensities from a\n";
		cout << "          single file.)\n";
		cout << "    -S    output CIFs are seperated. Will be ignored for [-q] or [-r].\n";
		cout << "          (Default: output intensities in a single CIF file.)\n";
		cout << "    -e    Specifies the total ends.(Default: only one end.) Data are \n";
		cout << "          processed independently for each ends.\n";
		cout << "    -c    specifies the begin and end cycle for each ends in the next\n";
		cout << "          arguments. Must be set after [-e]. Eg. [-c 1,101,102,109,110,210]\n";
		cout << "          specifies the cycles for the 3 ends of [-e 3].\n";
		cout << "    -i    specifies the input directory in the next arguments. Default:\n";
		cout << "          current folder.\n";
		cout << "    -o    specifies the output directory in the next arguments. Default:\n";
		cout << "          current folder.\n";
		cout << "    -m    specifies the .model file used for Phred-Score prediction. Details\n";
		cout << "          see the help of 3Dec-train.\n";
		cout << "    -n    does not correct ACC if -q or -t.\n";
		cout << "    -f    reduces iteration for latter blocks when estimating phasing. This\n";
		cout << "          will reduce calculation time while slightly reducing the accuracy.\n";
		cout << "    -p    specifies the processes used. Default: OPENMP default value.\n";
		cout << "    -z    output CIF file after correct color&phase crosstalk\n";
		cout << "    --inpath    the same as [-i].\n";
		cout << "    --loctype   the same as [-l].\n";
		cout << "    --outpath   the same as [-o].\n";
		cout << "    --version   prints 3Dec version.\n";
		cout << "    --inprefix  prefix for input\n";
		cout << "    --oprefix   prefix for output\n";
		cout << "    --insubfix  subfix for input\n";
		cout << "    --osubfix   subfix for output\n";
		cout << "        Arguments following the four commands specify the extra part of \n";
		cout << "        input and output CIF(fastq) files' names comparing with location files'\n";
		cout << "        names. The four arguments add prefix or subfix to the I/O files' names.\n\n";
		cout << "Examples:\n";
		cout << "  3Dec -i ./L001 -o ./output -q s_1_1101 s_1_12+\n";
		cout << "     This command reads location file s_1_1101.clocs in directory ./L001, then\n";
		cout << "     reads CIF file s_4_1101.cif in the same direcotory, and then does the base-\n";
		cout << "     calling and outputs s_4_1101.fastq in directory ./output. Then it searches\n";
		cout << "     the directory ./L001 for all files with the name pattern s_1_12*.clocs, and\n";
		cout << "     reads the cif files with the same tile name and write fastq files in ./output.\n\n";
		cout << "  3Dec -i ./L001 -o ./L001 -s -S -c 1,101 -t s+\n";
		cout << "     This command searches the directory ./L001 for all location files with the\n";
		cout << "     name in pattern s*.clocs, then for each location file sA.clocs, it reads\n";
		cout << "     seperated CIF files ./L001/C1.1/sA.cif, ./L001/C2.1/sA.cif, ... , \n";
		cout << "     ./L001/C101.1/sA.cif, and correct spatial crosstalk for them and then \n";
		cout << "     writes the corrected CIF files back (overwrite the original files).\n\n";
	}
	else
	{
		cout << "    Followed with its args: -l -e -c -i -o -m -p --inpath --loctype\n";
		cout << "                            --outpath --(in/o)prefix --insubfix\n";
		cout << "    Other options:          -L -s -S -n -f -z --version\n";
		cout << "Details of the options see the help documencation.\n";
	}
}




char * rmquote(char *str)
{
    int i=0;
    if(str[0]!='\"')
        return str;
    while(str[i+1]&&str[i+1]!='\"')
        str[i]=str[i+1],i++;
    str[i]=0;
    return str;
}

int patternfilter(const struct dirent *a)
{
    int n=(int)strlen(a->d_name);
    if(n<loclength)
        return 0;
    for(int i=0;i<loclength;i++)
        if(a->d_name[n-loclength+i]!=loctype[i])
            return 0;
    
    for(int i=0;i<patterns;i++)
    {
        int j=0,same=1,mylen=(int)strlen(a->d_name);
        while(cifpattern[i][j]||a->d_name[j])
        {
            if(cifpattern[i][j]=='+')
                return 1;
			if (cifpattern[i][j] == 0 && j + loclength == mylen)
				return 1;
            if(cifpattern[i][j]!=a->d_name[j])
            {
                same=0;
                break;
            }
            j++;
        }
        if(same)
            return 1;
    }
    return 0;
}
#ifdef _WIN32
int alphasort=0;


int scandir(const char* dirname,struct dirent ***namelist, int (*pattern)(const struct dirent *),int nnn)
{
    int n=0;
    DIR *mydir=opendir(dirname);
    if(!mydir) return 0;
    dirent *tmp;
    *namelist=(struct dirent **)malloc(sizeof(*namelist)*6000);
    do
    {
        tmp=readdir(mydir);
        if(!tmp) break;
        if(pattern(tmp))
		{
			(*namelist)[n]=(struct dirent *)malloc(sizeof(struct dirent));
			strcpy((*namelist)[n]->d_name,tmp->d_name);
            n++;
		}
    }
    while(tmp);
    return n;
}


#endif

int checkarg(int &i,int &j,char cc,int args, char *argv[])
{
    if(argv[i][j+1])
    {
        cerr<<"Error:Unexpected control character after ["<< cc <<"]:\\"<<(int)argv[i][j+1]<<"\n";

        return 0;
    }
    i++,j=0;
    if(i>args||argv[i][0]=='-')
    {
        cerr<<"Error:Control character ["<<cc<<"] needs followed arguments\n";
        return 0;
    }
    return 1;
}
int readmodel(Paras & par, char * name)
{
	FILE *in = fopen(name, "r");
	char test[20];
	fscanf(in, "%s", test);
	if (strcmp(test, "#MD3DECFQ"))
	{
		cerr << "Error: Cannot load the model file.\n";
		return -1;
	}
	for (int i = 0; i < tablelength; i++)
		fscanf(in, "%f", &par.table[i]), fgetc(in);
	return 0;
}

int analyseargs(int args, char * argv[],Paras &par)
{
    int i=1;
    int pp=0;
    while(i<args)
    {
        if(argv[i][0]=='-')
        {
            if(argv[i][1]=='-')
            {
               
                if(!strcmp(argv[i],"--loctype"))
                    i++,par.loctype=rmquote(argv[i]);
                else if(!strcmp(argv[i],"--inprefix"))
                    i++,par.inprefix=rmquote(argv[i]);
                else if(!strcmp(argv[i],"--insubfix"))
                    i++,par.insubfix=rmquote(argv[i]);
		        else if(!strcmp(argv[i],"--oprebfix"))
                    i++,par.oprefix=rmquote(argv[i]);
                else if(!strcmp(argv[i],"--osubfix"))
                    i++,par.osubfix=rmquote(argv[i]);
                else if(!strcmp(argv[i],"--inpath"))
                    i++,par.inpath=rmquote(argv[i]);
				else if (!strcmp(argv[i], "--outpath"))
					i++, par.insubfix = rmquote(argv[i]);
				else if (!strcmp(argv[i], "--version"))
				{
					version();
					exit(0);
				}
				else if (!strcmp(argv[i], "--help"))
				{
					putusage(1);
					exit(0);
				}
				else 
				{
					cerr<<"Error:Unrecognized command "<<argv[i]<<endl;
					return 0;
				}
            }
            else
            {
                int j=1,arging=1;
                while(argv[i][j]&&arging)
                {
                    switch(argv[i][j])
                    {
                        case 'L':
                            par.loctype="_pos.txt";
                            break;
                        case 'l':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            par.loctype=rmquote(argv[i]),arging=0;
                            break;
                        case 's':
                            par.insep=1;
                            break;
                        case 'S':
                            par.outsep=1;
                            break;
						case 'f':
							par.iters = -3;
							break;
                        case 'e':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            if(pp)
                            {
                                cerr<<"Error:Argument [e] must be set before [c]\n";
                                return 0;
                            }
                            sscanf(argv[i],"%d",&par.totalend),arging=0;
                            break;
                        case 'c':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            pp=1;
                            for(int k=0;k<par.totalend;k++)
                                sscanf(argv[i],"%d,%d,",&cycles[k][0],&cycles[k][1]),
                                cycles[k][1]=cycles[k][1]-cycles[k][0]+1,
                                cycles[k][0]--;
                            arging=0;
                            break;
                        case 'i':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            par.inpath=rmquote(argv[i]),arging=0;
                            break;
						case 'm':
							if (!checkarg(i, j, argv[i][j], args, argv))
								return 0;
							if (readmodel(par, rmquote(argv[i])) < 0)
								return 0;
							arging = 0;
							break;
                        case 'o':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            par.outpath=rmquote(argv[i]),arging=0;
                            break;
                        case 'p':
                            if(!checkarg(i,j,argv[i][j],args,argv))
                                return 0;
                            sscanf(argv[i],"%d",&pthread);
                            #ifdef _OPENMP
                            omp_set_num_threads(pthread);
                            #endif
                            arging=0;
                            break;
                        case 't':
                            par.ocif=1;
                            break;
                        case 'z':
                            par.oinnercif=1;
                            break;
                        case 'q':
                            par.ofastq=1;
                            break;
                        case 'r':
                            par.ofastq=1;
                            par.oiter=1;
                            break;
						case 'n':
							par.donothing = 1; 
							break;
                        default:
                            cerr<<"Error:Unrecognized control character ["<<argv[i][j]<<"]\n";
                            return 0;
                    }
                    j++;
                }
            
            }
            
        }
        else
            cifpattern[patterns]=argv[i],patterns++;
        i++;
    }
    if(!par.ofastq&&!par.oiter&&!par.ocif)
    {
        cerr<<"Error:At least one of [t], [q] or [r] should be chosen\n";
        return 0;
    }
    
    loclength=(int)par.loctype.length();
    loctype=par.loctype;
    return 1;
    
}


int main(int args, char *argv[])
{
    time_t start,end;

    start=time(NULL);
    string cnsubfix;
    Cif z1,z2;
    MatrixXf ll;
    Paras par;
    string str;
    int n;
    
    struct dirent **namelist;

    
    if(!analyseargs(args,argv,par))
    {
        putusage();
        return -1;
    }
	version();
    cnsubfix=par.osubfix;
    n=scandir(par.inpath.data(),&namelist,patternfilter,alphasort);
	cerr <<"The following "<< n << " tiles are selected:";
	for (int i = 0; i < n; i++)
	{
		if (i % 5 == 0)
			cerr << '\n';
		string name = namelist[i]->d_name;
		name = name.substr(0, name.length() - loclength);
		cerr << name << "  ";
		
	}
	cerr<<'\n';

#ifndef _WIN32
    mkdir((par.outpath+fs).c_str(),S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
#else
    _mkdir((par.outpath+fs).c_str());
#endif
    

    for(int i=0;i<n;i++)
    {
        string name=namelist[i]->d_name;
        name=name.substr(0,name.length()-loclength);
		
        for(int j=0;j<par.totalend;j++)
        {
            if(!cycles[j][1]&&par.insep)
            {
                cerr<<"Error: cycles must be specified for each ends for RunFolder inputs\n";
                return -1;
            }
            z1.read(par,name,0,cycles[j][1],0,cycles[j][0]);
            if(!par.oiter&&!par.ocif)
            {
                locread(par,name,ll);
                z1.basecalling(z1,par.iters);
                if (par.oinnercif){
                    	par.osubfix=cnsubfix+(string)"_inner";
                    	z1.write(par,name,cycles[j][0]);
                    }
				if (!par.donothing)
					correct_acc(z1,z1,ll);
                if(par.totalend>1)
                {
                    char num[10];
                    sprintf(num,"%d",j);
                    par.osubfix=cnsubfix+(string)"_end"+num;
                }
		
                z1.writefq(par,name);
		
                par.osubfix=cnsubfix;
                
            }
            else
            {
				if (par.ocif&&par.ofastq&&!par.oiter)
				{
					locread(par, name, ll);
					z1.basecalling(z1, par.iters);
					if (par.oinnercif){
                    	par.osubfix=cnsubfix+(string)"_inner";
                    	z1.write(par,name,cycles[j][0]);
                    }
					correct_acc(z1, z1, ll);
					if (par.totalend>1)
					{
						char num[10];
						sprintf(num, "%d", j);
						par.osubfix = cnsubfix + (string)"_end" + num;
					}
					cerr << "Writing constant Phred-score for alignment!\n";
					for (int i = 0; i < tablelength-1; i++)
						par.table[i] = 0;
					par.table[tablelength-1] = 5;
					z1.writefq(par, name);
					par.osubfix = cnsubfix;
				}
				else if	(!par.donothing)
				{
					locread(par, name, ll);
					z2.basecalling(z1,par.iters);
					correct_acc(z1, z2, ll);
				}
				if (par.ocif)
				{
					if (i == 0 && par.outsep)
					{
						for (int i = 0; i<z1.size[dcyc]; i++)
						{
							char num[10];
							sprintf(num, "C%d.1", cycles[j][0] + i + 1);
#ifndef _WIN32
							mkdir((par.outpath + fs + num + fs).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#else
							_mkdir((par.outpath + fs + num + fs).c_str());
#endif
						}
					}

                    if(!par.outsep&&par.totalend>1)
                    {
                        char num[10];
                        sprintf(num,"%d",j);
                        par.osubfix=cnsubfix+(string)"_end"+num;
                    }
                    z1.write(par,name,cycles[j][0]);
                    par.osubfix=cnsubfix;
                }
                z2.data.resize(0,0);
                if(par.oiter)
                {
                    z1.basecalling(z1,par.iters);
                    if (par.oinnercif){
                    	par.osubfix=cnsubfix+(string)"_inner";
                    	z1.write(par,name,cycles[j][0]);
                    }
                    if(par.totalend>1)
                    {
                        char num[10];
                        sprintf(num,"%d",j);
                        par.osubfix=cnsubfix+(string)"_end"+num;
                    }
                    z1.writefq(par,name);
                    par.osubfix=cnsubfix;
                }
                
            }
            
        }
        
        free(namelist[i]);
    }
    free(namelist);
    end=time(NULL);
    cerr<<"Total time:"<<difftime(end,start)<<endl;
    return 0;
    
}
