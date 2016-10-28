#include <stdio.h>
#include <iostream>
#include <inttypes.h>
using namespace std;
int main(int args,char *argv[]){
	FILE *in,*out;
	int32_t num;
	char dd;
	float *buffer;
	char begin[12]={1,0,0,0,0,0,128,63};
	if(args!=3){
		cerr<<"Error: #arguments should be 2!\n";
		cerr<<"Usage:\n    locs2txt fileA fileB\n        Read the locs file fileA, convert it into _pos.txt format, and output as fileB.\n"<<endl;
		return -1;
	}
	in=fopen(argv[1],"rb");
	for(int i=0;i<8;i++)
		if((dd=fgetc(in))!=begin[i]){
			cerr<<"Error: Cannot recognize this file!\n";
			cerr<<"Usage:\n    locs2txt fileA fileB\n        Read the locs file fileA, convert it into _pos.txt format, and output as fileB.\n"<<endl;
			return -1;
		}
	fread((void *) &num,4,1,in);
	out=fopen(argv[2],"wb");
	buffer=new float[2*num+1];
	fread((void *)buffer,4,2*num,in);
	for(int i=0;i<num;i++){
		fprintf(out,"%.2f %.2f\n",buffer[2*i],buffer[2*i+1]);
	}
	fclose(in);
	fclose(out);
	
	return 0;
}
