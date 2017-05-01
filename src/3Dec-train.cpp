/*
This file is part of 3Dec, an accurate base-calling software for sequences.

Copyright (c) 2015, Bo Wang, Academy of Mathematics and Systems Science,
Chinese Academy of Sciences, Beijing 100190, China


This Source Code Form is subject to the terms of the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International Public License.
If a copy of the licence was not distributed with this file, you can obtain
one at http://creativecommons.org/licenses/by-nc-sa/4.0/
*/

#ifdef _WIN32
#include "../include/liblinear/linear.h"
#else
#include <linear.h>
#endif
#include "classtype.h"
#include <math.h>
#include <stdio.h>

void putusage()
{
	cerr<<"\nUsage:\n";
	cerr<<"3Dec-train is used to train new logestic regression model for \n";
	cerr<<"Phred quality scores of 3Dec. Training new model requires an \n";
	cerr<<"entire tile, in which all reads should have a known reference.\n";
	cerr<<"Mismatches between short reads and reference will result in \n";
	cerr<<"underestimated quality scores.\n\n";
		
	cerr<<"Before the training, you should run 3Dec with options -q -t to\n";
	cerr<<"generate corrected intensity file \"cifname.cif\", then align \n";
	cerr<<"the generated .fastq file to the reference of the short reads \n";
	cerr<<"using a mapping software such as bowtie2 or BWA; the mapping \n";
	cerr<<"results should be stored in \"samname.sam\" in the same order.\n";
	cerr<<"Then the model file \"modelname\" used by 3Dec(with option -m) \n";
	cerr<<"can be generated by typing the command:\n";
	cerr<<"\t$ 3Dec-train cifname samname modelname\n";
}

void trainfq(Cif &in, string samname,string outname)
{
	// 1:9 10CR 10 SR2 11 SR2^2 12 SIG  13 1/ENDCYC 14-21 firstcycles
	// 22-69 idx[-1]-idx-subidx 70-133 idx[-3]idx[-2]idx[-1]                
	struct problem prob;
	struct parameter param;
	const int32_t nrfeature=70;
	const double loc[9] = { 0, 1, 2, 4, 6, 8, 11, 14, 17 };
	const double inter[8] = { 1, 1, 2, 2, 2, 3, 3,  3 };
	in.rotate(dchn);
	double mdif = 0,mdif2=0;
	mylong mychn = in.size[dchn];
	mylong mycyc = in.size[dcyc];
	mylong mycluster = in.size[dcluster];
	mylong unmatched=0;
	const mylong maxnode = 29000000;
	char *buffer = new char[mycyc*20 + 102];
	char *str1= new char [mycyc*3];
	char cigar[1000];
	int8_t *match = new int8_t[mycyc + 2];
	struct feature_node  **feature;
	float *d = new float[mychn];
	double *sr = new double[mycyc];
	int8_t *fid = new int8_t[mycyc]; 
	int8_t *sid = new int8_t[mycyc];
	const int8_t matchvalue=1;
	
	FILE *fin = fopen(samname.c_str(), "r");
	FILE *out = fopen(outname.c_str(), "w");
	if(fin==NULL||out==NULL)
	{
		cerr<<"Error: At least one of the files cannot be opened\n";
		putusage();
		return;
	}
	mylong n = mycyc*mycluster;
	n = n > maxnode ? maxnode : n;
	mylong m = 0;
	feature = new struct feature_node *[n];
#ifdef _linear_1
	int32_t *y;
	y = new int32_t[n];
#else
	double *y;
	y = new double[n];
	param.init_sol=NULL;
#endif
	mylong mtk=0;
	for (mylong i = 0; i < n; i++)
		feature[i] = new struct feature_node[12];
	for (mylong tk = 0; tk < mycluster; tk++)
	{
		
		int32_t i = 0,j=0,cj=0;
		int32_t inverse = 0;
		for (int32_t kk = 0; kk < mycyc; kk++)
			match[kk] =  0;
		char ci;

		////////////////////read sam
		buffer[0] = '@';
		while (buffer[0] == '@' || buffer[0] == '\n')
			fgets(buffer, mycyc * 20 + 100, fin);
		if (tk < ((float)mycluster/2-n/mycyc)/2)
			continue;
		sscanf(buffer, "%*s%d%s%*s%*s%s", &inverse, str1, cigar);
		inverse = (inverse & 16)>0;
		if (str1[0] == '*')
			continue;
		int32_t tt = 0, noind = 0;
		while (cigar[tt]) {
			if (cigar[tt] == 'D' || cigar[tt] == 'I')
			{
				noind = 1;
				break;
			}
			tt++;
		}
		if (noind)
			continue;

		i = 70;
		while (i < 1900 && !(buffer[i] == 'M' &&buffer[i + 1] == 'D'&&buffer[i + 2] == ':' \
			&&buffer[i + 3] == 'Z'&&buffer[i + 4] == ':'))
			i++;
		if (i > 1899)
			continue;
		i += 5;
		j = 0;

		do
		{
			ci = buffer[i];
			if (ci >= '0'&&ci <= '9')
				j = j * 10 + ci - '0';
			else if (ci == 'I')
			{
				for (int32_t k = 0; k < j; k++)
					match[cj] =-1, cj++;
				j = 0;
			}
			else if (ci == 'D')
				j = 0;
			else
			{
				for (int32_t kk = 0; kk < j; kk++)
					match[cj] =1+matchvalue, cj++;
				j = 0;
				if (ci == 'A' || ci == 'C' || ci == 'G' || ci == 'T')
					match[cj] = 1+!matchvalue, cj++,unmatched++;
				else if (ci == 'N')
					match[cj] = -1, cj++;
				else if (ci == '^')
					while ((buffer[i + 1] <= 'Z'&&buffer[i + 1] >= 'A') || (buffer[i + 1] <= 'z'&&buffer[i + 1] >= 'a'))
						i++, ci = buffer[i];
				else break;
			}
			i++;
		} while (ci != 0 && ci != ' '&& ci != '\t');
		
		/////////////////////////// end readsam
		/////////////////////////// calculate feature
		double  tsr=0, tsr2=0,sig=0;
		int8_t id, id2;
		float e, e2;

		for (mylong j = 0; j < mycyc; j++)
		{
			mylong tj = j + in.size[dcyc] * tk;
			for (int32_t kk = 0; kk < mychn; kk++)
				d[kk] = in.data.coeff(tj, kk)/4000;
			id = 0;
			e = d[0];
			for (int8_t kk = 1; kk < mychn; kk++)
				if (e < d[kk])
					e = d[kk], id = kk;
			id2 = 0;
			e2 = -100000;
			for (int8_t k = 0; k < mychn; k++)
				if (k != id&&e2 < d[k])
					id2 = k, e2 = d[k];
			
			fid[j] = id;
			sid[j] = id2;
			tsr += e - e2;
			sig += e;
			tsr2+=(e-e2)*(e-e2);
			sr[j] = e - e2;
		
		}

		tsr /= mycyc;
		sig /= mycyc;
		tsr2=tsr2/mycyc-tsr*tsr;
		
	
		for (mylong j = 0; j < mycyc; j++)
		if(match[inverse?mycyc-1-j:j]>0)
		{
			int32_t mm = 0;
			
			sr[j] *=5;
			
			if (sr[j]>19.5) sr[j] = 19.5;
			
			mtk++;
			mdif += sr[j];
			mdif2 += sr[j] * sr[j];

			int32_t cloc;
			if (match[j] == 0)
				continue;
			if (sr[j] >= 8)
				cloc = floor((sr[j] + 1) / 3 + 2);
			else if (sr[j] >= 2)
				cloc = floor((sr[j] / 2) + 1);
			else cloc = floor(sr[j]);
			if (cloc > 7)
				cloc = 7;
			float beta = 1 - (sr[j] - loc[cloc]) / inter[cloc];
			feature[m][mm].index = cloc + 1;
			feature[m][mm].value = beta;
			mm++;
			feature[m][mm].index = cloc + 2;
			feature[m][mm].value = 1-beta;
			mm++;
			
			if (tsr2 > 1)
				tsr2 = 1;
			if (sig > 3)
				sig = 3;
			if (sig < 0)
				sig = 0;
			float term = sig > 3 ? sig : 3;
			term = term > 30 ? 30 : term;
			feature[m][mm].index = 10;
			feature[m][mm].value = 1 / (0.02 + sqrt(tsr2));
			mm++;
			feature[m][mm].index = 11;
			feature[m][mm].value = log(0.02+tsr2);
			mm++;
			feature[m][mm].index = 12;
			feature[m][mm].value =1/ (0.01+ sig);
			mm++;
			feature[m][mm].index = 13;
			feature[m][mm].value = 1 / (double)(mycyc - j);
			mm++;
			feature[m][mm].index = 14;
			feature[m][mm].value = term;
			mm++;
			if (j <= 6)
			{
				feature[m][mm].index = 15+j;
				feature[m][mm].value = 1;
				mm++;
			}
			
			int8_t cc = fid[j] * 3 + sid[j] - (sid[j] > fid[j]);
			if(j>0)
			{
				cc+=fid[j - 1] * 12;
				feature[m][mm].index = 22 + cc;
				feature[m][mm].value = 1;
				mm++;
			}
		
		
			feature[m][mm].index = nrfeature;
			feature[m][mm].value = 1;
			mm++;

			feature[m][mm].index = -1;
			y[m] = match[inverse?mycyc-1-j:j];

			m++;
			if (m == n) break;
		}
		if (m == n) break;
		
	}
	cerr << "Average mean of F1:" << mdif / mtk << endl;
	cerr << "Average  SD  of F1:" << sqrt(mdif2 / mtk - (mdif / mtk) *(mdif / mtk)) << endl;
	prob.y=y;
	prob.x = feature;
	prob.bias = -1;
	prob.l = m;
	prob.n = nrfeature;
	param.solver_type = L1R_LR;
	param.eps = 0.01;
	param.C = 0.2;
	param.nr_weight = 0;
	param.weight_label=NULL;
	param.weight=NULL;

	cerr<<"Unmatched bases : "<<unmatched<<" Now training...\n";

	struct model *model = train(&prob, &param);
	
	
	cerr<<"done.\n";
	
	fprintf(out,"#MD3DECFQ\n");
	for (int32_t i = 0; i <  nrfeature; i++)
		fprintf(out,"%6.4lg,", (model->w)[i]*(model->label[0]==2?1:-1));
	
	fclose(fin);
	fclose(out);
	delete[]prob.y;
	for (mylong i = 0; i < n; i++)
		delete[](feature[i]);
	delete[] feature;
	delete[] fid;
	delete[] sid;
	delete[] sr;
	delete[] d;
	delete[] str1;
	delete model;
	delete[]match;
	delete[]buffer;
	return;
}

int main(int args, char * argv[])
{
	Cif in;
	Paras par;
	string str1, str2;
	if(args!=4)
	{
		cerr<<"Error: The number of input arguments should be 3\n";
		putusage();
		return -1;
	}
	str1 = argv[1];
	in.read(par,str1);
	str1 = argv[2];
	str1+= ".sam";
	str2 = argv[3];
	trainfq(in, str1,str2);
}
