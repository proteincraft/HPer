#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "mem.h"

#define LINELEN 81

void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+8*i);
}
void read_protein(PROTEIN *newpro, char *pn)
{
	FILE *pf;
	char line[LINELEN];
	char conform=' ';
	int MAXATMN=0;
	int ai=0;

	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s!\n",pn);
		exit(0);
	}
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ENDMDL",5)) break;
		if(!strncmp(line,"ATOM",4)||(!strncmp(line,"HETATM",6)&&!strncmp(line+17,"MSE",3))){
			 MAXATMN++;
		}
	}
	rewind(pf);
	snew(newpro->atx, MAXATMN);
	snew(newpro->atname, MAXATMN);
	snew(newpro->attype, MAXATMN);
	snew(newpro->is_surf, MAXATMN);
	snew(newpro->at_surf, MAXATMN);
	snew(newpro->is_target, MAXATMN);

	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)||(!strncmp(line,"HETATM",6)&&!strncmp(line+17,"MSE",3))){
			if(line[13]=='H'||line[12]=='H'){
				continue;
			}
			if(line[16]!=' '){
				if(conform==' ')
					conform=line[16];
				else if(line[16]!=conform)
					continue;
			}
			getxyz(line+30,newpro->atx[ai]);
			strncpy(newpro->atname[ai],line+13,13);
			newpro->atname[ai][13]='\0';
			ai++;
		}
	}
	newpro->atmN=ai;
	
	fclose(pf);
}
