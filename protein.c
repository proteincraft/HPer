#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "surface.h"
#include "mem.h"

float cal_surface(PROTEIN *pro,float water_rad)
{
	int ia,i;
	float surf;
	float *sol_rad;
	float *at_surf;

	snew(sol_rad, pro->atmN);
	snew(at_surf, pro->atmN);
	
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->atname[ia][0]=='C'){
			if(pro->atname[ia][1]=='A'||pro->atname[ia][1]=='B')
				sol_rad[ia]=1.87;
			else if(!strncmp(pro->atname[ia]+4,"HIS",3)||!strncmp(pro->atname[ia]+4,"PHE",3)||!strncmp(pro->atname[ia]+4,"TYR",3)||!strncmp(pro->atname[ia]+4,"TRP",3))
				sol_rad[ia]=1.76;
			else if(pro->atname[ia][1]==' ')
				sol_rad[ia]=1.76;
			else if(!strncmp(pro->atname[ia]+4,"ASP",3)||!strncmp(pro->atname[ia]+4,"ASN",3)||(pro->atname[ia][1]=='D'&&(!strncmp(pro->atname[ia]+4,"TYR",3)||!strncmp(pro->atname[ia]+4,"TRP",3))))
				sol_rad[ia]=1.76;
			else
				sol_rad[ia]=1.87;
		}
		else if(pro->atname[ia][0]=='N'){
			if(!strncmp(pro->atname[ia]+4,"LYS",3)&&pro->atname[ia][1]=='Z')
				sol_rad[ia]=1.50;
			else
				sol_rad[ia]=1.65;
		}
		else if(pro->atname[ia][0]=='O'){
			sol_rad[ia]=1.40;
		}
		else if(pro->atname[ia][0]=='S'){
			sol_rad[ia]=1.85;
		}
		else
			sol_rad[ia]=1.85;
	}

	surf=surfarea(pro->atmN, pro->atx,sol_rad,at_surf,SURF_CAL_DENSIT, water_rad);
	for(ia=0;ia<pro->atmN;ia++){
		pro->at_surf[ia]=at_surf[ia];
		if(at_surf[ia]>0.01)
			pro->is_surf[ia]=1;
		else
			pro->is_surf[ia]=0;
	}

	sfree(sol_rad);
	sfree(at_surf);

	return surf;
}
void get_attype(PROTEIN *pro)
{
	int ia;

	for(ia=0;ia<pro->atmN;ia++){
		if(strncmp(pro->atname[ia],"N  ",3)&&strncmp(pro->atname[ia],"CA ",3)&&strncmp(pro->atname[ia],"C  ",3)&&strncmp(pro->atname[ia],"O  ",3)&&strncmp(pro->atname[ia],"OXT",3)){
			//if(!strncmp(pro->atname[ia]+4,"TYR",3)&&!strncmp(pro->atname[ia],"OH ",3))
			//	printf("stop\n");
			if(!strncmp(pro->atname[ia]+4,"MSE",3)){
				pro->attype[ia]=ATM_HYDROPHOBIC;
				continue;
			}
		}
		if(!strncmp(pro->atname[ia]+4,"ASP",3)){
			if(!strncmp(pro->atname[ia],"OD1",3)||!strncmp(pro->atname[ia],"OD2",3))
				pro->attype[ia]= ATM_CHARGED;
			else if (!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]=ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"GLU",3)){
			if (!strncmp(pro->atname[ia],"OE1",3)||!strncmp(pro->atname[ia],"OE2",3))
				pro->attype[ia]= ATM_CHARGED;
			else if(!strncmp(pro->atname[ia],"CB ",3))
				;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"LYS",3)){
			if(!strncmp(pro->atname[ia],"NZ ",3))
				pro->attype[ia]= ATM_CHARGED;
			else if (!strncmp(pro->atname[ia],"CE ",3)||!strncmp(pro->atname[ia],"CD ",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"ARG",3)){
			if(!strncmp(pro->atname[ia],"NH1",3)||!strncmp(pro->atname[ia],"NH2",3)||!strncmp(pro->atname[ia],"NE ",3))
				pro->attype[ia]= ATM_CHARGED;
			else if (!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CD ",3)||!strncmp(pro->atname[ia],"CZ ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"HIS",3)){
			if(!strncmp(pro->atname[ia],"ND1",3)||!strncmp(pro->atname[ia],"NE2",3))
				pro->attype[ia]=ATM_POLAR;
			else if(!strncmp(pro->atname[ia],"CD2",3)||!strncmp(pro->atname[ia],"CE1",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]=ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"ASN",3)){
			if(!strncmp(pro->atname[ia],"OD1",3)||!strncmp(pro->atname[ia],"ND2",3))
				pro->attype[ia]= ATM_POLAR;
			else if(!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia] = ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"GLN",3)){
			if(!strncmp(pro->atname[ia],"OE1",3)||!strncmp(pro->atname[ia],"NE2",3))
				pro->attype[ia]= ATM_POLAR;
			else if (!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia] = ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"CYS",3)){
			if(!strncmp(pro->atname[ia],"SG ",3))
				pro->attype[ia]= ATM_POLAR;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"SER",3)){
			if(!strncmp(pro->atname[ia],"OG ",3))
				pro->attype[ia]= ATM_POLAR;
			else if(!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia] = ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"THR",3)){
			if(!strncmp(pro->atname[ia],"OG1",3))
				pro->attype[ia]= ATM_POLAR;
			else if (!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG2",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"TYR",3)){
			if(!strncmp(pro->atname[ia],"OH ",3))
				pro->attype[ia]= ATM_POLAR;
			else if(!strncmp(pro->atname[ia],"CE2",3)||!strncmp(pro->atname[ia],"CD2",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CD1",3)||!strncmp(pro->atname[ia],"CE1",3)||!strncmp(pro->atname[ia],"CZ ",3)||!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		
		}
		else if(!strncmp(pro->atname[ia]+4,"TRP",3)){
			if(!strncmp(pro->atname[ia],"NE1",3))
				pro->attype[ia]= ATM_POLAR;
			else if(!strncmp(pro->atname[ia],"CD1",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CD2",3)||!strncmp(pro->atname[ia],"CE3",3)||!strncmp(pro->atname[ia],"CZ3",3)||!strncmp(pro->atname[ia],"CH2",3)||!strncmp(pro->atname[ia],"CZ2",3)||!strncmp(pro->atname[ia],"CE2",3)||!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"ALA",3)){
			if(!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"PHE",3)){
			if(!strncmp(pro->atname[ia],"CD1",3)||!strncmp(pro->atname[ia],"CE1",3)||!strncmp(pro->atname[ia],"CZ ",3)||!strncmp(pro->atname[ia],"CE2",3)||!strncmp(pro->atname[ia],"CD2",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CB ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"ILE",3)){
			if(!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG1",3)||!strncmp(pro->atname[ia],"CG2",3)||!strncmp(pro->atname[ia],"CD1",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"LEU",3)){
			if(!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CD2",3)||!strncmp(pro->atname[ia],"CD1",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"VAL",3)){
			if(!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG1",3)||!strncmp(pro->atname[ia],"CG2",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"MET",3)){
			if(!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"CD ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else if(!strncmp(pro->atname[ia]+4,"PRO",3)){
			if(!strncmp(pro->atname[ia],"CB ",3)||!strncmp(pro->atname[ia],"CG ",3)||!strncmp(pro->atname[ia],"SD ",3)||!strncmp(pro->atname[ia],"CE ",3))
				pro->attype[ia]= ATM_HYDROPHOBIC;
			else
				pro->attype[ia]=ATM_OTHERS;
		}
		else
			pro->attype[ia]=ATM_OTHERS;
	}
}

void get_target(PROTEIN *pro, char *target_f)
{
	int trn=0;
	FILE *f;
	char (*Tresname)[10];
	char line[100];
	int ir,ia;
	int flag;
	snew(Tresname,pro->atmN);
	if((f=fopen(target_f,"r"))==NULL){
		printf("ERROR: Can not open the mutable residue file %s!\n",target_f);
		exit(0);
	}
	while(fgets(line, 100,f)){
		if(line[0]=='%') continue;
		if(trn>=pro->atmN) continue;
		strncpy(Tresname[trn],line,9);
		trn++;
	}
	fclose(f);
	for (ia=0;ia<pro->atmN;ia++){
	pro->is_target[ia]=0;
	}	
	for(ir=0;ir<trn;ir++){
		flag = 0;
		for (ia=0;ia<pro->atmN;ia++){
		//	printf("This is target residue:%s and atom:%s\n",Tresname[ir],pro->atname[ia]+4);
			if(!strncmp(pro->atname[ia]+4,Tresname[ir],9)){
//				printf("ATOM   2033  %s    %8.3f%8.3f%8.3f \n",pro->atname[ia],pro->atx[ia][0],pro->atx[ia][1],pro->atx[ia][2]);
				pro->is_target[ia]=1;
				flag += 1;
				//printf("This is flag:%f\n",flag);
			}
		}
		if(flag == 0){
			printf("ERROR: Can not find %s in the pdbfile!\n",Tresname[ir]);
			exit(0);
			}
		}

	sfree(Tresname);
}

void free_pro(PROTEIN *pro)
{
	sfree(pro->atx);
	sfree(pro->atname);
	sfree(pro->attype);
	sfree(pro->is_surf);
	sfree(pro->at_surf);
	sfree(pro->is_target);
}
