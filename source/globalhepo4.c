

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "paxis.h"
#include "wyang.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
	float center[3], distance[2];
	int pax;
	int len,i,rank;
	PROTEIN target;

	read_protein(&target, argv[1]);
	cal_surface(&target,WATERRAD);
	get_attype(&target);
	get_target(&target, argv[2]);
	pax=get_principle_axis(&target, center);
	MPI_Init(0,0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",11,11,center[0],center[1],center[2]);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",12,12,center[0]+normvec[pax][0]*6.0,center[1]+normvec[pax][1]*6.0,center[2]+normvec[pax][2]*6.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",13,13,center[0]-normvec[pax][0]*6.0,center[1]-normvec[pax][1]*6.0,center[2]-normvec[pax][2]*6.0);
	for(i=-79;i<=79;i++)
	{		
		;
//		printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",3000+i,3000+i,center[0]+normvec[pax][0]*i*1.43/8,center[1]+normvec[pax][1]*i*1.43/8,center[2]+normvec[pax][2]*i*1.43/8);
/*	float axv[3];
	axv[0]=-0.581713;
	axv[1]=-0.038127;
     	axv[2]=0.812500;
	cal_len(&target, axv, center, 3, distance);
	printf("This is calculated length %f\n",distance[0]-distance[1]);
*/
//
//	split_the_surface(&target, pax, center, distance);
	}
	}
	len=atoi(argv[3]);
	get_best_axis(&target,pax,center,len*LENPERAA);

//	get_wyang_score(&target, pax, center, len*LENPERAA);

	free_pro(&target);

	return EXIT_SUCCESS;
}
