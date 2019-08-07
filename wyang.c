#include "wyang.h"
#include "score.h"
#include <math.h>
#include <stdio.h>
#include "geometry.h"
#include "paxis.h"

void cal_len(PROTEIN *pro, float axv[3], float center[3], float len, float distance[2],float dy_center1[3],float dy_center2[3])
{
	int i;
	float sign;
	float length;
	int ia;
	float dis,dis0,ppdis;
	float d_max[4], d_min[4], d[4];
	float iv[3];
//	float axv[3];
//	float test;	
	d_max[0] = -100000;
	d_min[0] = 1000000;
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->is_target[ia] == 1){
			for(i=0;i<3;i++) iv[i] = pro->atx[ia][i] - center[i];
//			for(i=0;i<3;i++) axv[i]=normvec[pax][i];
//			test = (float)cos_angle(axv,iv);
//			printf("First cos:%f\n",test);
			
			sign = fabs(cos_angle(axv,iv))/cos_angle(axv,iv);
			dis=dis_point_line(pro->atx[ia],center,axv);
			ppdis=distance2(pro->atx[ia],center);
			d[0] = sign*(sqrt(ppdis - dis));
			d[1] = d[0]*axv[0]+center[0];
			d[2] = d[0]*axv[1]+center[1];
			d[3] = d[0]*axv[2]+center[2];
			for(i=0;i<4;i++){ 
			if(fabs(d[i]) >= 1000) {
				printf("ERROR:This is d[%d]: %e \n",i,d[i]);
				printf("ERROR:This is axv:%f %f %f and center:%f %f %f\n",axv[0],axv[1],axv[2],center[0],center[1],center[2]);
			}
			}
			if(d[0] <= d_min[0]) {
					for(i=0;i<4;i++) d_min[i] = *(d+i);
					
			}
			if(d[0] >= d_max[0]) {
					for(i=0;i<4;i++) d_max[i] = *(d+i);
			}
//			printf("This is d_min now: %f This is d_max now: %f\n",d_min, d_max);
		}
	}
	distance[0] = d_max[0];
	distance[1] = d_min[0];
	for (i=0;i<3;i++) {
//		if (fabs(d_max[i+1]) < 0.001|| fabs(d_min[i+1]) < 0.001){
//			d_max[i+1] = 0	;
//			d_min[i+1] = 0;
//		}
//		if(fabs(d_max[i+1]) >= 1000 || fabs(d_min[i+1]) >= 1000) {
//			printf("ERROR: This is d_max[%d]: %e and d_min[%d]: %e\n",i+1,d_max[i+1],i+1,d_min[i+1]);
//		}
		dy_center1[i] = (d_max[i+1] + d_min[i+1])/2;
	}
//	for(i=0;i<3;i++){
		//printf("This is dy_center1: %e %e %e\n",dy_center1[0],dy_center1[1],dy_center1[2]);
//		if(fabs(dy_center1[i]) >= 1000) {
//			printf("ERROR:This is dy_center1[%d]: %e \n",i,dy_center1[i]);
//		}
//	printf("This is step %f\n",floor((distance[0]-distance[1]-4*3.6*LENPERAA)/(3.6*LENPERAA)));
}


float cal_angle(float axv[3], float paxv[3])
{
	float result;
//	float i,j;

	//i=4;
	//j=sqrt(i);
//	printf("This is 2 and 2: %f %f\n",sqrt(i),j);
//	printf("This is axv %f %f %f and paxv %f %f %f.\n",axv[0],axv[1],axv[2],paxv[0],paxv[1],paxv[2]);
//	printf("This is cos of:%f\n",cos_angle(axv,paxv));
//	result = cos_angle(axv,paxv);
//	printf("This is cos :%f\n",result);
//	printf("");

//	printf("");
	result = sqrt(fabs(cos_angle(axv,paxv)));
//	printf("This is axv %f %f %f and paxv %f %f %f and cos %f.\n",axv[0],axv[1],axv[2],paxv[0],paxv[1],paxv[2], cos_angle(axv,paxv));
	//printf("This is cos, fabs, sqrt: %f %f %f\n",cos_angle(axv,paxv),fabs(cos_angle(axv,paxv)),result);
//	result = sqrt(fabs(cos_angle(axv,paxv)));
	return result;
}




void split_the_surface(PROTEIN *pro, int pax, float center[3], float distance[2])
{
	int i, step, nstep, ia, pv0, pv1, pv2,target_2, target_3;
	float start_point[3], end_point[3], axv[3], boundary[3], temp_v[3], dis_sum[3], center_1[3], center_2[3];
	
	step = floor((distance[0]-distance[1]-4*3.6*LENPERAA)/(3.6*LENPERAA));
	for(i=0;i<3;i++) axv[i]=normvec[pax][i];
	for(i=0;i<3;i++){
		start_point[i] = center[i] + axv[i]*distance[0];
		end_point[i] = center[i] + axv[i]*distance[1];
	}
	
	if( step >= 0){
		for(nstep=0;nstep <= step; nstep++){
			target_2 = 0;
			target_3 = 0;
			for(i=0;i<3;i++) boundary[i] = start_point[i] - axv[i]*LENPERAA*3.6;
			for(ia=0;ia<pro->atmN;ia++){
				if(pro->is_target[ia] == 1){
					for(i=0;i<3;i++) temp_v[i] = pro->atx[ia][i] - boundary[i];
					if(cos_angle(temp_v, axv) > 0 ){
						pro->is_target[ia] == 2;
						target_2++;
					}
					else{
						pro->is_target[ia] == 3;
						target_3++;
					}
				}
			}
			printf("In %d cycle, target_2 is %d, target_3 is %d\n",nstep, target_2, target_3);
			pv0 = get_principle_axis_1(pro, center_1, dis_sum);
			pv1 = get_principle_axis_2(pro, center_1, dis_sum);
			pv2 = get_principle_axis_3(pro, center_2, dis_sum);
			printf("These are three lines and distance: dis0 %f, dis1 %f, dis2 %f\n", dis_sum[0], dis_sum[1], dis_sum[2]);
			printf("These are three angles: 0&1 %f, 0&2 %f, 1&2 %f\n", cal_vec_angle(axv, normvec[pv1]), cal_vec_angle(axv, normvec[pv2]), cal_vec_angle(normvec[pv1], normvec[pv2]));
			if((dis_sum[0]>dis_sum[1]+dis_sum[2]) && (cal_vec_angle(normvec[pv1], normvec[pv2]) <30||cal_vec_angle(normvec[pv1], normvec[pv2]) >150)){
			;
		}
			else{
			;
		}	
		}
	}
	else{
	//do nothing
	;
	}	
}

float get_patial_score(PROTEIN *pro, int pax, float center[3], float len)
{
	
}

int get_principle_axis_1(PROTEIN *pro, float center[3], float dis_sum[3])
{
	
	int ia,i;
	int an=0;
	int iv,dmin_v;
	float d,dmin;
	float hc1[3],hc2[3];
	float count1,count2;
	float pax[3];
	int paxv;

	for(i=0;i<3;i++)
		center[i]=0.0;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->is_surf[ia]==0||pro->is_target[ia]==0 ) continue;
		for(i=0;i<3;i++)
			center[i]+=pro->atx[ia][i];
		an++;
	}
	for(i=0;i<3;i++)
		center[i]/=an;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_line(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
			dis_sum[0] = dmin;
		}
	}
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",1,1,center[0],center[1],center[2]);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",2,2,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",3,3,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",4,4,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",5,5,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);
	for(i=-79;i<=79;i++)
	;	
//		printf("pass");
	//	printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",2000+i,2000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++)
		pax[i]=normvec[dmin_v][i];
	paxv=dmin_v;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_plane(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
			dis_sum[1];
		}
	}

	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",6,6,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",7,7,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",8,8,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",9,9,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);

	for(i=-79;i<=10;i++)
		;
//		printf("pass");
//		printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",1000+i,1000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++){
		hc1[i]=center[i]+normvec[dmin_v][i]*9.0;
		hc2[i]=center[i]-normvec[dmin_v][i]*9.0;
	}
	count1=count2=0;
	for(ia=0;ia<pro->atmN;ia++){
		//if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
		d=dis_point_line(pro->atx[ia],hc1,pax);
		if(d<4.0*4.0) count1++;
		d=dis_point_line(pro->atx[ia],hc2,pax);
		if(d<4.0*4.0) count2++;
	}
	if(count1<count2){
		for(i=0;i<3;i++)
			center[i]=hc1[i];
	}
	else{
		for(i=0;i<3;i++)
			center[i]=hc2[i];
	}
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",10,10,center[0],center[1],center[2]);

	return 	paxv;

}

int get_principle_axis_2(PROTEIN *pro, float center[3], float dis_sum[3])
{
	
	int ia,i;
	int an=0;
	int iv,dmin_v;
	float d,dmin;
	float hc1[3],hc2[3];
	float count1,count2;
	float pax[3];
	int paxv;

	for(i=0;i<3;i++)
		center[i]=0.0;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->is_surf[ia]==0||pro->is_target[ia]==0 || pro->is_target[ia]==1 || pro->is_target[ia]==3) continue;
		for(i=0;i<3;i++)
			center[i]+=pro->atx[ia][i];
		an++;
	}
	for(i=0;i<3;i++)
		center[i]/=an;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_line(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
		}
	}
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",1,1,center[0],center[1],center[2]);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",2,2,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",3,3,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",4,4,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",5,5,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);
	for(i=-79;i<=79;i++)
	;	
//		printf("pass");
	//	printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",2000+i,2000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++)
		pax[i]=normvec[dmin_v][i];
	paxv=dmin_v;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_plane(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
			dis_sum[1];
		}
	}

	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",6,6,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",7,7,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",8,8,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",9,9,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);

	for(i=-79;i<=10;i++)
		;
//		printf("pass");
//		printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",1000+i,1000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++){
		hc1[i]=center[i]+normvec[dmin_v][i]*9.0;
		hc2[i]=center[i]-normvec[dmin_v][i]*9.0;
	}
	count1=count2=0;
	for(ia=0;ia<pro->atmN;ia++){
		//if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
		d=dis_point_line(pro->atx[ia],hc1,pax);
		if(d<4.0*4.0) count1++;
		d=dis_point_line(pro->atx[ia],hc2,pax);
		if(d<4.0*4.0) count2++;
	}
	if(count1<count2){
		for(i=0;i<3;i++)
			center[i]=hc1[i];
	}
	else{
		for(i=0;i<3;i++)
			center[i]=hc2[i];
	}
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",10,10,center[0],center[1],center[2]);

	return 	paxv;

}

int get_principle_axis_3(PROTEIN *pro, float center[3], float dis_sum[3])
{
	
	int ia,i;
	int an=0;
	int iv,dmin_v;
	float d,dmin;
	float hc1[3],hc2[3];
	float count1,count2;
	float pax[3];
	int paxv;

	for(i=0;i<3;i++)
		center[i]=0.0;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->is_surf[ia]==0||pro->is_target[ia]==0 || pro->is_target[ia]==1 || pro->is_target[ia]==2) continue;
		for(i=0;i<3;i++)
			center[i]+=pro->atx[ia][i];
		an++;
	}
	for(i=0;i<3;i++)
		center[i]/=an;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_line(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
			dis_sum[2] = dmin;
		}
	}
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",1,1,center[0],center[1],center[2]);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",2,2,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",3,3,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",4,4,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",5,5,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);
	for(i=-79;i<=79;i++)
	;	
//		printf("pass");
	//	printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",2000+i,2000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++)
		pax[i]=normvec[dmin_v][i];
	paxv=dmin_v;

	for(iv=0;iv<NORMVECNUM;iv++){
		d=0;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
			d+=dis_point_plane(pro->atx[ia],center,normvec[iv]);
		}

		if(iv==0||d<dmin){
			dmin=d;
			dmin_v=iv;
		}
	}

	//printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",6,6,center[0]+normvec[dmin_v][0]*3.0,center[1]+normvec[dmin_v][1]*3.0,center[2]+normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",7,7,center[0]+normvec[dmin_v][0]*6.0,center[1]+normvec[dmin_v][1]*6.0,center[2]+normvec[dmin_v][2]*6.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",8,8,center[0]-normvec[dmin_v][0]*3.0,center[1]-normvec[dmin_v][1]*3.0,center[2]-normvec[dmin_v][2]*3.0);
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",9,9,center[0]-normvec[dmin_v][0]*6.0,center[1]-normvec[dmin_v][1]*6.0,center[2]-normvec[dmin_v][2]*6.0);

	for(i=-79;i<=10;i++)
		;
//		printf("pass");
//		printf("HETATM %4d  O   HOH A%4d    %8.3f%8.3f%8.3f \n",1000+i,1000+i,center[0]+normvec[dmin_v][0]*i*1.43/8,center[1]+normvec[dmin_v][1]*i*1.43/8,center[2]+normvec[dmin_v][2]*i*1.43/8);

	for(i=0;i<3;i++){
		hc1[i]=center[i]+normvec[dmin_v][i]*9.0;
		hc2[i]=center[i]-normvec[dmin_v][i]*9.0;
	}
	count1=count2=0;
	for(ia=0;ia<pro->atmN;ia++){
		//if(pro->is_surf[ia]==0||pro->is_target[ia]==0) continue;
		d=dis_point_line(pro->atx[ia],hc1,pax);
		if(d<4.0*4.0) count1++;
		d=dis_point_line(pro->atx[ia],hc2,pax);
		if(d<4.0*4.0) count2++;
	}
	if(count1<count2){
		for(i=0;i<3;i++)
			center[i]=hc1[i];
	}
	else{
		for(i=0;i<3;i++)
			center[i]=hc2[i];
	}
//	printf("HETATM  %3d  O   HOH A %3d    %8.3f%8.3f%8.3f \n",10,10,center[0],center[1],center[2]);

	return 	paxv;

}

void get_wyang_score(PROTEIN *pro,int pax,float center[3],float len)
{
	
}
