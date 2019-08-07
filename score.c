#include "score.h"
#include <math.h>
#include <stdio.h>
#include "geometry.h"
#include "wyang.h"

float get_vector_score(PROTEIN *pro,float vc[3],float axv[3], float len,float *allscore,float pax[3],float other_score[10])
{
	int ia, i;
	float dis,dis0,ppdis;
	float score=0.0;
	float all_bsa = 0;
	float score2=0.0;
	float all_interface=0.0;
	float halflen;
	int out;
	float distance_length[2];
	float dy_center1[3],dy_center2[3];
	float angle_factor;
	float apolar_con=0.0,polar_con=0.0,charged_con=0.0;
	cal_len(pro, axv, vc, 3, distance_length,dy_center1,dy_center2);
	len=distance_length[0]-distance_length[1];
	other_score[0] = len;
	angle_factor=cal_angle(axv,pax);
	halflen=(len/2)*(len/2);
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->is_target[ia]==1){
			all_interface = all_interface + 1;
		}
		if(score<-50) return score/len;
		dis=dis_point_line(pro->atx[ia],vc,axv);
		ppdis=distance2(pro->atx[ia],vc);
		if(dis>11.5*11.5||ppdis>30.0*30.0) continue;
		if(halflen<ppdis-dis){
			dis0=dis;
			dis=ppdis+halflen-len*sqrt(ppdis-dis);
			if(dis0/dis<0.75)
				out=1;
			else out=0;
		}
		else out=0;
		if(dis>11.0*11.0) continue;
		else if(dis<6.25*6.25){
			score+=-300.0;
		}
		else if(dis<3*3){
			score2+=-100000;
		}
		else if(out==0&&pro->is_surf[ia]==1){
			if(dis<8.30*8.30&&dis>7.30*7.30){
				score+=pro->at_surf[ia]*2;
				all_bsa += pro->at_surf[ia]*2;
				if(pro->is_target[ia]==1)
					score2+=1;
				if(pro->attype[ia]==ATM_HYDROPHOBIC){
					apolar_con += pro->at_surf[ia]*2;
				}
				else if(pro->attype[ia]==ATM_POLAR)
					polar_con += pro->at_surf[ia]*2;
				else if(pro->attype[ia]==ATM_CHARGED)
					charged_con += pro->at_surf[ia]*2;
			}
			else if(dis<9.30*9.30){
				score+=pro->at_surf[ia]*2;
				all_bsa += pro->at_surf[ia]*2;
				if(pro->is_target[ia]==1)
					score2+=2.5;
				if(pro->attype[ia]==ATM_HYDROPHOBIC){
					apolar_con += pro->at_surf[ia]*2;
				}
				else if(pro->attype[ia]==ATM_POLAR)
					polar_con += pro->at_surf[ia]*2;
				else if(pro->attype[ia]==ATM_CHARGED)
					charged_con += pro->at_surf[ia]*2;
				else if(pro->attype[ia]!=ATM_COON)
					;
			}
			else if(dis<10.00*10.00){
				score+=pro->at_surf[ia]*2;
				all_bsa += pro->at_surf[ia]*2;
				if(pro->is_target[ia]==1)
					score2+=1;
				if(pro->attype[ia]==ATM_HYDROPHOBIC){
					apolar_con += pro->at_surf[ia]*2;
				}
				else if(pro->attype[ia]==ATM_POLAR)
					polar_con += pro->at_surf[ia]*2;
				else if(pro->attype[ia]==ATM_CHARGED)
					charged_con += pro->at_surf[ia]*2;
			}
			else{
				all_bsa += pro->at_surf[ia]*2;
				score+=pro->at_surf[ia]/4;
			}
		}
		
	}
	//REMARK   490-20-11-10   11.485,   0.681    41.110   472.165     0.701      2    1.768     162.056       2.297     819.936    2798.719       0.000
	other_score[1] = all_bsa;
	other_score[2] = score2/all_interface;
	score2 = score2/all_interface*angle_factor;
	*allscore=score2;
	other_score[3] = dy_center1[0];
	other_score[4] = dy_center1[1];
	other_score[5] = dy_center1[2];
	other_score[6] = polar_con;
	other_score[7] = apolar_con;
	other_score[8] = charged_con;
	return score/len;
}
