#include "geometry.h"
#include <math.h>

float distance2(float x1[3], float x2[3])
{
	return (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]);
}
void getdx(float x1[3],float x2[3],float dx[3])
{
	dx[0]=x1[0]-x2[0];
	dx[1]=x1[1]-x2[1];
	dx[2]=x1[2]-x2[2];
}
void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}
float iprod(float a[3],float b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

float dis_point_line(float p[3],float lp[3],float lv[3])
{
	float M0M1[3];
	float sMM[3];
	float d;

	getdx(p,lp,M0M1);
	oprod(lv,M0M1,sMM);
	d=(sMM[0]*sMM[0]+sMM[1]*sMM[1]+sMM[2]*sMM[2]);

	return d;
}
float dis_point_plane(float p[3],float lp[3],float lv[3])
{
	float M0M1[3];
	float d;

	getdx(p,lp,M0M1);
	d=iprod(M0M1,lv);

	return d*d;
}
float cos_angle(const float a[3],const float b[3])
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  float   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; 
  
  ip=ipa=ipb=0.0;
  for(m=0; m<3; m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
//printf("GEO:cos to return: %f\n",cos); 
  return cos;
}
float cal_vec_angle(float v1[3], float v2[3])
{
	float cosa;
	float a;

	cosa=cos_angle(v1,v2);
	a=acos(cosa);

	return TODEG(a);
}
