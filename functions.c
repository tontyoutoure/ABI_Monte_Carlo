#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>
#include"REF.h"
#define pi 3.141592653589793
#define sqrt2 1.4142135623730950488016887242097
/*The next two functions are for coordinate transformation between Cartesian coordinate and spherical coordiante, 
pc is an array of Cartesian coordinate, ps is an array of spherical polar coordinate */
/*
void sph2car(double *ps, double *pc)
{
	*pc=*ps*sin(*(ps+1))*cos(*(ps+2));
	*(pc+1)=*ps*sin(*(ps+1))*sin(*(ps+2));
	*(pc+2)=*ps*cos(*(ps+1));
}

void car2sph(double *pc, double *ps)
{
	double xy=hypot(pc[0],pc[1]);
	ps[0]=hypot(xy, pc[2]);
	if(!ps[0])
	{
		ps[1]=0;
		ps[2]=0;
	}
	else if(!xy)
	{
		ps[2]=0;
		if(pc[2]>0)
			ps[1]=0;
		else
			ps[1]=pi;			
	}
	else
	{
		ps[1]=acos(pc[2]/ps[0]);
		if(pc[1]>0)
			ps[2]=acos(pc[0]/xy);
		else
			ps[2]=acos((-pc[0])/xy)+pi;
	}
}
*/
/*The next two functions are for matrix multiplication, mm3 is for 3x3 times 3x3, mm31 is fore 3x3 times 1x3 
(of course the latter matrix is acturally a array, but seen as a 1x3 matrix, after all it's for coordinate rotation)*/
/*
inline void mm3(double *p1, double *p2, double *p3)
{
	p3[0]=p1[0]*p2[0]+p1[1]*p2[3]+p1[2]*p2[6];
	p3[1]=p1[0]*p2[1]+p1[1]*p2[4]+p1[2]*p2[7];
	p3[2]=p1[0]*p2[2]+p1[1]*p2[5]+p1[2]*p2[8];
	p3[3]=p1[3]*p2[0]+p1[4]*p2[3]+p1[5]*p2[6];
	p3[4]=p1[3]*p2[1]+p1[4]*p2[4]+p1[5]*p2[7];
	p3[5]=p1[3]*p2[2]+p1[4]*p2[5]+p1[5]*p2[8];
	p3[6]=p1[6]*p2[0]+p1[7]*p2[3]+p1[8]*p2[6];
	p3[7]=p1[6]*p2[1]+p1[7]*p2[4]+p1[8]*p2[7];
	p3[8]=p1[6]*p2[2]+p1[7]*p2[5]+p1[8]*p2[8];
}
inline void mm31(double *p1, double *p2, double *p3)
{
	p3[0]=p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
	p3[1]=p1[3]*p2[0]+p1[4]*p2[1]+p1[5]*p2[2];
	p3[2]=p1[6]*p2[0]+p1[7]*p2[1]+p1[8]*p2[2];
}
*/
/*Next function is for caculating the refraction of a ball*/


void refraction(double *ptheta1, double *pphi1, double *ptheta2, double *pphi2, const double n) 
{
	double rot[9], DIR[3], sc1[3], sc2[3], theta2_, theta3, theta4, arot[9], cache[3], xy;
	/*rot[], DIR[] are in Cartesian coordinate, sc1[] and sc2[] are in spherical coordinate*/
	double rot1[9]={cos(*pphi1), -sin(*pphi1), 0, sin(*pphi1), cos(*pphi1), 0, 0, 0, 1}, 
	       rot2[9]={cos(*ptheta1), 0, sin(*ptheta1), 0, 1, 0, -sin(*ptheta1), 0, cos(*ptheta1)};	
	mm3(rot1, rot2, rot);
	arot[0]=rot[0], arot[1]=rot[3], arot[2]=rot[6], arot[3]=rot[1], arot[4]=rot[4], arot[5]=rot[7], arot[6]=rot[2], arot[7]=rot[5], arot[8]=rot[8];
	cache[0]=sin(*ptheta2)*cos(*pphi2), cache[1]=sin(*ptheta2)*sin(*pphi2), cache[2]=cos(*ptheta2);

	/*cache is in Cartesian coordinate*/
	mm31(arot, cache, DIR);/*DIR is the direction vector in the rotated coordinate*/
	if (hypot(DIR[0],DIR[1])/n < 1) {
		car2sph(DIR, sc1);
		theta2_=asin(hypot(DIR[0],DIR[1])/n);
		theta3=2*(pi/2-theta2_);
		theta4=pi-sc1[1]+theta3;
		
		sc1[1]=theta3;
		sph2car(sc1, cache);
		mm31(rot, cache, DIR); 
		car2sph(DIR, sc2);
		*ptheta1=sc2[1],*pphi1=sc2[2];
		
		sc1[1]=theta4;
		sph2car(sc1, cache);
		mm31(rot, cache, DIR);
		car2sph(DIR, sc2);
		*ptheta2=sc2[1],*pphi2=sc2[2];	
	}
	else {
		DIR[2] = -DIR[2];
		mm31(rot, DIR, cache);
		car2sph(cache, sc2);
		*ptheta2=sc2[1],*pphi2=sc2[2];	
	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  readpho
 *  Description:  Read a 5 version of photon file for given numbers of photons
 * =====================================================================================
 */
int readpho(FILE *fppho, double *pho, const unsigned long number_of_photons) {
	unsigned long i;
	double data[5];
	for (i = 0; i < number_of_photons; i++)	{
		fread(data, sizeof(double), 5, fppho);
		pho[i*7] = data[0], pho[i*7+1] = data[1],
		pho[i*7+2] = 0,
		pho[i*7+3] = data[2], pho[i*7+4] = data[3],
		pho[i*7+5] = sqrt(1 - pow(data[2], 2) - pow(data[3], 2)),
		pho[i*7+6] = data[4];
	}
	return 0;
}		/* -----  end of function readpho  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  writepho
 *  Description:  Write given numbers of photons into a given file
 * =====================================================================================
 */
int writepho(FILE *fpout, double *pho, const unsigned long number_of_photons) {
	unsigned long i;
	double data[5];
	for (i = 0; i < number_of_photons; i++)	{
		data[0] = pho[i*7], data[1] = pho[i*7+1],
		data[2] = pho[i*7+3], data[3] = pho[i*7+4],
		data[4] = pho[i*7+6];
		fwrite(data, sizeof(double), 5, fpout);
	}
	fflush(fpout);
	return 0;
}		/* -----  end of function writepho  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  IndexPos
 *  Description:  Indexing the position array. The array to be indexed is delivered with
 *	*poss, the indexed array is posd. Every block has a size of size_per_block. The size
 *	of the array posd should be equal to xc*yc*size_per_block
 * =====================================================================================
 */

int IndexPos(double xs, double xe, unsigned long xc, double ys, double ye, unsigned long yc, double *poss, unsigned long size_s, double *posd, unsigned long *end_position, unsigned long size_per_block) {

	unsigned long N = xc*yc, i, x, y;
	double X, Y;

	for (i = 0; i < N; i++) {
		end_position[i] = 0;
	}
	for (i = 0; i < size_per_block*3*xc*yc; i++)
		posd[i] = nan("");

	for (i = 0; i < size_s; i++) {

		X = poss[i*3], Y = poss[i*3+1];
		if (X < xs || X > xe || Y < ys || Y > ye) {
			fprintf(stderr, "Size error in IndexPos()\n");
			exit(1);
		}
		
		x = floor((X-xs)/(xe-xs)*xc), y = floor((Y-ys)/(ye-ys)*yc);

		if (x == xc)
			x = xc-1;
		if (y == yc)
			y = yc-1;

		
		if(!isnan(posd[size_per_block*(y*xc+x)*3+end_position[y*xc+x]*3+2])) {
	//		fprintf(stderr, "Cache not large enough, index failed, %lf\n", posd[size_per_block*(y*xc+x)+end_position[y*xc+x]*3+2]);
			return 1;
		}

		posd[size_per_block*3*(y*xc+x)+end_position[y*xc+x]*3] = poss[i*3];
		posd[size_per_block*3*(y*xc+x)+end_position[y*xc+x]*3+1] = poss[i*3+1];
		posd[size_per_block*3*(y*xc+x)+end_position[y*xc+x]*3+2] = poss[i*3+2];
		end_position[y*xc+x]++;
		
	}
	
/*
	for (x = 0; x < xc; x++) {
		for (y = 0; y < yc; y++) 
			printf("%ld, ", end_position[x+y*xc]);
		printf("\n");
	}
*/
	return 0;
}		/* -----  end of function IndexPos  ----- */

