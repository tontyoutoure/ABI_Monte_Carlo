/*REF.h*/
#define pi 3.141592653589793
#define THREAD_MAX 100
#define sqrt2 1.4142135623730950488016887242097
#define debug 0

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>
#include<stdint.h>
uint32_t number_passes_tumor;

struct ref_dir{
	double geo_a;
	double geo_b;
	double geo_c;
	double geo_r;
	/*the geometry information of the tissue*/
	double ri_n;
	double ri_dL;
	char L;
	/*refraction index*/
};

struct ref_arg{
	unsigned long number_of_positions;
	unsigned long phos;
	unsigned long phoe;
	unsigned long index_count;
	unsigned long size_per_block;
	unsigned long *end_position;
	double *ref_pos;
	double *ref_pho;
	pthread_mutex_t *ref_pmutex;
	struct ref_dir ref_DIR;
};


#define mm3(p1, p2, p3)\
do {\
	(p3)[0]=(p1)[0]*(p2)[0]+(p1)[1]*(p2)[3]+(p1)[2]*(p2)[6];\
	(p3)[1]=(p1)[0]*(p2)[1]+(p1)[1]*(p2)[4]+(p1)[2]*(p2)[7];\
	(p3)[2]=(p1)[0]*(p2)[2]+(p1)[1]*(p2)[5]+(p1)[2]*(p2)[8];\
	(p3)[3]=(p1)[3]*(p2)[0]+(p1)[4]*(p2)[3]+(p1)[5]*(p2)[6];\
	(p3)[4]=(p1)[3]*(p2)[1]+(p1)[4]*(p2)[4]+(p1)[5]*(p2)[7];\
	(p3)[5]=(p1)[3]*(p2)[2]+(p1)[4]*(p2)[5]+(p1)[5]*(p2)[8];\
	(p3)[6]=(p1)[6]*(p2)[0]+(p1)[7]*(p2)[3]+(p1)[8]*(p2)[6];\
	(p3)[7]=(p1)[6]*(p2)[1]+(p1)[7]*(p2)[4]+(p1)[8]*(p2)[7];\
	(p3)[8]=(p1)[6]*(p2)[2]+(p1)[7]*(p2)[5]+(p1)[8]*(p2)[8];\
}\
while(0)

#define mm31(p1, p2, p3)\
do {\
	(p3)[0]=(p1)[0]*(p2)[0]+(p1)[1]*(p2)[1]+(p1)[2]*(p2)[2];\
	(p3)[1]=(p1)[3]*(p2)[0]+(p1)[4]*(p2)[1]+(p1)[5]*(p2)[2];\
	(p3)[2]=(p1)[6]*(p2)[0]+(p1)[7]*(p2)[1]+(p1)[8]*(p2)[2];\
}\
while(0)

#define sph2car(ps, pc)\
do {\
	(pc)[0] = (ps)[0]*sin((ps)[1])*cos((ps)[2]);\
	(pc)[1] = (ps)[0]*sin((ps)[1])*sin((ps)[2]);\
	(pc)[2] = (ps)[0]*cos((ps)[1]);\
}\
while(0)

#define car2sph(pc, ps)\
do {\
	xy = hypot((pc)[0], (pc)[1]);\
	(ps)[0]=hypot(xy, (pc)[2]);\
	if(!(ps)[0]) {\
		(ps)[1]=0;\
		(ps)[2]=0;\
	}\
	else if(!xy) {\
		(ps)[2]=0;\
		if((pc)[2]>0)\
			(ps)[1]=0;\
		else\
			(ps)[1]=pi;\
	}\
	else {\
		(ps)[1]=acos((pc)[2]/(ps)[0]);\
		if((pc)[1]>0)\
			(ps)[2]=acos((pc)[0]/xy);\
		else\
			(ps)[2]=acos((-(pc)[0])/xy)+pi;\
	}\
}\
while(0)


/*Global variables*/
extern time_t timer;
extern char str_log[FILENAME_MAX];
extern unsigned long count;


/* functions.c */
/*
void sph2car(double *ps, double *pc);
void car2sph(double *pc, double *ps);
*/
void refraction(double *ptheta1, double *pphi1, double *ptheta2, double *pphi2, const double n);
int readpho(FILE *fppho, double *pho, const unsigned long number_of_photons);
int writepho(FILE *fpout, double *pho, const unsigned long number_of_photons);
int IndexPos(double xs, double xe, unsigned long xc, double ys, double ye, unsigned long yc, double *poss, unsigned long size_s, double *posd, unsigned long *end_position, unsigned long size_per_block);

/* threads_fast.c */
void *photons(void *ARG);
