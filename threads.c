#define pi 3.141592653589793
#define SIZE_OF_CACHE 7200
#define SIZE_OF_INDEX 1200
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>
#include"REF.h"
#include<time.h>
#include<string.h>

static unsigned long FindIndex(const double xs, const double xe, const long xc, const double ys, const double ye, const long yc, const double x0, const double y0, const double r0, unsigned long *index, const unsigned long size_of_indices);
static unsigned long cou(pthread_mutex_t *mutex);
static inline int cache_compare(const void *a, const void *b);
static int clean_cache(double *cache, unsigned long size);
//static unsigned long copy_pos(double *cache, double *pos, unsigned long *pos_index, unsigned long *end_position, unsigned long size_per_block, unsigned long number_of_indices);

void *photons(void *ARG)/*parameter is the number of the balls*/
{	
	double d, d_, d1, d2, d3, sc1[3], sc2[3], x0, y0, r0, r2, X[3], X0[3], a, b, c, r, dis;
	unsigned long i, j, k, eM, index_count, size_per_block;	
	unsigned long *end_position;
	unsigned long pos_index[SIZE_OF_INDEX], number_of_indices;
	double *cache2, *pho, *pos;
	char is_using_dynamic_cache, is_within;
	double cache1[SIZE_OF_CACHE];
	double xy, n, dL;
	long stat1, stat3, stat4, stat5;
	double stat2;
	a=(((struct ref_arg *)ARG)->ref_DIR).geo_a;
	b=(((struct ref_arg *)ARG)->ref_DIR).geo_b;
	c=(((struct ref_arg *)ARG)->ref_DIR).geo_c;
	r=(((struct ref_arg *)ARG)->ref_DIR).geo_r;
	n=(((struct ref_arg *)ARG)->ref_DIR).ri_n;
	dL=(((struct ref_arg *)ARG)->ref_DIR).ri_dL;

	index_count=((struct ref_arg *)ARG)->index_count;
	size_per_block=((struct ref_arg *)ARG)->size_per_block;
	end_position=((struct ref_arg *)ARG)->end_position;

	pho=((struct ref_arg *)ARG)->ref_pho;
	pos=((struct ref_arg *)ARG)->ref_pos;

	eM=floor(c/r*40);
	
	cache2 = malloc(eM*sizeof(double)*3);
	r2=0.999*r;
	
	for(i=((struct ref_arg *)ARG)->phos;i<((struct ref_arg *)ARG)->phoe;i++)
	{	
		cou(((struct ref_arg *)ARG)->ref_pmutex);

		x0=pho[i*7], y0=pho[i*7+1];
		r0=r2;
		is_within=1;
		stat1 = 0; /*spheres that passes*/
		stat2 = 0; /*spatial diversion*/
		stat3 = -1; /*cache expasion times*/
		stat4 = 0; /*cache size*/
		stat5 = 0; /*if dynamic cache is applied*/
		while(is_within)
		{
			pho[i*7]=x0,pho[i*7+1]=y0,pho[i*7+2]=0;
			pho[i*7+3]=0,pho[i*7+4]=0,pho[i*7+5]=1;
			dis = 0;
			is_within = 0;
			r0 = r0+r2;
			stat2 = r0;
			stat3++;
			k = 0;
			number_of_indices = FindIndex(0, a+4*r, index_count, 0, b+4*r, index_count, x0, y0, r0, pos_index, SIZE_OF_INDEX);
					
			for (j = 0; j < number_of_indices; j++)
				k += end_position[pos_index[j]];
			if(k > SIZE_OF_CACHE/3) {
				is_using_dynamic_cache = 1;
				for (j = 0, k = 0; j < number_of_indices; j++) {
					memcpy(cache2+k*3, pos+pos_index[j]*3*size_per_block, 3*sizeof(double)*end_position[pos_index[j]]);
					k += end_position[pos_index[j]];
				}
//				k = copy_pos(cache2, pos, pos_index, end_position, size_per_block, number_of_indices);
				qsort(cache2, k, sizeof(double)*3, cache_compare);
			}
			else {
				is_using_dynamic_cache = 0;
				for (j = 0, k = 0; j < number_of_indices; j++) {
					memcpy(cache1+k*3, pos+pos_index[j]*3*size_per_block, 3*sizeof(double)*end_position[pos_index[j]]);
					k += end_position[pos_index[j]];
				}
//				k = copy_pos(cache1, pos, pos_index, end_position, size_per_block, number_of_indices);
				qsort(cache1, k, sizeof(double)*3, cache_compare);
			}
			stat4 = k;
			stat5 = is_using_dynamic_cache;
			
			if(is_using_dynamic_cache){
				for(j=0;j<k;j++){				
					d1=pho[i*7+3]*(cache2[j*3]-pho[i*7])+pho[i*7+4]*(cache2[j*3+1]-pho[i*7+1])+pho[i*7+5]*(cache2[j*3+2]-pho[i*7+2]);
				
					if(d1>0){					
						d2=hypot(hypot(pho[i*7]-cache2[j*3],pho[i*7+1]-cache2[j*3+1]),pho[i*7+2]-cache2[j*3+2]);
					
						d3=sqrt(d2*d2-d1*d1);
						if(d3<r2){			
							d=d1-sqrt(r2*r2-d3*d3);
							X[0]=pho[i*7]+pho[i*7+3]*d-cache2[j*3];
							X[1]=pho[i*7+1]+pho[i*7+4]*d-cache2[j*3+1];
							X[2]=pho[i*7+2]+pho[i*7+5]*d-cache2[j*3+2];
							memcpy(X0, X, sizeof(double)*3);
							car2sph(X,sc1);
							car2sph(pho+i*7+3,sc2);
						
							refraction(sc1+1,sc1+2,sc2+1,sc2+2,n);
							stat1++;
							sph2car(sc1,X);

							d_=sqrt(pow(X[0]-X0[0],2)+pow(X[1]-X0[1],2)+pow(X[2]-X0[2],2));
							if(dL > 0)
								dis += d;
							else
								dis += d_;

							sph2car(sc2,pho+i*7+3);
							pho[i*7]=cache2[j*3]+X[0];
							pho[i*7+1]=cache2[j*3+1]+X[1];
							pho[i*7+2]=cache2[j*3+2]+X[2];						
						
						}
					}
				}
				clean_cache(cache2, k*3);
			}
			else{
				for(j=0;j<k;j++){				
					d1=pho[i*7+3]*(cache1[j*3]-pho[i*7])+pho[i*7+4]*(cache1[j*3+1]-pho[i*7+1])+pho[i*7+5]*(cache1[j*3+2]-pho[i*7+2]);
				
					if(d1>0){					
						d2=hypot(hypot(pho[i*7]-cache1[j*3],pho[i*7+1]-cache1[j*3+1]),pho[i*7+2]-cache1[j*3+2]);
					
						d3=sqrt(d2*d2-d1*d1);
						if(d3<r2){			
							d=d1-sqrt(r2*r2-d3*d3);
							X[0]=pho[i*7]+pho[i*7+3]*d-cache1[j*3];
							X[1]=pho[i*7+1]+pho[i*7+4]*d-cache1[j*3+1];
							X[2]=pho[i*7+2]+pho[i*7+5]*d-cache1[j*3+2];
							memcpy(X0, X, sizeof(double)*3);
							car2sph(X,sc1);
							car2sph(pho+i*7+3,sc2);
						
							refraction(sc1+1,sc1+2,sc2+1,sc2+2,n);
							stat1++;
							sph2car(sc1,X);

							d_=sqrt(pow(X[0]-X0[0],2)+pow(X[1]-X0[1],2)+pow(X[2]-X0[2],2));
							if(dL > 0)
								dis += d;
							else
								dis += d_;

							sph2car(sc2,pho+i*7+3);
							pho[i*7]=cache1[j*3]+X[0];
							pho[i*7+1]=cache1[j*3+1]+X[1];
							pho[i*7+2]=cache1[j*3+2]+X[2];						
						}
					}
				}
				clean_cache(cache1, k*3);
			}
			
			
			is_within=0;
			if (dL > 0)
				dis+=(c-pho[i*7+2])/pho[i*7+5];
			pho[i*7]=(c-pho[i*7+2])*pho[i*7+3]/pho[i*7+5]+pho[i*7];
			pho[i*7+1]=(c-pho[i*7+2])*pho[i*7+4]/pho[i*7+5]+pho[i*7+1];
			if(hypot(x0-pho[i*7],y0-pho[i*7+1])>(r0-r2))
				is_within=1;
				
		}
		

		if (debug) {
			pho[i*7] = (double)stat1,
			pho[i*7+1] = stat2,
			pho[i*7+3] = (double)stat3,
			pho[i*7+4] = (double)stat4,
			pho[i*7+6] = (double)stat5;
		}
		else {
			pho[i*7]=pho[i*7]-2*r,pho[i*7+1]=pho[i*7+1]-2*r;
			pho[i*7+6]*=exp(-dis/fabs(dL));
			car2sph(pho+i*7+3,sc1);
			sc1[1]=asin(sin(sc1[1])/n);
			sph2car(sc1,pho+i*7+3);
		}

	}
	free(cache2);
	return NULL;
	
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  FindIndex
 *  Description:  find the index of spheres around a certain photon x0, y0. Returns the 
 *	amount of indices overall.
 * =====================================================================================
 */

unsigned long FindIndex(const double xs, const double xe, const long xc, const double ys, const double ye, const long yc, const double x0, const double y0, const double r0, unsigned long *index, const unsigned long size_of_indices) {
	const double X1 = (x0-r0-xs)/(xe-xs)*xc, X2 = (x0+r0-xs)/(xe-xs)*xc, Y1 = (y0-r0-xs)/(xe-xs)*xc, Y2 = (y0+r0-xs)/(xe-xs)*xc, R0 = r0/(xe-xs)*xc, X0 = (x0-xs)/(xe-xs)*xc, Y0 = (y0-ys)/(ye-ys)*yc;
	long i = 0, j, k;

	for (j = (long)X1; j <= (long)X2; j++) {
		if (i == size_of_indices)
			return 0;
		if(0 <= j && j < xc)
			index[i++] = j + (long)Y0*xc;
	}
	

	for (k = (long)Y1; k < (long)Y0; k++) {
		if (i == size_of_indices)
			return 0;
		if(0 <= k && k < yc)
			index[i++] = X0 + (long)k*xc;
	}

	for (k = (long)Y0+1; k <= (long)Y2; k++) {
		if (i == size_of_indices)
			return 0;
		if(0 <= k && k < yc)
			index[i++] = X0 + (long)k*xc;
	}


	for (j = (long)X1; j < (long)X0; j++) {
		for (k = (long)Y0 + 1; k <= (long)Y2; k++) {
			if (hypot(-X0+j+1, -Y0+k) < R0) {
				if (i == size_of_indices)
					return 0;

				if(0 <= k && k < yc && 0 <= j && j < xc)
					index[i++] = k*xc+j;	
			}
		}
	}
	
	for (j = (long)X0+1; j <= (long)X2; j++) {
		for (k = (long)Y0 + 1; k <= (long)Y2; k++) {
			if (hypot(-X0+j, -Y0+k) < R0) {
				if (i == size_of_indices)
					return 0;
				if(0 <= k && k < yc && 0 <= j && j < xc)
					index[i++] = k*xc+j;	
			}
		}
	}

	for (j = (long)X1; j < (long)X0; j++) {
		for (k = (long)Y1; k < (long)Y0; k++) {
			if (hypot(-X0+j+1, -Y0+k+1) < R0) {
				if (i == size_of_indices)
					return 0;
				if(0 <= k && k < yc && 0 <= j && j < xc)
					index[i++] = k*xc+j;	
			}
		}
	}
	
	for (j = (long)X0+1; j <= (long)X2; j++) {
		for (k = (long)Y1; k < (long)Y0; k++) {
			if (hypot(-X0+j, -Y0+k+1) < R0) {
				if (i == size_of_indices)
					return 0;
				if(0 <= k && k < yc && 0 <= j && j < xc)
					index[i++] = k*xc+j;	
			}
		}
	}
	

	return i;
}		/* -----  end of function find_index  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cou
 *  Description:  using for counting photons computed
 * =====================================================================================
 */
static unsigned long cou(pthread_mutex_t *mutex){
	FILE *flog;
	pthread_mutex_lock(mutex);
	if(++count%100000==0){
		flog=fopen(str_log,"a");
		fprintf(flog,"%ld phos has been computed, cost %ld seconds\n",count, time(0)-timer);
		timer=time(0);
		fclose(flog);
	/*	flog=fopen("cont","r");
		while(fgetc(flog)!='R'){
			fclose(flog);
			sleep(1200);
			flog=fopen("cont","r");
		}*/
	}
	pthread_mutex_unlock(mutex);
	return 0;
}		/* -----  end of function cou  ----- */
/*
	static unsigned long find(double *pos, double *cache1, double *cache2, const double x0, const double y0, const double r0, long k, const long number_of_positions, char p_is_using_dynamic_cache) {
		long j, l;
		*p_is_using_dynamic_cache = 0;
		for(j = 0; j < number_of_positions; j++){		
			if(hypot(pos[j*3]-x0,pos[j*3+1]-y0)<r0) {			
				cache1[k*3]=pos[j*3],cache1[k*3+1]=pos[j*3+1],cache1[k*3+2]=pos[j*3+2];
				k++;
			}
			if(k == 2399){
				for(l=0;l<2400;l++){
					cache2[l*3]=cache1[l*3],cache2[l*3+1]=cache1[l*3+1],cache2[l*3+2]=cache1[l*3+2];
				}
				for(j++;j<number_of_positions;j++)
					if(hypot(pos[j*3]-x0,pos[j*3+1]-y0)<r0) {
						cache2[k*3]=pos[j*3],cache2[k*3+1]=pos[j*3+1],cache2[k*3+2]=pos[j*3+2];
						k++;
					}
				*p_is_using_dynamic_cache=1;
				break;
			}
		}
		return k;
	}
*/

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  cache_compare
 *  Description:  
 * =====================================================================================
 */
static inline int cache_compare(const void *a, const void *b) {
	if (*((double *)a + 2) - *((double *)b + 2) > 0)
		return 1;
	else
		return -1;
}
		/* -----  end of function cache_compare  ----- */
static int clean_cache(double *cache, unsigned long size) {
	for (;size > 0; size--)
		cache[size] = 0;
	return 0;
}

/*
static unsigned long copy_pos(double *cache, double *pos, unsigned long *pos_index, unsigned long *end_position, unsigned long size_per_block, unsigned long number_of_indices) {
	unsigned long k, j;
	for (j = 0, k = 0; j < number_of_indices; j++) {
		memcpy(cache+k*3, pos+pos_index[j]*3*size_per_block, 3*sizeof(double)*end_position[pos_index[j]]);
		k += end_position[pos_index[j]];
	}
	return k;
}
*/
