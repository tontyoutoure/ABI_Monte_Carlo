#include<stdio.h>
#include<string.h>
#include<stdlib.h>
//#include<math.h>
#include<time.h>
#include"mt19937ar.c"
#include <gsl_cdf.h>

#define MODE_PLANE_PARALLEL 1
#define MODE_GAUSSIAN_DISTRIBUTE 2

int main(int argc, char *argv[]){
	long i, n;
	double xf, xc, yf, yc, sigma_x, sigma_y;
	int mode, argn=0;
	FILE *fp=NULL;
	double *pho;

	for(i=1;i<argc;i++){
		if(!strcmp(argv[i],"-x")){
			i++;
			if(!sscanf(argv[i],"%lf",&xf)){
				printf("x input error");
				return -1;
			}
			i++;
			if(!sscanf(argv[i],"%lf",&xc)){
				printf("x input error");
				return -1;
			}
			argn++;
		}
		else if(!strcmp(argv[i],"-y")){
			i++;
			if(!sscanf(argv[i],"%lf",&yf)){
				printf("y input error");
				return -1;
			}
			i++;
			if(!sscanf(argv[i],"%lf",&yc)){
				printf("y input error");
				return -1;
			}
			argn++;
		}
		else if(!strcmp(argv[i],"--mode-pp")){
			mode=MODE_PLANE_PARALLEL;
			argn++;
		}
		else if(!strcmp(argv[i],"--mode-gaussian")) {
			i++;
			if(!sscanf(argv[i],"%lf",&sigma_x)){
				printf("sigma input error");
				return -1;
			}
			i++;
			if(!sscanf(argv[i],"%lf",&sigma_y)){
				printf("sigma input error");
				return -1;
			}
			mode=MODE_GAUSSIAN_DISTRIBUTE;
			argn++;
		}
		else if(!strcmp(argv[i],"-n")){
			i++;
			if(!sscanf(argv[i],"%ld",&n)){
				printf("number of photons input error");
				return -1;
			}
			argn++;
		}
		else if(!strcmp(argv[i],"-o")){
			i++;
			fp=fopen(argv[i],"w");
			if(fp==NULL){
				printf("output file cannot open\n");
				return -1;
			}
			argn++;
		}
		else{
			printf("Input error\n");
			return -1;
		}
	}
	if(argn<5){
		printf("not enough parameter, plz check\n");
		return -2;
	}
	init_genrand(time(0));
	
	if(mode==MODE_PLANE_PARALLEL){
		pho=(double *)malloc(n*sizeof(double)*5);
		for(i=0;i<n;i++){
			pho[i*5]=genrand_res53()*(xc-xf)+xf;
			pho[i*5+1]=genrand_res53()*(yc-yf)+yf;
			pho[i*5+2]=0, pho[i*5+3]=0, pho[i*5+4]=1;
		}
		fwrite(pho,sizeof(double),5*n,fp);
		printf("%ld photons has been generated\n",n);
		fclose(fp);
	}
	else if (mode == MODE_GAUSSIAN_DISTRIBUTE) {
		pho=(double *)malloc(n*sizeof(double)*5);
		for (i = 0; i < n; i++) {
			pho[i*5]=genrand_res53()*(xc-xf)+xf;
			pho[i*5+1]=genrand_res53()*(yc-yf)+yf;
			pho[i*5+2] = gsl_cdf_gaussian_Pinv(genrand_res53(), sigma_x);
			pho[i*5+3] = gsl_cdf_gaussian_Pinv(genrand_res53(), sigma_y);
			pho[i*5+4]=1;
		}
		fwrite(pho,sizeof(double),5*n,fp);
		printf("%ld photons has been generated\n",n);
		fclose(fp);
	}


	return 0;
}
