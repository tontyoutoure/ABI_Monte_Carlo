#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>
#include"REF.h"
#include<string.h>
#include<time.h>
#include<argp.h>

static int REF(struct ref_dir DIR, int number_of_threads, double *pos, double *pho, unsigned long number_of_photons, unsigned long number_of_positions, unsigned long index_count, unsigned long size_per_block, unsigned long *end_position);
static int parse_opt(int key, char *arg, struct argp_state *state);

#define IS_A 1
#define IS_B 2
#define IS_C 4
#define IS_R 8
#define IS_I 16
#define IS_O 32
#define IS_T 64
#define IS_N 256
#define IS_H 512
#define IS_L 1024
#define IS_D 2048

static long input_check = 0;
const char *argp_program_version = "version 2.1";
const char *argp_program_bug_address = "tontyoutoure@gmail.com";
time_t timer;
char str_log[FILENAME_MAX];
unsigned long count;

struct argp_arg {
	double *a;
	double *b;
	double *c;
	double *r;
	double *n;
	double *dL;
	char *L;
	unsigned long *NoT;
	unsigned long *index_count;
	char *str_in;
	char *str_tissue;
	char *str_out;
	char *str_log;
};

static struct argp_option opt[] = {
	{NULL, 'a', "length", 0, "Set the length of a (*)"},
	{NULL, 'b', "length", 0, "Set the length of b (*)"},
	{NULL, 'c', "length", 0, "Set the length of c (*)"},
	{"radius", 'r', "length", 0, "Set the length of radius of little spheres (*)"},
	{"energy", 'e', "length", 0, "Specify the energy (*)"},
	{"attenuation-Length", 'L', "length", 0, "Set attenutation length directly"},
	{"refraction-index", 'n', "refractive index", 0, "Set refractive index directly"},
	{"input-file", 'i', "filename", 0, "Place the input file name (*)"},
	{"output-file", 'o', "filename", 0, "Place the output file name (*)"},
	{"log-file", 'l', "filename", 0, "Place the log file name (*)"},
	{"tissue-file", 't', "filename", 0, "Place the tissue file name (*)"},
	{"number-of-threads", 'h', "number", 0, "Place the number of threads (*)"},
	{"index-count", 'u', "number", 0, "Number of indexes of the looking table, each direction"},
	{"lenGth", 'G', NULL, 0, "Record length instead of attenuation factor"},
	{0}
};

int main(int argc, char *argv[]){
	char str_in[FILENAME_MAX], str_tissue[FILENAME_MAX], str_out[FILENAME_MAX];
	struct ref_dir DIR={0, 0, 0, 0, 0, 0, 0};
	unsigned long number_of_threads, number_of_photons, number_of_positions, index_count = 100, size_per_block;
	FILE *flog, *fppho, *fppos;
	double *pho, *pos, *posd;
	time_t timero; 
	unsigned long *end_position;
	double size_m = 2;

	struct argp argp = {opt, parse_opt, NULL, "A monte carlo program for the caculation of DEI of lung tissue\n"
	"Options with (*) should not be omitted"};
	struct argp_arg arg = {&DIR.geo_a, &DIR.geo_b, &DIR.geo_c, &DIR.geo_r, &DIR.ri_n, &DIR.ri_dL, &DIR.L, &number_of_threads, &index_count, str_in, str_tissue, str_out, str_log};

	/*initiating*/
	timer=time(0);
	timero = time(0);
	str_in[0]='\0';

	/*reading parameters*/
	
	argp_parse(&argp, argc, argv, 0, 0, &arg);

	if((input_check | (IS_A | IS_B | IS_C | IS_R | IS_O | IS_T | IS_N | IS_H | IS_L | IS_D)) != input_check){
		fprintf(stderr, "Maybe you lost a few inputs\n");
		return -2;
	}
	
	/*Read sphere positions*/

	fppos = fopen(str_tissue, "r");
	if(fppos == NULL) {
		fprintf(stderr, "Tissue file cannot be opened\n");
		return -2;
	}
	fseek(fppos, 0, SEEK_END);
	number_of_positions=ftell(fppos)/3/sizeof(double);
	rewind(fppos);
	pos = malloc(sizeof(double)*number_of_positions*3);
	fread(pos, sizeof(double), number_of_positions*3, fppos);
	fclose(fppos);
	
	/*Read input photons*/
	fppho = fopen(str_in, "r");
	if(fppho == NULL) {
		fprintf(stderr, "Input file cannot be opened\n");
		return -2;
	}

	fseek(fppho, 0, SEEK_END);
	number_of_photons = ftell(fppho)/5/sizeof(double);
	rewind(fppho);
	pho = malloc(sizeof(double)*7*number_of_photons);
	readpho(fppho, pho, number_of_photons);
	fclose(fppho);

	/*Originally, the program will generate a plane parallel light.
	 *This part of code is reserved for reference.
	 */

	/*
		else{
			pho=(double **)malloc(sizeof(double *)*number_of_photons);
			pho[0]=(double *)malloc(sizeof(double)*number_of_photons*7);
			for(i=1;i<number_of_photons;i++)
				pho[i]=pho[i-1]+7;
			for(i=0;i<number_of_photons;i++){
				pho[i][0]=genrand_res53()*DIR.geo_a+2*DIR.geo_r;
				pho[i][1]=genrand_res53()*DIR.geo_b+2*DIR.geo_r;
				pho[i][2]=0;
				pho[i][3]=0, pho[i][4]=0, pho[i][5]=1;
				pho[i][6]=1;
			}
		}
	 */
	end_position = malloc(sizeof(long)*index_count*index_count);
	
	while(1) {
		size_per_block = ceil((double)number_of_positions*size_m/index_count/index_count);
		posd = malloc(size_per_block*3*sizeof(double)*index_count*index_count);
		if (!IndexPos(0, DIR.geo_r*4+DIR.geo_a, index_count, 0, DIR.geo_r*4+DIR.geo_b, index_count, pos, number_of_positions, posd, end_position, size_per_block))
			break;
		free(posd);
		size_m += 1;
	}

	REF(DIR, number_of_threads, posd, pho, number_of_photons, number_of_positions, index_count, size_per_block, end_position);
	
	fppho=fopen(str_out,"a");//open a specific file for output
	writepho(fppho, pho, number_of_photons);
	fclose(fppho);
	printf("Data written to disk\n");
	
	free(pos);
	free(pho);
	free(posd);
	free(end_position);
	
	flog=fopen(str_log,"a");
	fprintf(flog,"%ld seconds has been passed\n",time(0)-timero);
	printf("%ld seconds has been passed\n",time(0)-timero);
	fclose(flog);

	return 0;
}

static int REF(struct ref_dir DIR, int number_of_threads, double *pos, double *pho, unsigned long number_of_photons, unsigned long number_of_positions, unsigned long index_count, unsigned long size_per_block, unsigned long *end_position){
	pthread_t tid[100];
	struct ref_arg ARG[THREAD_MAX];
	pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
	unsigned long i;
	
/*Initialization*/	
	for(i=0;i<number_of_threads-1;i++){
		(ARG+i)->number_of_positions=number_of_positions;
		(ARG+i)->phos=i*floor((double)number_of_photons/number_of_threads);
		(ARG+i)->phoe=(i+1)*floor((double)number_of_photons/number_of_threads);//this thread will caculate from pho[i*(...)] to pho[(i+1)*(...)-1]
		(ARG+i)->index_count = index_count;
		(ARG+i)->size_per_block = size_per_block;
		(ARG+i)->end_position = end_position;
		(ARG+i)->ref_DIR=DIR;
		(ARG+i)->ref_pos=pos;
		(ARG+i)->ref_pho=pho;
		(ARG+i)->ref_pmutex=&mutex;
		pthread_create(tid+i,NULL,photons,(void *)(ARG+i));
	}
	(ARG+i)->number_of_positions=number_of_positions;
	(ARG+i)->phos=i*floor((double)number_of_photons/number_of_threads);
	(ARG+i)->phoe=number_of_photons;
	(ARG+i)->index_count = index_count;
	(ARG+i)->size_per_block = size_per_block;
	(ARG+i)->end_position = end_position;
	(ARG+i)->ref_DIR=DIR;
	(ARG+i)->ref_pos=pos;
	(ARG+i)->ref_pho=pho;
	(ARG+i)->ref_pmutex=&mutex;
	pthread_create(tid+i,NULL,photons,(void *)(ARG+i));
	
	for(i=0;i<number_of_threads;i++)
		pthread_join(tid[i],NULL);
	printf("Simulation done\n");
	return 0;
}

static int parse_opt(int key, char *arg, struct argp_state *state) {
	FILE *fptest;
	struct argp_arg *input = state -> input;
	switch (key) {
		case 'a': 
			if(atof(arg) > 0) {
				*(input -> a) = atof(arg);
			}
			else {
				argp_failure(state, 1, 0, "Input error with a");
			}
			input_check |= IS_A;
			break;
		case 'b':
			if(atof(arg) > 0) {
				*(input -> b) = atof(arg);
			}
			else {
				argp_failure(state, 1, 0, "Input error with b");
			}
			input_check |= IS_B;
			break;
		case 'c':
			if(atof(arg) > 0) {
				*(input -> c) = atof(arg);
			}
			else {
				argp_failure(state, 1, 0, "Input error with c");
			}
			input_check |= IS_C;
			break;
		case 'r':
			if(atof(arg) > 0) {
				*(input -> r) = atof(arg);
			}
			else {
				argp_failure(state, 1, 0, "Input error with radius");
			}
			input_check |= IS_R;
			break;
		case 'i':
			if((fptest = fopen(arg, "r")) != NULL) {
				strcpy(input -> str_in, arg);
			}
			else {
				argp_failure(state, 1, 0, "Input file error");
			}
			input_check |= IS_I;
			break;
		case 'o':
			if((fptest = fopen(arg, "wx")) != NULL) {
				strcpy(input -> str_out, arg);
				fclose(fptest);
			}
			else {
				argp_failure(state, 1, 0, "Output file error");
			}
			input_check |= IS_O;
			break;
		case 'l':
			if((fptest = fopen(arg, "a")) != NULL) {
				strcpy(input -> str_log, arg);
				fclose(fptest);
			}
			else {
				argp_failure(state, 1, 0, "Log file error");
			}
			input_check |= IS_L;
			break;
		case 't':
			if((fptest = fopen(arg, "r")) != NULL) {
				strcpy(input -> str_tissue, arg);
				fclose(fptest);
			}
			else {
				argp_failure(state, 1, 0, "Tissue file error");
			}
			input_check |= IS_T;
			break;
		case 'h':
			if(atoi(arg) > 0 && atoi(arg) < THREAD_MAX) {
				*(input -> NoT) = atoi(arg);
			}
			else {
				argp_failure(state, 1, 0, "Input error with number of threads");
			}
			input_check |= IS_H;
			break;
		case 'L':
			*(input -> dL) = atof(arg);
			input_check |= IS_D;
			break;
		case 'G':
			*(input -> L) = 1;
			break;
		case 'n':
			*(input -> n) = atof(arg);
			input_check |= IS_N;
			break;
		case 'e':
			switch (atoi(arg)) {
				case 15:
					*(input -> n) = 1.000001025922213;
					*(input -> dL) = 0.005534;
					break;
				case 20:
					*(input -> n) = 1.000000660852943;
					*(input -> dL) = 0.014707;
					break;
				case 30:
					*(input -> n) = 1.000000256114134;
					*(input -> dL) = 0.024964;
					break;
				case 40:
					*(input -> n) = 1+1.496892e-07;
					*(input -> dL) = 3.527309e-02;
					break;
				case 60:
					*(input -> n) = 1.000000064028534;
					*(input -> dL) = 0.04639;
					break;
				case 80:
					*(input -> n) = 1+3.741268e-08;
					*(input -> dL) = 5.221423e-02;
					break;
				default: 
					argp_failure(state, 1, 0, "Input error with energy");
			}
			break;
		case 'u':
			if(atoi(arg) > 0) 
				*(input -> index_count) = atoi(arg);
			else
				argp_failure(state, 1, 0, "Input error with numebr of indexes for table lookup");
			break;
	}
	return 0;
}
