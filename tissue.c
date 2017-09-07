#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<unistd.h>
#include<string.h>
#include<pthread.h>
#include"mt19937ar.c"
#include"tissue.h"


double a,b,c,R,DX,DC,ce=0;
double *Cache;
long i,total=0;
int TryTime;
char StrPos[FILENAME_MAX],criteria_;
pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_t tid[NoT];


int main(){
	long CacheSize[2]={0}, j;
	time_t timer;
	FILE *FPPos;
	
	a=0.01, b=0.01, c=0.00022, R=0.0001, DX=R/1000;
	sprintf(StrPos,"%f_%f_%f.tissue",a,b,c);

	init_genrand(time(0));
	DC=MaxCacheSize*pow(R,3)/a/b/0.4;
	a+=4*R;
	b+=4*R;
	R+=4*DX;
	Cache=malloc(sizeof(double)*MaxCacheSize*3);
	time(&timer);
	
	FPPos=fopen(StrPos,"w");
	if(FPPos==NULL){
		printf("Input error, please try again\n");
		return -1;
	}//make sure the file exist
	fclose(FPPos);	
	

	while(ce<c){		
		for(i=0;i<CacheSize[1];i++)
			Cache[i*3+0]=Cache[(i+CacheSize[0])*3+0],
			Cache[i*3+1]=Cache[(i+CacheSize[0])*3+1],
			Cache[i*3+2]=Cache[(i+CacheSize[0])*3+2];
		CacheSize[0]=CacheSize[1],CacheSize[1]=0;
		ce+=DC;
		if(ce>c)
			ce=c;
		TryTime=0;
		criteria_=1;
		for(j=0;j<NoT;j++){
			pthread_create(tid+j,NULL,fallingstar,tid+j);
		}
		for(;criteria_;){
			sleep(2);
		}
		for(j=0;j<NoT;j++){
			pthread_join(tid[j],NULL);
		}
		CacheSize[1]=i-CacheSize[0];
		qsort(Cache,i,sizeof(double)*3,compare);
		FPPos=fopen(StrPos,"a");
				
		fwrite(Cache,sizeof(double),3*CacheSize[0],FPPos);
		fclose(FPPos);//sort the balls by the z coordinate of the ball and output
		
		
		
			
	}
	
	FPPos=fopen(StrPos,"a");
	fwrite(Cache+3*CacheSize[0],sizeof(double),3*CacheSize[1],FPPos);
	fclose(FPPos);
	return 0;
}

void *fallingstar(void *arg){
	double x[3],dist[MaxCacheSize],SC[6],CC[6],dtheta,dphi,middle,rot[9],rot1[9],rot2[9],arot[9],dx=DX,r=R;
	long MaxI,Ind[60],ib,j,k,l,m;
	char criteria;
	int mode;
	
	while(criteria_){
			ib=i;
			x[0]=genrand_res53()*(a-2*r)+r;
			x[1]=genrand_res53()*(b-2*r)+r;
			x[2]=ce+r;
			mode=1;
			
			
			/*	mode 1 is free falling mode
				mode 2 is sliding down a ball
				mode 3 is sliding down two balls
				mode 4 is sliding down a ball and x=0
				mode 5 is sliding down a ball and x=a
				mode 6 is sliding down a ball and y=0
				mode 7 is sliding down a ball and y=b
				mode 8 is sliding down x=0
				mode 9 is sliding down x=a
				mode 10 is sliding down y=0
				mode 11 is sliding down y=b
				mode 0 is stoping mode*/
			while(mode){
				switch(mode){
					case 1:

						for(j=0;j<ib;j++){
							dist[j]=DIST2(j);
						}
						for(j=0,MaxI=-1;j<ib;j++)
							if(dist[j]<4*r*r&&Cache[j*3+2]<x[2]){
								if(MaxI==-1)
									MaxI=j;
								else if(Cache[j*3+2]+sqrt(4*r*r-dist[j])>Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]))
									MaxI=j;									
							}
						if(MaxI==-1)
							x[2]=r,mode=0;
						else
							x[2]=Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]),j=MaxI,mode=2;
					
						break;
						
					case 2:
						CC[0*3+0]=x[0]-Cache[j*3+0],
						CC[0*3+1]=x[1]-Cache[j*3+1],
						CC[0*3+2]=x[2]-Cache[j*3+2];
						car2sph(CC,SC);
						for(k=0,l=0;k<ib;k++)
							if(pow((Cache[k*3+0]-Cache[j*3+0]),2)+pow((Cache[k*3+1]-Cache[j*3+1]),2)+pow((Cache[k*3+2]-Cache[j*3+2]),2)<25*r*r&&k!=j)
								Ind[l++]=k;
					
						dtheta=dx/2/r;
						while(SC[0*3+1]<pi/2){

							SC[0*3+1]+=dtheta;
							sph2car(SC,CC);
							x[0]=Cache[j*3+0]+CC[0*3+0],
							x[1]=Cache[j*3+1]+CC[0*3+1],
							x[2]=Cache[j*3+2]+CC[0*3+2];
						
							for(k=0;k<l;k++){
								if(DIST3(Ind[k])<4*r*r){
									mode=3,k=Ind[k];

									CC[1*3+0]=Cache[k*3+0]-Cache[j*3+0],
									CC[1*3+1]=Cache[k*3+1]-Cache[j*3+1],
									CC[1*3+2]=Cache[k*3+2]-Cache[j*3+2];
									car2sph(CC+3,SC+3);
									middle=CC[1*3+0]*cos(SC[0*3+2])+CC[1*3+1]*sin(SC[0*3+2]);									
									SC[0*3+1]=acos((CC[1*3+2]*pow(SC[1*3+0],2)/SC[0*3+0]+sqrt(4*middle*middle*(pow(CC[1*3+2],2)+middle*middle)-middle*middle*pow(SC[1*3+0],4)/pow(SC[0*3+0],2)))/(2*(pow(CC[1*3+2],2)+middle*middle)));
									sph2car(SC,CC);
									x[0]=Cache[j*3+0]+CC[0*3+0],
									x[1]=Cache[j*3+1]+CC[0*3+1],
									x[2]=Cache[j*3+2]+CC[0*3+2];

									break;
								}
							}
								
							if(mode==3){
								
								break;
								}
							if(x[0]<r){							
								SC[0*3+1]=asin((r-Cache[j*3+0])/2/r/cos(SC[0*3+2]));
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+0],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+2];
								mode=4;
								break;
							}
							if(x[0]>a-r){
								SC[0*3+1]=asin((a-r-Cache[j*3+0])/2/r/cos(SC[0*3+2]));
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+0],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+2];
								mode=5;
								break;
							}
							if(x[1]<r){							
								SC[0*3+1]=asin((r-Cache[j*3+1])/2/r/sin(SC[0*3+2]));
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+0],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+2];
								mode=6;
								break;
							}
							if(x[1]>b-r){							
								SC[0*3+1]=asin((b-r-Cache[j*3+1])/2/r/sin(SC[0*3+2]));
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+0],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+2];
								mode=7;
								break;
							}
						}
						if(mode==2)
							mode=1;		
											
						break;

					case 3:
						if(Cache[j*3+2]>Cache[k*3+2]){
							l=j,j=k,k=l;}
						CC[1*3+0]=Cache[k*3+0]-Cache[j*3+0],
						CC[1*3+1]=Cache[k*3+1]-Cache[j*3+1],
						CC[1*3+2]=Cache[k*3+2]-Cache[j*3+2];
						car2sph(CC+3,SC+3);
						for(l=0;l<9;l++)
							rot1[l]=0,rot2[l]=0;
						rot1[0]=cos(SC[1*3+2]),rot1[1]=-sin(SC[1*3+2]),rot1[3]=sin(SC[1*3+2]),rot1[4]=cos(SC[1*3+2]),rot1[8]=1;
						rot2[0]=cos(SC[1*3+1]),rot2[2]=sin(SC[1*3+1]),rot2[4]=1,rot2[6]=-sin(SC[1*3+1]),rot2[8]=cos(SC[1*3+1]);
						mm3(rot1,rot2,rot);
						arot[0]=rot[0], arot[1]=rot[3], arot[2]=rot[6], arot[3]=rot[1], arot[4]=rot[4], arot[5]=rot[7], arot[6]=rot[2], arot[7]=rot[5], arot[8]=rot[8];
						CC[1*3+0]=x[0]-Cache[j*3+0],
						CC[1*3+1]=x[1]-Cache[j*3+1],
						CC[1*3+2]=x[2]-Cache[j*3+2];
						mm31(arot,CC+3,CC);
						car2sph(CC,SC);
						dphi=dx/sin(SC[0*3+1])/SC[0*3+0];
						for(l=0,m=0;l<i;l++){
							if(((Cache[k*3+0]-Cache[l*3+0])*(Cache[k*3+0]-Cache[l*3+0])+(Cache[k*3+1]-Cache[l*3+1])*(Cache[k*3+1]-Cache[l*3+1])+(Cache[k*3+2]-Cache[l*3+2])*(Cache[k*3+2]-Cache[l*3+2])<25*r*r||(Cache[j*3+0]-Cache[l*3+0])*(Cache[j*3+0]-Cache[l*3+0])+(Cache[j*3+1]-Cache[l*3+1])*(Cache[j*3+1]-Cache[l*3+1])+(Cache[j*3+2]-Cache[l*3+2])*(Cache[j*3+2]-Cache[l*3+2])<25*r*r)&&l!=j&&l!=k)
							Ind[m++]=l;
						}
						if(SC[0*3+2]>pi)
							while(SC[0*3+2]<3*pi/2){
								SC[0*3+2]+=dphi;
								sph2car(SC,CC);
								mm31(rot,CC,CC+3);
								x[0]=Cache[j*3+0]+CC[1*3+0],
								x[1]=Cache[j*3+1]+CC[1*3+1],
								x[2]=Cache[j*3+2]+CC[1*3+2];
								if(x[0]<r||x[0]>a-r||x[1]<r||x[1]>b-r||x[2]<r){
									mode=0;break;
								}
								for(l=0;l<m;l++)								
									if((Cache[Ind[l]*3+0]-x[0])*(Cache[Ind[l]*3+0]-x[0])+(Cache[Ind[l]*3+1]-x[1])*(Cache[Ind[l]*3+1]-x[1])+(Cache[Ind[l]*3+2]-x[2])*(Cache[Ind[l]*3+2]-x[2])<4*r*r){
										mode=0;break;
									}
								if(mode==0)
									break;								
							}
						else
							while(SC[0*3+2]>pi/2){
								SC[0*3+2]-=dphi;
								sph2car(SC,CC);
								mm31(rot,CC,CC+3);
								x[0]=Cache[j*3+0]+CC[1*3+0],
								x[1]=Cache[j*3+1]+CC[1*3+1],
								x[2]=Cache[j*3+2]+CC[1*3+2];
								if(x[0]<r||x[0]>a-r||x[1]<r||x[1]>b-r||x[2]<r){
									mode=0;break;
								}
								for(l=0;l<m;l++)								
									if((Cache[Ind[l]*3+0]-x[0])*(Cache[Ind[l]*3+0]-x[0])+(Cache[Ind[l]*3+1]-x[1])*(Cache[Ind[l]*3+1]-x[1])+(Cache[Ind[l]*3+2]-x[2])*(Cache[Ind[l]*3+2]-x[2])<4*r*r){
										mode=0;break;
									}
								if(mode==0)
									break;
							}
						if(mode==3)
							mode=2;	
						break;
					case 4:
					
						CC[0*3+0]=x[2]-Cache[j*3+2],
						CC[0*3+1]=x[1]-Cache[j*3+1],
						CC[0*3+2]=Cache[j*3+0]-x[0];
						car2sph(CC,SC);
						dphi=dx/SC[0*3+0]/sin(SC[0*3+1]);
						for(k=0,l=0;k<i;k++)
							if((Cache[k*3+0]-Cache[j*3+0])*(Cache[k*3+0]-Cache[j*3+0])+(Cache[k*3+1]-Cache[j*3+1])*(Cache[k*3+1]-Cache[j*3+1])+(Cache[k*3+2]-Cache[j*3+2])*(Cache[k*3+2]-Cache[j*3+2])<25*r*r&&k!=j)
								Ind[l++]=k;
						if(SC[0*3+2]>pi){
							while(SC[0*3+2]>pi*3/2){
								SC[0*3+2]-=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]-CC[0*3+2],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[1]>b-r||x[1]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}							
						}
						else{
							while(SC[0*3+2]<pi/2){
								SC[0*3+2]+=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]-CC[0*3+2],
								x[1]=Cache[j*3+1]+CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[1]>b-r||x[1]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						if(mode==4)
								mode=8;
						break;
					case 5:
					
						CC[0*3+0]=x[2]-Cache[j*3+2],
						CC[0*3+1]=Cache[j*3+1]-x[1],
						CC[0*3+2]=x[0]-Cache[j*3+0];
						car2sph(CC,SC);
						dphi=dx/SC[0*3+0]/sin(SC[0*3+1]);
						for(k=0,l=0;k<i;k++)
							if((Cache[k*3+0]-Cache[j*3+0])*(Cache[k*3+0]-Cache[j*3+0])+(Cache[k*3+1]-Cache[j*3+1])*(Cache[k*3+1]-Cache[j*3+1])+(Cache[k*3+2]-Cache[j*3+2])*(Cache[k*3+2]-Cache[j*3+2])<25*r*r&&k!=j)
								Ind[l++]=k;
						if(SC[0*3+2]>pi){
							while(SC[0*3+2]>pi*3/2){
								SC[0*3+2]-=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+2],
								x[1]=Cache[j*3+1]-CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[1]>b-r||x[1]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						else{
							while(SC[0*3+2]<pi/2){
								SC[0*3+2]+=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+2],
								x[1]=Cache[j*3+1]-CC[0*3+1],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[1]>b-r||x[1]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						if(mode==5)
								mode=9;
						break;
					case 6:
				
						CC[0*3+0]=x[2]-Cache[j*3+2],
						CC[0*3+1]=Cache[j*3+0]-x[0],
						CC[0*3+2]=Cache[j*3+1]-x[1];
						car2sph(CC,SC);
						dphi=dx/SC[0*3+0]/sin(SC[0*3+1]);
						for(k=0,l=0;k<i;k++)
							if((Cache[k*3+0]-Cache[j*3+0])*(Cache[k*3+0]-Cache[j*3+0])+(Cache[k*3+1]-Cache[j*3+1])*(Cache[k*3+1]-Cache[j*3+1])+(Cache[k*3+2]-Cache[j*3+2])*(Cache[k*3+2]-Cache[j*3+2])<25*r*r&&k!=j)
								Ind[l++]=k;
						if(SC[0*3+2]>pi){
							while(SC[0*3+2]>pi*3/2){
								SC[0*3+2]-=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]-CC[0*3+1],
								x[1]=Cache[j*3+1]-CC[0*3+2],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[0]>a-r||x[0]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						else{
							while(SC[0*3+2]<pi/2){
								SC[0*3+2]+=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]-CC[0*3+1],
								x[1]=Cache[j*3+1]-CC[0*3+2],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[0]>a-r||x[0]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						if(mode==6)
								mode=10;
						break;
					case 7:
					
						CC[0*3+0]=x[2]-Cache[j*3+2],
						CC[0*3+1]=x[0]-Cache[j*3+0],
						CC[0*3+2]=x[1]-Cache[j*3+1];
						car2sph(CC,SC);
						dphi=dx/SC[0*3+0]/sin(SC[0*3+1]);
						for(k=0,l=0;k<i;k++)
							if((Cache[k*3+0]-Cache[j*3+0])*(Cache[k*3+0]-Cache[j*3+0])+(Cache[k*3+1]-Cache[j*3+1])*(Cache[k*3+1]-Cache[j*3+1])+(Cache[k*3+2]-Cache[j*3+2])*(Cache[k*3+2]-Cache[j*3+2])<25*r*r&&k!=j)
								Ind[l++]=k;
						if(SC[0*3+2]>pi){
							while(SC[0*3+2]>pi*3/2){
								SC[0*3+2]-=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+1],
								x[1]=Cache[j*3+1]+CC[0*3+2],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[0]>a-r||x[0]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						else{
							while(SC[0*3+2]<pi/2){
								SC[0*3+2]+=dphi;
								sph2car(SC,CC);
								x[0]=Cache[j*3+0]+CC[0*3+1],
								x[1]=Cache[j*3+1]+CC[0*3+2],
								x[2]=Cache[j*3+2]+CC[0*3+0];
								if(x[0]>a-r||x[0]<r){
									mode=0;
									break;
								}
								for(k=0;k<l;k++){
									if((Cache[Ind[k]*3+0]-x[0])*(Cache[Ind[k]*3+0]-x[0])+(Cache[Ind[k]*3+1]-x[1])*(Cache[Ind[k]*3+1]-x[1])+(Cache[Ind[k]*3+2]-x[2])*(Cache[Ind[k]*3+2]-x[2])<4*r*r){
										mode=0;
										break;
									}
								}
								if(mode==0)
									break;
							}
						}
						if(mode==7)
								mode=11;
						break;
					case 8:
					
						for(j=0;j<ib;j++){
							dist[j]=(Cache[j*3+0]-x[0])*(Cache[j*3+0]-x[0])+(Cache[j*3+1]-x[1])*(Cache[j*3+1]-x[1]);
						}
						for(j=0,MaxI=-1;j<ib;j++)
							if(dist[j]<4*r*r&&Cache[j*3+2]<x[2]){
								if(MaxI==-1)
									MaxI=j;
								else if(Cache[j*3+2]+sqrt(4*r*r-dist[j])>Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]))
									MaxI=j;									
							}
						if(MaxI==-1){
							mode=0,x[2]=r;
						}
						else{
							x[2]=Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]),j=MaxI,mode=4;							
						}
						break;
					case 9:
					
						for(j=0;j<ib;j++){
							dist[j]=(Cache[j*3+0]-x[0])*(Cache[j*3+0]-x[0])+(Cache[j*3+1]-x[1])*(Cache[j*3+1]-x[1]);
						}
						for(j=0,MaxI=-1;j<ib;j++)
							if(dist[j]<4*r*r&&Cache[j*3+2]<x[2]){
								if(MaxI==-1)
									MaxI=j;
								else if(Cache[j*3+2]+sqrt(4*r*r-dist[j])>Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]))
									MaxI=j;									
							}
						if(MaxI==-1){
							mode=0,x[2]=r;
						}
						else{
							x[2]=Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]),j=MaxI,mode=5;							
						}
						break;
					case 10:
				
						for(j=0;j<ib;j++){
							dist[j]=(Cache[j*3+0]-x[0])*(Cache[j*3+0]-x[0])+(Cache[j*3+1]-x[1])*(Cache[j*3+1]-x[1]);
						}
						for(j=0,MaxI=-1;j<ib;j++)
							if(dist[j]<4*r*r&&Cache[j*3+2]<x[2]){
								if(MaxI==-1)
									MaxI=j;
								else if(Cache[j*3+2]+sqrt(4*r*r-dist[j])>Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]))
									MaxI=j;									
							}
						if(MaxI==-1){
							mode=0,x[2]=r;
						}
						else{
							x[2]=Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]),j=MaxI,mode=6;							
						}
						break;
					case 11:
					
						for(j=0;j<ib;j++){
							dist[j]=(Cache[j*3+0]-x[0])*(Cache[j*3+0]-x[0])+(Cache[j*3+1]-x[1])*(Cache[j*3+1]-x[1]);
						}
						for(j=0,MaxI=-1;j<ib;j++)
							if(dist[j]<4*r*r&&Cache[j*3+2]<x[2]){
								if(MaxI==-1)
									MaxI=j;
								else if(Cache[j*3+2]+sqrt(4*r*r-dist[j])>Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]))
									MaxI=j;									
							}
						if(MaxI==-1){
							mode=0,x[2]=r;
						}
						else{
							x[2]=Cache[MaxI*3+2]+sqrt(4*r*r-dist[MaxI]),j=MaxI,mode=7;							
						}
						break;
					
				}
				
				
				
			}
			if(x[2]>ce-r){			
				TryTime++;
				if(TryTime>MaxTryTime)
					criteria_=0;	
			}				
			else{
			pthread_mutex_lock(&mutex);
			
			for(criteria=1,j=ib;j<i;j++)
				criteria*=(DIST3(j)>=5*r*r);
			if(criteria){
				Cache[i*3+0]=x[0],
				Cache[i*3+1]=x[1],
				Cache[i*3+2]=x[2];
				TryTime=0,i++;
				total++;
				printf("%ld\n",total);
			}
			
			pthread_mutex_unlock(&mutex);			
		}
		
	}
	return NULL;
}

int compare(const void *p1, const void *p2)
{
	return (*((double *)p1+2)-*((double *)p2+2)>0)-(*((double *)p1+2)-*((double *)p2+2)<0);
}
