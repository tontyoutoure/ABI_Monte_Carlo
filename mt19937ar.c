#include <stdio.h>
/* Period parameters */  
#define Nmt19937ar 624
#define Mmt19937ar 397
#define Mmt19937arATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_Mmt19937arASK 0x80000000UL /* most significant w-r bits */
#define LOWER_Mmt19937arASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt19937ar[Nmt19937ar]; /* the array for the state vector  */
static int mt19937ari=Nmt19937ar+1; /* mt19937ari==Nmt19937ar+1 means mt19937ar[Nmt19937ar] is not initialized */

/* initializes mt19937ar[Nmt19937ar] with a seed */
void init_genrand(unsigned long s)
{
    mt19937ar[0]= s & 0xffffffffUL;
    for (mt19937ari=1; mt19937ari<Nmt19937ar; mt19937ari++) {
        mt19937ar[mt19937ari] = 
	    (1812433253UL * (mt19937ar[mt19937ari-1] ^ (mt19937ar[mt19937ari-1] >> 30)) + mt19937ari); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, Mmt19937arSBs of the seed affect   */
        /* only Mmt19937arSBs of the array mt19937ar[].                        */
        /* 2002/01/09 modified by Mmt19937arakoto Mmt19937aratsumoto             */
        mt19937ar[mt19937ari] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (Nmt19937ar>key_length ? Nmt19937ar : key_length);
    for (; k; k--) {
        mt19937ar[i] = (mt19937ar[i] ^ ((mt19937ar[i-1] ^ (mt19937ar[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt19937ar[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=Nmt19937ar) { mt19937ar[0] = mt19937ar[Nmt19937ar-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=Nmt19937ar-1; k; k--) {
        mt19937ar[i] = (mt19937ar[i] ^ ((mt19937ar[i-1] ^ (mt19937ar[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt19937ar[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=Nmt19937ar) { mt19937ar[0] = mt19937ar[Nmt19937ar-1]; i=1; }
    }

    mt19937ar[0] = 0x80000000UL; /* Mmt19937arSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, Mmt19937arATRIX_A};
    /* mag01[x] = x * Mmt19937arATRIX_A  for x=0,1 */

    if (mt19937ari >= Nmt19937ar) { /* generate Nmt19937ar words at one time */
        int kk;

        if (mt19937ari == Nmt19937ar+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<Nmt19937ar-Mmt19937ar;kk++) {
            y = (mt19937ar[kk]&UPPER_Mmt19937arASK)|(mt19937ar[kk+1]&LOWER_Mmt19937arASK);
            mt19937ar[kk] = mt19937ar[kk+Mmt19937ar] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<Nmt19937ar-1;kk++) {
            y = (mt19937ar[kk]&UPPER_Mmt19937arASK)|(mt19937ar[kk+1]&LOWER_Mmt19937arASK);
            mt19937ar[kk] = mt19937ar[kk+(Mmt19937ar-Nmt19937ar)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt19937ar[Nmt19937ar-1]&UPPER_Mmt19937arASK)|(mt19937ar[0]&LOWER_Mmt19937arASK);
        mt19937ar[Nmt19937ar-1] = mt19937ar[Mmt19937ar-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mt19937ari = 0;
    }
  
    y = mt19937ar[mt19937ari++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


