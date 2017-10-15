#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include<time.h>

double randomNumber(int ubound, int lbound){
    double s;
    s = ((double)rand()/(RAND_MAX))*(ubound-lbound);
    return s;
}

void copyMatrix(double *a, double *b, int n){
    int i;
    for(i=0; i< n*n;i++){
        b[i]=a[i];
    }
}

double* assignMatVal(double *a, int n, int ubound, int lbound){
    int i;
    for(i=0;i<n*n;i++){
        a[i] = randomNumber(ubound, lbound);
    }
    return a;
}

void calcPerformance(clock_t start, clock_t end,int n){
    printf("Matrix Size = %d\n",n);
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    printf("Time taken in seconds = %f\n",time);
    double perfomance = 2*((double)pow(n,3))/(time*pow(10,9));
    printf("Performance in GFLOPS = %e\n",perfomance);
}

void calcExecutionTime(clock_t start, clock_t end,int n, char s[]){
    printf("Matrix Size = %d\n",n);
    printf("Matrix Config = %s\n", s);
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    printf("Time taken in seconds = %.10f\n",time);
}

void findCorrectness(double *a, double *b,int n){
    double error = 0.0000000;
    int i;
    for(i=0;i<n*n;i++){
        if(error < abs(a[i]-b[i]))
            error = abs(a[i]-b[i]);
    }
    printf("Error = %f\n",error);
    printf("\n");
}

void regCacheMat(double *a, double *b, double *c, int n, int B){
    int i,j,k,i1,j1,k1;
    clock_t start,end;
    start = clock();
    for (i = 0; i < n; i+=B){
        for (j = 0; j < n; j+=B){
            for (k = 0; k < n; k+=B){
     /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+B; i1+=2){
                    for (j1 = j; j1 < j+B; j1+=2){
                        register int t = i1*n+j1; register int tt = t+n;
                        register double c00 = c[t]; register double c01 = c[t+1];  register double c10 = c[tt]; register double c11 = c[tt+1];
                        for (k1 = k; k1 < k+B; k1+=2){
                            register int ta = i1*n+k1; register int tta = ta+n; register int tb = k1*n+j1; register int ttb = tb+n;
                            register double a00 = a[ta]; register double a10 = a[tta]; register double b00 = b[tb]; register double b01 = b[tb+1];

                            c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;

                            a00 = a[ta+1]; a10 = a[tta+1]; b00 = b[ttb]; b01 = b[ttb+1];

                            c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
                        }
                        c[t] = c00;
                        c[t+1] = c01;
                        c[tt] = c10;
                        c[tt+1] = c11;
                    }
                }
            }
        }
    }
    end = clock();
    printf("Alogrithm: dgemm3\n");
    printf("Block size: %d\n",B);
    calcPerformance(start,end,n);

}

void dgemm0(double *a, double *b, double *c, int n){
    clock_t start,end;
    int i,j,k;
    start = clock();
    for (i=0; i<n; i++ ){
        for (j=0; j<n; j++ ){
            for (k=0; k<n; k++ ){
                c[i*n+j]+=a[i*n+k]*b[k*n+j];
            }
        }
    }
    end = clock();
    printf("Algorithm: dgemm0\n");
    calcPerformance(start,end,n);
}

int main(){
    int i,n;
    double *arrA,*arrB,*arrC0, *arrC1, * arrC2;
    int arrN[] = {2,4,8,16,32,64,128,256,512,1024};
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    int len = sizeof(arrN)/sizeof(arrN[0]);
    n=2048;
    for(i=0;i<len;i++){
        int B = arrN[i];
		arrA = (double *)calloc(sizeof(double), n*n);
		arrB = (double *)calloc(sizeof(double), n*n);
		arrC0 = (double *)calloc(sizeof(double), n*n);
        arrC1 = (double *)calloc(sizeof(double), n*n);
        arrC2 = (double *)calloc(sizeof(double), n*n);
        arrA = assignMatVal(arrA,n,ubound,lbound);
        arrB = assignMatVal(arrB,n,ubound,lbound);
        arrC0 = assignMatVal(arrC0,n,ubound,lbound);
        copyMatrix(arrC0,arrC1,n);
        copyMatrix(arrC0,arrC2,n);
		dgemm0(arrA,arrB,arrC1,n);
		regCacheMat(arrA,arrB,arrC0,n,B);
		findCorrectness(arrC0,arrC1,n);
		free(arrA);
		free(arrB);
		free(arrC0);
        free(arrC1);
        free(arrC2);
    }
    return 0;
}
