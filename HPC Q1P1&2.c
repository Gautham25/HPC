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
    printf("Time taken = %.10f\n",time);
    double perfomance = 2*((double)pow(n,3))/(time*pow(10,9));
    printf("Performance in GFLOPS = %.10f\n",perfomance);
    printf("\n");
}

void calcExecutionTime(clock_t start, clock_t end,int n, char s[]){
    printf("Matrix Size = %d\n",n);
    printf("Matrix Config = %s\n", s);
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    printf("Time taken = %.10f\n",time);
}

void findCorrectness(double *a, double *b,int n){
    double error = 0.0000000;
    int i;
    for(i=0;i<n*n;i++){
        if(error < abs(a[i]-b[i]))
            error = abs(a[i]-b[i]);
    }
    printf("Error = %f\n",error);
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
    printf("\nAlgorithm: dgemm0\n");
    calcPerformance(start,end,n);
}

void dgemm1(double *a, double *b, double *c, int n){
    clock_t start,end;
    int i,j,k;
    start = clock();
    for (i=0; i<n; i++ ){
        for (j=0; j<n; j++ ){
            register double r = c[i*n+j];
            for (k=0; k<n; k++ ){
                r += a[i*n+k]*b[k*n+j];
            }
            c[i*n+j] = r;
        }
    }
    end = clock();
    printf("Algorithm: dgemm1\n");
    calcPerformance(start,end,n);
}

void dgemm2(double *a, double *b, double *c, int n){
    clock_t start,end;
    int i,j,k;
    start = clock();
    for (i=0; i<n; i+=2 ){
        for (j=0; j<n; j+=2 ){
            register int t = i*n+j;register int tt = t+n;
            register double c00=c[t], c01=c[t+1], c10=c[tt], c11=c[tt+1];
            for (k=0; k<n; k+=2 ){
                register int ta = i*n+k,tta = ta+n, tb = k*n+j, ttb = tb+n;
                register double a00 = a[ta], a01 = a[ta+1], a10 = a[tta], a11 = a[tta+1];
                register double b00 = b[tb], b01 = b[tb+1], b10 = b[ttb], b11 = b[ttb+1];
                c00 += a00*b00 + a01*b10;
                c01 += a00*b01 + a01*b11;
                c10 += a10*b00 + a11*b10;
                c11 += a10*b01 + a11*b11;
            }
            c[t] = c00;
            c[t+1] = c01;
            c[tt] = c10;
            c[tt+1] = c11;
        }
    }
    end = clock();
    printf("Algorithm: dgemm2\n");
    calcPerformance(start,end,n);
}

int main(){
    int i,B,n;
    double *arrA,*arrB,*arrC0, *arrC1, * arrC2;
    double arrN[] = {64,128,256,512,1024,2048};
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    int len = (sizeof(arrN)/sizeof(arrN[0]));
    for(i=0;i<len;i++){
		n = arrN[i];
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
		dgemm0(arrA,arrB,arrC0,n);
        dgemm1(arrA,arrB,arrC1,n);
		dgemm2(arrA,arrB,arrC2,n);
        printf("\nMatrix Config = dgemm0 & dgemm1\n");
        findCorrectness(arrC0,arrC1,arrN[i]);
        printf("\nMatrix Config = dgemm1 & dgemm2\n");
        findCorrectness(arrC1,arrC2,arrN[i]);
		free(arrA);
		free(arrB);
		free(arrC0);
        free(arrC1);
        free(arrC2);
    }
    return 0;
}
