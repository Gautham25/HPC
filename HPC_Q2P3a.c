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
    printf("Time taken in seconds = %.10f s\n",time);
    double perfomance = 2*((double)pow(n,3))/(time*pow(10,9));
    printf("Performance in GFLOPS = %.10f\n",perfomance);
    printf("\n");
}

void calcExecutionTime(clock_t start, clock_t end,int n, char s[]){
    printf("Matrix Size = %d\n",n);
    printf("Matrix Config = %s\n", s);
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    printf("Time taken in seconds = %.10f\n",time);
    double perfomance = 2*((double)pow(n,3))/(time*pow(10,9));
    printf("Performance in GFLOPS = %.10f\n",perfomance);
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

void matTripleijk(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
            register double r=c[i*n+j];
            for (k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j]=r;
    }
    end = clock();
    calcExecutionTime(start,end,n,"ijk");
}

void matTriplejik(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            register double r = c[i*n+j];
            for (k=0; k<n; k++)
                r += a[i*n+k] * b[k*n+j];
            c[i*n+j] = r;
        }
    }
    end = clock();
    calcExecutionTime(start,end,n,"jik");
}

void matTriplekij(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (k=0; k<n; k++)
        for (i=0; i<n; i++) {
            register double r=a[i*n+k];
            for (j=0; j<n; j++){
                c[i*n+j] += r * b[k*n+j];
            }

        }
    end = clock();
    calcExecutionTime(start,end,n,"kij");
}

void matTripleikj(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (i=0; i<n; i++)
        for (k=0; k<n; k++) {
            register double r= a[i*n+k];
            for (j=0; j<n; j++)
                c[i*n+j] += r * b[k*n+j];
        }
    end = clock();
    calcExecutionTime(start,end,n,"ikj");
}

void matTriplejki(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (j=0; j<n; j++)
        for (k=0; k<n; k++) {
            register double r=b[k*n+j];
            for (i=0; i<n; i++)
                c[i*n+j] += a[i*n+k] * r;
    }
    end = clock();
    calcExecutionTime(start,end,n,"jki");
}

void matTriplekji(double *a, double *b, double *c, int n){
    int i,j,k;
    double sum;
    clock_t start,end;
    start = clock();
    for (k=0; k<n; k++)
        for (j=0; j<n; j++) {
            register double r=b[k*n+j];
            for (i=0; i<n; i++)
                c[i*n+j] += a[i*n+k] * r;
        }
    end = clock();
    calcExecutionTime(start,end,n,"kji");
}

int main(){
    int i,n;
    double *arrA,*arrB,*arrC0, *arrC1, * arrC2;
    //double arrN[] = {64,128,256,512,1024,2048};
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    n=2048;

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

    //Triple loop matrix multiplications
    printf("\nTRIPLE LOOP MATRIX MULTIPLICATIONS\n");
    matTripleijk(arrA,arrB,arrC0,n);
    printf("\n");
    matTriplejik(arrA,arrB,arrC1,n);
    printf("\nMatrix Configs: ijk & jik\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matTriplekij(arrA,arrB,arrC0,n);
    printf("\nMatrix Configs: jik & kij\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC1);

    arrC1 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC1,n);
    matTripleikj(arrA,arrB,arrC1,n);
    printf("\nMatrix Configs: kij & ikj\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matTriplejki(arrA,arrB,arrC0,n);
    printf("\nMatrix Configs: ikj & jki\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC1);

    arrC1 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC1,n);
    matTriplekji(arrA,arrB,arrC1,n);
    printf("\nMatrix Configs: jki & kji\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matTripleikj(arrA,arrB,arrC0,n);
    printf("\nMatrix Configs: kji & ikj\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);
    free(arrC1);
    free(arrA);
    free(arrB);

	return 0;
}
