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

void matBlockijk(double *a, double *b, double *c, int n,int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (i = 0; i < n; i+=B)
        for (j = 0; j < n; j+=B)
            for (k = 0; k < n; k+=B)
    /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+B; i1++)
                    for (j1 = j; j1 < j+B; j1++) {
                        register double r=c[i1*n+j1];
                        for (k1 = k; k1 < k+B; k1++)
                            r += a[i1*n + k1]*b[k1*n + j1];
                        c[i1*n+j1]=r;
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block ijk");
}

void matBlockjik(double *a, double *b, double *c, int n,int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (j = 0; j < n; j+=B)
        for (i = 0; i < n; i+=B)
            for (k = 0; k < n; k+=B)
    /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+B; j1++)
                    for (i1 = i; i1 < i+B; i1++) {
                        register double r=c[i1*n+j1];
                        for (k1 = k; k1 < k+B; k1++)
                            r += a[i1*n + k1]*b[k1*n + j1];
                        c[i1*n+j1]=r;
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block jik");
}

void matBlockkij(double *a, double *b, double *c, int n,int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (k = 0; k < n; k+=B)
        for (i = 0; i < n; i+=B)
            for (j = 0; j < n; j+=B)
    /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+B; k1++)
                    for (i1 = i; i1 < i+B; i1++) {
                        register double r=a[i1*n+k1];
                        for (j1 = j; j1 < j+B; j1++)
                            c[i1*n+j1] += r*b[k1*n + j1];
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block kij");
}

void matBlockikj(double *a, double *b, double *c, int n,int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (i = 0; i < n; i+=B)
        for (k = 0; k < n; k+=B)
            for (j = 0; j < n; j+=B)
    /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+B; i1++)
                    for (k1 = k; k1 < k+B; k1++) {
                        register double r=a[i1*n+k1];
                        for (j1 = j; j1 < j+B; j1++)
                            c[i1*n+j1] += r*b[k1*n + j1];
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block ikj");
}

void matBlockjki(double *a, double *b, double *c, int n,int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (j = 0; j < n; j+=B)
        for (k = 0; k < n; k+=B)
            for (i = 0; i < n; i+=B)
    /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+B; j1++)
                    for (k1 = k; k1 < k+B; k1++) {
                        register double r=b[k1*n+j1];
                        for (i1 = i; i1 < i+B; i1++)
                            c[i1*n+j1] += a[i1*n+k1]*r;
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block jki");
}

void matBlockkji(double *a, double *b, double *c, int n, int B){
    int i,j,k,i1,k1,j1;
    double sum;
    clock_t start,end;
    start = clock();
    for (k = 0; k < n; k+=B)
        for (j = 0; j < n; j+=B)
            for (i = 0; i < n; i+=B)
    /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+B; k1++)
                    for (j1 = j; j1 < j+B; j1++) {
                        register double r=b[k1*n+j1];
                        for (i1 = i; i1 < i+B; i1++)
                            c[i1*n+j1] += a[i1*n+k1]*r;
                    }
    end = clock();
    calcExecutionTime(start,end,n,"Block kji");
}

int main(){
    int i,n;
    double *arrA,*arrB,*arrC0, *arrC1, * arrC2;
    //double arrN[] = {64,128,256,512,1024,2048};
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    n=2048;
    int B[]={16,32,64,128,256};

    //Blocked matrix multiplications
    for(i=0;i<(sizeof(B)/sizeof(B[0]));i++){
        arrA = (double *)calloc(sizeof(double), n*n);
        arrB = (double *)calloc(sizeof(double), n*n);
        arrC0 = (double *)calloc(sizeof(double), n*n);
        arrC1 = (double *)calloc(sizeof(double), n*n);
        arrC2 = (double *)calloc(sizeof(double), n*n);
        arrA = assignMatVal(arrA, n, ubound, lbound);
        arrB = assignMatVal(arrB, n, ubound, lbound);
        arrC0 = assignMatVal(arrC0, n, ubound, lbound);
        copyMatrix(arrC0,arrC1,n);
        copyMatrix(arrC0,arrC2,n);
        printf("\nBLOCKED MATRIX MULTIPLICATIONS\n");
        printf("\nBlock Size= %d\n",B[i]);
        matBlockijk(arrA,arrB,arrC0,n,8);
        printf("\n");
        matBlockjik(arrA,arrB,arrC1,n,B[i]);
        printf("\nMatrix Configs: Block ijk & jik\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrC0);

        arrC0 = (double *)calloc(sizeof(double), n*n);
        copyMatrix(arrC2,arrC0,n);
        matBlockkij(arrA,arrB,arrC0,n,B[i]);
        printf("\nMatrix Configs: Block jik & kij\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrC1);

        arrC1 = (double *)calloc(sizeof(double), n*n);
        copyMatrix(arrC2,arrC1,n);
        matBlockikj(arrA,arrB,arrC1,n,B[i]);
        printf("\nMatrix Configs: Block kij & ikj\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrC0);

        arrC0 = (double *)calloc(sizeof(double), n*n);
        copyMatrix(arrC2,arrC0,n);
        matBlockjki(arrA,arrB,arrC0,n,B[i]);
        printf("\nMatrix Configs: Block ikj & jki\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrC1);

        arrC1 = (double *)calloc(sizeof(double), n*n);
        copyMatrix(arrC2,arrC1,n);
        matBlockkji(arrA,arrB,arrC1,n,B[i]);
        printf("\nMatrix Configs: Block jki & kji\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrC0);

        arrC0 = (double *)calloc(sizeof(double), n*n);
        copyMatrix(arrC2,arrC0,n);
        matBlockikj(arrA,arrB,arrC0,n,B[i]);
        printf("\nMatrix Configs: Block kji & ikj\n");
        findCorrectness(arrC0,arrC1,n);
        printf("\n");
        free(arrA);
        free(arrB);
        free(arrC0);
        free(arrC1);
        free(arrC2);
    }
	return 0;
}
