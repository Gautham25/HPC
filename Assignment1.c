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
    printf("Time taken in seconds = %.10f\n",time);
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
    printf("Correctness = %f\n",error);
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

void dgemm1(double *a, double *b, double *c, int n){
    clock_t start,end;
    int i,j,k;
    start = clock();
    for (i=0; i<n; i++ ){
        for (j=0; j<n; j++ ){
            register double r = 0; // = c[i*n+j];
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
            register double c00=0.0, c01=0.0, c10=0.0, c11=0.0;
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

void regCacheMat(double *a, double *b, double *c, int n, int B){
    int i,j,k,i1,j1,k1;
    for (i = 0; i < n; i+=B){
        for (j = 0; j < n; j+=B){
            for (k = 0; k < n; k+=B){
     /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+B; i1++){
                    for (j1 = j; j1 < j+B; j1++){
                        register int t = i1*n+j1; register int tt = t+n;
                        register double c00 = c[t]; register double c01 = c[t+1];  register double c10 = c[tt]; register double c11 = c[tt+1];
                        for (k1 = k; k1 < k+B; k1++){
                            register int ta = i1*n+k1; register int tta = ta+n; register int tb = k1*n+j1; register int ttb = tb+n;
                            register double a00 = a[ta]; register double a01 = a[ta+1]; register double a10 = a[tta]; register double a11 = a[tta+1];
                            register double b00 = b[tb]; register double b01 = b[tb+1]; register double b10 = b[ttb]; register double b11 = b[ttb+1];
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
            }
        }
    }


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
    n=64;
    B=2;
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

    //Blocked matrix multiplications
    arrC0 = (double *)calloc(sizeof(double), n*n);
    arrC1 = (double *)calloc(sizeof(double), n*n);
    arrC2 = (double *)calloc(sizeof(double), n*n);
    arrC0 = assignMatVal(arrC0, n, ubound, lbound);
    copyMatrix(arrC0,arrC1,n);
    copyMatrix(arrC0,arrC2,n);
    printf("\nBLOCKED MATRIX MULTIPLICATIONS\n");
    matBlockijk(arrA,arrB,arrC0,n,B);
    printf("\n");
    matBlockjik(arrA,arrB,arrC1,n,B);
    printf("\nMatrix Configs: Block ijk & jik\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matBlockkij(arrA,arrB,arrC0,n,B);
    printf("\nMatrix Configs: Block jik & kij\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC1);

    arrC1 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC1,n);
    matBlockikj(arrA,arrB,arrC1,n,B);
    printf("\nMatrix Configs: Block kij & ikj\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matBlockjki(arrA,arrB,arrC0,n,B);
    printf("\nMatrix Configs: Block ikj & jki\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC1);

    arrC1 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC1,n);
    matBlockkji(arrA,arrB,arrC1,n,B);
    printf("\nMatrix Configs: Block jki & kji\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrC0);

    arrC0 = (double *)calloc(sizeof(double), n*n);
    copyMatrix(arrC2,arrC0,n);
    matBlockikj(arrA,arrB,arrC0,n,B);
    printf("\nMatrix Configs: Block kji & ikj\n");
    findCorrectness(arrC0,arrC1,n);
    printf("\n");
    free(arrA);
    free(arrB);
    free(arrC0);
    free(arrC1);
    free(arrC2);
	return 0;
}
