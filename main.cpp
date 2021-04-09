#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void RESHI(double *a, double *x, double *b,int N);
int Max(double *p,int nn, int M);
int Min(double *p,int w,int nn, int M);
void makeAandB(double *p,double *A, double *B, double t, double h,double a,int j,int M);
void PRINT(double *p,double t,double h, double s,int x, int y,int z,int c,int u, int g,int nn, int M, int T, double eps);
void makeDx(double *Dx, double *p, double h,int nn, int M);
void makeDt(double *Dt,double *p,double t,int nn, int M);



double diff(double a, double b, double p[],int nn, int M, int T, double eps,double *Dt)
{
    int i,j,k,*c,L,MAX,MIN,num=0,minDx,maxDx,minDt,maxDt;
    double x=0,*X,*A,*B,s,*Dx;
    double t, n1, M1, h;
    n1=nn;
    M1=M;
    t=T/(n1-1);
    h=1/(M1-1);
    for(int i=M*(nn-1); i<=M*nn-1; i++){p[i]=1;}

    for(int j=nn-2; j>=0; j--)
    {
        B=(double*)malloc((M-2)*sizeof(double));
        if(B==NULL){printf("2\n");return -1;}
        X=(double*)malloc((M-2)*sizeof(double));
        if(X==NULL){printf("3\n");return -1;}
        A=(double*)malloc((M-2)*(M-2)*sizeof(double));
        if(A==NULL){printf("5\n");return -1;}
        Dx=(double*)malloc(M*nn*sizeof(double));
        if(Dx==NULL){printf("6\n");return -1;}

        c=(int*)malloc((M-2)*sizeof(int));
        if(c==NULL){printf("7\n");return -1;}

        for(i=0;i<(M-2);i++){c[i]=i+1;}
        for(i=0;i<(M-2)*(M-2);i++){A[i]=0;}
        for(i=0;i<M*nn;i++){Dx[i]=100;}
        for(i=0;i<M*nn;i++){Dt[i]=100;}

        makeAandB(p,A,B,t,h,a,j,M);
        RESHI(A,X,B,M-2);
        L=0;
        for(i=1;i<=M-2;i++)
        {
            p[j*M+i]=X[i-1];
            if(fabs(p[j*M+i]-p[(j+1)*M+i])<eps){L++;}
        }
        p[j*M]=1;
        p[j*M+M-1] = 1;
        if(L==M-2){s=(nn-j)*t;num=j;j=-1;}
        makeDx(Dx,p,h,nn,M);
//        makeDt(Dt,p,t,nn,M);
    }
    /*for(int i=0; i<=nn*M-1; i++)
    {
        if (i%M==0) printf("\n");
        printf("%.8f   ", p[i]);
    }
    printf("\n");
    printf("\n");
    for(int i=0; i<=nn*M-1; i++)
    {
        if (i%M==0) printf("\n");
        printf("%.8f   ", Dx[i]);
    }
    printf("\n");
    printf("\n");
    for(int i=0; i<=nn*M-1; i++)
    {
        if (i%M==0) printf("\n");
        printf("%.8f   ", Dt[i]);
    }*/
//    MAX=Max(p,nn,M);
//    MIN=Min(p,num,nn,M);
//    minDx=Min(Dx,0,nn,M);
//    maxDx=Max(Dx,nn,M);
//    minDt=Min(Dt,0,nn,M);
//    maxDt=Max(Dt,nn,M);
//    PRINT(p,t,h,s,MIN,MAX,minDx,maxDx,minDt,maxDt,nn,M,T,eps);
}


int main ()
{
    int i=0, nn=10,M=10,nn1=10,M1=10,nn2=10,M2=10,T=1, u, minDt1, minDt2, minDt3;
    double  b = 0.1, a = 0.01, eps=0.00000001, *Dt,*Dt1,*Dt2,t1,n1,n2,t,n3,t2;
    for(int j=0; j<2;j++) {

        M = 10 * (j + 1);
        nn = M * M * 3;
        n1 = nn;
        t = T / (n1 - 1);
        double p[nn * M - 1];
        for (i = 0; i < nn * M; i++) { p[i] = 100; }
        Dt=(double*)malloc(M*nn*sizeof(double));
        if(Dt==NULL){printf("6\n");return -1;}
        diff(a, b, p, nn, M, T, eps,Dt);
        minDt1 = Min(Dt, 0, nn, M);

        M1 = 10 * (j + 2);
        nn1 = M1 * M1 * 3;
        n2 = nn1;
        t1 = T / (n2 - 1);
        double p1[nn1 * M1 - 1];
        for (i = 0; i < nn1 * M1; i++) { p1[i] = 100; }
        Dt1=(double*)malloc(M1*nn1*sizeof(double));
        if(Dt1==NULL){printf("61\n");return -1;}
        diff(a, b, p1, nn1, M1, T, eps,Dt1);
        minDt2 = Min(Dt1, 0, nn1, M1);

        M2 = 10 * (j + 3);
        nn2 = M2 * M2 * 3;
        n3 = nn2;
        t2 = T / (n3 - 1);
        double p2[nn2 * M2 - 1];
        for (i = 0; i < nn2 * M2; i++) { p2[i] = 100; }
        Dt2=(double*)malloc(M2*nn2*sizeof(double));
        if(Dt2==NULL){printf("62\n");return -1;}
        diff(a, b, p2, nn2, M2, T, eps,Dt2);
        minDt3 = Min(Dt2, 0, nn2, M2);


        M=M/2;
        M1=M1/2;
        M2=M2/2;
        printf("0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",
               (p[M * (nn / 10)] - p1[M1 * (nn1 / 10)]) / (p1[M1 * (nn1 / 10)] - p2[M2 * (nn2 / 10)]),
               (p[M * (nn / 10) + M / 4] - p1[M1 * (nn1 / 10) + M1 / 4]) /
               (p1[M1 * (nn1 / 10) + M1 / 4] - p2[M2 * (nn2 / 10) + M2 / 4]),
               (p[M * (nn / 10) + M / 2] - p1[M1 * (nn1 / 10) + M1 / 2]) /
               (p1[M1 * (nn1 / 10) + M1 / 2] - p2[M2 * (nn2 / 10) + M2 / 2]),
               (p[M * (nn / 10) + 3 * M / 4] - p1[M1 * (nn1 / 10) + 3 * M1 / 4]) /
               (p1[M1 * (nn1 / 10) + 3 * M1 / 4] - p2[M2 * (nn2 / 10) + 3 * M2 / 4]),
               (p[M * (nn / 10) + M - 1] - p1[M1 * (nn1 / 10) + M1 - 1]) /
               (p1[M1 * (nn1 / 10) + M1 - 1] - p2[M2 * (nn2 / 10) + M2 - 1]));
        printf("\\\\");
        printf("\n");
        printf("$t_{min D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f",
               (p[M * (minDt1 / M)] - p1[M1 * (minDt2 / M1)]) / (p1[M * (minDt2 / M)] - p2[M2 * (minDt3 / M2)]),
               (p[M * (minDt1 / M) + M / 4] - p1[M1 * (minDt2 / M1) + M1 / 4]) /
               (p1[M1 * (minDt2 / M1) + M1 / 4] - p2[M2 * (minDt3 / M2) + M2 / 4]),
               (p[M * (minDt1 / M) + M / 2] - p1[M1 * (minDt2 / M1) + M1 / 2]) /
               (p1[M1 * (minDt2 / M1) + M1 / 2] - p2[M2 * (minDt3 / M2) + M2 / 2]),
               (p[M * (minDt1 / M) + 3 * M / 4] - p1[M1 * (minDt2 / M1) + 3 * M1 / 4]) /
               (p1[M1 * (minDt2 / M1) + 3 * M1 / 4] - p2[M2 * (minDt3 / M2) + 3 * M2 / 4]),
               (p[M * (minDt1 / M) + M - 1] - p1[M1 * (minDt2 / M1) + M1 - 1]) /
               (p1[M1 * (minDt2 / M1) + M1 - 1] - p2[M2 * (minDt3 / M2) + M2 - 1]));
        printf("\\\\");
        printf("\n");
        printf("0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",
               (p[9 * M * (nn / 10)] - p1[9 * M1 * (nn1 / 10)]) / (p1[9 * M1 * (nn1 / 10)] - p2[9 * M2 * (nn2 / 10)]),
               (p[9 * M * (nn / 10) + M / 4] - p1[9 * M1 * (nn1 / 10) + M1 / 4]) /
               (p1[9 * M1 * (nn1 / 10) + M1 / 4] - p2[9 * M2 * (nn2 / 10) + M2 / 4]),
               (p[9 * M * (nn / 10) + M / 2] - p1[9 * M1 * (nn1 / 10) + M1 / 2]) /
               (p1[9 * M1 * (nn1 / 10) + M1 / 2] - p2[9 * M2 * (nn2 / 10) + M2 / 2]),
               (p[9 * M * (nn / 10) + 3 * M / 4] - p1[9 * M1 * (nn1 / 10) + 3 * M1 / 4]) /
               (p1[9 * M1 * (nn1 / 10) + 3 * M1 / 4] - p2[9 * M2 * (nn2 / 10) + 3 * M2 / 4]),
               (p[9 * M * (nn / 10) + M - 1] - p1[9 * M1 * (nn1 / 10) + M1 - 1]) /
               (p1[9 * M1 * (nn1 / 10) + M1 - 1] - p2[9 * M2 * (nn2 / 10) + M2 - 1]));
        printf("\\\\");
        printf("\n");
        M1=M1*2;
        M=M*2;
        M2=M2*2;
    }
    return 0;
}



int Max(double *p,int nn, int M)
{
    double s=0.0;
    int i,k;
    for(i=0;i<nn*M;i++)
    {
        if(p[i]>=s && !(p[i]>=100 && p[i]<=100)){s=p[i];k=i;}
    }
    return k;
}
int Min(double *p,int w,int nn, int M)
{
    double s=2.0;
    int i,k;
    for(i=0;i<nn*M;i++)
    {
        if(p[i]<s && i>=w*M){s=p[i];k=i;}
    }
    return k;
}

void RESHI(double *a, double *x, double *b,int N)
{
    double *s = new double[N];
    double *p = new double[N];
    double *y = new double[N];

    y[0]=a[0];
    s[0]=-a[1]/y[0];
    p[0]=b[0]/y[0];
    for(int i=1;i<N;i++)
    {
        y[i]=a[N*i+i]+a[N*i+i-1]*s[i-1];
        s[i]=-a[N*i+i+1]/y[i];
        p[i]=(b[i]-a[N*i+i-1]*p[i-1])/y[i];
    }
    y[N-1]=a[N*N-1]+a[N*N-2]*s[N-2];
    p[N-1]=(b[N-1]-a[N*N-2]*p[N-2])/y[N-1];
    x[N-1]=p[N-1];
    for(int i=N-2;i>=0;i--)
    {
        x[i]=s[i]*x[i+1]+p[i];
    }
    return;
}
void makeAandB(double *p,double *A, double *B, double t, double h,double a, int j, int M)
{
//    int i,q1=0;
    A[0] = 1/t+a/(h*h)+(h*h)/2-a/(2*h*h*h*(1/h+h/(2*a*t)));
    A[1] = -a/(2*h*h);
    A[(M-2)*(M-2)-2]=-a/(2*h*h);
    A[(M-2)*(M-2)-1]=1/t+a/(h*h)+(h*(M-2)*h*(M-2))/2-a/(2*h*h*h*(1/h+(h)/(2*a*t)+(h/(2*a))));
    B[0]=p[(j+1)*(M)+1]*(1/t-(h*h)/2)+a*(p[(j+1)*M]-2*p[(j+1)*M+1]+p[(j+1)*M+2])/(2*h*h) -2
         +a*((h*p[(j+1)*M])/(2*a*t)+(h)/(2*a))/(2*h*h*(-1/h-(h)/(2*a*t)));
    B[M-3]=p[(j+1)*(M)+M-2]*(1/t-(h*(M-2)*(M-2)*h)/2)+a*(p[(j+1)*M+M-3]-2*p[(j+1)*M+M-2]+p[(j+1)*M+M-1])/(2*h*h) -1
           +a*((h*p[(j+1)*M+M-1])/(2*a*t)-(h)/(2*a))/(2*h*h*(1/h+(h)/(2*a*t)+h/(2*a)));

    for(int i=1; i<M-3; i++)
    {
        A[(M-2)*i+i-1]=-a/(2*h*h);
        A[(M-2)*i+i]=1/t+a/(h*h)+(h*h*(i+1)*(i+1))/2;

        A[(M-2)*i+i+1]=-a/(2*h*h);
        B[i]=p[(j+1)*(M)+1+i]*(1/t-(h*(i+1)*(i+1)*h)/2)+a*(p[(j+1)*M+i]-2*p[(j+1)*M+1+i]+p[(j+1)*M+2+i])/(2*h*h) -1;
    }
    return;
}
void makeDx(double *Dx, double *p, double h,int nn, int M)
{
    int i,j;
    for(j=0;j<nn-1;j++)
    {
        for(i=j*M+1;i<=(j+1)*M-2;i++)
        {
            Dx[i]=(p[i+1]-p[i-1])/2/h;
        }
    }
    return;
}
void makeDt(double *Dt, double *p, double t,int nn, int M)
{
    int i,j;
    for(j=1;j<nn-1;j++)
    {
        for(i=j*M;i<(j+1)*M;i++)
        {
            Dt[i]=(p[i+M]-p[i-M])/2/t;
        }
    }
    return;
}
void PRINT(double *p,double t,double h, double s,int x, int y,int z,int c,int u,int g,int nn, int M, int T, double eps)
{
    printf("\n");
    printf("\n");

}


