#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void RESHI(double *a, double *x, double *b,int N);
int Max(double *p, int nn, int M);
int Min(double *p,int w, int nn, int M);
void makeAandB(double *p,double *A, double *B, double t, double h,double a, double b, int nn, int M);
void makeDt(double *Dt,double *p,double t, int nn, int M);


double
lambda (double x)
{
    if(x<=1){
        return x * (1 -  x) * (1 -  x);}
    else return 0;
}


double
diff (double a, double b, double p[], int nn, int M, int T)
{
    int  MAX, MIN, num = 0, minDx, maxDx, minDt, maxDt;
    double x = 0, s, *Dx, *Dt;
    double t, n1, M1, h;
    n1 = nn;
    M1 = M;
    t = T / (n1 - 1);
    h = 1 / (M1 - 1);
    for (int i = M  ; i <= M * 2 - 1; i++)
    {
        p[i] = lambda ((i-M)* h);
    }




    for (int j = nn - 2; j >= 0; j--)
    {
        double *A,*B,*X;
        B = (double *) malloc ((M - 2) * sizeof (double));
        X = (double *) malloc ((M - 2) * sizeof (double));
        A = (double *) malloc ((M - 2) * 3 * sizeof (double));

        makeAandB (p, A, B, t, h, a, b, nn,M);
        RESHI (A, X, B, M - 2);

        for (int i = 1; i <= M - 2; i++)
        {
            p[i] = X[i - 1];
        }
        for (int i =0; i <= M - 3; i++)
        {
            printf("%.8f   ", X[i]);
        }

        p[0] = (-p[1]/h+(h*p[(1)*M])/(2*a*t)+(h)/(2*a)-1)/((-1/h-(h)/(2*a*t))) ;
        p[M - 1] =(p[M - 2]/h+(h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(1/h+(h)/(2*a*t)+h/(2*a));
        for (int i = 0; i <= M - 1; i++){

            // printf("%.8f   ", p[i+M]);
            p[i+M]=p[i];
        }
        printf("\n");


    }


    for(int i=0; i<=2*M-1; i++)
    {
        if (i%M==0) printf("\n");
        printf("%.8f   ", p[i]);
    }
    printf("\n");
    printf("\n");/*
       for(int i=0; i<=nn*M-1; i++)
       {
       if (i%M==0) printf("\n");
       printf("%.8f   ", Dx[i]);
       }
       printf("\n");*/
    return 0;
}


int main ()
{
    int nn=10,M=10,T=1;
    double  b = 0.1, a = 0.01;
    for(int j=0; j<1;j++){
        M = 100*(j+1);
        nn=10;
        double *p;
        p = (double *)malloc((2*M)* sizeof (double));
        for(int i=0;i<2*M;i++){p[i]=100;}
        diff(a,b,p,nn,M,T);
        /*                   Dt = (double *) malloc (M * nn * sizeof (double));
          if (Dt == NULL)
          {
         printf ("6\n");
         return -1;
       }
       for (i = 0; i < M * nn; i++)
       {
         Dt[i] = 100;
       }
       makeDt (Dt, p, t,nn,M);
       minDt1=Min(Dt,0,nn,M);
       printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f", p[M*(nn/10)],p[M*(nn/10)+M/4],p[M*(nn/10)+M/2],p[M*(nn/10)+3*M/4],p[M*(nn/10)+M-1]) ;
       printf("\\\\");
       printf("\n");
       printf ( "$t_{min D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", p[M*(minDt1/M)],p[M*(minDt1/M)+M/4],p[M*(minDt1/M)+M/2],p[M*(minDt1/M)+3*M/4],p[M*(minDt1/M)+M-1]) ;
       printf("\\\\");
       printf("\n");
       printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f", p[M*9*(nn/10)],p[M*9*(nn/10)+M/4],p[M*9*(nn/10)+M/2],p[M*9*(nn/10)+3*M/4],p[M*9*(nn/10)+M-1]);
       printf("\\\\");
        printf("\n");
       }*/
        /*for(int j=0; j<2;j++){
            M = 20*(j+1);
            nn=M*M*3;
            n1=nn;
            t = T / (n1 - 1);
            double p[nn*M-1];
            for(i=0;i<nn*M;i++){p[i]=100;}
            diff(a,b,p,nn,M,T);
                    Dt = (double *) malloc (M * nn * sizeof (double));
           if (Dt == NULL)
           {
          printf ("6\n");
          return -1;
        }
        for (i = 0; i < M * nn; i++)
        {
          Dt[i] = 100;
        }
        makeDt (Dt, p, t,nn,M);
        minDt1=Min(Dt,0,nn,M);
            M1 = 20*(j+2);
            nn1= M1*M1*3;
            n2=nn1;
            t1 = T / (n2 - 1);
            double p1[nn1*M1-1];
            for(i=0;i<nn1*M1;i++){p1[i]=100;}
            diff(a,b,p1,nn1,M1,T);
                            Dt1 = (double *) malloc (M1 * nn1 * sizeof (double));
           if (Dt1 == NULL)
           {
          printf ("6\n");
          return -1;
        }
        for (i = 0; i < M1 * nn1; i++)
        {
          Dt1[i] = 100;
        }
            makeDt (Dt1, p1, t1,nn1,M1);
            minDt2=Min(Dt1,0,nn1,M1);
            M2 = 20*(j+3);
            nn2=M2*M2*3;
            n3=nn2;
            t2 = T / (n3 - 1);
            double p2[nn2*M2-1];
            for(i=0;i<nn2*M2;i++){p2[i]=100;}
            diff(a,b,p2,nn2,M2,T);
                              Dt2 = (double *) malloc (M2 * nn2 * sizeof (double));
           if (Dt2 == NULL)
           {
          printf ("6\n");
          return -1;
        }
        for (i = 0; i < M1 * nn1; i++)
        {
          Dt2[i] = 100;
        }
            makeDt (Dt2, p2, t2,nn2,M2);
            minDt3=Min(Dt2,0,nn2,M2);
        printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",(p[M*(nn/10)]-p1[M1*(nn1/10)]),(p[M*(nn/10)+M/4]-p1[M1*(nn1/10)+M1/4]),(p[M*(nn/10)+M/2]-p1[M1*(nn1/10)+M1/2]),(p[M*(nn/10)+3*M/4]-p1[M1*(nn1/10)+3*M1/4]),(p[M*(nn/10)+M-1]-p1[M1*(nn1/10)+M1-1])) ;
        printf("\\\\");
        printf("\n");
        printf ( "$t_{min D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", (p[M*(minDt1/M)]),(p[M*(minDt1/M)+M/4]-p1[M1*(minDt2/M1)+M1/4]),(p[M*(minDt1/M)+M/2]-p1[M1*(minDt2/M1)+M1/2]),(p[M*(minDt1/M)+3*M/4]-p1[M1*(minDt2/M1)+3*M1/4]),(p[M*(minDt1/M)+M-1]-p1[M1*(minDt2/M1)+M1-1])) ;
        printf("\\\\");
        printf("\n");
        printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",(p[9*M*(nn/10)]-p1[9*M1*(nn1/10)]),(p[9*M*(nn/10)+M/4]-p1[9*M1*(nn1/10)+M1/4]),(p[9*M*(nn/10)+M/2]-p1[9*M1*(nn1/10)+M1/2]),(p[9*M*(nn/10)+3*M/4]-p1[9*M1*(nn1/10)+3*M1/4]),(p[9*M*(nn/10)+M-1]-p1[9*M1*(nn1/10)+M1-1])) ;
        printf("\\\\");
        printf("\n");
        printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",(p[M*(nn/10)]-p1[M1*(nn1/10)])/(p1[M1*(nn1/10)]-p2[M2*(nn2/10)]),(p[M*(nn/10)+M/4]-p1[M1*(nn1/10)+M1/4])/(p1[M1*(nn1/10)+M1/4]-p2[M2*(nn2/10)+M2/4]),(p[M*(nn/10)+M/2]-p1[M1*(nn1/10)+M1/2])/(p1[M1*(nn1/10)+M1/2]-p2[M2*(nn2/10)+M2/2]),(p[M*(nn/10)+3*M/4]-p1[M1*(nn1/10)+3*M1/4])/(p1[M1*(nn1/10)+3*M1/4]-p2[M2*(nn2/10)+3*M2/4]),(p[M*(nn/10)+M-1]-p1[M1*(nn1/10)+M1-1])/(p1[M1*(nn1/10)+M1-1]-p2[M2*(nn2/10)+M2-1])) ;
        printf("\\\\");
        printf("\n");
        printf ( "$t_{min D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", (p[M*(minDt1/M)]-p1[M1*(minDt2/M1)])/(p1[M*(minDt2/M)]-p2[M2*(minDt3/M2)]),(p[M*(minDt1/M)+M/4]-p1[M1*(minDt2/M1)+M1/4])/(p1[M1*(minDt2/M1)+M1/4]-p2[M2*(minDt3/M2)+M2/4]),(p[M*(minDt1/M)+M/2]-p1[M1*(minDt2/M1)+M1/2])/(p1[M1*(minDt2/M1)+M1/2]-p2[M2*(minDt3/M2)+M2/2]),(p[M*(minDt1/M)+3*M/4]-p1[M1*(minDt2/M1)+3*M1/4])/(p1[M1*(minDt2/M1)+3*M1/4]-p2[M2*(minDt3/M2)+3*M2/4]),(p[M*(minDt1/M)+M-1]-p1[M1*(minDt2/M1)+M1-1])/(p1[M1*(minDt2/M1)+M1-1]-p2[M2*(minDt3/M2)+M2-1])) ;
        printf("\\\\");
        printf("\n");
      printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",(p[9*M*(nn/10)]-p1[9*M1*(nn1/10)])/(p1[9*M1*(nn1/10)]-p2[9*M2*(nn2/10)]),(p[9*M*(nn/10)+M/4]-p1[9*M1*(nn1/10)+M1/4])/(p1[9*M1*(nn1/10)+M1/4]-p2[9*M2*(nn2/10)+M2/4]),(p[9*M*(nn/10)+M/2]-p1[9*M1*(nn1/10)+M1/2])/(p1[9*M1*(nn1/10)+M1/2]-p2[9*M2*(nn2/10)+M2/2]),(p[9*M*(nn/10)+3*M/4]-p1[9*M1*(nn1/10)+3*M1/4])/(p1[9*M1*(nn1/10)+3*M1/4]-p2[9*M2*(nn2/10)+3*M2/4]),(p[9*M*(nn/10)+M-1]-p1[9*M1*(nn1/10)+M1-1])/(p1[9*M1*(nn1/10)+M1-1]-p2[9*M2*(nn2/10)+M2-1])) ;
        printf("\\\\");
        printf("\n");*/
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
int Min(double *p,int w, int nn, int M)
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
    double *s,*p,*y;
    s = (double *) malloc (N * sizeof (double));
    p = (double *) malloc (N* sizeof (double));
    y = (double *) malloc (N* sizeof (double));
    y[0]=a[0];
    s[0]=-a[1]/y[0];
    p[0]=b[0]/y[0];
    for(int i=1;i<N-1;i++)
    {
        y[i]=a[3*i]+a[3*i+2]*s[i-1];
        s[i]=-a[3*i+1]/y[i];
        p[i]=(b[i]-a[3*i+2]*p[i-1])/y[i];
    }
    y[N-1]=a[(N-1)*3]+a[(N-1)*3+2]*s[N-2];
    p[N-1]=(b[N-1]-a[(N-1)*3+2]*p[N-2])/y[N-1];
    x[N-1]=p[N-1];
    for(int i=N-2;i>=0;i--)
    {
        x[i]=s[i]*x[i+1]+p[i];
    }
    return;
}
void makeAandB (double *p, double *A, double *B, double t, double h, double a,
                double b, int nn, int M)
{
    int q1=0;
    A[0] = 1/t+a/(h*h)+(h*h)/2-a/(2*h*h*h*(1/h+h/(2*a*t)));
    A[1] = -a/(2*h*h);
    A[2] = 0;
    A[(M-2)*3-2] = 0;
    A[(M-2)*3-1]=-a/(2*h*h);
    A[(M-2)*3-3]=1/t+a/(h*h)+(h*(M-2)*h*(M-2))/2-a/(2*h*h*h*(1/h+(h)/(2*a*t)+(h/(2*a))));
    B[0]=p[(1)*(M)+1]*(1/t-(h*h)/2)+a*(p[(1)*M]-2*p[(1)*M+1]+p[(1)*M+2])/(2*h*h) -2
         +a*((h*p[(1)*M])/(2*a*t)+(h)/(2*a))/(2*h*h*(-1/h-(h)/(2*a*t)));
    B[M-3]=p[(1)*(M)+M-2]*(1/t-(h*(M-2)*(M-2)*h)/2)+a*(p[(1)*M+M-3]-2*p[(1)*M+M-2]+p[(1)*M+M-1])/(2*h*h) -1
           +a*((h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(2*h*h*(1/h+(h)/(2*a*t)+h/(2*a)));
    for (int i = 1; i < M - 3; i++)

    {

        A[3*i+1]=-a/(2*h*h);
        A[3*i]=1/t+a/(h*h)+(h*h*(i+1)*(i+1))/2;

        A[3*i+2]=-a/(2*h*h);
        B[i]=0;
    }
/*for (int i = 0; i <= (M - 2) * 3 - 1; i++)
    {
if (i % (M - 2) == 0)
	{
	  printf ("\n");
	}
	if (i % 3 == 0)
	{
	  printf ("%.8f   ", B[8]);
	  q1 = q1 + 1;
	}
printf ("%.8f   ", A[i]);
}
//printf ("\n");
//printf ("\n");
//printf ("\n");*/
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