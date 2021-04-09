#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int glavny(int i, double *a, int N);
void swapstolb(int n, int m, double *a,  int N,int *c);
void kdiagonali( double *a,double *b,  int N, int *c);
void reshi(double *a, double *x, double *b,int N,int *c);
int Max(double *p, int nn, int M, int T, double eps);
int Min(double *p,int w, int nn, int M, int T, double eps);
void makeAandB(double *p,double *A, double *B, double t, double h,double a, double b, int nn, int M, int T, double eps);
void makeDx(double *Dx, double *p, double h, int nn, int M, int T, double eps);
void makeDt(double *Dt,double *p,double t, int nn, int M, int T, double eps);


double
lambda (double x)
{
    if(x<=1){
        return x * (1 -  x) * (1 -  x);}

}




double
diff (double a, double b, double p[], int nn, int M, int T, double eps)
{
    int i, j, k, *c, L, MAX, MIN, num = 0, minDx, maxDx, minDt, maxDt;
    double x = 0, *X, *A, *B, s, *Dx, *Dt;
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
        B = (double *) malloc ((M - 2) * sizeof (double));
        if (B == NULL)
        {
            printf ("2\n");
            return -1;
        }
        X = (double *) malloc ((M - 2) * sizeof (double));
        if (X == NULL)
        {
            printf ("3\n");
            return -1;
        }
        A = (double *) malloc ((M - 2) * (M - 2) * sizeof (double));
        if (A == NULL)
        {
            printf ("5\n");
            return -1;
        }
        c = (int *) malloc ((M - 2) * sizeof (int));
        if (c == NULL)
        {
            printf ("7\n");
            return -1;
        }
        for (i = 0; i < (M - 2); i++)
        {
            c[i] = i + 1;
        }
        for (i = 0; i < (M - 2) * (M - 2); i++)
        {
            A[i] = 0;
        }
        makeAandB (p, A, B, t, h, a, b, nn,M,T,eps);
        reshi (A, X, B, M - 2, c);
        L = 0;
        for (i = 1; i <= M - 2; i++)
        {
            p[i] = X[i - 1];
            if (fabs (p[ i] - p[(1) * M + i]) < eps)
            {
                L++;
            }

        }
        p[0] = (-p[1]/h+(h*p[(1)*M])/(2*a*t)+(h)/(2*a)-1)/((-1/h-(h)/(2*a*t))) ;
        p[M - 1] =(p[M - 2]/h+(h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(1/h+(h)/(2*a*t)+h/(2*a));
        for (i = 0; i <= M - 1; i++){

            //printf("%.8f   ", p[i+M]);
            p[i+M]=p[i];
        }
        //printf("\n");


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
}


int main ()
{
    int i=0, nn=10,M=10,nn1=10,M1=10,nn2=10,M2=10,T=1, u, minDt1, minDt2, minDt3;
    double  b = 0.1, a = 0.01, eps=0.00000001, *Dt,*Dt1,*Dt2,t1,n1,n2,t,n3,t2;
    for(int j=0; j<1;j++){
        M = 100*(j+1);
        nn=M*M*3;
        double p[2*M-1];
        for(i=0;i<nn*M;i++){p[i]=100;}
        diff(a,b,p,nn,M,T,eps);
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
       makeDt (Dt, p, t,nn,M,T,eps);
       minDt1=Min(Dt,0,nn,M,T,eps);
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
            diff(a,b,p,nn,M,T,eps);
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
        makeDt (Dt, p, t,nn,M,T,eps);
        minDt1=Min(Dt,0,nn,M,T,eps);




            M1 = 20*(j+2);
            nn1= M1*M1*3;
            n2=nn1;
            t1 = T / (n2 - 1);
            double p1[nn1*M1-1];
            for(i=0;i<nn1*M1;i++){p1[i]=100;}
            diff(a,b,p1,nn1,M1,T,eps);
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
            makeDt (Dt1, p1, t1,nn1,M1,T,eps);
            minDt2=Min(Dt1,0,nn1,M1,T,eps);



            M2 = 20*(j+3);
            nn2=M2*M2*3;
            n3=nn2;
            t2 = T / (n3 - 1);
            double p2[nn2*M2-1];
            for(i=0;i<nn2*M2;i++){p2[i]=100;}
            diff(a,b,p2,nn2,M2,T,eps);
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
            makeDt (Dt2, p2, t2,nn2,M2,T,eps);
            minDt3=Min(Dt2,0,nn2,M2,T,eps);
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



int Max(double *p,int nn, int M,int T, double eps)
{
    double s=0.0;
    int i,k;
    for(i=0;i<nn*M;i++)
    {
        if(p[i]>=s && !(p[i]>=100 && p[i]<=100)){s=p[i];k=i;}
    }
    return k;
}
int Min(double *p,int w, int nn, int M,int T, double eps)
{
    double s=2.0;
    int i,k;
    for(i=0;i<nn*M;i++)
    {
        if(p[i]<s && i>=w*M){s=p[i];k=i;}
    }
    return k;
}
int glavny(int i, double *a, int N)
{
    int j,k;
    double b;
    k=i*N;
    b=fabs(a[i*N]);
    for(j=i*N;j<(i+1)*N;j++){if(fabs(a[j])>b){b=fabs(a[j]);k=j;}}
    k=k-i*N;
    return k;
}
void swapstolb(int n, int m, double *a,  int N,int *c)
{
    int i,k;
    double b;

    k=c[n];
    c[n]=c[m];
    c[m]=k;
    for(i=0;i<N;i++)
    {
        b=a[i*N+n];
        a[i*N+n]=a[i*N+m];
        a[i*N+m]=b;
    }
    return;
}
void kdiagonali( double *a,double *b,  int N, int *c)
{
    int k,i,j,p;
    double s,f;
    for(i=0;i<N-1;i++)
    {
        k=glavny(i,a,N);
        // if(fabs(a[i*N+k])>=0 && fabs(a[i*N+k])<=0){printf("vyrozdena %d\n",i);}
        s=a[i*N+k];
        if(fabs(s)<=1e-16){return;}
        else{
            swapstolb(i,k,a,N,c);
            for(j=i*N;j<(i+1)*N;j++){a[j]=a[j]/s;}
            b[i]=b[i]/s;
            for(j=i+1;j<N;j++)
            {
                f=a[j*N+i];
                for(p=j*N+i;p<(j+1)*N;p++)
                {
                    a[p]=a[p]-f*a[p-j*N+i*N];
                }
                b[j]=b[j]-f*b[i];
            }
        }

    }
    b[N-1]=b[N-1]/a[(N-1)*N+N-1];
    a[(N-1)*N+N-1]=1;
    return;
}
void reshi(double *a, double *x, double *b,int N,int *c)
{
    int i,j;
    kdiagonali(a,b,N,c);
    /*for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%lf  ",a[i*N+j]);
        }
        printf("%lf  ",b[i]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
*/
    x[N-1]=b[N-1];
    for(i=N-2;i>=0;i--)
    {
        x[i]=b[i];
        for(j=N-1;j>i;j--)
        {
            x[i]=x[i]-x[j]*a[i*N+j];
        }
    }
    return;
}
void
makeAandB (double *p, double *A, double *B, double t, double h, double a,
           double b, int nn, int M,int T, double eps)
{
    int q1 = 0;
    int i;
    A[0] = 1/t+a/(h*h)+(h*h)/2-a/(2*h*h*h*(1/h+h/(2*a*t)));
    A[1] = -a/(2*h*h);
    A[(M-2)*(M-2)-2]=-a/(2*h*h);
    A[(M-2)*(M-2)-1]=1/t+a/(h*h)+(h*(M-2)*h*(M-2))/2-a/(2*h*h*h*(1/h+(h)/(2*a*t)+(h/(2*a))));
    B[0]=p[(1)*(M)+1]*(1/t-(h*h)/2)+a*(p[(1)*M]-2*p[(1)*M+1]+p[(1)*M+2])/(2*h*h) -2
         +a*((h*p[(1)*M])/(2*a*t)+(h)/(2*a))/(2*h*h*(-1/h-(h)/(2*a*t)));
    B[M-3]=p[(1)*(M)+M-2]*(1/t-(h*(M-2)*(M-2)*h)/2)+a*(p[(1)*M+M-3]-2*p[(1)*M+M-2]+p[(1)*M+M-1])/(2*h*h) -1
           +a*((h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(2*h*h*(1/h+(h)/(2*a*t)+h/(2*a)));
    for (int i = 1; i < M - 3; i++)

    {

        A[(M-2)*i+i-1]=-a/(2*h*h);
        A[(M-2)*i+i]=1/t+a/(h*h)+(h*h*(i+1)*(i+1))/2;

        A[(M-2)*i+i+1]=-a/(2*h*h);
        B[i]=p[(1)*(M)+1+i]*(1/t-(h*(i+1)*(i+1)*h)/2)+a*(p[(1)*M+i]-2*p[(1)*M+1+i]+p[(1)*M+2+i])/(2*h*h) -1;
    }
/*for (int i = 0; i <= (M - 2) * (M - 2) - 1; i++)
    {
if (i % (M - 2) == 0)
	{
	  printf ("\n");
	}
	if (i % (M - 2) == 0)
	{
	  printf ("%.8f   ", B[q1]);
	  q1 = q1 + 1;
	}
printf ("%.8f   ", A[i]);



}*/
//printf ("\n");
//printf ("\n");
//printf ("\n");
    return;
}


void makeDx(double *Dx, double *p, double h,int nn, int M,int T, double eps)
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
void makeDt(double *Dt, double *p, double t,int nn, int M,int T, double eps)
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

