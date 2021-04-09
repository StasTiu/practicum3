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
diff (double a, double b, double p[], int nn, int M, int T,double strMax[],int numMax,double str01[],double str09[],double s)
{
    int  MAX, MIN, num = 0, minDx, maxDx, minDt, maxDt;
    double x = 0, *Dx, *Dt;
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

        p[0] = (-p[1]/h+(h*p[(1)*M])/(2*a*t)+(h)/(2*a)-1)/((-1/h-(h)/(2*a*t))) ;
        p[M - 1] =(p[M - 2]/h+(h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(1/h+(h)/(2*a*t)+h/(2*a));
        if(j==(nn-1)/10){
            for(int i=0;i<5;i++)
                str09[i]=p[(i*M)/5];
            // printf("%d////", j);
            // printf("\n");
        }
        if(j==9*(nn-1)/10){
            for(int i=0;i<5;i++)
                str01[i]=p[(i*M)/5];
            //printf("%d////", j);
        }
        for( int i=0;i<M-1;i++)

        {
            if(fabs(p[i]-p[M+i])>=fabs(s)){s=p[i]-p[M+i];numMax=j;
                //printf("%d////", numMax);
                //printf("%.8f   ", p[i]-p[M+i]);
                for(int k=0; k<5;k++){
                    strMax[k]=p[(k*M)/5];
                }
            }
        }
        for (int i = 0; i <= M - 1; i++){

            //printf("%.8f   ", p[i+M]);
            p[i+M]=p[i];
        }
        //printf("\n");


    }


    /*for(int i=0; i<=2*M-1; i++)
    {
        if (i%M==0) printf("\n");
        printf("%.8f   ", p[i]);
    }*/
    //printf("\n");
    // printf("\n");
    /*for(int i=0; i<=nn*M-1; i++)
    {
    if (i%M==0) printf("\n");
    printf("%.8f   ", Dx[i]);
    }
    printf("\n");*/
    return 0;
}


int main ()
{
    int nn=10,M=10,T=1,numMax,M1,nn1, numMax1, nn2, numMax2, M2;
    double  s=0,b = 0.1, a = 0.01,strMax[5],str01[5],str09[5],strMax1[5],str011[5],str091[5],n2,t2,str012[5],str092[5],n3,t3,strMax2[5];
    for(int j=0; j<3;j++){
        M = 50*(j+1);
        nn=M*M*3;
        double *p,*p1,*p2;
        p = (double *)malloc((2*M)* sizeof (double));
        diff(a,b,p,nn,M,T,strMax,numMax,str01,str09,s);
        /*printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",str01[0],str01[1],str01[2],str01[3],str01[4]);
        printf("\\\\");
        printf("\n");
        printf ( "$t_{Max D_t} %d$ & %.8f & %.8f  & %.8f & %.8f & %.8f",numMax, strMax[0],strMax[1],strMax[2],strMax[3],strMax[4]) ;
        printf("\\\\");
        printf("\n");
        printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",str09[0],str09[1],str09[2],str09[3],str09[4]);
        printf("\\\\");
        printf("\n");*/

        M1 = 100*(j+1);
        nn1=M1*M1*3;
        p1 = (double *)malloc((2*M1)* sizeof (double));
        diff(a,b,p1,nn1,M1,T,strMax1,numMax1,str011,str091,s);


        M2 = 150*(j+1);
        nn2=M2*M2*3;
        p2 = (double *)malloc((2*M2)* sizeof (double));
        diff(a,b,p2,nn2,M2,T,strMax2,numMax2,str012,str092,s);

        printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",str01[0]-str011[0],str01[1]-str011[1],
                 str01[2]-str011[2],str01[3]-str011[3],str01[4]-str011[4]);
        printf("\\\\");
        printf("\n");
        printf ( "$t_{Max D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", strMax[0]-strMax1[0],strMax[1]-strMax1[1],
                 strMax[2]-strMax1[2],strMax[3]-strMax1[3],strMax[4]-strMax1[4]) ;
        printf("\\\\");
        printf("\n");
        printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",str09[0]-str091[0],str09[1]-str091[1],
                 str09[2]-str091[2],str09[3]-str091[3],str09[4]-str091[4]);
        printf("\\\\");

        printf("\n");
        printf("\n");
        printf("\n");

        printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",str011[0]-str012[0],str011[1]-str012[1],str011[2]-str012[2],
                 str011[3]-str012[3],str011[4]-str012[4]);
        printf("\\\\");
        printf("\n");
        printf ( "$t_{Max D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", strMax1[0]-strMax2[0],strMax1[1]-strMax2[1],
                 strMax1[2]-strMax2[2],strMax1[3]-strMax2[3],strMax1[4]-strMax2[4]) ;
        printf("\\\\");
        printf("\n");
        printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",str091[0]-str092[0],str091[1]-str092[1],str091[2]-str092[2],
                 str091[3]-str092[3],str091[4]-str092[4]);
        printf("\\\\");

        printf("\n");
        printf("\n");
        printf("\n");

        printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",(str01[0]-str011[0])/(str011[0]-str012[0]),
                 (str01[1]-str011[1])/(str011[1]-str012[1]),(str01[2]-str011[2])/(str011[2]-str012[2]),
                 (str01[3]-str011[3])/(str011[3]-str012[3]),(str01[4]-str011[4])/(str011[4]-str012[4]));
        printf("\\\\");
        printf("\n");
        printf ( "$t_{Max D_t}$ & %.8f & %.8f  & %.8f & %.8f & %.8f", (strMax[0]-strMax1[0])/( strMax1[0]-strMax2[0]),
                 (strMax[1]-strMax1[1])/(strMax1[1]-strMax2[1]),(strMax[2]-strMax1[2])/(strMax1[2]-strMax2[2]),
                 (strMax[3]-strMax1[3])/(strMax1[3]-strMax2[3]),(strMax[4]-strMax1[4])/(strMax1[4]-strMax2[4]));
        printf("\\\\");
        printf("\n");
        printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",(str09[0]-str091[0])/(str091[0]-str092[0]),
                 (str09[1]-str091[1])/(str091[1]-str092[1]),(str09[2]-str091[2])/(str091[2]-str092[2]),
                 (str09[3]-str091[3])/(str091[3]-str092[3]),(str09[4]-str091[4])/(str091[4]-str092[4]));
        printf("\\\\");


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
    B[0]=p[(1)*(M)+1]*(1/t-(h*h)/2)+a*(p[(1)*M]-2*p[(1)*M+1]+p[(1)*M+2])/(2*h*h) -1
         +a*((h*p[(1)*M])/(2*a*t)+(h)/(2*a))/(2*h*h*(-1/h-(h)/(2*a*t)));
    B[M-3]=p[(1)*(M)+M-2]*(1/t-(h*(M-2)*(M-2)*h)/2)+a*(p[(1)*M+M-3]-2*p[(1)*M+M-2]+p[(1)*M+M-1])/(2*h*h) -1
           +a*((h*p[(1)*M+M-1])/(2*a*t)-(h)/(2*a))/(2*h*h*(1/h+(h)/(2*a*t)+h/(2*a)));
    for (int i = 1; i < M - 3; i++)

    {

        A[3*i+1]=-a/(2*h*h);
        A[3*i]=1/t+a/(h*h)+(h*h*(i+1)*(i+1))/2;

        A[3*i+2]=-a/(2*h*h);
        B[i]=p[(1)*(M)+1+i]*(1/t-(1+(h*(i+1)-1)*(h*(i+1)-1))/2)+a*(p[(1)*M+i]-2*p[(1)*M+1+i]+p[(1)*M+2+i])/(2*h*h)-1;
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