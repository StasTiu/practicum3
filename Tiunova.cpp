#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void trehdiagonMatrix(double *a, double *x, double *b, int N);
void makeAandB(double *p,double *A, double *B, double t, double h,double a,int M);
void setka (double a, double *p, int nn, int M, int T, double *strMax, int numMax, double *str01, double *str09, double s);


int main ()
{
    int nn=10,M=10,T=1,numMax,M1,nn1, numMax1, nn2, numMax2, M2;
    double  s=0, a = 0.01,strMax[5],str01[5],str09[5],strMax1[5],str011[5],str091[5],str012[5],str092[5],strMax2[5];
    M = 50;
    nn=M*M*3;
    double *p,*p1,*p2;
    p = (double *)malloc((2*M)* sizeof (double));
    setka(a,p,nn,M,T,strMax,numMax,str01,str09,s);
//first 3 tables, change M to 50,100,150 and comment other three tables

    /*printf ( "0.1 & %.8f & %.8f  & %.8f & %.8f & %.8f",str01[0],str01[1],str01[2],str01[3],str01[4]);
      printf("\\\\");
      printf("\n");
      printf ( "$t_{Max D_t} %d$ & %.8f & %.8f  & %.8f & %.8f & %.8f",numMax, strMax[0],strMax[1],strMax[2],strMax[3],strMax[4]) ;
      printf("\\\\");
      printf("\n");
      printf ( "0.9 & %.8f & %.8f  & %.8f & %.8f & %.8f",str09[0],str09[1],str09[2],str09[3],str09[4]);
      printf("\\\\");
      printf("\n");*/

    M1 = 100;
    nn1=M1*M1*3;
    p1 = (double *)malloc((2*M1)* sizeof (double));
    setka(a,p1,nn1,M1,T,strMax1,numMax1,str011,str091,s);


    M2 = 150;
    nn2=M2*M2*3;
    p2 = (double *)malloc((2*M2)* sizeof (double));
    setka(a,p2,nn2,M2,T,strMax2,numMax2,str012,str092,s);
//table 4 (don't forget to change M to 50)
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
//table 5
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
//table 6
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

    return 0;
}
void setka (double a, double *p, int nn, int M, int T, double *strMax, int numMax, double *str01, double *str09, double s)
{
    double t, n1, M1, h;
    n1=nn;
    M1=M;
    t=T/(n1-1);
    h=2/(M1-1);
    for(int i=M; i<=M*2-1; i++){p[i]=1;}
    for(int j=nn-2; j>=0; j--){
        double *A,*B,*X;
        B = (double *) malloc ((M - 2) * sizeof (double));
        X = (double *) malloc ((M - 2) * sizeof (double));
        A = (double *) malloc ((M - 2) * 3 * sizeof (double));

        makeAandB(p,A,B,t,h,a,M);
        trehdiagonMatrix(A,X,B,M-2);
        for(int i=1;i<=M-2;i++){
            p[i]=X[i-1];
        }
        p[0]=1;
        p[M-1] = 1;
        if(j==(nn-1)/10){
            for(int i=0;i<5;i++)
                str09[i]=p[(i*M)/10+(M/5)];
        }
        if(j==9*(nn-1)/10){
            for(int i=0;i<5;i++)
                str01[i]=p[(i*M)/10+(M/5)];
        }
        for( int i=0;i<M-1;i++){
            if(fabs(p[i]-p[M+i])>=fabs(s)){
                s=p[i]-p[M+i];
                numMax=j;
                for(int k=0; k<5;k++){
                    strMax[k]=p[(k*M)/10+(M/5)];
                }
            }
        }
        for (int i = 0; i <= M - 1; i++){
            p[i+M]=p[i];
        }
    }
}

//function solves SLU with 3-diagonal matrix
void trehdiagonMatrix(double *a, double *x, double *b, int N)
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

void makeAandB(double *p,double *A, double *B, double t, double h,double a, int M)
{
    A[0]=1/t+a/(h*h)+(1-h)*(1-h);
    A[1]=-a/(2*h*h);
    A[2] = 0;
    A[(M-2)*3-2] = 0;
    A[(M-2)*3-1]=-a/(2*h*h);
    A[(M-2)*3-3]=1/t+a/(h*h)+(1-h)*(1-h);

    B[0]=a/(2*h*h)+p[(1)*(M)+1]*(1/t-(1+(1-h)*(1-h)))+a*(p[(1)*M]-2*p[(1)*M+1]+p[(1)*M+2])/(2*h*h);
    B[M-3]= a/(2*h*h)+a*(p[(1)*M+M-3]-2*p[(1)*M+M-2]+p[(1)*M+M-1])/(2*h*h)+p[(1)*(M)+M-2]*(1/t-(1+(1-h)*(1-h)));

    for(int i=1; i<M-3; i++)
    {
        A[3*i+2]=-a/(2*h*h);
        A[3*i]=1/t+a/(h*h)+(1+(h*(i+1)-1)*(h*(i+1)-1))/2;
        A[3*i+1]=-a/(2*h*h);
        B[i]=p[(1)*(M)+1+i]*(1/t-(1+(h*(i+1)-1)*(h*(i+1)-1))/2)+a*(p[(1)*M+i]-2*p[(1)*M+1+i]+p[(1)*M+2+i])/(2*h*h);
    }
    return;
}