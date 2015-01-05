#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData
 { 
   double *V0, A[3], B[3], C[3];
   double *V0P, AP[3], BP[3], CP[3];
   double u[3];

 } IntegrandData;

int InnerIntegrand(unsigned ndim, const double *x, void *UserData,
                   unsigned fdim, double *fval)
{
  IntegrandData *Data=(IntegrandData *)UserData;
  double *u=Data->u;
  double *V0=Data->V0;
  double *A=Data->A;
  double *B=Data->B;
  double *C=Data->C;
  double *V0P=Data->V0P;
  double *AP=Data->AP;
  double *BP=Data->BP;
  double *CP=Data->CP;

  double Xi1=x[0];
  double Xi2=x[1];
  double Xi3=x[2];
  double Eta1=Xi1+u[0];
  double Eta2=Xi2+u[1];
  double Eta3=Xi3+u[2];

  double X[3], XP[3];
  for(int Mu=0; Mu<3; Mu++)
   { X[Mu] = V0[Mu] + Xi1*A[Mu] + Xi2*B[Mu] + Xi3*C[Mu];
     XP[Mu] = V0P[Mu] + Eta1*AP[Mu] + Eta2*BP[Mu] + Eta3*CP[Mu];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddInnerIntegral(double u1, double u2, double u3, void *UserData, double *I)
{
  IntegrandData *Data=(IntegrandData *)UserData;
  Data->u[0] = u1;
  Data->u[1] = u2;
  Data->u[2] = u3;

  double Lower[3], Upper[3];
  for(int i=0; i<3; i++)
   { if ( u[i] < 0.0 )
      { Lower[i] = -u[i];
        Upper[i] = 1.0;
      }
     else 
      { Lower[i] = 0.0;
        Upper[i] = 1.0-u[i];
      };
   };
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int OuterIntegrand(unsigned ndim, const double *x, void *UserData,
                   unsigned fdim, double *fval)
{

  double u1=x[0], Alpha2=x[1], Alpha3=x[2];
  double u2=u1*Alpha1, u3=u1*Alpha2*Alpha3;
  double Jacobian=Alpha1*u1*u1;

  memset(fval, 0, fdim*sizeof(double));
  for(int Sign1=1; Sign1>=-1; Sign1-=2)
   for(int Sign2=1; Sign2>=-1; Sign2-=2)
    for(int Sign3=1; Sign3>=-1; Sign3-=2)
     { 
       double su1 = Sign1*u1;
       double su2 = Sign2*u2;
       double su3 = Sign3*u3;

       AddInnerIntegral(su1, su2, su3, UserData, fval);
       AddInnerIntegral(su1, su3, su2, UserData, fval);
       AddInnerIntegral(su2, su1, su3, UserData, fval);
       AddInnerIntegral(su2, su3, su1, UserData, fval);
       AddInnerIntegral(su3, su1, su2, UserData, fval);
       AddInnerIntegral(su3, su2, su1, UserData, fval);
     };

  for(int nf=0; nf<fdim; nf++)
   fval[nf] *= Jacobian;

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
}
