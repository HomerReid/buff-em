#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData
 { 
   double A

 } IntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddInnerIntegral(double u1, double u2, double u3,
                      void *UserData, double *I)
{
 
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

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
}
