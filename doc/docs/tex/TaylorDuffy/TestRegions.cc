#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>

#if 0
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta=0.05;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Delta",    PA_DOUBLE,  1, 1, (void *)&Delta,        0,  "delta"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  double Volume=0.0;
  FILE *f=fopen("TestRegions.out","w");
  for(double u=0.0; u<=1.0; u+=Delta)
   for(double v=0.0; v<=1.0-u; v+=Delta)
    for(double w=0.0; w<=1.0-u-w; w+=Delta)
     for(double up=0.0; up<=1.0; up+=Delta)
      for(double vp=0.0; vp<=1.0-up; vp+=Delta)
       for(double wp=0.0; wp<=1.0-up-wp; wp+=Delta)
        {
          Volume += pow(Delta,6.0);
          double Xi1 = u-up;
          double Xi2 = v-vp;
          double Xi3 = w-wp;
          double Eta1 = u+up;
          double Eta2 = v+vp;
          double Eta3 = w+wp;
          fprintf(f,"%e %e %e %e %e %e \n",Xi1,Xi2,Xi3,Eta1,Eta2,Eta3);
        };
  fclose(f);
  printf("Volume=%e (should be %e)\n",Volume,1.0/36.0);

}
#endif

int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double DN=10000000;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"N",        PA_DOUBLE,     1, 1, (void *)&DN,            0,  "N"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  srand48(time(0));

  unsigned long N = (unsigned long)DN;
  unsigned long NRegion[12];
  memset(NRegion,0,12*sizeof(NRegion[0]));
  FILE *f=fopen("TestRegions.out","w");
  for(unsigned long n=0; n<N; n++)
   { 
     double Xi1=drand48();
     double Xi2=drand48();
     double Xi3=drand48();

     double Eta1=drand48();
     double Eta2=drand48();
     double Eta3=drand48();

     if (Xi2>Xi1) continue;
     if (Xi3>Xi2) continue;
     if (Eta2>Eta1) continue;
     if (Eta3>Eta2) continue;

     double u1 = Eta1 - Xi1;
     double u2 = Eta2 - Xi2;
     double u3 = Eta3 - Xi3;

     int Region=-1;
     if ( u3>=0.0 )
      { if ( u1>=0.0 && u2>=0.0 && u2>=u1)      Region=0;
        else if ( u1>=0.0 && u2>=0.0 && u2<=u1) Region=1;
        else if ( u1>=0.0 && u2<=0.0)           Region=2;
        else if ( u1<=0.0 && u2<=0.0 && u2<=u1) Region=3;
        else if ( u1<=0.0 && u2<=0.0 && u2>=u1) Region=4;
        else if ( u1<=0.0 && u2>=0.0)           Region=5;
      }
     else
      { if ( u1>=0.0 && u2>=0.0 && u2>=u1)      Region=6;
        else if ( u1>=0.0 && u2>=0.0 && u2<=u1) Region=7;
        else if ( u1>=0.0 && u2<=0.0)           Region=8;
        else if ( u1<=0.0 && u2<=0.0 && u2<=u1) Region=9;
        else if ( u1<=0.0 && u2<=0.0 && u2>=u1) Region=10;
        else if ( u1<=0.0 && u2>=0.0)           Region=11;
      };

     NRegion[Region]++;

     fprintf(f,"%i %e %e %e %e %e %e %e %e %e \n",Region,u1,u2,u3,Xi1,Xi2,Xi3,Eta1,Eta2,Eta3);
   };
  fclose(f);

  double TotalVolume=0.0;
  for(int nr=0; nr<12; nr++)
   { double Volume = ((double)NRegion[nr]) / DN;
     TotalVolume+=Volume;
     printf("Volume(%2i)=%.6f\n",nr+1,Volume);
   };
  printf("Total = %.6f (should be %e)\n",TotalVolume,1.0/36.0);

}

#if 0
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta=0.05;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"Delta",    PA_DOUBLE,  1, 1, (void *)&Delta,        0,  "delta"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  FILE *f=fopen("TestRegions.out","w");
  for(double Xi1=0.0; Xi1<=1.0; Xi1+=Delta)
   for(double Xi2=0.0; Xi2<=Xi1; Xi2+=Delta)
    for(double Alpha3=0.0; Alpha3<1.1; Alpha3+=0.25)
     for(double Eta1=0.0; Eta1<=1.0; Eta1+=Delta)
      for(double Eta2=0.0; Eta2<=Eta1; Eta2+=Delta)
       for(double Beta3=0.0; Beta3<1.1; Beta3+=0.25)
        { 
          double Xi3  = Alpha3*Xi2;
          double Eta3 = Beta3*Eta2;

          double u1=Eta1-Xi1;
          double u2=Eta2-Xi2;
          double u3=Eta3-Xi3;

          fprintf(f,"%e %e %e %e %e %e \n",u1,u2,u3,Xi1,Xi2,Xi3);
        };
  fclose(f);
}
#endif
