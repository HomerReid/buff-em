#################################################
#################################################
#################################################
function GetXiLimits(u1, u2, u3, N=1000000) 
  Xi1Max=0;
  Xi2Max=0;
  Xi3Max=0;

  Xi1Min=1;
  Xi2Min=1;
  Xi3Min=1;

  for n=1:N
    Eta1=rand();
    Eta2=Eta1*rand();
    Eta3=Eta2*rand();
    Xi1 = Eta1 - u1;
    Xi2 = Eta2 - u2;
    Xi3 = Eta3 - u3;

    if (Xi1<0.0 || Xi1>1.0 || Xi2<0.0 || Xi2>Xi1 || Xi3<0.0 || Xi3>Xi2) 
     continue;
    end

    Xi1Max = max(Xi1, Xi1Max);
    Xi2Max = max(Xi2, Xi2Max);
    Xi3Max = max(Xi3, Xi3Max);

    Xi1Min = min(Xi1, Xi1Min);
    Xi2Min = min(Xi2, Xi2Min);
    Xi3Min = min(Xi3, Xi3Min);
  end

  [Xi1Min Xi1Max; 
   Xi2Min Xi2Max; 
   Xi3Min Xi3Max]

end

##################################################
##################################################
##################################################
function GetSector(u1,u2,u3)

  if     (u1>0.0 && u2>0.0 && u3>0.0)
    Q1=1;
  elseif (u1>0.0 && u2>0.0 && u3<0.0)
    Q1=2;
  elseif (u1>0.0 && u2<0.0 && u3>0.0)
    Q1=3;
  elseif (u1>0.0 && u2<0.0 && u3<0.0)
    Q1=4;
  elseif (u1<0.0 && u2>0.0 && u3>0.0)
    Q1=5;
  elseif (u1<0.0 && u2>0.0 && u3<0.0)
    Q1=6;
  elseif (u1<0.0 && u2<0.0 && u3>0.0)
    Q1=7;
  elseif (u1<0.0 && u2<0.0 && u3<0.0)
    Q1=8;
  end

  if     (abs(u1)>abs(u2) && abs(u2)>abs(u3))
    Q2=1;
  elseif (abs(u1)>abs(u3) && abs(u3)>abs(u2))
    Q2=2;
  elseif (abs(u2)>abs(u1) && abs(u1)>abs(u3))
    Q2=3;
  elseif (abs(u2)>abs(u3) && abs(u3)>abs(u1))
    Q2=4;
  elseif (abs(u3)>abs(u1) && abs(u1)>abs(u2))
    Q2=5;
  elseif (abs(u3)>abs(u2) && abs(u2)>abs(u1))
    Q2=6;
  end

 Q1, Q2
end

#################################################
#################################################
#################################################
function DoSector(u1,u2,u3)

  (Q1, Q2) = GetSector(u1, u2, u3);
  L=GetXiLimits(u1,u2,u3);

  if ( (L[1,2] <= L[1,1]) || (L[2,2] <= L[2,1]) || (L[3,2] <= L[3,1]) )
    @printf("%i%i: 0 \n",Q1,Q2);
    return
  end

L[3,2]-=L[2,2]; L[2,2]-=L[1,2];

  Names=["Xi1Lower" "Xi1Upper";
         "Xi2Lower" "Xi2Upper";
         "Xi3Lower" "Xi3Upper"];

  Delta=0.1;
 
  LP=GetXiLimits(u1+Delta, u2, u3); LP[3,2]-=LP[2,2]; LP[2,2]-=LP[1,2];
  B=(LP - L) / Delta;

  LP=GetXiLimits(u1, u2+Delta, u3); LP[3,2]-=LP[2,2]; LP[2,2]-=LP[1,2];
  C=(LP - L) / Delta;

  LP=GetXiLimits(u1, u2, u3+Delta); LP[3,2]-=LP[2,2]; LP[2,2]-=LP[1,2];
  D=(LP - L) / Delta;

  for m=1:3, n=1:2
    @printf("%i%i %s: ",Q1,Q2,Names[m,n]);

    if -1.2 < B[m,n] < -0.8 
      @printf(" -u1 ");
    elseif 0.8 < B[m,n] < 1.2 
      @printf(" +u1 ");
    end

    if -1.2 < C[m,n] < -0.8 
      @printf(" -u2 ");
    elseif 0.8 < C[m,n] < 1.2 
      @printf(" +u2 ");
    end

    if -1.2 < D[m,n] < -0.8 
      @printf(" -u3 ");
    elseif 0.8 < D[m,n] < 1.2 
      @printf(" +u3 ");
    end
    @printf("\n");

  end

end

function DoAllSectors()

  u10=0.1; 
  u20=0.3; 
  u30=0.5; 

  for s1=-1:2:1, s2=-1:2:1, s3=-1:2:1
   u1=s1*u10;
   u2=s2*u20;
   u3=s3*u30;
   DoSector(u1, u2, u3);
   DoSector(u1, u3, u2);
   DoSector(u2, u1, u3);
   DoSector(u2, u3, u1);
   DoSector(u3, u1, u2);
   DoSector(u3, u2, u1);
  end 
  

end
