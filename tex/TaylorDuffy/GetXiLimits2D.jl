#################################################
#################################################
#################################################
function GetXiLimits(u1, u2, N=1000000) 
  Xi1Max=0;
  Xi2Max=0;

  Xi1Min=1;
  Xi2Min=1;

  for n=1:N
    Eta1=rand();
    Eta2=Eta1*rand();
    Xi1 = Eta1 - u1;
    Xi2 = Eta2 - u2;

    if (Xi1<0.0 || Xi1>1.0 || Xi2<0.0 || Xi2>Xi1)
     continue;
    end

    Xi1Max = max(Xi1, Xi1Max);
    Xi2Max = max(Xi2, Xi2Max);

    Xi1Min = min(Xi1, Xi1Min);
    Xi2Min = min(Xi2, Xi2Min);
  end

  [Xi1Min Xi1Max; 
   Xi2Min Xi2Max];

end

##################################################
##################################################
##################################################
function GetSector(u1,u2)

  if     (u1>0.0 && u2>0.0 && u1>u2)
    Sector=2;
  elseif (u1>0.0 && u2>0.0 && u1<u2)
    Sector=6;
  elseif (u1>0.0 && u2<0.0)
    Sector=4;
  elseif (u1<0.0 && u2<0.0 && u1<u2)
    Sector=1;
  elseif (u1<0.0 && u2<0.0 && u1>u2)
    Sector=5;
  elseif (u1<0.0 && u2>0.0)
    Sector=3;
  end

 Sector
end

#################################################
#################################################
#################################################
function DoSector(u1,u2)

  Q = GetSector(u1, u2)
  L=GetXiLimits(u1,u2);

  if ( (L[1,2] <= L[1,1]) || (L[2,2] <= L[2,1]))
    @printf("%i: 0 \n",Q);
    return
  end

  Names=["Xi1Lower" "Xi1Upper";
         "Xi2Lower" "Xi2Upper"];

  Delta=0.1;
  B=(GetXiLimits(u1+Delta, u2) - L) / Delta;
  C=(GetXiLimits(u1, u2+Delta) - L) / Delta;

  for m=1:2, n=1:2
    @printf("%i %s: ",Q,Names[m,n]);

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
    @printf("\n");

  end

end

function DoAllSectors()

  u10=0.25; 
  u20=0.5; 

  for s1=-1:2:1, s2=-1:2:1
   u1=s1*u10;
   u2=s2*u20;
   DoSector(u1, u2);
   DoSector(u2, u1);
  end 

end
