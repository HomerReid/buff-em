#################################################
#################################################
#################################################
function GetXiLimits(u1, u2, u3, f, N=100000) 
  Xi1Max=0;
  Xi2Max=0;
  Xi3Max=0;

  Xi1Min=1;
  Xi2Min=1;
  Xi3Min=1;

#  for n=1:N
#    Eta1=rand();
#    Eta2=Eta1*rand();
#    Eta3=Eta2*rand();
N=100
  for nx1=0:N, nx2=0:N, nx3=0:N

    Eta1=(nx1)/float(N);
    Eta2=Eta1*(nx2)/float(N);
    Eta3=Eta2*(nx3)/float(N);

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

    @printf(f,"%e %e %e\n",Xi1,Xi2,Xi3); 

  end


end

#################################################
#################################################
#################################################
function GetXiLimits2(u1, u2, u3, f)
  
  Xi1Upper=1 - maximum([ 0,u1 ]);
  Xi1Lower=maximum([0,-u1,-u2,-u3,u2-u1,u3-u1,u3-u2,u2-u3-u1]);
  DXi2Upper=maximum([0,-u1,u2-u1]);
  Xi2Lower=maximum([0,-u2,-u3,u3-u2]);
  DXi3Upper=maximum([0,-u2,u3-u2]);
  Xi3Lower=maximum([0,-u3]);
 
  N=10
  for nx1=0:N, nx2=0:N, nx3=0:N

    dx1=float(nx1)/float(N);
    dx2=float(nx2)/float(N);
    dx3=float(nx3)/float(N);
  
    Xi1 = Xi1Lower + dx1*(Xi1Upper-Xi1Lower); 

    Xi2Upper = Xi1 - DXi2Upper;
    Xi2 = Xi2Lower + dx2*(Xi2Upper-Xi2Lower);

    Xi3Upper = Xi2 - DXi3Upper;
    Xi3 = Xi3Lower + dx3*(Xi3Upper-Xi3Lower);

    if (Xi2>Xi1 || Xi3>Xi2) 
     continue; 
    end

    @printf(f,"%e %e %e \n",Xi1,Xi2,Xi3);

  end

end

#######################################################
#######################################################
#######################################################
function Compare(u1,u2,u3)

  f=open("/tmp/XiLimits1.out","w");
  GetXiLimits(u1, u2, u3, f);
  close(f);

  f=open("/tmp/XiLimits2.out","w");
  GetXiLimits2(u1, u2, u3, f);
  close(f);

end
