reload("GetXiLimits.jl");

##################################################
##################################################
##################################################
function P(Xi1, Xi2, Xi3, Eta1, Eta2, Eta3)
 1.0 + 2.1*Xi1 + 3.2*Xi2 + 0.9*Xi3 + 4.3*Eta1 + 5.4*Eta2 + 6.5*Eta3;
end

##################################################
##################################################
##################################################
function MCIntegrate(N)

  Buckets=zeros(8,6);

LogFile=open("/tmp/MC.dat","w");
  for n=1:N
  
    #################################################################
    #################################################################
    #################################################################
    Xi1=rand();
    Xi2=rand();
    if (Xi2>Xi1) continue; end
    Xi3=rand();
    if (Xi3>Xi2) continue; end
    Eta1=rand();
    Eta2=rand();
    if (Eta2>Eta1) continue; end
    Eta3=rand();
    if (Eta3>Eta2) continue; end
  
    u1=Eta1-Xi1;
    u2=Eta2-Xi2;
    u3=Eta3-Xi3;
  
    (Q1,Q2)=GetSector(u1,u2,u3);
    Buckets[Q1,Q2] += P(Xi1, Xi2, Xi3, Eta1, Eta2, Eta3);

@printf(LogFile,"%i%i %e %e %e %e %e %e %e %e %e\n",Q1,Q2,u1,u2,u3,Xi1,Xi2,Xi3,Eta1,Eta2,Eta3);
  end

  Total=0.0;
  for Q1=1:8, Q2=1:6
    Buckets[Q1,Q2] /= float(N);
    Total += Buckets[Q1, Q2];
    @printf("%i%i %e \n",Q1,Q2,Buckets[Q1,Q2]);
  end
  @printf("Total %e (%e)\n",Total,1.0/36.0);
  
end
