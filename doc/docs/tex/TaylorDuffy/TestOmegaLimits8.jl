function XiLimits(u1, u2, u3, N=1000000) 
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

  #  for ne1=0:1, ne2=0:1, ne3=0:1

  #  Eta1 =       float(ne1) / 1.0;
  #  Eta2 = Eta1*(float(ne2) / 1.0);
  #  Eta3 = Eta2*(float(ne3) / 1.0);

    Xi1 = Eta1 - u1;
    Xi2 = Eta2 - u2;
    Xi3 = Eta3 - u3;
#@printf("%e %e %e %e %e %e \n",Eta1,Eta2,Eta3,Xi1,Xi2,Xi3);

    if (Xi1<0.0 || Xi1>1.0) 
#     @printf("Bawonk 1 %e \n",Xi1);
     continue;
    end

    if (Xi2<0.0 || Xi2>Xi1) 
#     @printf("Bawonk 2 %e %e \n",Xi2,Xi1);
     continue;
    end

    if (Xi3<0.0 || Xi3>Xi2) 
#     @printf("Bawonk 3 %e %e \n",Xi3,Xi2);
     continue;
    end

    Xi1Max = max(Xi1, Xi1Max);
    Xi2Max = max(Xi2, Xi2Max);
    Xi3Max = max(Xi3, Xi3Max);

    Xi1Min = min(Xi1, Xi1Min);
    Xi2Min = min(Xi2, Xi2Min);
    Xi3Min = min(Xi3, Xi3Min);
  end

#  Xi1MaxShouldBe = 1.0 + max(0.0, u1)
#  XiMinShouldBe=max(0,

  [Xi1Min Xi1Max; 
   Xi2Min Xi2Max; 
   Xi3Min Xi3Max]

end
