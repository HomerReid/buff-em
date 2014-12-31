PointsPerRegion=zeros(48);
TotalPoints=10000000;

f=open("/tmp/TOL3.out","w");
for n=1:TotalPoints

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

 #################################################################
 #################################################################
 #################################################################
 if     (u1>0.0 && u2>0.0 && u3>0.0)
  Q1=0;
 elseif (u1>0.0 && u2>0.0 && u3<0.0)
  Q1=1;
 elseif (u1>0.0 && u2<0.0 && u3>0.0)
  Q1=2;
 elseif (u1>0.0 && u2<0.0 && u3<0.0)
  Q1=3;
 elseif (u1<0.0 && u2>0.0 && u3>0.0)
  Q1=4;
 elseif (u1<0.0 && u2>0.0 && u3<0.0)
  Q1=5;
 elseif (u1<0.0 && u2<0.0 && u3>0.0)
  Q1=6;
 elseif (u1<0.0 && u2<0.0 && u3<0.0)
  Q1=7;
 else
  @printf("Bawonkatage!\n");
 end

 if     (abs(u1)>abs(u2) && abs(u2)>abs(u3))
  Q2=0;
 elseif (abs(u1)>abs(u3) && abs(u3)>abs(u2))
  Q2=1;
 elseif (abs(u2)>abs(u1) && abs(u1)>abs(u3))
  Q2=2;
 elseif (abs(u2)>abs(u3) && abs(u3)>abs(u1))
  Q2=3;
 elseif (abs(u3)>abs(u1) && abs(u1)>abs(u2))
  Q2=4;
 elseif (abs(u3)>abs(u2) && abs(u2)>abs(u1))
  Q2=5;
 else
  @printf("Bawonkatage!\n");
 end

 Region = Q1*6 + Q2;
 PointsPerRegion[Region+1] += 1;
 @printf(f,"%i %i %e %e %e %e %e %e %e %e %e\n",
            Q1,Q2,u1,u2,u3,Xi1,Xi2,Xi3,Eta1,Eta2,Eta3);
end
close(f);

#################################################################
#################################################################
#################################################################
TotalVolume=0.0;
NumNonZero=0;
for Q1=0:7, Q2=0:5
  Region = Q1*6 + Q2;
  Volume = PointsPerRegion[Region+1] / TotalPoints;
  TotalVolume+=Volume;
  if Volume>0.0
    @printf("(%i,%i): %e\n",Q1,Q2,Volume);
    NumNonZero+=1;
  end
end
@printf("%i nonempty regions.\n",NumNonZero);
@printf("Total volume: %e (should be %e)\n",TotalVolume,1.0/36.0);
