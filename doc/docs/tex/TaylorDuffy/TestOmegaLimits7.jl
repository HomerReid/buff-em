N=3

for nXi1=0:N, nXi2=0:N, nXi3=0:N, nEta1=0:N, nEta2=0:N, nEta3=0:N

  Xi1=(float(nXi1)/float(N));
  Xi2=(float(nXi2)/float(N)) * Xi1;
  Xi3=(float(nXi2)/float(N)) * Xi2;

  Eta1=(float(nEta1)/float(N));
  Eta2=(float(nEta2)/float(N)) * Eta1;
  Eta3=(float(nEta2)/float(N)) * Eta2;

  u1=Eta1-Xi1;
  u2=Eta2-Xi2;
  u3=Eta3-Xi3;

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

  if     (u1>u2 && u2>u3)
    Q2=1;
  elseif (u1>u3 && u3>u2)
    Q2=2;
  elseif (u2>u1 && u1>u3)
    Q2=2;
  elseif (u2>u3 && u3>u1)
    Q2=3;
  elseif (u3>u1 && u1>u2)
    Q2=4;
  elseif (u3>u2 && u2>u1)
    Q2=5;
  end

  Xi1Max[Q1, Q2]=

end
