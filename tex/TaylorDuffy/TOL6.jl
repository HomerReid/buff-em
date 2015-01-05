N1=100;
N2=10;
N3=10;

function F(u1,u2,Q,f)

  Xi1Lower = max(0.0,    -u1, -u2, u2-u1);
  Xi1Upper = min(1.0, 1.0-u1);
if (Xi1Upper<=Xi1Lower) 
  return
end

  for n=1:N2
    Xi1 = Xi1Lower + (Xi1Upper-Xi1Lower)*rand();

    Xi2Lower = max(0.0,  -u2);
    Xi2Upper = Xi1 - max(0.0, u2-u1);

if (Xi2Upper<=Xi2Lower) 
  return
end
    for m=1:N3
       Xi2 = Xi2Lower + (Xi2Upper-Xi2Lower)*rand()
       @printf(f,"%i %e %e %e %e %e %e \n",Q,u1,u2,Xi1,Xi2,Xi1+u1,Xi2+u2);
    end
  end
end

f=open("/tmp/TOL6.out","w");
for n=1:N1

  u1 = rand();
  u2 = u1*rand();

  F(u1,u2,11,f)
  F(u2,u1,12,f)

  F(-u1,u2,21,f)
  F(-u2,u1,22,f)

  F(-u1,-u2,31,f)
  F(-u2,-u1,32,f)

  F(u1,-u2,41,f)
  F(u2,-u1,42,f)

end

close(f)
