f=open("BoundaryLines.out","w");
for x=0:0.01:1, BoundaryFlags=0:7

    Xi1  = ((BoundaryFlags & 1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags & 2)==0) ? 0 : Xi1;
    Eta1 = ((BoundaryFlags & 4)==0) ? 0 : 1;
    Eta2 = x;
    if (Eta2<=Eta1)
      @printf(f,"%i%i %e %e %e %e %e %e\n",BoundaryFlags,4,
               Eta1-Xi1, Eta2-Xi2, Xi1, Xi2, Eta1, Eta2); 
    end

    Xi1  = ((BoundaryFlags & 1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags & 2)==0) ? 0 : Xi1;
    Eta1 = x;
    Eta2 = ((BoundaryFlags & 4)==0) ? 0 : Eta1;
    @printf(f,"%i%i %e %e %e %e %e %e\n",BoundaryFlags,3,
            Eta1-Xi1, Eta2-Xi2, Xi1, Xi2, Eta1, Eta2); 

    Xi1  = ((BoundaryFlags & 1)==0) ? 0 : 1;
    Xi2  = x;
    Eta1 = ((BoundaryFlags & 2)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags & 4)==0) ? 0 : Eta1;
    if (Xi2<=Xi1)
      @printf(f,"%i%i %e %e %e %e %e %e\n",BoundaryFlags,2,
               Eta1-Xi1, Eta2-Xi2, Xi1, Xi2, Eta1, Eta2); 
    end

    Xi1  = x;
    Xi2  = ((BoundaryFlags & 1)==0) ? 0 : Xi1
    Eta1 = ((BoundaryFlags & 2)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags & 4)==0) ? 0 : Eta1;
    @printf(f,"%i%i %e %e %e %e %e %e\n",BoundaryFlags,1,
             Eta1-Xi1, Eta2-Xi2, Xi1, Xi2, Eta1, Eta2); 

end
close(f)
