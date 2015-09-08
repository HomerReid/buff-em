f=open("BoundaryLines.out","w");
for x=0:0.01:1, BoundaryFlags=0:31

    Xi1  = ((BoundaryFlags &  1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags &  2)==0) ? 0 : Xi1;
    Xi3  = ((BoundaryFlags &  4)==0) ? 0 : Xi2;
    Eta1 = ((BoundaryFlags &  8)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags & 16)==0) ? 0 : Eta1;
    Eta3 = x;
    if (Eta3<=Eta2)
      @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
               Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
               Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 
    end

    Xi1  = ((BoundaryFlags &  1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags &  2)==0) ? 0 : Xi1;
    Xi3  = ((BoundaryFlags &  4)==0) ? 0 : Xi2;
    Eta1 = ((BoundaryFlags &  8)==0) ? 0 : 1;
    Eta2 = x;
    Eta3 = ((BoundaryFlags & 16)==0) ? 0 : Eta2;
    if (Eta2<=Eta1)
      @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
               Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
               Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 
    end

    Xi1  = ((BoundaryFlags &  1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags &  2)==0) ? 0 : Xi1;
    Xi3  = ((BoundaryFlags &  4)==0) ? 0 : Xi2;
    Eta1 = x;
    Eta2 = ((BoundaryFlags &  8)==0) ? 0 : Eta1;
    Eta3 = ((BoundaryFlags & 16)==0) ? 0 : Eta2;
    @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
            Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
            Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 

    Xi1  = ((BoundaryFlags &  1)==0) ? 0 : 1;
    Xi2  = ((BoundaryFlags &  2)==0) ? 0 : Xi1;
    Xi3  = x
    Eta1 = ((BoundaryFlags &  4)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags &  8)==0) ? 0 : Eta1;
    Eta3 = ((BoundaryFlags & 16)==0) ? 0 : Eta2;
    if (Xi3<=Xi2)
      @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
               Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
               Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 
    end

    Xi1  = ((BoundaryFlags &  1)==0) ? 0 : 1;
    Xi2  = x
    Xi3  = ((BoundaryFlags &  2)==0) ? 0 : Xi2;
    Eta1 = ((BoundaryFlags &  4)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags &  8)==0) ? 0 : Eta1;
    Eta3 = ((BoundaryFlags & 16)==0) ? 0 : Eta2;
    if (Xi2<=Xi1)
      @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
               Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
               Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 
    end

    Xi1  = x
    Xi2  = ((BoundaryFlags &  1)==0) ? 0 : Xi1;
    Xi3  = ((BoundaryFlags &  2)==0) ? 0 : Xi2;
    Eta1 = ((BoundaryFlags &  4)==0) ? 0 : 1;
    Eta2 = ((BoundaryFlags &  8)==0) ? 0 : Eta1;
    Eta3 = ((BoundaryFlags & 16)==0) ? 0 : Eta2;
    @printf(f,"%i%i %e %e %e %e %e %e %e %e %e \n",BoundaryFlags,4,
            Eta1-Xi1, Eta2-Xi2, Eta3-Xi3, 
            Xi1, Xi2, Xi3, Eta1, Eta2, Eta3); 

end
