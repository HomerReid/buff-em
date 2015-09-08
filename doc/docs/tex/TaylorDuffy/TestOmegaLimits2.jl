Points0=zeros(1,9);

Points1=zeros(1,9);
Points2=zeros(1,9);
Points3=zeros(1,9);
Points4=zeros(1,9);
Points5=zeros(1,9);
Points6=zeros(1,9);

Points7=zeros(1,9);
Points8=zeros(1,9);
Points9=zeros(1,9);
Points10=zeros(1,9);
Points11=zeros(1,9);
Points12=zeros(1,9);

for Xi1=0.0:0.5:1, Xi2Mult=0.0:0.5:1, Xi3Mult=0.0:0.5:1, Eta1=0.0:0.5:1, Eta2Mult=0.0:0.5:1, Eta3Mult=0.0:0.5:1

 Xi2=Xi2Mult*Xi1;
 Xi3=Xi3Mult*Xi2;
 Eta2=Eta2Mult*Eta1;
 Eta3=Eta3Mult*Eta2;

 u1=Eta1-Xi1;
 u2=Eta2-Xi2;
 u3=Eta3-Xi3;

 Classified=0;
  
 if (Xi1==0.0 && Xi2==0.0 && Xi3==0.0 && Eta1==0.0 && Eta2==0.0 && Eta3==0.0)
   continue
 end

 Points0=[Points0; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3];

 if (u3>=0.0)
   if     (u1>=0.0 && u2>=0.0 && u2>=u1)
     Points1=[Points1; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1>=0.0 && u2>=0.0 && u2<=u1)
     Points2=[Points2; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1>=0.0 && u2<=0.0)
     Points3=[Points3; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2<=0.0 && u2<=u1)
     Points4=[Points4; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2<=0.0 && u2>=u1)
     Points5=[Points5; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2>=0.0)
     Points6=[Points6; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end
 end
 
 if (u3<=0.0)
   if     (u1>=0.0 && u2>=0.0 && u2>=u1)
     Points7=[Points7; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1>=0.0 && u2>=0.0 && u2<=u1)
     Points8=[Points8; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1>=0.0 && u2<=0.0)
     Points9=[Points9; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2<=0.0 && u2<=u1)
     Points10=[Points10; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2<=0.0 && u2>=u1)
     Points11=[Points11; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

   if     (u1<=0.0 && u2>=0.0)
     Points12=[Points12; u1 u2 u3 Xi1 Xi2 Xi3 Eta1 Eta2 Eta3]
     Classified=1;
   end

 end

   if (Classified==0)
     @printf("Bawonkatage! %e %e %e \n",u1,u2,u3);
   end

end
