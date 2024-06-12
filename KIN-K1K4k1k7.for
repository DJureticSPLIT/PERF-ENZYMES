       program     KIN
C
C Four-state-kinesin-1-kinetic-scheme (KIN-enzyme)
C Reference Hwang-2016. Table 2
C
C The first row of the output file KIN-K1K4k1k7.dat contains published parameters or parameters calculated
C from measured values. Catalytic efficiency k(cat)/K(M) and overall dissipation/RT parameters from the first
C row are coordinates we used to create the KIN point in Figure 1 (numbers in the row 30, columns A and B
C from our SI dataset (excell file)).
C
C We introduced random noise into forward microscopic rate constants k12 and k41.
C These rate constants (also labeled as k1 and k7) are, respectively, for substrate-enzyme
C association and for product release. Equilibrium constants K1 and K4 are then noisy.
C For fixed force, fixed temperature, and steady state condition the maximum
C in dissipation/RT and optimal catalytic efficiency can be identified in the output file
C KIN-K1K4k1k4.dat. We used that output file to create Figure 2B.
C
C
       real*8  thermale,rckonst1,rckonst2,
     &   xsec,sigma,enprodB1,enprodB2,enprodB3,enprodB4,
     &   k12,k21,k23,k32,k34,k43,k41,k14,ckonst1a(2000),
     &   keqk12,keqk23,keqk34,keqk41,sigmap,sigmab,
     &   sigmah,sigmaq,currB1,currB2,currB3,currB4,
     &   keqk3,keqk1,keqk2,keqk4,enprodT,
     &   sprob(4),entros(4),entrops,aX1,aX2,aX3,aX4,aX,
     &   kcat,KMich,Ktot,ENZeff,products,
     &   s1,s2,pi,gaussrand,substrates
C
C
       keqk3=0.0
       keqk1=0.0
       keqk4=0.0
       keqk2=0.0
       keqk12=0.0
       keqk23=0.0
       keqk34=0.0
       keqk41=0.0
       aX1=0.0
       aX2=0.0
       aX3=0.0
       aX4=0.0
       aX=0.0
       enprodB3=0.0
       enprodB4=0.0
       enprodB1=0.0
       enprodB2=0.0
       sigmap=0.0
       sigmab=0.0
       sigmah=0.0
       sigmaq=0.0
C
       emprodB1=0.0
       emprodB2=0.0
       k12=0.0
       k21=0.0
       k23=0.0
       k32=0.0
       k34=0.0
       k43=0.0
       k41=0.0
       k14=0.0
       substrates=0.0
       products=0.0
       currB1=0.0
       currB2=0.0
       currB3=0.0
       currB4=0.0
       thermale=0.0
       xsec=0.0
       sigma=0.0
       kcat=0.0
       KMich=0.0
       ENZeff=0.0
       Ktot=0.0
       mmm=0
       i=0
       entrops=0.0
C
       s1=0.0
       s2=0.0
       pi=0.0
       gaussrand=0.0
       nnn=0
       jjj=0
C
       do 25 nnn=1,2000
        ckonst1a(nnn)=0.0
25     continue
C
       substrates=0.00005
       products=0.0001
       k12=2300000.0*substrates
       k21=20.0
       k23=600.0
       k32=1.4
       k34=400.0
       k43=1.7
       k41=190.0
       k14=1200000.0*products
       keqk12=k12/k21
       keqk23=428.57143
       keqk34=235.29412
       keqk41=1.58333
C
       Ktot=keqk12*keqk23*keqk34*keqk41
C in moles
         thermale=2.47895703
C RT in kJ/mol for 25 C
         pi=3.141592653589793
C
        open (unit=9,status='unknown',file='KIN-K1K4k1k7.dat',
     & access='sequential',form='formatted')
       write(9,45)
45     format(2x,' step     k1          k2       k3       k4     ',
     & '  k5       k6       k7        k8         Ptot          Eff',
     & '           Ktot           kcat',
     & '           Km         P1           P2    ',
     & '        P3             P4             J     ',
     & '     S(I)          p1              p2       ',
     & '       p3           p4         ',
     & '    K4            K1            X')
C
        do 50 i=1,2000
C
C
        call random_number(s1)
        call random_number(s2)
C
        gaussrand=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1.0
CCCCC    Box-Muller transform for the random normal noise. SHIFT is present here.
C
       substrates=0.00005
       products=0.0001
       k12=2300000.0*substrates
       k21=20.0
       k23=600.0
       k32=1.4
       k34=400.0
       k43=1.7
       k41=190.0
       k14=1200000.0*products
       keqk12=k12/k21
       keqk23=428.57143
       keqk34=235.29412
       keqk41=1.58333
C
       Ktot=keqk12*keqk23*keqk34*keqk41
C
       if (i-1.eq.0) then
           k12=2300000.0*substrates
            else
            if ((gaussrand.gt.0).and.((i-1).gt.0)) then
         ckonst1a(i)=5.75*gaussrand
C
          keqk12=ckonst1a(i)
          k12=20.0*keqk12
C
          keqk41=Ktot/(keqk12*keqk23*keqk34)
C
          k41=k14*keqk41
            endif
         endif
C
       kcat=(k23*k34*k41)/((k23+k32)*(k43+k41)+k23*k34+k34*k41)
C
       Ktot=keqk12*keqk23*keqk34*keqk41
C
       KMich=(k21*(k32*(k43+k41)+k34*k41)+k23*k34*k41)*substrates/
     &       (k12*((k23+k32)*(k43+k41)+k23*k34+k34*k41))
C
        ENZeff=kcat/KMich
C
         sigma=0.0
         sigmab=0.0
         sigmah=0.0
         sigmaq=0.0
         sigmap=0.0
C  p=1,b=2,h=3,q=4
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
	     enprodB3=0.0
         enprodB4=0.0
	     enprodT=0.0
C
           do 55 jj=1,4
         sprob(jj)=0.0
         entros(jj)=0.0
55       continue
C
         sigmap=k43*k32*k21+k21*k34*k41+k41*k32*k21+k23*k34*k41
         sigmab=k34*k41*k12+k14*k43*k32+k12*k43*k32+k32*k41*k12
         sigmah=k41*k12*k23+k21*k14*k43+k14*k43*k23+k12*k23*k43
         sigmaq=k32*k21*k14+k12*k23*k34+k14*k23*k34+k21*k14*k34
C
C
         sigma=sigmab+sigmah+sigmaq+sigmap
C
         sigmab=sigmab/sigma
         sigmah=sigmah/sigma
         sigmaq=sigmaq/sigma
         sigmap=sigmap/sigma
C
         sprob(1)=sigmab
         sprob(2)=sigmah
         sprob(3)=sigmaq
         sprob(4)=sigmap
C
       do 700 mmmm=1,4
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
700    continue
C
C
       xsec=thermale*2.3025851*dlog10((k12*k23*k34*k41)/
     & (k21*k32*k43*k14))
Cn all currents are without N (N=1)
C
        currB1=(k12*sigmap-k21*sigmab)
        currB2=(k23*sigmab-k32*sigmah)
        currB3=(k34*sigmah-k43*sigmaq)
        currB4=(k41*sigmaq-k14*sigmap)
C
          aX1=2.3025851*dlog10((k12*sigmap)/
     &   (k21*sigmab))
C
          aX2=2.3025851*dlog10((k23*sigmab)/
     &   (k32*sigmah))
C
          aX3=2.3025851*dlog10((k34*sigmah)/
     &   (k43*sigmaq))
C
          aX4=2.3025851*dlog10((k41*sigmaq)/
     &   (k14*sigmap))
C
          aX=aX1+aX2+aX3+aX4
C in RT units
C
        enprodB1=aX1*currB1
        enprodB2=aX2*currB2
        enprodB3=aX3*currB3
        enprodB4=aX4*currB4
C
        enprodT=enprodB1+enprodB2+enprodB3+enprodB4
C
C These enprod values are entropy production divided by the gas constant R
C Thus, units are in inverse seconds
C
       write(9,'(i6,1x,f12.5,7(f8.1,1x),4(f15.5),1x,f8.5,13(f13.8,1x))')
     &  i,k12,k21,k23,k32,k34,k43,k41,k14,enprodT,ENZeff,Ktot,kcat,
     &  KMich,enprodB1,enprodB2,enprodB3,enprodB4,currB1,entrops,
     &  sigmap,sigmab,sigmah,sigmaq,keqk41,keqk12,aX
C
 50    continue
C
       close(unit=9)
C
       stop
       end
C Author: Davor Juretic, year:2024
C
C run time around 2 seconds
