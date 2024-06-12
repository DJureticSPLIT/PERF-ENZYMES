       program     GI
C
C Two-state-glucose-isomerase-enzyme-kinetic-scheme
C
CCCC  Values of rate constants are from Converti-1998 paper.
C     Second order rate constants k1s and k4s
C     are expressed in inverse moles and inverse seconds.
C
C
C The first row of the output file GI-K1K2k1k3.dat contains published parameters or parameters calculated
C from measured values. Catalytic efficiency k(cat)/K(M) and overall dissipation/RT parameters from the first
C row are coordinates we used to create the GI point in Figure 1 (numbers in the row 2, columns A and B
C from our SI dataset (excell file).
C
C We introduced random noise into forward microscopic rate constants k1 and k3.
C These rate constants are, respectively, for substrate-enzyme
C association and for product release. Equilibrium constants K1 and K2 are then noisy.
C For fixed force, fixed temperature, and steady state condition the maximum
C in dissipation/RT and optimal catalytic efficiency can be identified in the output file GI-K1K2k1k3.
C We used that output file to create Figure 2A
C
       real*8  thermale,k1,k2,k3,k4,k1s,k4s,
     &   xsec,sigma,enprodB1,enprodB2,aX1,aX2,aXtot,
     &   ckonst1a(10000),gaussrand1, currB1,currB2,
     &   keqk12L,keqk21D,sigmap,sigmab,enprodT,
     &   sprob(2),entros(2),entrops,substrates,products,
     &   kcat,Km,Ktot,enzeff,s1,s2,pi
C
C
       aX1=0.0
       aX2=0.0
C
       aXtot=0.0
C
       enprodB1=0.0
       enprodB2=0.0
C
       sigmap=0.0
       sigmab=0.0
C
C
       k1=0.0
       k2=0.0
       k3=0.0
       k4=0.0
       k1s=0.0
       k4s=0.0
C
       keqk12L=0.0
       keqk21D=0.0
C
       currB1=0.0
       currB2=0.0
C
       thermale=0.0
       xsec=0.0
       sigma=0.0
       kcat=0.0
       Km=0.0
       enzeff=0.0
       Ktot=0.0
       mmm=0
       i=0
       entrops=0.0
C
       substrates=0.0
       products=0.0
C
       s1=0.0
       s2=0.0
       pi=0.0
C
       gaussrand1=0.0
C
       nnn=0
       jjj=0
C
       do 25 nnn=1,10000
        ckonst1a(nnn)=0.0
25     continue
C
       substrates=0.5
C
        products=0.05
C
         thermale=2.8115355
C Corresponds to 65 C. Units are in kJ/mol
C
        k1s=0.063
C
        k1=k1s*substrates
C
        k2=0.021
C
        k3=0.029
C
        k4s=0.082
C
        k4=k4s*products
C
         kcat=k3
C
       keqk12L=k1/k2
C
       keqk21D=k3/k4
C
       Ktot=keqk12L*keqk21D
C
         Km=(k2+k3)/k1s
         enzeff=kcat/Km
C
         sigma=0.0
         sigmab=0.0
         sigmap=0.0
C 2-th state probability  p=1,b=2
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
C
	     enprodT=0.0
C
           do 75 jj=1,2
         sprob(jj)=0.0
         entros(jj)=0.0
75       continue
C
C
            open (unit=9,status='unknown',file='GI-K1K2k1k3.dat',
     & access='sequential',form='formatted')
       write(9,45)
45      format(2x,'   i     k1      k2       k3    ',
     &  '   k4     S(I)    p1    p2        P1   ',
     &  '    P2      Ptot    enzeff     kcat  ',
     &  '    Km         J        xsec  ',
     &  '    Xtot       X1         X2     Ktot   ',
     &  '      S         P       S+P')
C
C
        substrates=0.5
C
        products=0.05
C
C
        k1s=0.063
C
        k1=k1s*substrates
C
        k2=0.021
C
        k3=0.029
C
        k4s=0.082
C
        k4=k4s*products
C
C
         keqk12L=k1/k2
C
         keqk21D=k3/k4
C
         Ktot=keqk12L*keqk21D
C
         thermale=2.8115355
C Corresponds to 65 C. Units are in kJ/mol
C
         pi=3.141592653589793
C
        do 50 i=1,10000
C
        call random_number(s1)
        call random_number(s2)
C
        gaussrand1=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1
C
        if (i-1.eq.0) then
          k1s=0.063
          k1=k1s*substrates
            else
C
            if ((gaussrand1.gt.0).and.((i-1).gt.0)) then
C
         ckonst1a(i)=1.5*gaussrand1
C
          keqk12L=ckonst1a(i)
          k1=k2*keqk12L
          k1s=k1/substrates
C
          keqk21D=Ktot/keqk12L
C
          k3=k4*keqk21D
            endif
         endif
C
C
        kcat=k3
        Km=(k2+k3)/k1s
        enzeff=kcat/Km
C
         sigma=0.0
         sigmab=0.0
         sigmap=0.0
C 2-th state probability  p=1,b=2
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
C
	     enprodT=0.0
C
           do 55 jj=1,2
         sprob(jj)=0.0
         entros(jj)=0.0
55       continue
C
         sigmap=k2+k3
C
         sigmab=k1+k4
C
         sigma=sigmap+sigmab
C
         sigmab=sigmab/sigma
C
         sigmap=sigmap/sigma
C
         sigma=sigmap+sigmab
C
         sprob(1)=sigmab
C
         sprob(2)=sigmap
C
       do 700 mmmm=1,2
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
700    continue
C
C
         xsec=thermale*2.3025851*dlog10((k1*k3)/(k2*k4))
Cn all currents are without N (N=1)
Cn
        currB1=(k1*sigmap-k2*sigmab)
        currB2=(k3*sigmab-k4*sigmap)
C
          aX1=2.3025851*dlog10((k1*sigmap)/
     &   (k2*sigmab))
Cn
          aX2=2.3025851*dlog10((k3*sigmab)/
     &   (k4*sigmap))
Cn
C
          aXtot=aX1+aX2
C
Cn all "entropy productions" are TP expressions in RT units
Cn
        enprodB1=aX1*currB1
        enprodB2=aX2*currB2
C
        enprodT=enprodB1+enprodB2
C
C
        write(9,'(i6,4(f8.4,1x),3(f6.4,1x),4(f9.4),11(f10.4))')
     &  i,k1,k2,k3,k4,entrops,sigmap,sigmab,enprodB1,enprodB2,
     &  enprodT,enzeff,kcat,Km,currB1,xsec,aXtot,aX1,aX2,Ktot,
     &  substrates,products,substrates+products
C
 50    continue
C
       close(unit=9)
C
       stop
       end
C Author: Davor Juretic, year:2024
C
C The run time should be around 2 sec.
C
