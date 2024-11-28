
       program     RTEM
C
C  Three-state reversible cycle kinetic model.
C  Experimental values are used for beta-lactamase RTEM at 0.12s, as in
C  Christensen's 1990 paper with remaining kinetic constants derived from
C  kinetic modeling in Juretic et al 2019 lactamase paper.
C
C
C The output of the program RTEM-SimII.for are RTEM-SIexper.dat and RTEM-SimII.dat.
C The first row with numbers from the RTEM-SimII.dat must give the same parameters
C already contained in the single row with numbers from the RTEM-SIexper.dat.
C Remaining rows (10000 rows) serve to find optimal catalytic efficiency
C and optimal catalytic constant from maximal dissipation.
C Maximal dissipation follows from stochastic variations in enzyme-substrate association
C rate constant k1=ckonst1 and enzyme-product dissociation rate constant k5=ck3
C when overall force X/RT does not change.
C
       real*8  ckonst1a(10000),
     &   thermale,rckonst1,xsec,sigma,enprodB1,rckonst2,
     &   enprodB2,ckonst1,ckonst2,ck3,rck3,sigmap,sigmab,
     &   sigmah,sigmaq,currB1,currB2,currB3,enprodB3,
     &   abh,ahq,aqp,ahp,aXtot,substrates,products,
     &   keqk3,keqk1,keqk2,Ktot, s1,s2,pi,gaussrand1,
     &   kcat,Kmm,specificity,sprob(3),entros(3),entrops,Ptot
C
       s1=0.0
       s2=0.0
       pi=0.0
       gaussrand1=0.0
       keqk3=0.0
       keqk1=0.0
       keqk2=0.0
       Ktot=0.0
       ahp=0.0
       aqp=0.0
       ahq=0.0
       aXtot=0.0
       enprodB3=0.0
       enprodB1=0.0
       enprodB2=0.0
       enprodT=0.0
       Ptot=0.0
	   kcat=0.0
	   Kmm=0.0
	   specificity=0.0
	   substrates=0.0
       products=0.0
C
       sigmap=0.0
       sigmab=0.0
       sigmah=0.0
       sigmaq=0.0
C
       ckonst1=0.0
       ckonst2=0.0
       currB1=0.0
       currB2=0.0
       currB3=0.0
C
       thermale=0.0
C
       rckonst1=0.0
       rckonst2=0.0
       ck3=0.0
       rck3=0.0
C
       xsec=0.0

       sigma=0.0
C
       mmm=0
       i=0
       entrops=0.0
C
       nnn=0
       jjj=0
       do 25 j=1,10000
       ckonst1a(j)=0.0
25     continue
C
C       substrate and product conc in moles for
C         RTEM beta lactamase at 0.12s
         substrates=0.001390
	     products=0.000110
C
	   ckonst1=123*1000000.0*substrates
       ckonst2=2800.0
       ck3=1500.0
       rckonst1=11800.0
       rckonst2=6.0
       rck3=40*1000000.0*products
C
	     keqk1=14.489
	     keqk2=466.667
	     keqk3=0.3409
C
           thermale=2.4361
C
        open(unit=15,status='unknown',file='RTEM-SIexper.dat',
     & access='sequential',form='formatted')
	    write(15,275)
275      format(2x,'         k1              k2                   k3 ',
     &  '                  k4               k5                 k6     ',
     &  '          Ptot              Ptot             Eff            ',
     &  '           kcat               Km               J   ',
     &  '                   P1                 P2                 P3 ',
     &  '                 X               X/RT                S(i)  ',
     &  '              p1                p2                  p3')
C
C
        kcat=(ckonst2*ck3)/(ckonst2+ck3+rckonst2)
	     Kmm=(ckonst2*ck3+rckonst1*rckonst2+rckonst1*ck3)*substrates/
     &       (ckonst1*(ckonst2+rckonst2+ck3))
	     specificity=kcat/Kmm
C
         sigma=0.0
         sigmab=0.0
         sigmah=0.0
         sigmaq=0.0
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
	     enprodB3=0.0
	     enprodT=0.0
         do 45 jj=1,3
         sprob(jj)=0.0
         entros(jj)=0.0
45       continue
C
       sigmab=rckonst2*rckonst1+ckonst2*ck3+rckonst1*ck3
Cn
       sigmah=ckonst1*rckonst2+rck3*rckonst2+ck3*ckonst1
Cn
       sigmaq= rckonst1*rck3+ckonst1*ckonst2+rck3*ckonst2
Cn
         sigma=sigmab+sigmah+sigmaq
C
         sigmab=sigmab/sigma
         sigmah=sigmah/sigma
         sigmaq=sigmaq/sigma
C
         sprob(1)=sigmab
         sprob(2)=sigmah
         sprob(3)=sigmaq
C
       do 770 mmmm=1,3
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
770    continue
C
       xsec=thermale*2.3025851*dlog10(ckonst1*ckonst2*ck3/
     & (rckonst1*rckonst2*rck3))
Cn all currents are without N (N=1)
Cn
        currB1=(ckonst1*sigmab-rckonst1*sigmah)
        currB2=(ckonst2*sigmah-rckonst2*sigmaq)
        currB3=(ck3*sigmaq-rck3*sigmab)
Cn
        abh=thermale*2.3025851*dlog10((ckonst1*sigmab)/
     &  (rckonst1*sigmah))
Cn
         ahq=thermale*2.3025851*dlog10((ckonst2*sigmah)/
     &  (rckonst2*sigmaq))
Cn
         aqp=thermale*2.3025851*dlog10((ck3*sigmaq)/
     &  (rck3*sigmab))
C
          aXtot=(abh+ahq+aqp)/thermale
C
Cn all "(entropy productions/R" are in units of inverse seconds
Cn
        enprodB1=currB1*2.3025851*dlog10((ckonst1*sigmab)/
     &   (rckonst1*sigmah))
        enprodB2=currB2*2.3025851*dlog10((ckonst2*sigmah)/
     &   (rckonst2*sigmaq))
        enprodB3=currB3*2.3025851*dlog10((ck3*sigmaq)/
     &   (rck3*sigmab))
C
C
         enprodT=currB2*2.3025851*dlog10((ckonst1*ckonst2*ck3)/
     &   (rckonst1*rckonst2*rck3))
C entropy production divided by R is in units of inverse seconds
C
         Ptot=enprodB1+enprodB2+enprodB3
C

	     write(15,'(21(f18.6,1x))')
     &  ckonst1,rckonst1,ckonst2,rckonst2,ck3,rck3,Ptot,enprodT,
     &  specificity,kcat,Kmm,currB3,enprodB1,enprodB2,enprodB3,
     &  xsec,aXtot,entrops,sigmab,sigmah,sigmaq
         close(unit=15)
C
        open (unit=7,status='unknown',file='RTEM-SimII.dat',
     & access='sequential',form='formatted')
       write(7,40)
40     format(2x,'   i         k1          k2            k3      ',
     &  '        k4           k5           k6           S(I) ',
     &  '             P1             P2              P3 ',
     &  '            Ptot             J                   Km',
     &  '                kcat               Ptot          Eff        ',
     &  '               p1                 p2                 p3',
     &  '                K1                  K2                 K3   ',
     &  '             X/RT             Ktot ')
C
       keqk3=0.0
       keqk1=0.0
       keqk2=0.0
       kcat=0.0
       Kmm=0.0
       specificity=0.0
C
          substrates=0.001390
	     products=0.000110
C
	     ckonst1=123*1000000.0*substrates
       ckonst2=2800.0
       ck3=1500.0
       rckonst1=11800.0
       rckonst2=6.0
       rck3=40*1000000.0*products
C
           thermale=2.4361
C for 20 C in kJ pr mol

	     keqk1=ckonst1/rckonst1
CC
	     keqk2=ckonst2/rckonst2
CC
	     keqk3=ck3/rck3
C
          Ktot=keqk1*keqk2*keqk3
C
         pi=3.141592653589793
C
        do 50 i=1,10000
C
        call random_number(s1)
        call random_number(s2)
C
        gaussrand1=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1.0
C
           if (i-1.eq.0) then
         ckonst1=123*1000000.0*substrates
C
            else
            if ((gaussrand1.gt.0).and.((i-1).gt.0)) then
         ckonst1a(i)=14.489*gaussrand1
C
          keqk1=ckonst1a(i)
          ckonst1=11800.0*keqk1
C
          keqk3=Ktot/(keqk1*keqk2)
C
          ck3=rck3*keqk3
            endif
         endif
C
C calculate k(cat), K(M), and k(cat)/K(M)
C
         kcat=(ckonst2*ck3)/(ckonst2+ck3+rckonst2)
	     Kmm=(ckonst2*ck3+rckonst1*rckonst2+rckonst1*ck3)*substrates/
     &       (ckonst1*(ckonst2+rckonst2+ck3))
	     specificity=kcat/Kmm
C
C
         sigma=0.0
         sigmab=0.0
         sigmah=0.0
         sigmaq=0.0
         entrops=0.0
C
         do 55 jj=1,3
         sprob(jj)=0.0
         entros(jj)=0.0
55       continue
C
       sigmab=rckonst2*rckonst1+ckonst2*ck3+rckonst1*ck3
Cn
       sigmah=ckonst1*rckonst2+rck3*rckonst2+ck3*ckonst1
Cn
       sigmaq= rckonst1*rck3+ckonst1*ckonst2+rck3*ckonst2
Cn
         sigma=sigmab+sigmah+sigmaq
C
         sigmab=sigmab/sigma
         sigmah=sigmah/sigma
         sigmaq=sigmaq/sigma
C
         sprob(1)=sigmab
         sprob(2)=sigmah
         sprob(3)=sigmaq
C
       do 700 mmmm=1,3
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
700    continue
C
       xsec=thermale*2.3025851*dlog10(ckonst1*ckonst2*ck3/
     & (rckonst1*rckonst2*rck3))
Cn all currents are without N (N=1)
Cn
        currB1=(ckonst1*sigmab-rckonst1*sigmah)
        currB2=(ckonst2*sigmah-rckonst2*sigmaq)
        currB3=(ck3*sigmaq-rck3*sigmab)
Cn
         abh=thermale*2.3025851*dlog10((ckonst1*sigmab)/
     &  (rckonst1*sigmah))
Cn
         ahq=thermale*2.3025851*dlog10((ckonst2*sigmah)/
     &  (rckonst2*sigmaq))
Cn
         aqp=thermale*2.3025851*dlog10((ck3*sigmaq)/
     &  (rck3*sigmab))
C
         aXtot=(abh+ahq+aqp)/thermale
C
Cn all "(entropy productions/R" are in units of inverse seconds
Cn
        enprodB1=currB1*2.3025851*dlog10((ckonst1*sigmab)/
     &   (rckonst1*sigmah))
        enprodB2=currB2*2.3025851*dlog10((ckonst2*sigmah)/
     &   (rckonst2*sigmaq))
        enprodB3=currB3*2.3025851*dlog10((ck3*sigmaq)/
     &   (rck3*sigmab))
C
C
         enprodT=currB2*2.3025851*dlog10((ckonst1*ckonst2*ck3)/
     &   (rckonst1*rckonst2*rck3))
C entropy production divided by R is in units of inverse seconds
C
         Ptot=enprodB1+enprodB2+enprodB3
C
       write(7,'(i5,1x,6(f12.1,2x),(f10.7,1x),5(f17.5),12(f18.6,1x))')
     &  i,ckonst1,rckonst1,ckonst2,rckonst2,ck3,rck3,entrops,
     &  enprodB1,enprodB2,enprodB3,enprodT,currB3,Kmm,kcat,
     &  Ptot,specificity,sigmab,sigmah,sigmaq,keqk1,keqk2,keqk3,
     &  aXtot,Ktot
C
 50    continue
C
215     continue
C
250     continue
C
          close(unit=7)
C
       stop
       end
C Author: Davor Juretic, 2024
C Run time about 2 sec

