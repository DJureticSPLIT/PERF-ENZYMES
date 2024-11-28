
       program     PC1
C
C  Experimental values are used for beta-lactamase PC1 at 0.12s, as in
C  Christensen's 1990 paper with remaining kinetic parameters derived from
C  kinetic modeling in Juretic et al. 2019 lactamase paper.
C
C
C The output of the program PC1-SimII.for are PC1-SIexper.dat and PC1-SimII.dat.
C The first row with numbers from the PC1-SimII.dat must give the same parameters
C already contained in the single row with numbers from the PC1-SIexper.dat.
C Remaining rows (10000 rows) serve to find optimal catalytic efficiency
C and optimal catalytic constant from maximal dissipation.
C Maximal dissipation follows from stochastic variations in enzyme-substrate association
C rate constant k1=ckonst1 and enzyme-product dissociation rate constant k5=ck3
C when overall force X/RT does not change.
C
       real*8  ckonst1a(10000),ckonst2a(10000),
     &   thermale,rckonst1,xsec,sigma,enprodB1,rckonst2,
     &   enprodB2,ckonst1,ckonst2,ck3,rck3,sigmap,sigmab,
     &   sigmah,sigmaq,currB1,currB2,currB3,enprodB3,
     &   abh,ahq,aqp,ahp,aXtot,substrates,products,
     &   ckonst3a(10000),ckonst4a(10000),ckonst5a(10000),
     &   ckonst6a(10000),keqk3,keqk1,keqk2,Ktot,
     &   kcat,Kmm,specificity,sprob(3),entros(3),entrops,Ptot,
     &   s1,s2,pi,gaussrand1,gaussrand2,gaussrand3,gaussrand4,
     &   gaussrand5,gaussrand6
C
       s1=0.0
       s2=0.0
       pi=0.0
       gaussrand1=0.0
       gaussrand2=0.0
       gaussrand3=0.0
       gaussrand4=0.0
       gaussrand5=0.0
       gaussrand6=0.0
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
       ckonst2a(j)=0.0
       ckonst3a(j)=0.0
       ckonst4a(j)=0.0
       ckonst6a(j)=0.0
       ckonst5a(j)=0.0
C
25     continue
C
C       substrate and product conc in moles for
C         beta lactamase I at 0.12s
         substrates=0.001492
	     products=0.000008
C
	     ckonst1=22*1000000.0*substrates
        ckonst2=173.0
       ck3=96.0
       rckonst1=196.0
       rckonst2=4.0
       rck3=1000000.0*products
C
	     keqk1=167.4694
	     keqk2=43.25
	     keqk3=12.0
C
           thermale=2.4361
C
        open(unit=15,status='unknown',file='PC1-SIexper.dat',
     & access='sequential',form='formatted')
	    write(15,275)
275      format(2x,'       k1            k2             k3         ',
     &  '        k4              k5              k6             Ptot',
     &  '            enprT             Eff             kcat',
     &  '              Km              J              P1 ',
     &  '            P2               P3               X ',
     &  '            X/RT            S(I)             p1 ',
     &  '             p2               p3')

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
C
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
C
C
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

	     write(15,'(21(f15.6,1x))')
     &  ckonst1,rckonst1,ckonst2,rckonst2,ck3,rck3,Ptot,enprodT,
     &  specificity,kcat,Kmm,currB3,enprodB1,enprodB2,enprodB3,
     &  xsec,aXtot,entrops,sigmab,sigmah,sigmaq
         close(unit=15)
C
        open (unit=7,status='unknown',file='PC1-SimII.dat',
     & access='sequential',form='formatted')
       write(7,40)
40     format(2x,'  i      k1            k2          k3           k4  ',
     &  '       k5           k6       S(I)        P1          P2      ',
     &  '   P3           Ptot             J               Km',
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
        substrates=0.001492
	     products=0.000008
C
	     ckonst1=22*1000000.0*substrates
        ckonst2=173.0
       ck3=96.0
       rckonst1=196.0
       rckonst2=4.0
       rck3=1000000.0*products
C
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
        gaussrand1=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1
C
           if (i-1.eq.0) then
         ckonst1=22*1000000.0*substrates
C
            else
            if ((gaussrand1.gt.0).and.((i-1).gt.0)) then
         ckonst1a(i)=167.4694*gaussrand1
C
          keqk1=ckonst1a(i)
          ckonst1=196.0*keqk1
          k1s=ckonst1/substrates
C
          keqk3=Ktot/(keqk1*keqk2)
C
          ck3=rck3*keqk3
            endif
         endif
C
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
	     enprodB1=0.0
	     enprodB2=0.0
	     enprodB3=0.0
	     enprodT=0.0
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
C
C
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
C
       write(7,'(i5,6(f12.1),(f10.6),3(f11.3),2(f16.6),12(f18.6,1x))')
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

