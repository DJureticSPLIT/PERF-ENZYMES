       program     HsKYNase
C
C Four-state human enzyme kynureninase variant 93D9 kinetic scheme (HsKYNase-93D9)
C
C Rate constants are from Table 3 of Karamitros-2024 paper.
C The substrate is L-kynurenine at the concentration of 0.000075 M. First product is the anthranilic acid.
C The second product is L-alanine. The choice of product concentrations is 0.00075 M.
C The substrate is L-kynurenine at low concentration of because higher concentrations are inhibitory.
C
C The output of the program HsKYN93D9-SI.for are three data files:
C a) KYN93D9SI-exper.dat containing the row #46 parameters from the SI excel dataset (s=sum of probabilities).
C b) KYN93D9SimI-K1K2K3K4c.dat containing the 2000 rows with the Simulation I results. The first row
C repeats KYN93D9SI parameters. Rows with negative random numbers (rand) are replaced
C with the KYN93D9SI-exper.dat row.
C c) KYN93D9SimII-K1K4k1k7.dat containing the 2000 rows with the Simulation II results. The first row
C repeats KYN93D9SI parameters. Rows with negative random numbers (300 rows) are replaced with the first
C previous row with a positive random number.
C The output parameters are ordered as follows: microscopic rate constants, force X/RT, force X,
C total equilibrium constant Ktot, net flux J, total dissipation Ptot, enzyme efficiency Eff,
C catalytic constant k(cat), partial dissipations P1, P2, P3, and P4, substrate S and product P concentration,
C equilibrium constant K1, Michaelis-Menten constant Km, information entropy S(I),
C state probabilities p1, p2, p3, and p4, and random number rand.
C The first row's k(cat) value should be similar to the 0.75 inverse seconds (Karamitros-2024).
C
C The k(cat) value from the first row should be similar to the 0.68 inverse seconds (Karamitros-2024).
C
       real*8  ckonst1a(2000),ckonst2a(2000),
     &   thermale,rckonst1,rckonst2,
     &   xsec,sigma,enprodB1,enprodB2,enprodB3,enprodB4,
     &   effic,k12,k21,k23,k32,k34,k43,k41,k14,
     &   keqk12,keqk23,keqk34,keqk41,
     &   sigmap,sigmab,sigmah,sigmaq,
     &   currB1,currB2,currB3,currB4,
     &   ckonst3a(2000),ckonst4a(2000),
     &   enprodT,s1,s2,gaussrand1(2000),pi,gaussrand(2000),
     &   sprob(4),entros(4),entrops,ENZeff,k12s,
     &   aX1,aX2,aX3,aX4,kcat,aXtot,Ktot,KMich,substrates,products,
     &   Jplus,Jminus,Jratio,Kplus,Kminus,Kratio
C
       s1=0.0
       s2=0.0
       pi=0.0
       keqk12=0.0
       keqk23=0.0
       keqk34=0.0
       keqk41=0.0
       aX1=0.0
       aX2=0.0
       aX3=0.0
       aX4=0.0
       enprodB3=0.0
       enprodB4=0.0
       enprodB1=0.0
       enprodB2=0.0
       enprodT=0.0
C
       sigmap=0.0
       sigmab=0.0
       sigmah=0.0
       sigmaq=0.0
C
       k12=0.0
       k21=0.0
       k23=0.0
       k32=0.0
       k34=0.0
       k43=0.0
       k41=0.0
       k14=0.0
       k12s=0.0
       currB1=0.0
       currB2=0.0
       currB3=0.0
       currB4=0.0
C
       thermale=0.0
C
       rckonst1=0.0
       rckonst2=0.0
C
       sigma=0.0
C
       kcat=0.0
       aXtot=0.0
       Ktot=0.0
       KMich=0.0
       ENZeff=0.0
       substrates=0.0
       products=0.0
       Jplus=0.0
       Jminus=0.0
       Jratio=0.0
       Kplus=0.0
       Kminus=0.0
       Kratio=0.0
       mmm=0
       entrops=0.0
       nnn=0
       jjj=0
C
C
       do 85 nnn=1,2000
       ckonst1a(nnn)=0.0
       ckonst2a(nnn)=0.0
       ckonst3a(nnn)=0.0
       ckonst4a(nnn)=0.0
       gaussrand(nnn)=0.0
       gaussrand1(nnn)=0.0
85     continue
C
C
            open (unit=11,status='unknown',file='KYN93D9SI-exper.dat',
     & access='sequential',form='formatted')
          write (11,115)
115     format(2x,'     k1              k2           k3            k4 ',
     & '          k5            k6             k7             k8      ',
     & '       X/RT          X           Ktot          J ',
     & '            Ptot        Eff              kcat           P1 ',
     & '          P2              P3          P4            S  ',
     & '         P             s              Km            S(I)',
     & '         p1          p2          p3           p4')
C
C
            products=0.00075
            substrates=0.000075
C
          k12=240000.0*substrates
          k21=4.0
          k23=0.79
          k32=0.34
          k34=11.5
          k43=1000.0*products
          k41=11.0
          k14=20000.0*products
C
	     keqk12=k12/k21
	     keqk23=k23/k32
	     keqk34=k34/k43
	     keqk41=k41/k14
C
       Jplus=k12*k23*k34*k41
       Jminus=k21*k32*k43*k14
       Jratio=Jplus/Jminus
       Kplus=Jplus/substrates
CC                 division with the substrate conc. in moles
       Kminus=Jminus/products
CC                 division with the product conc. in moles
       Kratio=Kplus/Kminus
C
          thermale=2.57873058
C RT in kJ/mol for 37 C
C
         pi=3.141592653589793
C
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
C 4-th state probability added p=1,b=2,h=3,q=4
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
	     enprodB3=0.0
         enprodB4=0.0
	     enprodT=0.0
CC
           do 75 jj=1,4
         sprob(jj)=0.0
         entros(jj)=0.0
75       continue
C
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
         sigma=sigmab+sigmah+sigmaq+sigmap
C
         sprob(1)=sigmab
         sprob(2)=sigmah
         sprob(3)=sigmaq
         sprob(4)=sigmap
C
        do 750 mmmm=1,4
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
750     continue
C
C
        xsec=thermale*2.3025851*dlog10((k12*k23*k34*k41)/
     & (k21*k32*k43*k14))
Cn all currents are without N (N=1)
Cn
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
          aXtot=aX1+aX2+aX3+aX4
C
Cn All "entropy productions" are XJ/RT expressions
Cn
        enprodB1=aX1*currB1
        enprodB2=aX2*currB2
        enprodB3=aX3*currB3
        enprodB4=aX4*currB4
C
        enprodT=enprodB1+enprodB2+enprodB3+enprodB4
C
C These enprod expressions are entropy production divided by the gas const R
        write(11,'(19(1x,f13.7),2x,3(f12.7,2x),1x,6(f12.7,1x))')
     & k12,k21,k23,k32,k34,k43,k41,k14,aXtot,xsec,Ktot,currB1,
     & enprodT,Enzeff,kcat,enprodB1,enprodB2,enprodB3,enprodB4,
     & substrates,products,sigma,KMich,entrops,
     & sigmab,sigmah,sigmaq,sigmap
C
C
        open (unit=13,status='unknown',file='KYN93D9SimII-K1K4k1k7.dat',
     & access='sequential',form='formatted')
       write(13,65)
65     format(2x,'  i          k1              k2               k3    ',
     &  '           k4           k5               k6               k7 ',
     &  '             k8             X/RT             X        ',
     &  '    Ktot              J               Ptot          Eff      ',
     &  '         kcat            P1               P2              P3',
     &  '                 P4            S                P           ',
     &  '        K1            K(M)            S(I)           p1     ',
     &  '         p2               p3             p4             rand')
C
C
C
            products=0.00075
            substrates=0.000075
C
          k12=240000.0*substrates
          k21=4.0
          k23=0.79
          k32=0.34
          k34=11.5
          k43=1000.0*products
          k41=11.0
          k14=20000.0*products
C
	     keqk12=k12/k21
	     keqk23=k23/k32
	     keqk34=k34/k43
	     keqk41=k41/k14
C
          thermale=2.57873058
C RT in kJ/mol for 37 C
         pi=3.141592653589793
C
         Ktot=keqk12*keqk23*keqk34*keqk41
C
        do 50 i=1,2000
C
        call random_number(s1)
        call random_number(s2)
C
        gaussrand1(i)=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1.0
C
C
C
         if (i-1.eq.0) then
         k12=240000.0*substrates
         k12s=240000.0
            else
            if ((gaussrand1(i).gt.0).and.((i-1).gt.0)) then
         ckonst1a(i)=4.5*gaussrand1(i)
C
          keqk12=ckonst1a(i)
          k12=4.0*keqk12
          k12s=k12/substrates
C
          keqk41=Ktot/(keqk12*keqk23*keqk34)
C
          k41=k14*keqk41
            endif
         endif
C
       kcat=(k23*k34*k41)/((k23+k32)*(k43+k41)+k23*k34+k34*k41)
C
C
       KMich=(k21*(k32*(k43+k41)+k34*k41)+k23*k34*k41)*substrates/
     &       (k12*((k23+k32)*(k43+k41)+k23*k34+k34*k41))
C
        ENZeff=kcat/KMich
C
C
         sigma=0.0
         sigmab=0.0
         sigmah=0.0
         sigmaq=0.0
         sigmap=0.0
C p=1,b=2,h=3,q=4
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
         sigma=sigmab+sigmah+sigmaq+sigmap
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
Cn
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
          aXtot=aX1+aX2+aX3+aX4
C
        enprodB1=aX1*currB1
        enprodB2=aX2*currB2
        enprodB3=aX3*currB3
        enprodB4=aX4*currB4
C
        enprodT=enprodB1+enprodB2+enprodB3+enprodB4
C
C These enprod expressions are entropy production divided by the gas const R
C
       write(13,'(i5,2x,19(1x,f15.7),2x,4(f15.9,2x),1x,6(f14.10,1x))')
     & i,k12,k21,k23,k32,k34,k43,k41,k14,aXtot,xsec,Ktot,currB1,
     & enprodT,Enzeff,kcat,enprodB1,enprodB2,enprodB3,enprodB4,
     & substrates,products,keqk12,KMich,
     & entrops,sigmab,sigmah,sigmaq,sigmap,gaussrand1(i)
C
50    continue
C
       close (unit=13)
C
C
        open (unit=9,status='unknown',file='KYN93D9SimI-K1K2K3K4c.dat',
     & access='sequential',form='formatted')
       write(9,95)
95     format(2x,' step       k1            k2           k3      ',
     & '      k4          k5           k6           k7            k8 ',
     & '         Ptot             Eff           kcat      ',
     & ' Km            P1         P2          P3           P4     ',
     & '      J             S(I)        p1            p2        ',
     & '     p3          p4          Ktot          X            rand')
C
C
            products=0.00075
            substrates=0.000075
C
C in moles
C      Kinetic constants without [S] and [P] are from Toney-2019
          k12=240000.0*substrates
          k21=4.0
          k23=0.79
          k32=0.34
          k34=11.5
          k43=1000.0*products
          k41=11.0
          k14=20000.0*products
C
C
       keqk23=k23/k32
       keqk34=k34/k43
       keqk41=k41/k14
   	   keqk12=k12/k21
C
       Ktot=keqk12*keqk23*keqk34*keqk41
C
C
         do 60 i=1,2000
       do 185 nnn=1,2000
       ckonst1a(nnn)=0.0
       ckonst2a(nnn)=0.0
       ckonst3a(nnn)=0.0
       ckonst4a(nnn)=0.0
       gaussrand(nnn)=0.0
185     continue
        call random_number(s1)
        call random_number(s2)
C
        gaussrand(i)=sqrt(-2*2.3025851*dlog10(s1))*cos(2*pi*s2)+1.0
C
C
         if (i-1.eq.0) then
          k12=240000.0*substrates
          k21=4.0
          k23=0.79
          k32=0.34
          k34=11.5
          k43=1000.0*products
          k41=11.0
          k14=20000.0*products
C
         keqk23=k23/k32
         keqk34=k34/k43
         keqk41=k41/k14
   	     keqk12=k12/k21
C
         Ktot=keqk12*keqk23*keqk34*keqk41
C
            else
            if ((gaussrand(i).gt.0).and.((i-1).gt.0)) then
            products=0.00075
            substrates=0.000075
          k12=240000.0*substrates
          k21=4.0
          k23=0.79
          k32=0.34
          k34=11.5
          k43=1000.0*products
          k41=11.0
          k14=20000.0*products
C
          ckonst1a(i)=k12*gaussrand(i)
          k12=ckonst1a(i)
          k21=k12/keqk12
C
          ckonst2a(i)=k23*gaussrand(i)
          k23=ckonst2a(i)
          k32=k23/keqk23
C
          ckonst3a(i)=k34*gaussrand(i)
          k34=ckonst3a(i)
          k43=k34/keqk34
C
          ckonst4a(i)=k41*gaussrand(i)
          k41=ckonst4a(i)
          k14=k41/keqk41
C
            endif
         endif
C
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
C 4-th state probability added p=1,b=2,h=3,q=4
         entrops=0.0
	     enprodB1=0.0
	     enprodB2=0.0
	     enprodB3=0.0
         enprodB4=0.0
	     enprodT=0.0
C
           do 105 jj=1,4
         sprob(jj)=0.0
         entros(jj)=0.0
105       continue
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
       do 800 mmmm=1,4
          entros(mmmm)=-sprob(mmmm)*2.3*dlog10(sprob(mmmm))
          entrops=entrops+entros(mmmm)
800    continue
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
       write(9,'(i6,8(f12.4,1x),3(f14.4,1x),1x,f9.6,1x,13(f12.5,1x))')
     &  i,k12,k21,k23,k32,k34,k43,k41,k14,enprodT,ENZeff,kcat,KMich,
     &  enprodB1,enprodB2,enprodB3,enprodB4,
     &  currB1,entrops,sigmap,sigmab,sigmah,sigmaq,Ktot,aX,gaussrand(i)
C
 60    continue
C
       close(unit=9)
C
250     continue
C
C
       stop
       end
C Author: Davor Juretic, year:2024
C
C The run time should be around 2.0 sec.

