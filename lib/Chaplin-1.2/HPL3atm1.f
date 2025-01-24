      double complex function HPL3atm1(n1,n2,n3)
      implicit none
      integer n1,n2,n3,j
      double complex ris,myi
      double precision pi, zeta2, zeta3,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

      if(j.ge.10) then
         select case (j)
         case (10)
            ris = zeta3
         case (11)
            ris = 2*zeta3 - myi*pi*zeta2
         case(12)
            ris = pi**2*ll2/4d0 - zeta3
         case(13)
            ris = -zeta3
         case(14)
            ris = -myi*pi*zeta2
         case(15)
            ris = -3*zeta3/4d0
         case(16)
            ris = -pi**2*ll2/4d0 + 13d0*zeta3/8d0
         case(17)
            ris = 3*zeta3/2d0 - myi*pi**3/12d0
         case(18)
            ris = zeta3/8d0
         case(19)
            ris = 0.5d0*zeta2*ll2 - ll2**3/6d0 
     &           - 7d0*zeta3/8d0
         case(20)
            ris = 0.5d0*zeta2*ll2 + myi*pi*(0.5d0*zeta2 
     &           - 0.5d0*ll2**2) - zeta3
         case(21)
            ris = -0.5d0*zeta2*ll2 + ll2**3/6d0 
     &           + zeta3/4d0
         case(22)
            ris = zeta2*ll2 - 5d0*zeta3/8d0
         case(23)
            ris = 0.5d0*myi*pi*zeta2 + pi**2*ll2/2d0 
     &           - 3*zeta3/4d0
         case(24)
            ris = 0.5d0*zeta2*ll2 - zeta3/4d0
         case(25)
            ris = ll2**3/6d0 - zeta3/8d0
         case(26)
            ris = -0.5d0*zeta2*ll2 + zeta3/8d0 
     &           + 0.5d0*myi*pi*ll2**2
         case(27)
            ris = -1d0/6d0*ll2**3
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3atm1=ris
      return
      end function
