      double complex function HPL2atm1(n1,n2)
      implicit none
      integer n1,n2,j
      double complex ris,myi
      double precision pi,ll2

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if(j.gt.3) then
         select case (j)
         case(4)
            ris = -pi**2/6d0
         case(5)
            ris = -pi**2/2d0
         case(6)
            ris = -pi**2/12d0
         case(7)
            ris = pi**2/12d0 - ll2**2/2d0
         case(8)
            ris = pi**2/12d0 - myi*pi*ll2
         case(9)
            ris = ll2**2/2d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL2atm1=ris
      return
      end function
