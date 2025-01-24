      double complex function HPL2at1(n1,n2)
      implicit none
      integer n1,n2,j
      double complex ris,myi
      double precision pi,ll2

      pi=3.1415926535897932385D0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

      if ((j.eq.7).or.(j.eq.9)) then
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL2: "
         print*, "HPL2(",n1,",",n2
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      else
         select case (j)
         case(1)
            ris = ll2**2/2d0
         case(2)
            ris = -pi**2/12d0
         case(3)
            ris = pi**2/12d0 - ll2**2/2d0
         case(4)
            ris = pi**2/12d0
         case(5)
            ris = 0d0
         case(6)
            ris = pi**2/6d0
         case(8)
            ris = -pi**2/6d0
         end select
      endif
      HPL2at1=ris
      return
      end function
