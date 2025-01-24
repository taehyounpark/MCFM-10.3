      double complex function HPL3at1(n1,n2,n3)
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

      if(j.le.18.or.j.eq.23) then
         select case (j)
         case (1)
            ris =ll2**3d0/6d0
         case (2)
            ris =-(pi**2d0*ll2)/12d0 + zeta3/8d0
         case (3)
            ris =-ll2**3d0/6d0 + zeta3/8d0
         case (4)
            ris =(pi**2d0*ll2)/12d0 - zeta3/4d0
         case (5)
            ris =(3d0*zeta3)/4d0
         case (6)
            ris =(pi**2d0*ll2)/6d0 - (5d0*zeta3)/8d0
         case (7)
            ris =(pi**2d0*ll2)/12d0 
     &           - ll2**3d0/6d0-zeta3/4d0
         case (8)
            ris =(pi**2d0*ll2)/12d0 - zeta3
         case (9)
            ris =-(pi**2*ll2)/12d0+ll2**3/6d0
     &           +(7d0*zeta3)/8d0
         case (10)
            ris =zeta3/8d0
         case (11)
            ris =(-3d0*zeta3)/2d0
         case (12)
            ris =-(pi**2d0*ll2)/4d0 + (13d0*zeta3)/8d0
         case (13)
            ris =(3d0*zeta3)/4d0
         case (14)
            ris =0d0
         case (15)
            ris =zeta3
         case (16)
            ris =(pi**2d0*ll2)/4d0 - zeta3
         case (17)
            ris =-2d0*zeta3
         case (18)
            ris =zeta3
         case (23)
            ris =zeta3
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL3: "
         print*, "HPL3(",n1,",",n2,",",n3
     &        ,",1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL3at1=ris
      return
      end function
