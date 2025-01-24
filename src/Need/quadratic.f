!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module quadratic
      implicit none
      public solvequadratic,solvequadratic_r,solvequadratic_c

      interface solvequadratic
      module procedure solvequadratic,solvequadratic_r,solvequadratic_c
      end interface

      contains 
      
      subroutine solvequadratic(a,b,c,lamp,lamm)
      use types
      real(dp):: a,b,c,q,rt,lamp,lamm
      rt=b**2-4*a*c
      if (rt .ge. 0._dp) then
         rt=sqrt(rt)
         if (b .ge. 0) then
            q=-0.5d0*(b+rt)
         else
            q=-0.5d0*(b-rt)
         endif
      else
         write(6,*) 'no real roots'
         stop
      endif
      lamp=q/a
      lamm=c/q
      end subroutine solvequadratic
      
      subroutine solvequadratic_r(a,b,c,lamp,lamm)
      use types
      real(dp):: a,b,c
      complex(dp):: q,rt,lamp,lamm
      rt=cmplx(b**2-4*a*c,kind=dp)
      rt=sqrt(rt)
      if (b .ge. 0) then
         q=-0.5d0*(b+rt)
      else
         q=-0.5d0*(b-rt)
      endif
      lamp=q/a
      lamm=c/q
      end subroutine solvequadratic_r

      subroutine solvequadratic_c(a,b,c,lamp,lamm)
      use types
      complex(dp):: a,b,c,q,rt,lamp,lamm
      rt=b**2-4*a*c
      rt=sqrt(rt)
      if (real(conjg(b)*rt) .gt. 0._dp) then
      q=-0.5_dp*(b+rt)
      else
      q=-0.5_dp*(b-rt)
      endif
      lamp=q/a
      lamm=c/q
      end subroutine solvequadratic_c

      end module quadratic
