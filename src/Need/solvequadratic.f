!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine solvequadratic(a,b,c,lamp,lamm)
      implicit none 
      include 'types.f'
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
      
