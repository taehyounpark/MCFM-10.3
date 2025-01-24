!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function easy(p1sq,p2sq,p3sq,p4sq,
     & s12,s23,m1sq,m2sq,m3sq,m4sq)
c      Two opposite external offshellnesses
      implicit none
      include'types.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
      complex(dp):: easy,Lsm1_2me
      if (
     &       (p2sq  /=  0d0)
     & .and. (p4sq  /=  0d0)
     & .and. (s12  /=  0d0)
     & .and. (s23  /=  0d0)
     & .and. (p1sq  ==  0d0)
     & .and. (p3sq  ==  0d0)
     & .and. (m1sq  ==  0d0)
     & .and. (m2sq  ==  0d0)
     & .and. (m3sq  ==  0d0)
     & .and. (m4sq  ==  0d0)) then
      easy=2d0*Lsm1_2me(s12,s23,p2sq,p4sq)/(s23*s12-p4sq*p2sq)
      else
      write(6,*) 'Unimplemented configuration of arguments for easy'
      write(6,*) p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
      stop
      endif
      return
      end
