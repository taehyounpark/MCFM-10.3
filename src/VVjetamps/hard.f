!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function hard(p1sq,p2sq,p3sq,p4sq,
     & s12,s23,m1sq,m2sq,m3sq,m4sq)
      implicit none
      include'types.f'
c     Two adjacent external offshellnesses
c     calculated for a limited order of arguments
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
      complex(dp):: hard,Lsm1_2mht
      if (
     &       (p3sq  /=  0d0)
     & .and. (p4sq  /=  0d0)
     & .and. (s12  /=  0d0)
     & .and. (s23  /=  0d0)
     & .and. (p1sq  ==  0d0)
     & .and. (p2sq  ==  0d0)
     & .and. (m1sq  ==  0d0)
     & .and. (m2sq  ==  0d0)
     & .and. (m3sq  ==  0d0)
     & .and. (m4sq  ==  0d0)) then
      hard=2d0*Lsm1_2mht(s12,s23,p3sq,p4sq)/(s23*s12)
      elseif (
     &       (p1sq  /=  0d0)
     & .and. (p4sq  /=  0d0)
     & .and. (s12  /=  0d0)
     & .and. (s23  /=  0d0)
     & .and. (p2sq  ==  0d0)
     & .and. (p3sq  ==  0d0)
     & .and. (m1sq  ==  0d0)
     & .and. (m2sq  ==  0d0)
     & .and. (m3sq  ==  0d0)
     & .and. (m4sq  ==  0d0)) then
      hard=2d0*Lsm1_2mht(s23,s12,p1sq,p4sq)/(s23*s12)
      else
      write(6,*) 'Unimplemented configuration of arguments for hard'
      write(6,*) p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
      stop
      endif
      return
      end
