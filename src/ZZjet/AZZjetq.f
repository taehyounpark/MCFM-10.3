!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjetq(p,j1,j2,j3,j4,j5,j6,j7,A7qv,A7qvmLL,A7qvmLR)
c--- compute amplitude A^a_7 and A^q_7 virtual, including both
c--- leading and subleading pieces, for both positive
c--- and negative gluon helicities, and for all quark and lepton helicities
c--- labelling is A7xxx(gluehel, Zlephel)
c--- Note: includes both the usual strong coupling renormalization
c--- and also the finite renormalization in the dred scheme
c--- leading order amplitudes are also returned
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'WWjetlabels.f'
      include 'masses.f'
      integer j1,j2,j3,j4,j5,j6,j7,k1,k2,k3,k4,k5,k6,i12,i34,i56
      real(dp):: p(mxpart,4),mq
      complex(dp):: A7qv(2,2,2,2),A7qvmLL(2,2,2,2),A7qvmLR(2,2,2,2),
     & A7q(0:2),scintsmq(4,60,-2:0)
      character, parameter:: ah(2) = (/ 'm', 'p' /)
      character(len=4):: help,helm

      call spinoru(7,p,za,zb)

      call WWjetcomputescalars(j1,j2,j3,j4,j5,j6,j7,scints)

      mq=mt

c uncomment these for printing out KCheck comparison
c (also need to use kinpoint_new.f in setupmcfmmom)
c      mq = sqrt(7d0/11d0)
c      verbose=.true.
c      call spinorz(7,p,za,zb)

      call ZZjetqloopcomputescalars(mq,j1,j2,j3,j4,j5,j6,j7,scintsmq)

      do i12=1,2
      do i34=1,2
      do i56=1,2

      if (i12 == 1) then
        k1=j1
        k2=j2
      else
        k1=j2
        k2=j1
      endif
      if (i34 == 1) then
        k3=j3
        k4=j4
      else
        k3=j4
        k4=j3
      endif
      if (i56 == 1) then
        k5=j5
        k6=j6
      else
        k5=j6
        k6=j5
      endif

      help=ah(i12)//ah(i34)//ah(i56)//'p'
      helm=ah(i12)//ah(i34)//ah(i56)//'m'

c--- positive helicity gluon
c--- down-type amplitudes
      call AZZjet_qloop_new(help,.false.,k1,k2,k3,k4,k5,k6,j7,za,zb,mq,scints,scintsmq,A7q)
      A7qv(i12,i34,i56,2)=+A7q(0)
      A7qvmLL(i12,i34,i56,2)=+A7q(1)
      A7qvmLR(i12,i34,i56,2)=+A7q(2)

c--- negative helicity gluon
c--- down-type amplitudes
c (These calls are equivalent but clearly better here to never use swapVV)
c      call AZZjet_qloop_new(helm,.true.,k2,k1,k6,k5,k4,k3,j7,zb,za,scints,A7q)
      call AZZjet_qloop_new(helm,.false.,k2,k1,k4,k3,k6,k5,j7,zb,za,mq,scints,scintsmq,A7q)
      A7qv(i12,i34,i56,1)=-A7q(0)
      A7qvmLL(i12,i34,i56,1)=-A7q(1)
      A7qvmLR(i12,i34,i56,1)=-A7q(2)

      enddo
      enddo
      enddo

      return
      end

