!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function WWjetcheckpiDpjk(p)
c---- function that returns true if pi.pjk is too small for any
c---- combination of (pi,pj,pk)=(1,2,5,6)
c---- "too small" is set by parameter "tiny" below (default = 1d-6)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      logical:: WWjetcheckpiDpjk
      integer:: n,ifailed
      real(dp):: p(mxpart,4),sij,sijk,sijkl
      real(dp):: piDp34,piDp56,p34Dp56,BoxG,BGa,BGb,BGc,piDpj,
     & pjDp34,pjDp56,BoxG34,BoxG56
c--- 3-mass Box Gram singularities
      integer, parameter :: ib(3)=(/1,2,7/)
c--- 2-mass Box Gram singularities
      integer, parameter :: i2(3)=(/1,1,2/)
      integer, parameter :: j2(3)=(/2,7,7/)
c--- singularities of the form sijk-sij
      integer, parameter :: ii(9)=(/5,3,1,2,5,3,1,5,3/)
      integer, parameter :: jj(9)=(/6,4,7,7,6,4,2,6,4/)
      integer, parameter :: kk(9)=(/1,2,2,1,7,7,7,2,1/)
c--- singularities of the form sijkl-sijk
      integer, parameter :: i4(12)=(/2,1,2,1,1,2,1,2,7,7,7,7/)
      integer, parameter :: j4(12)=(/3,5,3,5,3,5,3,5,5,3,3,5/)
      integer, parameter :: k4(12)=(/4,6,4,6,4,6,4,6,6,4,4,6/)
      integer, parameter :: l4(12)=(/7,7,1,2,7,7,2,1,1,1,2,2/)
c      integer, parameter :: ifailed/0/
c      save ifailed
       real(dp), parameter :: tiny=1d-6

      ifailed=1

      WWjetcheckpiDpjk=.false.

      call dotem(7,p,s)

c--- three-mass box Gram determinants
      p34Dp56=(s(3,5)+s(3,6)+s(4,5)+s(4,6))/2d0
      do n=1,3
        piDp34=(s(ib(n),3)+s(ib(n),4))/2d0
        piDp56=(s(ib(n),5)+s(ib(n),6))/2d0
        BGa=2d0*piDp34*piDp56*p34Dp56
        BGb=s(3,4)*piDp56**2
        BGc=s(5,6)*piDp34**2
        BoxG=abs(BGa-BGb-BGc)
c        write(6,*) 'n,BoxG=',n,BoxG,BoxG/max(abs(BGa),abs(BGb),abs(BGc))
        if (BoxG/max(abs(BGa),abs(BGb),abs(BGc))  <  tiny) then
          ifailed=ifailed+1
c           write(6,*) '>>>>> WARNING:  discarded ',ifailed,
c     &      ' phase space points due to small 3m box Gram <<<<<'
           WWjetcheckpiDpjk=.true.
           return
        endif
      enddo

c--- two-mass box Gram determinants
      do n=1,3
        piDp34=(s(i2(n),3)+s(i2(n),4))/2d0
        piDp56=(s(i2(n),5)+s(i2(n),6))/2d0
        pjDp34=(s(j2(n),3)+s(j2(n),4))/2d0
        pjDp56=(s(j2(n),5)+s(j2(n),6))/2d0
        piDpj=s(i2(n),j2(n))/2d0
        BoxG34=abs(2d0*piDp34*pjDp34-piDpj*s(3,4))
        BoxG56=abs(2d0*piDp56*pjDp56-piDpj*s(5,6))
c        write(6,*) 'n,BoxG=',n,BoxG,BoxG/max(abs(BGa),abs(BGb),abs(BGc))
        if ((BoxG34/abs(piDpj*s(3,4))  <  tiny)
     &  .or.(BoxG56/abs(piDpj*s(5,6))  <  tiny)) then
          ifailed=ifailed+1
c           write(6,*) '>>>>> WARNING:  discarded ',ifailed,
c     &      ' phase space points due to small 2m box Gram <<<<<'
           WWjetcheckpiDpjk=.true.
           return
        endif
      enddo

      do n=1,9
        sij=s(ii(n),jj(n))
        sijk=sij+s(ii(n),kk(n))+s(jj(n),kk(n))
c        write(6,*) ii(n),jj(n),kk(n),abs((sijk-sij)/sij)
        if (abs((sijk-sij)/sij)  <  tiny) then
          ifailed=ifailed+1
c          write(6,*) 'small pk.pij: k,i,j=',kk(n),ii(n),jj(n)
c           write(6,*) '>>>>> WARNING:  discarded ',ifailed,
c     &      ' phase space points due to small pk.pij <<<<<'
           WWjetcheckpiDpjk=.true.
           return
        endif
      enddo

      do n=1,12
        sijk=s(i4(n),j4(n))+s(i4(n),k4(n))+s(j4(n),k4(n))
        sijkl=sijk+s(i4(n),l4(n))+s(j4(n),l4(n))+s(k4(n),l4(n))
c        write(6,*) i4(n),j4(n),k4(n),l4(n),abs((sijkl-sijk)/sijk)
        if (abs((sijkl-sijk)/sijk)  <  tiny) then
          ifailed=ifailed+1
c          write(6,*) 'small pk.pij: k,i,j=',kk(n),ii(n),jj(n)
c           write(6,*) '>>>>> WARNING:  discarded ',ifailed,
c     &      ' phase space points due to small pk.pij <<<<<'
           WWjetcheckpiDpjk=.true.
           return
        endif
      enddo

      return
      end

