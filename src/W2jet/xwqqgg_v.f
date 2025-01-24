!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine xwqqgg_v(mqqb)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author J.M.Campbell, February 2000                               *
c     Returns the interference of the tree and loop                    *
c     amplitudes for the process                                       *
c     0---> q(p1)+g(p2)+g(p3)+qbar(p4)+l(p5)+a(p6)                     *
c                                                                      *
c     The value of COLOURCHOICE determines which colour structures     *
c     are included in the calculation                                  *
c***********************************************************************

      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      include 'heldefs.f'
      integer:: j,lh,h2,h3,hq,h(2:3)
      real(dp):: mqqb,fac
      complex(dp):: m(2),ml1(2),ml2(2),ml3,ml4(2)
      complex(dp):: a6treeg1,a61g1lc,a61g1slc,a61g1nf,a63g1
      integer,parameter::
     & i1(2)=(/1,4/),i2(2)=(/2,3/),i3(2)=(/3,2/),
     & i4(2)=(/4,1/),i5(2)=(/6,5/),i6(2)=(/5,6/)
      integer,parameter:: st1(2,2)=reshape(
     & (/hqpgmgmqbm ,hqpgmgpqbm ,hqpgpgmqbm ,hqpgpgpqbm /),(/2,2/))
      integer,parameter:: st2(2,2)=reshape(
     & (/hqpqbmgpgp ,hqpqbmgpgm ,hqpqbmgmgp ,hqpqbmgmgm /),(/2,2/))
      integer,parameter:: st3(2,2)=reshape(
     & (/hqpqbmgmgm ,hqpqbmgmgp ,hqpqbmgpgm ,hqpqbmgpgp /),(/2,2/))


c ---final matrix element squared is needed for left-handed
c--  quark and lepton line helicity

      fac=avegg*8._dp*gsq**2*esq**2*cf*xn**3*ason2pi
c--- no extra factor here since colour algebra is already done in (2.12)

      hq=1
      lh=1
      mqqb=0._dp

      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
        do j=1,2
        if (hq == 1) then
        m(j)=  a6treeg1(st1(3-h(i2(j)),3-h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        if (colourchoice <= 2) then
        ml1(j)= a61g1lc(st1(3-h(i2(j)),3-h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        endif
        if ((colourchoice >= 2) .or. (colourchoice == 0)) then
        ml2(j)=a61g1slc(st2(3-h(i2(j)),3-h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        endif
        if (colourchoice == 0) then
        ml4(j)=a61g1nf(st1(3-h(i2(j)),3-h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
        endif
        else
        m(j)=  a6treeg1(st1(h(i2(j)),h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        if (colourchoice <= 2) then
        ml1(j)= a61g1lc(st1(h(i2(j)),h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        endif
        if ((colourchoice >= 2) .or. (colourchoice == 0)) then
        ml2(j)=a61g1slc(st2(h(i2(j)),h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        endif
        if (colourchoice == 0) then
        ml4(j)=a61g1nf(st1(h(i2(j)),h(i3(j))),
     &     i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
        endif
        endif
        enddo

        if ((colourchoice == 2) .or. (colourchoice == 0)) then
        if (hq == 1) then
        ml3=a63g1(st3(3-h2,3-h3),1,4,2,3,i6(lh),i5(lh),zb,za)
        else
        ml3=a63g1(st3(h2,h3),1,4,2,3,i5(lh),i6(lh),za,zb)
        endif
        endif

      if     (colourchoice == 1) then
        mqqb=mqqb+fac*(real(conjg(m(1))*ml1(1),dp)
     &                +real(conjg(m(2))*ml1(2),dp))
      elseif (colourchoice == 2) then
        mqqb=mqqb+fac*(
     &    real(conjg(m(1))*(
     &     -(ml1(1)+ml2(1)+ml1(2)-ml3)/xnsq
     &     -(ml2(1)+ml2(2))/xnsq),dp)
     &   +real(conjg(m(2))*(
     &     -(ml1(2)+ml2(2)+ml1(1)-ml3)/xnsq
     &     -(ml2(1)+ml2(2))/xnsq),dp))
      elseif (colourchoice == 3) then
        mqqb=mqqb+fac*(1._dp+xnsq)/xnsq**2*
     &                (real(conjg(m(1))*(ml2(1)+ml2(2)),dp)
     &                +real(conjg(m(2))*(ml2(1)+ml2(2)),dp))
      else
        mqqb=mqqb+fac*(
     &    real(conjg(m(1))*(
     &      ml1(1)
     &     -(ml1(1)+ml2(1)+ml1(2)-ml3)/xnsq
     &     +(ml2(1)+ml2(2))/xnsq**2
     &     +ml4(1)/xn
     &     -(ml4(1)+ml4(2))/xn**3),dp)
     &   +real(conjg(m(2))*(
     &      ml1(2)
     &     -(ml1(2)+ml2(2)+ml1(1)-ml3)/xnsq
     &     +(ml2(1)+ml2(2))/xnsq**2
     &     +ml4(2)/xn
     &     -(ml4(1)+ml4(2))/xn**3),dp))
      endif

      enddo
      enddo


      return
      end

