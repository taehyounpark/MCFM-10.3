!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine CSsq(order,Mhsq,musq,coeff)
      implicit none
c   Wilson Coefficients for gg->H
c   Taken from Ahrens,Becher,Neubert,Yang, 0809.4283v3
c   Also allows for timelike value of musq and corresponding complex alphas
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nfl.f'
      include 'scet_const.f'
      integer,intent(in)::order
      real(dp),intent(in)::Mhsq,musq
      real(dp),intent(out)::coeff(2)
      complex(dp)::L,lnrat
      complex(dp)::Cs(2)
      complex(dp)::ason2pi
      complex(dp)::getason2pic

      L=lnrat(-Mhsq,musq)
      
      Cs(:)=zip
c----Eq.(17) of 0809.4283v3
      if (order > 0) then
      Cs(1)=CA*(-L**2+zeta2)
      endif
      if (order > 1) then
      Cs(2)=
     & +CA**2*(1/2._dp*L**4+11/9._dp*L**3+(-67/9._dp+zeta2)*L**2
     & +L*(80/27._dp-22/3._dp*zeta2-2*zeta3)
     & +5105/162._dp+67/6._dp*zeta2+1/2._dp*zeta2**2-143/9._dp*zeta3)
     & +TR*nfl*CF*(4*L-67/3._dp+16*zeta3)
     & +TR*nfl*CA*(
     & -4/9._dp*L**3+20/9._dp*L**2+(104/27._dp+8/3._dp*zeta2)*L
     & -1832/81._dp-10/3._dp*zeta2-92/9._dp*zeta3)
      endif

! allow for complex alpha-s (timelike musq)
      ason2pi=getason2pic(musq)

! apply factors of ason2pi
      Cs(1)=Cs(1)*ason2pi
      Cs(2)=Cs(2)*ason2pi**2
       
      coeff(1)=real(Cs(1),dp)
      coeff(2)=(2*real(Cs(2),dp)+real(Cs(1)*conjg(Cs(1)),dp))/four

      return
      end

