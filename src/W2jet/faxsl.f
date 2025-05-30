!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Faxsl(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Faxsl
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'heldefs.f'
      integer st
      complex(dp):: L1
      real(dp):: t,mtsq
      mtsq=mt**2
      if(st==hqpqbmgpgm) then
      Faxsl=
     &(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j4),-s(j5,j6))/(2._dp*s(j5,j6)))*zb(
     &j1,j3)*
     &(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))*zb(j3,j6))/
     &(s(j5,j6)*zb(j1,j4)*zb(j2,j4))+
     &(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j3),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     &j2,j4)*za(j4,j5)*
     &(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6)))/
     &(s(j5,j6)*za(j1,j3)*za(j2,j3))
      elseif(st==hqpqbmgpgp) then
      Faxsl=
     &(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j4),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     &j2,j5)*
     &(-(za(j1,j2)*zb(j1,j3))-za(j2,j4)*zb(j3,j4))*zb(j3,j6))/
     &(s(j5,j6)*za(j1,j4)*za(j2,j4))+
     &(2._dp*(1/(24._dp*mtsq)-L1(-t(j1,j2,j3),-s(j5,j6))/(2._dp*s(j5,j6)))*za(
     &j2,j5)*
     &(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*zb(j4,j6))/
     &(s(j5,j6)*za(j1,j3)*za(j2,j3))
      endif
      return
      end
