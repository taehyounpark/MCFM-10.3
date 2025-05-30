!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gagaAGTYassemble(ss,tt,uu,Apie,Qpie)
c======C.Williams August 2015
c======routine to reproduce eq. 4.6 of AGTY, without overall
c=====factor of N_c, which lives with the LO in our units (factor of 2 is included)

c=====Apie is the piece which goes like sum over quark charges
c=====Qpie is the piece which goes like the LO quark charge

c=====note that neither Apie nor Qpie is dressed with EM charges in this
c===== routine (for consistency)

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scale.f'
      include 'scet_const.f'
      real(dp), intent(in) :: ss,tt,uu
      real(dp), intent(out) :: Apie,Qpie
      real(dp):: ddilog,Li3,Li4

c=======basis functions from AGTY
      real(dp) :: AGTYAs,AGTYBs,AGTYD2s
      real(dp) :: AGTYE3s

      real(dp) :: Bigx,Bigy,Bigs,BigT,BigU,x,y,z
      real(dp) Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li2z,Li3z,Li4z,Li4zinv

      include 'scet_beta.f'

c----- define various pieces
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq)
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)

      x=-tt/ss
      y=-uu/ss
      z=-uu/tt

      Li2x=ddilog(x)
      Li3x=Li3(x)
      Li4x=Li4(x)

      Li2y=ddilog(y)
      Li3y=Li3(y)
      Li4y=Li4(y)

      Li2z=ddilog(z)
      Li3z=Li3(z)
      Li4z=Li4(z)

      Li4zinv=Li4(1._dp/z)


c======== now build the functions

      Apie=two*cf*tr*AGTYAs(ss,tt,uu,BigX,BigY,Li2x,Li2y,Li3x,Li3y,
     & Li4x,Li4y,Li4z,Li4zinv)

c======= B pieces (CF**2)
      Qpie=two*Cf**2*AGTYBs(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

c======= D2 pieces (CF*CA)
      Qpie=Qpie+two*CF*CA*AGTYD2s(ss,tt,uu,BigX,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

c======= E3 pieces (NF*CF)
      Qpie=Qpie+two*nf*Cf*AGTYE3s(tt,uu,BigX,Bigy,Bigs)

      return
      end
