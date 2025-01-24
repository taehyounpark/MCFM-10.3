!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AGTYassemble_new(ss,tt,uu,F1x1,F2x0,Remain)
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
      real(dp), intent(in) :: ss,tt,uu
      real(dp):: ddilog,Li3,Li4

      real(dp), intent(out) :: F2x0,F1x1,Remain
c=======basis functions from AGTY
      real(dp) :: AGTYAs,AGTYBs,AGTYD2s,AGTYE3s,AGTYG1s

      real(dp) :: Bigx,Bigy,Bigs,x,y,z
      real(dp) Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z,Li4zinv


c----- define various pieces
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq)

      x=-tt/ss
      y=-uu/ss
      z=-uu/tt

      Li2x=ddilog(x)
      Li3x=Li3(x)
      Li4x=Li4(x)

      Li2y=ddilog(y)
      Li3y=Li3(y)
      Li4y=Li4(y)

      Li4z=Li4(z)

      Li4zinv=Li4(1._dp/z)


c======== now build the functions

c====== this bit is special and doesnt scale with the LO quark charge.
      Remain=two*cf*tr*AGTYAs(ss,tt,uu,BigX,BigY,Li2x,Li2y,Li3x,Li3y,
     & Li4x,Li4y,Li4z,Li4zinv)

c======= B piecs (CF**2)
      F2x0=two*Cf**2*AGTYBs(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

c======= D2 piecs (CF*CA)
      F2x0=F2x0
     & +two*CF*CA*AGTYD2s(ss,tt,uu,BigX,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)

c======= E3 pieces (NF*CF)
      F2x0=F2x0+two*nf*Cf*AGTYE3s(tt,uu,BigX,Bigy,Bigs)

      F1x1=CF**2*AGTYG1s(tt,uu,BigX,BigY)


      return
      end
