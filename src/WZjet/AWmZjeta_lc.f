!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWmZjeta_lc(swap12,swapVV,p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz,scints,coeff,Alo,A7a_lc)
c--- MCFM notation
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'WWjetlabels.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      complex(dp):: Alo,A7treea,A7a_lc,lnrat
      complex(dp):: tri,tri_sl
      complex(dp):: xd7x1x56,xd1x27x34,xd7x2x34,xd2x17x56,xd1x7x2,
     & bub17,bub156,bub234,bub127,bub56
      real(dp):: p(mxpart,4)
      logical:: swap12,swapVV,swapz

c--- call standard MCFM routine to get LO
      Alo=A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)


c Note interchange of 34 and 65 in calls below, to reflect W+W- vs W-W+

      coeff(4,d7x1x34)=xd7x1x56(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d27x1x34)=xd1x27x34(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d7x2x56)=xd7x2x34(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d17x2x56)=xd2x17x56(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x7x2)=xd1x7x2(j1,j2,j6,j5,j3,j4,j7,za,zb)

c--- Triangles (set up subleading color here too)
      call tri56x34(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c56x34)=tri
      coeff(3,c56x34sl)=tri_sl

      call tri56x17(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c17x34)=tri
      coeff(3,c17x34sl)=tri_sl

      call tri34x27(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c27x56)=tri
      coeff(3,c27x56sl)=tri_sl

      call tri12x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c12x56sl)=tri_sl

      call tri12x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c12x34sl)=tri_sl

c--- Leading color bubbles
      coeff(2,b17)=bub17(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b134)=bub156(j1,j2,j6,j5,j4,j3,j7,za,zb)
      coeff(2,b256)=bub234(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff(2,b127)=bub127(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b34)=bub56(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
c--- remaining bubble from pole structure
      coeff(2,b56)=-1.5_dp*Alo-coeff(2,b17)-coeff(2,b134)
     & -coeff(2,b256)-coeff(2,b127)-coeff(2,b34)

c Now assemble total
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        A7a_lc=
     &  +scints(4,d7x1x34, 0)*coeff(4,d7x1x34)
     &  +scints(4,d27x1x34,0)*coeff(4,d27x1x34)
     &  +scints(4,d7x2x56, 0)*coeff(4,d7x2x56)
     &  +scints(4,d17x2x56,0)*coeff(4,d17x2x56)
     &  +scints(4,d1x7x2,  0)*coeff(4,d1x7x2)
     &  +scints(3,c56x34,0)*coeff(3,c56x34)
     &  +scints(3,c17x34,0)*coeff(3,c17x34)
     &  +scints(3,c27x56,0)*coeff(3,c27x56)
     &  +scints(2,b17 ,0)*coeff(2,b17)
     &  +scints(2,b134,0)*coeff(2,b134)
     &  +scints(2,b256,0)*coeff(2,b256)
     &  +scints(2,b127,0)*coeff(2,b127)
     &  +scints(2,b56 ,0)*coeff(2,b56)
     &  +scints(2,b34 ,0)*coeff(2,b34)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        A7a_lc=
     &  +scints(4,d7x1x56, 0)*coeff(4,d7x1x34)
     &  +scints(4,d1x27x34,0)*coeff(4,d27x1x34)
     &  +scints(4,d7x2x34, 0)*coeff(4,d7x2x56)
     &  +scints(4,d2x17x56,0)*coeff(4,d17x2x56)
     &  +scints(4,d1x7x2,  0)*coeff(4,d1x7x2)
     &  +scints(3,c56x34,0)*coeff(3,c56x34)
     &  +scints(3,c56x17,0)*coeff(3,c17x34)
     &  +scints(3,c34x27,0)*coeff(3,c27x56)
     &  +scints(2,b17 ,0)*coeff(2,b17)
     &  +scints(2,b156,0)*coeff(2,b134)
     &  +scints(2,b234,0)*coeff(2,b256)
     &  +scints(2,b127,0)*coeff(2,b127)
     &  +scints(2,b34 ,0)*coeff(2,b56)
     &  +scints(2,b56 ,0)*coeff(2,b34)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        A7a_lc=
     &  +scints(4,d7x2x34, 0)*coeff(4,d7x1x34)
     &  +scints(4,d2x17x56,0)*coeff(4,d27x1x34)
     &  +scints(4,d7x1x56, 0)*coeff(4,d7x2x56)
     &  +scints(4,d1x27x34,0)*coeff(4,d17x2x56)
     &  +scints(4,d1x7x2,  0)*coeff(4,d1x7x2)
     &  +scints(3,c56x34,0)*coeff(3,c56x34)
     &  +scints(3,c34x27,0)*coeff(3,c17x34)
     &  +scints(3,c56x17,0)*coeff(3,c27x56)
     &  +scints(2,b27 ,0)*coeff(2,b17)
     &  +scints(2,b234,0)*coeff(2,b134)
     &  +scints(2,b156,0)*coeff(2,b256)
     &  +scints(2,b127,0)*coeff(2,b127)
     &  +scints(2,b56 ,0)*coeff(2,b56)
     &  +scints(2,b34 ,0)*coeff(2,b34)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        A7a_lc=
     &  +scints(4,d7x2x56, 0)*coeff(4,d7x1x34)
     &  +scints(4,d17x2x56,0)*coeff(4,d27x1x34)
     &  +scints(4,d7x1x34, 0)*coeff(4,d7x2x56)
     &  +scints(4,d27x1x34,0)*coeff(4,d17x2x56)
     &  +scints(4,d1x7x2,  0)*coeff(4,d1x7x2)
     &  +scints(3,c56x34,0)*coeff(3,c56x34)
     &  +scints(3,c27x56,0)*coeff(3,c17x34)
     &  +scints(3,c17x34,0)*coeff(3,c27x56)
     &  +scints(2,b27 ,0)*coeff(2,b17)
     &  +scints(2,b256,0)*coeff(2,b134)
     &  +scints(2,b134,0)*coeff(2,b256)
     &  +scints(2,b127,0)*coeff(2,b127)
     &  +scints(2,b34 ,0)*coeff(2,b56)
     &  +scints(2,b56 ,0)*coeff(2,b34)
     &  +coeff(0,irat)
      endif

c Contribution from known pole terms
      A7a_lc=A7a_lc+Alo*(
     & +epinv*epinv2*(-2._dp)
     & +epinv*(-1.5_dp-(lnrat(musq,-s(j1,j7))+lnrat(musq,-s(j2,j7))))
     & -0.5_dp*(lnrat(musq,-s(j1,j7))**2+lnrat(musq,-s(j2,j7))**2))

      if (swap12 .eqv. .true.) then
        Alo=-Alo
        A7a_lc=-A7a_lc
      endif

      return
      end
