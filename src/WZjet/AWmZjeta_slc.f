!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWmZjeta_slc(swap12,swapVV,p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz,scints,coeff,Alo,A7a_slc)
c--- MCFM notation
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c--- correspoding to Tania notation
c---   u(6) ubar(1) nu(3) e+(2) e-(5) nubar(4) g(7)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'WWjetlabels.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer j1,j2,j6,j5,j3,j4,j7
      complex(dp):: Alo,A7treea,A7a_slc,lnrat
      complex(dp)::
     & xd1x27x34sl,xd1x2x34sl,xd1x2x7sl,xd1x56x7sl,xd1x7x56sl,
     & xd2x1x56sl,xd2x1x7sl,xd2x34x7sl,xd2x7x34sl,xd17x2x34sl,
     & xd7x34x12sl,xd7x12x34sl,xd7x12x56sl,
     & bub347sl,bub127sl,bub567sl,bub156sl,bub234sl,bub17sl,
     & bub12sl,bub56sl
      real(dp):: p(mxpart,4)
      logical:: swap12,swapVV,swapz

      Alo=A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)

c      s156=s(j1,j5)+s(j1,j6)+s(j5,j6)
c      s127=s(j1,j2)+s(j1,j7)+s(j2,j7)

      if (swap12 .neqv. swapVV) then
c        call parser('runWmZ','SLau',helname,KCD,KCC,KCB,rattot,tot)
      else
c        call parser('runWmZ','SLad',helname,KCD,KCC,KCB,rattot,tot)
      endif
c      include 'runSL.f'

      coeff(4,d27x1x34sl)=xd1x27x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x2x56sl)=xd1x2x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x2x7sl)=xd1x2x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x34x7sl)=xd1x56x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x7x34sl)=xd1x7x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x1x34sl)=xd2x1x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x1x7sl)=xd2x1x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x56x7sl)=xd2x34x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x7x56sl)=xd2x7x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d17x2x56sl)=xd17x2x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d7x34x12sl)=xd7x34x12sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d7x12x56sl)=xd7x12x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d7x12x34sl)=xd7x12x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)

c--- Subleading color bubbles
      coeff(2,b567sl)=bub567sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b127sl)=bub127sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b347sl)=bub347sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b134sl)=bub156sl(j1,j2,j6,j5,j4,j3,j7,za,zb)+coeff(2,b134)
      coeff(2,b256sl)=bub234sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)+coeff(2,b256)
      coeff(2,b17sl)=bub17sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)+coeff(2,b17)
      coeff(2,b12sl)=bub12sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b34sl)=bub56sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
c--- remaining bubble from pole structure
      coeff(2,b56sl)=-1.5d0*Alo-coeff(2,b134sl)-coeff(2,b256sl)
     & -coeff(2,b34sl)-coeff(2,b17sl)-coeff(2,b127sl)-coeff(2,b347sl)
     & -coeff(2,b567sl)-coeff(2,b12sl)

c Now assemble total
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        A7a_slc=
     &  +scints(4,d2x56x7sl ,0)*coeff(4,d2x56x7sl )
     &  +scints(4,d2x7x56sl ,0)*coeff(4,d2x7x56sl )
     &  +scints(4,d17x2x56sl,0)*coeff(4,d17x2x56sl)
     &  +scints(4,d1x34x7sl ,0)*coeff(4,d1x34x7sl )
     &  +scints(4,d1x7x34sl ,0)*coeff(4,d1x7x34sl )
     &  +scints(4,d27x1x34sl,0)*coeff(4,d27x1x34sl)
     &  +scints(4,d2x1x7sl  ,0)*coeff(4,d2x1x7sl  )
     &  +scints(4,d1x2x7sl  ,0)*coeff(4,d1x2x7sl  )
     &  +scints(4,d2x1x34sl ,0)*coeff(4,d2x1x34sl )
     &  +scints(4,d7x34x12sl,0)*coeff(4,d7x34x12sl)
     &  +scints(4,d1x2x56sl ,0)*coeff(4,d1x2x56sl )
     &  +scints(4,d7x12x34sl,0)*coeff(4,d7x12x34sl)
     &  +scints(4,d7x12x56sl,0)*coeff(4,d7x12x56sl)
     &  +scints(3,c56x34sl,0)*coeff(3,c56x34sl)
     &  +scints(3,c17x34sl,0)*coeff(3,c17x34sl)
     &  +scints(3,c27x56sl,0)*coeff(3,c27x56sl)
     &  +scints(3,c12x34sl,0)*coeff(3,c12x34sl)
     &  +scints(3,c12x56sl,0)*coeff(3,c12x56sl)
     &  +scints(2,b347sl,0)*coeff(2,b347sl)
     &  +scints(2,b127sl,0)*coeff(2,b127sl)
     &  +scints(2,b567sl,0)*coeff(2,b567sl)
     &  +scints(2,b134sl,0)*coeff(2,b134sl)
     &  +scints(2,b256sl,0)*coeff(2,b256sl)
     &  +scints(2,b17sl ,0)*coeff(2,b17sl)
     &  +scints(2,b12sl ,0)*coeff(2,b12sl)
     &  +scints(2,b56sl ,0)*coeff(2,b56sl)
     &  +scints(2,b34sl ,0)*coeff(2,b34sl)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        A7a_slc=
     &  +scints(4,d2x34x7sl ,0)*coeff(4,d2x56x7sl )
     &  +scints(4,d2x7x34sl ,0)*coeff(4,d2x7x56sl )
     &  +scints(4,d17x2x34sl,0)*coeff(4,d17x2x56sl)
     &  +scints(4,d1x56x7sl ,0)*coeff(4,d1x34x7sl )
     &  +scints(4,d1x7x56sl ,0)*coeff(4,d1x7x34sl )
     &  +scints(4,d1x27x34sl,0)*coeff(4,d27x1x34sl)
     &  +scints(4,d2x1x7sl  ,0)*coeff(4,d2x1x7sl  )
     &  +scints(4,d1x2x7sl  ,0)*coeff(4,d1x2x7sl  )
     &  +scints(4,d2x1x56sl ,0)*coeff(4,d2x1x34sl )
     &  +scints(4,d7x34x12sl,0)*coeff(4,d7x34x12sl)
     &  +scints(4,d1x2x34sl ,0)*coeff(4,d1x2x56sl )
     &  +scints(4,d7x12x56sl,0)*coeff(4,d7x12x34sl)
     &  +scints(4,d7x12x34sl,0)*coeff(4,d7x12x56sl)
     &  +scints(3,c56x34sl,0)*coeff(3,c56x34sl)
     &  +scints(3,c56x17sl,0)*coeff(3,c17x34sl)
     &  +scints(3,c34x27sl,0)*coeff(3,c27x56sl)
     &  +scints(3,c12x56sl,0)*coeff(3,c12x34sl)
     &  +scints(3,c12x34sl,0)*coeff(3,c12x56sl)
     &  +scints(2,b567sl,0)*coeff(2,b347sl)
     &  +scints(2,b127sl,0)*coeff(2,b127sl)
     &  +scints(2,b347sl,0)*coeff(2,b567sl)
     &  +scints(2,b156sl,0)*coeff(2,b134sl)
     &  +scints(2,b234sl,0)*coeff(2,b256sl)
     &  +scints(2,b17sl ,0)*coeff(2,b17sl)
     &  +scints(2,b12sl ,0)*coeff(2,b12sl)
     &  +scints(2,b34sl ,0)*coeff(2,b56sl)
     &  +scints(2,b56sl ,0)*coeff(2,b34sl)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        A7a_slc=
     &  +scints(4,d1x56x7sl ,0)*coeff(4,d2x56x7sl )
     &  +scints(4,d1x7x56sl ,0)*coeff(4,d2x7x56sl )
     &  +scints(4,d1x27x34sl,0)*coeff(4,d17x2x56sl)
     &  +scints(4,d2x34x7sl ,0)*coeff(4,d1x34x7sl )
     &  +scints(4,d2x7x34sl ,0)*coeff(4,d1x7x34sl )
     &  +scints(4,d17x2x34sl,0)*coeff(4,d27x1x34sl)
     &  +scints(4,d1x2x7sl  ,0)*coeff(4,d2x1x7sl  )
     &  +scints(4,d2x1x7sl  ,0)*coeff(4,d1x2x7sl  )
     &  +scints(4,d1x2x34sl ,0)*coeff(4,d2x1x34sl )
     &  +scints(4,d7x34x12sl,0)*coeff(4,d7x34x12sl)
     &  +scints(4,d2x1x56sl ,0)*coeff(4,d1x2x56sl )
     &  +scints(4,d7x12x34sl,0)*coeff(4,d7x12x34sl)
     &  +scints(4,d7x12x56sl,0)*coeff(4,d7x12x56sl)
     &  +scints(3,c56x34sl,0)*coeff(3,c56x34sl)
     &  +scints(3,c34x27sl,0)*coeff(3,c17x34sl)
     &  +scints(3,c17x56sl,0)*coeff(3,c27x56sl)
     &  +scints(3,c12x34sl,0)*coeff(3,c12x34sl)
     &  +scints(3,c12x56sl,0)*coeff(3,c12x56sl)
     &  +scints(2,b347sl,0)*coeff(2,b347sl)
     &  +scints(2,b127sl,0)*coeff(2,b127sl)
     &  +scints(2,b567sl,0)*coeff(2,b567sl)
     &  +scints(2,b234sl,0)*coeff(2,b134sl)
     &  +scints(2,b156sl,0)*coeff(2,b256sl)
     &  +scints(2,b27sl ,0)*coeff(2,b17sl)
     &  +scints(2,b12sl ,0)*coeff(2,b12sl)
     &  +scints(2,b56sl ,0)*coeff(2,b56sl)
     &  +scints(2,b34sl ,0)*coeff(2,b34sl)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        A7a_slc=
     &  +scints(4,d1x34x7sl ,0)*coeff(4,d2x56x7sl )
     &  +scints(4,d1x7x34sl ,0)*coeff(4,d2x7x56sl )
     &  +scints(4,d27x1x34sl,0)*coeff(4,d17x2x56sl)
     &  +scints(4,d2x56x7sl ,0)*coeff(4,d1x34x7sl )
     &  +scints(4,d2x7x56sl ,0)*coeff(4,d1x7x34sl )
     &  +scints(4,d17x2x56sl,0)*coeff(4,d27x1x34sl)
     &  +scints(4,d1x2x7sl  ,0)*coeff(4,d2x1x7sl  )
     &  +scints(4,d2x1x7sl  ,0)*coeff(4,d1x2x7sl  )
     &  +scints(4,d1x2x56sl ,0)*coeff(4,d2x1x34sl )
     &  +scints(4,d7x34x12sl,0)*coeff(4,d7x34x12sl)
     &  +scints(4,d2x1x34sl ,0)*coeff(4,d1x2x56sl )
     &  +scints(4,d7x12x56sl,0)*coeff(4,d7x12x34sl)
     &  +scints(4,d7x12x34sl,0)*coeff(4,d7x12x56sl)
     &  +scints(3,c56x34sl,0)*coeff(3,c56x34sl)
     &  +scints(3,c27x56sl,0)*coeff(3,c17x34sl)
     &  +scints(3,c17x34sl,0)*coeff(3,c27x56sl)
     &  +scints(3,c12x56sl,0)*coeff(3,c12x34sl)
     &  +scints(3,c12x34sl,0)*coeff(3,c12x56sl)
     &  +scints(2,b567sl,0)*coeff(2,b347sl)
     &  +scints(2,b127sl,0)*coeff(2,b127sl)
     &  +scints(2,b347sl,0)*coeff(2,b567sl)
     &  +scints(2,b256sl,0)*coeff(2,b134sl)
     &  +scints(2,b134sl,0)*coeff(2,b256sl)
     &  +scints(2,b27sl ,0)*coeff(2,b17sl)
     &  +scints(2,b12sl ,0)*coeff(2,b12sl)
     &  +scints(2,b34sl ,0)*coeff(2,b56sl)
     &  +scints(2,b56sl ,0)*coeff(2,b34sl)
     &  +coeff(0,iratsl)
      endif

c Contribution from known pole terms
      A7a_slc=A7a_slc+Alo*(
     & +epinv*epinv2*(-1._dp)
     & +epinv*(-1.5_dp-lnrat(musq,-s(j1,j2)))
     & -0.5_dp*lnrat(musq,-s(j1,j2))**2)

      if (swap12 .eqv. .true.) then
        Alo=-Alo
        A7a_slc=-A7a_slc
      endif

      return
      end


