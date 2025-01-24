!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjeta_slc(helname,swap12,swapVV,p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz,scints,coeff,Alo,A7a_slc)
c--- MCFM notation
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c--- corresponding to Tania notation
c---   u(6) ubar(1) nu(3) e+(2) e-(5) nubar(4) g(7)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'WWjetlabels.f'
      include 'KCdef.f'
      include 'verbose.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j6,j5,j3,j4,j7
      complex(dp):: Alo,Aloswap,A7treea,A7a_slc
      complex(dp)::coeff3456(10:14),coeff6543(10:14),
     & xd1x27x34sl,xd1x2x34sl,xd1x2x7sl,xd1x56x7sl,xd1x7x56sl,
     & xd2x1x56sl,xd2x1x7sl,xd2x34x7sl,xd2x7x34sl,xd17x2x34sl,
     & xd7x34x12sl,xd7x12x34sl,xd7x12x56sl,
     & bub347sl,bub127sl,bub567sl,bub156sl,bub234sl,bub17sl,
     & bub12sl,bub56sl,bub17,b34extra,b56extra,lnrat
      real(dp):: p(mxpart,4)
      character(len=4):: helname
      logical:: swap12,swapVV,swapz

      Alo=A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)
      Aloswap=A7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)


c      if (verbose) call parser('runZZ','SLa',helname,KCD,KCC,KCB,rattot,tot)

c      include 'runZZSLa.f'
c      include 'runZZSLa.txt'

      coeff(4,d27x1x34sl)=xd1x27x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x27x34sl)=xd1x27x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d1x2x56sl)=xd1x2x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x2x34sl)=xd1x2x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d1x34x7sl)=xd1x56x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x56x7sl)=xd1x56x7sl(j1,j2,j3,j4,j6,j5,j7,za,zb)
      coeff(4,d1x2x7sl)=xd1x2x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
     &                 +xd1x2x7sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d2x1x7sl)=xd2x1x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
     &                 +xd2x1x7sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d1x7x34sl)=xd1x7x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x7x56sl)=xd1x7x56sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d2x1x34sl)=xd2x1x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x1x56sl)=xd2x1x56sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d2x56x7sl)=xd2x34x7sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x34x7sl)=xd2x34x7sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d2x7x56sl)=xd2x7x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d2x7x34sl)=xd2x7x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d17x2x56sl)=xd17x2x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d17x2x34sl)=xd17x2x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d7x34x12sl)=
     & +xd7x34x12sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
     & +xd7x34x12sl(j1,j2,j3,j4,j6,j5,j7,za,zb)
      coeff(4,d12x7x34sl)=
     &  xd7x12x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb)
     & +xd7x12x56sl(j1,j2,j3,j4,j6,j5,j7,za,zb)
      coeff(4,d7x12x34sl)=
     &  xd7x12x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb)
     & +xd7x12x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb)

      if (verbose) then
      write(6,*)
      write(6,*) 'SL:Comparison to KCheck'
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
      write(6,*)'d1x27x34sl ',KCD(d1x27x34sl )/coeff(4,d1x27x34sl )
      write(6,*)'d27x1x34sl',KCD(d27x1x34sl)/coeff(4,d27x1x34sl)
      write(6,*)'d1x2x34sl  ',KCD(d1x2x34sl  )/coeff(4,d1x2x34sl  )
      write(6,*)'d1x2x56sl',KCD(d1x2x56sl)/coeff(4,d1x2x56sl)
      write(6,*)'d1x56x7sl  ',KCD(d1x56x7sl  )/coeff(4,d1x56x7sl  )
      write(6,*)'d1x34x7sl',KCD(d1x34x7sl)/coeff(4,d1x34x7sl)

      write(6,*)'d1x2x7sl',KCD(d1x2x7sl)/coeff(4,d1x2x7sl)
      write(6,*)'d2x1x7sl',KCD(d2x1x7sl)/coeff(4,d2x1x7sl)
      write(6,*)'d1x7x56sl  ',KCD(d1x7x56sl  )/coeff(4,d1x7x56sl  )
      write(6,*)'d1x7x34sl',KCD(d1x7x34sl)/coeff(4,d1x7x34sl)
      write(6,*)'d2x1x56sl',KCD(d2x1x56sl)/coeff(4,d2x1x56sl)
      write(6,*)'d2x1x34sl',KCD(d2x1x34sl)/coeff(4,d2x1x34sl)
      write(6,*)'d2x34x7sl  ',KCD(d2x34x7sl  )/coeff(4,d2x34x7sl  )
      write(6,*)'d2x56x7sl  ',KCD(d2x56x7sl  )/coeff(4,d2x56x7sl  )
      write(6,*)'d2x7x34sl  ',KCD(d2x7x34sl  )/coeff(4,d2x7x34sl  )
      write(6,*)'d2x7x56sl  ',KCD(d2x7x56sl  )/coeff(4,d2x7x56sl  )
      write(6,*)'d17x2x34sl ',KCD(d17x2x34sl )/coeff(4,d17x2x34sl )
      write(6,*)'d17x2x56sl',KCD(d17x2x56sl)/coeff(4,d17x2x56sl)
      write(6,*)'d7x34x12sl',KCD(d7x34x12sl)/coeff(4,d7x34x12sl)

      write(6,*) 'd7x12x34sl',KCD(d7x12x34sl)/coeff(4,d7x12x34sl)
      write(6,*) 'd12x7x34sl',KCD(d12x7x34sl)/coeff(4,d12x7x34sl)

      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
      write(6,*)'d1x27x34sl ',KCD(d1x27x34sl )/coeff(4,d27x1x34sl )
      write(6,*)'d27x1x34sl',KCD(d27x1x34sl)/coeff(4,d1x27x34sl)
      write(6,*)'d1x2x34sl  ',KCD(d1x2x34sl  )/coeff(4,d1x2x56sl  )
      write(6,*)'d1x2x56sl',KCD(d1x2x56sl)/coeff(4,d1x2x34sl)
      write(6,*)'d1x56x7sl  ',KCD(d1x56x7sl  )/coeff(4,d1x34x7sl  )
      write(6,*)'d1x34x7sl',KCD(d1x34x7sl)/coeff(4,d1x56x7sl)

      write(6,*)'d1x2x7sl',KCD(d1x2x7sl)/coeff(4,d1x2x7sl)
      write(6,*)'d2x1x7sl',KCD(d2x1x7sl)/coeff(4,d2x1x7sl)
      write(6,*)'d1x7x56sl  ',KCD(d1x7x56sl  )/coeff(4,d1x7x34sl  )
      write(6,*)'d1x7x34sl',KCD(d1x7x34sl)/coeff(4,d1x7x56sl)
      write(6,*)'d2x1x56sl',KCD(d2x1x56sl)/coeff(4,d2x1x34sl)
      write(6,*)'d2x1x34sl',KCD(d2x1x34sl)/coeff(4,d2x1x56sl)
      write(6,*)'d2x34x7sl  ',KCD(d2x34x7sl  )/coeff(4,d2x56x7sl  )
      write(6,*)'d2x56x7sl  ',KCD(d2x56x7sl  )/coeff(4,d2x34x7sl  )
      write(6,*)'d2x7x34sl  ',KCD(d2x7x34sl  )/coeff(4,d2x7x56sl  )
      write(6,*)'d2x7x56sl  ',KCD(d2x7x56sl  )/coeff(4,d2x7x34sl  )
      write(6,*)'d17x2x34sl ',KCD(d17x2x34sl )/coeff(4,d17x2x56sl )
      write(6,*)'d17x2x56sl',KCD(d17x2x56sl)/coeff(4,d17x2x34sl)
      write(6,*)'d7x34x12sl',KCD(d7x34x12sl)/coeff(4,d7x34x12sl)

      write(6,*) 'd7x12x34sl',KCD(d7x12x34sl)/coeff(4,d12x7x34sl)
      write(6,*) 'd12x7x34sl',KCD(d12x7x34sl)/coeff(4,d7x12x34sl)

      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
      write(6,*)'d1x27x34sl ',KCD(d1x27x34sl )/coeff(4,d17x2x56sl )
      write(6,*)'d27x1x34sl',KCD(d27x1x34sl)/coeff(4,d17x2x34sl)
      write(6,*)'d1x2x34sl  ',KCD(d1x2x34sl  )/coeff(4,d2x1x34sl  )
      write(6,*)'d1x2x56sl',KCD(d1x2x56sl)/coeff(4,d2x1x56sl)
      write(6,*)'d1x56x7sl  ',KCD(d1x56x7sl  )/coeff(4,d2x56x7sl  )
      write(6,*)'d1x34x7sl',KCD(d1x34x7sl)/coeff(4,d2x34x7sl)

      write(6,*)'d1x2x7sl',KCD(d1x2x7sl)/coeff(4,d2x1x7sl)
      write(6,*)'d2x1x7sl',KCD(d2x1x7sl)/coeff(4,d1x2x7sl)
      write(6,*)'d1x7x56sl  ',KCD(d1x7x56sl  )/coeff(4,d2x7x56sl  )
      write(6,*)'d1x7x34sl',KCD(d1x7x34sl)/coeff(4,d2x7x34sl)
      write(6,*)'d2x1x56sl',KCD(d2x1x56sl)/coeff(4,d1x2x56sl)
      write(6,*)'d2x1x34sl',KCD(d2x1x34sl)/coeff(4,d1x2x34sl)
      write(6,*)'d2x34x7sl  ',KCD(d2x34x7sl  )/coeff(4,d1x34x7sl  )
      write(6,*)'d2x56x7sl  ',KCD(d2x56x7sl  )/coeff(4,d1x56x7sl  )
      write(6,*)'d2x7x34sl  ',KCD(d2x7x34sl  )/coeff(4,d1x7x34sl  )
      write(6,*)'d2x7x56sl  ',KCD(d2x7x56sl  )/coeff(4,d1x7x56sl  )
      write(6,*)'d17x2x34sl ',KCD(d17x2x34sl )/coeff(4,d27x1x34sl )
      write(6,*)'d17x2x56sl',KCD(d17x2x56sl)/coeff(4,d1x27x34sl)
      write(6,*)'d7x34x12sl',KCD(d7x34x12sl)/coeff(4,d7x34x12sl)

      write(6,*) 'd7x12x34sl',KCD(d7x12x34sl)/coeff(4,d7x12x34sl)
      write(6,*) 'd12x7x34sl',KCD(d12x7x34sl)/coeff(4,d12x7x34sl)

      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
      write(6,*)'d1x27x34sl ',KCD(d1x27x34sl )/coeff(4,d17x2x34sl )
      write(6,*)'d27x1x34sl',KCD(d27x1x34sl)/coeff(4,d17x2x56sl)
      write(6,*)'d1x2x34sl  ',KCD(d1x2x34sl  )/coeff(4,d2x1x56sl  )
      write(6,*)'d1x2x56sl',KCD(d1x2x56sl)/coeff(4,d2x1x34sl)
      write(6,*)'d1x56x7sl  ',KCD(d1x56x7sl  )/coeff(4,d2x34x7sl  )
      write(6,*)'d1x34x7sl',KCD(d1x34x7sl)/coeff(4,d2x56x7sl)

      write(6,*)'d1x2x7sl',KCD(d1x2x7sl)/coeff(4,d2x1x7sl)
      write(6,*)'d2x1x7sl',KCD(d2x1x7sl)/coeff(4,d1x2x7sl)
      write(6,*)'d1x7x56sl  ',KCD(d1x7x56sl  )/coeff(4,d2x7x34sl  )
      write(6,*)'d1x7x34sl',KCD(d1x7x34sl)/coeff(4,d2x7x56sl)
      write(6,*)'d2x1x56sl',KCD(d2x1x56sl)/coeff(4,d1x2x34sl)
      write(6,*)'d2x1x34sl',KCD(d2x1x34sl)/coeff(4,d1x2x56sl)
      write(6,*)'d2x34x7sl  ',KCD(d2x34x7sl  )/coeff(4,d1x56x7sl  )
      write(6,*)'d2x56x7sl  ',KCD(d2x56x7sl  )/coeff(4,d1x34x7sl  )
      write(6,*)'d2x7x34sl  ',KCD(d2x7x34sl  )/coeff(4,d1x7x56sl  )
      write(6,*)'d2x7x56sl  ',KCD(d2x7x56sl  )/coeff(4,d1x7x34sl  )
      write(6,*)'d17x2x34sl ',KCD(d17x2x34sl )/coeff(4,d1x27x34sl )
      write(6,*)'d17x2x56sl',KCD(d17x2x56sl)/coeff(4,d27x1x34sl)
      write(6,*)'d7x34x12sl',KCD(d7x34x12sl)/coeff(4,d7x34x12sl)

      write(6,*) 'd7x12x34sl',KCD(d7x12x34sl)/coeff(4,d12x7x34sl)
      write(6,*) 'd12x7x34sl',KCD(d12x7x34sl)/coeff(4,d7x12x34sl)

      endif

      write(6,*)

      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'c56x17',KCC(c56x17sl)/coeff(3,c17x56sl)
        write(6,*)'c17x34',KCC(c17x34sl)/coeff(3,c17x34sl)
        write(6,*)'c27x56',KCC(c27x56sl)/coeff(3,c27x56sl)
        write(6,*)'c34x27',KCC(c34x27sl)/coeff(3,c34x27sl)
        write(6,*)'c56x34',KCC(c56x34sl)/coeff(3,c34x56sl)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'c56x17',KCC(c56x17sl)/coeff(3,c17x34sl)
        write(6,*)'c17x34',KCC(c17x34sl)/coeff(3,c56x17sl)
        write(6,*)'c27x56',KCC(c27x56sl)/coeff(3,c34x27sl)
        write(6,*)'c34x27',KCC(c34x27sl)/coeff(3,c27x56sl)
        write(6,*)'c56x34',KCC(c56x34sl)/coeff(3,c56x34sl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'c56x17',KCC(c27x56sl)/coeff(3,c56x17sl)
        write(6,*)'c17x34',KCC(c34x27sl)/coeff(3,c17x34sl)
        write(6,*)'c27x56',KCC(c56x17sl)/coeff(3,c27x56sl)
        write(6,*)'c34x27',KCC(c17x34sl)/coeff(3,c34x27sl)
        write(6,*)'c56x34',KCC(c56x34sl)/coeff(3,c56x34sl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'c56x17',KCC(c27x56sl)/coeff(3,c17x34sl)
        write(6,*)'c17x34',KCC(c34x27sl)/coeff(3,c56x17sl)
        write(6,*)'c27x56',KCC(c56x17sl)/coeff(3,c34x27sl)
        write(6,*)'c34x27',KCC(c17x34sl)/coeff(3,c27x56sl)
        write(6,*)'c56x34',KCC(c56x34sl)/coeff(3,c56x34sl)
      endif

 !     write(6,*) 'Results from leading colour'

c      write(6,*) 'c56x17sl',coeff(3,c56x17sl)/KCC(c17x56sl)
c      write(6,*) 'c34x17sl',coeff(3,c17x34sl)/KCC(c17x34sl)
c      write(6,*) 'c27x56sl',coeff(3,c27x56sl)/KCC(c27x56sl)
c      write(6,*) 'c34x27sl',coeff(3,c34x27sl)/KCC(c34x27sl)
c      write(6,*) 'c56x34sl',coeff(3,c56x34sl)/KCC(c34x56sl)

c      write(6,*) 'c12x34sl',coeff(3,c12x34sl)/KCC(c12x34sl)
c      write(6,*) 'c12x56sl',coeff(3,c12x56sl)/KCC(c12x56sl)
      write(6,*)
      endif
c--- Subleading color bubbles
      coeff(2,b134sl)=bub156sl(j1,j2,j6,j5,j4,j3,j7,za,zb)+coeff(2,b134)
      coeff(2,b156sl)=bub156sl(j1,j2,j3,j4,j5,j6,j7,za,zb)+coeff(2,b156)

      coeff(2,b256sl)=bub234sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)+coeff(2,b256)
      coeff(2,b234sl)=bub234sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)+coeff(2,b234)

      coeff(2,b34sl)=bub56sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff(2,b56sl)=bub56sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)

      coeff3456(b567sl)=bub567sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b567sl)=bub567sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff3456(b347sl)=bub347sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b347sl)=bub347sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff3456(b12sl)=bub12sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
c     b12sl antisymmetric under flip (3456)<-->(6543)
      coeff6543(b12sl)=-coeff3456(b12sl)
      coeff3456(b127sl)=bub127sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b127sl)=bub127sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff3456(b17sl)=bub17sl(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b17sl)=bub17sl(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)

      coeff(2,b567sl)=coeff3456(b567sl)+coeff6543(b347sl)
      coeff(2,b347sl)=coeff6543(b567sl)+coeff3456(b347sl)
      coeff(2,b127sl)=coeff3456(b127sl)+coeff6543(b127sl)
      coeff(2,b17sl)=coeff3456(b17sl)+coeff6543(b17sl)+coeff(2,b17)
c     b12sl not needed in this case: antisymmetric under flip (3456)<-->(6543)
      coeff(2,b12sl)=czip

c--- remaining bubble bits from pole structure
      b56extra=-1.5d0*Alo
     & -coeff(2,b134sl)-coeff(2,b256sl)-coeff(2,b34sl)
     & -coeff3456(b12sl)
     & -coeff3456(b17sl)-bub17(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
     & -coeff3456(b127sl)-coeff3456(b347sl)-coeff3456(b567sl)
      b34extra=-1.5d0*Aloswap
     & -coeff(2,b156sl)-coeff(2,b234sl)-coeff(2,b56sl)
     & -coeff6543(b12sl)
     & -coeff6543(b17sl)-bub17(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
     & -coeff6543(b127sl)-coeff6543(b347sl)-coeff6543(b567sl)
      coeff(2,b56sl)=coeff(2,b56sl)+b56extra
      coeff(2,b34sl)=coeff(2,b34sl)+b34extra

      if (verbose) then
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        write(6,*) 'b127sl',coeff(2,b127sl)/KCB(b127sl)
        write(6,*) 'b156sl',coeff(2,b156sl)/KCB(b156sl)
        write(6,*) 'b134sl',coeff(2,b134sl)/KCB(b134sl)
        write(6,*) 'b56sl ',coeff(2,b56sl) /KCB(b56sl )
        write(6,*) 'b34sl ',coeff(2,b34sl) /KCB(b34sl )
        write(6,*) 'b234sl',coeff(2,b234sl)/KCB(b234sl)
        write(6,*) 'b256sl',coeff(2,b256sl)/KCB(b256sl)
        write(6,*) 'b567sl',coeff(2,b567sl)/KCB(b567sl)
        write(6,*) 'b347sl',coeff(2,b347sl)/KCB(b347sl)
        write(6,*) 'b17sl ',coeff(2,b17sl) /KCB(b17sl )
        write(6,*) 'b12sl ',coeff(2,b12sl) ,KCB(b12sl )
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        write(6,*) 'b127sl',coeff(2,b127sl)/KCB(b127sl)
        write(6,*) 'b156sl',coeff(2,b156sl)/KCB(b134sl)
        write(6,*) 'b134sl',coeff(2,b134sl)/KCB(b156sl)
        write(6,*) 'b56sl ',coeff(2,b56sl) /KCB(b34sl )
        write(6,*) 'b34sl ',coeff(2,b34sl) /KCB(b56sl )
        write(6,*) 'b234sl',coeff(2,b234sl)/KCB(b256sl)
        write(6,*) 'b256sl',coeff(2,b256sl)/KCB(b234sl)
        write(6,*) 'b567sl',coeff(2,b567sl)/KCB(b347sl)
        write(6,*) 'b347sl',coeff(2,b347sl)/KCB(b567sl)
        write(6,*) 'b17sl ',coeff(2,b17sl) /KCB(b17sl )
        write(6,*) 'b12sl ',coeff(2,b12sl) ,KCB(b12sl )
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        write(6,*) 'b127sl',coeff(2,b127sl)/KCB(b127sl)
        write(6,*) 'b156sl',coeff(2,b156sl)/KCB(b256sl)
        write(6,*) 'b134sl',coeff(2,b134sl)/KCB(b234sl)
        write(6,*) 'b56sl ',coeff(2,b56sl) /KCB(b56sl )
        write(6,*) 'b34sl ',coeff(2,b34sl) /KCB(b34sl )
        write(6,*) 'b234sl',coeff(2,b234sl)/KCB(b134sl)
        write(6,*) 'b256sl',coeff(2,b256sl)/KCB(b156sl)
        write(6,*) 'b567sl',coeff(2,b567sl)/KCB(b567sl)
        write(6,*) 'b347sl',coeff(2,b347sl)/KCB(b347sl)
        write(6,*) 'b17sl ',coeff(2,b17sl) /KCB(b27sl )
        write(6,*) 'b12sl ',coeff(2,b12sl) ,KCB(b12sl )
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        write(6,*) 'b127sl',coeff(2,b127sl)/KCB(b127sl)
        write(6,*) 'b156sl',coeff(2,b156sl)/KCB(b234sl)
        write(6,*) 'b134sl',coeff(2,b134sl)/KCB(b256sl)
        write(6,*) 'b56sl ',coeff(2,b56sl) /KCB(b34sl )
        write(6,*) 'b34sl ',coeff(2,b34sl) /KCB(b56sl )
        write(6,*) 'b234sl',coeff(2,b234sl)/KCB(b156sl)
        write(6,*) 'b256sl',coeff(2,b256sl)/KCB(b134sl)
        write(6,*) 'b567sl',coeff(2,b567sl)/KCB(b347sl)
        write(6,*) 'b347sl',coeff(2,b347sl)/KCB(b567sl)
        write(6,*) 'b17sl ',coeff(2,b17sl) /KCB(b27sl )
        write(6,*) 'b12sl ',coeff(2,b12sl) ,KCB(b12sl )
      endif
      endif

      if (verbose) then
        write(6,*)
        write(6,*) 'rational',coeff(0,iratsl)/rattot
      endif

c Now assemble total
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        A7a_slc=
     &  +coeff(4,d1x27x34sl)*scints(4,d1x27x34sl,0)
     &  +coeff(4,d27x1x34sl)*scints(4,d27x1x34sl,0)
     &  +coeff(4,d1x2x34sl )*scints(4,d1x2x34sl ,0)
     &  +coeff(4,d1x2x56sl )*scints(4,d1x2x56sl ,0)
     &  +coeff(4,d1x56x7sl )*scints(4,d1x56x7sl ,0)
     &  +coeff(4,d1x34x7sl )*scints(4,d1x34x7sl ,0)
     &  +coeff(4,d1x2x7sl  )*scints(4,d1x2x7sl  ,0)
     &  +coeff(4,d2x1x7sl  )*scints(4,d2x1x7sl  ,0)
     &  +coeff(4,d1x7x56sl )*scints(4,d1x7x56sl ,0)
     &  +coeff(4,d1x7x34sl )*scints(4,d1x7x34sl ,0)
     &  +coeff(4,d2x1x56sl )*scints(4,d2x1x56sl ,0)
     &  +coeff(4,d2x1x34sl )*scints(4,d2x1x34sl ,0)
     &  +coeff(4,d2x34x7sl )*scints(4,d2x34x7sl ,0)
     &  +coeff(4,d2x56x7sl )*scints(4,d2x56x7sl ,0)
     &  +coeff(4,d2x7x34sl )*scints(4,d2x7x34sl ,0)
     &  +coeff(4,d2x7x56sl )*scints(4,d2x7x56sl ,0)
     &  +coeff(4,d17x2x34sl)*scints(4,d17x2x34sl,0)
     &  +coeff(4,d17x2x56sl)*scints(4,d17x2x56sl,0)
     &  +coeff(4,d7x34x12sl)*scints(4,d7x34x12sl,0)
     &  +coeff(4,d7x12x34sl)*scints(4,d7x12x34sl,0)
     &  +coeff(4,d12x7x34sl)*scints(4,d12x7x34sl,0)
        A7a_slc=A7a_slc
     &  +coeff(3,c17x56sl)*scints(3,c17x56sl,0)
     &  +coeff(3,c17x34sl)*scints(3,c17x34sl,0)
     &  +coeff(3,c27x56sl)*scints(3,c27x56sl,0)
     &  +coeff(3,c34x27sl)*scints(3,c34x27sl,0)
     &  +coeff(3,c34x56sl)*scints(3,c34x56sl,0)
     &  +coeff(3,c12x34sl)*scints(3,c12x34sl,0)
     &  +coeff(3,c12x56sl)*scints(3,c12x56sl,0)
     &  +coeff(2,b127sl)*scints(2,b127sl,0)
     &  +coeff(2,b156sl)*scints(2,b156sl,0)
     &  +coeff(2,b134sl)*scints(2,b134sl,0)
     &  +coeff(2,b56sl) *scints(2,b56sl,0)
     &  +coeff(2,b34sl) *scints(2,b34sl,0)
     &  +coeff(2,b234sl)*scints(2,b234sl,0)
     &  +coeff(2,b256sl)*scints(2,b256sl,0)
     &  +coeff(2,b567sl)*scints(2,b567sl,0)
     &  +coeff(2,b347sl)*scints(2,b347sl,0)
     &  +coeff(2,b17sl) *scints(2,b17sl,0)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
      A7a_slc=
     &  +coeff(4,d1x27x34sl)*scints(4,d27x1x34sl,0)
     &  +coeff(4,d27x1x34sl)*scints(4,d1x27x34sl,0)
     &  +coeff(4,d1x2x34sl )*scints(4,d1x2x56sl ,0)
     &  +coeff(4,d1x2x56sl )*scints(4,d1x2x34sl ,0)
     &  +coeff(4,d1x56x7sl )*scints(4,d1x34x7sl ,0)
     &  +coeff(4,d1x34x7sl )*scints(4,d1x56x7sl ,0)
     &  +coeff(4,d1x2x7sl  )*scints(4,d1x2x7sl  ,0)
     &  +coeff(4,d2x1x7sl  )*scints(4,d2x1x7sl  ,0)
     &  +coeff(4,d1x7x56sl )*scints(4,d1x7x34sl ,0)
     &  +coeff(4,d1x7x34sl )*scints(4,d1x7x56sl ,0)
     &  +coeff(4,d2x1x56sl )*scints(4,d2x1x34sl ,0)
     &  +coeff(4,d2x1x34sl )*scints(4,d2x1x56sl ,0)
     &  +coeff(4,d2x34x7sl )*scints(4,d2x56x7sl ,0)
     &  +coeff(4,d2x56x7sl )*scints(4,d2x34x7sl ,0)
     &  +coeff(4,d2x7x34sl )*scints(4,d2x7x56sl ,0)
     &  +coeff(4,d2x7x56sl )*scints(4,d2x7x34sl ,0)
     &  +coeff(4,d17x2x34sl)*scints(4,d17x2x56sl,0)
     &  +coeff(4,d17x2x56sl)*scints(4,d17x2x34sl,0)
     &  +coeff(4,d7x34x12sl)*scints(4,d7x34x12sl,0)
     &  +coeff(4,d7x12x34sl)*scints(4,d7x12x56sl,0)
     &  +coeff(4,d12x7x34sl)*scints(4,d7x12x34sl,0)
      A7a_slc=A7a_slc
     &  +coeff(3,c17x56sl)*scints(3,c17x34sl,0)
     &  +coeff(3,c17x34sl)*scints(3,c17x56sl,0)
     &  +coeff(3,c27x56sl)*scints(3,c34x27sl,0)
     &  +coeff(3,c34x27sl)*scints(3,c27x56sl,0)
     &  +coeff(3,c34x56sl)*scints(3,c56x34sl,0)
     &  +coeff(3,c12x34sl)*scints(3,c12x56sl,0)
     &  +coeff(3,c12x56sl)*scints(3,c12x34sl,0)
      A7a_slc=A7a_slc
     &  +coeff(2,b127sl)*scints(2,b127sl,0)
     &  +coeff(2,b156sl)*scints(2,b134sl,0)
     &  +coeff(2,b134sl)*scints(2,b156sl,0)
     &  +coeff(2,b56sl) *scints(2,b34sl,0)
     &  +coeff(2,b34sl) *scints(2,b56sl,0)
     &  +coeff(2,b234sl)*scints(2,b256sl,0)
     &  +coeff(2,b256sl)*scints(2,b234sl,0)
     &  +coeff(2,b567sl)*scints(2,b347sl,0)
     &  +coeff(2,b347sl)*scints(2,b567sl,0)
     &  +coeff(2,b17sl) *scints(2,b17sl,0)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        A7a_slc=
     &  +coeff(4,d1x27x34sl)*scints(4,d17x2x56sl,0)
     &  +coeff(4,d27x1x34sl)*scints(4,d17x2x34sl,0)
     &  +coeff(4,d1x2x34sl )*scints(4,d2x1x34sl ,0)
     &  +coeff(4,d1x2x56sl )*scints(4,d2x1x56sl ,0)
     &  +coeff(4,d1x56x7sl )*scints(4,d2x56x7sl ,0)
     &  +coeff(4,d1x34x7sl )*scints(4,d2x34x7sl ,0)
     &  +coeff(4,d1x2x7sl  )*scints(4,d2x1x7sl  ,0)
     &  +coeff(4,d2x1x7sl  )*scints(4,d1x2x7sl  ,0)
     &  +coeff(4,d1x7x56sl )*scints(4,d2x7x56sl ,0)
     &  +coeff(4,d1x7x34sl )*scints(4,d2x7x34sl ,0)
     &  +coeff(4,d2x1x56sl )*scints(4,d1x2x56sl ,0)
     &  +coeff(4,d2x1x34sl )*scints(4,d1x2x34sl ,0)
     &  +coeff(4,d2x34x7sl )*scints(4,d1x34x7sl ,0)
     &  +coeff(4,d2x56x7sl )*scints(4,d1x56x7sl ,0)
     &  +coeff(4,d2x7x34sl )*scints(4,d1x7x34sl ,0)
     &  +coeff(4,d2x7x56sl )*scints(4,d1x7x56sl ,0)
     &  +coeff(4,d17x2x34sl)*scints(4,d27x1x34sl,0)
     &  +coeff(4,d17x2x56sl)*scints(4,d1x27x34sl,0)
     &  +coeff(4,d7x34x12sl)*scints(4,d7x34x12sl,0)
     &  +coeff(4,d7x12x34sl)*scints(4,d7x12x34sl,0)
     &  +coeff(4,d12x7x34sl)*scints(4,d12x7x34sl,0)
      A7a_slc=A7a_slc
     &  +coeff(3,c17x56sl)*scints(3,c27x56sl,0)
     &  +coeff(3,c17x34sl)*scints(3,c34x27sl,0)
     &  +coeff(3,c27x56sl)*scints(3,c17x56sl,0)
     &  +coeff(3,c34x27sl)*scints(3,c17x34sl,0)
     &  +coeff(3,c34x56sl)*scints(3,c34x56sl,0)
     &  +coeff(3,c12x34sl)*scints(3,c12x34sl,0)
     &  +coeff(3,c12x56sl)*scints(3,c12x56sl,0)
      A7a_slc=A7a_slc
     &  +coeff(2,b127sl)*scints(2,b127sl,0)
     &  +coeff(2,b156sl)*scints(2,b256sl,0)
     &  +coeff(2,b134sl)*scints(2,b234sl,0)
     &  +coeff(2,b56sl) *scints(2,b56sl,0)
     &  +coeff(2,b34sl) *scints(2,b34sl,0)
     &  +coeff(2,b234sl)*scints(2,b134sl,0)
     &  +coeff(2,b256sl)*scints(2,b156sl,0)
     &  +coeff(2,b567sl)*scints(2,b567sl,0)
     &  +coeff(2,b347sl)*scints(2,b347sl,0)
     &  +coeff(2,b17sl) *scints(2,b27sl,0)
     &  +coeff(0,iratsl)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        A7a_slc=
     &  +coeff(4,d1x27x34sl)*scints(4,d17x2x34sl,0)
     &  +coeff(4,d27x1x34sl)*scints(4,d17x2x56sl,0)
     &  +coeff(4,d1x2x34sl )*scints(4,d2x1x56sl ,0)
     &  +coeff(4,d1x2x56sl )*scints(4,d2x1x34sl ,0)
     &  +coeff(4,d1x56x7sl )*scints(4,d2x34x7sl ,0)
     &  +coeff(4,d1x34x7sl )*scints(4,d2x56x7sl ,0)
     &  +coeff(4,d1x2x7sl  )*scints(4,d2x1x7sl  ,0)
     &  +coeff(4,d2x1x7sl  )*scints(4,d1x2x7sl  ,0)
     &  +coeff(4,d1x7x56sl )*scints(4,d2x7x34sl ,0)
     &  +coeff(4,d1x7x34sl )*scints(4,d2x7x56sl ,0)
     &  +coeff(4,d2x1x56sl )*scints(4,d1x2x34sl ,0)
     &  +coeff(4,d2x1x34sl )*scints(4,d1x2x56sl ,0)
     &  +coeff(4,d2x34x7sl )*scints(4,d1x56x7sl ,0)
     &  +coeff(4,d2x56x7sl )*scints(4,d1x34x7sl ,0)
     &  +coeff(4,d2x7x34sl )*scints(4,d1x7x56sl ,0)
     &  +coeff(4,d2x7x56sl )*scints(4,d1x7x34sl ,0)
     &  +coeff(4,d17x2x34sl)*scints(4,d1x27x34sl,0)
     &  +coeff(4,d17x2x56sl)*scints(4,d27x1x34sl,0)
     &  +coeff(4,d7x34x12sl)*scints(4,d7x34x12sl,0)
     &  +coeff(4,d7x12x34sl)*scints(4,d7x12x56sl,0)
     &  +coeff(4,d12x7x34sl)*scints(4,d7x12x34sl,0)
        A7a_slc=A7a_slc
     &  +coeff(3,c17x56sl)*scints(3,c34x27sl,0)
     &  +coeff(3,c17x34sl)*scints(3,c27x56sl,0)
     &  +coeff(3,c27x56sl)*scints(3,c17x34sl,0)
     &  +coeff(3,c34x27sl)*scints(3,c17x56sl,0)
     &  +coeff(3,c34x56sl)*scints(3,c56x34sl,0)
     &  +coeff(3,c12x34sl)*scints(3,c12x56sl,0)
     &  +coeff(3,c12x56sl)*scints(3,c12x34sl,0)
        A7a_slc=A7a_slc
     &  +coeff(2,b127sl)*scints(2,b127sl,0)
     &  +coeff(2,b156sl)*scints(2,b234sl,0)
     &  +coeff(2,b134sl)*scints(2,b256sl,0)
     &  +coeff(2,b56sl) *scints(2,b34sl,0)
     &  +coeff(2,b34sl) *scints(2,b56sl,0)
     &  +coeff(2,b234sl)*scints(2,b156sl,0)
     &  +coeff(2,b256sl)*scints(2,b134sl,0)
     &  +coeff(2,b567sl)*scints(2,b347sl,0)
     &  +coeff(2,b347sl)*scints(2,b567sl,0)
     &  +coeff(2,b17sl) *scints(2,b27sl,0)
     &  +coeff(0,iratsl)
      endif

      Alo=Alo+Aloswap
c Contribution from known pole terms
      A7a_slc=A7a_slc+Alo*(
     & +epinv*epinv2*(-1._dp)
     & +epinv*(-1.5_dp-lnrat(musq,-s(j1,j2)))
     & -0.5_dp*lnrat(musq,-s(j1,j2))**2)

      if (swapVV .eqv. .true.) then
        Alo=-Alo
        A7a_slc=-A7a_slc
      endif

      if (verbose) then
        write(6,*)
        write(6,*) 'total',A7a_slc/tot
      endif

      return
      end
