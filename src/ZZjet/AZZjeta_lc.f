!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjeta_lc(helname,swap12,swapVV,p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz,scints,coeff,Alo,A7a_lc)
c--- MCFM notation
c---   u(1) ubar(2) e-(3) e+(4) mu-(5) mu+(6) g(7)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'WWjetlabels.f'
      include 'KCdef.f'
      include 'verbose.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j6,j5,j7
      complex(dp):: Alo,Aloswap,A7treea,A7a_lc,
     & coeff3456(2),coeff6543(2),tri,tri_sl,lnrat
      complex(dp):: xd7x1x56,xd1x27x34,xd7x2x34,xd2x17x56,xd1x7x2,
     & bub17,bub156,bub234,bub127,bub56,b34extra,b56extra
      real(dp):: p(mxpart,4)
      character(len=4):: helname
      logical:: swap12,swapVV,swapz

c--- call standard MCFM routine to get LO
      Alo=A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)
      Aloswap=A7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)

c      if (verbose) call parser('runZZ','LCa',helname,KCD,KCC,KCB,rattot,tot)
c      include 'runZZLCa.f'
c      include 'runZZLCa.txt'

      coeff(4,d7x1x56)=xd7x1x56(j1,j2,j3,j4,j6,j5,j7,za,zb)
      coeff(4,d7x1x34)=xd7x1x56(j1,j2,j6,j5,j3,j4,j7,za,zb)

      coeff(4,d27x1x34)=xd1x27x34(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d1x27x34)=xd1x27x34(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d7x2x56)=xd7x2x34(j1,j2,j6,j5,j3,j4,j7,za,zb)
      coeff(4,d7x2x34)=xd7x2x34(j1,j2,j3,j4,j6,j5,j7,za,zb)

      coeff(4,d2x17x56)=xd2x17x56(j1,j2,j3,j4,j6,j5,j7,za,zb)
      coeff(4,d17x2x56)=xd2x17x56(j1,j2,j6,j5,j3,j4,j7,za,zb)

      coeff(4,d1x7x2)=(xd1x7x2(j1,j2,j3,j4,j6,j5,j7,za,zb)
     &                +xd1x7x2(j1,j2,j6,j5,j3,j4,j7,za,zb))

      if (verbose) then
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'d7x1x56 ',coeff(4,d7x1x56)/KCD(d7x1x56)
        write(6,*)'d7x1x34 ',coeff(4,d7x1x34)/KCD(d7x1x34)
        write(6,*)'d27x1x34',coeff(4,d27x1x34)/KCD(d27x1x34)
        write(6,*)'d1x27x34',coeff(4,d1x27x34)/KCD(d1x27x34)
        write(6,*)'d7x2x56 ',coeff(4,d7x2x56)/KCD(d7x2x56)
        write(6,*)'d7x2x34 ',coeff(4,d7x2x34)/KCD(d7x2x34)
        write(6,*)'d2x17x56',coeff(4,d2x17x56)/KCD(d2x17x56)
        write(6,*)'d17x2x56',coeff(4,d17x2x56)/KCD(d17x2x56)
        write(6,*)'d1x7x2  ',coeff(4,d1x7x2)/KCD(d1x7x2)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'d7x1x56 ',coeff(4,d7x1x56)/KCD(d7x1x34)
        write(6,*)'d7x1x34 ',coeff(4,d7x1x34)/KCD(d7x1x56)
        write(6,*)'d27x1x34',coeff(4,d27x1x34)/KCD(d1x27x34)
        write(6,*)'d1x27x34',coeff(4,d1x27x34)/KCD(d27x1x34)
        write(6,*)'d7x2x56 ',coeff(4,d7x2x56)/KCD(d7x2x34)
        write(6,*)'d7x2x34 ',coeff(4,d7x2x34)/KCD(d7x2x56)
        write(6,*)'d2x17x56',coeff(4,d2x17x56)/KCD(d17x2x56)
        write(6,*)'d17x2x56',coeff(4,d17x2x56)/KCD(d2x17x56)
        write(6,*)'d1x7x2  ',coeff(4,d1x7x2)/KCD(d1x7x2)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'d7x1x56 ',coeff(4,d7x1x56)/KCD(d7x2x56)
        write(6,*)'d7x1x34 ',coeff(4,d7x1x34)/KCD(d7x2x34)
        write(6,*)'d27x1x34',coeff(4,d27x1x34)/KCD(d2x17x56)
        write(6,*)'d1x27x34',coeff(4,d1x27x34)/KCD(d17x2x56)
        write(6,*)'d7x2x56 ',coeff(4,d7x2x56)/KCD(d7x1x56)
        write(6,*)'d7x2x34 ',coeff(4,d7x2x34)/KCD(d7x1x34)
        write(6,*)'d2x17x56',coeff(4,d2x17x56)/KCD(d27x1x34)
        write(6,*)'d17x2x56',coeff(4,d17x2x56)/KCD(d1x27x34)
        write(6,*)'d1x7x2  ',coeff(4,d1x7x2)/KCD(d1x7x2)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'d7x1x56 ',coeff(4,d7x1x56)/KCD(d7x2x34)
        write(6,*)'d7x1x34 ',coeff(4,d7x1x34)/KCD(d7x2x56)
        write(6,*)'d27x1x34',coeff(4,d27x1x34)/KCD(d17x2x56)
        write(6,*)'d1x27x34',coeff(4,d1x27x34)/KCD(d2x17x56)
        write(6,*)'d7x2x56 ',coeff(4,d7x2x56)/KCD(d7x1x34)
        write(6,*)'d7x2x34 ',coeff(4,d7x2x34)/KCD(d7x1x56)
        write(6,*)'d2x17x56',coeff(4,d2x17x56)/KCD(d1x27x34)
        write(6,*)'d17x2x56',coeff(4,d17x2x56)/KCD(d27x1x34)
        write(6,*)'d1x7x2  ',coeff(4,d1x7x2)/KCD(d1x7x2)
      endif
      endif

c--- Triangles (set up subleading color here too)
      call tri56x34(j1,j2,j3,j4,j6,j5,j7,za,zb,tri,tri_sl)
      coeff(3,c56x34)=tri
      coeff(3,c56x34sl)=tri_sl
      call tri56x34(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c56x34)=coeff(3,c56x34)+tri
      coeff(3,c56x34sl)=coeff(3,c56x34sl)+tri_sl

      call tri56x17(j1,j2,j3,j4,j6,j5,j7,za,zb,tri,tri_sl)
      coeff(3,c56x17)=tri
      coeff(3,c56x17sl)=tri_sl

      call tri56x17(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c17x34)=tri
      coeff(3,c17x34sl)=tri_sl

      call tri34x27(j1,j2,j3,j4,j6,j5,j7,za,zb,tri,tri_sl)
      coeff(3,c34x27)=tri
      coeff(3,c34x27sl)=tri_sl

      call tri34x27(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c27x56)=tri
      coeff(3,c27x56sl)=tri_sl

      call tri12x34sl(j1,j2,j3,j4,j6,j5,j7,za,zb,tri,tri_sl)
      coeff(3,c12x34sl)=tri_sl
      call tri12x56sl(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c12x34sl)=coeff(3,c12x34sl)+tri_sl

      call tri12x56sl(j1,j2,j3,j4,j6,j5,j7,za,zb,tri,tri_sl)
      coeff(3,c12x56sl)=tri_sl
      call tri12x34sl(j1,j2,j6,j5,j3,j4,j7,za,zb,tri,tri_sl)
      coeff(3,c12x56sl)=coeff(3,c12x56sl)+tri_sl

      if (verbose) then
      write(6,*)
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'c56x17',KCC(c56x17)/coeff(3,c56x17)
        write(6,*)'c17x34',KCC(c17x34)/coeff(3,c17x34)
        write(6,*)'c27x56',KCC(c27x56)/coeff(3,c27x56)
        write(6,*)'c34x27',KCC(c34x27)/coeff(3,c34x27)
        write(6,*)'c56x34',KCC(c56x34)/coeff(3,c56x34)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'c56x17',KCC(c56x17)/coeff(3,c17x34)
        write(6,*)'c17x34',KCC(c17x34)/coeff(3,c56x17)
        write(6,*)'c27x56',KCC(c27x56)/coeff(3,c34x27)
        write(6,*)'c34x27',KCC(c34x27)/coeff(3,c27x56)
        write(6,*)'c56x34',KCC(c56x34)/coeff(3,c56x34)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        write(6,*)'c56x17',KCC(c27x56)/coeff(3,c56x17)
        write(6,*)'c17x34',KCC(c34x27)/coeff(3,c17x34)
        write(6,*)'c27x56',KCC(c56x17)/coeff(3,c27x56)
        write(6,*)'c34x27',KCC(c17x34)/coeff(3,c34x27)
        write(6,*)'c56x34',KCC(c56x34)/coeff(3,c56x34)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        write(6,*)'c56x17',KCC(c27x56)/coeff(3,c17x34)
        write(6,*)'c17x34',KCC(c34x27)/coeff(3,c56x17)
        write(6,*)'c27x56',KCC(c56x17)/coeff(3,c34x27)
        write(6,*)'c34x27',KCC(c17x34)/coeff(3,c27x56)
        write(6,*)'c56x34',KCC(c56x34)/coeff(3,c56x34)
      endif

c      write(6,*)
c      write(6,*)'c7x134',KCC(c7x134),coeff(3,c7x134)
c      write(6,*)'c1x7',KCC(c1x7),coeff(3,c1x7)
c      write(6,*)'c1x56',KCC(c1x56),coeff(3,c1x56)
c      write(6,*)'c1x34',KCC(c1x34),coeff(3,c1x34)
c      write(6,*)'c2x7',KCC(c2x7),coeff(3,c2x7)
c      write(6,*)'c2x56',KCC(c2x56),coeff(3,c2x56)
c      write(6,*)'c2x34',KCC(c2x34),coeff(3,c2x34)
c      write(6,*)'c2x17',KCC(c2x17),coeff(3,c2x17)
      endif

c--- Leading color bubbles
      coeff(2,b134)=bub156(j1,j2,j6,j5,j4,j3,j7,za,zb)
      coeff(2,b156)=bub156(j1,j2,j3,j4,j5,j6,j7,za,zb)
      coeff(2,b256)=bub234(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff(2,b234)=bub234(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff(2,b34)=bub56(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff(2,b56)=bub56(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff3456(b127)=bub127(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b127)=bub127(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
c      write(6,*) 'coeff3456(b127)+coeff6543(b127)',coeff3456(b127)+coeff6543(b127)
c      tmp=2*bub127_sym(j1,j2,j3,j4,j5,j6,j7,za,zb)
c      write(6,*) '2*bub127_sym(p,j1,j2,j6,j5,j4,j3,j7)',tmp
c      write(6,*) 'coeff3456(b127)+coeff6543(b127)',(coeff3456(b127)+coeff6543(b127))/tmp
c      write(6,*) 'pause in AZZjeta_lc.f'
c      read(5,*)
      coeff3456(b17)=bub17(p,j1,j2,j3,j4,j5,j6,j7,za,zb,swapz)
      coeff6543(b17)=bub17(p,j1,j2,j6,j5,j4,j3,j7,za,zb,swapz)
      coeff(2,b127)=coeff3456(b127)+coeff6543(b127)
      coeff(2,b17)=coeff3456(b17)+coeff6543(b17)

c--- remaining bubble bits from pole structure
      b56extra=-1.5d0*Alo-coeff(2,b134)
     & -coeff(2,b256)-coeff(2,b34)
     & -coeff3456(b17)-coeff3456(b127)
      b34extra=-1.5d0*Aloswap-coeff(2,b156)
     & -coeff(2,b234)-coeff(2,b56)
     & -coeff6543(b17)-coeff6543(b127)
      coeff(2,b56)=coeff(2,b56)+b56extra
      coeff(2,b34)=coeff(2,b34)+b34extra
      if (verbose) then
      write(6,*)
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        write(6,*) 'b156',coeff(2,b156)/KCB(b156)
        write(6,*) 'b134',coeff(2,b134)/KCB(b134)
        write(6,*) 'b234',coeff(2,b234)/KCB(b234)
        write(6,*) 'b256',coeff(2,b256)/KCB(b256)
        write(6,*) 'b56',coeff(2,b56)/KCB(b56)
        write(6,*) 'b34',coeff(2,b34)/KCB(b34)
        write(6,*) 'b127',coeff(2,b127)/KCB(b127)
        write(6,*) 'b17',coeff(2,b17)/KCB(b17)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        write(6,*) 'b156',coeff(2,b156)/KCB(b134)
        write(6,*) 'b134',coeff(2,b134)/KCB(b156)
        write(6,*) 'b234',coeff(2,b234)/KCB(b256)
        write(6,*) 'b256',coeff(2,b256)/KCB(b234)
        write(6,*) 'b56',coeff(2,b56)/KCB(b34)
        write(6,*) 'b34',coeff(2,b34)/KCB(b56)
        write(6,*) 'b127',coeff(2,b127)/KCB(b127)
        write(6,*) 'b17',coeff(2,b17)/KCB(b17)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        write(6,*) 'b156',coeff(2,b156)/KCB(b256)
        write(6,*) 'b134',coeff(2,b134)/KCB(b234)
        write(6,*) 'b234',coeff(2,b234)/KCB(b134)
        write(6,*) 'b256',coeff(2,b256)/KCB(b156)
        write(6,*) 'b56',coeff(2,b56)/KCB(b56)
        write(6,*) 'b34',coeff(2,b34)/KCB(b34)
        write(6,*) 'b127',coeff(2,b127)/KCB(b127)
        write(6,*) 'b17',coeff(2,b17)/KCB(b27)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        write(6,*) 'b156',coeff(2,b156)/KCB(b234)
        write(6,*) 'b134',coeff(2,b134)/KCB(b256)
        write(6,*) 'b234',coeff(2,b234)/KCB(b156)
        write(6,*) 'b256',coeff(2,b256)/KCB(b134)
        write(6,*) 'b56',coeff(2,b56)/KCB(b34)
        write(6,*) 'b34',coeff(2,b34)/KCB(b56)
        write(6,*) 'b127',coeff(2,b127)/KCB(b127)
        write(6,*) 'b17',coeff(2,b17)/KCB(b27)
      endif
      endif

      if (verbose) then
        write(6,*)
        write(6,*) 'rational',coeff(0,irat)/rattot
      endif

c Now assemble total
      if     ((swap12 .eqv. .false.) .and. (swapVV .eqv. .false.)) then
        A7a_lc=
     &  +coeff(4,d7x1x56) *scints(4,d7x1x56,0)
     &  +coeff(4,d7x1x34) *scints(4,d7x1x34,0)
     &  +coeff(4,d27x1x34)*scints(4,d27x1x34,0)
     &  +coeff(4,d1x27x34)*scints(4,d1x27x34,0)
     &  +coeff(4,d7x2x56) *scints(4,d7x2x56,0)
     &  +coeff(4,d7x2x34) *scints(4,d7x2x34,0)
     &  +coeff(4,d2x17x56)*scints(4,d2x17x56,0)
     &  +coeff(4,d17x2x56)*scints(4,d17x2x56,0)
     &  +coeff(4,d1x7x2)  *scints(4,d1x7x2,0)
     &  +coeff(3,c56x17)*scints(3,c56x17,0)
     &  +coeff(3,c17x34)*scints(3,c17x34,0)
     &  +coeff(3,c27x56)*scints(3,c27x56,0)
     &  +coeff(3,c34x27)*scints(3,c34x27,0)
     &  +coeff(3,c56x34)*scints(3,c56x34,0)
     &  +coeff(2,b156)*scints(2,b156,0)
     &  +coeff(2,b134)*scints(2,b134,0)
     &  +coeff(2,b234)*scints(2,b234,0)
     &  +coeff(2,b256)*scints(2,b256,0)
     &  +coeff(2,b56) *scints(2,b56,0)
     &  +coeff(2,b34) *scints(2,b34,0)
     &  +coeff(2,b127)*scints(2,b127,0)
     &  +coeff(2,b17) *scints(2,b17,0)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .false.) .and. (swapVV .eqv. .true.)) then
        A7a_lc=
     &  +coeff(4,d7x1x56) *scints(4,d7x1x34,0)
     &  +coeff(4,d7x1x34) *scints(4,d7x1x56,0)
     &  +coeff(4,d27x1x34)*scints(4,d1x27x34,0)
     &  +coeff(4,d1x27x34)*scints(4,d27x1x34,0)
     &  +coeff(4,d7x2x56) *scints(4,d7x2x34,0)
     &  +coeff(4,d7x2x34) *scints(4,d7x2x56,0)
     &  +coeff(4,d2x17x56)*scints(4,d17x2x56,0)
     &  +coeff(4,d17x2x56)*scints(4,d2x17x56,0)
     &  +coeff(4,d1x7x2)  *scints(4,d1x7x2,0)
     &  +coeff(3,c17x34)*scints(3,c56x17,0)
     &  +coeff(3,c56x17)*scints(3,c17x34,0)
     &  +coeff(3,c34x27)*scints(3,c27x56,0)
     &  +coeff(3,c27x56)*scints(3,c34x27,0)
     &  +coeff(3,c56x34)*scints(3,c56x34,0)
     &  +coeff(2,b156)*scints(2,b134,0)
     &  +coeff(2,b134)*scints(2,b156,0)
     &  +coeff(2,b234)*scints(2,b256,0)
     &  +coeff(2,b256)*scints(2,b234,0)
     &  +coeff(2,b56) *scints(2,b34,0)
     &  +coeff(2,b34) *scints(2,b56,0)
     &  +coeff(2,b127)*scints(2,b127,0)
     &  +coeff(2,b17) *scints(2,b17,0)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .false.)) then
        A7a_lc=
     &  +coeff(4,d7x1x56) *scints(4,d7x2x56,0)
     &  +coeff(4,d7x1x34) *scints(4,d7x2x34,0)
     &  +coeff(4,d27x1x34)*scints(4,d2x17x56,0)
     &  +coeff(4,d1x27x34)*scints(4,d17x2x56,0)
     &  +coeff(4,d7x2x56) *scints(4,d7x1x56,0)
     &  +coeff(4,d7x2x34) *scints(4,d7x1x34,0)
     &  +coeff(4,d2x17x56)*scints(4,d27x1x34,0)
     &  +coeff(4,d17x2x56)*scints(4,d1x27x34,0)
     &  +coeff(4,d1x7x2)  *scints(4,d1x7x2,0)
     &  +coeff(3,c56x17)*scints(3,c27x56,0)
     &  +coeff(3,c17x34)*scints(3,c34x27,0)
     &  +coeff(3,c27x56)*scints(3,c56x17,0)
     &  +coeff(3,c34x27)*scints(3,c17x34,0)
     &  +coeff(3,c56x34)*scints(3,c56x34,0)
     &  +coeff(2,b156)*scints(2,b256,0)
     &  +coeff(2,b134)*scints(2,b234,0)
     &  +coeff(2,b234)*scints(2,b134,0)
     &  +coeff(2,b256)*scints(2,b156,0)
     &  +coeff(2,b56) *scints(2,b56,0)
     &  +coeff(2,b34) *scints(2,b34,0)
     &  +coeff(2,b127)*scints(2,b127,0)
     &  +coeff(2,b17) *scints(2,b27,0)
     &  +coeff(0,irat)
      elseif ((swap12 .eqv. .true.) .and. (swapVV .eqv. .true.)) then
        A7a_lc=
     &  +coeff(4,d7x1x56) *scints(4,d7x2x34,0)
     &  +coeff(4,d7x1x34) *scints(4,d7x2x56,0)
     &  +coeff(4,d27x1x34)*scints(4,d17x2x56,0)
     &  +coeff(4,d1x27x34)*scints(4,d2x17x56,0)
     &  +coeff(4,d7x2x56) *scints(4,d7x1x34,0)
     &  +coeff(4,d7x2x34) *scints(4,d7x1x56,0)
     &  +coeff(4,d2x17x56)*scints(4,d1x27x34,0)
     &  +coeff(4,d17x2x56)*scints(4,d27x1x34,0)
     &  +coeff(4,d1x7x2)  *scints(4,d1x7x2,0)
     &  +coeff(3,c17x34)*scints(3,c27x56,0)
     &  +coeff(3,c56x17)*scints(3,c34x27,0)
     &  +coeff(3,c34x27)*scints(3,c56x17,0)
     &  +coeff(3,c27x56)*scints(3,c17x34,0)
     &  +coeff(3,c56x34)*scints(3,c56x34,0)
     &  +coeff(2,b156)*scints(2,b234,0)
     &  +coeff(2,b134)*scints(2,b256,0)
     &  +coeff(2,b234)*scints(2,b156,0)
     &  +coeff(2,b256)*scints(2,b134,0)
     &  +coeff(2,b56) *scints(2,b34,0)
     &  +coeff(2,b34) *scints(2,b56,0)
     &  +coeff(2,b127)*scints(2,b127,0)
     &  +coeff(2,b17) *scints(2,b27,0)
     &  +coeff(0,irat)
      endif

      Alo=Alo+Aloswap
c Contribution from known pole terms
      A7a_lc=A7a_lc+Alo*(
     & +epinv*epinv2*(-2._dp)
     & +epinv*(-1.5_dp-(lnrat(musq,-s(j1,j7))+lnrat(musq,-s(j2,j7))))
     & -0.5_dp*(lnrat(musq,-s(j1,j7))**2+lnrat(musq,-s(j2,j7))**2))

      if (swapVV .eqv. .true.) then
        Alo=-Alo
        A7a_lc=-A7a_lc
      endif

      if (verbose) then
        write(6,*)
        write(6,*) 'total',A7a_lc/tot
      endif

      return
      end
