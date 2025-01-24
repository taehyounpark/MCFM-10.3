!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c Note: this routine is now in DKS notation, i.e. W-(l3 a4) W+(a5 l6)
      subroutine AWWjetsr_amps(helname,swap12,j1,j2,j3,j4,j5,j6,j7,
     & Qid,hq,za,zb,swapz,Alo,A7b_lc,A7b_slc)
c--- compute amplitude A^b_7 virtual, leading and subleading color pieces,
c--- for positive gluon helicity
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j5,j6,j4,j3,j7,Qid,hq
      complex(dp):: Alo,A7b_slc,A7b_lc
      logical:: swap12,swapz
      character(len=2):: helname
c---- positive-helicity gluon (standard case)
      call AZWWsr2(Alo,A7b_slc,j2,j1,j7,j3,j4,j6,j5,Qid,hq,za,zb)
      call AZWWsr1(Alo,A7b_lc,j2,j7,j1,j3,j4,j6,j5,Qid,hq,za,zb)
      Alo=-Alo
      A7b_lc=-A7b_lc
      A7b_slc=-A7b_slc

c Extra minus sign to get amplitude for RH quark from LH one
      if (swap12) then
        Alo=-Alo
        A7b_lc=-A7b_lc
        A7b_slc=-A7b_slc
      endif

      return
      end

