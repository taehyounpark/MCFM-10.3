!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjeta_rat(j1,j2,j3,j4,j5,j6,j7,za,zb,coeff,Alo)
c---- Fills appropriate entries in coeff with the rational parts
c---- present at leading and subleading color

c--- MCFM notation
c---   u(1) ubar(2) nu(3) e+(4) e-(5) nubar(6) g(7)
c--- correspoding to Tania notation
c---   u(6) ubar(1) nu(3) e+(2) e-(5) nubar(4) g(7)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp):: Alo,rat_lc_save,rat_slc_save


      call AWWjeta_rat(j1,j2,j3,j4,j5,j6,j7,za,zb,coeff,Alo)
      rat_lc_save=coeff(0,irat)
      rat_slc_save=coeff(0,iratsl)
      call AWWjeta_rat(j1,j2,j6,j5,j4,j3,j7,za,zb,coeff,Alo)

      coeff(0,irat)  =coeff(0,irat)+rat_lc_save
      coeff(0,iratsl)=coeff(0,iratsl)+rat_slc_save

      return
      end

