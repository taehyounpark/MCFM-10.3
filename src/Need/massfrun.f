!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c-----------------------------------------------------------------------------

      function massfrun(mf,scale,asmz,nloop)
      implicit none
      include 'types.f'
      real(dp):: massfrun

c-----------------------------------------------------------------------------

c       This function returns the 'nloop' value of a MSbar fermion mass
c       at a given scale.

c       INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale
c              scale = scale at which the running mass is evaluated
c              asmz  = AS(MZ) : this is passed to alphas(scale,asmz,2)
c              nloop = # of loops in the evolutionC

c       COMMON BLOCKS: COMMON/QMASS/CMASS,BMASS,TMASS
c                      contains the MS-bar masses of the heavy quarks.

c       EXTERNAL:      real(dp):: alphas(scale,asmz,2)

c-----------------------------------------------------------------------------


      include 'masses.f'

c     ARGUMENTS

      real(dp)::  mf,scale,asmz
      integer::           nloop

c     LOCAL

      real(dp)::  beta0, beta1,gamma0,gamma1
      real(dp)::  A1,as,asmf,l2
      integer::  nf

c     EXTERNAL

      real(dp)::  alphas
      external          alphas

c     COMMON

c      real *8      cmass,bmass,tmass
c      COMMON/QMASS/CMASS,BMASS,TMASS

c     CONSTANTS

      real(dp)::  One, Two, Three, Pi
      parameter( One = 1._dp, Two = 2._dp, Three = 3._dp )
      parameter( Pi = 3.14159265358979323846_dp)
cc



      if ( mf>mt ) then
         nf = 6
      else
         nf = 5
      end if

      beta0 = ( 11._dp - Two/Three *nf )/4._dp
      beta1 = ( 102._dp  - 38._dp/Three*nf )/16._dp
      gamma0= one
      gamma1= ( 202._dp/3._dp  - 20._dp/9._dp*nf )/16._dp
      A1    = -beta1*gamma0/beta0**2+gamma1/beta0
      as    = alphas(scale,asmz,nloop)
      asmf  = alphas(mf   ,asmz,nloop)
      l2    = (one+A1*as/Pi)/(one+A1*asmf/Pi)

      massfrun = mf * (as/asmf)**(gamma0/beta0)

      if(nloop==2) massfrun=massfrun*l2
ccc
      return
      end


