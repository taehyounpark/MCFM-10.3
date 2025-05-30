!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function ampsq_2gam2g(i1,i2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      real(dp):: ampsq_2gam2g
c--- Matrix element squared for the process
c---       0  -->  qb(j5) + q(j6) + g(i1) + g(i2) + gam(j3) + gam(j4)
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---

      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,h5,icol
      complex(dp):: amp(2,2,2,2,2,2),
     & amp_2gam2g_mppppm,amp_2gam2g_pmpppm,
     & amp_2gam2g_ppmppm,amp_2gam2g_pppmpm,amp_2gam2g_mmpppm,
     & amp_2gam2g_pmmppm,amp_2gam2g_mpmppm

c--- loop over color ordering of gluons
      do icol=1,2

      if (icol == 1) then
        j1=i1
        j2=i2
      else
        j1=i2
        j2=i1
      endif

c--- ordering of labels in amp is as follows:
c---  (gluon j1, gluon j2, photon j3, photon j4, antiquark j5)
c--- for consistency with function names (and quark j6=3-j5)


c--- trivial amplitude
      amp(icol,2,2,2,2,2)=czip
c--- basic MHV amplitudes
      amp(icol,1,2,2,2,2)=amp_2gam2g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,1,2,2,2)=amp_2gam2g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,2,1,2,2)=amp_2gam2g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,2,2,1,2)=amp_2gam2g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)

c--- non-MHV amplitudes
      amp(icol,1,1,2,2,2)=amp_2gam2g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,1,1,2,2)=amp_2gam2g_pmmppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,1,2,1,2,2)=amp_2gam2g_mpmppm(j1,j2,j3,j4,j5,j6,za,zb)

c--- parity and charge conjugation
c--- (complex conjugate, interchange q, qb and switch gluon order)
      amp(icol,1,1,1,1,2)=czip

      amp(icol,2,1,1,1,2)=-amp_2gam2g_pmpppm(j2,j1,j3,j4,j6,j5,zb,za)
      amp(icol,1,2,1,1,2)=-amp_2gam2g_mppppm(j2,j1,j3,j4,j6,j5,zb,za)
      amp(icol,1,1,2,1,2)=-amp_2gam2g_ppmppm(j2,j1,j3,j4,j6,j5,zb,za)
      amp(icol,1,1,1,2,2)=-amp_2gam2g_pppmpm(j2,j1,j3,j4,j6,j5,zb,za)

      amp(icol,2,2,1,1,2)=-amp_2gam2g_mmpppm(j2,j1,j3,j4,j6,j5,zb,za)
      amp(icol,1,2,2,1,2)=-amp_2gam2g_mpmppm(j2,j1,j3,j4,j6,j5,zb,za)
      amp(icol,2,1,2,1,2)=-amp_2gam2g_pmmppm(j2,j1,j3,j4,j6,j5,zb,za)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(icol,h1,h2,h3,h4,1)=amp(icol,3-h1,3-h2,3-h3,3-h4,2)
      enddo
      enddo
      enddo
      enddo

      enddo

c--- note: obvious redundancy in this routine, but might be
c--- worth checking relations for use in virtual

      ampsq_2gam2g=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
c      write(6,*) h1*10000+h2*1000+h3*100+h4*10+h5,
c     & amp(icol,h1,h2,h3,h4,h5),
c     & abs(amp(1,h1,h2,h3,h4,h5)+amp(2,h2,h1,h3,h4,h5))
      ampsq_2gam2g=ampsq_2gam2g
     & +abs(amp(1,h1,h2,h3,h4,h5))**2
     & +abs(amp(2,h2,h1,h3,h4,h5))**2
     & -1._dp/xn**2*abs(amp(1,h1,h2,h3,h4,h5)+amp(2,h2,h1,h3,h4,h5))**2
      enddo
      enddo
      enddo
      enddo
      enddo

c      write(6,*) 'ampsq_2gam2g',ampsq_2gam2g*(-9._dp)
c      pause

      return
      end

