!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c     =--- nf loops -+-+

      function amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_del1_nf
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: ALO,Vpole,lnrat

      ALO=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))

c      write(6,*) ALO
c---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)

      amp_qqbQQbga_del1_nf=Vpole*ALO

c      write(6,*) 'dell 1 -+-+'
c      write(6,*) 'LO',ALO
c      write(6,*) 'total = ',amp_qqbQQbga_del1_nf*im
      return
      end


      function amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_del2_nf
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: ALO,Vpole,lnrat

      ALO=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))

c---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)
      amp_qqbQQbga_del2_nf=-Vpole*ALO

c      write(6,*) 'dell 2 -+-+'
c      write(6,*) 'LO',ALO
c      write(6,*) 'total = ',amp_qqbQQbga_del2_nf*im
c      write(6,*)
      return
      end
c=====nf loops -++-

      function amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del1_nf
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: ALO,Vpole,lnrat

      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))

c---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)

      amp_qqbQQbga_mhvalt_del1_nf=Vpole*ALO

c      write(6,*) 'dell 1 --++'

c      write(6,*) 'LO',ALO
c      write(6,*) 'total = ',amp_qqbQQbga_mhvalt_del1_nf*im
c      write(6,*)
      return
      end


      function amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del2_nf
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp) :: ALO,Vpole,lnrat

      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))

c---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)

      amp_qqbQQbga_mhvalt_del2_nf=-Vpole*ALO
c      write(6,*) 'dell 2 -++-'

c      write(6,*) 'LO',ALO
c      write(6,*) 'total = ',amp_qqbQQbga_mhvalt_del2_nf*im
c      write(6,*)
      return
      end

