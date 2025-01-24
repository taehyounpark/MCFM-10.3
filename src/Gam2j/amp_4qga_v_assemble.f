!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- CW May 16
c----- this routine fills the various virtual amplitudes for the following process

c-----q(i1)+qb(i2)+Q(i3)+QB(i4)+gamma(i5)

c-----the array is 3 dim and corresponds to the helicities of q(i1), Q(i3) and gamma(i5)

c---- THE PHOTON COUPLES TO THE  q(i1),q(ib) line, no charges are added.

c---  routine returns 2 arrays:  lc prop to del_1  = delta(i1,i3)*delta(i2,i4)
c---                          :  slc prop to del_2 = delta(i1,i2)*delta(i3,i4)

c Removed unnecessary calls, for speed, JC, 5/20/21
c (original routine _original below)

      subroutine amp_qqbQQbga_v_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del1,amp_del2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'blha.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_del1(2,2,2),amp_del2(2,2,2)
      complex(dp) :: amp_del1_lc(2,2,2),amp_del1_slc(2,2,2)
     &     ,amp_del1_nf(2,2,2)
      complex(dp) :: amp_del2_lc(2,2,2),amp_del2_slc(2,2,2)
     &     ,amp_del2_nf(2,2,2)

      complex(dp):: amp_qqbQQbga_del1_lc,amp_qqbQQbga_del1_slc
      complex(dp):: amp_qqbQQbga_del2_lc,amp_qqbQQbga_del2_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del1_lc
     &     ,amp_qqbQQbga_mhvalt_del1_slc
      complex(dp):: amp_qqbQQbga_del2_nf
      complex(dp):: amp_qqbQQbga_mhvalt_del1_nf
     &     ,amp_qqbQQbga_del1_nf
      integer h1,h2,h3


c------lc del1 pieces

      amp_del1_lc(1,1,2) = amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_lc(1,2,2) = amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_lc(2,2,1) = -amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_lc(2,1,1) = -amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_lc(1,1,1) = amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,2,2) = -amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_lc(1,2,1) = amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,1,2) = -amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,za,zb)


c------slc del1 pieces

      amp_del1_slc(1,1,2) = amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_slc(1,2,2) = amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_slc(2,2,1) = -amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_slc(2,1,1) = -amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_slc(1,1,1) = amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,2,2) = -amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_slc(1,2,1) = amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,1,2) = -amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,za,zb)


c------nf del1 pieces

      amp_del1_nf(1,1,2) = amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,za,zb)
      amp_del1_nf(1,2,2) = amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del1_nf(2,2,1) = -amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,zb,za)
      amp_del1_nf(2,1,1) = -amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del1_nf(1,1,1) = amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,2,2) = -amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,za,zb)

      amp_del1_nf(1,2,1) = amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,1,2) = -amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,za,zb)

      if ((useblha == 0).or.(blhatype == 2)) then

c These only contribute if h1=h2, so other contributions are correct but commented out here
c------lc del2 pieces

      amp_del2_lc(1,1,2) = amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,za,zb)
c      amp_del2_lc(1,2,2) = amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_lc(2,2,1) = -amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,zb,za)
c      amp_del2_lc(2,1,1) = -amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_lc(1,1,1) = amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_lc(2,2,2) = -amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,za,zb)

c      amp_del2_lc(1,2,1) = amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,zb,za)
c      amp_del2_lc(2,1,2) = -amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,za,zb)


c------slc del2 pieces

      amp_del2_slc(1,1,2) = amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,za,zb)
c      amp_del2_slc(1,2,2) = amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_slc(2,2,1) = -amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,zb,za)
c      amp_del2_slc(2,1,1) = -amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_slc(1,1,1) = amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_slc(2,2,2) = -amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,za,zb)

c      amp_del2_slc(1,2,1) = amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,zb,za)
c      amp_del2_slc(2,1,2) = -amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,za,zb)


c------nf del2 pieces

      amp_del2_nf(1,1,2) = amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,za,zb)
c      amp_del2_nf(1,2,2) = amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del2_nf(2,2,1) = -amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,zb,za)
c      amp_del2_nf(2,1,1) = -amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del2_nf(1,1,1) = amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del2_nf(2,2,2) = -amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,za,zb)

c      amp_del2_nf(1,2,1) = amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,zb,za)
c      amp_del2_nf(2,1,2) = -amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,za,zb)

      endif


c=======final assemble
      amp_del1(:,:,:)=czip
      amp_del2(:,:,:)=czip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               amp_del1(h1,h2,h3) = xn*amp_del1_lc(h1,h2,h3)
     &      +1._dp/xn*amp_del1_slc(h1,h2,h3)+nf/xn*amp_del1_nf(h1,h2,h3)
               amp_del2(h1,h2,h3) = +amp_del2_lc(h1,h2,h3)
     &              +1._dp/xn**2*amp_del2_slc(h1,h2,h3)
     &              +nf/xn**2*amp_del2_nf(h1,h2,h3)
            enddo
         enddo
      enddo

      return
      end


      subroutine amp_qqbQQbga_v_assemble_original(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del1,amp_del2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_del1(2,2,2),amp_del2(2,2,2)
      complex(dp) :: amp_del1_lc(2,2,2),amp_del1_slc(2,2,2)
     &     ,amp_del1_nf(2,2,2)
      complex(dp) :: amp_del2_lc(2,2,2),amp_del2_slc(2,2,2)
     &     ,amp_del2_nf(2,2,2)

      complex(dp):: amp_qqbQQbga_del1_lc,amp_qqbQQbga_del1_slc
      complex(dp):: amp_qqbQQbga_del2_lc,amp_qqbQQbga_del2_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del1_lc
     &     ,amp_qqbQQbga_mhvalt_del1_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del2_lc
     &     ,amp_qqbQQbga_mhvalt_del2_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del2_nf
     &     ,amp_qqbQQbga_del2_nf
      complex(dp):: amp_qqbQQbga_mhvalt_del1_nf
     &     ,amp_qqbQQbga_del1_nf
      integer h1,h2,h3


c------lc del1 pieces

      amp_del1_lc(1,1,2) = amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_lc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_lc(2,2,1) = -amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_lc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_lc(1,1,1) = amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,2,2) = -amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_lc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,za,zb)


c------slc del1 pieces

      amp_del1_slc(1,1,2) = amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_slc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_slc(2,2,1) = -amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_slc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_slc(1,1,1) = amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,2,2) = -amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_slc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,za,zb)


c------nf del1 pieces

      amp_del1_nf(1,1,2) = amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,za,zb)
      amp_del1_nf(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del1_nf(2,2,1) = -amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,zb,za)
      amp_del1_nf(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del1_nf(1,1,1) = amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,2,2) = -amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,za,zb)

      amp_del1_nf(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,za,zb)



c------lc del2 pieces

      amp_del2_lc(1,1,2) = amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,za,zb)
      amp_del2_lc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_lc(2,2,1) = -amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,zb,za)
      amp_del2_lc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_lc(1,1,1) = amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_lc(2,2,2) = -amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,za,zb)

      amp_del2_lc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_lc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,za,zb)


c------slc del2 pieces

      amp_del2_slc(1,1,2) = amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,za,zb)
      amp_del2_slc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_slc(2,2,1) = -amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,zb,za)
      amp_del2_slc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_slc(1,1,1) = amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_slc(2,2,2) = -amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,za,zb)

      amp_del2_slc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_slc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,za,zb)


c------nf del2 pieces

      amp_del2_nf(1,1,2) = amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,za,zb)
      amp_del2_nf(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del2_nf(2,2,1) = -amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,zb,za)
      amp_del2_nf(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del2_nf(1,1,1) = amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del2_nf(2,2,2) = -amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,za,zb)

      amp_del2_nf(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del2_nf(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,za,zb)


c=======final assemble
      amp_del1(:,:,:)=czip
      amp_del2(:,:,:)=czip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               amp_del1(h1,h2,h3) = xn*amp_del1_lc(h1,h2,h3)
     &      +1._dp/xn*amp_del1_slc(h1,h2,h3)+nf/xn*amp_del1_nf(h1,h2,h3)
               amp_del2(h1,h2,h3) = +amp_del2_lc(h1,h2,h3)
     &              +1._dp/xn**2*amp_del2_slc(h1,h2,h3)
     &              +nf/xn**2*amp_del2_nf(h1,h2,h3)
            enddo
         enddo
      enddo

c------test with the Kirill Check
      if(1 == 2) then
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* xn lc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_lc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* 1/xn slc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_slc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* nf  pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_nf(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)


                  write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* 1 lc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_lc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* 1/xn**2 slc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_slc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* nf/xn  pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_nf(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         stop
      endif

      return
      end
