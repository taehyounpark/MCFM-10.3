!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine ga_qqbQQb_v(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del1_a,amp_del2_a,amp_del1_b,amp_del2_b)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_del1_a(2,2,2),amp_del2_a(2,2,2)
      complex(dp) :: amp_del1_b(2,2,2),amp_del2_b(2,2,2)
      complex(dp) :: phase,temp(2)
      integer h1,h2,h3

c----- option for extracted phase,
      phase=-cone

      call amp_qqbQQbga_v_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del1_a,amp_del2_a)
      call amp_qqbQQbga_v_assemble(i3,i4,i1,i2,i5,za,zb
     &     ,amp_del1_b,amp_del2_b)

      temp(:)=amp_del1_b(1,2,:)
      amp_del1_b(1,2,:)=amp_del1_b(2,1,:)
      amp_del1_b(2,1,:)=temp
      temp(:)=amp_del2_b(1,2,:)
      amp_del2_b(1,2,:)=amp_del2_b(2,1,:)
      amp_del2_b(2,1,:)=temp


      amp_del1_a(:,:,:)=phase*amp_del1_a(:,:,:)
      amp_del2_a(:,:,:)=phase*amp_del2_a(:,:,:)
      amp_del1_b(:,:,:)=phase*amp_del1_b(:,:,:)
      amp_del2_b(:,:,:)=phase*amp_del2_b(:,:,:)

      if(1 == 2) then
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* ALL DEL 1 QA pieces                   *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_a(h1,h2,h3)
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* ALL DEL 1 QB pieces                   *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_b(h1,h2,h3)
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
c         pause
      endif

      return
      end




