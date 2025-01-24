!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!


      subroutine qqb_dm_monojet_v_Axamps(p,i1,i2,i3,i4,i5,qgqb)
      implicit none
      include 'types.f'
c----- combined colour and Born intefered amplitudes as a function of quark line helicity
c------ default is q(i1)+g(i2)+qb(i3)+x(i4)+x(i5)
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
c      include 'zprods_decl.f'
      real(dp)::qgqb(2)
      integer::i1,i2,i3,i4,i5
      real(dp)::p(mxpart,4)
      integer::h1,h2,h3,h4
      complex(dp)::amp_tree(2,2,2,2),amp_lc(2,2,2,2),
     & amp_slc(2,2,2,2)
c      real(dp)::s45
c------ tree
      call qqb_dm_monojet_Axamps(p,i1,i2,i3,i4,i5,amp_tree)
c------ leading colour
       call qqb_dm_monojet_lc_Axamps(p,i1,i2,i3,i4,i5,amp_lc)
c------ subleading colour (swap i2 and i3)
      call qqb_dm_monojet_slc_Axamps(p,i1,i3,i2,i4,i5,amp_slc)

c------ interfere with born as function of quark line helicity


c      call spinoru(5,p,za,zb)
c      write(6,*) A51(i1,i2,i3,i4,i5,za,zb)
c      write(6,*) amp_lc(2,2,1,2)/(cone*s45)
c      pause
      qgqb(:)=0d0
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
c                  if(h3+h4==3) then
c        write(6,*) h1,h2,h3,h4,amp_tree(h1,h2,h3,h4),amp_lc(h1,h2,h3,h4)
c                  endif
                  qgqb(h1)=qgqb(h1)
     &                 +ason2pi*Dble(conjg(amp_tree(h1,h2,h3,h4))*
     &                 (amp_lc(h1,h2,h3,h4)
     &                 +amp_slc(h1,h2,h3,h4)/xnsq))
               enddo
            enddo
         enddo
      enddo

c      write(6,*)
c      write(6,*) i1,i2,i3,i4,i5




c      do h1=1,2
c         do h3=1,2
c            do h4=1,2
c               if(h3 /= h4) then
c                  write(6,*) h1,h3,h4,
c     &             '1',amp_tree(h1,1,h3,h4)/(s45*cone)
c     &  ,amp_lc(h1,1,h3,h4)/(s45*cone),
c     &  Dble(conjg(amp_tree(h1,1,h3,h4))*amp_lc(h1,1,h3,h4))/s45**2,
c     &            '2',amp_tree(h1,2,h3,h4)/(s45*cone)
c     &  ,amp_lc(h1,2,h3,h4)/(s45*cone),
c     &  Dble(conjg(amp_tree(h1,2,h3,h4))*amp_lc(h1,1,h3,h4))/s45**2

c     &    ason2pi*Dble(conjg(amp_tree(h1,1,h3,h4))*amp_lc(h1,1,h3,h4))
c     &                 /s45**2+
c     &     ason2pi*Dble(conjg(amp_tree(h1,2,h3,h4))*amp_lc(h1,2,h3,h4))
c     &                 /s45**2
c               endif
c            enddo
c         enddo
c      enddo

      return
      end


