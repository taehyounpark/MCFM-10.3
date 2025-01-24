!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine ZZmassivebub(j1,j2,j3,j4,j5,j6,za,zb,mt,bub,totrat)
      implicit none
      include 'types.f'
c----Bubble coefficients extracted from BDK 11.5, 11.8
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'ggZZintegrals.f'
      include 'blabels.f'
      include 'docheck.f'
c      logical:: ggZZunstable
      integer:: j,j1,j3,j4,j2,j5,j6,h1,h2,h3,h4
      complex(dp):: b(2,2,2,2,7),bub(2,2,2,2,-2:0),totrat(2,2,2,2)
      real(dp):: mt
c      common/ggZZunstable/ggZZunstable

      call ZZmassivebubcore(j1,j2,j3,j4,j5,j6,zb,za,mt,b)

c--- obtain h2=1 results by c.c.
      do h3=1,2
      do h4=1,2
      b(2,1,h3,h4,:)=b(1,2,3-h3,3-h4,:)
      b(1,1,h3,h4,:)=b(2,2,3-h3,3-h4,:)
      enddo
      enddo

      call ZZmassivebubcore(j1,j2,j3,j4,j5,j6,za,zb,mt,b)

c--- compare with obtaining h2=1 coefficients by c.c.
c      do h3=1,2
c      do h4=1,2
c      write(6,*) b(2,1,h3,h4,:),conjg(b(1,2,3-h3,3-h4,:))
c      write(6,*) b(1,1,h3,h4,:),conjg(b(2,2,3-h3,3-h4,:))
c      enddo
c      enddo
c      write(6,*)

c--- copy rational parts to another array to return separately
      totrat(:,:,:,:)=b(:,:,:,:,rat)

      if (docheck) then
      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b34),'( s34;mtsq,mtsq)')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b34),'( s34;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b34),'( s34;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b34),'( s34;mtsq,mtsq)')

      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b134),'(s134;mtsq,mtsq)')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b134),'(s134;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b134),'(s134;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b134),'(s134;mtsq,mtsq)')

      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b56),'( s56;mtsq,mtsq)')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b56),'( s56;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b56),'( s56;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b56),'( s56;mtsq,mtsq)')

      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b234),'(s156;mtsq,mtsq)')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b234),'(s156;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b234),'(s156;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b234),'(s156;mtsq,mtsq)')

      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,b12),'( s12;mtsq,mtsq)')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,b12),'( s12;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,b12),'( s12;mtsq,mtsq)')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,b12),'( s12;mtsq,mtsq)')

      call ggZZparsecheck('LL',-1,-1,b(1,1,:,:,rat),'xe,res/lo  1')
      call ggZZparsecheck('LL',-1,+1,b(1,2,:,:,rat),'xe,res/lo  1')
      call ggZZparsecheck('LL',+1,-1,b(2,1,:,:,rat),'xe,res/lo  1')
      call ggZZparsecheck('LL',+1,+1,b(2,2,:,:,rat),'xe,res/lo  1')
      endif

c--- Check that coefficients sum to zero to test numerical stability
c--- note: only do for (1,2,1,1) and (2,2,1,1) for simplicity
c--- and multiply by s12 to remove dimensions
c      testmm=(b(1,1,1,1,1)+b(1,1,1,1,2)+b(1,1,1,1,3)
c     &       +b(1,1,1,1,4)+b(1,1,1,1,5))*s12
c      testmp=(b(1,2,1,1,1)+b(1,2,1,1,2)+b(1,2,1,1,3)
c     &       +b(1,2,1,1,4)+b(1,2,1,1,5))*s12
c      if ((abs(testmm) > 1d-4) .or. (abs(testmp) > 1d-4)) then
c        write(6,*) 'testmm=',testmm
c        write(6,*) 'testmp=',testmp
c        ggZZunstable=.true.
c      else
c        ggZZunstable=.false.
c      endif

      bub(:,:,:,:,:)=czip

c--- Note: Bint(6) does not appear in this calculation
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,7
      if (j  /=  6) then
      bub(h1,h2,h3,h4,0)=bub(h1,h2,h3,h4,0)+b(h1,h2,h3,h4,j)*B0(j)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo

      if (docheck) then
        do j=1,7
        do h3=1,2
        do h4=1,2
        call ggZZcapture('bubmp'//char(ichar('0')+j),h3,h4,1,2,3,4,5,6,
     &   b(1,2,h3,h4,j),
     &   czip,czip)
        call ggZZcapture('bubpp'//char(ichar('0')+j),h3,h4,1,2,3,4,5,6,
     &   b(2,2,h3,h4,j),
     &   czip,czip)
        enddo
        enddo
        enddo
      endif

c      if (docheck) then
c        do j=1,7
c        write(6,*) 'j,Bint(j)',j,Bint(j,0)*b(2,1,j)
c        enddo
c      endif

      return
      end



      subroutine ZZmassivebubcore(j1,j2,j3,j4,j5,j6,za,zb,mt,b)
      implicit none
      include 'types.f'
c----Bubble coefficients extracted from BDK 11.5, 11.8
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'blabels.f'
c      logical:: ggZZunstable
      character(len=9):: st,st1
      integer:: j,j1,j3,j4,j2,j5,j6,k1,k3,k4,k2,k5,k6,h1,h2,h3,h4
      complex(dp):: b(2,2,2,2,7),bcoeff(8),tmp
      real(dp):: mt
c      common/ggZZunstable/ggZZunstable
c--- QCDLoop already initialized from call to massivebox

      k1=j1
      k2=j2

      do h1=1,2
c      do h2=1,2
c--- only compute h2=1: will obtain others by c.c.
      h2=2
      do h3=1,2
      do h4=1,2

      if (h3 ==1) then
          k3=j3
          k4=j4
      elseif (h3 ==2) then
          k3=j4
          k4=j3
      endif
      if (h4 ==1) then
          k5=j5
          k6=j6
      elseif (h4 ==2) then
          k5=j6
          k6=j5
      endif

      if ((h1==1) .and. (h2==1)) st='q+qb-g-g-'
      if ((h1==1) .and. (h2==2)) st='q+qb-g-g+'
      if ((h1==2) .and. (h2==1)) st='q+qb-g+g-'
      if ((h1==2) .and. (h2==2)) st='q+qb-g+g+'

      if (st == 'q+qb-g-g-') then
          st1='q+qb-g+g+'
          call ZZmbc(st1,k2,k1,k4,k3,k6,k5,zb,za,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
          bcoeff(b134)=bcoeff(b234)
          bcoeff(b234)=tmp
c--- (-,-) amplitudes
c      write(6,*) '--'
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b1)',
c     & bcoeff(b1),abs(bcoeff(b1))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b1)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
c      stop

      elseif (st=='q+qb-g-g+') then
          call ZZmbc(st,k1,k2,k3,k4,k5,k6,za,zb,bcoeff)

c--- new:(-,+) amplitudes
c      write(6,*) '-+'
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b1)',
c     & bcoeff(b1),abs(bcoeff(b1))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b1)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
c      stop

      elseif (st=='q+qb-g+g+') then
          call ZZmbc(st,k1,k2,k3,k4,k5,k6,za,zb,bcoeff)
c--- (+,+) amplitudes
c      write(6,*) '++'
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b1)',
c     & bcoeff(b1),abs(bcoeff(b1))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b1)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
c      stop

      elseif (st=='q+qb-g+g-') then
          st1='q+qb-g+g-'
          call ZZmbc(st1,k1,k2,k4,k3,k6,k5,zb,za,bcoeff)
c--- (+,-) amplitudes
c      write(6,*) 'old:+-'
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b1)',
c     & bcoeff(b1),abs(bcoeff(b1))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b1)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
c      pause

      endif

      do j=1,7
      b(h1,h2,h3,h4,j)=bcoeff(j)
      enddo

      enddo

      enddo
c      enddo
      enddo

      return
      end

