!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_gamgam_frag(p,msq)
      implicit none
      include 'types.f'
c=====
c-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+g(p4)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'frag.f'
      include 'zcouple_cms.f'

      integer:: j,k,i
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qa,aq,qg,gq,ag,ga
      real(dp):: D(0:5),fsq
      common/D/D
!$omp threadprivate(/D/)

      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0._dp
         if     (fragset == 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))
         elseif (fragset == 'BFGsetII') then
            call get_frag(z_frag,fsq,2,i,D(i))
         elseif (fragset == 'GdRG__LO') then
            call GGdR_frag(z_frag,i,D(i),0)
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop
         endif
      enddo



      call dotem(3,p,s)
      fac=4._dp*V*gsq*abs(zesq)

      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa
      qg=-fac*aveqg*(s(1,3)/s(1,2)+s(1,2)/s(1,3))
      ag=qg
      gq=-fac*aveqg*(s(1,2)/s(2,3)+s(2,3)/s(1,2))
      ga=gq

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qa
      if ((j > 0) .and. (k < 0)) then
          if (j == -k) msq(j,k)=Q(j)**2*qa*D(0)
c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) msq(j,k)=Q(k)**2*aq*D(0)
c--qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*qg*D(abs(j))
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=Q(j)**2*ag*D(abs(j))
c--gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=Q(k)**2*gq*D(abs(k))
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=Q(k)**2*ga*D(abs(k))
      endif

      enddo
      enddo


      return
      end



