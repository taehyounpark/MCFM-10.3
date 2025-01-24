!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zaj_frag(p,msq)
      implicit none
      include 'types.f'
 !------- FRAGMENTATION CONTRIBUTION, Always fragment p5 (=j) , allow phase space integration to symmetrize.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'frag.f'
      integer:: j,k,i
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: D(0:5),fsq
      common/D/D
      real(dp):: mqq_z2j(0:2,fn:nf,fn:nf)
      real(dp):: msqx_z2j(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: msqx_cs_z2j(0:2,-nf:nf,-nf:nf)
      real(dp):: msq_z2j(-nf:nf,-nf:nf)
      integer:: m,icol
      real(dp):: qcd_z2j(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
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


c====== NEED TO SETUP equivlient of qqii_ii etc for
c====== this process
c====== Underlying QCD matrix element
      call qqb_z2jetx(p,msq_z2j,mqq_z2j,msqx_z2j,msqx_cs_z2j)
c----- not naming follows z2jet process, we will work with msqx_z2j


c====== We do not need individual colour compenents so define new matrices which are sums over colour

      do i=-nf,nf
         do j=-nf,nf
            do k=-nf,nf
               do m=-nf,nf
                  qcd_z2j(i,j,k,m)=0._dp
                  do icol=0,2
               qcd_z2j(i,j,k,m)=qcd_z2j(i,j,k,m)+msqx_z2j(icol,i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       enddo



       do j=-nf,nf
       do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qq
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then

c---- p5 = j
            msq(j,k)=two*qcd_z2j(j,j,j,j)*D(abs(j))
          else
            msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))+qcd_z2j(j,k,k,j)*D(k)
          endif

c--qa
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then

             msq(j,k)=two*qcd_z2j(j,k,j,k)*D(j)
             do m=1,5
                if(m /= j) then
                   msq(j,k)=msq(j,k)+two*qcd_z2j(j,k,m,-m)*D(m)
                endif
             enddo
             msq(j,k)=msq(j,k)+qcd_z2j(j,k,0,0)*D(0)*two
          else
c-------p3=j
             msq(j,k)=qcd_z2j(j,k,j,k)*D(j)+qcd_z2j(j,k,k,j)*D(abs(k))
          endif

c--aa
       elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
             msq(j,k)=two*qcd_z2j(j,k,j,k)*D(abs(j))
          else
             msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))
     &            +qcd_z2j(j,k,k,j)*D(abs(k))
          endif

c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=two*(qcd_z2j(j,k,j,k))*D(abs(j))
            do m=1,5
               if(m /= -j) then
                  msq(j,k)=msq(j,k)+two*qcd_z2j(j,k,m,-m)*D(m)
               endif
            enddo
            msq(j,k)=msq(j,k)+two*(qcd_z2j(j,k,0,0)*D(0))
          else
            msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))+qcd_z2j(j,k,k,j)*D(k)
          endif

c--qg_qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qcd_z2j(j,k,j,0)*D(j)+qcd_z2j(j,k,0,j)*D(0)
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=qcd_z2j(j,k,j,0)*D(abs(j))+qcd_z2j(j,k,0,j)*D(0)
c--gq_gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=qcd_z2j(j,k,k,0)*D(k)+qcd_z2j(j,k,0,k)*D(0)
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=qcd_z2j(j,k,k,0)*D(abs(k))+qcd_z2j(j,k,0,k)*D(0)
c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=qcd_z2j(j,k,0,0)*D(0)
            do m=1,5
               msq(j,k)=msq(j,k)+qcd_z2j(0,0,m,-m)*D(m)
            enddo
      endif

      enddo
      enddo


      return
      end



