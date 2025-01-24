!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine vv_amps(p,msq,order)
      use vv_msq
      use handyG
      use mod_vvamp
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'nf.f'
      include 'zcouple_cms.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'newpoint.f'
      include 'nproc.f'
      include 'nflav.f'
      integer:: i,j,order,fmax,ztype1,ztype2
      integer,parameter:: dntype=1,uptype=2
      integer,parameter:: eltype=1,nutype=2
      real(dp)::p(mxpart,4),msq(-nf:nf,-nf:nf,0:2),s12,t,u,p3sq,p4sq
      complex(dp)::qqbAj(0:2,4,10)
      logical::qqb,reset0,reset1,reset2

c Convert MCFM virtual result from dred into tH-V (equivalent to CDR) scheme
c      msqvMCFM(:)=msqvMCFM(:)-ason2pi*CF*(1d0-zeta2)*msqMCFM(:)

      msq(:,:,:)=zip

c    p all outgoing 0 -> 1+2+3+4+5+6
c Signs changed in subsequent routines, qqbVVMLLL and qqbVVMRLL
c so that crossing p1->-p1,p2->-p2 no longer required
      call spinoru(6,p,za,zb)

c assign correct lepton types for Z decays
      ztype1=eltype
      if (nproc == 82) then
        ztype2=nutype
      else
        ztype2=eltype
      endif

      reset0=.true.
      reset1=.true.
      reset2=.true.

      do i=-2,2
c Only fill these arrays twice (qqb and qbq) per call
        if (i == -2) then
          if (newpoint) then
            call clearcache
            s12=s(1,2)
            t=s(1,3)+s(1,4)+s(3,4)
            u=s(2,3)+s(2,4)+s(3,4)
            p3sq=s(3,4)
            p4sq=s(5,6)
            call qqbAjfill(s12,u,t,p3sq,p4sq,order,nflav,qqbAj)   ! qb-q case
          endif
          qqb=.false.
        elseif (i == 1) then
          if (newpoint) then
            call clearcache
            s12=s(1,2)
            t=s(1,3)+s(1,4)+s(3,4)
            u=s(2,3)+s(2,4)+s(3,4)
            p3sq=s(3,4)
            p4sq=s(5,6)
            call qqbAjfill(s12,t,u,p3sq,p4sq,order,nflav,qqbAj)   ! q-qb case
          endif
          qqb=.true.
        endif

      do j=-2,2
          if (j == 0) cycle
          if ((kcase==kWWqqbr .or. kcase==kZZlept) .and. i /= -j) cycle
          if (kcase==kWZbbar .and. (i+j /= nwz)) cycle
c no bottom-initiated contribution for WW
          if (kcase==kWWqqbr .and. ((abs(i)==5) .or. (abs(j)==5))) cycle

          if (i < 0) then       ! ub-u contributions
                qqb=.false.
              if (mod(i,2) == 0) then
                msq(i,j,0)=msq_tree(qqb,uptype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset0)
                if (order >= 1) then
                msq(i,j,1)=msq_onelp(qqb,uptype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset1)
                endif
                if (order >= 2) then
                msq(i,j,2)=msq_twolp(qqb,uptype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset2)
                endif

              else               ! db-d contributions
                msq(i,j,0)=msq_tree(qqb,dntype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset0)
                if (order >= 1) then
                msq(i,j,1)=msq_onelp(qqb,dntype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset1)
                endif
                if (order >= 2) then
                msq(i,j,2)=msq_twolp(qqb,dntype,ztype1,ztype2,2,1,3,4,5,6,qqbAj,reset2)
                endif
              endif
          elseif (i > 0) then
c u-ub contributions
              qqb=.true.
              if ( ((mod(i,2) == 0) .and. (kcase /= kWZbbar))
     &        .or. ((mod(i,2) == 1) .and. (kcase == kWZbbar)) ) then
                msq(i,j,0)=msq_tree(qqb,uptype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset0)
                if (order >= 1) then
                msq(i,j,1)=msq_onelp(qqb,uptype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset1)
                endif
                if (order >= 2) then
                msq(i,j,2)=msq_twolp(qqb,uptype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset2)
                endif
              else
c d-db contributions
                msq(i,j,0)=msq_tree(qqb,dntype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset0)
                if (order >= 1) then
                msq(i,j,1)=msq_onelp(qqb,dntype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset1)
                endif
                if (order >= 2) then
                msq(i,j,2)=msq_twolp(qqb,dntype,ztype1,ztype2,1,2,3,4,5,6,qqbAj,reset2)
                endif
              endif
            endif
        enddo
      enddo

c Extend to other flavours
      if (kcase == kWWqqbr) then
        fmax=4
      else
        fmax=5
      endif

      if (kcase == kWZbbar) then
        do j=3,fmax,2
        msq(j,-2,:)=msq(1,-2,:)
        msq(-2,j,:)=msq(-2,1,:)
        msq(2,-j,:)=msq(2,-1,:)
        msq(-j,2,:)=msq(-1,2,:)
        enddo
        do j=1,fmax,2
        msq(j,-4,:)=msq(1,-2,:)
        msq(-4,j,:)=msq(-2,1,:)
        msq(4,-j,:)=msq(2,-1,:)
        msq(-j,4,:)=msq(-1,2,:)
        enddo
      else
        do j=3,fmax
        msq(+j,-j,:)=msq(+(j-2),-(j-2),:)
        msq(-j,+j,:)=msq(-(j-2),+(j-2),:)
        enddo
      endif

c Apply overall couplings
      do i=-nf,nf
      do j=-nf,nf
      if (kcase==kWZbbar) then
        msq(i,j,:)=4*abs(zesq**4/zxw**2)*msq(i,j,:)*Vsq(i,j)*aveqq*xn
      elseif (kcase==kWWqqbr) then
        msq(i,j,:)=abs(zesq/zxw)**4*msq(i,j,:)*aveqq*xn
      elseif (kcase==kZZlept) then
        msq(i,j,:)=16*abs(zesq)**4*msq(i,j,:)*aveqq*xn
      endif
      enddo
      enddo

      newpoint = .false.

      return
      end

