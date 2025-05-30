!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monophot_v_Samps(p,i1,i2,i3,i4,i5,qgqb)
      implicit none
      include 'types.f'

c----- combined colour and Born intefered amplitudes as a function of quark line helicity
c------ default is q(i1)+g(i2)+qb(i3)+x(i4)+x(i5)
      include 'mxpart.f'
      include 'qcdcouple.f'
c      include 'zprods_decl.f'
      real(dp):: qgqb(2)
      integer:: i1,i2,i3,i4,i5
      real(dp):: p(mxpart,4)
      integer:: h1,h2,h3,h4
      complex(dp):: amp_tree(2,2,2,2),amp_slc(2,2,2,2)

c------ tree
      call qqb_dm_monojet_Samps(p,i1,i2,i3,i4,i5,amp_tree)
c------ subleading colour (swap i2 and i3)
      call qqb_dm_monophot_Samps_fill(p,i1,i3,i2,i4,i5,amp_slc)

      qgqb(:)=0d0
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  qgqb(h1)=qgqb(h1)
     &                 +ason2pi*Dble(conjg(amp_tree(h1,h2,h3,h4))*
     &                 amp_slc(h1,h2,h3,h4))
               enddo
            enddo
         enddo
      enddo


      return
      end


      subroutine qqb_dm_monophot_Samps_fill(p,i1,i2,i3,i4,i5,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'epinv.f'
c----- SUB LEADING COLOR
c----- fills amplitude for q qb g chi,chib
      complex(dp):: amp(2,2,2,2)
      complex(dp):: amp_tree(2,2,2,2)
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
      integer:: h1,h2,h3,h4
      complex(dp):: Lsm1
      real(dp):: s12,s23,s123,s34
      complex(dp):: lnrat
      real(dp):: s(mxpart,mxpart)
      complex(dp):: vfac,ffac,rat(2,2),amp_dec(2,2)
      integer:: j,k
      real(dp):: bp,beta
      real(dp):: uvsub
      real(dp):: s13

      amp(:,:,:,:)=czip
      if(xmass>1d-8) then
c--------- generate massless phase space
      call gen_masslessvecs(p,q,i4,i5)
c--------- generate spinors
      call spinoru(5,q,za,zb)
      else
c-------- massless dm can use usual spinoru
         call spinoru(5,p,za,zb)

      endif
      s(:,:)=0d0
      do j=1,5
         do k=1,5
            s(j,k)=Dble(za(j,k)*zb(k,j))
         enddo
      enddo

c------ basis integrals
c      l12=lnrat(musq,-s(i1,i2))
c      l23=lnrat(musq,-s(i2,i3))

      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s12=s(i1,i2) ! s12 = s_qqb
      s13=s(i1,i3)
      s23=s(i2,i3)
c------ univerisal v and f pieces
      vfac=-epinv**2-(epinv*lnrat(musq,-s12))
     &-(0.5d0*lnrat(musq,-s12)**2)
      ffac=-Lsm1(-s23,-s123,-s12,-s123)-Lsm1(-s13,-s123,-s12,-s123)
      uvsub=-1.5d0*epinv
c------helicity dependent rational pieces

      rat(1,1)=0.5d0*(s13+s23)/(zb(i3,i2)*zb(i1,i3))
      rat(2,2)=0.5d0*(s13+s23)/(za(i3,i2)*za(i1,i3))
      rat(1,2)=czip
      rat(2,1)=czip


c====== scalar decay
      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

      call dm_scal_decay(i4,i5,za,zb,bp,amp_dec)
c======= tree (note swap on gluon quark in these amplitudes)
      call qqb_dm_monojet_Samps(p,i1,i3,i2,i4,i5,amp_tree)

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
               amp(h1,h2,h3,h4)=(vfac+ffac+uvsub)*amp_tree(h1,h2,h3,h4)
     &                 +rat(h1,h2)*amp_dec(h3,h4)
               enddo
            enddo
         enddo
      enddo

      return
      end
