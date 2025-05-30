!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monophot_v_Axamps(p,i1,i2,i3,i4,i5,qgqb)
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
      call qqb_dm_monojet_Axamps(p,i1,i2,i3,i4,i5,amp_tree)
c------ subleading colour (swap i2 and i3)
      call qqb_dm_monophot_Axamps_fill(p,i1,i3,i2,i4,i5,amp_slc)

      qgqb(:)=0d0
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  qgqb(h1)=qgqb(h1)
     &                 -ason2pi*Dble(conjg(amp_tree(h1,h2,h3,h4))*
     &                 amp_slc(h1,h2,h3,h4))
               enddo
            enddo
         enddo
      enddo

      return
      end


      subroutine qqb_dm_monophot_Axamps_fill(p,i1,i2,i3,i4,i5,amp)
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
      complex(dp):: amp(2,2,2,2),Vamp(2,2,2,2),Famp(2,2,2,2)
      complex(dp):: amp_tree(2,2,2,2)
      real(dp):: p(mxpart,4),q(mxpart,4)
      integer:: i1,i2,i3,i4,i5
      integer:: h1,h2,h3,h4
      complex(dp):: l12,l45,L0,Lsm1,L1
      complex(dp):: lnrat
      real(dp):: s(mxpart,mxpart)
      complex(dp):: vfac
      integer:: j,k
      real(dp):: s12,s23,s13,s45
      real(dp):: bp,beta,s34

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
      s(:,:)=zip
      do j=1,5
         do k=1,5
            s(j,k)=Dble(za(j,k)*zb(k,j))
         enddo
      enddo
c beta = sqrt(1-4xmass**2/s34)
      s34=Dble(za(i4,i5)*zb(i5,i4))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)


      s12=s(i1,i2)
      s23=s(i2,i3)
      s13=s(i1,i3)
      s45=s(i4,i5)

c------ basis integrals
      l12=lnrat(musq,-s(i1,i2))
      l45=lnrat(musq,-s(i4,i5))

c----- use same notation as Born
c----- split into Vamp and Famp with amp=Vamp+Famp

c------ Vamp piece universial
      vfac=
     &-(epinv**2+epinv*l12+0.5d0*l12**2)
     & -2d0*(epinv+l45)-4d0
     & +0.5d0*(epinv+l45)+0.5d0

c------ note have swapped order of i3 and i2 to calculate proper tree
c------ need minus sign too
      call qqb_dm_monojet_Axamps(p,i1,i3,i2,i4,i5,amp_tree)
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  Vamp(h1,h2,h3,h4)=-amp_tree(h1,h2,h3,h4)*vfac
               enddo
            enddo
         enddo
      enddo

c      write(6,*) 'amp_tree ',amp_tree(2,2,1,2)/s(i4,i5)

c------ Fpieces
c------ Fpieces

c====== Helicity conserving

c------ Fpieces

c====== Helicity conserving

      Famp(1,1,1,2)= -(za(i1,i3)*za(i1,i4)*zb(i5,i1))/(2.*za(i1,i2)*
     &     zb(i3,i1)) +
     &  (bp*za(i1,i3)*za(i1,i4)*zb(i5,i1))/(za(i1,i2)*zb(i3,i1)) -
     &  (za(i1,i3)**2*za(i2,i4)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*zb(i3,i1)) +
     &  (bp*za(i1,i3)**2*za(i2,i4)*zb(i5,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i3,i1)) +
     &  (L0(-s45,-s12)*za(i1,i3)*za(i1,i4)*zb(i2,i1)*zb(i5,i1))/
     &   (s12*zb(i3,i1)) - (2*bp*L0(-s45,-s12)*za(i1,i3)*za(i1,i4)*
     &     zb(i2,i1)*zb(i5,i1))/(s12*zb(i3,i1)) -
     &  (L1(-s45,-s23)*za(i1,i3)**2*za(i3,i4)*zb(i3,i2)*zb(i5,i1))/
     &   (2.*s23**2) + (bp*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i4)*
     &     zb(i3,i2)*zb(i5,i1))/s23**2 +
     &  (L0(-s45,-s12)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*zb(i3,i1)**2) -
     &  (2*bp*L0(-s45,-s12)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*zb(i3,i1)**2) +
     &  (L0(-s45,-s23)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s23*zb(i3,i1)**2) -
     &  (2*bp*L0(-s45,-s23)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s23*zb(i3,i1)**2) +
     &  (L0(-s45,-s23)*za(i1,i3)*za(i3,i4)*zb(i3,i2)*zb(i5,i1))/
     &   (s23*zb(i3,i1)) - (2*bp*L0(-s45,-s23)*za(i1,i3)*za(i3,i4)*
     &     zb(i3,i2)*zb(i5,i1))/(s23*zb(i3,i1)) -
     &  (L1(-s45,-s23)*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(2.*s23**2*zb(i3,i1)) +
     &  (bp*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s23**2*zb(i3,i1)) -
     &  (2*L0(-s23,-s45)*za(i1,i3)*za(i3,i4)*zb(i5,i2))/s45 +
     &  (4*bp*L0(-s23,-s45)*za(i1,i3)*za(i3,i4)*zb(i5,i2))/s45 -
     &  (za(i1,i4)*za(i2,i3)*zb(i5,i2))/(2.*za(i1,i2)*zb(i3,i1)) +
     &  (bp*za(i1,i4)*za(i2,i3)*zb(i5,i2))/(za(i1,i2)*zb(i3,i1)) -
     &  (za(i1,i3)*za(i2,i4)*zb(i5,i2))/(2.*za(i1,i2)*zb(i3,i1)) +
     &  (bp*za(i1,i3)*za(i2,i4)*zb(i5,i2))/(za(i1,i2)*zb(i3,i1)) +
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i3,i4)*zb(i5,i2))/zb(i3,i1) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i4)*zb(i5,i2))/zb(i3,i1)-
     &  (2*L0(-s23,-s45)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i5,i2))/
     &   (s45*zb(i3,i1)) + (4*bp*L0(-s23,-s45)*za(i1,i3)*za(i2,i4)*
     &     zb(i2,i1)*zb(i5,i2))/(s45*zb(i3,i1)) -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)*zb(i5,i2))/
     &   (zb(i3,i1)*zb(i3,i2)) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)*zb(i5,i2))/
     &   (zb(i3,i1)*zb(i3,i2)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)**2*zb(i5,i3))/
     &   zb(i3,i1)**3 + (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*
     &     zb(i2,i1)**2*zb(i5,i3))/zb(i3,i1)**3 +
     &  (s45*L1(-s45,-s12)*za(i1,i3)*za(i3,i4)*zb(i2,i1)*zb(i5,i3))/
     &   (s12**2*zb(i3,i1)) -
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i3)*za(i3,i4)*zb(i2,i1)*
     &     zb(i5,i3))/
     &   (s12**2*zb(i3,i1)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i2,i1)**2*zb(i5,i3))/
     &   (zb(i3,i1)**2*zb(i3,i2)) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i2,i1)**2*
     &     zb(i5,i3))/(zb(i3,i1)**2*zb(i3,i2))
      Famp(1,2,1,2)= (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i5,i1))/
     &   za(i2,i3)**3 - (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*
     &     za(i3,i4)*zb(i5,i1))/za(i2,i3)**3 -
     &  (2*L0(-s13,-s45)*za(i1,i2)*za(i1,i4)*zb(i3,i2)*zb(i5,i1))/
     &   (s45*za(i2,i3)) + (4*bp*L0(-s13,-s45)*za(i1,i2)*za(i1,i4)*
     &     zb(i3,i2)*zb(i5,i1))/(s45*za(i2,i3)) +
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*za(i2,i3)**2) -
     &  (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*za(i2,i3)**2) +
     &  (L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i1))/(s13*za(i2,i3)**2) -
     &  (2*bp*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i5,i1))/(s13*za(i2,i3)**2) -
     &  (za(i1,i4)*zb(i3,i2)*zb(i5,i1))/(2.*za(i2,i3)*zb(i2,i1)) +
     &  (bp*za(i1,i4)*zb(i3,i2)*zb(i5,i1))/(za(i2,i3)*zb(i2,i1)) -
     &  (L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)**2*
     &     zb(i5,i1))/(2.*s13**2*za(i2,i3)) +
     &  (bp*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)**2*
     &     zb(i5,i1))/(s13**2*za(i2,i3)) -
     &  (za(i2,i4)*zb(i3,i2)**2*zb(i5,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i3,i1)) +
     &  (bp*za(i2,i4)*zb(i3,i2)**2*zb(i5,i1))/
     &   (za(i2,i3)*zb(i2,i1)*zb(i3,i1)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i4)*zb(i5,i2))/
     &   (za(i1,i3)*za(i2,i3)) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i4)*zb(i5,i2))/
     &   (za(i1,i3)*za(i2,i3)) +
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i4)*zb(i5,i2))/
     &   (za(i1,i3)*za(i2,i3)**2) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i5,i2))/(za(i1,i3)*za(i2,i3)**2) -
     &  (za(i1,i4)*zb(i3,i1)*zb(i5,i2))/(2.*za(i2,i3)*zb(i2,i1)) +
     &  (bp*za(i1,i4)*zb(i3,i1)*zb(i5,i2))/(za(i2,i3)*zb(i2,i1)) +
     &  (L0(-s45,-s12)*za(i1,i2)*za(i2,i4)*zb(i3,i2)*zb(i5,i2))/
     &   (s12*za(i2,i3)) - (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i2,i4)*
     &     zb(i3,i2)*zb(i5,i2))/(s12*za(i2,i3)) -
     &  (za(i2,i4)*zb(i3,i2)*zb(i5,i2))/(2.*za(i2,i3)*zb(i2,i1)) +
     &  (bp*za(i2,i4)*zb(i3,i2)*zb(i5,i2))/(za(i2,i3)*zb(i2,i1)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i5,i3))/za(i2,i3) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i5,i3))/za(i2,i3) +
     &  (2*L0(-s13,-s45)*za(i1,i4)*zb(i3,i2)*zb(i5,i3))/s45 -
     &  (4*bp*L0(-s13,-s45)*za(i1,i4)*zb(i3,i2)*zb(i5,i3))/s45 -
     &  (L0(-s45,-s13)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*zb(i5,i3))/
     &   (s13*za(i2,i3)) + (2*bp*L0(-s45,-s13)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)*zb(i5,i3))/(s13*za(i2,i3)) +
     &  (s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*zb(i3,i2)*zb(i5,i3))/
     &   (s12**2*za(i2,i3)) -
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*zb(i3,i2)*
     &     zb(i5,i3))/
     &   (s12**2*za(i2,i3)) +
     &  (L1(-s45,-s13)*za(i1,i3)*za(i2,i4)*zb(i3,i2)**2*zb(i5,i3))/
     &   (2.*s13**2) - (bp*L1(-s45,-s13)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)**2*zb(i5,i3))/s13**2
      Famp(2,1,1,2)= (-2*L0(-s13,-s45)*za(i2,i3)*za(i3,i4)*
     &     zb(i5,i1))/s45 +
     &  (4*bp*L0(-s13,-s45)*za(i2,i3)*za(i3,i4)*
     &     zb(i5,i1))/s45 +
     &  (za(i1,i4)*za(i2,i3)*zb(i5,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)) -
     &  (bp*za(i1,i4)*za(i2,i3)*zb(i5,i1))/
     &   (za(i1,i2)*zb(i3,i2)) +
     &  (za(i1,i3)*za(i2,i4)*zb(i5,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)) -
     &  (bp*za(i1,i3)*za(i2,i4)*zb(i5,i1))/
     &   (za(i1,i2)*zb(i3,i2)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i3,i4)*
     &     zb(i5,i1))/zb(i3,i2) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i4)*
     &     zb(i5,i1))/zb(i3,i2) +
     &  (2*L0(-s13,-s45)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     &     zb(i5,i1))/(s45*zb(i3,i2)) -
     &  (4*bp*L0(-s13,-s45)*za(i1,i4)*za(i2,i3)*
     &     zb(i2,i1)*zb(i5,i1))/(s45*zb(i3,i2)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)*
     &     zb(i5,i1))/(zb(i3,i1)*zb(i3,i2)) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*
     &     zb(i2,i1)*zb(i5,i1))/(zb(i3,i1)*zb(i3,i2))
     &   - (L1(-s45,-s13)*za(i2,i3)**2*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i2))/(2.*s13**2) +
     &  (bp*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i2))/s13**2 -
     &  (L0(-s45,-s12)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     &     zb(i3,i1)*zb(i5,i2))/(s12*zb(i3,i2)**2) +
     &  (2*bp*L0(-s45,-s12)*za(i1,i4)*za(i2,i3)*
     &     zb(i2,i1)*zb(i3,i1)*zb(i5,i2))/
     &   (s12*zb(i3,i2)**2) -
     &  (L0(-s45,-s13)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     &     zb(i3,i1)*zb(i5,i2))/(s13*zb(i3,i2)**2) +
     &  (2*bp*L0(-s45,-s13)*za(i1,i4)*za(i2,i3)*
     &     zb(i2,i1)*zb(i3,i1)*zb(i5,i2))/
     &   (s13*zb(i3,i2)**2) +
     &  (za(i1,i4)*za(i2,i3)**2*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*zb(i3,i2)) -
     &  (bp*za(i1,i4)*za(i2,i3)**2*zb(i5,i2))/
     &   (za(i1,i2)*za(i1,i3)*zb(i3,i2)) +
     &  (za(i2,i3)*za(i2,i4)*zb(i5,i2))/
     &   (2.*za(i1,i2)*zb(i3,i2)) -
     &  (bp*za(i2,i3)*za(i2,i4)*zb(i5,i2))/
     &   (za(i1,i2)*zb(i3,i2)) -
     &  (L0(-s45,-s12)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*
     &     zb(i5,i2))/(s12*zb(i3,i2)) +
     &  (2*bp*L0(-s45,-s12)*za(i2,i3)*za(i2,i4)*
     &     zb(i2,i1)*zb(i5,i2))/(s12*zb(i3,i2)) +
     &  (L0(-s45,-s13)*za(i2,i3)*za(i3,i4)*zb(i3,i1)*
     &     zb(i5,i2))/(s13*zb(i3,i2)) -
     &  (2*bp*L0(-s45,-s13)*za(i2,i3)*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i2))/(s13*zb(i3,i2)) +
     &  (L1(-s45,-s13)*za(i1,i4)*za(i2,i3)**2*
     &     zb(i2,i1)*zb(i3,i1)*zb(i5,i2))/
     &   (2.*s13**2*zb(i3,i2)) -
     &  (bp*L1(-s45,-s13)*za(i1,i4)*za(i2,i3)**2*
     &     zb(i2,i1)*zb(i3,i1)*zb(i5,i2))/
     &   (s13**2*zb(i3,i2)) -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*
     &     zb(i2,i1)**2*zb(i5,i3))/zb(i3,i2)**3 +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*
     &     zb(i2,i1)**2*zb(i5,i3))/zb(i3,i2)**3 -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*
     &     zb(i2,i1)**2*zb(i5,i3))/
     &   (zb(i3,i1)*zb(i3,i2)**2) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*
     &     zb(i2,i1)**2*zb(i5,i3))/
     &   (zb(i3,i1)*zb(i3,i2)**2) -
     &  (s45*L1(-s45,-s12)*za(i2,i3)*za(i3,i4)*
     &     zb(i2,i1)*zb(i5,i3))/(s12**2*zb(i3,i2)) +
     &  (2*bp*s45*L1(-s45,-s12)*za(i2,i3)*za(i3,i4)*
     &     zb(i2,i1)*zb(i5,i3))/(s12**2*zb(i3,i2))
      Famp(2,2,1,2)=(Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i4)*
     &     zb(i5,i1))/
     &   (za(i1,i3)*za(i2,i3)) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i4)*
     &     zb(i5,i1))/(za(i1,i3)*za(i2,i3)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*zb(i5,i1))/
     &   (za(i1,i3)**2*za(i2,i3)) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i5,i1))/(za(i1,i3)**2*za(i2,i3)) -
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*zb(i3,i1)*zb(i5,i1))/
     &   (s12*za(i1,i3)) + (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*
     &     zb(i3,i1)*zb(i5,i1))/(s12*za(i1,i3)) +
     &  (za(i1,i4)*zb(i3,i1)*zb(i5,i1))/(2.*za(i1,i3)*zb(i2,i1)) -
     &  (bp*za(i1,i4)*zb(i3,i1)*zb(i5,i1))/(za(i1,i3)*zb(i2,i1)) +
     &  (za(i2,i4)*zb(i3,i2)*zb(i5,i1))/(2.*za(i1,i3)*zb(i2,i1)) -
     &  (bp*za(i2,i4)*zb(i3,i2)*zb(i5,i1))/(za(i1,i3)*zb(i2,i1)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*zb(i5,i2))/
     &   za(i1,i3)**3 - (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*
     &     za(i3,i4)*zb(i5,i2))/za(i1,i3)**3 -
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*za(i1,i3)**2) +
     &  (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*za(i1,i3)**2) -
     &  (L0(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s23*za(i1,i3)**2) +
     &  (2*bp*L0(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s23*za(i1,i3)**2) +
     &  (2*L0(-s23,-s45)*za(i1,i2)*za(i2,i4)*zb(i3,i1)*zb(i5,i2))/
     &   (s45*za(i1,i3)) - (4*bp*L0(-s23,-s45)*za(i1,i2)*za(i2,i4)*
     &     zb(i3,i1)*zb(i5,i2))/(s45*za(i1,i3)) +
     &  (za(i2,i4)*zb(i3,i1)*zb(i5,i2))/(2.*za(i1,i3)*zb(i2,i1)) -
     &  (bp*za(i2,i4)*zb(i3,i1)*zb(i5,i2))/(za(i1,i3)*zb(i2,i1)) +
     &  (L1(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i5,i2))/(2.*s23**2*za(i1,i3)) -
     &  (bp*L1(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i5,i2))/(s23**2*za(i1,i3)) +
     &  (za(i1,i4)*zb(i3,i1)**2*zb(i5,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i3,i2)) -
     &  (bp*za(i1,i4)*zb(i3,i1)**2*zb(i5,i2))/
     &   (za(i1,i3)*zb(i2,i1)*zb(i3,i2)) -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i5,i3))/za(i1,i3) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i5,i3))/
     &   za(i1,i3) - (L0(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i3))/(s23*za(i1,i3)) +
     &  (2*bp*L0(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &   (s23*za(i1,i3)) + (2*L0(-s23,-s45)*za(i2,i4)*zb(i3,i1)*
     &     zb(i5,i3))/s45 - (4*bp*L0(-s23,-s45)*za(i2,i4)*zb(i3,i1)*
     &     zb(i5,i3))/s45 - (s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*
     &     zb(i3,i1)*zb(i5,i3))/(s12**2*za(i1,i3)) +
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*zb(i3,i1)*
     &     zb(i5,i3))/(s12**2*za(i1,i3)) +
     &  (L1(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)**2*zb(i5,i3))/
     &   (2.*s23**2) - (bp*L1(-s45,-s23)*za(i1,i4)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i5,i3))/s23**2

      Famp(1,1,2,1)=(za(i1,i3)*za(i1,i5)*zb(i4,i1))/(2.*za(i1,i2)*
     &     zb(i3,i1)) -
     &  (bp*za(i1,i3)*za(i1,i5)*zb(i4,i1))/(za(i1,i2)*zb(i3,i1)) +
     &  (za(i1,i3)**2*za(i2,i5)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*zb(i3,i1)) -
     &  (bp*za(i1,i3)**2*za(i2,i5)*zb(i4,i1))/
     &   (za(i1,i2)*za(i2,i3)*zb(i3,i1)) -
     &  (L0(-s45,-s12)*za(i1,i3)*za(i1,i5)*zb(i2,i1)*zb(i4,i1))/
     &   (s12*zb(i3,i1)) + (2*bp*L0(-s45,-s12)*za(i1,i3)*za(i1,i5)*
     &     zb(i2,i1)*zb(i4,i1))/(s12*zb(i3,i1)) +
     &  (L1(-s45,-s23)*za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*zb(i4,i1))/
     &   (2.*s23**2) - (bp*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i5)*
     &     zb(i3,i2)*zb(i4,i1))/s23**2 -
     &  (L0(-s45,-s12)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*zb(i3,i1)**2) +
     &  (2*bp*L0(-s45,-s12)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*zb(i3,i1)**2) -
     &  (L0(-s45,-s23)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s23*zb(i3,i1)**2) +
     &  (2*bp*L0(-s45,-s23)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s23*zb(i3,i1)**2) -
     &  (L0(-s45,-s23)*za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i4,i1))/
     &   (s23*zb(i3,i1)) + (2*bp*L0(-s45,-s23)*za(i1,i3)*za(i3,i5)*
     &     zb(i3,i2)*zb(i4,i1))/(s23*zb(i3,i1)) +
     &  (L1(-s45,-s23)*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(2.*s23**2*zb(i3,i1)) -
     &  (bp*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s23**2*zb(i3,i1)) +
     &  (2*L0(-s23,-s45)*za(i1,i3)*za(i3,i5)*zb(i4,i2))/s45 -
     &  (4*bp*L0(-s23,-s45)*za(i1,i3)*za(i3,i5)*zb(i4,i2))/s45 +
     &  (za(i1,i5)*za(i2,i3)*zb(i4,i2))/(2.*za(i1,i2)*zb(i3,i1)) -
     &  (bp*za(i1,i5)*za(i2,i3)*zb(i4,i2))/(za(i1,i2)*zb(i3,i1)) +
     &  (za(i1,i3)*za(i2,i5)*zb(i4,i2))/(2.*za(i1,i2)*zb(i3,i1)) -
     &  (bp*za(i1,i3)*za(i2,i5)*zb(i4,i2))/(za(i1,i2)*zb(i3,i1)) -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i3,i5)*zb(i4,i2))/zb(i3,i1) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i5)*zb(i4,i2))/
     &   zb(i3,i1) + (2*L0(-s23,-s45)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*
     &     zb(i4,i2))/(s45*zb(i3,i1)) -
     &  (4*bp*L0(-s23,-s45)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i4,i2))/
     &   (s45*zb(i3,i1)) + (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*
     &     zb(i2,i1)*zb(i4,i2))/(zb(i3,i1)*zb(i3,i2)) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)*
     &     zb(i4,i2))/(zb(i3,i1)*zb(i3,i2)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)**2*zb(i4,i3))/
     &   zb(i3,i1)**3 - (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*
     &     zb(i2,i1)**2*zb(i4,i3))/zb(i3,i1)**3 -
     &  (s45*L1(-s45,-s12)*za(i1,i3)*za(i3,i5)*zb(i2,i1)*zb(i4,i3))/
     &   (s12**2*zb(i3,i1)) +
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*zb(i3,i1)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i2,i1)**2*zb(i4,i3))/
     &   (zb(i3,i1)**2*zb(i3,i2)) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i1)**2*zb(i3,i2))
      Famp(1,2,2,1)=-((Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &       zb(i4,i1))/za(i2,i3)**3) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i4,i1))/za(i2,i3)**3 +
     &  (2*L0(-s13,-s45)*za(i1,i2)*za(i1,i5)*zb(i3,i2)*zb(i4,i1))/
     &   (s45*za(i2,i3)) - (4*bp*L0(-s13,-s45)*za(i1,i2)*za(i1,i5)*
     &     zb(i3,i2)*zb(i4,i1))/(s45*za(i2,i3)) -
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*za(i2,i3)**2) +
     &  (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*za(i2,i3)**2) -
     &  (L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i1))/(s13*za(i2,i3)**2) +
     &  (2*bp*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i1))/(s13*za(i2,i3)**2) +
     &  (za(i1,i5)*zb(i3,i2)*zb(i4,i1))/(2.*za(i2,i3)*zb(i2,i1)) -
     &  (bp*za(i1,i5)*zb(i3,i2)*zb(i4,i1))/(za(i2,i3)*zb(i2,i1)) +
     &  (L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)**2*
     &     zb(i4,i1))/(2.*s13**2*za(i2,i3)) -
     &  (bp*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)**2*
     &     zb(i4,i1))/(s13**2*za(i2,i3)) +
     &  (za(i2,i5)*zb(i3,i2)**2*zb(i4,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i3,i1)) -
     &  (bp*za(i2,i5)*zb(i3,i2)**2*zb(i4,i1))/
     &   (za(i2,i3)*zb(i2,i1)*zb(i3,i1)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i5)*zb(i4,i2))/
     &   (za(i1,i3)*za(i2,i3)) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i5)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)) -
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*zb(i4,i2))/
     &   (za(i1,i3)*za(i2,i3)**2) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)**2) +
     &  (za(i1,i5)*zb(i3,i1)*zb(i4,i2))/(2.*za(i2,i3)*zb(i2,i1)) -
     &  (bp*za(i1,i5)*zb(i3,i1)*zb(i4,i2))/(za(i2,i3)*zb(i2,i1)) -
     &  (L0(-s45,-s12)*za(i1,i2)*za(i2,i5)*zb(i3,i2)*zb(i4,i2))/
     &   (s12*za(i2,i3)) + (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i2,i5)*
     &     zb(i3,i2)*zb(i4,i2))/(s12*za(i2,i3)) +
     &  (za(i2,i5)*zb(i3,i2)*zb(i4,i2))/(2.*za(i2,i3)*zb(i2,i1)) -
     &  (bp*za(i2,i5)*zb(i3,i2)*zb(i4,i2))/(za(i2,i3)*zb(i2,i1)) +
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i4,i3))/za(i2,i3) -
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i4,i3))/
     &   za(i2,i3) - (2*L0(-s13,-s45)*za(i1,i5)*zb(i3,i2)*zb(i4,i3))/
     &   s45 + (4*bp*L0(-s13,-s45)*za(i1,i5)*zb(i3,i2)*zb(i4,i3))/
     &   s45 + (L0(-s45,-s13)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i4,i3))/(s13*za(i2,i3)) -
     &  (2*bp*L0(-s45,-s13)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*zb(i4,i3))/
     &   (s13*za(i2,i3)) - (s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*
     &     zb(i3,i2)*zb(i4,i3))/(s12**2*za(i2,i3)) +
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*zb(i3,i2)*
     &     zb(i4,i3))/(s12**2*za(i2,i3)) -
     &  (L1(-s45,-s13)*za(i1,i3)*za(i2,i5)*zb(i3,i2)**2*zb(i4,i3))/
     &   (2.*s13**2) + (bp*L1(-s45,-s13)*za(i1,i3)*za(i2,i5)*
     &     zb(i3,i2)**2*zb(i4,i3))/s13**2
      Famp(2,1,2,1)= (2*L0(-s13,-s45)*za(i2,i3)*za(i3,i5)*zb(i4,i1))/s45
     &      -
     &  (4*bp*L0(-s13,-s45)*za(i2,i3)*za(i3,i5)*zb(i4,i1))/s45 -
     &  (za(i1,i5)*za(i2,i3)*zb(i4,i1))/(2.*za(i1,i2)*zb(i3,i2)) +
     &  (bp*za(i1,i5)*za(i2,i3)*zb(i4,i1))/(za(i1,i2)*zb(i3,i2)) -
     &  (za(i1,i3)*za(i2,i5)*zb(i4,i1))/(2.*za(i1,i2)*zb(i3,i2)) +
     &  (bp*za(i1,i3)*za(i2,i5)*zb(i4,i1))/(za(i1,i2)*zb(i3,i2)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i3,i5)*zb(i4,i1))/zb(i3,i2) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i5)*zb(i4,i1))/
     &   zb(i3,i2) - (2*L0(-s13,-s45)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*
     &     zb(i4,i1))/(s45*zb(i3,i2)) +
     &  (4*bp*L0(-s13,-s45)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i4,i1))/
     &   (s45*zb(i3,i2)) - (Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*
     &     zb(i2,i1)*zb(i4,i1))/(zb(i3,i1)*zb(i3,i2)) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)*
     &     zb(i4,i1))/(zb(i3,i1)*zb(i3,i2)) +
     &  (L1(-s45,-s13)*za(i2,i3)**2*za(i3,i5)*zb(i3,i1)*zb(i4,i2))/
     &   (2.*s13**2) - (bp*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i5)*
     &     zb(i3,i1)*zb(i4,i2))/s13**2 +
     &  (L0(-s45,-s12)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*zb(i3,i2)**2) -
     &  (2*bp*L0(-s45,-s12)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*zb(i3,i2)**2) +
     &  (L0(-s45,-s13)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s13*zb(i3,i2)**2) -
     &  (2*bp*L0(-s45,-s13)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s13*zb(i3,i2)**2) -
     &  (za(i1,i5)*za(i2,i3)**2*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*zb(i3,i2)) +
     &  (bp*za(i1,i5)*za(i2,i3)**2*zb(i4,i2))/
     &   (za(i1,i2)*za(i1,i3)*zb(i3,i2)) -
     &  (za(i2,i3)*za(i2,i5)*zb(i4,i2))/(2.*za(i1,i2)*zb(i3,i2)) +
     &  (bp*za(i2,i3)*za(i2,i5)*zb(i4,i2))/(za(i1,i2)*zb(i3,i2)) +
     &  (L0(-s45,-s12)*za(i2,i3)*za(i2,i5)*zb(i2,i1)*zb(i4,i2))/
     &   (s12*zb(i3,i2)) - (2*bp*L0(-s45,-s12)*za(i2,i3)*za(i2,i5)*
     &     zb(i2,i1)*zb(i4,i2))/(s12*zb(i3,i2)) -
     &  (L0(-s45,-s13)*za(i2,i3)*za(i3,i5)*zb(i3,i1)*zb(i4,i2))/
     &   (s13*zb(i3,i2)) + (2*bp*L0(-s45,-s13)*za(i2,i3)*za(i3,i5)*
     &     zb(i3,i1)*zb(i4,i2))/(s13*zb(i3,i2)) -
     &  (L1(-s45,-s13)*za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(2.*s13**2*zb(i3,i2)) +
     &  (bp*L1(-s45,-s13)*za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s13**2*zb(i3,i2)) +
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)**2*zb(i4,i3))/
     &   zb(i3,i2)**3 - (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*
     &     zb(i2,i1)**2*zb(i4,i3))/zb(i3,i2)**3 +
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i2,i1)**2*zb(i4,i3))/
     &   (zb(i3,i1)*zb(i3,i2)**2) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i1)*zb(i3,i2)**2) +
     &  (s45*L1(-s45,-s12)*za(i2,i3)*za(i3,i5)*zb(i2,i1)*zb(i4,i3))/
     &   (s12**2*zb(i3,i2)) -
     &  (2*bp*s45*L1(-s45,-s12)*za(i2,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*zb(i3,i2))
      Famp(2,2,2,1)= -((Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i5)*
     &     zb(i4,i1))/
     &     (za(i1,i3)*za(i2,i3))) +
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i5)*
     &     zb(i4,i1))/(za(i1,i3)*za(i2,i3)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*zb(i4,i1))/
     &   (za(i1,i3)**2*za(i2,i3)) +
     &  (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i4,i1))/(za(i1,i3)**2*za(i2,i3)) +
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*zb(i3,i1)*zb(i4,i1))/
     &   (s12*za(i1,i3)) - (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*
     &     zb(i3,i1)*zb(i4,i1))/(s12*za(i1,i3)) -
     &  (za(i1,i5)*zb(i3,i1)*zb(i4,i1))/(2.*za(i1,i3)*zb(i2,i1)) +
     &  (bp*za(i1,i5)*zb(i3,i1)*zb(i4,i1))/(za(i1,i3)*zb(i2,i1)) -
     &  (za(i2,i5)*zb(i3,i2)*zb(i4,i1))/(2.*za(i1,i3)*zb(i2,i1)) +
     &  (bp*za(i2,i5)*zb(i3,i2)*zb(i4,i1))/(za(i1,i3)*zb(i2,i1)) -
     &  (Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*zb(i4,i2))/
     &   za(i1,i3)**3 + (2*bp*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*
     &     za(i3,i5)*zb(i4,i2))/za(i1,i3)**3 +
     &  (L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*za(i1,i3)**2) -
     &  (2*bp*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*za(i1,i3)**2) +
     &  (L0(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s23*za(i1,i3)**2) -
     &  (2*bp*L0(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s23*za(i1,i3)**2) -
     &  (2*L0(-s23,-s45)*za(i1,i2)*za(i2,i5)*zb(i3,i1)*zb(i4,i2))/
     &   (s45*za(i1,i3)) + (4*bp*L0(-s23,-s45)*za(i1,i2)*za(i2,i5)*
     &     zb(i3,i1)*zb(i4,i2))/(s45*za(i1,i3)) -
     &  (za(i2,i5)*zb(i3,i1)*zb(i4,i2))/(2.*za(i1,i3)*zb(i2,i1)) +
     &  (bp*za(i2,i5)*zb(i3,i1)*zb(i4,i2))/(za(i1,i3)*zb(i2,i1)) -
     &  (L1(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i4,i2))/(2.*s23**2*za(i1,i3)) +
     &  (bp*L1(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i4,i2))/(s23**2*za(i1,i3)) -
     &  (za(i1,i5)*zb(i3,i1)**2*zb(i4,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i3,i2)) +
     &  (bp*za(i1,i5)*zb(i3,i1)**2*zb(i4,i2))/
     &   (za(i1,i3)*zb(i2,i1)*zb(i3,i2)) +
     &  (Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i4,i3))/za(i1,i3) -
     &  (2*bp*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i4,i3))/
     &   za(i1,i3) + (L0(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i3))/(s23*za(i1,i3)) -
     &  (2*bp*L0(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i4,i3))/
     &   (s23*za(i1,i3)) - (2*L0(-s23,-s45)*za(i2,i5)*zb(i3,i1)*
     &     zb(i4,i3))/s45 + (4*bp*L0(-s23,-s45)*za(i2,i5)*zb(i3,i1)*
     &     zb(i4,i3))/s45 + (s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*
     &     zb(i3,i1)*zb(i4,i3))/(s12**2*za(i1,i3)) -
     &  (2*bp*s45*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*zb(i3,i1)*
     &     zb(i4,i3))/(s12**2*za(i1,i3)) -
     &  (L1(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)**2*zb(i4,i3))/
     &   (2.*s23**2) + (bp*L1(-s45,-s23)*za(i1,i5)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i4,i3))/s23**2
c===== Hel violating

      Famp(1,1,1,1)=(xmass*za(i1,i3)*za(i1,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)**2*za(i2,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*zb(i3,i1)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i1,i4)*zb(i2,i1)*zb(i4,i1))/
     &   (s12*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(2.*s23**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*zb(i3,i1)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s23*zb(i3,i1)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i3,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (s23*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i3,i2)*zb(i4,i1))/(2.*s23**2*zb(i3,i1)*zb(i5,i4)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i3,i4)*zb(i4,i2))/
     &   (s45*zb(i5,i4)) + (xmass*za(i1,i4)*za(i2,i3)*zb(i4,i2))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)*za(i2,i4)*zb(i4,i2))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i4)*zb(i4,i2))/
     &   (zb(i3,i1)*zb(i5,i4)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i2))/(s45*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)*
     &     zb(i4,i2))/(zb(i3,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i1)**3*zb(i5,i4)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i3)*za(i3,i4)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i1)**2*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)**2*za(i2,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*zb(i3,i1)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i1,i5)*zb(i2,i1)*zb(i5,i1))/
     &   (s12*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(2.*s23**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*zb(i3,i1)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s23*zb(i3,i1)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s23*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i3,i2)*zb(i5,i1))/(2.*s23**2*zb(i3,i1)*zb(i5,i4)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i3,i5)*zb(i5,i2))/
     &   (s45*zb(i5,i4)) + (xmass*za(i1,i5)*za(i2,i3)*zb(i5,i2))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*za(i1,i3)*za(i2,i5)*zb(i5,i2))/
     &   (2.*za(i1,i2)*zb(i3,i1)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i5)*zb(i5,i2))/
     &   (zb(i3,i1)*zb(i5,i4)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i2))/(s45*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)*
     &     zb(i5,i2))/(zb(i3,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(zb(i3,i1)**3*zb(i5,i4)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i5,i3))/(s12**2*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(zb(i3,i1)**2*zb(i3,i2)*zb(i5,i4))
      Famp(1,2,1,1)=-((xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*
     &     za(i3,i4)*
     &       zb(i4,i1))/(za(i2,i3)**3*zb(i5,i4))) +
     &  (2*xmass*L0(-s13,-s45)*za(i1,i2)*za(i1,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s45*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*za(i2,i3)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s13*za(i2,i3)**2*zb(i5,i4)) +
     &  (xmass*za(i1,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)**2*zb(i4,i1))/(2.*s13**2*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i2,i4)*zb(i3,i2)**2*zb(i4,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i4)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)**2*zb(i5,i4)) +
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i2,i4)*zb(i3,i2)*zb(i4,i2))/
     &   (s12*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i2,i4)*zb(i3,i2)*zb(i4,i2))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i4,i3))/
     &   (za(i2,i3)*zb(i5,i4)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i4)*zb(i3,i2)*zb(i4,i3))/
     &   (s45*zb(i5,i4)) + (xmass*L0(-s45,-s13)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)*zb(i4,i3))/(s13*za(i2,i3)*zb(i5,i4)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*zb(i3,i2)*
     &     zb(i4,i3))/(s12**2*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i3)*za(i2,i4)*zb(i3,i2)**2*
     &     zb(i4,i3))/(2.*s13**2*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i1))/(za(i2,i3)**3*zb(i5,i4)) +
     &  (2*xmass*L0(-s13,-s45)*za(i1,i2)*za(i1,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s45*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*za(i2,i3)**2*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s13*za(i2,i3)**2*zb(i5,i4)) +
     &  (xmass*za(i1,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*
     &     zb(i3,i2)**2*zb(i5,i1))/(2.*s13**2*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i2,i5)*zb(i3,i2)**2*zb(i5,i1))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i5)*
     &     zb(i5,i2))/(za(i1,i3)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i3)*za(i2,i3)**2*zb(i5,i4)) +
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i2,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (s12*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*za(i2,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (2.*za(i2,i3)*zb(i2,i1)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i5,i3))/
     &   (za(i2,i3)*zb(i5,i4)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i5)*zb(i3,i2)*zb(i5,i3))/
     &   (s45*zb(i5,i4)) + (xmass*L0(-s45,-s13)*za(i1,i3)*za(i2,i5)*
     &     zb(i3,i2)*zb(i5,i3))/(s13*za(i2,i3)*zb(i5,i4)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*zb(i3,i2)*
     &     zb(i5,i3))/(s12**2*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i3)*za(i2,i5)*zb(i3,i2)**2*
     &     zb(i5,i3))/(2.*s13**2*zb(i5,i4))
      Famp(2,1,1,1)=(2*xmass*L0(-s13,-s45)*za(i2,i3)*za(i3,i4)*
     &     zb(i4,i1))/
     &   (s45*zb(i5,i4)) - (xmass*za(i1,i4)*za(i2,i3)*zb(i4,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i1,i3)*za(i2,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i4)*zb(i4,i1))/
     &   (zb(i3,i2)*zb(i5,i4)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     &     zb(i4,i1))/(s45*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i1))/(zb(i3,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i4,i2))/(2.*s13**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*zb(i3,i2)**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s13)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s13*zb(i3,i2)**2*zb(i5,i4)) -
     &  (xmass*za(i1,i4)*za(i2,i3)**2*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i2,i3)*za(i2,i4)*zb(i4,i2))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2))/
     &   (s12*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s13)*za(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (s13*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i4)*za(i2,i3)**2*zb(i2,i1)*
     &     zb(i3,i1)*zb(i4,i2))/(2.*s13**2*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i2)**3*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(zb(i3,i1)*zb(i3,i2)**2*zb(i5,i4)) +
     &  (s45*xmass*L1(-s45,-s12)*za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*zb(i3,i2)*zb(i5,i4)) +
     &  (2*xmass*L0(-s13,-s45)*za(i2,i3)*za(i3,i5)*zb(i5,i1))/
     &   (s45*zb(i5,i4)) - (xmass*za(i1,i5)*za(i2,i3)*zb(i5,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i1,i3)*za(i2,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i5)*zb(i5,i1))/
     &   (zb(i3,i2)*zb(i5,i4)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*
     &     zb(i5,i1))/(s45*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i1))/(zb(i3,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i5,i2))/(2.*s13**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*zb(i3,i2)**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s13)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i5,i2))/(s13*zb(i3,i2)**2*zb(i5,i4)) -
     &  (xmass*za(i1,i5)*za(i2,i3)**2*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i5,i2))/
     &   (2.*za(i1,i2)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i2,i3)*za(i2,i5)*zb(i2,i1)*zb(i5,i2))/
     &   (s12*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L0(-s45,-s13)*za(i2,i3)*za(i3,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (s13*zb(i3,i2)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*
     &     zb(i3,i1)*zb(i5,i2))/(2.*s13**2*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(zb(i3,i2)**3*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(zb(i3,i1)*zb(i3,i2)**2*zb(i5,i4)) +
     &  (s45*xmass*L1(-s45,-s12)*za(i2,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i5,i3))/(s12**2*zb(i3,i2)*zb(i5,i4))
      Famp(2,2,1,1)=-((xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*
     &     za(i2,i4)*
     &       zb(i4,i1))/(za(i1,i3)*za(i2,i3)*zb(i5,i4))) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i1))/(za(i1,i3)**2*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (s12*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*za(i2,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i3)**3*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*za(i1,i3)**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s23*za(i1,i3)**2*zb(i5,i4)) -
     &  (2*xmass*L0(-s23,-s45)*za(i1,i2)*za(i2,i4)*zb(i3,i1)*
     &     zb(i4,i2))/(s45*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i2,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i4,i2))/(2.*s23**2*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i4)*zb(i3,i1)**2*zb(i4,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i4,i3))/
     &   (za(i1,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*zb(i4,i3))/
     &   (s23*za(i1,i3)*zb(i5,i4)) -
     &  (2*xmass*L0(-s23,-s45)*za(i2,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (s45*zb(i5,i4)) + (s45*xmass*L1(-s45,-s12)*za(i1,i2)*
     &     za(i3,i4)*zb(i3,i1)*zb(i4,i3))/(s12**2*za(i1,i3)*zb(i5,i4))
     &    - (xmass*L1(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i4,i3))/(2.*s23**2*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i5)*
     &     zb(i5,i1))/(za(i1,i3)*za(i2,i3)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i1))/(za(i1,i3)**2*za(i2,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (s12*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i3)**3*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*za(i1,i3)**2*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s23*za(i1,i3)**2*zb(i5,i4)) -
     &  (2*xmass*L0(-s23,-s45)*za(i1,i2)*za(i2,i5)*zb(i3,i1)*
     &     zb(i5,i2))/(s45*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i2,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i5,i4)) -
     &  (xmass*L1(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i5,i2))/(2.*s23**2*za(i1,i3)*zb(i5,i4)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)**2*zb(i5,i2))/
     &   (2.*za(i1,i3)*zb(i2,i1)*zb(i3,i2)*zb(i5,i4)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i5,i3))/
     &   (za(i1,i3)*zb(i5,i4)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &   (s23*za(i1,i3)*zb(i5,i4)) -
     &  (2*xmass*L0(-s23,-s45)*za(i2,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (s45*zb(i5,i4)) + (s45*xmass*L1(-s45,-s12)*za(i1,i2)*
     &     za(i3,i5)*zb(i3,i1)*zb(i5,i3))/(s12**2*za(i1,i3)*zb(i5,i4))
     &    - (xmass*L1(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i5,i3))/(2.*s23**2*zb(i5,i4))

      Famp(1,1,2,2)=(xmass*za(i1,i3)*za(i1,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*za(i1,i3)**2*za(i2,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*za(i4,i5)*zb(i3,i1)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i1,i4)*zb(i2,i1)*zb(i4,i1))/
     &   (s12*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(2.*s23**2*za(i4,i5)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*za(i4,i5)*zb(i3,i1)**2) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i4,i1))/(s23*za(i4,i5)*zb(i3,i1)**2) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i3,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (s23*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*
     &     zb(i3,i2)*zb(i4,i1))/(2.*s23**2*za(i4,i5)*zb(i3,i1)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i3,i4)*zb(i4,i2))/
     &   (s45*za(i4,i5)) + (xmass*za(i1,i4)*za(i2,i3)*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*za(i1,i3)*za(i2,i4)*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i4)*zb(i4,i2))/
     &   (za(i4,i5)*zb(i3,i1)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i2))/(s45*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)*
     &     zb(i4,i2))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(za(i4,i5)*zb(i3,i1)**3) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i3)*za(i3,i4)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(za(i4,i5)*zb(i3,i1)**2*zb(i3,i2)) +
     &  (xmass*za(i1,i3)*za(i1,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*za(i1,i3)**2*za(i2,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i2,i3)*za(i4,i5)*zb(i3,i1)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i1,i5)*zb(i2,i1)*zb(i5,i1))/
     &   (s12*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i3,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(2.*s23**2*za(i4,i5)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*za(i4,i5)*zb(i3,i1)**2) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*zb(i3,i2)*
     &     zb(i5,i1))/(s23*za(i4,i5)*zb(i3,i1)**2) -
     &  (xmass*L0(-s45,-s23)*za(i1,i3)*za(i3,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (s23*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*L1(-s45,-s23)*za(i1,i3)**2*za(i2,i5)*zb(i2,i1)*
     &     zb(i3,i2)*zb(i5,i1))/(2.*s23**2*za(i4,i5)*zb(i3,i1)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i3,i5)*zb(i5,i2))/
     &   (s45*za(i4,i5)) + (xmass*za(i1,i5)*za(i2,i3)*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*za(i1,i3)*za(i2,i5)*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i1)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i3,i5)*zb(i5,i2))/
     &   (za(i4,i5)*zb(i3,i1)) +
     &  (2*xmass*L0(-s23,-s45)*za(i1,i3)*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i2))/(s45*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)*
     &     zb(i5,i2))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(za(i4,i5)*zb(i3,i1)**3) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i5,i3))/(s12**2*za(i4,i5)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(za(i4,i5)*zb(i3,i1)**2*zb(i3,i2))
      Famp(1,2,2,2)=-((xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*
     &     za(i3,i4)*
     &       zb(i4,i1))/(za(i2,i3)**3*za(i4,i5))) +
     &  (2*xmass*L0(-s13,-s45)*za(i1,i2)*za(i1,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s45*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s12*za(i2,i3)**2*za(i4,i5)) -
     &  (xmass*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*zb(i3,i2)*
     &     zb(i4,i1))/(s13*za(i2,i3)**2*za(i4,i5)) +
     &  (xmass*za(i1,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) +
     &  (xmass*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)**2*zb(i4,i1))/(2.*s13**2*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i2,i4)*zb(i3,i2)**2*zb(i4,i1))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i4)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i3)*za(i2,i3)**2*za(i4,i5)) +
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i2,i4)*zb(i3,i2)*zb(i4,i2))/
     &   (s12*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i2,i4)*zb(i3,i2)*zb(i4,i2))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i4)*zb(i4,i3))/
     &   (za(i2,i3)*za(i4,i5)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i4)*zb(i3,i2)*zb(i4,i3))/
     &   (s45*za(i4,i5)) + (xmass*L0(-s45,-s13)*za(i1,i3)*za(i2,i4)*
     &     zb(i3,i2)*zb(i4,i3))/(s13*za(i2,i3)*za(i4,i5)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i2)*za(i3,i4)*zb(i3,i2)*
     &     zb(i4,i3))/(s12**2*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i3)*za(i2,i4)*zb(i3,i2)**2*
     &     zb(i4,i3))/(2.*s13**2*za(i4,i5)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i1))/(za(i2,i3)**3*za(i4,i5)) +
     &  (2*xmass*L0(-s13,-s45)*za(i1,i2)*za(i1,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s45*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s12*za(i2,i3)**2*za(i4,i5)) -
     &  (xmass*L0(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*zb(i3,i2)*
     &     zb(i5,i1))/(s13*za(i2,i3)**2*za(i4,i5)) +
     &  (xmass*za(i1,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) +
     &  (xmass*L1(-s45,-s13)*za(i1,i2)*za(i1,i3)*za(i2,i5)*
     &     zb(i3,i2)**2*zb(i5,i1))/(2.*s13**2*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i2,i5)*zb(i3,i2)**2*zb(i5,i1))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)*za(i1,i5)*
     &     zb(i5,i2))/(za(i1,i3)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i3)*za(i2,i3)**2*za(i4,i5)) +
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i2,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (s12*za(i2,i3)*za(i4,i5)) +
     &  (xmass*za(i2,i5)*zb(i3,i2)*zb(i5,i2))/
     &   (2.*za(i2,i3)*za(i4,i5)*zb(i2,i1)) +
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i5)*zb(i5,i3))/
     &   (za(i2,i3)*za(i4,i5)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i5)*zb(i3,i2)*zb(i5,i3))/
     &   (s45*za(i4,i5)) + (xmass*L0(-s45,-s13)*za(i1,i3)*za(i2,i5)*
     &     zb(i3,i2)*zb(i5,i3))/(s13*za(i2,i3)*za(i4,i5)) -
     &  (s45*xmass*L1(-s45,-s12)*za(i1,i2)*za(i3,i5)*zb(i3,i2)*
     &     zb(i5,i3))/(s12**2*za(i2,i3)*za(i4,i5)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i3)*za(i2,i5)*zb(i3,i2)**2*
     &     zb(i5,i3))/(2.*s13**2*za(i4,i5))
      Famp(2,1,2,2)= (2*xmass*L0(-s13,-s45)*za(i2,i3)*za(i3,i4)*
     &     zb(i4,i1))/
     &   (s45*za(i4,i5)) - (xmass*za(i1,i4)*za(i2,i3)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*za(i1,i3)*za(i2,i4)*zb(i4,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i4)*zb(i4,i1))/
     &   (za(i4,i5)*zb(i3,i2)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*
     &     zb(i4,i1))/(s45*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i4)*zb(i2,i1)*
     &     zb(i4,i1))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)) +
     &  (xmass*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i4)*zb(i3,i1)*
     &     zb(i4,i2))/(2.*s13**2*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*za(i4,i5)*zb(i3,i2)**2) +
     &  (xmass*L0(-s45,-s13)*za(i1,i4)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i4,i2))/(s13*za(i4,i5)*zb(i3,i2)**2) -
     &  (xmass*za(i1,i4)*za(i2,i3)**2*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*za(i2,i3)*za(i2,i4)*zb(i4,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) +
     &  (xmass*L0(-s45,-s12)*za(i2,i3)*za(i2,i4)*zb(i2,i1)*zb(i4,i2))/
     &   (s12*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*L0(-s45,-s13)*za(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (s13*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i4)*za(i2,i3)**2*zb(i2,i1)*
     &     zb(i3,i1)*zb(i4,i2))/(2.*s13**2*za(i4,i5)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(za(i4,i5)*zb(i3,i2)**3) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i2,i1)**2*
     &     zb(i4,i3))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)**2) +
     &  (s45*xmass*L1(-s45,-s12)*za(i2,i3)*za(i3,i4)*zb(i2,i1)*
     &     zb(i4,i3))/(s12**2*za(i4,i5)*zb(i3,i2)) +
     &  (2*xmass*L0(-s13,-s45)*za(i2,i3)*za(i3,i5)*zb(i5,i1))/
     &   (s45*za(i4,i5)) - (xmass*za(i1,i5)*za(i2,i3)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*za(i1,i3)*za(i2,i5)*zb(i5,i1))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i3,i5)*zb(i5,i1))/
     &   (za(i4,i5)*zb(i3,i2)) -
     &  (2*xmass*L0(-s13,-s45)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*
     &     zb(i5,i1))/(s45*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i2,i5)*zb(i2,i1)*
     &     zb(i5,i1))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)) +
     &  (xmass*L1(-s45,-s13)*za(i2,i3)**2*za(i3,i5)*zb(i3,i1)*
     &     zb(i5,i2))/(2.*s13**2*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*za(i4,i5)*zb(i3,i2)**2) +
     &  (xmass*L0(-s45,-s13)*za(i1,i5)*za(i2,i3)*zb(i2,i1)*zb(i3,i1)*
     &     zb(i5,i2))/(s13*za(i4,i5)*zb(i3,i2)**2) -
     &  (xmass*za(i1,i5)*za(i2,i3)**2*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i1,i3)*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*za(i2,i3)*za(i2,i5)*zb(i5,i2))/
     &   (2.*za(i1,i2)*za(i4,i5)*zb(i3,i2)) +
     &  (xmass*L0(-s45,-s12)*za(i2,i3)*za(i2,i5)*zb(i2,i1)*zb(i5,i2))/
     &   (s12*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*L0(-s45,-s13)*za(i2,i3)*za(i3,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (s13*za(i4,i5)*zb(i3,i2)) -
     &  (xmass*L1(-s45,-s13)*za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*
     &     zb(i3,i1)*zb(i5,i2))/(2.*s13**2*za(i4,i5)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(za(i4,i5)*zb(i3,i2)**3) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i2,i1)**2*
     &     zb(i5,i3))/(za(i4,i5)*zb(i3,i1)*zb(i3,i2)**2) +
     &  (s45*xmass*L1(-s45,-s12)*za(i2,i3)*za(i3,i5)*zb(i2,i1)*
     &     zb(i5,i3))/(s12**2*za(i4,i5)*zb(i3,i2))
      Famp(2,2,2,2)=-((xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*
     &     za(i2,i4)*
     &       zb(i4,i1))/(za(i1,i3)*za(i2,i3)*za(i4,i5))) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i1))/(za(i1,i3)**2*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (s12*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i4)*zb(i3,i1)*zb(i4,i1))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*za(i2,i4)*zb(i3,i2)*zb(i4,i1))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i4)*
     &     zb(i4,i2))/(za(i1,i3)**3*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s12*za(i1,i3)**2*za(i4,i5)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*
     &     zb(i4,i2))/(s23*za(i1,i3)**2*za(i4,i5)) -
     &  (2*xmass*L0(-s23,-s45)*za(i1,i2)*za(i2,i4)*zb(i3,i1)*
     &     zb(i4,i2))/(s45*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i2,i4)*zb(i3,i1)*zb(i4,i2))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*L1(-s45,-s23)*za(i1,i2)*za(i1,i4)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i4,i2))/(2.*s23**2*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i4)*zb(i3,i1)**2*zb(i4,i2))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i4)*zb(i4,i3))/
     &   (za(i1,i3)*za(i4,i5)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)*zb(i4,i3))/
     &   (s23*za(i1,i3)*za(i4,i5)) -
     &  (2*xmass*L0(-s23,-s45)*za(i2,i4)*zb(i3,i1)*zb(i4,i3))/
     &   (s45*za(i4,i5)) + (s45*xmass*L1(-s45,-s12)*za(i1,i2)*
     &     za(i3,i4)*zb(i3,i1)*zb(i4,i3))/(s12**2*za(i1,i3)*za(i4,i5))
     &    - (xmass*L1(-s45,-s23)*za(i1,i4)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i4,i3))/(2.*s23**2*za(i4,i5)) -
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i1,i2)*za(i2,i5)*
     &     zb(i5,i1))/(za(i1,i3)*za(i2,i3)*za(i4,i5)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i1))/(za(i1,i3)**2*za(i2,i3)*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (s12*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)*zb(i5,i1))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*za(i2,i5)*zb(i3,i2)*zb(i5,i1))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*Lsm1(-s12,-s45,-s23,-s45)*za(i1,i2)**2*za(i3,i5)*
     &     zb(i5,i2))/(za(i1,i3)**3*za(i4,i5)) +
     &  (xmass*L0(-s45,-s12)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s12*za(i1,i3)**2*za(i4,i5)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*
     &     zb(i5,i2))/(s23*za(i1,i3)**2*za(i4,i5)) -
     &  (2*xmass*L0(-s23,-s45)*za(i1,i2)*za(i2,i5)*zb(i3,i1)*
     &     zb(i5,i2))/(s45*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i2,i5)*zb(i3,i1)*zb(i5,i2))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)) -
     &  (xmass*L1(-s45,-s23)*za(i1,i2)*za(i1,i5)*za(i2,i3)*
     &     zb(i3,i1)**2*zb(i5,i2))/(2.*s23**2*za(i1,i3)*za(i4,i5)) -
     &  (xmass*za(i1,i5)*zb(i3,i1)**2*zb(i5,i2))/
     &   (2.*za(i1,i3)*za(i4,i5)*zb(i2,i1)*zb(i3,i2)) +
     &  (xmass*Lsm1(-s12,-s45,-s13,-s45)*za(i2,i5)*zb(i5,i3))/
     &   (za(i1,i3)*za(i4,i5)) +
     &  (xmass*L0(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)*zb(i5,i3))/
     &   (s23*za(i1,i3)*za(i4,i5)) -
     &  (2*xmass*L0(-s23,-s45)*za(i2,i5)*zb(i3,i1)*zb(i5,i3))/
     &   (s45*za(i4,i5)) + (s45*xmass*L1(-s45,-s12)*za(i1,i2)*
     &     za(i3,i5)*zb(i3,i1)*zb(i5,i3))/(s12**2*za(i1,i3)*za(i4,i5))
     &    - (xmass*L1(-s45,-s23)*za(i1,i5)*za(i2,i3)*zb(i3,i1)**2*
     &     zb(i5,i3))/(2.*s23**2*za(i4,i5))





      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp(h1,h2,h3,h4)=Vamp(h1,h2,h3,h4)+Famp(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo

c      do h1=1,2
c         do h2=1,2
c            do h3=1,2
c               do h4=1,2
c              write(6,*) h1,h2,h3,h4,Vamp(h1,h2,h3,h4),Famp(h1,h2,h3,h4)
c               enddo
c            enddo
c         enddo
c      enddo

c     call qqb_dm_monophot_Vamps_fill(p,i1,i2,i3,i4,i5,amp)
c      pause


c      call writeout(p)
c      write(6,*) 'my tree ',amp_tree(2,2,1,2)/s(i4,i5)
c      call A52(i1,i2,i3,i4,i5,za,zb)
c      write(6,*) ' my imp ', amp(2,2,1,2)/s(i4,i5)
c      pause

      return
      end


