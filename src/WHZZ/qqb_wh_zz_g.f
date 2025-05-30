!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_wh_zz_g(P,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p9)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> Z(e^-(p5),e^+(p6)) Z(mu^-(p7),mu^+(p8))
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p9)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> Z(e^-(p5),e^+(p6)) Z(mu^-(p7),mu^+(p8))
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ckm.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radi_zz
      real(dp):: qqbWHg,qbqWHg,qgWHq,gqWHq,gqbWHqb,qbgWHqb

      msq(:,:)=0._dp

      call dotem(9,p,s)

      qqbWHg=aveqq*radi_zz(1,2,9,5,6,7,8,3,4)
      qbqWHg=aveqq*radi_zz(2,1,9,5,6,7,8,3,4)
      qgWHq=-radi_zz(1,9,2,5,6,7,8,3,4)*aveqg
      gqWHq=-radi_zz(2,9,1,5,6,7,8,3,4)*aveqg

      gqbWHqb=-radi_zz(9,2,1,5,6,7,8,3,4)*aveqg
      qbgWHqb=-radi_zz(9,1,2,5,6,7,8,3,4)*aveqg

c      write(6,*) 'qqbWHg',qqbWHg
c      write(6,*) 'qbqWHg',qbqWHg
c      write(6,*) 'qbgWHqb',qbgWHqb
c      write(6,*) 'gqbWHqb',gqbWHqb
c      write(6,*) 'qgWHq',qgWHq
c      write(6,*) 'gqWHq',gqWHq


      do j=-nf,nf
      do k=-nf,nf

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWHg
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWHg
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqb
      endif

      enddo
      enddo
      return

      end


      function radi_zz(j1,j2,j3,j4,j5,j6,j7,j8,j9)
      implicit none
      include 'types.f'
      real(dp):: radi_zz
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,j8,j9
      real(dp):: s4567,s12,s13,s23,s123,prop
      real(dp):: fac,hdecay

      s4567=s(j4,j5)+s(j4,j6)+s(j4,j7)+s(j5,j6)+s(j5,j7)+s(j6,j7)
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 W propagators
      prop=       ((s123-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(j8,j9)-wmass**2)**2+(wmass*wwidth)**2)

      fac=2._dp*cf*xn*gsq*gwsq**3*wmass**2/prop

      hdecay=gwsq**3*zmass**2*4._dp*xw**2/(one-xw)*
     & ( ((l1*l2)**2+(r1*r2)**2)*s(j4,j6)*s(j5,j7)
     &  +((r1*l2)**2+(r2*l1)**2)*s(j4,j7)*s(j5,j6))
      hdecay=hdecay/((s(j4,j5)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(j6,j7)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s4567-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay
c-old
c      radi_zz=s12/s13/s23
c     & *(2._dp*s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)+s(j2,j8)*s(j3,j9))
c     & +(s(j1,j9)*s(j2,j8)+s(j2,j8)*s(j3,j9)-s(j1,j8)*s(j1,j9))/s13
c     & +(s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)-s(j2,j8)*s(j2,j9))/s23


      radi_zz=
     & (s(j1,j9)*((s12+s13)*(s(j2,j8)+s(j3,j8))-s23*s(j1,j8))
     & +s(j2,j8)*((s12+s23)*(s(j1,j9)+s(j3,j9))-s13*s(j2,j9)))/(s13*s23)

      radi_zz=fac*radi_zz
      return
      end

