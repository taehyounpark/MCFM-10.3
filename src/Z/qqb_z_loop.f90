!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      ! F_ns: non-singlet form factors at 1,2,3 loop
      ! Fvec_s: singlet vector form factors at 1,2,3 loop
      ! Fax_s: singlet axial-vector form factors at 1,2,3 loop

module qqb_z_hardfun
    use types
    implicit none

    public :: qqb_z_loop

    contains

!---- Hard function multi-loop matrix elements for Z production
!----averaged over initial colours and spins
!     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))

      ! results in msqall returned as coefficients in series (alphas/4/pi)^k

      subroutine qqb_z_loop(p, musq, msqall, withaxsinglet)
          use types
      implicit none
      
      include 'src/Inc/constants.f'
      include 'src/Inc/nf.f'
      include 'src/Inc/mxpart.f'
      include 'src/Inc/cplx.h'
      include 'src/Inc/masses.f'
      include 'src/Inc/ewcouple.f'
      include 'src/Inc/zcouple.f'
      include 'src/Inc/ewcharge.f'
      include 'src/Inc/zprods_decl.f'
      integer:: j,k
      real(dp), intent(in) :: p(mxpart,4), musq
      logical, intent(in), optional :: withaxsinglet

      logical :: withaxsinglet_use

      ! loop corrections to vector and axial-vector form-factors with
      ! non-singlet (ns) and singlet (s) contributions each
      complex(dp) :: F_ns(3), Fvec_s(3), Fax_s(3)
      real(dp), intent(out) :: msqall(-nf:nf,-nf:nf,0:3)
      real(dp):: s,fac,s34
      complex(dp):: prop,qqb,qbq,Zprop

      complex(dp) :: amptree(4), amp1loop(4), amp2loop(4), amp3loop(4)
      real(dp) :: prefac
      integer :: i1,i2
      
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      s34=s(3,4)

      F_ns(:) = 0._dp
      Fvec_s(:) = 0._dp
      Fax_s(:) = 0._dp
      msqall(:,:,:) = 0._dp

      if (present(withaxsinglet)) then
          withaxsinglet_use = withaxsinglet
      else
          withaxsinglet_use = .true.
      endif

      call qqcoeff3(s(1,2), musq, F_ns, Fvec_s(3))
      ! DEBUG to disable vector singlet contribution
      !Fvec_s(:) = 0._dp

      ! comment out to disable axial singlet contributions
      if (withaxsinglet_use .eqv. .false.) then
          Fax_s(:) = 0._dp
      else
          call qqcoeff3_ax(s(1,2), musq, Fax_s)
      endif

      Zprop = (s34-zmass**2)+im*zmass*zwidth

      call spinoru(4,p,za,zb)

      ! helicity amplitudes:
      ! (1) = -+-+ PL (quark part), PL (leptons part)
      ! (2) = +--+ PR, PL

      ! (3) = +-+- PR, PR
      ! (4) = -++- PL, PR

      if (s(1,2) > 4*mt**2) then
          write (*,*) "Warning: Q^2 > 4*mt^2, 3-loop amplitude no longer valid"
      endif

      ! prefactor for squared amplitudes
      prefac = aveqq*esq**2*xn

      do k=-nf,nf
          if (k == 0) then
            continue
          else
            if (k > 0) then
                j = k
                i1 = 1
                i2 = 2
            else
                ! charge conjugation for qbar q case
                j = -k
                i1 = 2
                i2 = 1
            endif

            amptree(1) = 2._dp*za(i2,3)*zb(4,i1)*(Q(j)*q1/s34+L(j)*l1/Zprop)
            amptree(2) = 2._dp*za(i1,3)*zb(4,i2)*(Q(j)*q1/s34+R(j)*l1/Zprop)
            amptree(3) = 2._dp*za(i1,4)*zb(3,i2)*(Q(j)*q1/s34+R(j)*r1/Zprop)
            amptree(4) = 2._dp*za(i2,4)*zb(i1,3)*(Q(j)*q1/s34+L(j)*r1/Zprop)

            ! no singlet contribution at one-loop
            amp1loop(1) = 2._dp*za(i2,3)*zb(4,i1)*(Q(j)*q1/s34*F_ns(1) + L(j)*F_ns(1)*l1/Zprop)
            amp1loop(2) = 2._dp*za(i1,3)*zb(4,i2)*(Q(j)*q1/s34*F_ns(1) + R(j)*F_ns(1)*l1/Zprop)
            amp1loop(3) = 2._dp*za(i1,4)*zb(3,i2)*(Q(j)*q1/s34*F_ns(1) + R(j)*F_ns(1)*r1/Zprop)
            amp1loop(4) = 2._dp*za(i2,4)*zb(i1,3)*(Q(j)*q1/s34*F_ns(1) + L(j)*F_ns(1)*r1/Zprop)

            ! at two-loops there's an axial singlet contribution
            amp2loop(1) = 2._dp*za(i2,3)*zb(4,i1)*( &
                            F_ns(2)*Q(j)*q1/s34 + &
                            F_ns(2)*L(j)*l1/Zprop + & ! non-singlet contribution
                            Fax_s(2)*(L(5)-R(5))/2._dp * l1/Zprop & ! third generation axial singlet
                             ! Fax_s = Fax_s(bottom) - Fax_s(top)
                          )
            amp2loop(2) = 2._dp*za(i1,3)*zb(4,i2)*( &
                            F_ns(2)*Q(j)*q1/s34 + &
                            F_ns(2)*R(j)*l1/Zprop + &
                            Fax_s(2)*(R(5)-L(5))/2._dp * l1/Zprop &
                          )
            amp2loop(3) = 2._dp*za(i1,4)*zb(3,i2)*( &
                            F_ns(2)*Q(j)*q1/s34 + &
                            F_ns(2)*R(j)*r1/Zprop + &
                            Fax_s(2)*(R(5)-L(5))/2._dp * r1/Zprop &
                          )
            amp2loop(4) = 2._dp*za(i2,4)*zb(i1,3)*( &
                            F_ns(2)*Q(j)*q1/s34 + &
                            F_ns(2)*L(j)*r1/Zprop + &
                            Fax_s(2)*(L(5)-R(5))/2._dp * r1/Zprop &
                          )

            ! at three-loops there's both vector and axial-vector singlet
            ! contributions. For the vector part there is also a top-loop, but
            ! it is power suppressed and a tiny contribution of a small contribution itself
            amp3loop(1) = 2._dp*za(i2,3)*zb(4,i1)*( &
                            F_ns(3)*Q(j)*q1/s34 + & ! photon non-singlet contribution
                            Fvec_s(3)*sum(Q(1:nf))*q1/s34 + & ! photon singlet contribution
                            F_ns(3)*L(j)*l1/Zprop + & ! Z non-singlet contribution
                            Fvec_s(3)*sum(L(1:nf)+R(1:nf))/2._dp*l1/Zprop + & ! Z vector singlet
                            Fax_s(3)*(L(5)-R(5))/2._dp * l1/Zprop & ! third generation axial singlet
                          )
            amp3loop(2) = 2._dp*za(i1,3)*zb(4,i2)*( &
                            F_ns(3)*Q(j)*q1/s34 + &
                            Fvec_s(3)*sum(Q(1:nf))*q1/s34 + & ! photon singlet contribution
                            F_ns(3)*R(j)*l1/Zprop + &
                            Fvec_s(3)*sum(L(1:nf)+R(1:nf))/2._dp*l1/Zprop + & ! Z vector singlet
                            Fax_s(3)*(R(5)-L(5))/2._dp * l1/Zprop &
                          )
            amp3loop(3) = 2._dp*za(i1,4)*zb(3,i2)*( &
                            F_ns(3)*Q(j)*q1/s34 + &
                            Fvec_s(3)*sum(Q(1:nf))*q1/s34 + & ! photon singlet contribution
                            F_ns(3)*R(j)*r1/Zprop + &
                            Fvec_s(3)*sum(L(1:nf)+R(1:nf))/2._dp*r1/Zprop + & ! Z vector singlet
                            Fax_s(3)*(R(5)-L(5))/2._dp * r1/Zprop &
                          )
            amp3loop(4) = 2._dp*za(i2,4)*zb(i1,3)*( &
                            F_ns(3)*Q(j)*q1/s34 + &
                            Fvec_s(3)*sum(Q(1:nf))*q1/s34 + & ! photon singlet contribution
                            F_ns(3)*L(j)*r1/Zprop + &
                            Fvec_s(3)*sum(L(1:nf)+R(1:nf))/2._dp*r1/Zprop + & ! Z vector singlet
                            Fax_s(3)*(L(5)-R(5))/2._dp * r1/Zprop &
                          )

            msqall(k,-k,0) = prefac*sum(abs(amptree(:))**2)
            msqall(k,-k,1) = prefac*2._dp*sum(real(amptree(:)*conjg(amp1loop(:))))
            msqall(k,-k,2) = prefac*sum(real(amp1loop(:)*conjg(amp1loop(:)))) + &
                            prefac*2._dp*sum(real(amptree(:)*conjg(amp2loop(:))))
            msqall(k,-k,3) = prefac*2._dp*sum(real(amp1loop(:)*conjg(amp2loop(:)))) + &
                            prefac*2._dp*sum(real(amptree(:)*conjg(amp3loop(:))))

          endif
      enddo

      return
      end

end module
