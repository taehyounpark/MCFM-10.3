!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c a64ax and a65ax have had all top contributions removed, so that they
c  corresonpd to only the massless loops
c Also, a65ax has been un-symmetrized to remove the 3<->4 exchange,
c  a65ax = a65axunsym(3,4) + a65axunsym(4,3)

c color structure for gg case:

c       + [T(C4,C3,i1,i2)-d_(C3,C4)*d_(i1,i2)/N] * (
c        + A64ax(4,3)
c        )
c       + [T(C3,C4,i1,i2)-d_(C3,C4)*d_(i1,i2)/N] * (
c        + A64ax(3,4)
c        )
c       + d_(C3,C4)*d_(i1,i2) * (
c        + A65ax*N^-1
c        );

c  color-structure for g(3) photon(4) case:
c   note terms proportional to quark charge in the loop (will be _q3)
c   and one proportional to charge of external quark line (will be _q4)

c       + T(C3,i1,i2)*Qloop * (
c        + A64ax(3,4)
c        + A64ax(4,3)
c        - A65axunsym(3,4)
c        - A65axunsym(4,3)
c        )
c       + T(C3,i1,i2)*Qextl * (
c        + A65axunsym(3,4)
c        );

      module zaj_virtamps_m
      use helicities
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'scale.f'!impure
      include 'epinv.f'!impure
      include 'epinv2.f'!impure

      contains

c     virt helicity amplitude (q;1), couplings stripped,
c       receiving contributions from different leading color
c       and subleading color primitive amps, each of it
c       having a pole, and cc, sc part (see BDK paper)
      subroutine zaj_virtamp_q1(i1,i2,i3,i4,i5,i6,za,zb,amps)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_lc(2,2,2,2), amps_slc(2,2,2,2)

      call zaj_virt_a6_lc(i1,i2,i3,i4,i5,i6,za,zb,amps_lc)

      call zaj_virt_a6_slc(i1,i2,i3,i4,i5,i6,za,zb,amps_slc)

      amps = -Nc*amps_lc - 1._dp/Nc*amps_slc

      end subroutine

c     virt helicity amplitude (q;2), couplings stripped
      subroutine zaj_virtamp_q2(i1,i2,i3,i4,i5,i6,za,zb,amps)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_fvs(2,2,2,2), amps_fvf(2,2,2,2)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a64v_fvs_pp, zaj_virt_a64v_fvs_pm,
     &     zaj_virt_a64v_fvs_mp, amps_fvs)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a64v_fvf_pp, zaj_virt_a64v_fvf_pm,
     &     zaj_virt_a64v_fvf_mp, amps_fvf)

      ! apparently there is a 3<->4 symmetry,
      ! such that a factor 2 is sufficient
      amps = -2*(amps_fvs + amps_fvf)

      end subroutine

      subroutine zaj_virtamp_q3(i1,i2,i3,i4,i5,i6,za,zb,amps)
      use helicities
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_s1(2,2,2,2), amps_s2(2,2,2,2)
      complex(dp) :: amps_s3(2,2,2,2), amps_s4(2,2,2,2)

       call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &    zaj_virt_a64ax_pp, zaj_virt_a64ax_pm,
     &    zaj_virt_a64ax_mp, amps_s1)

       call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &    zaj_virt_a64ax_pp, zaj_virt_a64ax_mp,
     &    zaj_virt_a64ax_pm, amps_s2)

       call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &    zaj_virt_a65ax_pp, zaj_virt_a65ax_pm,
     &    zaj_virt_a65ax_mp, amps_s3)

       call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &    zaj_virt_a65ax_pp, zaj_virt_a65ax_mp,
     &    zaj_virt_a65ax_pm, amps_s4)

      amps = amps_s1 + amps_s2 - amps_s3 - amps_s4

      ! Additional minus sign required for axial part, c.f. xzqqgg_v.f
      amps(helMinus,:,:,:) =  - amps(helMinus,:,:,:)

      end subroutine

      subroutine zaj_virtamp_q4(i1,i2,i3,i4,i5,i6,za,zb,amps)
      use helicities
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_s3(2,2,2,2)

       call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &    zaj_virt_a65ax_pp, zaj_virt_a65ax_pm,
     &    zaj_virt_a65ax_mp, amps_s3)

      amps = amps_s3

      ! Additional minus sign required for axial part, c.f. xzqqgg_v.f
      amps(helMinus,:,:,:) =  - amps(helMinus,:,:,:)

      end subroutine

      ! assembly of q;1 leading color (lc) piece
      subroutine zaj_virt_a6_lc(i1,i2,i3,i4,i5,i6,za,zb,amps)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_pole(2,2,2,2), amps_fcc(2,2,2,2), amps_fsc(2,2,2,2)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_pole_pp, zaj_virt_a6_lc_pole_pm,
     &     zaj_virt_a6_lc_pole_mp, amps_pole)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_fcc_pp, zaj_virt_a6_lc_fcc_pm,
     &     zaj_virt_a6_lc_fcc_mp, amps_fcc)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_lc_fsc_pp, zaj_virt_a6_lc_fsc_pm,
     &     zaj_virt_a6_lc_fsc_mp, amps_fsc)

      amps = amps_pole + amps_fcc + amps_fsc

      end subroutine

      subroutine zaj_virt_a6_slc(i1,i2,i3,i4,i5,i6,za,zb,amps)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      ! helicities for particles q, l, g, gam
      complex(dp), intent(out) :: amps(2,2,2,2)
      complex(dp) :: amps_pole(2,2,2,2), amps_fcc(2,2,2,2), amps_fsc(2,2,2,2)
      complex(dp) :: amps_pole_fl(2,2,2,2), amps_fcc_fl(2,2,2,2), amps_fsc_fl(2,2,2,2)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_pole_pp, zaj_virt_a6_slc_pole_pm,
     &     zaj_virt_a6_slc_pole_mp, amps_pole)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fcc_pp, zaj_virt_a6_slc_fcc_pm,
     &     zaj_virt_a6_slc_fcc_mp, amps_fcc)

      call zaj_crossings(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fsc_pp, zaj_virt_a6_slc_fsc_pm,
     &     zaj_virt_a6_slc_fsc_mp, amps_fsc)

c       symmetrization configuration
      call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_pole_pp, zaj_virt_a6_slc_pole_mp,
     &     zaj_virt_a6_slc_pole_pm, amps_pole_fl)

      call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fcc_pp, zaj_virt_a6_slc_fcc_mp,
     &     zaj_virt_a6_slc_fcc_pm, amps_fcc_fl)

      call zaj_crossings(i1,i2,i4,i3,i5,i6,za,zb,
     &     zaj_virt_a6_slc_fsc_pp, zaj_virt_a6_slc_fsc_mp,
     &     zaj_virt_a6_slc_fsc_pm, amps_fsc_fl)

      amps = amps_pole + amps_fcc + amps_fsc
     &      + amps_pole_fl + amps_fcc_fl + amps_fsc_fl

      end subroutine




c     primitive amplitudes for q;1 subleading color piece
      function zaj_virt_a6_slc_pole_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_pole_pp
      complex(dp) :: zaj_virt_a6_slc_tree_pp
      complex(dp) :: Vcc, Vsc

      complex(dp) :: xl12, xl56
      complex(dp) :: lnrat

      real(dp) :: s,t
      s(i1,i2)=real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

      xl12=Lnrat(musq,real(-s(i1,i2),dp))
      xl56=Lnrat(musq,real(-s(i5,i6),dp))

      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      zaj_virt_a6_slc_tree_pp =
     &    (-za(i2,i5)**2)/(za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))

      zaj_virt_a6_slc_pole_pp = zaj_virt_a6_slc_tree_pp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_pole_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_pole_pm
      complex(dp) :: zaj_virt_a6_slc_tree_pm
      complex(dp) :: Vcc, Vsc

      complex(dp) :: xl12, xl56
      complex(dp) :: lnrat

      real(dp) :: s,t
      s(i1,i2) = real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

      xl12=Lnrat(musq,real(-s(i1,i2),dp))
      xl56=Lnrat(musq,real(-s(i5,i6),dp))

      zaj_virt_a6_slc_tree_pm =
     &    -(-((za(i2,i4)*za(i2,i5)*zb(i1,i3)*zb(i1,i6))/
     &    (s(i3,i4)*s(i5,i6)*za(i2,i3)*zb(i1,i4)))+
     &    (za(i2,i5)*zb(i1,i3)**2*(-(za(i1,i4)*zb(i1,i6))-za(i3,i4)*zb(i3,i6
     &    )))/
     &    (s(i3,i4)*s(i5,i6)*zb(i1,i4)*t(i1,i3,i4))+
     &    (za(i2,i4)**2*zb(i1,i6)*(-(za(i2,i5)*zb(i2,i3))+za(i4,i5)*zb(i3,i4
     &    )))/
     &    (s(i3,i4)*s(i5,i6)*za(i2,i3)*t(i2,i3,i4)))

      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      zaj_virt_a6_slc_pole_pm = zaj_virt_a6_slc_tree_pm*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_pole_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_pole_mp
      complex(dp) :: zaj_virt_a6_slc_tree_mp
      complex(dp) :: Vcc, Vsc

      complex(dp) :: xl12, xl56
      complex(dp) :: lnrat

      real(dp) :: s,t
      s(i1,i2) = real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

      xl12=Lnrat(musq,real(-s(i1,i2),dp))
      xl56=Lnrat(musq,real(-s(i5,i6),dp))

      zaj_virt_a6_slc_tree_mp =
     &    (((-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*s(i5,i6)*za(i1,i4)*zb(i2,i3))-
     &    (za(i1,i3)*za(i2,i5)*zb(i1,i4)*
     &    (-(za(i1,i3)*zb(i1,i6))+za(i3,i4)*zb(i4,i6)))/
     &    (s(i3,i4)*s(i5,i6)*za(i1,i4)*t(i1,i3,i4))-
     &    (za(i2,i3)*zb(i1,i6)*zb(i2,i4)*
     &    (-(za(i2,i5)*zb(i2,i4))-za(i3,i5)*zb(i3,i4)))/
     &    (s(i3,i4)*s(i5,i6)*zb(i2,i3)*t(i2,i3,i4)))

      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      zaj_virt_a6_slc_pole_mp = zaj_virt_a6_slc_tree_mp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_slc_fcc_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_fcc_pp

      include 'heldefs.f'
      complex(dp) :: fcc

      zaj_virt_a6_slc_fcc_pp=Fcc(hqpqbmgpgp,i1,i2,i3,i4,i5,i6,za,zb)


      end function

      function zaj_virt_a6_slc_fcc_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_fcc_pm

      include 'heldefs.f'
      complex(dp) :: Fcc

      zaj_virt_a6_slc_fcc_pm=Fcc(hqpqbmgpgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_slc_fcc_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_fcc_mp
      include 'heldefs.f'
      complex(dp) :: Fcc

      zaj_virt_a6_slc_fcc_mp=Fcc(hqpqbmgmgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_slc_fsc_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_fsc_pp
      include 'heldefs.f'
      complex(dp) :: Fsc

      zaj_virt_a6_slc_fsc_pp=Fsc(hqpqbmgpgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_slc_fsc_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_slc_fsc_pm
      include 'heldefs.f'
      complex(dp) :: Fsc

      zaj_virt_a6_slc_fsc_pm=Fsc(hqpqbmgpgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_slc_fsc_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_slc_fsc_mp
      include 'heldefs.f'
      complex(dp) :: Fsc

      zaj_virt_a6_slc_fsc_mp=Fsc(hqpqbmgmgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

c     primitive amplitudes for q;1 leading color piece
c     !!  note that the argument ordering is changed to match
c     !!  q,qb,g,a,lb,l for the leading color pieces
      function zaj_virt_a6_lc_pole_pm(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_tree_pm
      complex(dp) :: zaj_virt_a6_lc_pole_pm
      complex(dp) :: Vcc, Vsc

      complex(dp) :: xl12, xl23, xl56
      complex(dp) :: lnrat

      real(dp) :: s,t
      s(i1,i2)=real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)


      xl12=Lnrat(musq,real(-s(i1,i2),dp))
      xl23=Lnrat(musq,real(-s(i2,i3),dp))
      xl56=Lnrat(musq,real(-s(i5,i6),dp))

      zaj_virt_a6_lc_tree_pm =
     &    (((za(i3,i5)*zb(i1,i3)+za(i4,i5)*zb(i1,i4))*
     &    (-(za(i1,i3)*zb(i1,i6))-za(i2,i3)*zb(i2,i6)))/
     &    (s(i5,i6)*za(i1,i2)*za(i2,i3)*zb(i1,i4)*zb(i3,i4))+
     &    (za(i3,i5)*zb(i1,i2)*(-(za(i1,i4)*zb(i1,i6))-za(i2,i4)*zb(i2,i6)))
     &    /
     &    (s(i5,i6)*za(i1,i2)*zb(i1,i4)*t(i1,i2,i4))-
     &    (za(i3,i4)*zb(i1,i6)*(za(i3,i5)*zb(i2,i3)+za(i4,i5)*zb(i2,i4)))/
     &    (s(i5,i6)*za(i2,i3)*zb(i3,i4)*t(i2,i3,i4)))

      Vcc=-four
     &   -(+(epinv*epinv2+xl12*epinv+half*xl12**2)
     &     +(epinv*epinv2+xl23*epinv+half*xl23**2))-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      zaj_virt_a6_lc_pole_pm = zaj_virt_a6_lc_tree_pm*(Vcc+Vsc)

      end function

      ! this helper routine makes it simpler to use the same
      ! crossing code from the LO amplitudes
      function zaj_virt_a6_lc_pole_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_pole_mp

      zaj_virt_a6_lc_pole_mp =
     &    zaj_virt_a6_lc_pole_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_pole_pp(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_tree_pp
      complex(dp) :: zaj_virt_a6_lc_pole_pp
      complex(dp) :: Vcc, Vsc

      complex(dp) :: xl12, xl23, xl56
      complex(dp) :: lnrat

      real(dp) :: s,t
      s(i1,i2)=real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)


      xl12=Lnrat(musq,real(-s(i1,i2),dp))
      xl23=Lnrat(musq,real(-s(i2,i3),dp))
      xl56=Lnrat(musq,real(-s(i5,i6),dp))

      zaj_virt_a6_lc_tree_pp =
     &    (-za(i1,i3)*za(i3,i5)**2)/
     &    (za(i1,i2)*za(i1,i4)*za(i2,i3)*za(i3,i4)*za(i5,i6))

      Vcc=-four
     &   -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   -(epinv*epinv2+xl23*epinv+half*xl23**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      zaj_virt_a6_lc_pole_pp = zaj_virt_a6_lc_tree_pp*(Vcc+Vsc)

      end function

      function zaj_virt_a6_lc_fcc_pp(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_fcc_pp

      include 'heldefs.f'
      complex(dp) :: Fcc

      zaj_virt_a6_lc_fcc_pp=Fcc(hqpgpqbmgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_lc_fcc_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_lc_fcc_mp

      zaj_virt_a6_lc_fcc_mp = zaj_virt_a6_lc_fcc_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_fcc_pm(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_fcc_pm

      include 'heldefs.f'
      complex(dp) :: Fcc

      zaj_virt_a6_lc_fcc_pm=Fcc(hqpgpqbmgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_lc_fsc_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a6_lc_fsc_mp

      zaj_virt_a6_lc_fsc_mp = zaj_virt_a6_lc_fsc_pm(i2,i1,i3,i4,i6,i5,zb,za)

      end function

      function zaj_virt_a6_lc_fsc_pm(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_fsc_pm
      include 'heldefs.f'
      complex(dp) :: Fsc

      zaj_virt_a6_lc_fsc_pm=Fsc(hqpgpqbmgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a6_lc_fsc_pp(i1,i3,i2,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a6_lc_fsc_pp
      include 'heldefs.f'
      complex(dp) :: Fsc

      zaj_virt_a6_lc_fsc_pp=Fsc(hqpgpqbmgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvf_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a64v_fvf_pp
      include 'heldefs.f'
      complex(dp) :: Fvf

      zaj_virt_a64v_fvf_pp=Fvf(hqpqbmgpgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvf_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a64v_fvf_mp

      zaj_virt_a64v_fvf_mp = zaj_virt_a64v_fvf_pm(i1,i2,i4,i3,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvf_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a64v_fvf_pm
      include 'heldefs.f'
      complex(dp) :: Fvf

      zaj_virt_a64v_fvf_pm=Fvf(hqpqbmgpgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvs_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a64v_fvs_pp
      include 'heldefs.f'
      complex(dp) :: Fvs

      zaj_virt_a64v_fvs_pp=Fvs(hqpqbmgpgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvs_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a64v_fvs_mp

      zaj_virt_a64v_fvs_mp = zaj_virt_a64v_fvs_pm(i1,i2,i4,i3,i5,i6,za,zb)

      end function

      function zaj_virt_a64v_fvs_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      complex(dp) :: zaj_virt_a64v_fvs_pm
      include 'heldefs.f'
      complex(dp) :: Fvs

      zaj_virt_a64v_fvs_pm=Fvs(hqpqbmgpgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64ax_pp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a64ax_pp
      include 'heldefs.f'
      complex(dp) :: a64ax

      zaj_virt_a64ax_pp=a64ax(hqpqbmgpgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64ax_pm(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a64ax_pm
      include 'heldefs.f'
      complex(dp) :: a64ax

      zaj_virt_a64ax_pm=a64ax(hqpqbmgpgm,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a64ax_mp(i1,i2,i3,i4,i5,i6,za,zb)
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a64ax_mp
      include 'heldefs.f'
      complex(dp) :: a64ax

      zaj_virt_a64ax_mp=a64ax(hqpqbmgmgp,i1,i2,i3,i4,i5,i6,za,zb)

      end function

      function zaj_virt_a65ax_pp(i1,i2,i3,i4,i5,i6,za,zb)
      include 'masses.f'
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a65ax_pp
      complex(dp) :: L1
      real(dp) :: s,t

      s(i1,i2) = real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

c Only include the exch34 term from (11.4) since the other one corresponds
c to the diagram with the photon attached inside the loop
      zaj_virt_a65ax_pp=
     & 2._dp*za(i2,i5)*zb(i3,i6)
     & *(za(i2,i1)*zb(i1,i3)+za(i2,i4)*zb(i4,i3))
     &  /(za(i1,i4)*za(i2,i4)*s(i5,i6))
     &  *(-0.5_dp*L1(-t(i1,i2,i4),-s(i5,i6))/s(i5,i6))

      end function

      function zaj_virt_a65ax_pm(i1,i2,i3,i4,i5,i6,za,zb)
      include 'masses.f'
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a65ax_pm
      complex(dp) :: L1
      real(dp) :: s,t

      s(i1,i2) = real(za(i1,i2)*zb(i2,i1),kind=dp)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

c Only include the flip2 term from (11.4) since the other one corresponds
c to the diagram with the photon attached inside the loop
      zaj_virt_a65ax_pm=
     &  +2._dp*zb(i1,i3)*zb(i3,i6)*(zb(i1,i2)*za(i2,i5)+zb(i1,i4)*za(i4,i5))
     &  /(zb(i2,i4)*zb(i1,i4)*s(i6,i5))
     &  *(-0.5_dp*L1(-t(i2,i1,i4),-s(i6,i5))/s(i6,i5))

      end function

      function zaj_virt_a65ax_mp(i1,i2,i3,i4,i5,i6,za,zb)
      include 'masses.f'
      integer, intent(in) :: i1,i2,i3,i4,i5,i6
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp) :: zaj_virt_a65ax_mp
      complex(dp) :: L1
      real(dp) :: s,t

      s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)

c This is the regular term from (11.4) with 3<->4, c.f. (11.11);
c the other term corresponds to the diagram with the photon attached inside the loop
      zaj_virt_a65ax_mp=
     &   2._dp*za(i2,i3)*za(i3,i5)*(za(i2,i1)*zb(i1,i6)+za(i2,i4)*zb(i4,i6))
     &  /(za(i1,i4)*za(i2,i4)*s(i5,i6))
     &  *(-0.5_dp*L1(-t(i1,i2,i4),-s(i5,i6))/s(i5,i6))

      end function

      end module
