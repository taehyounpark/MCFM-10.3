!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function amp_qqagg_ql(i1,qh,i2,h2,i3,h3,i4,h4,i5,lh,j6,j7,za,zb)
      !use zajj_treeamps_m, only: zajj_tree_qqgg_fsr_ppp, zajj_tree_qqgg_anomZZ_ppp
      !use zaj_treeamps_m, only: zaj_tree_fsr_pp, zaj_tree_anomZZ_pp
      !use mcfm_debug
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      complex(dp):: amp_qqagg_ql
      integer, intent(in) :: i1,i2,i3,i4,i5,j6,j7
      integer, intent(in) :: qh,h2,h3,h4,lh
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

      integer :: g1,g2,g3,g4,lg,i6,i7
      complex(dp) :: s,t
      complex(dp) :: A7h1,A7h2,A7h3,A7h4,A7h5,A7h6,A7h7,A7h8
      !real(dp) :: s267
      !complex(dp) :: ampSeven, ampSix

      s(i1,i2) = za(i1,i2)*zb(i2,i1)
      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


c---A(qh,h2,h3,h4,lh)
c---h=1 LH
c---h=2 RH
      g1=qh
      g2=h2
      g3=h3
      g4=h4
      lg=lh
      i6=j6
      i7=j7

      if (g1*lg==2)then
      i7=j6
      i6=j7
      lg=3-lg
      endif

      if (g1==2) then

      elseif (g1==1) then
      g1=2
      g2=3-g2
      g3=3-g3
      g4=3-g4
      lg=3-lg

      endif

      !s267 = s(2,6) + s(6,7) + s(7,2)

c-----A7h1
c-----A(2,2,2,2,2)
      if (
     & (g1==2).and.(g2==2).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,2,2,2)
      A7h1=czip
      A7h1=za(i5,i6)**2
     &/(za(i1,i3)*za(i3,i4)*za(i4,i5)*za(i2,i6)*za(i2,i7))
      amp_qqagg_ql = A7h1

c     ampSeven = zajj_tree_qqgg_fsr_ppp(i5,i1,i6,i7,i2,i3,i4,za,zb)
c     ampSix = zaj_tree_fsr_pp(i1,i5,i4,i2,i6,i7,za,zb)*s267
c     write (*,*) "ppp", ampSeven/(A7h1*s267)
c     write (*,*) "i1,i2,i3,i4,i5,i6,i7", i1,i2,i3,i4,i5,i6,i7
c       call check_softfact(ampSix,ampSeven,za,zb, i1, i3, i4, 2)

c     write (*,*) "ANOM: "
c     ampSeven = zajj_tree_qqgg_anomZZ_ppp(i5,i1,i6,i7,i2,i3,i4,za,zb,(1d0,0d0),(1d0,0d0))
c     ampSix = zaj_tree_anomZZ_pp(i1,i5,i4,i2,i6,i7,za,zb,.false.)*s267
c     call check_softfact(ampSix,ampSeven,za,zb, i1, i3, i4, 2)
c     write (*,*) ""

c-----A7h2
c-----A(2,2,2,1,2) =
      elseif
     &((g1==2).and.(g2==2).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,2,1,2)
      A7h2=czip
      A7h2=1/(s(i3,i4)*t(i2,i6,i7)*za(i2,i6)*za(i2,i7))*(
     & za(i4,i1)*zb(i1,i3)*( za(i4,i1)*zb(i1,i2)*za(i2,i6)
     & +za(i4,i1)*zb(i1,i7)*za(i7,i6)+za(i4,i3)*zb(i3,i2)*za(i2,i6)
     & +za(i4,i3)*zb(i3,i7)*za(i7,i6) )*za(i5,i6)/za(i1,i3)/t(i1,i3,i4)
     &-za(i4,i5)*zb(i5,i3)*( za(i6,i4)*zb(i4,i3)+za(i6,i5)*zb(i5,i3) )*
     & ( za(i6,i2)*zb(i2,i1)+za(i6,i7)*zb(i7,i1) )/zb(i4,i5)/t(i3,i4,i5)
     &-( za(i4,i1)*zb(i1,i2)*za(i2,i6)+za(i4,i1)*zb(i1,i7)*za(i7,i6)
     &  +za(i4,i3)*zb(i3,i2)*za(i2,i6)+za(i4,i3)*zb(i3,i7)*za(i7,i6) )*
     & ( za(i6,i4)*zb(i4,i3)+za(i6,i5)*zb(i5,i3) )/za(i1,i3)/zb(i4,i5)
     &)
      amp_qqagg_ql = A7h2
      !write (*,*) "ppm", zajj_tree_qqgg_fsr_ppm(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h2*s267)
c-----A7h3
c-----A(2,2,1,2,2) =
      elseif
     &((g1==2).and.(g2==2).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,2,1,2,2)
      A7h3=czip
      A7h3=1/(s(i3,i4)*t(i2,i6,i7)*za(i2,i6)*za(i2,i7))*(
     &-zb(i1,i4)**2*( za(i3,i1)*zb(i1,i2)*za(i2,i6)
     & +za(i3,i1)*zb(i1,i7)*za(i7,i6)+za(i3,i4)*zb(i4,i2)*za(i2,i6)
     & +za(i3,i4)*zb(i4,i7)*za(i7,i6) )*za(i5,i6)/zb(i1,i3)/t(i1,i3,i4)
     &+za(i3,i5)**2*( za(i6,i3)*zb(i3,i4)+za(i6,i5)*zb(i5,i4) )*
     & ( za(i6,i2)*zb(i2,i1)+za(i6,i7)*zb(i7,i1) )/za(i4,i5)/t(i3,i4,i5)
     &+zb(i1,i4)*za(i3,i5)*za(i5,i6)*( za(i6,i2)*zb(i2,i1)
     & +za(i6,i7)*zb(i7,i1) )/zb(i1,i3)/za(i4,i5) )
      amp_qqagg_ql = A7h3
      !write (*,*) "pmp", zajj_tree_qqgg_fsr_pmp(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h3*s267)
c-----A7h4
c-----A(2,1,2,2,2)
      elseif
     &((g1==2).and.(g2==1).and.(g3==2).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,2,2,2)
      A7h4=czip
      A7h4=(za(i5,i2)*zb(i2,i7)+za(i5,i6)*zb(i6,i7))**2
     &/(t(i2,i6,i7)*za(i1,i3)*za(i3,i4)*za(i4,i5)*zb(i2,i6)*zb(i2,i7))
      amp_qqagg_ql = A7h4
      !write (*,*) "mpp", zajj_tree_qqgg_fsr_mpp(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h4*s267)
c-----A7h5
c-----A(2,2,1,1,2) =
      elseif
     &((g1==2).and.(g2==2).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,2,1,1,2)
      A7h5=czip
      A7h5=-(za(i6,i2)*zb(i2,i1)+za(i6,i7)*zb(i7,i1))**2
     &/(t(i2,i6,i7)*zb(i1,i3)*zb(i3,i4)*zb(i4,i5)*za(i2,i6)*za(i2,i7))
      amp_qqagg_ql = A7h5
      !write (*,*) "pmm", zajj_tree_qqgg_fsr_pmm(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h5*s267)
c-----A7h6
c-----A(2,1,2,1,2) =
      elseif
     &((g1==2).and.(g2==1).and.(g3==2).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,2,1,2)
      A7h6=czip
      A7h6=1/(s(i3,i4)*t(i2,i6,i7)*zb(i2,i6)*zb(i2,i7))*(
     & za(i4,i1)*zb(i1,i3)*( za(i4,i1)*zb(i1,i7)+za(i4,i3)*zb(i3,i7) )*
     & ( za(i5,i2)*zb(i2,i7)+za(i5,i6)*zb(i6,i7) )/za(i1,i3)/t(i1,i3,i4)
     &-za(i4,i5)*zb(i5,i3)*( zb(i3,i4)*za(i4,i2)*zb(i2,i7)
     & + zb(i3,i4)*za(i4,i6)*zb(i6,i7) + zb(i3,i5)*za(i5,i2)*zb(i2,i7)
     & + zb(i3,i5)*za(i5,i6)*zb(i6,i7) )*zb(i1,i7)/zb(i4,i5)/t(i3,i4,i5)
     &-( zb(i3,i4)*za(i4,i2)*zb(i2,i7) + zb(i3,i4)*za(i4,i6)*zb(i6,i7)
     &  +zb(i3,i5)*za(i5,i2)*zb(i2,i7) + zb(i3,i5)*za(i5,i6)*zb(i6,i7))*
     & ( za(i4,i1)*zb(i1,i7)+za(i4,i3)*zb(i3,i7) )/za(i1,i3)/zb(i4,i5))
      amp_qqagg_ql = A7h6
      !write (*,*) "mpm", zajj_tree_qqgg_fsr_mpm(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h6*s267)
c-----A7h7
c-----A(2,1,1,2,2) =
      elseif
     &((g1==2).and.(g2==1).and.(g3==1).and.(g4==2).and.(lg==2)
     & ) then
c-----A(2,1,1,2,2)
      A7h7=czip
      A7h7=1/(s(i3,i4)*t(i2,i6,i7)*zb(i2,i6)*zb(i2,i7))*(
     &-zb(i1,i4)**2*( za(i3,i1)*zb(i1,i7)+za(i3,i4)*zb(i4,i7) )*
     & ( za(i5,i2)*zb(i2,i7)+za(i5,i6)*zb(i6,i7) )/zb(i1,i3)/t(i1,i3,i4)
     &+za(i3,i5)**2*zb(i1,i7)*( zb(i4,i3)*za(i3,i2)*zb(i2,i7)
     & +zb(i4,i3)*za(i3,i6)*zb(i6,i7)+zb(i4,i5)*za(i5,i2)*zb(i2,i7)
     & +zb(i4,i5)*za(i5,i6)*zb(i6,i7) )/za(i4,i5)/t(i3,i4,i5)
     &-zb(i1,i4)*za(i3,i5)*( za(i5,i2)*zb(i2,i7)+za(i5,i6)*zb(i6,i7) )*
     & zb(i7,i1)/zb(i1,i3)/za(i4,i5) )
      amp_qqagg_ql = A7h7
      !write (*,*) "mmp", zajj_tree_qqgg_fsr_mmp(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h7*s267)
c-----A7h8
c-----A(2,1,1,1,2)
      elseif
     &((g1==2).and.(g2==1).and.(g3==1).and.(g4==1).and.(lg==2)
     & ) then
c-----A(2,1,1,1,2)
      A7h8=czip
      A7h8=-zb(i1,i7)**2
     &/(zb(i1,i3)*zb(i3,i4)*zb(i4,i5)*zb(i2,i6)*zb(i2,i7))
      amp_qqagg_ql = A7h8
      !write (*,*) "mmm", zajj_tree_qqgg_fsr_mmm(i5,i1,i6,i7,i2,i3,i4,za,zb)/(A7h8*s267)
c----
      else
      amp_qqagg_ql = czip
      endif
c      write(*,*) qh,h2,h3,h4,lh
c      write(*,*) amp_qqagg_ql
c-----done
      return
      end

