!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tri12x56sl(i1,i2,i3,i4,i5,i6,i7,za,zb,tot,tot_sl)
      implicit none
      include'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      real(dp):: p56Dp12
      integer i,i1,i2,i3,i4,i5,i6,i7,j,k
      real(dp):: gam(2),S1,S2,
     & al01(2),al02(2),al11(2),al12(2),al21(2),al22(2),
     & a01,a02,a11,a12,a21,a22,D(2)
      complex(dp):: tot,res_sl,tot_sl,xab1(7,7),xab2(7,7),zab2
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      S1=s(i5,i6)
      S2=s(i1,i2)
      p56Dp12=0.5_dp*(s(i5,i1)+s(i5,i2)+s(i6,i1)+s(i6,i2))
      call threemtrisoln(S1,S2,p56Dp12,
     & gam,al01,al02,al11,al12,al21,al22,D)

      tot_sl=czip
      do i=1,2
      a01=al01(i)
      a02=al02(i)
      a11=al11(i)
      a12=al12(i)
      a21=al21(i)
      a22=al22(i)
      do j=1,7
      do k=1,7
      xab1(j,k)=
     & (+zab2(j,i5,i6,k)+S1/gam(i)*zab2(j,i1,i2,k))*gam(i)**2/D(i)
      xab2(j,k)=
     & (-zab2(j,i1,i2,k)-S2/gam(i)*zab2(j,i5,i6,k))*gam(i)**2/D(i)
      enddo
      enddo

      res_sl =  + zb(i3,i4) * (
     &     - xab1(i1,i2)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i6)/xab2(i7,i2)*a22
     &     + xab1(i1,i2)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i6)/xab2(i7,i2)*a02
     &     + xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*a12
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)**2*xab2(i2,i6)*a21
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i6)*a21
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i6)*a11
     &     + xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i5,i6)*a01)
      res_sl = res_sl + zb(i3,i4) * (
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i3,i6)*a11
     &     + xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i1,i6)*a21)

      tot_sl=tot_sl+res_sl
      enddo
      tot_sl=-tot_sl/(2d0*s(i3,i4)*s(i5,i6))

      tot=czip ! no leading color piece

      return
      end

