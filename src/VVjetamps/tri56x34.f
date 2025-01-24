!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tri56x34(i1,i2,i3,i4,i5,i6,i7,za,zb,tot,tot_sl)
      implicit none
      include'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i,i1,i2,i3,i4,i5,i6,i7,j,k
      real(dp):: p34Dp56,gam(2),S1,S2,
     & al01(2),al02(2),al11(2),al12(2),al21(2),al22(2),
     & a01,a02,a11,a12,a21,D(2)
      complex(dp):: res,tot,res_sl,tot_sl,
     & xab1(7,7),xab2(7,7),iza,zab2
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      iza(i1,i2)=1._dp/za(i1,i2)
      S1=s(i3,i4)
      S2=s(i5,i6)
      p34Dp56=0.5_dp*(s(i5,i3)+s(i5,i4)+s(i6,i3)+s(i6,i4))
      call threemtrisoln(S1,S2,p34Dp56,
     & gam,al01,al02,al11,al12,al21,al22,D)

      tot=czip
      tot_sl=czip
      do i=1,2
      a01=al01(i)
      a02=al02(i)
      a11=al11(i)
      a12=al12(i)
      a21=al21(i)

      do j=1,7
      do k=1,7
      xab1(j,k)=
     & (+zab2(j,i3,i4,k)+S1/gam(i)*zab2(j,i5,i6,k))*gam(i)**2/D(i)
      xab2(j,k)=
     & (-zab2(j,i5,i6,k)-S2/gam(i)*zab2(j,i3,i4,k))*gam(i)**2/D(i)
      enddo
      enddo

      res =  + iza(i1,i7)*iza(i2,i7) * (
     &     - za(i1,i2)*xab1(i1,i4)*xab1(i3,i6)*xab1(i5,i7)/xab1(i2,i7)
     &    *a02*a12
     &     + za(i1,i2)*xab1(i1,i4)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)
     &    **2*xab2(i2,i6)*a11*a12
     &     - za(i1,i2)*xab1(i1,i4)*xab1(i3,i7)/xab1(i2,i7)*xab2(i5,i6)
     &    *a12*a21
     &     - za(i1,i2)*xab1(i1,i4)*xab1(i5,i7)/xab1(i2,i7)*xab2(i3,i6)
     &    *a12*a01
     &     + za(i1,i2)*xab1(i1,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)
     &    **2*xab2(i2,i4)*a02*a11
     &     - za(i1,i2)*xab1(i1,i6)*xab1(i3,i7)/xab1(i2,i7)*xab2(i5,i4)
     &    *a02*a21
     &     - za(i1,i2)*xab1(i1,i6)*xab1(i5,i7)/xab1(i2,i7)*xab2(i3,i4)
     &    *a02*a01
     &     )
      res = res + iza(i1,i7)*iza(i2,i7) * (
     &     - za(i1,i2)*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)
     &    **3*xab2(i2,i4)*xab2(i2,i6)*a11**2
     &     + za(i1,i2)*xab1(i1,i7)*xab1(i3,i7)/xab1(i2,i7)**2*xab2(i2,
     &    i4)*xab2(i5,i6)*a11*a21
     &     + za(i1,i2)*xab1(i1,i7)*xab1(i5,i7)/xab1(i2,i7)**2*xab2(i2,
     &    i4)*xab2(i3,i6)*a11*a01
     &     - za(i1,i2)*xab1(i1,i7)/xab1(i2,i7)*xab2(i3,i4)*xab2(i5,i6)
     &    *a21*a01
     &     + za(i1,i2)*xab1(i3,i4)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,i6)
     &    *a12*a21
     &     - 2*za(i1,i2)*xab1(i3,i4)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,
     &    i6)*a11*a12
     &     + za(i1,i2)*xab1(i3,i6)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,i4)
     &    *a02*a21
     &     )
      res = res + iza(i1,i7)*iza(i2,i7) * (
     &     - 2*za(i1,i2)*xab1(i3,i6)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,
     &    i4)*a02*a11
     &     - za(i1,i2)*xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i2,i7
     &    )*xab2(i1,i4)*xab2(i1,i6)*a21**2
     &     + 2*za(i1,i2)*xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i2,
     &    i7)*xab2(i1,i4)*xab2(i1,i6)*a11*a21
     &     - za(i1,i2)*xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i2,i7
     &    )*xab2(i1,i4)*xab2(i1,i6)*a11**2
     &     - za(i1,i2)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**2*xab2(i1,
     &    i4)*xab2(i2,i6)*a11*a21
     &     + 2*za(i1,i2)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**2*xab2(
     &    i1,i4)*xab2(i2,i6)*a11**2
     &     + za(i1,i2)*xab1(i3,i7)/xab1(i2,i7)*xab2(i1,i4)*xab2(i5,i6)
     &    *a21**2
     &     )
      res = res + iza(i1,i7)*iza(i2,i7) * (
     &     - 2*za(i1,i2)*xab1(i3,i7)/xab1(i2,i7)*xab2(i1,i4)*xab2(i5,
     &    i6)*a11*a21
     &     + za(i1,i2)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,i4)*xab2(i3,i6)
     &    *a21*a01
     &     - 2*za(i1,i2)*xab1(i5,i7)/xab1(i2,i7)*xab2(i1,i4)*xab2(i3,
     &    i6)*a11*a01
     &     )

      res_sl =  + iza(i1,i7) * (
     &     - xab1(i1,i4)*xab1(i1,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,
     &    i7)/xab1(i7,i7)*a02*a12
     &     + xab1(i1,i4)*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,
     &    i7)**2/xab1(i7,i7)*xab2(i2,i6)*a11*a12
     &     + xab1(i1,i4)*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,
     &    i7)/xab1(i7,i7)**2*xab2(i7,i6)*a12*a21
     &     - xab1(i1,i4)*xab1(i1,i7)*xab1(i3,i7)/xab1(i2,i7)/xab1(i7,
     &    i7)*xab2(i5,i6)*a12*a21
     &     - xab1(i1,i4)*xab1(i1,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,
     &    i7)*xab2(i3,i6)*a12*a01
     &     - 2*xab1(i1,i4)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i1,i6)*a11*a12
     &     + xab1(i1,i6)*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,
     &    i7)**2/xab1(i7,i7)*xab2(i2,i4)*a02*a11
     &     )
      res_sl = res_sl + iza(i1,i7) * (
     &     + xab1(i1,i6)*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,
     &    i7)/xab1(i7,i7)**2*xab2(i7,i4)*a02*a21
     &     - xab1(i1,i6)*xab1(i1,i7)*xab1(i3,i7)/xab1(i2,i7)/xab1(i7,
     &    i7)*xab2(i5,i4)*a02*a21
     &     - xab1(i1,i6)*xab1(i1,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,
     &    i7)*xab2(i3,i4)*a02*a01
     &     - 2*xab1(i1,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i1,i4)*a02*a11
     &     - xab1(i1,i7)**2*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**3
     &    /xab1(i7,i7)*xab2(i2,i4)*xab2(i2,i6)*a11**2
     &     - xab1(i1,i7)**2*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**2/
     &    xab1(i7,i7)**2*xab2(i2,i4)*xab2(i7,i6)*a11*a21
     &     - xab1(i1,i7)**2*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i4)*xab2(i7,i6)*a21**2
     &     )
      res_sl = res_sl + iza(i1,i7) * (
     &     + xab1(i1,i7)**2*xab1(i3,i7)/xab1(i2,i7)**2/xab1(i7,i7)*
     &    xab2(i2,i4)*xab2(i5,i6)*a11*a21
     &     + xab1(i1,i7)**2*xab1(i3,i7)/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i5,i4)*xab2(i7,i6)*a21**2
     &     + xab1(i1,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2/xab1(i7,i7)*
     &    xab2(i2,i4)*xab2(i3,i6)*a11*a01
     &     + xab1(i1,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i3,i4)*xab2(i7,i6)*a21*a01
     &     - xab1(i1,i7)**2/xab1(i2,i7)/xab1(i7,i7)*xab2(i3,i4)*xab2(
     &    i5,i6)*a21*a01
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)*xab2(i1,i4)*xab2(i2,i6)*a11**2
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i1,i4)*xab2(i7,i6)*a11*a21
     &     )
      res_sl = res_sl + iza(i1,i7) * (
     &     - 2*xab1(i1,i7)*xab1(i3,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i1,i4)*xab2(i5,i6)*a11*a21
     &     - 2*xab1(i1,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i1,i4)*xab2(i3,i6)*a11*a01
     &     - xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(i1,
     &    i4)*xab2(i1,i6)*a11**2
     &     )
      res_sl = res_sl + iza(i2,i7) * (
     &     + xab1(i1,i4)*xab1(i3,i6)*xab1(i5,i7)/xab1(i7,i7)*a02*a12
     &     - xab1(i1,i4)*xab1(i3,i7)*xab1(i5,i7)/xab1(i7,i7)**2*xab2(
     &    i7,i6)*a11*a12
     &     + xab1(i1,i4)*xab1(i3,i7)/xab1(i7,i7)*xab2(i5,i6)*a12*a21
     &     + xab1(i1,i4)*xab1(i5,i7)/xab1(i7,i7)*xab2(i3,i6)*a12*a01
     &     - xab1(i1,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i7,i7)**2*xab2(
     &    i7,i4)*a02*a11
     &     + xab1(i1,i6)*xab1(i3,i7)/xab1(i7,i7)*xab2(i5,i4)*a02*a21
     &     + xab1(i1,i6)*xab1(i5,i7)/xab1(i7,i7)*xab2(i3,i4)*a02*a01
     &     + xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i7,i7)**3*xab2(
     &    i7,i4)*xab2(i7,i6)*a11**2
     &     - xab1(i1,i7)*xab1(i3,i7)/xab1(i7,i7)**2*xab2(i5,i4)*xab2(
     &    i7,i6)*a11*a21
     &     - xab1(i1,i7)*xab1(i5,i7)/xab1(i7,i7)**2*xab2(i3,i4)*xab2(
     &    i7,i6)*a11*a01
     &     )
      res_sl = res_sl + iza(i2,i7) * (
     &     + xab1(i1,i7)/xab1(i7,i7)*xab2(i3,i4)*xab2(i5,i6)*a21*a01
     &     - xab1(i3,i4)*xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i6)*a12*a21
     &     + 2*xab1(i3,i4)*xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i6)*a11*a12
     &     - xab1(i3,i6)*xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i4)*a02*a21
     &     + 2*xab1(i3,i6)*xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i4)*a02*a11
     &     + xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i7,i7)*xab2(i1,
     &    i4)*xab2(i1,i6)*a21**2
     &     - 2*xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i7,i7)*xab2(
     &    i1,i4)*xab2(i1,i6)*a11*a21
     &     + xab1(i3,i7)*xab1(i5,i7)/xab1(i1,i7)/xab1(i7,i7)*xab2(i1,
     &    i4)*xab2(i1,i6)*a11**2
     &     + xab1(i3,i7)*xab1(i5,i7)/xab1(i7,i7)**2*xab2(i1,i4)*xab2(
     &    i7,i6)*a11*a21
     &     - 2*xab1(i3,i7)*xab1(i5,i7)/xab1(i7,i7)**2*xab2(i1,i4)*
     &    xab2(i7,i6)*a11**2
     &     )
      res_sl = res_sl + iza(i2,i7) * (
     &     - xab1(i3,i7)/xab1(i7,i7)*xab2(i1,i4)*xab2(i5,i6)*a21**2
     &     + 2*xab1(i3,i7)/xab1(i7,i7)*xab2(i1,i4)*xab2(i5,i6)*a11*a21
     &     - xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i4)*xab2(i3,i6)*a21*a01
     &     + 2*xab1(i5,i7)/xab1(i7,i7)*xab2(i1,i4)*xab2(i3,i6)*a11*a01
     &     )

      tot=tot+res
      tot_sl=tot_sl+res_sl
      enddo
      tot=tot/(2d0*s(i3,i4)*s(i5,i6))
      tot_sl=-tot_sl/(2d0*s(i3,i4)*s(i5,i6))
      return
      end




