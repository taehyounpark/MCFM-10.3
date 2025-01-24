!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tri56x17(i1,i2,i3,i4,i5,i6,i7,za,zb,tot,tot_sl)
      implicit none
      include'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i,i1,i2,i3,i4,i5,i6,i7,j,k
      real(dp):: p56Dp17,gam(2),S1,S2,S12,
     & al01(2),al02(2),al11(2),al12(2),al21(2),al22(2),
     & a01,a02,a11,a12,a21,a22,D(2),ga,gammS1,a22ma02
      complex(dp):: res,tot,res_sl,tot_sl,
     & xab1(7,7),xab2(7,7),zab2
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      S1=s(i5,i6)
      S2=s(i1,i7)
      p56Dp17=0.5_dp*(s(i5,i1)+s(i5,i7)+s(i6,i1)+s(i6,i7))
      call threemtrisoln(S1,S2,p56Dp17,
     & gam,al01,al02,al11,al12,al21,al22,D)

      tot=czip
      tot_sl=czip
      do i=1,2
      ga=gam(i)
      gammS1=ga-S1
      a01=al01(i)
      a02=al02(i)
      a11=al11(i)
      a12=al12(i)
      a21=al21(i)
      a22=al22(i)
      a22ma02=a22-a02
      S12=-2*(ga/D(i))**2
     & *((ga**2+S1*S2)*p56Dp17+ga*(s(i1,i7)*S1+s(i5,i6)*S2))
      do j=1,7
      do k=1,7
      xab1(j,k)=
     & (+zab2(j,i5,i6,k)+S1/ga*zab2(j,i1,i7,k))*ga**2/D(i)
      xab2(j,k)=
     & (-zab2(j,i1,i7,k)-S2/ga*zab2(j,i5,i6,k))*ga**2/D(i)
      enddo
      enddo

      res =  + zb(i3,i4)*gammS1**(-1)*a22ma02**(-1) * (
     &     + xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)*xab2(i2,i7)*a12*a21*ga*S12**(-1)
     &     + xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*a12*a22*ga*S12**(-1)
     &     - 2*xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)
     &    /xab1(i7,i7)*a02*a12*ga*S12**(-1)
     &     + xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i7)*a12*a21*ga*S12**(-1)
     &     - xab1(i1,i6)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i5,i7)*a12*a01*ga*S12**(-1)
     &     - 2*xab1(i1,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i3,i7)*a11*a12*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**3
     &    /xab1(i7,i7)*xab2(i2,i6)*xab2(i2,i7)*a21**2*ga*S12**(-1)
     &     )
      res = res + zb(i3,i4)*gammS1**(-1)*a22ma02**(-1) * (
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)**2*xab2(i2,i6)*xab2(i7,i7)*a21**2*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)*xab2(i2,i6)*a21*a22*ga*S12**(-1)
     &     + 2*xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)*xab2(i2,i6)*a02*a21*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i6)*xab2(i7,i7)*a21**2*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i6)*a21*a22*ga*S12**(-1)
     &     + 2*xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)
     &    /xab1(i7,i7)**2*xab2(i7,i6)*a02*a21*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i7,i6)/xab2(i7,i7)*a22**2*ga*S12**(-1)
     &     )
      res = res + zb(i3,i4)*gammS1**(-1)*a22ma02**(-1) * (
     &     + 2*xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)
     &    /xab1(i7,i7)*xab2(i7,i6)/xab2(i7,i7)*a02*a22*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i7,i6)/xab2(i7,i7)*a02**2*ga*S12**(-1)
     &     + xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)**2/xab1(i7,i7)*
     &    xab2(i2,i6)*xab2(i5,i7)*a21*a01*ga*S12**(-1)
     &     + xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i5,i6)*xab2(i7,i7)*a21*a01*ga*S12**(-1)
     &     + xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i5,i6)*a22*a01*ga*S12**(-1)
     &     - 2*xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)*
     &    xab2(i5,i6)*a02*a01*ga*S12**(-1)
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)*xab2(i2,i6)*xab2(i3,i7)*a11*a21*ga*S12**(-1)
     &     )
      res = res + zb(i3,i4)*gammS1**(-1)*a22ma02**(-1) * (
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i3,i6)*xab2(i7,i7)*a11*a21*ga*S12**(-1)
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i3,i6)*a11*a22*ga*S12**(-1)
     &     - 4*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)*xab2(i3,i6)*a02*a11*ga*S12**(-1)
     &     - 2*xab1(i1,i7)*xab1(i3,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i3,i6)*xab2(i5,i7)*a11*a01*ga*S12**(-1)
     &     - xab1(i1,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(i3,
     &    i6)*xab2(i3,i7)*a11**2*ga*S12**(-1)
     &     - xab1(i3,i6)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,
     &    i7)*xab2(i1,i7)*a12*a21*ga*S12**(-1)
     &     + xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2/xab1(i7,i7)*
     &    xab2(i1,i6)*xab2(i2,i7)*a21**2*ga*S12**(-1)
     &     )
      res = res + zb(i3,i4)*gammS1**(-1)*a22ma02**(-1) * (
     &     + xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i1,i6)*xab2(i7,i7)*a21**2*ga*S12**(-1)
     &     + xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i1,i6)*a21*a22*ga*S12**(-1)
     &     - 2*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*
     &    xab2(i1,i6)*a02*a21*ga*S12**(-1)
     &     - xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)*xab2(i1,i6)*xab2(
     &    i5,i7)*a21*a01*ga*S12**(-1)
     &     - 2*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)*xab2(
     &    i1,i6)*xab2(i3,i7)*a11*a21*ga*S12**(-1)
     &     )

      res_sl =  + zb(i3,i4)*gammS1**(-1) * (
     &     + xab1(i1,i6)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*a12*ga
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)**2
     &    /xab1(i7,i7)**2*xab2(i2,i6)*a21*ga
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i6)*a01*ga
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**3*xab2(i7,i6)*a21*ga
     &     - xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i6)/xab2(i7,i7)*a22*ga
     &     + xab1(i1,i7)*xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i7,i6)/xab2(i7,i7)*a02*ga
     &     + xab1(i1,i7)*xab1(i3,i7)**2/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i5,i6)*a01*ga
     &     )
      res_sl = res_sl + zb(i3,i4)*gammS1**(-1) * (
     &     + 2*xab1(i1,i7)*xab1(i3,i7)*xab1(i5,i7)/xab1(i2,i7)/xab1(
     &    i7,i7)**2*xab2(i3,i6)*a11*ga
     &     + xab1(i3,i7)**2*xab1(i5,i7)/xab1(i2,i7)/xab1(i7,i7)**2*
     &    xab2(i1,i6)*a21*ga
     &     )
      tot=tot+res
      tot_sl=tot_sl+res_sl
      enddo
      tot=-tot/(2d0*s(i3,i4)*s(i5,i6))
      tot_sl=-tot_sl/(2d0*s(i3,i4)*s(i5,i6))
      tot_sl=tot_sl+tot
      return
      end




