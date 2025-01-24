!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine tri12x34sl(j1,j2,j3,j4,j5,j6,j7,za,zb,tot,tot_sl)
      implicit none
      include'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i,j1,j2,j3,j4,j5,j6,j7,j,k
      real(dp):: p12Dp34,gam(2),S1,S2,
     & al01(2),al02(2),al11(2),al12(2),al21(2),al22(2),
     & a01,a02,a11,a12,a21,D(2)
      complex(dp):: tot,res_sl,tot_sl,
     & xab1(7,7),xab2(7,7),zab2
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      S1=s(j1,j2)
      S2=s(j3,j4)
      p12Dp34=0.5_dp*(s(j1,j3)+s(j1,j4)+s(j2,j3)+s(j2,j4))
      call threemtrisoln(S1,S2,p12Dp34,
     & gam,al01,al02,al11,al12,al21,al22,D)

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
     & (+zab2(j,j1,j2,k)+S1/gam(i)*zab2(j,j3,j4,k))*gam(i)**2/D(i)
      xab2(j,k)=
     & (-zab2(j,j3,j4,k)-S2/gam(i)*zab2(j,j1,j2,k))*gam(i)**2/D(i)
      enddo
      enddo
      res_sl=zb(j5,j6)*xab1(j5,j7)**2*xab1(j3,j7)*(
     & +xab1(j1,j2)*xab2(j7,j4)/xab2(j7,j2)*a12
     & -xab1(j1,j2)*xab2(j7,j4)/xab2(j7,j2)*a02
     & +xab1(j1,j4)*a02
     & -xab1(j1,j7)/xab1(j2,j7)*xab2(j2,j4)*a01
     & +xab1(j1,j7)/xab1(j3,j7)*xab2(j3,j4)*a21
     & +2*xab1(j1,j7)/xab1(j5,j7)*xab2(j5,j4)*a11
     & -xab1(j1,j7)/xab1(j7,j7)*xab2(j7,j4)*a21
     & -xab1(j1,j7)/xab1(j7,j7)*xab2(j7,j4)*a11
     & +xab2(j1,j4)*a01)/(xab1(j2,j7)*xab1(j7,j7)**2)

      tot_sl=tot_sl+res_sl
      enddo
      tot_sl=-tot_sl/(2d0*s(j3,j4)*s(j5,j6))
      tot=czip ! no leading color piece
      return
      end




