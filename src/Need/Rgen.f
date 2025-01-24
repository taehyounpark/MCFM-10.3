!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Rgen(pjet,i,p,j)
      implicit none
      include 'types.f'
      real(dp):: Rgen
c--- Calculate the angular separation between pjet(i) and p(j)
c--- This routine is a generalization of R.f: Rgen(p,i,p,j) == R(p,i,j)
      include 'mxpart.f'
      real(dp):: pjet(mxpart,4),p(mxpart,4),
     & r1,r2,dely,delphi,ei,ej,biti,bitj
      integer:: i,j
      real(dp), parameter:: tiny=1.e-9_dp

      ei=sqrt(pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
      ej=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)

      biti=pjet(i,3)/ei
      bitj=p(j,3)/ej
      if  ((abs(1._dp+biti) < tiny) .or. (abs(1._dp-biti) < tiny)
     & .or.(abs(1._dp+bitj) < tiny) .or. (abs(1._dp-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
        dely=100._dp
      else
        r1=(1._dp+biti)*(1._dp-bitj)/((1._dp+bitj)*(1._dp-biti))
        dely=0.5_dp*dlog(r1)
      endif

      r2= (pjet(i,1)*p(j,1)+pjet(i,2)*p(j,2))
     &     /sqrt((pjet(i,1)**2+pjet(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      Rgen=sqrt(dely**2+delphi**2)

      return
      end

