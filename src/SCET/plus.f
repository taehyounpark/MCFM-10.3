!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function plus(pz,fx,L0,p1,fx1,L01,z,jaco)
c     pz function evaluated as z
c     fx parton evaluated as xb/z
c     singular denominator i.e. 1/(1-z) or log(1-z)/(1-z)
c     p1 function evaluated at z=1
c     fx1 parton evaluated as xb
c     L01=minus integrand of singular denominator from 0 to xb
c     i.e. log(1-xb) or 1/2*log^2(1-xb)
c     value of z;
c     jaco=jacobian
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: plus
      real(dp), intent(in)::pz,p1,fx,fx1,z,L0,L01,jaco
      if(abs(one-z) < 1.e-5_dp) then
      plus=L01*p1*fx1
      else
      plus=(pz*fx/z-p1*fx1)*L0*jaco+L01*p1*fx1
      endif
      return
      end
