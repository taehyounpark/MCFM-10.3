!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function fc002(msq,taugs,qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: fc002
      real(dp):: qsq,msq,qsqhat,taugs,qgs,cgs2,rat,rat2
      complex(dp):: lnrat,I3me

      qsqhat=qsq-msq
      qgs=qsqhat-taugs
      cgs2=(qsqhat-taugs)**2-four*taugs*msq
      rat=qsqhat/qsq
      rat2=rat/qsq
      fc002=
     & -1.D0/12.D0*qgs*cgs2**(-1)*msq
     & +1.D0/6.D0*cgs2**(-1)*msq*taugs
     & +1.D0/12.D0*cgs2**(-1)*msq**2*taugs*qsq**(-1)
     & +1.D0/12.D0*cgs2**(-1)*msq**2
     & -1.D0/12.D0*cgs2**(-1)*msq**3*qsq**(-1)
     & -I3me(msq,taugs,qsq)*cgs2**(-2)*msq**2*taugs**3
     & +1.D0/2.D0*lnrat(-taugs,msq)*qgs*cgs2**(-2)*msq*taugs**2
     & -1.D0/12.D0*lnrat(-taugs,msq)*qgs*cgs2**(-1)*taugs
     & -1.D0/2.D0*lnrat(-qsqhat,msq)*qgs*cgs2**(-2)*msq*taugs**2*rat
     & +1.D0/12.D0*lnrat(-qsqhat,msq)*qgs*cgs2**(-1)*msq*rat
     & +1.D0/12.D0*lnrat(-qsqhat,msq)*qgs*cgs2**(-1)*taugs*rat

      fc002=fc002
     & -lnrat(-qsqhat,msq)*cgs2**(-2)*msq**2*taugs**2*rat
     & -1.D0/12.D0*lnrat(-qsqhat,msq)*cgs2**(-1)*msq**2*taugs*rat2
     & -1.D0/12.D0*lnrat(-qsqhat,msq)*cgs2**(-1)*msq**2*rat
     & +1.D0/12.D0*lnrat(-qsqhat,msq)*cgs2**(-1)*msq**3*rat2
      return
      end

      function fc00(msq,taugs,qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp):: fc00

      real(dp):: qsq,msq,qsqhat,taugs,qgs,cgs2
      complex(dp):: lnrat,I3me

      qsqhat=qsq-msq
      qgs=qsqhat-taugs
      cgs2=(qsqhat-taugs)**2-four*taugs*msq

      fc00 =
     &  - 1.D0/2.D0*I3me(msq,taugs,qsq)*cgs2**(-1)*msq*taugs**2
     &  + 1.D0/4.D0*lnrat( - taugs,msq)*qgs*cgs2**(-1)*taugs
     &  - 1.D0/4.D0*lnrat( - qsqhat,msq)*qgs*cgs2**(-1)*taugs*qsq**(-1)
     & *qsqhat
     &  - 1.D0/2.D0*lnrat( - qsqhat,msq)*cgs2**(-1)*msq*taugs*qsq**(-1)
     & *qsqhat

      return
      end
