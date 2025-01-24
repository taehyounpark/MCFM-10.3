!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
! expansions of S and a up to 1-loop can be found in Eq. (B3) of 1109.6027
! expansions up to 2-loop can be found in Eqs. (93) and (94) of hep-ph/0607228 (Becher, Neubert, Pecjak)
! expansion up to 3-loop obtain from direct solution of Eq. (39) using Mathematica

      r = alphasMu/alphasMuH
      lnr=Log(r)
      aGamma0 = (Gamma0*lnr)/(2.*beta0)
      aGamma1 = (alphasMuH*(-(beta1*Gamma0) + beta0*Gamma1)*(r - 1))/
     &   (8.*beta0**2*Pi)
      aGamma2 = (alphasMuH**2*(beta1**2*Gamma0 - beta0*beta1*Gamma1 + 
     &       beta0*(-(beta2*Gamma0) + beta0*Gamma2))*(r - 1)*(1 + r))/
     &   (64.*beta0**3*Pi**2)
      aGamma3 = -(alphasMuH**3*(beta1**3*Gamma0 - beta0*beta1**2*Gamma1+ 
     &        beta0*beta1*(-2*beta2*Gamma0 + beta0*Gamma2) + 
     &        beta0**2*(beta3*Gamma0 + beta2*Gamma1 - beta0*Gamma3))*
     &      (-1 + r**3))/(384.*beta0**4*Pi**3)

      agammai0 = (gammai0*lnr)/(2.*beta0)
      agammai1 = (alphasMuH*(-beta1*gammai0 + beta0*gammai1)*(r - 1))/
     &   (8.*beta0**2*Pi)
      agammai2 = (alphasMuH**2*(beta1**2*gammai0 - beta0*beta1*gammai1 + 
     &       beta0*(-(beta2*gammai0) + beta0*gammai2))*(r - 1)*(1 + r))/
     &   (64.*beta0**3*Pi**2)
      agammai3 =
     & -(alphasMuH**3*(beta1**3*gammai0 - beta0*beta1**2*gammai1 + 
     &        beta0*beta1*(-2*beta2*gammai0 + beta0*gammai2) + 
     &        beta0**2*(beta3*gammai0 + beta2*gammai1 - beta0*gammai3))*
     &      (-1 + r**3))/(384.*beta0**4*Pi**3)

       
       S0 = (2*(beta1*Gamma0 - beta0*Gamma1)*(r - 1) + 
     &     lnr*(-2*beta1*Gamma0 + 2*beta0*Gamma1 + beta1*Gamma0*lnr) + 
     &     (8*beta0*Gamma0*Pi*(-1 + r - r*lnr))/(alphasMuH*r))/
     &   (8.*beta0**3)
       S1 = -(alphasMuH*((r - 1)*(-(beta0*beta1*Gamma1*(-3 + r)) + 
     &           beta1**2*Gamma0*(r - 1) + 
     &           beta0*(beta0*Gamma2*(r - 1) - beta2*Gamma0*(1 + r))) + 
     &        2*(beta0*beta2*Gamma0 + beta1**2*Gamma0*(r - 1) - 
     &           beta0*beta1*Gamma1*r)*lnr))/(32.*beta0**4*Pi)
       S2 = (alphasMuH**2*((r - 1)*
     &        (-4*beta1**3*Gamma0 + beta0*beta1**2*Gamma1*(5 + r) + 
     &          beta0*beta1*(-3*beta0*Gamma2*(1 + r) + 
     &             beta2*Gamma0*(5 + r)) + 
     &          beta0**2*(-(beta3*Gamma0*(1 + r)) + 
     &             2*(-2*beta2*Gamma1 + beta0*Gamma3*(1 + r)))) + 
     &      2*(-((beta1**3-2*beta0*beta1*beta2+beta0**2*beta3)*Gamma0) + 
     &          beta1*(beta1**2*Gamma0 - beta0*beta1*Gamma1 + 
     &             beta0*(-(beta2*Gamma0) + beta0*Gamma2))*r**2)*lnr))/
     &   (256.*beta0**5*Pi**2) 
       S3 = (alphasMuH**3*(-3*beta1**4*Gamma0*(r-1)**2*(7+r*(8 + 3*r))+ 
     &       beta0*beta1**3*Gamma1*(r - 1)*(-25 + r*(-7 + r*(11+9*r))) - 
     &       beta0**2*beta1*(r - 1)*
     &        (2*beta3*Gamma0*(-7 + 9*r)*(1 + r + r**2) + 
     &          2*beta2*Gamma1*(-20 + r*(1 + r)*(-2 + 9*r)) + 
     &          beta0*Gamma3*(7 + r*(7 + (7 - 9*r)*r))) + 
     &       beta0*beta1**2*(-1 + r)*
     &        (beta0*Gamma2*(13 + r*(13 - r*(5 + 9*r))) + 
     &          beta2*Gamma0*(-41 + r*(-5 + r*(31 + 27*r)))) + 
     &       beta0**2*(9*beta0*beta2*Gamma2*(-1 + r**2)**2 - 
     &          beta2**2*Gamma0*(5 + r**2*(-18 + r*(4 + 9*r))) + 
     &          beta0*(beta4*Gamma0*(-1 + r**3*(-8 + 9*r)) - 
     &             3*(beta3*Gamma1*(-5 + 6*r + 2*r**3 - 3*r**4) + 
     &                beta0*Gamma4*(1 + r**3*(-4 + 3*r))))) - 
     &       12*(-((beta1**4 - 3*beta0*beta1**2*beta2 + 
     &               2*beta0**2*beta1*beta3 + 
     &               beta0**2*(beta2**2 - beta0*beta4))*Gamma0) + 
     &          beta1*(beta1**3*Gamma0 - beta0*beta1**2*Gamma1 + 
     &             beta0*beta1*(-2*beta2*Gamma0 + beta0*Gamma2) + 
     &       beta0**2*(beta3*Gamma0 + beta2*Gamma1 - beta0*Gamma3))*r**3
     &          )*lnr))/(9216.*beta0**6*Pi**3)

      if (order < 3) then
        agammai3=0._dp
        aGamma3=0._dp
        S3=0._dp
       endif
      if (order < 2) then
        agammai2=0._dp
        aGamma2=0._dp
        S2=0._dp
       endif
      if (order < 1) then
        agammai1=0._dp
        aGamma1=0._dp
        S1=0._dp
       endif

! expansions
      S=S0+S1+S2+S3
      aGamma=aGamma0+aGamma1+aGamma2+aGamma3
      agammai=agammai0+agammai1+agammai2+agammai3

      
