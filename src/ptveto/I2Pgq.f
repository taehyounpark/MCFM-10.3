!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2P(gq,rglr)= + 224.D0/27.D0*CF*x**(-1)*nfl + 160.D0/9.D0*CF*
     &    x**(-1)*lnmuonpt*nfl + 40.D0/9.D0*CF*x**(-1)*lnomx*nfl + 32.D0
     &    /3.D0*CF*x**(-1)*lnomx*lnmuonpt*nfl + 4.D0/3.D0*CF*x**(-1)*
     &    lnomx**2*nfl - 172.D0/27.D0*CF*nfl - 56.D0/9.D0*CF*lnmuonpt*
     &    nfl - 32.D0/9.D0*CF*lnomx*nfl - 16.D0/3.D0*CF*lnomx*lnmuonpt*
     &    nfl - 2.D0/3.D0*CF*lnomx**2*nfl - 52.D0/27.D0*CF*omx*nfl - 
     &    104.D0/9.D0*CF*omx*lnmuonpt*nfl - 8.D0/9.D0*CF*omx*lnomx*nfl
     &     - 16.D0/3.D0*CF*omx*lnomx*lnmuonpt*nfl - 2.D0/3.D0*CF*omx*
     &    lnomx**2*nfl - 3160.D0/27.D0*CF*CA*x**(-1) - 8*CF*CA*x**(-1)*
     &    Li3mx - 8*CF*CA*x**(-1)*lnmuonpt - 248.D0/3.D0*CF*CA*x**(-1)*
     &    lnmuonpt**2 - 8.D0/3.D0*CF*CA*x**(-1)*lnopx**3 - 32*CF*CA*
     &    x**(-1)*Li2mx*lnmuonpt - 304.D0/9.D0*CF*CA*x**(-1)*lnomx - 
     &    176.D0/3.D0*CF*CA*x**(-1)*lnomx*lnmuonpt + 32*CF*CA*x**(-1)*
     &    lnomx*lnmuonpt**2 - 22.D0/3.D0*CF*CA*x**(-1)*lnomx**2 - 16*CF
     &    *CA*x**(-1)*lnomx**2*lnmuonpt - 4.D0/3.D0*CF*CA*x**(-1)*
     &    lnomx**3
      I2P(gq,rglr) = I2P(gq,rglr) - 32*CF*CA*x**(-1)*lnx*lnmuonpt**2 - 
     &    32*CF*CA*x**(-1)*lnx*lnopx*lnmuonpt + 8*CF*CA*x**(-1)*lnx*
     &    lnopx**2 + 8*CF*CA*x**(-1)*lnx*Li2mx - 88.D0/3.D0*CF*CA*
     &    x**(-1)*lnx*lnomx + 32*CF*CA*x**(-1)*lnx*lnomx*lnmuonpt + 4*
     &    CF*CA*x**(-1)*lnx*lnomx**2 - 4*CF*CA*x**(-1)*lnx**2*lnopx + 4
     &    *CF*CA*x**(-1)*lnx**2*lnomx + 16*CF*CA*x**(-1)*Li3xonopx - 40
     &    *CF*CA*x**(-1)*Li3x - 88.D0/3.D0*CF*CA*x**(-1)*Li2x + 24*CF*
     &    CA*x**(-1)*Li2x*lnx + 100*CF*CA - 12*CF*CA*Li3mx - 668.D0/9.D0
     &    *CF*CA*lnmuonpt + 248.D0/3.D0*CF*CA*lnmuonpt**2 - 4*CF*CA*
     &    lnopx**3 + 4*CF*CA*Li2mx - 48*CF*CA*Li2mx*lnmuonpt + 218.D0/9.
     &    D0*CF*CA*lnomx + 40.D0/3.D0*CF*CA*lnomx*lnmuonpt - 16*CF*CA*
     &    lnomx*lnmuonpt**2 + 17.D0/3.D0*CF*CA*lnomx**2 + 8*CF*CA*
     &    lnomx**2*lnmuonpt + 2.D0/3.D0*CF*CA*lnomx**3 - 662.D0/9.D0*CF
     &    *CA*lnx + 472.D0/3.D0*CF*CA*lnx*lnmuonpt - 64*CF*CA*lnx*
     &    lnmuonpt**2 + 4*CF*CA*lnx*lnopx - 48*CF*CA*lnx*lnopx*lnmuonpt
     &     + 12*CF*CA*lnx*lnopx**2
      I2P(gq,rglr) = I2P(gq,rglr) + 12*CF*CA*lnx*Li2mx + 76.D0/3.D0*CF*
     &    CA*lnx*lnomx - 16*CF*CA*lnx*lnomx*lnmuonpt - 2*CF*CA*lnx*
     &    lnomx**2 + 53.D0/3.D0*CF*CA*lnx**2 - 24*CF*CA*lnx**2*lnmuonpt
     &     - 6*CF*CA*lnx**2*lnopx - 2*CF*CA*lnx**2*lnomx - 2*CF*CA*
     &    lnx**3 + 24*CF*CA*Li3xonopx - 12*CF*CA*Li3x + 88.D0/3.D0*CF*
     &    CA*Li2x + 4*CF*CA*Li2x*lnx - 16.D0/3.D0*CF*CA*omx + 4*CF*CA*
     &    omx*Li3mx + 868.D0/9.D0*CF*CA*omx*lnmuonpt - 88.D0/3.D0*CF*CA
     &    *omx*lnmuonpt**2 + 4.D0/3.D0*CF*CA*omx*lnopx**3 - 4*CF*CA*omx
     &    *Li2mx + 16*CF*CA*omx*Li2mx*lnmuonpt + 86.D0/9.D0*CF*CA*omx*
     &    lnomx + 136.D0/3.D0*CF*CA*omx*lnomx*lnmuonpt - 16*CF*CA*omx*
     &    lnomx*lnmuonpt**2 + 5.D0/3.D0*CF*CA*omx*lnomx**2 + 8*CF*CA*
     &    omx*lnomx**2*lnmuonpt + 2.D0/3.D0*CF*CA*omx*lnomx**3 + 340.D0/
     &    9.D0*CF*CA*omx*lnx - 248.D0/3.D0*CF*CA*omx*lnx*lnmuonpt + 32*
     &    CF*CA*omx*lnx*lnmuonpt**2 - 4*CF*CA*omx*lnx*lnopx + 16*CF*CA*
     &    omx*lnx*lnopx*lnmuonpt - 4*CF*CA*omx*lnx*lnopx**2 - 4*CF*CA*
     &    omx*lnx*Li2mx
      I2P(gq,rglr) = I2P(gq,rglr) + 4.D0/3.D0*CF*CA*omx*lnx*lnomx - 16*
     &    CF*CA*omx*lnx*lnomx*lnmuonpt - 2*CF*CA*omx*lnx*lnomx**2 - 25.D
     &    0/3.D0*CF*CA*omx*lnx**2 + 8*CF*CA*omx*lnx**2*lnmuonpt + 2*CF*
     &    CA*omx*lnx**2*lnopx - 2*CF*CA*omx*lnx**2*lnomx + 2.D0/3.D0*CF
     &    *CA*omx*lnx**3 - 8*CF*CA*omx*Li3xonopx + 20*CF*CA*omx*Li3x - 
     &    8.D0/3.D0*CF*CA*omx*Li2x - 12*CF*CA*omx*Li2x*lnx + 608.D0/27.D
     &    0*CF*CA*omx**2 - 352.D0/9.D0*CF*CA*omx**2*lnmuonpt + 32.D0/3.D
     &    0*CF*CA*omx**2*lnmuonpt**2 - 176.D0/9.D0*CF*CA*omx**2*lnx + 
     &    64.D0/3.D0*CF*CA*omx**2*lnx*lnmuonpt + 16.D0/3.D0*CF*CA*
     &    omx**2*lnx*lnomx + 8.D0/3.D0*CF*CA*omx**2*lnx**2 + 16.D0/3.D0
     &    *CF*CA*omx**2*Li2x + 44.D0/9.D0*CF*CA*pisq*x**(-1) + 4.D0/3.D0
     &    *CF*CA*pisq*x**(-1)*lnmuonpt + 4.D0/3.D0*CF*CA*pisq*x**(-1)*
     &    lnopx - 47.D0/9.D0*CF*CA*pisq - 6*CF*CA*pisq*lnmuonpt + 2*CF*
     &    CA*pisq*lnopx + 7.D0/9.D0*CF*CA*pisq*omx - 2.D0/3.D0*CF*CA*
     &    pisq*omx*lnmuonpt - 2.D0/3.D0*CF*CA*pisq*omx*lnopx - 8.D0/9.D0
     &    *CF*CA*pisq*omx**2
      I2P(gq,rglr) = I2P(gq,rglr) + 48*CF*CA*zeta3*x**(-1) - 32*CF*CA*
     &    zeta3 - 24*CF*CA*zeta3*omx + 32*CF**2*x**(-1)*lnomx + 48*
     &    CF**2*x**(-1)*lnomx*lnmuonpt + 32*CF**2*x**(-1)*lnomx*
     &    lnmuonpt**2 + 6*CF**2*x**(-1)*lnomx**2 + 16*CF**2*x**(-1)*
     &    lnomx**2*lnmuonpt + 4.D0/3.D0*CF**2*x**(-1)*lnomx**3 + 9*
     &    CF**2 + 36*CF**2*lnmuonpt + 12*CF**2*lnmuonpt**2 - 22*CF**2*
     &    lnomx - 24*CF**2*lnomx*lnmuonpt - 16*CF**2*lnomx*lnmuonpt**2
     &     - 5*CF**2*lnomx**2 - 8*CF**2*lnomx**2*lnmuonpt - 2.D0/3.D0*
     &    CF**2*lnomx**3 - 10*CF**2*lnx - 36*CF**2*lnx*lnmuonpt + 8*
     &    CF**2*lnx*lnmuonpt**2 - 7.D0/2.D0*CF**2*lnx**2 + 4*CF**2*
     &    lnx**2*lnmuonpt + 1.D0/3.D0*CF**2*lnx**3 + CF**2*omx - 24*
     &    CF**2*omx*lnmuonpt + 4*CF**2*omx*lnmuonpt**2 - 10*CF**2*omx*
     &    lnomx - 24*CF**2*omx*lnomx*lnmuonpt - 16*CF**2*omx*lnomx*
     &    lnmuonpt**2 - CF**2*omx*lnomx**2 - 8*CF**2*omx*lnomx**2*
     &    lnmuonpt - 2.D0/3.D0*CF**2*omx*lnomx**3 - 5*CF**2*omx*lnx + 
     &    20*CF**2*omx*lnx*lnmuonpt
      I2P(gq,rglr) = I2P(gq,rglr) + 8*CF**2*omx*lnx*lnmuonpt**2 + 3.D0/
     &    2.D0*CF**2*omx*lnx**2 + 4*CF**2*omx*lnx**2*lnmuonpt + 1.D0/3.D
     &    0*CF**2*omx*lnx**3

