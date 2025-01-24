!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2P(gg,delt)= - 164.D0/81.D0*CA*nfl + 16.D0/3.D0*CA*lnmuonpt*nfl
     &     + 5.D0/18.D0*CA*pisq*nfl - 2.D0/9.D0*CA*zeta3*nfl + 1214.D0/
     &    81.D0*CA**2 - 64.D0/3.D0*CA**2*lnmuonpt - 67.D0/36.D0*CA**2*
     &    pisq - 16.D0/3.D0*CA**2*pisq*lnmuonpt**2 + 5.D0/72.D0*CA**2*
     &    pisq**2 + 11.D0/9.D0*CA**2*zeta3 - 56*CA**2*zeta3*lnmuonpt + 
     &    4*CF*lnmuonpt*nfl

      I2P(gg,plus)= + 112.D0/27.D0*CA*nfl + 80.D0/9.D0*CA*lnmuonpt*nfl
     &     - 16.D0/3.D0*CA*lnmuonpt**2*nfl - 808.D0/27.D0*CA**2 - 536.D0
     &    /9.D0*CA**2*lnmuonpt + 88.D0/3.D0*CA**2*lnmuonpt**2 + 4*CA**2
     &    *pisq*lnmuonpt + 28*CA**2*zeta3

      I2P(gg,lpls)= + 64*CA**2*lnmuonpt**2

      I2P(gg,rglr)= + 242.D0/27.D0*CA*x**(-1)*nfl + 184.D0/9.D0*CA*
     &    x**(-1)*lnmuonpt*nfl - 16.D0/3.D0*CA*x**(-1)*lnmuonpt**2*nfl
     &     - 130.D0/9.D0*CA*nfl - 88.D0/3.D0*CA*lnmuonpt*nfl + 32.D0/3.D
     &    0*CA*lnmuonpt**2*nfl - 2.D0/3.D0*CA*lnomx*nfl + 46.D0/9.D0*CA
     &    *lnx*nfl + 32.D0/3.D0*CA*lnx*lnmuonpt*nfl + 4.D0/3.D0*CA*
     &    lnx**2*nfl + 112.D0/9.D0*CA*omx*nfl + 24*CA*omx*lnmuonpt*nfl
     &     - 16.D0/3.D0*CA*omx*lnmuonpt**2*nfl + 2.D0/3.D0*CA*omx*lnomx
     &    *nfl - 20.D0/9.D0*CA*omx*lnx*nfl - 16.D0/3.D0*CA*omx*lnx*
     &    lnmuonpt*nfl - 2.D0/3.D0*CA*omx*lnx**2*nfl - 278.D0/27.D0*CA*
     &    omx**2*nfl - 184.D0/9.D0*CA*omx**2*lnmuonpt*nfl + 16.D0/3.D0*
     &    CA*omx**2*lnmuonpt**2*nfl - 32*CA**2*omx**(-1)*lnx*
     &    lnmuonpt**2 + 32*CA**2*omx**(-1)*lnx*lnomx*lnmuonpt + 4*CA**2
     &    *omx**(-1)*lnx*lnomx**2 - 8*CA**2*omx**(-1)*lnx**2*lnmuonpt
     &     + 4*CA**2*omx**(-1)*lnx**2*lnomx - 2.D0/3.D0*CA**2*omx**(-1)
     &    *lnx**3 - 24*CA**2*omx**(-1)*Li3x + 16*CA**2*omx**(-1)*Li2x*
     &    lnx
      I2P(gg,rglr) = I2P(gg,rglr) - 3160.D0/27.D0*CA**2*x**(-1) - 8*
     &    CA**2*x**(-1)*Li3mx - 88*CA**2*x**(-1)*lnmuonpt**2 - 8.D0/3.D0
     &    *CA**2*x**(-1)*lnopx**3 - 32*CA**2*x**(-1)*Li2mx*lnmuonpt + 
     &    64*CA**2*x**(-1)*lnomx*lnmuonpt**2 - 32*CA**2*x**(-1)*lnx*
     &    lnmuonpt**2 - 32*CA**2*x**(-1)*lnx*lnopx*lnmuonpt + 8*CA**2*
     &    x**(-1)*lnx*lnopx**2 + 8*CA**2*x**(-1)*lnx*Li2mx - 88.D0/3.D0
     &    *CA**2*x**(-1)*lnx*lnomx + 32*CA**2*x**(-1)*lnx*lnomx*
     &    lnmuonpt + 4*CA**2*x**(-1)*lnx*lnomx**2 - 4*CA**2*x**(-1)*
     &    lnx**2*lnopx + 4*CA**2*x**(-1)*lnx**2*lnomx + 16*CA**2*
     &    x**(-1)*Li3xonopx - 40*CA**2*x**(-1)*Li3x - 88.D0/3.D0*CA**2*
     &    x**(-1)*Li2x + 24*CA**2*x**(-1)*Li2x*lnx + 8*CA**2*opx**(-1)*
     &    Li3mx + 8.D0/3.D0*CA**2*opx**(-1)*lnopx**3 + 32*CA**2*
     &    opx**(-1)*Li2mx*lnmuonpt + 32*CA**2*opx**(-1)*lnx*lnopx*
     &    lnmuonpt - 8*CA**2*opx**(-1)*lnx*lnopx**2 - 8*CA**2*opx**(-1)
     &    *lnx*Li2mx - 8*CA**2*opx**(-1)*lnx**2*lnmuonpt + 4*CA**2*
     &    opx**(-1)*lnx**2*lnopx
      I2P(gg,rglr) = I2P(gg,rglr) - 2.D0/3.D0*CA**2*opx**(-1)*lnx**3 - 
     &    16*CA**2*opx**(-1)*Li3xonopx + 16*CA**2*opx**(-1)*Li3x - 8*
     &    CA**2*opx**(-1)*Li2x*lnx + 4004.D0/27.D0*CA**2 - 32*CA**2*
     &    Li3mx + 536.D0/9.D0*CA**2*lnmuonpt + 176.D0/3.D0*CA**2*
     &    lnmuonpt**2 - 32.D0/3.D0*CA**2*lnopx**3 - 128*CA**2*Li2mx*
     &    lnmuonpt + 2.D0/3.D0*CA**2*lnomx - 128*CA**2*lnomx*
     &    lnmuonpt**2 - 154*CA**2*lnx + 464.D0/3.D0*CA**2*lnx*lnmuonpt
     &     - 64*CA**2*lnx*lnmuonpt**2 - 128*CA**2*lnx*lnopx*lnmuonpt + 
     &    32*CA**2*lnx*lnopx**2 + 32*CA**2*lnx*Li2mx + 88.D0/3.D0*CA**2
     &    *lnx*lnomx - 64*CA**2*lnx*lnomx*lnmuonpt - 8*CA**2*lnx*
     &    lnomx**2 + 58.D0/3.D0*CA**2*lnx**2 - 16*CA**2*lnx**2*lnmuonpt
     &     - 16*CA**2*lnx**2*lnopx - 8*CA**2*lnx**2*lnomx - 4.D0/3.D0*
     &    CA**2*lnx**3 + 64*CA**2*Li3xonopx - 16*CA**2*Li3x + 88.D0/3.D0
     &    *CA**2*Li2x - 3640.D0/27.D0*CA**2*omx + 24*CA**2*omx*Li3mx - 
     &    436.D0/9.D0*CA**2*omx*lnmuonpt - 328.D0/3.D0*CA**2*omx*
     &    lnmuonpt**2
      I2P(gg,rglr) = I2P(gg,rglr) + 8*CA**2*omx*lnopx**3 + 96*CA**2*omx
     &    *Li2mx*lnmuonpt - 2.D0/3.D0*CA**2*omx*lnomx + 64*CA**2*omx*
     &    lnomx*lnmuonpt**2 + 407.D0/3.D0*CA**2*omx*lnx - 616.D0/3.D0*
     &    CA**2*omx*lnx*lnmuonpt + 32*CA**2*omx*lnx*lnmuonpt**2 + 96*
     &    CA**2*omx*lnx*lnopx*lnmuonpt - 24*CA**2*omx*lnx*lnopx**2 - 24
     &    *CA**2*omx*lnx*Li2mx - 80.D0/3.D0*CA**2*omx*lnx*lnomx + 32*
     &    CA**2*omx*lnx*lnomx*lnmuonpt + 4*CA**2*omx*lnx*lnomx**2 - 77.D
     &    0/3.D0*CA**2*omx*lnx**2 + 12*CA**2*omx*lnx**2*lnopx + 4*CA**2
     &    *omx*lnx**2*lnomx - 48*CA**2*omx*Li3xonopx + 24*CA**2*omx*
     &    Li3x - 80.D0/3.D0*CA**2*omx*Li2x - 8*CA**2*omx*Li2x*lnx + 
     &    3340.D0/27.D0*CA**2*omx**2 - 8*CA**2*omx**2*Li3mx + 88*CA**2*
     &    omx**2*lnmuonpt**2 - 8.D0/3.D0*CA**2*omx**2*lnopx**3 - 32*
     &    CA**2*omx**2*Li2mx*lnmuonpt - 64*CA**2*omx**2*lnomx*
     &    lnmuonpt**2 - 536.D0/9.D0*CA**2*omx**2*lnx + 352.D0/3.D0*
     &    CA**2*omx**2*lnx*lnmuonpt + 32*CA**2*omx**2*lnx*lnmuonpt**2
     &     - 32*CA**2*omx**2*lnx*lnopx*lnmuonpt
      I2P(gg,rglr) = I2P(gg,rglr) + 8*CA**2*omx**2*lnx*lnopx**2 + 8*
     &    CA**2*omx**2*lnx*Li2mx + 88.D0/3.D0*CA**2*omx**2*lnx*lnomx - 
     &    32*CA**2*omx**2*lnx*lnomx*lnmuonpt - 4*CA**2*omx**2*lnx*
     &    lnomx**2 + 44.D0/3.D0*CA**2*omx**2*lnx**2 + 16*CA**2*omx**2*
     &    lnx**2*lnmuonpt - 4*CA**2*omx**2*lnx**2*lnopx - 4*CA**2*
     &    omx**2*lnx**2*lnomx + 4.D0/3.D0*CA**2*omx**2*lnx**3 + 16*
     &    CA**2*omx**2*Li3xonopx + 8*CA**2*omx**2*Li3x + 88.D0/3.D0*
     &    CA**2*omx**2*Li2x - 8*CA**2*omx**2*Li2x*lnx + 44.D0/9.D0*
     &    CA**2*pisq*x**(-1) + 4.D0/3.D0*CA**2*pisq*x**(-1)*lnmuonpt + 
     &    4.D0/3.D0*CA**2*pisq*x**(-1)*lnopx + 8.D0/3.D0*CA**2*pisq*
     &    opx**(-1)*lnmuonpt - 4.D0/3.D0*CA**2*pisq*opx**(-1)*lnopx - 
     &    44.D0/9.D0*CA**2*pisq - 56.D0/3.D0*CA**2*pisq*lnmuonpt + 16.D0
     &    /3.D0*CA**2*pisq*lnopx + 40.D0/9.D0*CA**2*pisq*omx + 12*CA**2
     &    *pisq*omx*lnmuonpt - 4*CA**2*pisq*omx*lnopx - 44.D0/9.D0*
     &    CA**2*pisq*omx**2 - 20.D0/3.D0*CA**2*pisq*omx**2*lnmuonpt + 4.
     &    D0/3.D0*CA**2*pisq*omx**2*lnopx
      I2P(gg,rglr) = I2P(gg,rglr) + 24*CA**2*zeta3*omx**(-1) + 48*CA**2
     &    *zeta3*x**(-1) + 4*CA**2*zeta3*opx**(-1) - 120*CA**2*zeta3 + 
     &    64*CA**2*zeta3*omx - 56*CA**2*zeta3*omx**2 - 4.D0/3.D0*CF*
     &    x**(-1)*nfl - 16.D0/3.D0*CF*x**(-1)*lnmuonpt*nfl + 32.D0/3.D0
     &    *CF*x**(-1)*lnmuonpt**2*nfl + 4.D0/3.D0*CF*nfl + 16.D0/3.D0*
     &    CF*lnmuonpt*nfl - 32.D0/3.D0*CF*lnmuonpt**2*nfl + 24*CF*lnx*
     &    nfl + 48*CF*lnx*lnmuonpt*nfl + 32*CF*lnx*lnmuonpt**2*nfl + 4*
     &    CF*lnx**2*nfl + 16*CF*lnx**2*lnmuonpt*nfl + 4.D0/3.D0*CF*
     &    lnx**3*nfl + 88.D0/3.D0*CF*omx*nfl + 184.D0/3.D0*CF*omx*
     &    lnmuonpt*nfl + 88.D0/3.D0*CF*omx*lnmuonpt**2*nfl - 12*CF*omx*
     &    lnx*nfl - 24*CF*omx*lnx*lnmuonpt*nfl - 16*CF*omx*lnx*
     &    lnmuonpt**2*nfl - CF*omx*lnx**2*nfl - 8*CF*omx*lnx**2*
     &    lnmuonpt*nfl - 2.D0/3.D0*CF*omx*lnx**3*nfl + 4.D0/3.D0*CF*
     &    omx**2*nfl - 32.D0/3.D0*CF*omx**2*lnmuonpt*nfl - 32.D0/3.D0*
     &    CF*omx**2*lnmuonpt**2*nfl

