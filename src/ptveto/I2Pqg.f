!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2P(qg,rglr)= + 172.D0/27.D0*CA*x**(-1) - 104.D0/9.D0*CA*x**(-1)*
     &    lnmuonpt + 16.D0/3.D0*CA*x**(-1)*lnmuonpt**2 + 8.D0/3.D0*CA*
     &    x**(-1)*lnx*lnomx + 8.D0/3.D0*CA*x**(-1)*Li2x - 226.D0/27.D0*
     &    CA + 10*CA*Li3mx + 68.D0/9.D0*CA*lnmuonpt - 16.D0/3.D0*CA*
     &    lnmuonpt**2 + 10.D0/3.D0*CA*lnopx**3 + 8*CA*Li2mx + 40*CA*
     &    Li2mx*lnmuonpt - CA*lnomx + 8*CA*lnomx*lnmuonpt**2 + 4*CA*
     &    lnomx**2*lnmuonpt + 1.D0/3.D0*CA*lnomx**3 + 118.D0/9.D0*CA*
     &    lnx - 188.D0/3.D0*CA*lnx*lnmuonpt + 40*CA*lnx*lnmuonpt**2 + 8
     &    *CA*lnx*lnopx + 40*CA*lnx*lnopx*lnmuonpt - 10*CA*lnx*lnopx**2
     &     - 10*CA*lnx*Li2mx - 8.D0/3.D0*CA*lnx*lnomx - 35.D0/6.D0*CA*
     &    lnx**2 + 12*CA*lnx**2*lnmuonpt + 5*CA*lnx**2*lnopx + CA*
     &    lnx**3 - 20*CA*Li3xonopx + 16*CA*Li3x + 2*CA*Li3omx - 2*CA*
     &    Li2omx*lnomx - 8.D0/3.D0*CA*Li2x - 8*CA*Li2x*lnx + 209.D0/27.D
     &    0*CA*omx - 12*CA*omx*Li3mx - 388.D0/9.D0*CA*omx*lnmuonpt + 
     &    152.D0/3.D0*CA*omx*lnmuonpt**2 - 4*CA*omx*lnopx**3 - 12*CA*
     &    omx*Li2mx
      I2P(qg,rglr) = I2P(qg,rglr) - 48*CA*omx*Li2mx*lnmuonpt + 5*CA*omx
     &    *lnomx - 16*CA*omx*lnomx*lnmuonpt**2 - 2*CA*omx*lnomx**2 - 8*
     &    CA*omx*lnomx**2*lnmuonpt - 2.D0/3.D0*CA*omx*lnomx**3 - 212.D0/
     &    9.D0*CA*omx*lnx + 352.D0/3.D0*CA*omx*lnx*lnmuonpt - 32*CA*omx
     &    *lnx*lnmuonpt**2 - 12*CA*omx*lnx*lnopx - 48*CA*omx*lnx*lnopx*
     &    lnmuonpt + 12*CA*omx*lnx*lnopx**2 + 12*CA*omx*lnx*Li2mx + 40.D
     &    0/3.D0*CA*omx*lnx*lnomx + 38.D0/3.D0*CA*omx*lnx**2 - 8*CA*omx
     &    *lnx**2*lnmuonpt - 6*CA*omx*lnx**2*lnopx - 2.D0/3.D0*CA*omx*
     &    lnx**3 + 24*CA*omx*Li3xonopx - 16*CA*omx*Li3x - 4*CA*omx*
     &    Li3omx + 4*CA*omx*Li2omx*lnomx + 40.D0/3.D0*CA*omx*Li2x + 8*
     &    CA*omx*Li2x*lnx - 298.D0/27.D0*CA*omx**2 + 4*CA*omx**2*Li3mx
     &     + 464.D0/9.D0*CA*omx**2*lnmuonpt - 124.D0/3.D0*CA*omx**2*
     &    lnmuonpt**2 + 4.D0/3.D0*CA*omx**2*lnopx**3 + 4*CA*omx**2*
     &    Li2mx + 16*CA*omx**2*Li2mx*lnmuonpt - 4*CA*omx**2*lnomx + 16*
     &    CA*omx**2*lnomx*lnmuonpt**2 + 2*CA*omx**2*lnomx**2 + 8*CA*
     &    omx**2*lnomx**2*lnmuonpt
      I2P(qg,rglr) = I2P(qg,rglr) + 2.D0/3.D0*CA*omx**2*lnomx**3 + 136.D
     &    0/9.D0*CA*omx**2*lnx - 176.D0/3.D0*CA*omx**2*lnx*lnmuonpt + 4
     &    *CA*omx**2*lnx*lnopx + 16*CA*omx**2*lnx*lnopx*lnmuonpt - 4*CA
     &    *omx**2*lnx*lnopx**2 - 4*CA*omx**2*lnx*Li2mx - 44.D0/3.D0*CA*
     &    omx**2*lnx*lnomx - 22.D0/3.D0*CA*omx**2*lnx**2 + 2*CA*omx**2*
     &    lnx**2*lnopx - 8*CA*omx**2*Li3xonopx + 4*CA*omx**2*Li3omx - 4
     &    *CA*omx**2*Li2omx*lnomx - 44.D0/3.D0*CA*omx**2*Li2x - 4.D0/9.D
     &    0*CA*pisq*x**(-1) + 10.D0/9.D0*CA*pisq + 8.D0/3.D0*CA*pisq*
     &    lnmuonpt - 5.D0/3.D0*CA*pisq*lnopx - 26.D0/9.D0*CA*pisq*omx
     &     - 8.D0/3.D0*CA*pisq*omx*lnmuonpt + 2*CA*pisq*omx*lnopx + 22.D
     &    0/9.D0*CA*pisq*omx**2 - 2.D0/3.D0*CA*pisq*omx**2*lnopx + 2*CA
     &    *zeta3 - 4*CA*zeta3*omx**2 - 5*CF - 10*CF*lnmuonpt + 6*CF*
     &    lnmuonpt**2 + CF*lnomx + 8*CF*lnomx*lnmuonpt**2 - 4*CF*
     &    lnomx**2*lnmuonpt - 1.D0/3.D0*CF*lnomx**3 + 15.D0/2.D0*CF*lnx
     &     - 2*CF*lnx*lnmuonpt - 12*CF*lnx*lnmuonpt**2 + 8*CF*lnx*lnomx
     &    *lnmuonpt
      I2P(qg,rglr) = I2P(qg,rglr) + CF*lnx*lnomx**2 + 5.D0/4.D0*CF*
     &    lnx**2 - 6*CF*lnx**2*lnmuonpt + CF*lnx**2*lnomx - 1.D0/2.D0*
     &    CF*lnx**3 - 2*CF*Li3x - 2*CF*Li3omx + 2*CF*Li2omx*lnomx + 2*
     &    CF*Li2x*lnx + 69.D0/2.D0*CF*omx + 34*CF*omx*lnmuonpt - 8*CF*
     &    omx*lnmuonpt**2 - 5*CF*omx*lnomx - 16*CF*omx*lnomx*lnmuonpt
     &     - 16*CF*omx*lnomx*lnmuonpt**2 + 2*CF*omx*lnomx**2 + 8*CF*omx
     &    *lnomx**2*lnmuonpt + 2.D0/3.D0*CF*omx*lnomx**3 + 1.D0/2.D0*CF
     &    *omx*lnx + 16*CF*omx*lnx*lnmuonpt + 24*CF*omx*lnx*lnmuonpt**2
     &     - 4*CF*omx*lnx*lnomx - 16*CF*omx*lnx*lnomx*lnmuonpt - 2*CF*
     &    omx*lnx*lnomx**2 + CF*omx*lnx**2 + 12*CF*omx*lnx**2*lnmuonpt
     &     - 2*CF*omx*lnx**2*lnomx + CF*omx*lnx**3 + 4*CF*omx*Li3x + 4*
     &    CF*omx*Li3omx - 4*CF*omx*Li2omx*lnomx - 4*CF*omx*Li2x*lnx - 
     &    36*CF*omx**2 - 44*CF*omx**2*lnmuonpt + 4*CF*omx**2*lnomx + 16
     &    *CF*omx**2*lnomx*lnmuonpt + 16*CF*omx**2*lnomx*lnmuonpt**2 - 
     &    2*CF*omx**2*lnomx**2 - 8*CF*omx**2*lnomx**2*lnmuonpt - 2.D0/3.
     &    D0*CF*omx**2*lnomx**3
      I2P(qg,rglr) = I2P(qg,rglr) - 4*CF*omx**2*lnx - 16*CF*omx**2*lnx*
     &    lnmuonpt - 16*CF*omx**2*lnx*lnmuonpt**2 + 4*CF*omx**2*lnx*
     &    lnomx + 16*CF*omx**2*lnx*lnomx*lnmuonpt + 2*CF*omx**2*lnx*
     &    lnomx**2 - 2*CF*omx**2*lnx**2 - 8*CF*omx**2*lnx**2*lnmuonpt
     &     + 2*CF*omx**2*lnx**2*lnomx - 2.D0/3.D0*CF*omx**2*lnx**3 - 4*
     &    CF*omx**2*Li3x - 4*CF*omx**2*Li3omx + 4*CF*omx**2*Li2omx*
     &    lnomx + 4*CF*omx**2*Li2x*lnx + 5.D0/3.D0*CF*pisq*lnmuonpt - 
     &    CF*pisq*omx - 10.D0/3.D0*CF*pisq*omx*lnmuonpt + CF*pisq*
     &    omx**2 + 10.D0/3.D0*CF*pisq*omx**2*lnmuonpt + 16*CF*zeta3 - 
     &    32*CF*zeta3*omx + 32*CF*zeta3*omx**2

