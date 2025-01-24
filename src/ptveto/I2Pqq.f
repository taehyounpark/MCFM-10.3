!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2P(qq,delt)= - 164.D0/81.D0*CF*nfl + 2.D0/3.D0*CF*lnmuonpt*nfl
     &     + 4*CF*lnmuonpt**2*nfl + 5.D0/18.D0*CF*pisq*nfl + 10.D0/9.D0
     &    *CF*pisq*lnmuonpt*nfl - 2.D0/9.D0*CF*zeta3*nfl + 1214.D0/81.D0
     &    *CF*CA - 17.D0/3.D0*CF*CA*lnmuonpt - 22*CF*CA*lnmuonpt**2 - 
     &    67.D0/36.D0*CF*CA*pisq - 55.D0/9.D0*CF*CA*pisq*lnmuonpt + 1.D0
     &    /18.D0*CF*CA*pisq**2 + 11.D0/9.D0*CF*CA*zeta3 + 24*CF*CA*
     &    zeta3*lnmuonpt - 3*CF**2*lnmuonpt + 18*CF**2*lnmuonpt**2 + 5*
     &    CF**2*pisq*lnmuonpt - 16.D0/3.D0*CF**2*pisq*lnmuonpt**2 + 1.D0
     &    /72.D0*CF**2*pisq**2 - 80*CF**2*zeta3*lnmuonpt

      I2P(qq,plus)= + 112.D0/27.D0*CF*nfl + 80.D0/9.D0*CF*lnmuonpt*nfl
     &     + 16.D0/3.D0*CF*lnmuonpt**2*nfl - 808.D0/27.D0*CF*CA - 536.D0
     &    /9.D0*CF*CA*lnmuonpt - 88.D0/3.D0*CF*CA*lnmuonpt**2 + 8.D0/3.D
     &    0*CF*CA*pisq*lnmuonpt + 28*CF*CA*zeta3 + 48*CF**2*lnmuonpt**2
     &     + 4.D0/3.D0*CF**2*pisq*lnmuonpt

      I2P(qq,lpls)= + 64*CF**2*lnmuonpt**2

      I2P(qq,rglr)= + 20.D0/9.D0*CF*omx**(-1)*lnx*nfl + 16.D0/3.D0*CF*
     &    omx**(-1)*lnx*lnmuonpt*nfl + 2.D0/3.D0*CF*omx**(-1)*lnx**2*
     &    nfl + 172.D0/27.D0*CF*x**(-1) - 104.D0/9.D0*CF*x**(-1)*
     &    lnmuonpt + 16.D0/3.D0*CF*x**(-1)*lnmuonpt**2 + 8.D0/3.D0*CF*
     &    x**(-1)*lnx*lnomx + 8.D0/3.D0*CF*x**(-1)*Li2x - 172.D0/27.D0*
     &    CF - 112.D0/27.D0*CF*nfl + 104.D0/9.D0*CF*lnmuonpt - 80.D0/9.D
     &    0*CF*lnmuonpt*nfl - 16.D0/3.D0*CF*lnmuonpt**2 - 16.D0/3.D0*CF
     &    *lnmuonpt**2*nfl + 46.D0/9.D0*CF*lnx - 20.D0/9.D0*CF*lnx*nfl
     &     - 80.D0/3.D0*CF*lnx*lnmuonpt - 16.D0/3.D0*CF*lnx*lnmuonpt*
     &    nfl + 16*CF*lnx*lnmuonpt**2 - 8.D0/3.D0*CF*lnx*lnomx - 7.D0/3.
     &    D0*CF*lnx**2 - 2.D0/3.D0*CF*lnx**2*nfl + 8*CF*lnx**2*lnmuonpt
     &     + 2.D0/3.D0*CF*lnx**3 - 8.D0/3.D0*CF*Li2x - 7.D0/27.D0*CF*
     &    omx + 38.D0/27.D0*CF*omx*nfl - 136.D0/9.D0*CF*omx*lnmuonpt + 
     &    64.D0/9.D0*CF*omx*lnmuonpt*nfl + 44.D0/3.D0*CF*omx*
     &    lnmuonpt**2 + 8.D0/3.D0*CF*omx*lnmuonpt**2*nfl - 68.D0/9.D0*
     &    CF*omx*lnx
      I2P(qq,rglr) = I2P(qq,rglr) + 10.D0/9.D0*CF*omx*lnx*nfl + 100.D0/
     &    3.D0*CF*omx*lnx*lnmuonpt + 8.D0/3.D0*CF*omx*lnx*lnmuonpt*nfl
     &     - 8*CF*omx*lnx*lnmuonpt**2 + 4.D0/3.D0*CF*omx*lnx*lnomx + 19.
     &    D0/6.D0*CF*omx*lnx**2 + 1.D0/3.D0*CF*omx*lnx**2*nfl - 4*CF*
     &    omx*lnx**2*lnmuonpt - 1.D0/3.D0*CF*omx*lnx**3 + 4.D0/3.D0*CF*
     &    omx*Li2x - 136.D0/27.D0*CF*omx**2 + 176.D0/9.D0*CF*omx**2*
     &    lnmuonpt - 16.D0/3.D0*CF*omx**2*lnmuonpt**2 + 64.D0/9.D0*CF*
     &    omx**2*lnx - 32.D0/3.D0*CF*omx**2*lnx*lnmuonpt - 8.D0/3.D0*CF
     &    *omx**2*lnx*lnomx - 4.D0/3.D0*CF*omx**2*lnx**2 - 8.D0/3.D0*CF
     &    *omx**2*Li2x - 4.D0/9.D0*CF*pisq*x**(-1) + 4.D0/9.D0*CF*pisq
     &     - 2.D0/9.D0*CF*pisq*omx + 4.D0/9.D0*CF*pisq*omx**2 - 152.D0/
     &    9.D0*CF*CA*omx**(-1)*lnx - 88.D0/3.D0*CF*CA*omx**(-1)*lnx*
     &    lnmuonpt - 11.D0/3.D0*CF*CA*omx**(-1)*lnx**2 - 8*CF*CA*
     &    omx**(-1)*lnx**2*lnmuonpt - 2.D0/3.D0*CF*CA*omx**(-1)*lnx**3
     &     + 16*CF*CA*omx**(-1)*Li3x - 8*CF*CA*omx**(-1)*Li3omx + 8*CF*
     &    CA*omx**(-1)*Li2omx*lnomx
      I2P(qq,rglr) = I2P(qq,rglr) - 8*CF*CA*omx**(-1)*Li2x*lnx + 808.D0/
     &    27.D0*CF*CA + 536.D0/9.D0*CF*CA*lnmuonpt + 88.D0/3.D0*CF*CA*
     &    lnmuonpt**2 + 2*CF*CA*lnomx + 260.D0/9.D0*CF*CA*lnx + 40.D0/3.
     &    D0*CF*CA*lnx*lnmuonpt + 5.D0/3.D0*CF*CA*lnx**2 + 8*CF*CA*
     &    lnx**2*lnmuonpt + 2.D0/3.D0*CF*CA*lnx**3 - 16*CF*CA*Li3x + 8*
     &    CF*CA*Li3omx - 8*CF*CA*Li2omx*lnomx + 8*CF*CA*Li2x*lnx - 8.D0/
     &    27.D0*CF*CA*omx - 616.D0/9.D0*CF*CA*omx*lnmuonpt - 44.D0/3.D0
     &    *CF*CA*omx*lnmuonpt**2 - 2*CF*CA*omx*lnomx - 166.D0/9.D0*CF*
     &    CA*omx*lnx - 20.D0/3.D0*CF*CA*omx*lnx*lnmuonpt + 4*CF*CA*omx*
     &    lnx*lnomx + 1.D0/6.D0*CF*CA*omx*lnx**2 - 4*CF*CA*omx*lnx**2*
     &    lnmuonpt - 1.D0/3.D0*CF*CA*omx*lnx**3 + 8*CF*CA*omx*Li3x - 4*
     &    CF*CA*omx*Li3omx + 4*CF*CA*omx*Li2omx*lnomx + 4*CF*CA*omx*
     &    Li2x - 4*CF*CA*omx*Li2x*lnx - 8.D0/3.D0*CF*CA*pisq*lnmuonpt
     &     - CF*CA*pisq*omx + 4.D0/3.D0*CF*CA*pisq*omx*lnmuonpt - 16*CF
     &    *CA*zeta3*omx**(-1) - 12*CF*CA*zeta3 + 6*CF*CA*zeta3*omx + 16
     &    *CF**2*omx**(-1)*lnx
      I2P(qq,rglr) = I2P(qq,rglr) + 24*CF**2*omx**(-1)*lnx*lnmuonpt - 
     &    32*CF**2*omx**(-1)*lnx*lnmuonpt**2 + 32*CF**2*omx**(-1)*lnx*
     &    lnomx*lnmuonpt + 4*CF**2*omx**(-1)*lnx*lnomx**2 + 3*CF**2*
     &    omx**(-1)*lnx**2 + 4*CF**2*omx**(-1)*lnx**2*lnomx - 40*CF**2*
     &    omx**(-1)*Li3x + 8*CF**2*omx**(-1)*Li3omx - 8*CF**2*omx**(-1)
     &    *Li2omx*lnomx + 24*CF**2*omx**(-1)*Li2x*lnx - 48*CF**2*
     &    lnmuonpt**2 - 2*CF**2*lnomx - 64*CF**2*lnomx*lnmuonpt**2 - 38
     &    *CF**2*lnx + 16*CF**2*lnx*lnmuonpt + 48*CF**2*lnx*lnmuonpt**2
     &     - 32*CF**2*lnx*lnomx*lnmuonpt - 4*CF**2*lnx*lnomx**2 + 2*
     &    CF**2*lnx**2 + 8*CF**2*lnx**2*lnmuonpt - 4*CF**2*lnx**2*lnomx
     &     + 2.D0/3.D0*CF**2*lnx**3 + 40*CF**2*Li3x - 8*CF**2*Li3omx + 
     &    8*CF**2*Li2omx*lnomx - 24*CF**2*Li2x*lnx - 22*CF**2*omx + 44*
     &    CF**2*omx*lnmuonpt + 8*CF**2*omx*lnmuonpt**2 + 2*CF**2*omx*
     &    lnomx - 16*CF**2*omx*lnomx*lnmuonpt + 32*CF**2*omx*lnomx*
     &    lnmuonpt**2 + 32*CF**2*omx*lnx - 8*CF**2*omx*lnx*lnmuonpt - 
     &    24*CF**2*omx*lnx*lnmuonpt**2
      I2P(qq,rglr) = I2P(qq,rglr) - 12*CF**2*omx*lnx*lnomx + 16*CF**2*
     &    omx*lnx*lnomx*lnmuonpt + 2*CF**2*omx*lnx*lnomx**2 - 2*CF**2*
     &    omx*lnx**2 - 4*CF**2*omx*lnx**2*lnmuonpt + 2*CF**2*omx*lnx**2
     &    *lnomx - 1.D0/3.D0*CF**2*omx*lnx**3 - 20*CF**2*omx*Li3x + 4*
     &    CF**2*omx*Li3omx - 4*CF**2*omx*Li2omx*lnomx - 8*CF**2*omx*
     &    Li2x + 12*CF**2*omx*Li2x*lnx - 4.D0/3.D0*CF**2*pisq*lnmuonpt
     &     + CF**2*pisq*omx + 2.D0/3.D0*CF**2*pisq*omx*lnmuonpt + 40*
     &    CF**2*zeta3*omx**(-1) - 40*CF**2*zeta3 + 20*CF**2*zeta3*omx

