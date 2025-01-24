!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2P(qbq,rglr)= + 172.D0/27.D0*CF*x**(-1) - 104.D0/9.D0*CF*x**(-1)
     &    *lnmuonpt + 16.D0/3.D0*CF*x**(-1)*lnmuonpt**2 + 8.D0/3.D0*CF*
     &    x**(-1)*lnx*lnomx + 8.D0/3.D0*CF*x**(-1)*Li2x - 172.D0/27.D0*
     &    CF + 104.D0/9.D0*CF*lnmuonpt - 16.D0/3.D0*CF*lnmuonpt**2 + 46.
     &    D0/9.D0*CF*lnx - 80.D0/3.D0*CF*lnx*lnmuonpt + 16*CF*lnx*
     &    lnmuonpt**2 - 8.D0/3.D0*CF*lnx*lnomx - 7.D0/3.D0*CF*lnx**2 + 
     &    8*CF*lnx**2*lnmuonpt + 2.D0/3.D0*CF*lnx**3 - 8.D0/3.D0*CF*
     &    Li2x - 7.D0/27.D0*CF*omx - 136.D0/9.D0*CF*omx*lnmuonpt + 44.D0
     &    /3.D0*CF*omx*lnmuonpt**2 - 68.D0/9.D0*CF*omx*lnx + 100.D0/3.D0
     &    *CF*omx*lnx*lnmuonpt - 8*CF*omx*lnx*lnmuonpt**2 + 4.D0/3.D0*
     &    CF*omx*lnx*lnomx + 19.D0/6.D0*CF*omx*lnx**2 - 4*CF*omx*lnx**2
     &    *lnmuonpt - 1.D0/3.D0*CF*omx*lnx**3 + 4.D0/3.D0*CF*omx*Li2x
     &     - 136.D0/27.D0*CF*omx**2 + 176.D0/9.D0*CF*omx**2*lnmuonpt - 
     &    16.D0/3.D0*CF*omx**2*lnmuonpt**2 + 64.D0/9.D0*CF*omx**2*lnx
     &     - 32.D0/3.D0*CF*omx**2*lnx*lnmuonpt - 8.D0/3.D0*CF*omx**2*
     &    lnx*lnomx
      I2P(qbq,rglr) = I2P(qbq,rglr) - 4.D0/3.D0*CF*omx**2*lnx**2 - 8.D0/
     &    3.D0*CF*omx**2*Li2x - 4.D0/9.D0*CF*pisq*x**(-1) + 4.D0/9.D0*
     &    CF*pisq - 2.D0/9.D0*CF*pisq*omx + 4.D0/9.D0*CF*pisq*omx**2 - 
     &    8*CF*CA*opx**(-1)*Li3mx - 8.D0/3.D0*CF*CA*opx**(-1)*lnopx**3
     &     - 32*CF*CA*opx**(-1)*Li2mx*lnmuonpt - 32*CF*CA*opx**(-1)*lnx
     &    *lnopx*lnmuonpt + 8*CF*CA*opx**(-1)*lnx*lnopx**2 + 8*CF*CA*
     &    opx**(-1)*lnx*Li2mx + 8*CF*CA*opx**(-1)*lnx**2*lnmuonpt - 4*
     &    CF*CA*opx**(-1)*lnx**2*lnopx + 2.D0/3.D0*CF*CA*opx**(-1)*
     &    lnx**3 + 16*CF*CA*opx**(-1)*Li3xonopx - 16*CF*CA*opx**(-1)*
     &    Li3x + 8*CF*CA*opx**(-1)*Li2x*lnx + 8*CF*CA*Li2mx - 14*CF*CA*
     &    lnx + 16*CF*CA*lnx*lnmuonpt + 8*CF*CA*lnx*lnopx - 15*CF*CA*
     &    omx + 4*CF*CA*omx*Li3mx + 16*CF*CA*omx*lnmuonpt + 4.D0/3.D0*
     &    CF*CA*omx*lnopx**3 - 4*CF*CA*omx*Li2mx + 16*CF*CA*omx*Li2mx*
     &    lnmuonpt + 11*CF*CA*omx*lnx - 8*CF*CA*omx*lnx*lnmuonpt - 4*CF
     &    *CA*omx*lnx*lnopx + 16*CF*CA*omx*lnx*lnopx*lnmuonpt - 4*CF*CA
     &    *omx*lnx*lnopx**2
      I2P(qbq,rglr) = I2P(qbq,rglr) - 4*CF*CA*omx*lnx*Li2mx - 4*CF*CA*
     &    omx*lnx*lnomx - 4*CF*CA*omx*lnx**2*lnmuonpt + 2*CF*CA*omx*
     &    lnx**2*lnopx - 1.D0/3.D0*CF*CA*omx*lnx**3 - 8*CF*CA*omx*
     &    Li3xonopx + 8*CF*CA*omx*Li3x - 4*CF*CA*omx*Li2x - 4*CF*CA*omx
     &    *Li2x*lnx - 8.D0/3.D0*CF*CA*pisq*opx**(-1)*lnmuonpt + 4.D0/3.D
     &    0*CF*CA*pisq*opx**(-1)*lnopx + 2.D0/3.D0*CF*CA*pisq + 1.D0/3.D
     &    0*CF*CA*pisq*omx + 4.D0/3.D0*CF*CA*pisq*omx*lnmuonpt - 2.D0/3.
     &    D0*CF*CA*pisq*omx*lnopx - 4*CF*CA*zeta3*opx**(-1) + 2*CF*CA*
     &    zeta3*omx + 16*CF**2*opx**(-1)*Li3mx + 16.D0/3.D0*CF**2*
     &    opx**(-1)*lnopx**3 + 64*CF**2*opx**(-1)*Li2mx*lnmuonpt + 64*
     &    CF**2*opx**(-1)*lnx*lnopx*lnmuonpt - 16*CF**2*opx**(-1)*lnx*
     &    lnopx**2 - 16*CF**2*opx**(-1)*lnx*Li2mx - 16*CF**2*opx**(-1)*
     &    lnx**2*lnmuonpt + 8*CF**2*opx**(-1)*lnx**2*lnopx - 4.D0/3.D0*
     &    CF**2*opx**(-1)*lnx**3 - 32*CF**2*opx**(-1)*Li3xonopx + 32*
     &    CF**2*opx**(-1)*Li3x - 16*CF**2*opx**(-1)*Li2x*lnx - 16*CF**2
     &    *Li2mx
      I2P(qbq,rglr) = I2P(qbq,rglr) + 28*CF**2*lnx - 32*CF**2*lnx*
     &    lnmuonpt - 16*CF**2*lnx*lnopx + 30*CF**2*omx - 8*CF**2*omx*
     &    Li3mx - 32*CF**2*omx*lnmuonpt - 8.D0/3.D0*CF**2*omx*lnopx**3
     &     + 8*CF**2*omx*Li2mx - 32*CF**2*omx*Li2mx*lnmuonpt - 22*CF**2
     &    *omx*lnx + 16*CF**2*omx*lnx*lnmuonpt + 8*CF**2*omx*lnx*lnopx
     &     - 32*CF**2*omx*lnx*lnopx*lnmuonpt + 8*CF**2*omx*lnx*lnopx**2
     &     + 8*CF**2*omx*lnx*Li2mx + 8*CF**2*omx*lnx*lnomx + 8*CF**2*
     &    omx*lnx**2*lnmuonpt - 4*CF**2*omx*lnx**2*lnopx + 2.D0/3.D0*
     &    CF**2*omx*lnx**3 + 16*CF**2*omx*Li3xonopx - 16*CF**2*omx*Li3x
     &     + 8*CF**2*omx*Li2x + 8*CF**2*omx*Li2x*lnx + 16.D0/3.D0*CF**2
     &    *pisq*opx**(-1)*lnmuonpt - 8.D0/3.D0*CF**2*pisq*opx**(-1)*
     &    lnopx - 4.D0/3.D0*CF**2*pisq - 2.D0/3.D0*CF**2*pisq*omx - 8.D0
     &    /3.D0*CF**2*pisq*omx*lnmuonpt + 4.D0/3.D0*CF**2*pisq*omx*
     &    lnopx + 8*CF**2*zeta3*opx**(-1) - 4*CF**2*zeta3*omx

