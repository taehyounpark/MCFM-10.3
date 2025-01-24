!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      if (x < xmin) then
      I2R(qpq,rglr)= + R8 * (  - 40*CF*omx**(-7)*Li2mx + 74.D0/3.D0*CF*
     &    omx**(-7)*lnx + 40*CF*omx**(-7)*Li2omx - 46.D0/3.D0*CF*
     &    omx**(-6) + 20*CF*omx**(-6)*lnopxon2 + 240*CF*omx**(-6)*Li2mx
     &     - 148*CF*omx**(-6)*lnx - 240*CF*omx**(-6)*Li2omx + 283.D0/3.D
     &    0*CF*omx**(-5) - 104*CF*omx**(-5)*lnopxon2 - 1928.D0/3.D0*CF*
     &    omx**(-5)*Li2mx + 3568.D0/9.D0*CF*omx**(-5)*lnx + 1928.D0/3.D0
     &    *CF*omx**(-5)*Li2omx - 4577.D0/18.D0*CF*omx**(-4) + 238*CF*
     &    omx**(-4)*lnopxon2 + 3034.D0/3.D0*CF*omx**(-4)*Li2mx - 11239.D
     &    0/18.D0*CF*omx**(-4)*lnx - 3034.D0/3.D0*CF*omx**(-4)*Li2omx
     &     + 395*CF*omx**(-3) - 944.D0/3.D0*CF*omx**(-3)*lnopxon2 - 
     &    9295.D0/9.D0*CF*omx**(-3)*Li2mx + 1919.D0/3.D0*CF*omx**(-3)*
     &    lnx + 9295.D0/9.D0*CF*omx**(-3)*Li2omx - 167839.D0/432.D0*CF*
     &    omx**(-2) + 791.D0/3.D0*CF*omx**(-2)*lnopxon2 + 2129.D0/3.D0*
     &    CF*omx**(-2)*Li2mx - 3983.D0/9.D0*CF*omx**(-2)*lnx - 2129.D0/
     &    3.D0*CF*omx**(-2)*Li2omx + 107741.D0/432.D0*CF*omx**(-1) - 
     &    863.D0/6.D0*CF*omx**(-1)*lnopxon2 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R8 * (  - 5875.D0/18.D0*CF*
     &    omx**(-1)*Li2mx + 3745.D0/18.D0*CF*omx**(-1)*lnx + 5875.D0/18.
     &    D0*CF*omx**(-1)*Li2omx + 5.D0/24.D0*CF*opx**(-6) - 13.D0/32.D0
     &    *CF*opx**(-5) + 55.D0/192.D0*CF*opx**(-4) - 83.D0/1152.D0*CF*
     &    opx**(-3) - 127.D0/2304.D0*CF*opx**(-2) - 491.D0/2304.D0*CF*
     &    opx**(-1) - 5103029123.D0/50803200.D0*CF + 157.D0/3.D0*CF*
     &    lnopxon2 + 94*CF*Li2mx - 610.D0/9.D0*CF*lnx - 94*CF*Li2omx + 
     &    1764498929.D0/90316800.D0*CF*omx - 127.D0/9.D0*CF*omx*
     &    lnopxon2 - 235.D0/18.D0*CF*omx*Li2mx + 33.D0/2.D0*CF*omx*lnx
     &     + 235.D0/18.D0*CF*omx*Li2omx + 699712933.D0/812851200.D0*CF*
     &    omx**2 + 47.D0/18.D0*CF*omx**2*lnopxon2 - 1.D0/9.D0*CF*omx**2
     &    *Li2mx - 47.D0/18.D0*CF*omx**2*lnx + 1.D0/9.D0*CF*omx**2*
     &    Li2omx - 131597267.D0/812851200.D0*CF*omx**3 - 20521273.D0/
     &    203212800.D0*CF*omx**4 - 2982737.D0/67737600.D0*CF*omx**5 - 
     &    10.D0/3.D0*CF*pisq*omx**(-7) + 20*CF*pisq*omx**(-6) - 482.D0/
     &    9.D0*CF*pisq*omx**(-5) )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R8 * ( 1517.D0/18.D0*CF*pisq*
     &    omx**(-4) - 9295.D0/108.D0*CF*pisq*omx**(-3) + 2129.D0/36.D0*
     &    CF*pisq*omx**(-2) - 5875.D0/216.D0*CF*pisq*omx**(-1) + 47.D0/
     &    6.D0*CF*pisq - 235.D0/216.D0*CF*pisq*omx - 1.D0/108.D0*CF*
     &    pisq*omx**2 + 40*CF*ln2*omx**(-6) - 220*CF*ln2*omx**(-5) + 
     &    536*CF*ln2*omx**(-4) - 760*CF*ln2*omx**(-3) + 2068.D0/3.D0*CF
     &    *ln2*omx**(-2) - 7399.D0/18.D0*CF*ln2*omx**(-1) + 297163.D0/
     &    1890.D0*CF*ln2 - 254069.D0/7560.D0*CF*ln2*omx + 12763.D0/7560.
     &    D0*CF*ln2*omx**2 + 487.D0/1890.D0*CF*ln2*omx**3 + 37.D0/315.D0
     &    *CF*ln2*omx**4 + 17.D0/420.D0*CF*ln2*omx**5 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R6 * (  - 16*CF*omx**(-5)*Li2mx
     &     + 28.D0/3.D0*CF*omx**(-5)*lnx + 16*CF*omx**(-5)*Li2omx - 20.D
     &    0/3.D0*CF*omx**(-4) + 8*CF*omx**(-4)*lnopxon2 + 72*CF*
     &    omx**(-4)*Li2mx - 42*CF*omx**(-4)*lnx - 72*CF*omx**(-4)*
     &    Li2omx + 92.D0/3.D0*CF*omx**(-3) - 88.D0/3.D0*CF*omx**(-3)*
     &    lnopxon2 - 140*CF*omx**(-3)*Li2mx + 82*CF*omx**(-3)*lnx + 140
     &    *CF*omx**(-3)*Li2omx - 176.D0/3.D0*CF*omx**(-2) + 136.D0/3.D0
     &    *CF*omx**(-2)*lnopxon2 + 460.D0/3.D0*CF*omx**(-2)*Li2mx - 88*
     &    CF*omx**(-2)*lnx - 460.D0/3.D0*CF*omx**(-2)*Li2omx + 190.D0/3.
     &    D0*CF*omx**(-1) - 116.D0/3.D0*CF*omx**(-1)*lnopxon2 - 302.D0/
     &    3.D0*CF*omx**(-1)*Li2mx + 190.D0/3.D0*CF*omx**(-1)*lnx + 302.D
     &    0/3.D0*CF*omx**(-1)*Li2omx + 1.D0/2.D0*CF*opx**(-4) - 7.D0/24.
     &    D0*CF*opx**(-3) + 7.D0/24.D0*CF*opx**(-2) + 5.D0/12.D0*CF*
     &    opx**(-1) - 461551.D0/14400.D0*CF + 80.D0/3.D0*CF*lnopxon2 + 
     &    32*CF*Li2mx - 124.D0/3.D0*CF*lnx - 32*CF*Li2omx - 9244189.D0/
     &    1411200.D0*CF*omx )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R6 * (  - 50.D0/3.D0*CF*omx*
     &    lnopxon2 + 10.D0/3.D0*CF*omx*Li2mx + 64.D0/3.D0*CF*omx*lnx - 
     &    10.D0/3.D0*CF*omx*Li2omx + 1618037.D0/156800.D0*CF*omx**2 + 
     &    14.D0/3.D0*CF*omx**2*lnopxon2 - 4*CF*omx**2*Li2mx - 14.D0/3.D0
     &    *CF*omx**2*lnx + 4*CF*omx**2*Li2omx - 175019.D0/151200.D0*CF*
     &    omx**3 - 6712.D0/33075.D0*CF*omx**4 - 4.D0/3.D0*CF*pisq*
     &    omx**(-5) + 6*CF*pisq*omx**(-4) - 35.D0/3.D0*CF*pisq*
     &    omx**(-3) + 115.D0/9.D0*CF*pisq*omx**(-2) - 151.D0/18.D0*CF*
     &    pisq*omx**(-1) + 8.D0/3.D0*CF*pisq + 5.D0/18.D0*CF*pisq*omx
     &     - 1.D0/3.D0*CF*pisq*omx**2 + 16*CF*ln2*omx**(-4) - 64*CF*ln2
     &    *omx**(-3) + 328.D0/3.D0*CF*ln2*omx**(-2) - 310.D0/3.D0*CF*
     &    ln2*omx**(-1) + 838.D0/15.D0*CF*ln2 - 1999.D0/210.D0*CF*ln2*
     &    omx - 409.D0/70.D0*CF*ln2*omx**2 + 62.D0/45.D0*CF*ln2*omx**3
     &     + 62.D0/315.D0*CF*ln2*omx**4 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R4 * (  - 8*CF*omx**(-3)*Li2mx + 
     &    4*CF*omx**(-3)*lnx + 8*CF*omx**(-3)*Li2omx - 4*CF*omx**(-2)
     &     + 4*CF*omx**(-2)*lnopxon2 + 24*CF*omx**(-2)*Li2mx - 12*CF*
     &    omx**(-2)*lnx - 24*CF*omx**(-2)*Li2omx + 12*CF*omx**(-1) - 8*
     &    CF*omx**(-1)*lnopxon2 - 28*CF*omx**(-1)*Li2mx + 28*CF*
     &    omx**(-1)*lnx + 28*CF*omx**(-1)*Li2omx + 2*CF*opx**(-2) + 9.D0
     &    /2.D0*CF*opx**(-1) - 53.D0/36.D0*CF + 16*CF*lnopxon2 - 36*CF*
     &    lnx - 199111.D0/7200.D0*CF*omx - 12*CF*omx*lnopxon2 + 20*CF*
     &    omx*Li2mx + 16*CF*omx*lnx - 20*CF*omx*Li2omx + 105341.D0/7200.
     &    D0*CF*omx**2 - 8*CF*omx**2*Li2mx + 8*CF*omx**2*Li2omx - 7141.D
     &    0/7200.D0*CF*omx**3 - 2.D0/3.D0*CF*pisq*omx**(-3) + 2*CF*pisq
     &    *omx**(-2) - 7.D0/3.D0*CF*pisq*omx**(-1) + 5.D0/3.D0*CF*pisq*
     &    omx - 2.D0/3.D0*CF*pisq*omx**2 + 8*CF*ln2*omx**(-2) - 20*CF*
     &    ln2*omx**(-1) + 56.D0/3.D0*CF*ln2 + 121.D0/15.D0*CF*ln2*omx
     &     - 221.D0/15.D0*CF*ln2*omx**2 + 16.D0/15.D0*CF*ln2*omx**3 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R2 * (  - 8*CF*omx**(-1)*Li2mx + 
     &    8*CF*omx**(-1)*Li2omx - 8*CF + 8*CF*lnx + 49.D0/18.D0*CF*omx
     &     + 8*CF*omx*lnopxon2 + 8*CF*omx*Li2mx - 16*CF*omx*lnx - 8*CF*
     &    omx*Li2omx - 8*CF*omx**2*lnopxon2 + 8*CF*omx**2*lnx - 2.D0/3.D
     &    0*CF*pisq*omx**(-1) + 2.D0/3.D0*CF*pisq*omx + 8*CF*ln2 + 16.D0
     &    /3.D0*CF*ln2*omx - 8*CF*ln2*omx**2 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + picube * (  - 4*CF*omx**(-1)*
     &    BCQQBARSx + 4*CF*omx**(-1)*BCQQBARS1 + 4*CF*BCQQBARSx )
      else
      I2R(qpq,rglr)= + R8 * (  - 44803.D0/38707200.D0*CF*omx - 3563.D0/
     &    16588800.D0*CF*omx**2 + 16515817.D0/406425600.D0*CF*omx**3 - 
     &    3774762301.D0/84304281600.D0*CF*omx**4 - 25999035719.D0/
     &    1180259942400.D0*CF*omx**5 + 12948910797859.D0/
     &    1196783581593600.D0*CF*omx**6 + 1262074977101.D0/
     &    199463930265600.D0*CF*omx**7 + 65048978167.D0/15542643916800.D
     &    0*CF*omx**8 + 462519857.D0/155831195520.D0*CF*omx**9 - 
     &    8562227.D0/26127360.D0*CF*omx**10 + 15649297.D0/52254720.D0*
     &    CF*omx**11 + 1.D0/1080.D0*CF*ln2*omx - 1.D0/1080.D0*CF*ln2*
     &    omx**2 - 241.D0/7560.D0*CF*ln2*omx**3 + 257.D0/5940.D0*CF*ln2
     &    *omx**4 + 569.D0/41580.D0*CF*ln2*omx**5 - 38243.D0/3243240.D0
     &    *CF*ln2*omx**6 - 12713.D0/2162160.D0*CF*ln2*omx**7 - 67.D0/
     &    21060.D0*CF*ln2*omx**8 - 785.D0/432432.D0*CF*ln2*omx**9 - 949.
     &    D0/648.D0*CF*ln2*omx**10 - 1.D0/81.D0*CF*ln2*omx**11 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R6 * ( 473.D0/52920.D0*CF*omx + 
     &    155089.D0/4233600.D0*CF*omx**2 - 1604599.D0/8467200.D0*CF*
     &    omx**3 + 649541.D0/20321280.D0*CF*omx**4 + 5964449.D0/
     &    67737600.D0*CF*omx**5 + 953588653.D0/24588748800.D0*CF*omx**6
     &     + 1271522809.D0/73766246400.D0*CF*omx**7 + 126487248037.D0/
     &    19946393026560.D0*CF*omx**8 + 31162391441.D0/66487976755200.D0
     &    *CF*omx**9 - 16279.D0/9720.D0*CF*omx**10 + 1269757.D0/1451520.
     &    D0*CF*omx**11 - 1.D0/126.D0*CF*ln2*omx - 11.D0/630.D0*CF*ln2*
     &    omx**2 + 103.D0/630.D0*CF*ln2*omx**3 - 1.D0/21.D0*CF*ln2*
     &    omx**4 - 23.D0/315.D0*CF*ln2*omx**5 - 169.D0/6930.D0*CF*ln2*
     &    omx**6 - 299.D0/41580.D0*CF*ln2*omx**7 - 10.D0/27027.D0*CF*
     &    ln2*omx**8 + 443.D0/180180.D0*CF*ln2*omx**9 - 7.D0/54.D0*CF*
     &    ln2*omx**10 - 4.D0/9.D0*CF*ln2*omx**11 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R4 * (  - 211.D0/7200.D0*CF*omx
     &     - 631.D0/1440.D0*CF*omx**2 + 1901.D0/3600.D0*CF*omx**3 + 
     &    117097.D0/705600.D0*CF*omx**4 - 17903.D0/282240.D0*CF*omx**5
     &     - 32327.D0/290304.D0*CF*omx**6 - 19541.D0/172800.D0*CF*
     &    omx**7 - 1251729527.D0/12294374400.D0*CF*omx**8 - 719907829.D0
     &    /8196249600.D0*CF*omx**9 - 3983363.D0/1451520.D0*CF*omx**10
     &     + 518587.D0/725760.D0*CF*omx**11 + 1.D0/15.D0*CF*ln2*omx + 1.
     &    D0/3.D0*CF*ln2*omx**2 - 7.D0/15.D0*CF*ln2*omx**3 - 8.D0/105.D0
     &    *CF*ln2*omx**4 + 2.D0/21.D0*CF*ln2*omx**5 + 1.D0/9.D0*CF*ln2*
     &    omx**6 + 1.D0/10.D0*CF*ln2*omx**7 + 296.D0/3465.D0*CF*ln2*
     &    omx**8 + 167.D0/2310.D0*CF*ln2*omx**9 + 11.D0/9.D0*CF*ln2*
     &    omx**10 - 8.D0/9.D0*CF*ln2*omx**11 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + R2 * (  - 23.D0/18.D0*CF*omx + 23.
     &    D0/9.D0*CF*omx**2 - 11.D0/6.D0*CF*omx**3 - 2059.D0/1800.D0*CF
     &    *omx**4 - 541.D0/720.D0*CF*omx**5 - 30391.D0/58800.D0*CF*
     &    omx**6 - 373.D0/1008.D0*CF*omx**7 - 1393619.D0/5080320.D0*CF*
     &    omx**8 - 169031.D0/806400.D0*CF*omx**9 + 11867.D0/181440.D0*
     &    CF*omx**10 - 511.D0/576.D0*CF*omx**11 + 4.D0/3.D0*CF*ln2*omx
     &     - 8.D0/3.D0*CF*ln2*omx**2 + 2*CF*ln2*omx**3 + 16.D0/15.D0*CF
     &    *ln2*omx**4 + 2.D0/3.D0*CF*ln2*omx**5 + 16.D0/35.D0*CF*ln2*
     &    omx**6 + 1.D0/3.D0*CF*ln2*omx**7 + 16.D0/63.D0*CF*ln2*omx**8
     &     + 1.D0/5.D0*CF*ln2*omx**9 + 8.D0/9.D0*CF*ln2*omx**10 )
      I2R(qpq,rglr) = I2R(qpq,rglr) + picube * (  - 4*CF*omx**(-1)*
     &    BCQQBARSx + 4*CF*omx**(-1)*BCQQBARS1 + 4*CF*BCQQBARSx )
      endif
