!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2R(qq,delt)= + R10 * (  - 149.D0/3870720.D0*CF**2 )
      I2R(qq,delt) = I2R(qq,delt) + R8 * (  - 30571386962947.D0/
     &    327168056033280000.D0*CF*nfl + 2692634479.D0/35407798272000.D0
     &    *CF*ln2*nfl - 4001.D0/59875200.D0*CF*ln2**2*nfl - 
     &    153495372888217.D0/163584028016640000.D0*CF*CA + 1.D0/518400.D
     &    0*CF*CA*zeta3 + 48471314251.D0/70815596544000.D0*CF*CA*ln2 - 
     &    28529.D0/119750400.D0*CF*CA*ln2**2 - 49.D0/69120.D0*CF**2 )
      I2R(qq,delt) = I2R(qq,delt) + R6 * ( 651413779.D0/3840721920000.D0
     &    *CF*nfl - 1272767.D0/3048192000.D0*CF*ln2*nfl + 11.D0/12600.D0
     &    *CF*ln2**2*nfl + 100116883541.D0/3840721920000.D0*CF*CA - 
     &    25763593.D0/3048192000.D0*CF*CA*ln2 + 23.D0/4200.D0*CF*CA*
     &    ln2**2 - 11.D0/144.D0*CF**2 )
      I2R(qq,delt) = I2R(qq,delt) + R4 * ( 13532213.D0/5334336000.D0*CF
     &    *nfl - 57311.D0/6350400.D0*CF*ln2*nfl - 53.D0/3780.D0*CF*
     &    ln2**2*nfl - 1108793639.D0/5334336000.D0*CF*CA + 1.D0/72.D0*
     &    CF*CA*zeta3 + 1706191.D0/12700800.D0*CF*CA*ln2 - 587.D0/7560.D
     &    0*CF*CA*ln2**2 + 2*CF**2 - 2*CF**2*ln2 )
      I2R(qq,delt) = I2R(qq,delt) + R2 * (  - 13033.D0/81000.D0*CF*nfl
     &     + 881.D0/5400.D0*CF*ln2*nfl + 7.D0/45.D0*CF*ln2**2*nfl + 
     &    368533.D0/81000.D0*CF*CA - CF*CA*zeta3 - 23381.D0/5400.D0*CF*
     &    CA*ln2 - 52.D0/45.D0*CF*CA*ln2**2 + 4*CF**2*zeta3 )
      I2R(qq,delt) = I2R(qq,delt) + lnR * ( 326.D0/27.D0*CF*nfl - 116.D0
     &    /9.D0*CF*ln2*nfl - 16.D0/3.D0*CF*ln2**2*nfl - 1622.D0/27.D0*
     &    CF*CA + 8*CF*CA*zeta3 + 548.D0/9.D0*CF*CA*ln2 + 88.D0/3.D0*CF
     &    *CA*ln2**2 )
      I2R(qq,delt) = I2R(qq,delt) - 325.D0/9.D0*CF*nfl + 29.D0/54.D0*CF
     &    *pisq*nfl + 640.D0/27.D0*CF*ln2*nfl + 4.D0/9.D0*CF*ln2*pisq*
     &    nfl + 20*CF*ln2**2*nfl + 8.D0/9.D0*CF*ln2**2*ln64*nfl + 16.D0/
     &    9.D0*CF*ln2**3*nfl + 1621.D0/9.D0*CF*CA + 128*CF*CA*Li4xhalf
     &     - 137.D0/54.D0*CF*CA*pisq - 47.D0/30.D0*CF*CA*pisq**2 - 3232.
     &    D0/27.D0*CF*CA*ln2 - 22.D0/9.D0*CF*CA*ln2*pisq + 112*CF*CA*
     &    ln2*zeta3 - 92*CF*CA*ln2**2 - 16.D0/3.D0*CF*CA*ln2**2*pisq - 
     &    352.D0/9.D0*CF*CA*ln2**3 + 16.D0/3.D0*CF*CA*ln2**4 + 32*CF**2
     &    *zeta3*lnmuonpt

      I2R(qq,plus)= + R8 * (  - 5664846191.D0/35407798272000.D0*CF*nfl
     &     + 4001.D0/29937600.D0*CF*ln2*nfl - 50937246539.D0/
     &    70815596544000.D0*CF*CA - 1.D0/777600.D0*CF*CA*pisq + 28529.D0
     &    /59875200.D0*CF*CA*ln2 )
      I2R(qq,plus) = I2R(qq,plus) + R6 * ( 7001023.D0/3048192000.D0*CF*
     &    nfl - 11.D0/6300.D0*CF*ln2*nfl + 74801417.D0/3048192000.D0*CF
     &    *CA - 23.D0/2100.D0*CF*CA*ln2 )
      I2R(qq,plus) = I2R(qq,plus) + R4 * (  - 168401.D0/6350400.D0*CF*
     &    nfl + 53.D0/1890.D0*CF*ln2*nfl - 9383279.D0/12700800.D0*CF*CA
     &     - 1.D0/108.D0*CF*CA*pisq + 587.D0/3780.D0*CF*CA*ln2 + 2*
     &    CF**2 )
      I2R(qq,plus) = I2R(qq,plus) + R2 * ( 3071.D0/5400.D0*CF*nfl - 14.D
     &    0/45.D0*CF*ln2*nfl + 1429.D0/5400.D0*CF*CA + 2.D0/3.D0*CF*CA*
     &    pisq + 104.D0/45.D0*CF*CA*ln2 - 8.D0/3.D0*CF**2*pisq )
      I2R(qq,plus) = I2R(qq,plus) + picube * (  - 4*CF*BCQQBARS1 - 4*CF
     &    *nfl*BCQQCFTF1 - 8*CF*CA*BCQQCACF1 - 8*CF**2*BCQQCFCF1 )
      I2R(qq,plus) = I2R(qq,plus) + lnR * (  - 92.D0/9.D0*CF*nfl + 32.D0
     &    /3.D0*CF*ln2*nfl + 524.D0/9.D0*CF*CA - 16.D0/3.D0*CF*CA*pisq
     &     - 176.D0/3.D0*CF*CA*ln2 )
      I2R(qq,plus) = I2R(qq,plus) - 26.D0/9.D0*CF*nfl + 92.D0/9.D0*CF*
     &    lnpi*nfl + 116.D0/9.D0*CF*ln2*nfl - 32.D0/3.D0*CF*ln2*lnpi*
     &    nfl - 32.D0/3.D0*CF*ln2**2*nfl + 26.D0/9.D0*CF*CA - 524.D0/9.D
     &    0*CF*CA*lnpi + 16.D0/3.D0*CF*CA*pisq*lnpi - 548.D0/9.D0*CF*CA
     &    *ln2 + 176.D0/3.D0*CF*CA*ln2*lnpi + 16.D0/3.D0*CF*CA*ln2*pisq
     &     + 176.D0/3.D0*CF*CA*ln2**2

      if (x < xmin) then
      I2R(qq,rglr)= + R8 * (  - 40*CF*omx**(-7)*Li2mx + 74.D0/3.D0*CF*
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
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 5875.D0/18.D0*CF*
     &    omx**(-1)*Li2mx + 3745.D0/18.D0*CF*omx**(-1)*lnx + 5875.D0/18.
     &    D0*CF*omx**(-1)*Li2omx + 5.D0/24.D0*CF*opx**(-6) - 13.D0/32.D0
     &    *CF*opx**(-5) + 55.D0/192.D0*CF*opx**(-4) - 83.D0/1152.D0*CF*
     &    opx**(-3) - 127.D0/2304.D0*CF*opx**(-2) - 491.D0/2304.D0*CF*
     &    opx**(-1) - 5103029123.D0/50803200.D0*CF + 37432579841.D0/
     &    35407798272000.D0*CF*nfl + 157.D0/3.D0*CF*lnopxon2 + 94*CF*
     &    Li2mx - 610.D0/9.D0*CF*lnx - 94*CF*Li2omx + 1764498929.D0/
     &    90316800.D0*CF*omx - 219702602291.D0/70815596544000.D0*CF*omx
     &    *nfl - 127.D0/9.D0*CF*omx*lnopxon2 - 235.D0/18.D0*CF*omx*
     &    Li2mx + 33.D0/2.D0*CF*omx*lnx + 235.D0/18.D0*CF*omx*Li2omx + 
     &    699712933.D0/812851200.D0*CF*omx**2 + 768879803.D0/
     &    98354995200.D0*CF*omx**2*nfl + 47.D0/18.D0*CF*omx**2*lnopxon2
     &     - 1.D0/9.D0*CF*omx**2*Li2mx - 47.D0/18.D0*CF*omx**2*lnx + 1.D
     &    0/9.D0*CF*omx**2*Li2omx - 131597267.D0/812851200.D0*CF*omx**3
     &     - 4159161163.D0/295064985600.D0*CF*omx**3*nfl )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 20521273.D0/203212800.D0*
     &    CF*omx**4 + 83202043.D0/6147187200.D0*CF*omx**4*nfl - 2982737.
     &    D0/67737600.D0*CF*omx**5 - 2062751.D0/409812480.D0*CF*omx**5*
     &    nfl - 10.D0/3.D0*CF*pisq*omx**(-7) + 20*CF*pisq*omx**(-6) - 
     &    482.D0/9.D0*CF*pisq*omx**(-5) + 1517.D0/18.D0*CF*pisq*
     &    omx**(-4) - 9295.D0/108.D0*CF*pisq*omx**(-3) + 2129.D0/36.D0*
     &    CF*pisq*omx**(-2) - 5875.D0/216.D0*CF*pisq*omx**(-1) + 47.D0/
     &    6.D0*CF*pisq - 235.D0/216.D0*CF*pisq*omx - 1.D0/108.D0*CF*
     &    pisq*omx**2 + 40*CF*ln2*omx**(-6) - 220*CF*ln2*omx**(-5) + 
     &    536*CF*ln2*omx**(-4) - 760*CF*ln2*omx**(-3) + 2068.D0/3.D0*CF
     &    *ln2*omx**(-2) - 7399.D0/18.D0*CF*ln2*omx**(-1) + 297163.D0/
     &    1890.D0*CF*ln2 - 24821.D0/29937600.D0*CF*ln2*nfl - 254069.D0/
     &    7560.D0*CF*ln2*omx + 139241.D0/59875200.D0*CF*ln2*omx*nfl + 
     &    12763.D0/7560.D0*CF*ln2*omx**2 - 241.D0/41580.D0*CF*ln2*
     &    omx**2*nfl + 487.D0/1890.D0*CF*ln2*omx**3 + 113.D0/10395.D0*
     &    CF*ln2*omx**3*nfl )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * ( 37.D0/315.D0*CF*ln2*omx**4
     &     - 38.D0/3465.D0*CF*ln2*omx**4*nfl + 17.D0/420.D0*CF*ln2*
     &    omx**5 + 1.D0/231.D0*CF*ln2*omx**5*nfl + 140*CF*CA*omx**(-9)*
     &    Li2mx - 533.D0/6.D0*CF*CA*omx**(-9)*lnx - 140*CF*CA*omx**(-9)
     &    *Li2omx + 307.D0/6.D0*CF*CA*omx**(-8) - 70*CF*CA*omx**(-8)*
     &    lnopxon2 - 980*CF*CA*omx**(-8)*Li2mx + 3731.D0/6.D0*CF*CA*
     &    omx**(-8)*lnx + 980*CF*CA*omx**(-8)*Li2omx - 4411.D0/12.D0*CF
     &    *CA*omx**(-7) + 435*CF*CA*omx**(-7)*lnopxon2 + 3060*CF*CA*
     &    omx**(-7)*Li2mx - 11647.D0/6.D0*CF*CA*omx**(-7)*lnx - 3060*CF
     &    *CA*omx**(-7)*Li2omx + 20945.D0/18.D0*CF*CA*omx**(-6) - 3565.D
     &    0/3.D0*CF*CA*omx**(-6)*lnopxon2 - 16840.D0/3.D0*CF*CA*
     &    omx**(-6)*Li2mx + 64063.D0/18.D0*CF*CA*omx**(-6)*lnx + 16840.D
     &    0/3.D0*CF*CA*omx**(-6)*Li2omx - 76805.D0/36.D0*CF*CA*
     &    omx**(-5) + 11239.D0/6.D0*CF*CA*omx**(-5)*lnopxon2 + 1205999.D
     &    0/180.D0*CF*CA*omx**(-5)*Li2mx - 9168679.D0/2160.D0*CF*CA*
     &    omx**(-5)*lnx )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 1205999.D0/180.D0*CF*CA*
     &    omx**(-5)*Li2omx + 10806163.D0/4320.D0*CF*CA*omx**(-4) - 
     &    675899.D0/360.D0*CF*CA*omx**(-4)*lnopxon2 - 326179.D0/60.D0*
     &    CF*CA*omx**(-4)*Li2mx + 7433033.D0/2160.D0*CF*CA*omx**(-4)*
     &    lnx + 326179.D0/60.D0*CF*CA*omx**(-4)*Li2omx - 16843579.D0/
     &    8640.D0*CF*CA*omx**(-3) + 178631.D0/144.D0*CF*CA*omx**(-3)*
     &    lnopxon2 + 1638799.D0/540.D0*CF*CA*omx**(-3)*Li2mx - 8295319.D
     &    0/4320.D0*CF*CA*omx**(-3)*lnx - 1638799.D0/540.D0*CF*CA*
     &    omx**(-3)*Li2omx + 26257819.D0/25920.D0*CF*CA*omx**(-2) - 
     &    77735.D0/144.D0*CF*CA*omx**(-2)*lnopxon2 - 620363.D0/540.D0*
     &    CF*CA*omx**(-2)*Li2mx + 262039.D0/360.D0*CF*CA*omx**(-2)*lnx
     &     + 620363.D0/540.D0*CF*CA*omx**(-2)*Li2omx - 17621993.D0/
     &    51840.D0*CF*CA*omx**(-1) + 42911.D0/288.D0*CF*CA*omx**(-1)*
     &    lnopxon2 + 67399.D0/240.D0*CF*CA*omx**(-1)*Li2mx - 258071.D0/
     &    1440.D0*CF*CA*omx**(-1)*lnx - 67399.D0/240.D0*CF*CA*omx**(-1)
     &    *Li2omx )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 5.D0/768.D0*CF*CA*
     &    opx**(-6) + 3.D0/256.D0*CF*CA*opx**(-5) - 13.D0/1536.D0*CF*CA
     &    *opx**(-4) + 11.D0/3072.D0*CF*CA*opx**(-3) - 1117.D0/276480.D0
     &    *CF*CA*opx**(-2) - 701.D0/138240.D0*CF*CA*opx**(-1) + 
     &    136322637038569.D0/2023302758400.D0*CF*CA - 34079.D0/1440.D0*
     &    CF*CA*lnopxon2 - 85319.D0/2160.D0*CF*CA*Li2mx + 27929.D0/1080.
     &    D0*CF*CA*lnx + 85319.D0/2160.D0*CF*CA*Li2omx - 
     &    126103337593253.D0/20233027584000.D0*CF*CA*omx + 5.D0/3.D0*CF
     &    *CA*omx*lnopxon2 + 7.D0/3.D0*CF*CA*omx*Li2mx - 5.D0/3.D0*CF*
     &    CA*omx*lnx - 7.D0/3.D0*CF*CA*omx*Li2omx + 1267757111.D0/
     &    10928332800.D0*CF*CA*omx**2 - 37547707.D0/3688312320.D0*CF*CA
     &    *omx**3 + 60744253.D0/12294374400.D0*CF*CA*omx**4 - 1929101.D0
     &    /234178560.D0*CF*CA*omx**5 + 35.D0/3.D0*CF*CA*pisq*omx**(-9)
     &     - 245.D0/3.D0*CF*CA*pisq*omx**(-8) + 255*CF*CA*pisq*
     &    omx**(-7) - 4210.D0/9.D0*CF*CA*pisq*omx**(-6) + 1205999.D0/
     &    2160.D0*CF*CA*pisq*omx**(-5) )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 326179.D0/720.D0*CF*CA*
     &    pisq*omx**(-4) + 1638799.D0/6480.D0*CF*CA*pisq*omx**(-3) - 
     &    620363.D0/6480.D0*CF*CA*pisq*omx**(-2) + 67399.D0/2880.D0*CF*
     &    CA*pisq*omx**(-1) - 2559569.D0/777600.D0*CF*CA*pisq + 302399.D
     &    0/1555200.D0*CF*CA*pisq*omx - 140*CF*CA*ln2*omx**(-8) + 910*
     &    CF*CA*ln2*omx**(-7) - 7850.D0/3.D0*CF*CA*ln2*omx**(-6) + 4375
     &    *CF*CA*ln2*omx**(-5) - 845339.D0/180.D0*CF*CA*ln2*omx**(-4)
     &     + 242287.D0/72.D0*CF*CA*ln2*omx**(-3) - 116303.D0/72.D0*CF*
     &    CA*ln2*omx**(-2) + 8059.D0/16.D0*CF*CA*ln2*omx**(-1) - 
     &    800547899.D0/8553600.D0*CF*CA*ln2 + 140871551.D0/17107200.D0*
     &    CF*CA*ln2*omx - 1312.D0/10395.D0*CF*CA*ln2*omx**2 + 17.D0/
     &    2079.D0*CF*CA*ln2*omx**3 - 23.D0/6930.D0*CF*CA*ln2*omx**4 + 1.
     &    D0/132.D0*CF*CA*ln2*omx**5 - 280*CF**2*omx**(-9)*Li2mx + 533.D
     &    0/3.D0*CF**2*omx**(-9)*lnx + 280*CF**2*omx**(-9)*Li2omx - 307.
     &    D0/3.D0*CF**2*omx**(-8) + 140*CF**2*omx**(-8)*lnopxon2 + 1960
     &    *CF**2*omx**(-8)*Li2mx )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 3731.D0/3.D0*CF**2*
     &    omx**(-8)*lnx - 1960*CF**2*omx**(-8)*Li2omx + 4411.D0/6.D0*
     &    CF**2*omx**(-7) - 870*CF**2*omx**(-7)*lnopxon2 - 6120*CF**2*
     &    omx**(-7)*Li2mx + 11647.D0/3.D0*CF**2*omx**(-7)*lnx + 6120*
     &    CF**2*omx**(-7)*Li2omx - 20945.D0/9.D0*CF**2*omx**(-6) + 7130.
     &    D0/3.D0*CF**2*omx**(-6)*lnopxon2 + 33680.D0/3.D0*CF**2*
     &    omx**(-6)*Li2mx - 64063.D0/9.D0*CF**2*omx**(-6)*lnx - 33680.D0
     &    /3.D0*CF**2*omx**(-6)*Li2omx + 76805.D0/18.D0*CF**2*omx**(-5)
     &     - 11239.D0/3.D0*CF**2*omx**(-5)*lnopxon2 - 1205999.D0/90.D0*
     &    CF**2*omx**(-5)*Li2mx + 9168679.D0/1080.D0*CF**2*omx**(-5)*
     &    lnx + 1205999.D0/90.D0*CF**2*omx**(-5)*Li2omx - 10806163.D0/
     &    2160.D0*CF**2*omx**(-4) + 675899.D0/180.D0*CF**2*omx**(-4)*
     &    lnopxon2 + 326179.D0/30.D0*CF**2*omx**(-4)*Li2mx - 7433033.D0/
     &    1080.D0*CF**2*omx**(-4)*lnx - 326179.D0/30.D0*CF**2*omx**(-4)
     &    *Li2omx + 16843579.D0/4320.D0*CF**2*omx**(-3) - 178631.D0/72.D
     &    0*CF**2*omx**(-3)*lnopxon2 )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 1638799.D0/270.D0*CF**2*
     &    omx**(-3)*Li2mx + 8295319.D0/2160.D0*CF**2*omx**(-3)*lnx + 
     &    1638799.D0/270.D0*CF**2*omx**(-3)*Li2omx - 26257819.D0/12960.D
     &    0*CF**2*omx**(-2) + 77735.D0/72.D0*CF**2*omx**(-2)*lnopxon2
     &     + 620363.D0/270.D0*CF**2*omx**(-2)*Li2mx - 262039.D0/180.D0*
     &    CF**2*omx**(-2)*lnx - 620363.D0/270.D0*CF**2*omx**(-2)*Li2omx
     &     + 17621993.D0/25920.D0*CF**2*omx**(-1) - 42911.D0/144.D0*
     &    CF**2*omx**(-1)*lnopxon2 - 67399.D0/120.D0*CF**2*omx**(-1)*
     &    Li2mx + 258071.D0/720.D0*CF**2*omx**(-1)*lnx + 67399.D0/120.D0
     &    *CF**2*omx**(-1)*Li2omx + 5.D0/384.D0*CF**2*opx**(-6) - 3.D0/
     &    128.D0*CF**2*opx**(-5) + 13.D0/768.D0*CF**2*opx**(-4) - 11.D0/
     &    1536.D0*CF**2*opx**(-3) + 1117.D0/138240.D0*CF**2*opx**(-2)
     &     + 701.D0/69120.D0*CF**2*opx**(-1) - 3285928125457.D0/
     &    24385536000.D0*CF**2 + 34079.D0/720.D0*CF**2*lnopxon2 + 85319.
     &    D0/1080.D0*CF**2*Li2mx - 27929.D0/540.D0*CF**2*lnx - 85319.D0/
     &    1080.D0*CF**2*Li2omx )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * ( 5439776039.D0/435456000.D0*
     &    CF**2*omx - 10.D0/3.D0*CF**2*omx*lnopxon2 - 14.D0/3.D0*CF**2*
     &    omx*Li2mx + 10.D0/3.D0*CF**2*omx*lnx + 14.D0/3.D0*CF**2*omx*
     &    Li2omx - 10364567.D0/33868800.D0*CF**2*omx**2 + 10090177.D0/
     &    135475200.D0*CF**2*omx**3 - 207637.D0/20321280.D0*CF**2*
     &    omx**4 + 86623.D0/11289600.D0*CF**2*omx**5 - 70.D0/3.D0*CF**2
     &    *pisq*omx**(-9) + 490.D0/3.D0*CF**2*pisq*omx**(-8) - 510*
     &    CF**2*pisq*omx**(-7) + 8420.D0/9.D0*CF**2*pisq*omx**(-6) - 
     &    1205999.D0/1080.D0*CF**2*pisq*omx**(-5) + 326179.D0/360.D0*
     &    CF**2*pisq*omx**(-4) - 1638799.D0/3240.D0*CF**2*pisq*
     &    omx**(-3) + 620363.D0/3240.D0*CF**2*pisq*omx**(-2) - 67399.D0/
     &    1440.D0*CF**2*pisq*omx**(-1) + 85319.D0/12960.D0*CF**2*pisq
     &     - 7.D0/18.D0*CF**2*pisq*omx + 280*CF**2*ln2*omx**(-8) - 1820
     &    *CF**2*ln2*omx**(-7) + 15700.D0/3.D0*CF**2*ln2*omx**(-6) - 
     &    8750*CF**2*ln2*omx**(-5) + 845339.D0/90.D0*CF**2*ln2*
     &    omx**(-4) )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 242287.D0/36.D0*CF**2*ln2
     &    *omx**(-3) + 116303.D0/36.D0*CF**2*ln2*omx**(-2) - 8059.D0/8.D
     &    0*CF**2*ln2*omx**(-1) + 42452423.D0/226800.D0*CF**2*ln2 - 
     &    1068241.D0/64800.D0*CF**2*ln2*omx + 583.D0/1890.D0*CF**2*ln2*
     &    omx**2 - 473.D0/7560.D0*CF**2*ln2*omx**3 + 1.D0/126.D0*CF**2*
     &    ln2*omx**4 - 1.D0/140.D0*CF**2*ln2*omx**5 )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 16*CF*omx**(-5)*Li2mx + 
     &    28.D0/3.D0*CF*omx**(-5)*lnx + 16*CF*omx**(-5)*Li2omx - 20.D0/
     &    3.D0*CF*omx**(-4) + 8*CF*omx**(-4)*lnopxon2 + 72*CF*omx**(-4)
     &    *Li2mx - 42*CF*omx**(-4)*lnx - 72*CF*omx**(-4)*Li2omx + 92.D0/
     &    3.D0*CF*omx**(-3) - 88.D0/3.D0*CF*omx**(-3)*lnopxon2 - 140*CF
     &    *omx**(-3)*Li2mx + 82*CF*omx**(-3)*lnx + 140*CF*omx**(-3)*
     &    Li2omx - 176.D0/3.D0*CF*omx**(-2) + 136.D0/3.D0*CF*omx**(-2)*
     &    lnopxon2 + 460.D0/3.D0*CF*omx**(-2)*Li2mx - 88*CF*omx**(-2)*
     &    lnx - 460.D0/3.D0*CF*omx**(-2)*Li2omx + 190.D0/3.D0*CF*
     &    omx**(-1) - 116.D0/3.D0*CF*omx**(-1)*lnopxon2 - 302.D0/3.D0*
     &    CF*omx**(-1)*Li2mx + 190.D0/3.D0*CF*omx**(-1)*lnx + 302.D0/3.D
     &    0*CF*omx**(-1)*Li2omx + 1.D0/2.D0*CF*opx**(-4) - 7.D0/24.D0*
     &    CF*opx**(-3) + 7.D0/24.D0*CF*opx**(-2) + 5.D0/12.D0*CF*
     &    opx**(-1) - 461551.D0/14400.D0*CF - 31946443.D0/3048192000.D0
     &    *CF*nfl + 80.D0/3.D0*CF*lnopxon2 + 32*CF*Li2mx - 124.D0/3.D0*
     &    CF*lnx )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 32*CF*Li2omx - 9244189.D0/
     &    1411200.D0*CF*omx + 8250923.D0/304819200.D0*CF*omx*nfl - 50.D0
     &    /3.D0*CF*omx*lnopxon2 + 10.D0/3.D0*CF*omx*Li2mx + 64.D0/3.D0*
     &    CF*omx*lnx - 10.D0/3.D0*CF*omx*Li2omx + 1618037.D0/156800.D0*
     &    CF*omx**2 - 7688197.D0/152409600.D0*CF*omx**2*nfl + 14.D0/3.D0
     &    *CF*omx**2*lnopxon2 - 4*CF*omx**2*Li2mx - 14.D0/3.D0*CF*
     &    omx**2*lnx + 4*CF*omx**2*Li2omx - 175019.D0/151200.D0*CF*
     &    omx**3 + 2108497.D0/38102400.D0*CF*omx**3*nfl - 6712.D0/33075.
     &    D0*CF*omx**4 - 91859.D0/3810240.D0*CF*omx**4*nfl - 4.D0/3.D0*
     &    CF*pisq*omx**(-5) + 6*CF*pisq*omx**(-4) - 35.D0/3.D0*CF*pisq*
     &    omx**(-3) + 115.D0/9.D0*CF*pisq*omx**(-2) - 151.D0/18.D0*CF*
     &    pisq*omx**(-1) + 8.D0/3.D0*CF*pisq + 5.D0/18.D0*CF*pisq*omx
     &     - 1.D0/3.D0*CF*pisq*omx**2 + 16*CF*ln2*omx**(-4) - 64*CF*ln2
     &    *omx**(-3) + 328.D0/3.D0*CF*ln2*omx**(-2) - 310.D0/3.D0*CF*
     &    ln2*omx**(-1) + 838.D0/15.D0*CF*ln2 + 17.D0/2100.D0*CF*ln2*
     &    nfl )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 1999.D0/210.D0*CF*ln2*omx
     &     - 19.D0/945.D0*CF*ln2*omx*nfl - 409.D0/70.D0*CF*ln2*omx**2
     &     + 37.D0/945.D0*CF*ln2*omx**2*nfl + 62.D0/45.D0*CF*ln2*omx**3
     &     - 43.D0/945.D0*CF*ln2*omx**3*nfl + 62.D0/315.D0*CF*ln2*
     &    omx**4 + 4.D0/189.D0*CF*ln2*omx**4*nfl + 160.D0/3.D0*CF*CA*
     &    omx**(-7)*Li2mx - 296.D0/9.D0*CF*CA*omx**(-7)*lnx - 160.D0/3.D
     &    0*CF*CA*omx**(-7)*Li2omx + 184.D0/9.D0*CF*CA*omx**(-6) - 80.D0
     &    /3.D0*CF*CA*omx**(-6)*lnopxon2 - 880.D0/3.D0*CF*CA*omx**(-6)*
     &    Li2mx + 1628.D0/9.D0*CF*CA*omx**(-6)*lnx + 880.D0/3.D0*CF*CA*
     &    omx**(-6)*Li2omx - 1040.D0/9.D0*CF*CA*omx**(-5) + 376.D0/3.D0
     &    *CF*CA*omx**(-5)*lnopxon2 + 2092.D0/3.D0*CF*CA*omx**(-5)*
     &    Li2mx - 3869.D0/9.D0*CF*CA*omx**(-5)*lnx - 2092.D0/3.D0*CF*CA
     &    *omx**(-5)*Li2omx + 7469.D0/27.D0*CF*CA*omx**(-4) - 2234.D0/9.
     &    D0*CF*CA*omx**(-4)*lnopxon2 - 936*CF*CA*omx**(-4)*Li2mx + 
     &    1730.D0/3.D0*CF*CA*omx**(-4)*lnx + 936*CF*CA*omx**(-4)*Li2omx
     &     - 19687.D0/54.D0*CF*CA*omx**(-3) )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * ( 2417.D0/9.D0*CF*CA*omx**(-3)
     &    *lnopxon2 + 2326.D0/3.D0*CF*CA*omx**(-3)*Li2mx - 8597.D0/18.D0
     &    *CF*CA*omx**(-3)*lnx - 2326.D0/3.D0*CF*CA*omx**(-3)*Li2omx + 
     &    7733.D0/27.D0*CF*CA*omx**(-2) - 1541.D0/9.D0*CF*CA*omx**(-2)*
     &    lnopxon2 - 1214.D0/3.D0*CF*CA*omx**(-2)*Li2mx + 4501.D0/18.D0
     &    *CF*CA*omx**(-2)*lnx + 1214.D0/3.D0*CF*CA*omx**(-2)*Li2omx - 
     &    537.D0/4.D0*CF*CA*omx**(-1) + 1157.D0/18.D0*CF*CA*omx**(-1)*
     &    lnopxon2 + 386.D0/3.D0*CF*CA*omx**(-1)*Li2mx - 1457.D0/18.D0*
     &    CF*CA*omx**(-1)*lnx - 386.D0/3.D0*CF*CA*omx**(-1)*Li2omx - 1.D
     &    0/48.D0*CF*CA*opx**(-4) + 1.D0/144.D0*CF*CA*opx**(-3) - 11.D0/
     &    576.D0*CF*CA*opx**(-2) - 1.D0/18.D0*CF*CA*opx**(-1) + 
     &    104508748843.D0/3048192000.D0*CF*CA - 241.D0/18.D0*CF*CA*
     &    lnopxon2 - 22*CF*CA*Li2mx + 271.D0/18.D0*CF*CA*lnx + 22*CF*CA
     &    *Li2omx - 38774789.D0/12192768.D0*CF*CA*omx + 4.D0/3.D0*CF*CA
     &    *omx*lnopxon2 + 4.D0/3.D0*CF*CA*omx*Li2mx - 4.D0/3.D0*CF*CA*
     &    omx*lnx )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 4.D0/3.D0*CF*CA*omx*
     &    Li2omx - 2493689.D0/21772800.D0*CF*CA*omx**2 + 102367.D0/
     &    7620480.D0*CF*CA*omx**3 - 168493.D0/3810240.D0*CF*CA*omx**4
     &     + 40.D0/9.D0*CF*CA*pisq*omx**(-7) - 220.D0/9.D0*CF*CA*pisq*
     &    omx**(-6) + 523.D0/9.D0*CF*CA*pisq*omx**(-5) - 78*CF*CA*pisq*
     &    omx**(-4) + 1163.D0/18.D0*CF*CA*pisq*omx**(-3) - 607.D0/18.D0
     &    *CF*CA*pisq*omx**(-2) + 193.D0/18.D0*CF*CA*pisq*omx**(-1) - 
     &    11.D0/6.D0*CF*CA*pisq + 1.D0/9.D0*CF*CA*pisq*omx - 160.D0/3.D0
     &    *CF*CA*ln2*omx**(-6) + 800.D0/3.D0*CF*CA*ln2*omx**(-5) - 5116.
     &    D0/9.D0*CF*CA*ln2*omx**(-4) + 6046.D0/9.D0*CF*CA*ln2*
     &    omx**(-3) - 4294.D0/9.D0*CF*CA*ln2*omx**(-2) + 1840.D0/9.D0*
     &    CF*CA*ln2*omx**(-1) - 102667.D0/2100.D0*CF*CA*ln2 + 1801.D0/
     &    378.D0*CF*CA*ln2*omx + 14.D0/135.D0*CF*CA*ln2*omx**2 - 1.D0/
     &    189.D0*CF*CA*ln2*omx**3 + 8.D0/189.D0*CF*CA*ln2*omx**4 - 320.D
     &    0/3.D0*CF**2*omx**(-7)*Li2mx + 592.D0/9.D0*CF**2*omx**(-7)*
     &    lnx )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * ( 320.D0/3.D0*CF**2*omx**(-7)*
     &    Li2omx - 368.D0/9.D0*CF**2*omx**(-6) + 160.D0/3.D0*CF**2*
     &    omx**(-6)*lnopxon2 + 1760.D0/3.D0*CF**2*omx**(-6)*Li2mx - 
     &    3256.D0/9.D0*CF**2*omx**(-6)*lnx - 1760.D0/3.D0*CF**2*
     &    omx**(-6)*Li2omx + 2080.D0/9.D0*CF**2*omx**(-5) - 752.D0/3.D0
     &    *CF**2*omx**(-5)*lnopxon2 - 4184.D0/3.D0*CF**2*omx**(-5)*
     &    Li2mx + 7738.D0/9.D0*CF**2*omx**(-5)*lnx + 4184.D0/3.D0*CF**2
     &    *omx**(-5)*Li2omx - 14938.D0/27.D0*CF**2*omx**(-4) + 4468.D0/
     &    9.D0*CF**2*omx**(-4)*lnopxon2 + 1872*CF**2*omx**(-4)*Li2mx - 
     &    3460.D0/3.D0*CF**2*omx**(-4)*lnx - 1872*CF**2*omx**(-4)*
     &    Li2omx + 19687.D0/27.D0*CF**2*omx**(-3) - 4834.D0/9.D0*CF**2*
     &    omx**(-3)*lnopxon2 - 4652.D0/3.D0*CF**2*omx**(-3)*Li2mx + 
     &    8597.D0/9.D0*CF**2*omx**(-3)*lnx + 4652.D0/3.D0*CF**2*
     &    omx**(-3)*Li2omx - 15466.D0/27.D0*CF**2*omx**(-2) + 3082.D0/9.
     &    D0*CF**2*omx**(-2)*lnopxon2 + 2428.D0/3.D0*CF**2*omx**(-2)*
     &    Li2mx )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 4501.D0/9.D0*CF**2*
     &    omx**(-2)*lnx - 2428.D0/3.D0*CF**2*omx**(-2)*Li2omx + 537.D0/
     &    2.D0*CF**2*omx**(-1) - 1157.D0/9.D0*CF**2*omx**(-1)*lnopxon2
     &     - 772.D0/3.D0*CF**2*omx**(-1)*Li2mx + 1457.D0/9.D0*CF**2*
     &    omx**(-1)*lnx + 772.D0/3.D0*CF**2*omx**(-1)*Li2omx + 1.D0/24.D
     &    0*CF**2*opx**(-4) - 1.D0/72.D0*CF**2*opx**(-3) + 11.D0/288.D0
     &    *CF**2*opx**(-2) + 1.D0/9.D0*CF**2*opx**(-1) - 290232017.D0/
     &    4233600.D0*CF**2 + 241.D0/9.D0*CF**2*lnopxon2 + 44*CF**2*
     &    Li2mx - 271.D0/9.D0*CF**2*lnx - 44*CF**2*Li2omx + 514637.D0/
     &    84672.D0*CF**2*omx - 8.D0/3.D0*CF**2*omx*lnopxon2 - 8.D0/3.D0
     &    *CF**2*omx*Li2mx + 8.D0/3.D0*CF**2*omx*lnx + 8.D0/3.D0*CF**2*
     &    omx*Li2omx + 2225543.D0/4233600.D0*CF**2*omx**2 - 16111.D0/
     &    264600.D0*CF**2*omx**3 + 1751.D0/37800.D0*CF**2*omx**4 - 80.D0
     &    /9.D0*CF**2*pisq*omx**(-7) + 440.D0/9.D0*CF**2*pisq*omx**(-6)
     &     - 1046.D0/9.D0*CF**2*pisq*omx**(-5) + 156*CF**2*pisq*
     &    omx**(-4) )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 1163.D0/9.D0*CF**2*pisq*
     &    omx**(-3) + 607.D0/9.D0*CF**2*pisq*omx**(-2) - 193.D0/9.D0*
     &    CF**2*pisq*omx**(-1) + 11.D0/3.D0*CF**2*pisq - 2.D0/9.D0*
     &    CF**2*pisq*omx + 320.D0/3.D0*CF**2*ln2*omx**(-6) - 1600.D0/3.D
     &    0*CF**2*ln2*omx**(-5) + 10232.D0/9.D0*CF**2*ln2*omx**(-4) - 
     &    12092.D0/9.D0*CF**2*ln2*omx**(-3) + 8588.D0/9.D0*CF**2*ln2*
     &    omx**(-2) - 3680.D0/9.D0*CF**2*ln2*omx**(-1) + 30799.D0/315.D0
     &    *CF**2*ln2 - 587.D0/63.D0*CF**2*ln2*omx - 146.D0/315.D0*CF**2
     &    *ln2*omx**2 + 29.D0/630.D0*CF**2*ln2*omx**3 - 2.D0/45.D0*
     &    CF**2*ln2*omx**4 )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * (  - 8*CF*omx**(-3)*Li2mx + 4*
     &    CF*omx**(-3)*lnx + 8*CF*omx**(-3)*Li2omx - 4*CF*omx**(-2) + 4
     &    *CF*omx**(-2)*lnopxon2 + 24*CF*omx**(-2)*Li2mx - 12*CF*
     &    omx**(-2)*lnx - 24*CF*omx**(-2)*Li2omx + 12*CF*omx**(-1) - 8*
     &    CF*omx**(-1)*lnopxon2 - 28*CF*omx**(-1)*Li2mx + 28*CF*
     &    omx**(-1)*lnx + 28*CF*omx**(-1)*Li2omx + 2*CF*opx**(-2) + 9.D0
     &    /2.D0*CF*opx**(-1) - 53.D0/36.D0*CF + 371677.D0/3175200.D0*CF
     &    *nfl + 16*CF*lnopxon2 - 36*CF*lnx - 199111.D0/7200.D0*CF*omx
     &     - 2735987.D0/12700800.D0*CF*omx*nfl - 12*CF*omx*lnopxon2 + 
     &    20*CF*omx*Li2mx + 16*CF*omx*lnx - 20*CF*omx*Li2omx + 105341.D0
     &    /7200.D0*CF*omx**2 + 87469.D0/352800.D0*CF*omx**2*nfl - 8*CF*
     &    omx**2*Li2mx + 8*CF*omx**2*Li2omx - 7141.D0/7200.D0*CF*omx**3
     &     - 14863.D0/117600.D0*CF*omx**3*nfl - 2.D0/3.D0*CF*pisq*
     &    omx**(-3) + 2*CF*pisq*omx**(-2) - 7.D0/3.D0*CF*pisq*omx**(-1)
     &     + 5.D0/3.D0*CF*pisq*omx - 2.D0/3.D0*CF*pisq*omx**2 + 8*CF*
     &    ln2*omx**(-2) )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * (  - 20*CF*ln2*omx**(-1) + 56.D
     &    0/3.D0*CF*ln2 - 167.D0/1890.D0*CF*ln2*nfl + 121.D0/15.D0*CF*
     &    ln2*omx + 641.D0/3780.D0*CF*ln2*omx*nfl - 221.D0/15.D0*CF*ln2
     &    *omx**2 - 22.D0/105.D0*CF*ln2*omx**2*nfl + 16.D0/15.D0*CF*ln2
     &    *omx**3 + 4.D0/35.D0*CF*ln2*omx**3*nfl + 24*CF*CA*omx**(-5)*
     &    Li2mx - 14*CF*CA*omx**(-5)*lnx - 24*CF*CA*omx**(-5)*Li2omx + 
     &    10*CF*CA*omx**(-4) - 12*CF*CA*omx**(-4)*lnopxon2 - 96*CF*CA*
     &    omx**(-4)*Li2mx + 56*CF*CA*omx**(-4)*lnx + 96*CF*CA*omx**(-4)
     &    *Li2omx - 41*CF*CA*omx**(-3) + 38*CF*CA*omx**(-3)*lnopxon2 + 
     &    470.D0/3.D0*CF*CA*omx**(-3)*Li2mx - 274.D0/3.D0*CF*CA*
     &    omx**(-3)*lnx - 470.D0/3.D0*CF*CA*omx**(-3)*Li2omx + 395.D0/6.
     &    D0*CF*CA*omx**(-2) - 139.D0/3.D0*CF*CA*omx**(-2)*lnopxon2 - 
     &    400.D0/3.D0*CF*CA*omx**(-2)*Li2mx + 78*CF*CA*omx**(-2)*lnx + 
     &    400.D0/3.D0*CF*CA*omx**(-2)*Li2omx - 619.D0/12.D0*CF*CA*
     &    omx**(-1) + 83.D0/3.D0*CF*CA*omx**(-1)*lnopxon2 + 61*CF*CA*
     &    omx**(-1)*Li2mx )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * (  - 229.D0/6.D0*CF*CA*
     &    omx**(-1)*lnx - 61*CF*CA*omx**(-1)*Li2omx - 1.D0/8.D0*CF*CA*
     &    opx**(-2) - 3.D0/8.D0*CF*CA*opx**(-1) + 9937055.D0/508032.D0*
     &    CF*CA - 19.D0/2.D0*CF*CA*lnopxon2 - 37.D0/3.D0*CF*CA*Li2mx + 
     &    35.D0/3.D0*CF*CA*lnx + 37.D0/3.D0*CF*CA*Li2omx - 36061259.D0/
     &    25401600.D0*CF*CA*omx + 13.D0/6.D0*CF*CA*omx*lnopxon2 - 13.D0/
     &    6.D0*CF*CA*omx*lnx - 2111.D0/352800.D0*CF*CA*omx**2 - 3301.D0/
     &    11760.D0*CF*CA*omx**3 + 2*CF*CA*pisq*omx**(-5) - 8*CF*CA*pisq
     &    *omx**(-4) + 235.D0/18.D0*CF*CA*pisq*omx**(-3) - 100.D0/9.D0*
     &    CF*CA*pisq*omx**(-2) + 61.D0/12.D0*CF*CA*pisq*omx**(-1) - 55.D
     &    0/54.D0*CF*CA*pisq - 1.D0/216.D0*CF*CA*pisq*omx - 24*CF*CA*
     &    ln2*omx**(-4) + 84*CF*CA*ln2*omx**(-3) - 350.D0/3.D0*CF*CA*
     &    ln2*omx**(-2) + 81*CF*CA*ln2*omx**(-1) - 21433.D0/756.D0*CF*
     &    CA*ln2 + 25787.D0/7560.D0*CF*CA*ln2*omx + 8.D0/105.D0*CF*CA*
     &    ln2*omx**2 + 2.D0/7.D0*CF*CA*ln2*omx**3 - 48*CF**2*omx**(-5)*
     &    Li2mx )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * ( 28*CF**2*omx**(-5)*lnx + 48*
     &    CF**2*omx**(-5)*Li2omx - 20*CF**2*omx**(-4) + 24*CF**2*
     &    omx**(-4)*lnopxon2 + 192*CF**2*omx**(-4)*Li2mx - 112*CF**2*
     &    omx**(-4)*lnx - 192*CF**2*omx**(-4)*Li2omx + 82*CF**2*
     &    omx**(-3) - 76*CF**2*omx**(-3)*lnopxon2 - 940.D0/3.D0*CF**2*
     &    omx**(-3)*Li2mx + 548.D0/3.D0*CF**2*omx**(-3)*lnx + 940.D0/3.D
     &    0*CF**2*omx**(-3)*Li2omx - 395.D0/3.D0*CF**2*omx**(-2) + 278.D
     &    0/3.D0*CF**2*omx**(-2)*lnopxon2 + 800.D0/3.D0*CF**2*omx**(-2)
     &    *Li2mx - 156*CF**2*omx**(-2)*lnx - 800.D0/3.D0*CF**2*
     &    omx**(-2)*Li2omx + 619.D0/6.D0*CF**2*omx**(-1) - 166.D0/3.D0*
     &    CF**2*omx**(-1)*lnopxon2 - 122*CF**2*omx**(-1)*Li2mx + 229.D0/
     &    3.D0*CF**2*omx**(-1)*lnx + 122*CF**2*omx**(-1)*Li2omx + 1.D0/
     &    4.D0*CF**2*opx**(-2) + 3.D0/4.D0*CF**2*opx**(-1) - 881821.D0/
     &    21600.D0*CF**2 + 19*CF**2*lnopxon2 + 74.D0/3.D0*CF**2*Li2mx
     &     - 70.D0/3.D0*CF**2*lnx - 74.D0/3.D0*CF**2*Li2omx + 26809.D0/
     &    5400.D0*CF**2*omx )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * (  - 13.D0/3.D0*CF**2*omx*
     &    lnopxon2 + 13.D0/3.D0*CF**2*omx*lnx - 353.D0/900.D0*CF**2*
     &    omx**2 + 239.D0/720.D0*CF**2*omx**3 - 4*CF**2*pisq*omx**(-5)
     &     + 16*CF**2*pisq*omx**(-4) - 235.D0/9.D0*CF**2*pisq*omx**(-3)
     &     + 200.D0/9.D0*CF**2*pisq*omx**(-2) - 61.D0/6.D0*CF**2*pisq*
     &    omx**(-1) + 37.D0/18.D0*CF**2*pisq + 48*CF**2*ln2*omx**(-4)
     &     - 168*CF**2*ln2*omx**(-3) + 700.D0/3.D0*CF**2*ln2*omx**(-2)
     &     - 162*CF**2*ln2*omx**(-1) + 2581.D0/45.D0*CF**2*ln2 - 767.D0/
     &    90.D0*CF**2*ln2*omx + 4.D0/15.D0*CF**2*ln2*omx**2 - 1.D0/3.D0
     &    *CF**2*ln2*omx**3 )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * (  - 8*CF*omx**(-1)*Li2mx + 8*
     &    CF*omx**(-1)*Li2omx - 8*CF - 305.D0/216.D0*CF*nfl + 8*CF*lnx
     &     + 49.D0/18.D0*CF*omx + 107.D0/72.D0*CF*omx*nfl + 8*CF*omx*
     &    lnopxon2 + 8*CF*omx*Li2mx - 16*CF*omx*lnx - 8*CF*omx*Li2omx
     &     - 253.D0/300.D0*CF*omx**2*nfl - 8*CF*omx**2*lnopxon2 + 8*CF*
     &    omx**2*lnx - 2.D0/3.D0*CF*pisq*omx**(-1) + 2.D0/3.D0*CF*pisq*
     &    omx + 8*CF*ln2 + 10.D0/9.D0*CF*ln2*nfl + 16.D0/3.D0*CF*ln2*
     &    omx - 4.D0/3.D0*CF*ln2*omx*nfl - 8*CF*ln2*omx**2 + 4.D0/5.D0*
     &    CF*ln2*omx**2*nfl + 16*CF*CA*omx**(-3)*Li2mx - 8*CF*CA*
     &    omx**(-3)*lnx - 16*CF*CA*omx**(-3)*Li2omx + 8*CF*CA*omx**(-2)
     &     - 8*CF*CA*omx**(-2)*lnopxon2 - 40*CF*CA*omx**(-2)*Li2mx + 24
     &    *CF*CA*omx**(-2)*lnx + 40*CF*CA*omx**(-2)*Li2omx - 16*CF*CA*
     &    omx**(-1) + 16*CF*CA*omx**(-1)*lnopxon2 + 32*CF*CA*omx**(-1)*
     &    Li2mx - 28*CF*CA*omx**(-1)*lnx - 32*CF*CA*omx**(-1)*Li2omx + 
     &    41.D0/216.D0*CF*CA - 12*CF*CA*lnopxon2 - 4*CF*CA*Li2mx + 16*
     &    CF*CA*lnx )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * ( 4*CF*CA*Li2omx + 607.D0/72.D0
     &    *CF*CA*omx + 4*CF*CA*omx*lnopxon2 - 4*CF*CA*omx*Li2mx - 4*CF*
     &    CA*omx*lnx + 4*CF*CA*omx*Li2omx - 847.D0/300.D0*CF*CA*omx**2
     &     + 4.D0/3.D0*CF*CA*pisq*omx**(-3) - 10.D0/3.D0*CF*CA*pisq*
     &    omx**(-2) + 8.D0/3.D0*CF*CA*pisq*omx**(-1) - CF*CA*pisq - 16*
     &    CF*CA*ln2*omx**(-2) + 32*CF*CA*ln2*omx**(-1) - 130.D0/9.D0*CF
     &    *CA*ln2 - 11.D0/3.D0*CF*CA*ln2*omx + 16.D0/5.D0*CF*CA*ln2*
     &    omx**2 - 32*CF**2*omx**(-3)*Li2mx + 16*CF**2*omx**(-3)*lnx + 
     &    32*CF**2*omx**(-3)*Li2omx - 16*CF**2*omx**(-2) + 16*CF**2*
     &    omx**(-2)*lnopxon2 + 80*CF**2*omx**(-2)*Li2mx - 48*CF**2*
     &    omx**(-2)*lnx - 80*CF**2*omx**(-2)*Li2omx + 32*CF**2*
     &    omx**(-1) - 32*CF**2*omx**(-1)*lnopxon2 - 64*CF**2*omx**(-1)*
     &    Li2mx + 56*CF**2*omx**(-1)*lnx + 64*CF**2*omx**(-1)*Li2omx + 
     &    59.D0/9.D0*CF**2 + 24*CF**2*lnopxon2 + 8*CF**2*Li2mx - 32*
     &    CF**2*lnx - 8*CF**2*Li2omx - 118.D0/9.D0*CF**2*omx - 8*CF**2*
     &    omx*lnopxon2 )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * ( 8*CF**2*omx*Li2mx + 8*CF**2*
     &    omx*lnx - 8*CF**2*omx*Li2omx + 11.D0/3.D0*CF**2*omx**2 - 8.D0/
     &    3.D0*CF**2*pisq*omx**(-3) + 20.D0/3.D0*CF**2*pisq*omx**(-2)
     &     - 16.D0/3.D0*CF**2*pisq*omx**(-1) + 10.D0/3.D0*CF**2*pisq - 
     &    2.D0/3.D0*CF**2*pisq*omx + 32*CF**2*ln2*omx**(-2) - 64*CF**2*
     &    ln2*omx**(-1) + 52.D0/3.D0*CF**2*ln2 + 16.D0/3.D0*CF**2*ln2*
     &    omx - 4*CF**2*ln2*omx**2 )
      I2R(qq,rglr) = I2R(qq,rglr) + picube * (  - 4*CF*omx**(-1)*
     &    BCQQBARSx + 4*CF*omx**(-1)*BCQQBARS1 - 4*CF*omx**(-1)*nfl*
     &    BCQQCFTFx + 4*CF*omx**(-1)*nfl*BCQQCFTF1 + 4*CF*BCQQBARSx - 8
     &    *CF*CA*omx**(-1)*BCQQCACFx + 8*CF*CA*omx**(-1)*BCQQCACF1 - 8*
     &    CF**2*omx**(-1)*BCQQCFCFx + 8*CF**2*omx**(-1)*BCQQCFCF1 + 8*
     &    CF**2*BCQQCFCFx )
      I2R(qq,rglr) = I2R(qq,rglr) + lnR * ( 92.D0/9.D0*CF*nfl - 46.D0/9.
     &    D0*CF*omx*nfl - 32.D0/3.D0*CF*ln2*nfl + 16.D0/3.D0*CF*ln2*omx
     &    *nfl - 524.D0/9.D0*CF*CA + 262.D0/9.D0*CF*CA*omx + 16.D0/3.D0
     &    *CF*CA*pisq - 8.D0/3.D0*CF*CA*pisq*omx + 176.D0/3.D0*CF*CA*
     &    ln2 - 88.D0/3.D0*CF*CA*ln2*omx )
      I2R(qq,rglr) = I2R(qq,rglr) + 26.D0/9.D0*CF*nfl - 92.D0/9.D0*CF*
     &    lnpi*nfl + 46.D0/9.D0*CF*lnpi*omx*nfl - 116.D0/9.D0*CF*ln2*
     &    nfl + 46.D0/9.D0*CF*ln2*omx*nfl + 32.D0/3.D0*CF*ln2*lnpi*nfl
     &     - 16.D0/3.D0*CF*ln2*lnpi*omx*nfl + 32.D0/3.D0*CF*ln2**2*nfl
     &     - 16.D0/3.D0*CF*ln2**2*omx*nfl - 26.D0/9.D0*CF*CA + 524.D0/9.
     &    D0*CF*CA*lnpi - 262.D0/9.D0*CF*CA*lnpi*omx - 16.D0/3.D0*CF*CA
     &    *pisq*lnpi + 8.D0/3.D0*CF*CA*pisq*lnpi*omx + 548.D0/9.D0*CF*
     &    CA*ln2 - 262.D0/9.D0*CF*CA*ln2*omx - 176.D0/3.D0*CF*CA*ln2*
     &    lnpi + 88.D0/3.D0*CF*CA*ln2*lnpi*omx - 16.D0/3.D0*CF*CA*ln2*
     &    pisq + 8.D0/3.D0*CF*CA*ln2*pisq*omx - 176.D0/3.D0*CF*CA*
     &    ln2**2 + 88.D0/3.D0*CF*CA*ln2**2*omx
      else
      I2R(qq,rglr)= + R8 * ( 37432579841.D0/35407798272000.D0*CF*nfl - 
     &    44803.D0/38707200.D0*CF*omx - 219702602291.D0/70815596544000.D
     &    0*CF*omx*nfl - 3563.D0/16588800.D0*CF*omx**2 + 768879803.D0/
     &    98354995200.D0*CF*omx**2*nfl + 16515817.D0/406425600.D0*CF*
     &    omx**3 - 4159161163.D0/295064985600.D0*CF*omx**3*nfl - 
     &    3774762301.D0/84304281600.D0*CF*omx**4 + 83202043.D0/
     &    6147187200.D0*CF*omx**4*nfl - 25999035719.D0/1180259942400.D0
     &    *CF*omx**5 - 2062751.D0/409812480.D0*CF*omx**5*nfl + 
     &    12948910797859.D0/1196783581593600.D0*CF*omx**6 + 
     &    1262074977101.D0/199463930265600.D0*CF*omx**7 + 65048978167.D0
     &    /15542643916800.D0*CF*omx**8 + 462519857.D0/155831195520.D0*
     &    CF*omx**9 - 8562227.D0/26127360.D0*CF*omx**10 + 15649297.D0/
     &    52254720.D0*CF*omx**11 - 24821.D0/29937600.D0*CF*ln2*nfl + 1.D
     &    0/1080.D0*CF*ln2*omx + 139241.D0/59875200.D0*CF*ln2*omx*nfl
     &     - 1.D0/1080.D0*CF*ln2*omx**2 - 241.D0/41580.D0*CF*ln2*omx**2
     &    *nfl )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 241.D0/7560.D0*CF*ln2*
     &    omx**3 + 113.D0/10395.D0*CF*ln2*omx**3*nfl + 257.D0/5940.D0*
     &    CF*ln2*omx**4 - 38.D0/3465.D0*CF*ln2*omx**4*nfl + 569.D0/
     &    41580.D0*CF*ln2*omx**5 + 1.D0/231.D0*CF*ln2*omx**5*nfl - 
     &    38243.D0/3243240.D0*CF*ln2*omx**6 - 12713.D0/2162160.D0*CF*
     &    ln2*omx**7 - 67.D0/21060.D0*CF*ln2*omx**8 - 785.D0/432432.D0*
     &    CF*ln2*omx**9 - 949.D0/648.D0*CF*ln2*omx**10 - 1.D0/81.D0*CF*
     &    ln2*omx**11 + 668290504547.D0/70815596544000.D0*CF*CA - 
     &    4432606176083.D0/141631193088000.D0*CF*CA*omx + 12400520203.D0
     &    /357654528000.D0*CF*CA*omx**2 - 42948144601.D0/3934199808000.D
     &    0*CF*CA*omx**3 + 8944028266087.D0/1495979476992000.D0*CF*CA*
     &    omx**4 - 12249080663411.D0/1595711442124800.D0*CF*CA*omx**5
     &     + 9751918291669.D0/35903507447808000.D0*CF*CA*omx**6 + 
     &    18306337797031.D0/143614029791232000.D0*CF*CA*omx**7 + 
     &    2033425369528529.D0/36892848541925376000.D0*CF*CA*omx**8 + 
     &    6399316516853537.D0/332035636877328384000.D0*CF*CA*omx**9 )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 58687.D0/2488320.D0*CF*CA
     &    *omx**10 + 1.D0/777600.D0*CF*CA*pisq - 1.D0/1555200.D0*CF*CA*
     &    pisq*omx - 361901.D0/59875200.D0*CF*CA*ln2 + 516253.D0/
     &    23950080.D0*CF*CA*ln2*omx - 143.D0/5600.D0*CF*CA*ln2*omx**2
     &     + 27949.D0/3326400.D0*CF*CA*ln2*omx**3 - 603503.D0/129729600.
     &    D0*CF*CA*ln2*omx**4 + 14921.D0/2162160.D0*CF*CA*ln2*omx**5 - 
     &    32573.D0/97297200.D0*CF*CA*ln2*omx**6 - 140149.D0/778377600.D0
     &    *CF*CA*ln2*omx**7 - 51577.D0/490089600.D0*CF*CA*ln2*omx**8 - 
     &    218357.D0/3308104800.D0*CF*CA*ln2*omx**9 + 7.D0/27.D0*CF*CA*
     &    ln2*omx**10 - 374402423.D0/24385536000.D0*CF**2 + 121429859.D0
     &    /1354752000.D0*CF**2*omx - 281987213633.D0/1967099904000.D0*
     &    CF**2*omx**2 + 448219881523.D0/5901299712000.D0*CF**2*omx**3
     &     - 1313630254921.D0/106855676928000.D0*CF**2*omx**4 + 
     &    5225830601843.D0/797855721062400.D0*CF**2*omx**5 - 
     &    9751918291669.D0/17951753723904000.D0*CF**2*omx**6 - 
     &    18306337797031.D0/71807014895616000.D0*CF**2*omx**7 )
      I2R(qq,rglr) = I2R(qq,rglr) + R8 * (  - 2033425369528529.D0/
     &    18446424270962688000.D0*CF**2*omx**8 - 6399316516853537.D0/
     &    166017818438664192000.D0*CF**2*omx**9 + 58687.D0/1244160.D0*
     &    CF**2*omx**10 + 1867.D0/226800.D0*CF**2*ln2 - 8927.D0/151200.D
     &    0*CF**2*ln2*omx + 89071.D0/831600.D0*CF**2*ln2*omx**2 - 
     &    104809.D0/1663200.D0*CF**2*ln2*omx**3 + 98249.D0/9266400.D0*
     &    CF**2*ln2*omx**4 - 6263.D0/1081080.D0*CF**2*ln2*omx**5 + 
     &    32573.D0/48648600.D0*CF**2*ln2*omx**6 + 140149.D0/389188800.D0
     &    *CF**2*ln2*omx**7 + 51577.D0/245044800.D0*CF**2*ln2*omx**8 + 
     &    218357.D0/1654052400.D0*CF**2*ln2*omx**9 - 14.D0/27.D0*CF**2*
     &    ln2*omx**10 )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 31946443.D0/3048192000.D0
     &    *CF*nfl + 473.D0/52920.D0*CF*omx + 8250923.D0/304819200.D0*CF
     &    *omx*nfl + 155089.D0/4233600.D0*CF*omx**2 - 7688197.D0/
     &    152409600.D0*CF*omx**2*nfl - 1604599.D0/8467200.D0*CF*omx**3
     &     + 2108497.D0/38102400.D0*CF*omx**3*nfl + 649541.D0/20321280.D
     &    0*CF*omx**4 - 91859.D0/3810240.D0*CF*omx**4*nfl + 5964449.D0/
     &    67737600.D0*CF*omx**5 + 953588653.D0/24588748800.D0*CF*omx**6
     &     + 1271522809.D0/73766246400.D0*CF*omx**7 + 126487248037.D0/
     &    19946393026560.D0*CF*omx**8 + 31162391441.D0/66487976755200.D0
     &    *CF*omx**9 - 16279.D0/9720.D0*CF*omx**10 + 1269757.D0/1451520.
     &    D0*CF*omx**11 + 17.D0/2100.D0*CF*ln2*nfl - 1.D0/126.D0*CF*ln2
     &    *omx - 19.D0/945.D0*CF*ln2*omx*nfl - 11.D0/630.D0*CF*ln2*
     &    omx**2 + 37.D0/945.D0*CF*ln2*omx**2*nfl + 103.D0/630.D0*CF*
     &    ln2*omx**3 - 43.D0/945.D0*CF*ln2*omx**3*nfl - 1.D0/21.D0*CF*
     &    ln2*omx**4 + 4.D0/189.D0*CF*ln2*omx**4*nfl - 23.D0/315.D0*CF*
     &    ln2*omx**5 )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 169.D0/6930.D0*CF*ln2*
     &    omx**6 - 299.D0/41580.D0*CF*ln2*omx**7 - 10.D0/27027.D0*CF*
     &    ln2*omx**8 + 443.D0/180180.D0*CF*ln2*omx**9 - 7.D0/54.D0*CF*
     &    ln2*omx**10 - 4.D0/9.D0*CF*ln2*omx**11 - 56938211.D0/
     &    435456000.D0*CF*CA + 70230283.D0/304819200.D0*CF*CA*omx - 
     &    37263887.D0/304819200.D0*CF*CA*omx**2 + 3918967.D0/152409600.D
     &    0*CF*CA*omx**3 - 952188059.D0/24588748800.D0*CF*CA*omx**4 + 
     &    1128139.D0/447068160.D0*CF*CA*omx**5 + 26487992161.D0/
     &    19946393026560.D0*CF*CA*omx**6 + 105065687.D0/127209139200.D0
     &    *CF*CA*omx**7 + 465984941.D0/779155977600.D0*CF*CA*omx**8 + 
     &    2453900509.D0/5114459750400.D0*CF*CA*omx**9 + 125273.D0/
     &    4354560.D0*CF*CA*omx**10 + 77.D0/900.D0*CF*CA*ln2 - 164.D0/
     &    945.D0*CF*CA*ln2*omx + 197.D0/1890.D0*CF*CA*ln2*omx**2 - 22.D0
     &    /945.D0*CF*CA*ln2*omx**3 + 79.D0/2310.D0*CF*CA*ln2*omx**4 - 1.
     &    D0/252.D0*CF*CA*ln2*omx**5 - 235.D0/108108.D0*CF*CA*ln2*
     &    omx**6 )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 101.D0/77220.D0*CF*CA*ln2
     &    *omx**7 - 457.D0/540540.D0*CF*CA*ln2*omx**8 - 2.D0/3465.D0*CF
     &    *CA*ln2*omx**9 + 4.D0/27.D0*CF*CA*ln2*omx**10 + 235439.D0/
     &    846720.D0*CF**2 - 393241.D0/529200.D0*CF**2*omx + 82471789.D0/
     &    152409600.D0*CF**2*omx**2 - 1302319.D0/15240960.D0*CF**2*
     &    omx**3 + 1303067441.D0/36883123200.D0*CF**2*omx**4 - 1128139.D
     &    0/223534080.D0*CF**2*omx**5 - 26487992161.D0/9973196513280.D0
     &    *CF**2*omx**6 - 105065687.D0/63604569600.D0*CF**2*omx**7 - 
     &    465984941.D0/389577988800.D0*CF**2*omx**8 - 2453900509.D0/
     &    2557229875200.D0*CF**2*omx**9 - 125273.D0/2177280.D0*CF**2*
     &    omx**10 - 11.D0/63.D0*CF**2*ln2 + 176.D0/315.D0*CF**2*ln2*omx
     &     - 439.D0/945.D0*CF**2*ln2*omx**2 + 31.D0/378.D0*CF**2*ln2*
     &    omx**3 - 293.D0/10395.D0*CF**2*ln2*omx**4 + 1.D0/126.D0*CF**2
     &    *ln2*omx**5 + 235.D0/54054.D0*CF**2*ln2*omx**6 + 101.D0/38610.
     &    D0*CF**2*ln2*omx**7 + 457.D0/270270.D0*CF**2*ln2*omx**8 + 4.D0
     &    /3465.D0*CF**2*ln2*omx**9 )
      I2R(qq,rglr) = I2R(qq,rglr) + R6 * (  - 8.D0/27.D0*CF**2*ln2*
     &    omx**10 )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * ( 371677.D0/3175200.D0*CF*nfl
     &     - 211.D0/7200.D0*CF*omx - 2735987.D0/12700800.D0*CF*omx*nfl
     &     - 631.D0/1440.D0*CF*omx**2 + 87469.D0/352800.D0*CF*omx**2*
     &    nfl + 1901.D0/3600.D0*CF*omx**3 - 14863.D0/117600.D0*CF*
     &    omx**3*nfl + 117097.D0/705600.D0*CF*omx**4 - 17903.D0/282240.D
     &    0*CF*omx**5 - 32327.D0/290304.D0*CF*omx**6 - 19541.D0/172800.D
     &    0*CF*omx**7 - 1251729527.D0/12294374400.D0*CF*omx**8 - 
     &    719907829.D0/8196249600.D0*CF*omx**9 - 3983363.D0/1451520.D0*
     &    CF*omx**10 + 518587.D0/725760.D0*CF*omx**11 - 167.D0/1890.D0*
     &    CF*ln2*nfl + 1.D0/15.D0*CF*ln2*omx + 641.D0/3780.D0*CF*ln2*
     &    omx*nfl + 1.D0/3.D0*CF*ln2*omx**2 - 22.D0/105.D0*CF*ln2*
     &    omx**2*nfl - 7.D0/15.D0*CF*ln2*omx**3 + 4.D0/35.D0*CF*ln2*
     &    omx**3*nfl - 8.D0/105.D0*CF*ln2*omx**4 + 2.D0/21.D0*CF*ln2*
     &    omx**5 + 1.D0/9.D0*CF*ln2*omx**6 + 1.D0/10.D0*CF*ln2*omx**7
     &     + 296.D0/3465.D0*CF*ln2*omx**8 + 167.D0/2310.D0*CF*ln2*
     &    omx**9 )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * ( 11.D0/9.D0*CF*ln2*omx**10 - 
     &    8.D0/9.D0*CF*ln2*omx**11 + 21593027.D0/12700800.D0*CF*CA - 
     &    7191319.D0/5080320.D0*CF*CA*omx + 418501.D0/1411200.D0*CF*CA*
     &    omx**2 - 48191.D0/282240.D0*CF*CA*omx**3 + 17161.D0/313600.D0
     &    *CF*CA*omx**4 + 1107689.D0/33868800.D0*CF*CA*omx**5 + 
     &    801904249.D0/36883123200.D0*CF*CA*omx**6 + 227148409.D0/
     &    14753249280.D0*CF*CA*omx**7 + 5389808999.D0/474914119680.D0*
     &    CF*CA*omx**8 + 429027099437.D0/49865982566400.D0*CF*CA*omx**9
     &     + 6643.D0/27648.D0*CF*CA*omx**10 + 1.D0/108.D0*CF*CA*pisq - 
     &    1.D0/216.D0*CF*CA*pisq*omx - 3761.D0/3780.D0*CF*CA*ln2 + 1663.
     &    D0/1512.D0*CF*CA*ln2*omx - 23.D0/70.D0*CF*CA*ln2*omx**2 + 13.D
     &    0/84.D0*CF*CA*ln2*omx**3 - 73.D0/1260.D0*CF*CA*ln2*omx**4 - 
     &    13.D0/420.D0*CF*CA*ln2*omx**5 - 389.D0/20790.D0*CF*CA*ln2*
     &    omx**6 - 205.D0/16632.D0*CF*CA*ln2*omx**7 - 89.D0/10296.D0*CF
     &    *CA*ln2*omx**8 - 859.D0/135135.D0*CF*CA*ln2*omx**9 - 110279.D0
     &    /21600.D0*CF**2 )
      I2R(qq,rglr) = I2R(qq,rglr) + R4 * ( 17843.D0/3600.D0*CF**2*omx
     &     - 703697.D0/705600.D0*CF**2*omx**2 + 15811.D0/141120.D0*
     &    CF**2*omx**3 - 17161.D0/156800.D0*CF**2*omx**4 - 1107689.D0/
     &    16934400.D0*CF**2*omx**5 - 801904249.D0/18441561600.D0*CF**2*
     &    omx**6 - 227148409.D0/7376624640.D0*CF**2*omx**7 - 5389808999.
     &    D0/237457059840.D0*CF**2*omx**8 - 429027099437.D0/
     &    24932991283200.D0*CF**2*omx**9 - 6643.D0/13824.D0*CF**2*
     &    omx**10 + 119.D0/45.D0*CF**2*ln2 - 39.D0/10.D0*CF**2*ln2*omx
     &     + 113.D0/105.D0*CF**2*ln2*omx**2 - 1.D0/14.D0*CF**2*ln2*
     &    omx**3 + 73.D0/630.D0*CF**2*ln2*omx**4 + 13.D0/210.D0*CF**2*
     &    ln2*omx**5 + 389.D0/10395.D0*CF**2*ln2*omx**6 + 205.D0/8316.D0
     &    *CF**2*ln2*omx**7 + 89.D0/5148.D0*CF**2*ln2*omx**8 + 1718.D0/
     &    135135.D0*CF**2*ln2*omx**9 )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * (  - 305.D0/216.D0*CF*nfl - 23.
     &    D0/18.D0*CF*omx + 107.D0/72.D0*CF*omx*nfl + 23.D0/9.D0*CF*
     &    omx**2 - 253.D0/300.D0*CF*omx**2*nfl - 11.D0/6.D0*CF*omx**3
     &     - 2059.D0/1800.D0*CF*omx**4 - 541.D0/720.D0*CF*omx**5 - 
     &    30391.D0/58800.D0*CF*omx**6 - 373.D0/1008.D0*CF*omx**7 - 
     &    1393619.D0/5080320.D0*CF*omx**8 - 169031.D0/806400.D0*CF*
     &    omx**9 + 11867.D0/181440.D0*CF*omx**10 - 511.D0/576.D0*CF*
     &    omx**11 + 10.D0/9.D0*CF*ln2*nfl + 4.D0/3.D0*CF*ln2*omx - 4.D0/
     &    3.D0*CF*ln2*omx*nfl - 8.D0/3.D0*CF*ln2*omx**2 + 4.D0/5.D0*CF*
     &    ln2*omx**2*nfl + 2*CF*ln2*omx**3 + 16.D0/15.D0*CF*ln2*omx**4
     &     + 2.D0/3.D0*CF*ln2*omx**5 + 16.D0/35.D0*CF*ln2*omx**6 + 1.D0/
     &    3.D0*CF*ln2*omx**7 + 16.D0/63.D0*CF*ln2*omx**8 + 1.D0/5.D0*CF
     &    *ln2*omx**9 + 8.D0/9.D0*CF*ln2*omx**10 - 1135.D0/216.D0*CF*CA
     &     + 15.D0/8.D0*CF*CA*omx - 1459.D0/1800.D0*CF*CA*omx**2 + 13.D0
     &    /18.D0*CF*CA*omx**3 + 21211.D0/58800.D0*CF*CA*omx**4 + 9281.D0
     &    /44100.D0*CF*CA*omx**5 )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * ( 3432551.D0/25401600.D0*CF*CA
     &    *omx**6 + 1178879.D0/12700800.D0*CF*CA*omx**7 + 41241581.D0/
     &    614718720.D0*CF*CA*omx**8 + 44329601.D0/878169600.D0*CF*CA*
     &    omx**9 + 1162447.D0/1451520.D0*CF*CA*omx**10 - 2.D0/3.D0*CF*
     &    CA*pisq + 1.D0/3.D0*CF*CA*pisq*omx + 26.D0/9.D0*CF*CA*ln2 - 
     &    CF*CA*ln2*omx + 16.D0/15.D0*CF*CA*ln2*omx**2 - 2.D0/3.D0*CF*
     &    CA*ln2*omx**3 - 11.D0/35.D0*CF*CA*ln2*omx**4 - 19.D0/105.D0*
     &    CF*CA*ln2*omx**5 - 37.D0/315.D0*CF*CA*ln2*omx**6 - 26.D0/315.D
     &    0*CF*CA*ln2*omx**7 - 85.D0/1386.D0*CF*CA*ln2*omx**8 - 47.D0/
     &    990.D0*CF*CA*ln2*omx**9 - 4.D0/9.D0*CF*CA*ln2*omx**10 + 157.D0
     &    /9.D0*CF**2 - 323.D0/900.D0*CF**2*omx**2 - 13.D0/9.D0*CF**2*
     &    omx**3 - 21211.D0/29400.D0*CF**2*omx**4 - 9281.D0/22050.D0*
     &    CF**2*omx**5 - 3432551.D0/12700800.D0*CF**2*omx**6 - 1178879.D
     &    0/6350400.D0*CF**2*omx**7 - 41241581.D0/307359360.D0*CF**2*
     &    omx**8 - 44329601.D0/439084800.D0*CF**2*omx**9 - 1162447.D0/
     &    725760.D0*CF**2*omx**10 )
      I2R(qq,rglr) = I2R(qq,rglr) + R2 * ( 8.D0/3.D0*CF**2*pisq - 4.D0/
     &    3.D0*CF**2*pisq*omx - 52.D0/3.D0*CF**2*ln2 + 4.D0/15.D0*CF**2
     &    *ln2*omx**2 + 4.D0/3.D0*CF**2*ln2*omx**3 + 22.D0/35.D0*CF**2*
     &    ln2*omx**4 + 38.D0/105.D0*CF**2*ln2*omx**5 + 74.D0/315.D0*
     &    CF**2*ln2*omx**6 + 52.D0/315.D0*CF**2*ln2*omx**7 + 85.D0/693.D
     &    0*CF**2*ln2*omx**8 + 47.D0/495.D0*CF**2*ln2*omx**9 + 8.D0/9.D0
     &    *CF**2*ln2*omx**10 )
      I2R(qq,rglr) = I2R(qq,rglr) + picube * (  - 4*CF*omx**(-1)*
     &    BCQQBARSx + 4*CF*omx**(-1)*BCQQBARS1 - 4*CF*omx**(-1)*nfl*
     &    BCQQCFTFx + 4*CF*omx**(-1)*nfl*BCQQCFTF1 + 4*CF*BCQQBARSx - 8
     &    *CF*CA*omx**(-1)*BCQQCACFx + 8*CF*CA*omx**(-1)*BCQQCACF1 - 8*
     &    CF**2*omx**(-1)*BCQQCFCFx + 8*CF**2*omx**(-1)*BCQQCFCF1 + 8*
     &    CF**2*BCQQCFCFx )
      I2R(qq,rglr) = I2R(qq,rglr) + lnR * ( 92.D0/9.D0*CF*nfl - 46.D0/9.
     &    D0*CF*omx*nfl - 32.D0/3.D0*CF*ln2*nfl + 16.D0/3.D0*CF*ln2*omx
     &    *nfl - 524.D0/9.D0*CF*CA + 262.D0/9.D0*CF*CA*omx + 16.D0/3.D0
     &    *CF*CA*pisq - 8.D0/3.D0*CF*CA*pisq*omx + 176.D0/3.D0*CF*CA*
     &    ln2 - 88.D0/3.D0*CF*CA*ln2*omx )
      I2R(qq,rglr) = I2R(qq,rglr) + 26.D0/9.D0*CF*nfl - 92.D0/9.D0*CF*
     &    lnpi*nfl + 46.D0/9.D0*CF*lnpi*omx*nfl - 116.D0/9.D0*CF*ln2*
     &    nfl + 46.D0/9.D0*CF*ln2*omx*nfl + 32.D0/3.D0*CF*ln2*lnpi*nfl
     &     - 16.D0/3.D0*CF*ln2*lnpi*omx*nfl + 32.D0/3.D0*CF*ln2**2*nfl
     &     - 16.D0/3.D0*CF*ln2**2*omx*nfl - 26.D0/9.D0*CF*CA + 524.D0/9.
     &    D0*CF*CA*lnpi - 262.D0/9.D0*CF*CA*lnpi*omx - 16.D0/3.D0*CF*CA
     &    *pisq*lnpi + 8.D0/3.D0*CF*CA*pisq*lnpi*omx + 548.D0/9.D0*CF*
     &    CA*ln2 - 262.D0/9.D0*CF*CA*ln2*omx - 176.D0/3.D0*CF*CA*ln2*
     &    lnpi + 88.D0/3.D0*CF*CA*ln2*lnpi*omx - 16.D0/3.D0*CF*CA*ln2*
     &    pisq + 8.D0/3.D0*CF*CA*ln2*pisq*omx - 176.D0/3.D0*CF*CA*
     &    ln2**2 + 88.D0/3.D0*CF*CA*ln2**2*omx
      endif
