
c     ALGORITHM 490 COLLECTED ALGORITHMS FROM ACM.
c     ALGORITHM APPEARED IN COMM. ACM, VOL. 18, NO. 4N,
c     P. 200.
      DOUBLE PRECISION FUNCTION DILOG(X)
c REAL PART OF THE DILOGARITHM FUNCTION FOR A REAL
c ARGUMENT. REF. NO. 1=L. LEWIN, *DILOGARITHMS +
c ASSOCIATED FUNCTIONS*
c                      (MAC-DONALD, LONDON, 1958).
c NUMERICAL CONSTANTS USED ARE C(N)=(N(N+1)(N+2))**2
c FOR N=1 TO 30, (PI**2)/3=3.289868...,
c (PI**2)/6=1.644394..., AND ZERO OF DILOG ON THE
c POSITIVE REAL AXIS, X0=12.59517...
      DOUBLE PRECISION A, B, BY, C, C1, C2, C3, C4,
     & DX, DY, TEST, W, X, X0, Y, Z
      DIMENSION C(30)
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7),
     & C(8), C(9), C(10), C(11), C(12), C(13),
     & C(14), C(15), C(16), C(17), C(18), C(19),
     & C(20), C(21), C(22), C(23), C(24), C(25),
     & C(26), C(27), C(28), C(29), C(30)
     & /36.D0,576.D0,36.D2,144.D2,441.D2,112896.D0,
     & 254016.D0,5184.D2,9801.D2,17424.D2,2944656.D0,
     & 4769856.D0,74529.D2,112896.D2,166464.D2,
     & 23970816.D0,33802596.D0,467856.D2,636804.D2,
     & 853776.D2,112911876.D0,147476736.D0,19044.D4,
     & 24336.D4,3080025.D2,386358336.D0,480661776.D0,
     & 5934096.D2,7273809.D2,8856576.D2/
      IF (X.GT.12.6D0) GO TO 10
      IF (X.GE.12.59D0) GO TO 100
      IF (X.GE.2.D0) GO TO 10
      IF (X.GT.1.D0) GO TO 20
      IF (X.EQ.1.D0) GO TO 30
      IF (X.GT..5D0) GO TO 40
      IF (X.GT.1.D-2) GO TO 50
      IF (X.LT.-1.D0) GO TO 60
      IF (X.LT.-1.D-2) GO TO 70
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(1).
      DILOG = X*(1.D0+X*(.25D0+X*(1.D0/9.D0+X*
     & (625.D-4+X*(4.D-2+X*(1.D0/36.D0+X*(1.D0/
     & 49.D0+X/64.D0)))))))
      RETURN
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(6),
c AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   10 Y = 1.D0/X
      BY = -1.D0 - Y*(4.D0+Y)
      DILOG = 3.28986813369645287D0 -
     & .5D0*DLOG(X)**2 + (Y*(4.D0+5.75D0*Y)+3.D0*
     & (1.D0+Y)*(1.D0-Y)*DLOG(1.D0-Y))/BY
      IF (DILOG+4.D0*Y.EQ.DILOG) RETURN
      GO TO 80
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7) WITH
c X=1/X + EQ(6), AND DESCRIPTION OF THIS ALGORITHM,
c EQ(4).
   20 Y = 1.D0 - 1.D0/X
      DX = DLOG(X)
      BY = 1.D0 + Y*(4.D0+Y)
      DILOG = 1.64493406684822643D0 +
     & DX*(.5D0*DX-DLOG(X-1.D0)) +
     & (Y*(4.D0+5.75D0*Y)-3.D0*(1.D0+Y)*DX/X)/BY
      GO TO 80
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(2).
   30 DILOG = 1.64493406684822643D0
      RETURN
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(7),
c AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   40 Y = 1.D0 - X
      DX = DLOG(X)
      BY = -1.D0 - Y*(4.D0+Y)
      DILOG = 1.64493406684822643D0 - DX*DLOG(Y) +
     & (Y*(4.D0+5.75D0*Y)+3.D0*(1.D0+Y)*DX*X)/BY
      GO TO 80
c DILOG COMPUTED FROM DESCRIPTION OF THIS ALGORITHM,
c EQ(4)
   50 Y = X
      BY = 1.D0 + Y*(4.D0+Y)
      DILOG = (Y*(4.D0+5.75D0*Y)+3.D0*(1.D0+Y)*
     & (1.D0-Y)*DLOG(1.D0-Y))/BY
      GO TO 80
c DILOG COMPUTED FROM REF. NO. 1, P.245, EQ(12) WITH
c X=-X, AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   60 Y = 1.D0/(1.D0-X)
      DX = DLOG(-X)
      DY = DLOG(Y)
      BY = 1.D0 + Y*(4.D0+Y)
      DILOG = -1.64493406684822643D0 +
     & .5D0*DY*(DY+2.D0*DX) + (Y*(4.D0+5.75D0*Y)
     & +3.D0*(1.D0+Y)*(1.D0-Y)*(DX+DY))/BY
      IF (DILOG+4.D0*Y.EQ.DILOG) RETURN
      GO TO 80
c DILOG COMPUTED FROM REF. NO. 1, P.244, EQ(8),
c AND DESCRIPTION OF THIS ALGORITHM, EQ(4).
   70 Y = X/(X-1.D0)
      DX = DLOG(1.D0-X)
      BY = -1.D0 - Y*(4.D0+Y)
      DILOG = (Y*(4.D0+5.75D0*Y)-3.D0*(1.D0+Y)*
     & (1.D0-Y)*DX)/BY - .5D0*DX*DX
   80 B = 4.D0*Y*Y/BY
      DO 90 N=1,30
        B = B*Y
        A = B/C(N)
        TEST = DILOG
        DILOG = DILOG + A
        IF (DILOG.EQ.TEST) RETURN
   90 CONTINUE
      RETURN
c DILOG COMPUTED FROM TAYLOR SERIES ABOUT ZERO OF
c DILOG, X0.
  100 X0 = 12.5951703698450184D0
      Y = X/X0 - 1.D0
      Z = 1.D0/11.5951703698450184D0
      W = Y*Z
      C1 = (3.D0*X0-2.D0)/6.D0
      C2 = ((11.D0*X0-15.D0)*X0+6.D0)/24.D0
      C3 = (((5.D1*X0-104.D0)*X0+84.D0)*X0-24.D0)/
     & 12.D1
      C4 = ((((274.D0*X0-77.D1)*X0+94.D1)*X0-54.D1)*
     & X0+12.D1)/72.D1
      DILOG = Y*(1.D0-Y*(.5D0-Y*(1.D0/3.D0-Y*
     & (.25D0-Y*(.2D0-Y/6.D0)))))*DLOG(Z) -
     & W*X0*Y*(.5D0-W*(C1-W*(C2-W*(C3-W*C4))))
      RETURN
      END
