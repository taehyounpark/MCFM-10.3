!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c###[ ffxdb0:
	subroutine ffxdb0(cdb0,cdb0p,xp,xma,xmb,ier)
c**#[*comment:***********************************************************
c									*
c	Calculates the the derivative of the two-point function with	*
c	respect to p2 and the same times p2 (one is always well-defined)*
c									*
c	Input:	xp	(real)	  k2, in B&D metric			*
c		xma	(real)	  mass2					*
c		xmb	(real)	  mass2					*
c									*
c	Output:	cdb0	(complex) dB0/dxp				*
c		cdb0p	(complex) xp*dB0/dxp				*
c		ier	(integer) # of digits lost, if >=100: error	*
c									*
c	Calls:	ffxdba							*
c									*
c**#]*comment:***********************************************************
c  #[ declarations:
	implicit none

c	arguments

	integer ier
	DOUBLE COMPLEX cdb0,cdb0p
	DOUBLE PRECISION xp,xma,xmb

c	local variables

	integer ier0
	DOUBLE PRECISION dmamb,dmap,dmbp

c	common blocks

	include 'ff.h'

c  #] declarations:
c  #[ check input:
	if ( lwrite ) then
	    print *,'ffxdb0: input:'
	    print *,'xma,xmb,xp,ier = ',xma,xmb,xp,ier
	endif
	if ( ltest ) then
	    if ( xma  <  0 .or. xmb  <  0 ) then
		print *,'ffxdb0: error: xma,b < 0: ',xma,xmb
		stop
	    endif
	endif
c  #] check input:
c  #[ get differences:
	ier0 = 0
	dmamb = xma - xmb
	dmap = xma - xp
	dmbp = xmb - xp
	if ( lwarn ) then
	    if ( abs(dmamb)  <  xloss*abs(xma) .and. xma  /=  xmb )
     &		call ffwarn(97,ier0,dmamb,xma)
	    if ( abs(dmap)  <  xloss*abs(xp) .and. xp  /=  xma )
     &		call ffwarn(98,ier0,dmap,xp)
	    if ( abs(dmbp)  <  xloss*abs(xp) .and. xp  /=  xmb )
     &		call ffwarn(99,ier0,dmbp,xp)
	endif
c  #] get differences:
c  #[ calculations:
	call ffxdbp(cdb0,cdb0p,xp,xma,xmb,dmap,dmbp,dmamb,ier)
	if ( lwrite ) print *,'B0'' = ',cdb0,cdb0p,ier
c  #] calculations:
c###] ffxdb0:
	end
c###[ ffxdbp:
	subroutine ffxdbp(cdb0,cdb0p,xp,xma,xmb,dmap,dmbp,dmamb,ier)
c**#[*comment:***********************************************************
c									*
c	calculates the derivatives of the two-point function		*
c	Veltman) for all possible cases: masses equal, unequal,		*
c	equal to zero.							*
c									*
c	Input:	xp	(real) p.p, in B&D metric			*
c		xma	(real) mass2,					*
c		xmb	(real) mass2,					*
c		dm[ab]p	(real) xm[ab] - xp				*
c		dmamb	(real) xma - xmb				*
c									*
c	Output:	cdb0	(complex) B0' = dB0/dxp				*
c		cdb0p	(complex) xp*dB0/dxp				*
c		ier	(integer) 0=ok,>0=numerical problems,>100=error	*
c									*
c	Calls:	ffxdbp.							*
c									*
c**#]*comment:***********************************************************
c  #[ declarations:
	implicit none

c	arguments

	integer ier
	DOUBLE COMPLEX cdb0,cdb0p
	DOUBLE PRECISION xp,xma,xmb,dmap,dmbp,dmamb

c	local variables

	integer i,initeq,jsign,initir
	DOUBLE PRECISION ax,ffbnd,
     &		xprceq,bdeq01,bdeq05,bdeq11,bdeq17,bdeq25,
     &		xprcn3,bdn301,bdn305,bdn310,bdn315,
     &		xprcn5,bdn501,bdn505,bdn510,bdn515,
     &		xprec0,bdn001,bdn005,bdn010,bdn015,bdn020
	DOUBLE PRECISION xcheck,xm,dmp,xm1,xm2,dm1m2,dm1p,
     &		dm2p,s,s1,s1a,s1b,s1p,s2,s2a,s2b,s2p,x,y,som,
     &		xlam,slam,xlogmm,alpha,alph1,xnoe,xpneq(30),
     &		xx,dfflo1,dfflo3,d1,d2,diff,h,a,b,c,d,beta,
     &		betm2n,xmax,s1c,s1d,s1e,s1f,s3
	DOUBLE COMPLEX cc,zxfflg
	save initeq,xpneq,initir,
     &		xprceq,bdeq01,bdeq05,bdeq11,bdeq17,bdeq25,
     &		xprcn3,bdn301,bdn305,bdn310,bdn315,
     &		xprcn5,bdn501,bdn505,bdn510,bdn515,
     &		xprec0,bdn001,bdn005,bdn010,bdn015,bdn020

c	common blocks

	include 'ff.h'
	DOUBLE PRECISION delta
	common /ffcut/ delta

c	data

	data xprceq /-1./
	data xprec0 /-1./
	data xprcn3 /-1./
	data xprcn5 /-1./
	data initeq /0/

c  #] declarations:
!$omp threadprivate(initeq,xpneq,initir,
!$omp&		xprceq,bdeq01,bdeq05,bdeq11,bdeq17,bdeq25,
!$omp&		xprcn3,bdn301,bdn305,bdn310,bdn315,
!$omp&		xprcn5,bdn501,bdn505,bdn510,bdn515,
!$omp&		xprec0,bdn001,bdn005,bdn010,bdn015,bdn020,/ffcut/)

c  #[ check input:
	if (ltest) then
	    xcheck = xma - xmb - dmamb
	    if ( abs(xcheck)  >  precx*max(abs(xma),abs(xmb),abs(
     &			dmamb))/xloss ) then
		print *,'ffxdbp: input not OK, dmamb <> xma-xmb',xcheck
	    endif
	    xcheck = -xp + xma - dmap
	    if ( abs(xcheck)  >  precx*max(abs(xp),abs(xma),abs(
     &			dmap))/xloss ) then
		print *,'ffxdbp: input not OK, dmap <> xma - xp',xcheck
	    endif
	    xcheck = -xp + xmb - dmbp
	    if ( abs(xcheck)  >  precx*max(abs(xp),abs(xmb),abs(
     &			dmbp))/xloss ) then
		print *,'ffxdbp: input not OK, dmbp <> xmb - xp',xcheck
	    endif
	endif
c  #] check input:
c  #[ which case:

c	sort according to the type of masscombination encountered:
c	100: both masses zero, 200: one equal to zero, 300: both equal
c	400: rest.

	if ( xma  ==  0 ) then
		if ( xmb  ==  0 ) then
			goto 100
		endif
		xm = xmb
		dmp = dmbp
		goto 200
	endif
	if ( xmb  ==  0 ) then
		xm = xma
		dmp = dmap
		goto 200
	elseif ( dmamb  ==  0 ) then
		xm = xma
		dmp = dmap
		goto 300
	elseif ( xma  >  xmb ) then
		xm2 = xma
		xm1 = xmb
		dm1m2 = -dmamb
		dm1p = dmbp
		dm2p = dmap
	else
		xm1 = xma
		xm2 = xmb
		dm1m2 = dmamb
		dm1p = dmap
		dm2p = dmbp
	endif
	goto 400
c  #] which case:
c  #[ both masses equal to zero:
  100	continue
	if ( xp /= 0 ) cdb0 = -1/xp
	cdb0p = -1
	return
c  #] both masses equal to zero:
c  #[ one mass equal to zero:
  200	continue

c	special case xp = 0

	if ( xp  ==  0 ) then
	    cdb0p = 0
	    cdb0 = 1/(2*xm)
	    goto 990

c	special case xp = xm

	elseif ( dmp == 0 ) then
	    if ( lsmug ) then
		if ( DBLE(cmipj(1,3)) < DBLE(cmipj(2,3)) ) then
		    cdb0p = -1 - log(cmipj(1,3)*DBLE(1/xm))
		else
		    cdb0p = -1 - log(cmipj(2,3)*DBLE(1/xm))
		endif
	    else
		if ( initir == 0 ) then
		    initir = 1
		    print *,'ffxdb0: IR divergent B0'', using cutoff ',
     &		    	delta
		endif
		if ( delta == 0 ) then
		    call fferr(74,ier)
		    cdb0p = 0
		else
		    cdb0p = -1 + log(xm/delta)/2
		endif
	    endif
	    cdb0 = cdb0p*(1/DBLE(xp))
	    goto 990
	endif

c	Normal case:

	x = xp/xm
	ax = abs(x)
	if ( ax  <  xloss ) then
c 	#[ Taylor expansion:
	    if ( xprec0  /=  precx ) then
		xprec0 = precx
		bdn001 = ffbnd(2,1,xninv)
		bdn005 = ffbnd(2,5,xninv)
		bdn010 = ffbnd(2,10,xninv)
		bdn015 = ffbnd(2,15,xninv)
		bdn020 = ffbnd(2,20,xninv)
	    endif
	    if ( lwarn .and. ax  >  bdn020 ) then
	    	call ffwarn(15,ier,precx,xninv(21)*ax**20)
	    endif
	    if ( ax  >  bdn015 ) then
		som = x*(xninv(17) + x*(xninv(18) + x*(xninv(19) +
     &		      x*(xninv(20) + x*(xninv(21) )))))
	    else
		som = 0
	    endif
	    if ( ax  >  bdn010 ) then
		som = x*(xninv(12) + x*(xninv(13) + x*(xninv(14) +
     &		      x*(xninv(15) + x*(xninv(16) + som )))))
	    endif
	    if ( ax  >  bdn005 ) then
		som = x*(xninv(7) + x*(xninv(8) + x*(xninv(9) +
     &		      x*(xninv(10) + x*(xninv(11) + som )))))
	    endif
	    if ( ax  >  bdn001 ) then
		som = x*(xninv(3) + x*(xninv(4) + x*(xninv(5) +
     &		      x*(xninv(6) + som ))))
	    endif
	    cdb0p = x*(xninv(2) + som)
	    if ( lwrite ) then
		print *,'cdb0p = ',cdb0p
		print *,'verg   ',-1 - xm/xp*dfflo1(x,ier),1
	    endif
c 	#] Taylor expansion:
	else
c 	#[ short formula:
	    s = log(abs(dmp/xm))
	    cdb0p = -(1 + s*xm/xp)
	    if ( xp > xm ) cdb0p = cdb0p+DCMPLX(DBLE(0),DBLE(xm/xp*pi))
c 	#] short formula:
	endif
	cdb0 = cdb0p*(1/DBLE(xp))
	goto 990
c  #] one mass equal to zero:
c  #[ both masses equal:
  300	continue

c	Both masses are equal.	Not only this speeds up things, some
c	cancellations have to be avoided as well.

c	first a special case

	if ( abs(xp)  <  8*xloss*xm ) then
c -#[	    taylor expansion:

c	    a Taylor expansion seems appropriate as the result will go
c	    as k^2 but seems to go as 1/k !!

c--#[	    data and bounds:
	    if ( initeq  ==  0 ) then
		initeq = 1
		xpneq(1) = x1/6
		do 1 i=2,30
		    xpneq(i) = - xpneq(i-1)*DBLE(i)/DBLE(2*(2*i+1))
    1		continue
	    endif
	    if (xprceq  /=  precx ) then

c		calculate the boundaries for the number of terms to be
c		included in the taylorexpansion

		xprceq = precx
		bdeq01 = ffbnd(1,1,xpneq)
		bdeq05 = ffbnd(1,5,xpneq)
		bdeq11 = ffbnd(1,11,xpneq)
		bdeq17 = ffbnd(1,17,xpneq)
		bdeq25 = ffbnd(1,25,xpneq)
	    endif
c--#]	    data and bounds:
	    x = -xp/xm
	    ax = abs(x)
	    if ( lwarn .and. ax  >  bdeq25 ) then
		call ffwarn(15,ier,precx,abs(xpneq(25))*ax**25)
	    endif
	    if ( ax  >  bdeq17 ) then
		som = x*(xpneq(18) + x*(xpneq(19) + x*(xpneq(20) +
     &		x*(xpneq(21) + x*(xpneq(22) + x*(xpneq(23) +
     &		x*(xpneq(24) + x*(xpneq(25) ))))))))
	    else
		som = 0
	    endif
	    if ( ax  >  bdeq11 ) then
		som = x*(xpneq(12) + x*(xpneq(13) + x*(xpneq(14) +
     &		x*(xpneq(15) + x*(xpneq(16) + x*(xpneq(17) + som ))))
     &		))
	    endif
	    if ( ax  >  bdeq05 ) then
		som = x*(xpneq(6) + x*(xpneq(7) + x*(xpneq(8) + x*(
     &		xpneq(9) + x*(xpneq(10) + x*(xpneq(11) + som ))))))
	    endif
	    if ( ax  >  bdeq01 ) then
		som = x*(xpneq(2) + x*(xpneq(3) + x*(xpneq(4) + x*(
     &		xpneq(5) + som ))))
	    endif
	    cdb0p = -x*(xpneq(1)+som)
	    if (lwrite) then
		print *,'ffxdbp: m1 = m2, Taylor expansion in ',x
		print *,'cdb0p = ',cdb0p
	    endif
	    if ( xp /= 0 ) then
		cdb0 = cdb0p*(1/DBLE(xp))
	    else
		cdb0 = xpneq(1)/xm
	    endif
	    goto 990
c -#]	    taylor expansion:
	endif
c -#[	normal case:

c	normal case

	call ffxlmb(xlam,-xp,-xm,-xm,dmp,dmp,x0,ier)
	if ( xlam  ==  0 ) then
	    call fferr(86,ier)
	    return
	elseif ( xlam  >  0 ) then
c	    cases 1,2 and 4
	    slam = sqrt(xlam)
	    s2a = dmp + xm
	    s2 = s2a + slam
	    if ( abs(s2)  >  xloss*slam ) then
c		looks fine
		jsign = 1
	    else
		s2 = s2a - slam
		jsign = -1
	    endif
	    ax = abs(s2/(2*xm))
	    if ( ax  <  xalogm ) then
		if ( lwarn ) call ffwarn(16,ier,ax,xalogm)
		s = 0
	    elseif( ax-1  <  .1 .and. s2  >  0 ) then
c		In this case a quicker and more accurate way is to
c		calculate log(1-x).
		s2 = (xp - slam)
c		the following line is superfluous.
		if ( lwarn .and. abs(s2)  <  xloss*slam )
     &			call ffwarn(17,ier,s2,slam)
		s = 2*xm/slam*dfflo1(s2/(2*xm),ier)
	    else
c		finally the normal case
		s = 2*xm/slam*log(ax)
		if ( jsign  ==  -1 ) s = -s
	    endif
	    if ( xp  >  2*xm ) then
c		in this case ( xlam>0, so xp>(2*m)^2) ) there also
c		is an imaginary part
		y = pi*2*xm/slam
	    else
		y = 0
	    endif
	else
c	    the root is complex (k^2 between 0 and (2*m1)^2)
	    slam = sqrt(-xlam)
	    s = 4*xm/slam*atan2(xp,slam)
	    y = 0
	endif
	if (lwrite) print *,'s =   ',s
	xx = s - 1
	if ( lwarn .and. abs(xx) < xloss ) call ffwarn(18,ier,xx,x1)
	cdb0p = DCMPLX(DBLE(xx),DBLE(y))
	cdb0 = cdb0p*(1/DBLE(xp))
	goto 990
c -#]	normal case:

c  #] both masses equal:
c  #[ unequal nonzero masses:
c -#[	get log(xm2/xm1):
  400	continue
	x = xm2/xm1
	if ( 1  <  xalogm*x ) then
	    call fferr(8,ier)
	    xlogmm = 0
	elseif ( abs(x-1)  <  xloss ) then
	    xlogmm = dfflo1(dm1m2/xm1,ier)
	else
	    xlogmm = log(x)
	endif
c -#]	get log(xm2/xm1):
c -#[	xp = 0:

c	first a special case

	if ( xp  ==  0 ) then

c	    repaired 19-nov-1993, see b2.frm

	    s1 = xm1*xm2*xlogmm/dm1m2**3
	    s2 = (xm1+xm2)/(2*dm1m2**2)
	    s = s1 + s2
	    if ( abs(s)  <  xloss**2*s2 ) then

c		second try

		h = dfflo3(dm1m2/xm1,ier)
		s1 = -xm1*h/dm1m2**2
		s2 = 1/(2*xm1)
		s3 = xm1**2*h/dm1m2**3
		s = s1 + s2 + s3
		if ( abs(s)  <  xloss*max(abs(s2),abs(s3)) ) then
		    call ffwarn(228,ier,s,s2)
		endif
	    endif
	    cdb0 = s
	    cdb0p = 0
	    goto 990
	endif
c -#]	xp = 0:
c -#[	normal case:

c	proceeding with the normal case

	call ffxlmb(xlam,-xp,-xm2,-xm1,dm2p,dm1p,dm1m2,ier)
	diff = xlam + xp*(dm2p+xm1)
	if ( lwrite ) print *,'diff = ',diff
	if ( abs(diff)  <  xloss*xlam ) then
	    h = dm1m2**2 - xp*(xm1+xm2)
	    if ( lwrite ) print *,'diff+= ',h
	    if ( abs(h)  <  xloss*dm1m2**2 ) then
		if ( dm1m2**2  <  abs(xlam) ) diff = h
		if ( lwarn ) then
		    call ffwarn(221,ier,diff,min(dm1m2**2,abs(xlam)))
		endif
	    endif
	endif
	if ( xlam  ==  0 ) then
	    call fferr(86,ier)
	    return
	elseif ( xlam  >  0 ) then
c	    cases k^2 < -(m2+m1)^2 or k^2 > -(m2-m1)^2:
c--#[	    first try:
c	    first try the normal way
	    slam = sqrt(xlam)
	    s2a = dm2p + xm1
	    s2 = s2a + slam
	    if ( abs(s2)  >  xloss*slam ) then
c		looks fine
		jsign = 1
	    else
		s2 = s2a - slam
		jsign = -1
	    endif
	    s2 = s2**2/(4*xm1*xm2)
	    if ( abs(s2)  <  xalogm ) then
		call fferr(9,ier)
		s2 = 0
	    elseif ( abs(s2-1)  <  xloss ) then
		if ( jsign == 1 ) then
		    if (lwrite) print *,'s2 ',-diff/(2*slam*xp)*log(s2)
		    s2 = -slam*(s2a+slam)/(2*xm1*xm2)
		    s2 = -diff/(2*slam*xp)*dfflo1(s2,ier)
		else
		    ier = ier + 50
		    print *,'ffxdb0: untested: s2 better in first try'
		    if (lwrite) print *,'s2 ',+diff/(2*slam*xp)*log(s2)
		    s2 = +slam*(s2a-slam)/(2*xm1*xm2)
		    s2 = +diff/(2*slam*xp)*dfflo1(s2,ier)
		endif
		if ( lwrite ) print *,'s2+ ',s2,jsign
	    else
		s2 = -diff/(2*slam*xp)*log(s2)
		if ( jsign  ==  -1 ) s2 = -s2
	    endif
	    s1 = -dm1m2*xlogmm/(2*xp)
	    xx = s1+s2-1
	    if (lwrite) then
		print *,'ffxdbp: lam>0, first try, xx  = ',xx,s1,s2,-1
	    endif
c--#]	    first try:
	    if ( abs(xx)  <  xloss**2*max(abs(s1),abs(s2)) ) then
c--#[		second try:
c		this is unacceptable, try a better solution
		s1a = diff + slam*dm1m2
		if (lwrite) print *,'s1 = ',-s1a/(2*xp*slam),diff/
     &			(2*xp*slam)
		if ( abs(s1a)  >  xloss*diff ) then
c		    this works
		    s1 = -s1a/(2*xp*slam)
		else
c		    by division a more accurate form can be found
		    s1 = -2*xm1*xm2*xp/(slam*(diff - slam*dm1m2))
		    if (lwrite) print *,'s1+= ',s1
		endif
		s = s1
		s1 = s1*xlogmm
		if ( abs(xp)  <  xm2 ) then
		    s2a = xp - dm1m2
		else
		    s2a = xm2 - dm1p
		endif
		s2 = s2a - slam
		if (lwrite) print *,'s2 = ',s2/(2*xm2),slam/(2*xm2)
		if ( abs(s2)  >  xloss*slam ) then
c		    at least reasonable
		    s2 = s2 / (2*xm2)
		else
c		    division again
		    s2 = (2*xp) / (s2a+slam)
		    if (lwrite) print *,'s2+= ',s2
		endif
		if ( abs(s2)  <  .1 ) then
c		    choose a quick way to get the logarithm
		    s2 = dfflo1(s2,ier)
		elseif ( s2 == 1 ) then
		    print *,'ffxdbp: error: arg log would be 0!'
		    print *,'        xp,xma,xmb = ',xp,xma,xmb
		    goto 600
		else
		    h = abs(1-s2)
		    s2 = zxfflg(h,0,c0,ier)
		endif
		s2 = -diff/(slam*xp)*s2
		xx = s1 + s2 - 1
		if (lwrite) then
		    print *,'ffxdbp: lam>0, 2nd try, xx  = ',xx,s1,s2,-1
		endif
c--#]		second try:
		if ( abs(xx)  <  xloss**2*max(abs(s1),abs(s2)) ) then
c--#[		    third try:
c		    (we accept two times xloss because that's the same
c		    as in this try)
c		    A Taylor expansion might work.  We expand
c		    inside the logs. Only do the necessary work.

c		#[ split up 1:
		    xnoe = s2a+slam
		    a = 1
		    b = 2/xnoe-1/xp
		    c = -4/(xp*xnoe)
		    d = sqrt((2/xnoe)**2 + 1/xp**2)
		    call ffroot(d1,d2,a,b,c,d,ier)
		    if ( xp > 0 ) then
			beta = d2
		    else
			beta = d1
		    endif
		    alpha = beta*diff/slam
		    alph1 = 1-alpha
		    if ( alph1  <  xloss ) then
			s1a = 4*xp**2*xm1*xm2/(slam*dm1m2*(diff-slam*
     &				dm1m2))
			s1b = -diff/slam*4*xm1*xp/(dm1m2*xnoe*(2*xp-
     &				xnoe))
			b = -1/xp
			c = -(2/xnoe)**2
			call ffroot(d1,d2,a,b,c,d,ier)
			if ( xp > 0 ) then
			    betm2n = d2
			else
			    betm2n = d1
			endif
			d1 = s1a + s1b - diff/slam*betm2n
			if ( lwrite ) then
			    print *,'alph1    = ',d1,s1a,s1b,-diff/slam*
     &				betm2n
			    print *,'verg       ',1-alpha
			endif
			xmax = max(abs(s1a),abs(s1b))
			if ( xmax  <  1 ) then
			    alph1 = d1
			else
			    xmax = 1
			endif
			if ( lwarn .and. abs(alph1) < xloss*xmax ) then
			    call ffwarn(222,ier,alph1,xmax)
			    if ( lwrite ) print *,'d1,s1a,s2b,... = ',
     &				d1,s1a,s1b,diff/slam*betm2n
			endif
		    else
			betm2n = beta - 2/xnoe
		    endif
		    if ( lwrite ) then
			print *,'     s1 - alph1 = ',s1-alph1
			print *,'     s2 - alpha = ',s2-alpha
		    endif
c		#] split up 1:
c		#[ s2:

c		    first s2:

  490		    continue
		    s2p = s2 - alpha
		    if ( abs(s2p)  <  xloss*abs(s2) ) then
c -#[			bounds:
c			determine the boundaries for 1,5,10,15 terms
			if ( xprcn5  /=  precx ) then
			    xprcn5 = precx
			    bdn501 = ffbnd(3,1,xinfac)
			    bdn505 = ffbnd(3,5,xinfac)
			    bdn510 = ffbnd(3,10,xinfac)
			    bdn515 = ffbnd(3,15,xinfac)
			endif
c -#]			bounds:
			x = beta*xp
			ax = abs(x)
			if ( lwarn .and. ax  >  bdn515 ) then
c			    do not do the Taylor expansion
			    call ffwarn(23,ier,s2p,s2)
			    goto 495
			endif
			if ( ax  >  bdn510 ) then
			    s2a = x*(xinfac(13) + x*(xinfac(14) + x*(
     &				     xinfac(15) + x*(xinfac(16) + x*(
     &				     xinfac(17))))))
			else
			    s2a = 0
			endif
			if ( ax  >  bdn505 ) then
			    s2a = x*(xinfac(8) + x*(xinfac(9) + x*(
     &				    xinfac(10) + x*(xinfac(11) + x*(
     &				    xinfac(12) + s2a)))))
			endif
			if ( ax  >  bdn501 ) then
			    s2a = x*(xinfac(4) + x*(xinfac(5) + x*(
     &				     xinfac(6) + x*(xinfac(7) + s2a))))
			endif
			s2a = x**3*(xinfac(3)+s2a)
			s2b = 2*xp/xnoe*(s2a + x**2/2)
			s2p = s2b - s2a
			if ( lwarn .and. abs(s2p) < xloss*abs(s2a) )
     &				call ffwarn(24,ier,s2p,s2a)
			s2p = -diff/(xp*slam)*dfflo1(s2p,ier)
			if (lwrite) then
			    print *,'ffxdbp: Taylor expansion of s2-a'
			    print *,'	     in x = ',x
			    print *,'	     gives s2p = ',s2p
			endif
		    endif
c		#] s2:
c		#[ s1:

c		    next s1:

  495		    continue
		    s1p = s1 - alph1
		    if ( abs(s1p)  <  xloss*abs(s1) ) then
c -#[			bounds:
c			determine the boundaries for 1,5,10,15 terms
			if ( xprcn3  /=  precx ) then
			    xprcn3 = precx
			    bdn301 = ffbnd(3,1,xinfac)
			    bdn305 = ffbnd(3,5,xinfac)
			    bdn310 = ffbnd(3,10,xinfac)
			    bdn315 = ffbnd(3,15,xinfac)
			endif
c -#]			bounds:

			x = slam*(diff-slam*dm1m2)*alph1/(2*xp*xm1*xm2)
			h = (2*xp*(xm1+xm2) - xp**2)/(slam-dm1m2)
			ax = abs(x)
			if ( lwarn .and. ax  >  bdn315 ) then
c			    do not do the Taylor expansion
			    call ffwarn(21,ier,s1p,s1)
			    goto 500
			endif

c			see form job gets1.frm

			s1b = diff*(diff-slam*dm1m2)*betm2n/(2*xp*xm1*
     &				xm2)
			s1c = 1/(xm1*xnoe*(2*xp-xnoe))*(
     &				xp*( 4*xp*xm2 + 2*dm1m2**2/xm2*(xp-h) +
     &				2*dm1m2*(3*xp-h) - 8*dm1m2**2 )
     &				- 2*dm1m2**3/xm2*(3*xp-h)
     &				+ 4*dm1m2**4/xm2
     &				)
			if ( lwrite ) then
			    print *,'s1c was ',-2*xp/dm1m2 + 2*diff*
     &				(diff-slam*dm1m2)/(xm2*dm1m2*xnoe*(2*xp-
     &				xnoe)) + dm1m2/xm1
			    print *,'  en is ',s1c
			    print *,'s1b+s1c was ',dm1m2/xm1-x
			    print *,'      en is ',s1b+s1c
			endif
			s1d = x*dm1m2/xm1
			s1e = -x**2/2
			if ( ax  >  bdn310 ) then
			    s1a = x*(xinfac(13) + x*(xinfac(14) + x*(
     &				     xinfac(15) + x*(xinfac(16) + x*(
     &				     xinfac(17))))))
			else
			    s1a = 0
			endif
			if ( ax  >  bdn305 ) then
			    s1a = x*(xinfac(8) + x*(xinfac(9) + x*(
     &				   xinfac(10) + x*(xinfac(11) + x*(
     &				   xinfac(12) + s1a)))))
			endif
			if ( ax  >  bdn301 ) then
			    s1a = x*(xinfac(4) + x*(xinfac(5) + x*(
     &				     xinfac(6) + x*(xinfac(7) + s1a))))
			endif
			s1a = -x**3 *(xinfac(3) + s1a)
			s1f = dm1m2/xm1*(x**2/2 - s1a)
			s1p = s1e + s1d + s1c + s1b + s1a + s1f
			xmax = max(abs(s1a),abs(s1b),abs(s1c),abs(s1d),
     &				abs(s1e))
			if ( lwarn .and. abs(s1p) < xloss*xmax ) then
			    call ffwarn(223,ier,s1p,xmax)
			    if ( lwrite )
     &				print *,'s1p,s1e,s1d,s1c,s1b,s1a,s1f = '
     &				,s1p,s1e,s1d,s1c,s1b,s1a,s1f
			endif
			s1p = s*dfflo1(s1p,ier)
			if (lwrite) then
			    print *,'s1a = ',s1a
			    print *,'s1b = ',s1b
			    print *,'s1c = ',s1c
			    print *,'s1d = ',s1d
			    print *,'s1e = ',s1e
			    print *,'s1f = ',s1f
			    print *,'s   = ',s
			    print *,'ffxdbp: Taylor exp. of s1-(1-a)'
			    print *,'        in x = ',x
			    print *,'        gives s1p = ',s1p
			    print *,'        verg        ',s*log(xm2/xm1
     &				*exp(x))
			endif
		    endif
c		#] s1:

c		    finally ...

  500		    continue
		    xx = s1p + s2p
		    if ( lwarn .and. abs(xx)  <  xloss*abs(s1p) ) then
			call ffwarn(25,ier,xx,s1p)
		    endif
c--#]		    third try:
		endif
	    endif
  600	    continue
	    if ( xp  >  xm1+xm2 ) then
c--#[		imaginary part:
c		in this case ( xlam>0, so xp>(m1+m2)^2) ) there also
c		is an imaginary part
		y = -pi*diff/(slam*xp)
	    else
		y = 0
c--#]		imaginary part:
	    endif
	 else
c	    the root is complex (k^2 between -(m1+m2)^2 and -(m2-m1)^2)
c--#[	    first try:
	    slam = sqrt(-xlam)
	    xnoe = dm2p + xm1
	    s1 = -(dm1m2/(2*xp))*xlogmm
	    s2 = -diff/(slam*xp)*atan2(slam,xnoe)
	    xx = s1 + s2 - 1
	    if (lwrite) then
		print *,'ffxdbp: lam<0, first try, xx  = ',xx,s1,s2,-1
c		alpha = -xlam/(2*xp*xnoe)
c		alph1 = -(xp**2-dm1m2**2)/(2*xp*xnoe)
c		print *,'	 alpha = ',alpha
c		print *,'	 s1 = ',s1,' - 2alph1 = ',s1-2*alph1
c		print *,'	 s2 = ',s2,' - 2alpha = ',s2-2*alpha
	    endif
c--#]	    first try:
	    if ( lwarn .and. abs(xx) < xloss**2*max(abs(s1),abs(s2)) )
     &			then
		call ffwarn(224,ier,xx,max(abs(s1),abs(s2)))
	    endif
	    y = 0
	endif
  590	continue
	cdb0p = DCMPLX(DBLE(xx),DBLE(y))
	cdb0 = cdb0p*(1/DBLE(xp))
	goto 990
c -#]	normal case:
c  #] unequal nonzero masses:
c  #[ debug:
  990	continue
	if (lwrite) then
	    print *,'cdb0  = ',cdb0,cdb0p
	endif
c  #] debug:
c###] ffxdbp:
	end
