!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gen4(r,p,wt4,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'kprocess.f'
      include 'phasemin.f'
      include 'x1x2.f'
      include 'ipsgen.f'
      include 'breit.f'
      include 'zerowidth.f'
      include 'limits.f'
      include 'leptcuts.f'
      integer:: nu,icount
      real(dp):: r(mxdim)
      real(dp):: wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p(mxpart,4)
      real(dp):: pswt,xjac,taubound
      real(dp):: tau,x1mx2,surd
      real(dp):: wtprop,zaprop,rmass,rwidth,s34
      include 'energy.f'
      data icount/1/
      save icount
!$omp threadprivate(icount)

c--- statement function
      wtprop(s34,rmass,rwidth)=(s34-rmass**2)**2+(rmass*rwidth)**2
      zaprop(s34,rmass,rwidth)=s34*((s34-rmass**2)**2+(rmass*rwidth)**2)

      wt4=0._dp

      taubound=zip
      if (n2 == 1) then
        if (zerowidth) then
          taubound=taubound+mass2
        else
          taubound=taubound+sqrt(wsqmin)
        endif
      endif
      if (n3 == 1) then
        if (zerowidth) then
          taubound=taubound+mass3
        else
          taubound=taubound+sqrt(bbsqmin)
        endif
      endif
      if ((n2 == 1) .and. (n3 == 1)) then
        taubound=max(taubound,m3456min)
      endif
      if  ((kcase == kWgajet) .or. (kcase == kZgajet) .or. (kcase == kWgaj_a)
     &.or. (kcase == kWga_ew) .or. (kcase==kWgajew)) then
        taubound=taubound+gammptmin
      endif

      taubound=max((taubound/sqrts)**2,taumin)

      xjac=1._dp
      call pick(2,tau,taubound,1._dp,r(9),xjac)

c      lntaum=log(taumin)
c      tau=exp(lntaum*(one-r(9)))
c      xjac=-lntaum*tau

c      tau=(one-taumin)*r(9)**2+taumin
c      xjac=2*r(9)*(one-taumin)

      x1mx2=two*r(10)-one
      surd=sqrt(x1mx2**2+four*tau)

      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)
c      write(*,*) r(1),r(9),r(10)
c      write(*,*) xx(1),xx(2)

      xjac=xjac*two/surd

c--- for comparison with C. Oleari's e+e- --> QQbg calculation
c      if (runstring(1:5) == 'carlo') then
c        xx(1)=1._dp
c        xx(2)=1._dp
c        xjac=1._dp
c      endif

      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      if  (    (kcase==kt_bbar)
     &    .or. (kcase==kbq_tpq)) then
        call phase41(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      elseif ( (kcase==kqqttbb)
     &    .or. (kcase==kqqttgg))  then
        call phase4m(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      elseif ( (kcase==kZgamma) .or. (kcase==kWgamma)
     &    .or. (kcase==kZgajet) .or. (kcase==kWgajet) .or. (kcase==kWgajew)
     &    .or. (kcase==kWgaj_a) .or. (kcase==kWga_ew))  then
        if (ipsgen == 1) then
          call phase4Vgam(r,p1,p2,p3,p4,p5,p6,pswt,*999)
        elseif (ipsgen == 2) then
          call phase4Vgam_fsr(r,p1,p2,p3,p4,p5,p6,pswt,*999)
        else
            print *, "unsupported ipsgen ..."
            call exit(1)
        endif

      elseif (kcase==kvlchk4)  then
        call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      else
        call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0._dp
      enddo

c      if (interference) then
c        if (icount == 1) then
c          bw34_56=.true.
c          icount=icount-1
c        else
c          bw34_56=.false.
c          do nu=1,4
c            p(4,nu)=p6(nu)
c            p(6,nu)=p4(nu)
c          enddo
c          icount=icount+1
c        endif
c      endif

      wt4=xjac*pswt
      if (debug) write(6,*) 'wt4 in gen4',wt4
      return

 999  return 1
      end

