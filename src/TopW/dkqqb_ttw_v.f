!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine dkqqb_ttw_v(pin,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     April, 2012.                                                     *
c     calculate the Interference matrix element squared                *
c     to one-loop in semileptonic decay of the top and anti-top        *
c     for the process                                                  *
c     My notation                                                      *
c     q(-p1) +qbar(-p2) -> t(nu(p3)+e+(p4)+b(p5))                      *
c                         +tbar(bbar(p6)+e-(p7)+nubar(p8))             *
c                         +W(l(p9)+a(p10))                             *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'plabel.f'

      integer:: j,k,nu,j1,j2,hb,hc
      real(dp):: msq(-nf:nf,-nf:nf),pin(mxpart,4),p(mxpart,4),
     & s34,fac,qqb,qbq,q(mxpart,4),
     & beta,betasq,bp,bm,s56,dot
      complex(dp):: mqqbprod(2,2),mqbqprod(2,2),prop,
     & mtotqqb(2,2),mtotqbq(2,2),mtotqqbv(2,2),mtotqbqv(2,2),
     & mtop(2,2),manti(2,2),mtopv(2,2),mantiv(2,2)

c-----Zero matrix element
      msq(:,:)=0._dp

c--setup amplitudes for top decay
      call tdecayrod(pin,3,4,5,6,7,8,0,mtop)
      call adecayrod(pin,3,4,5,6,7,8,0,manti)
      call tdecayrod_v(pin,3,4,5,6,7,8,mtopv)
      call adecayrod_v(pin,3,4,5,6,7,8,mantiv)

      do nu=1,4
      p(1,nu)=pin(1,nu)
      p(2,nu)=pin(2,nu)
      p(3,nu)=pin(9,nu)
      p(4,nu)=pin(10,nu)
      p(5,nu)=pin(3,nu)+pin(4,nu)+pin(5,nu)
      p(6,nu)=pin(6,nu)+pin(7,nu)+pin(8,nu)
      enddo

c calculate betap
      s56=2._dp*mt**2+2._dp*dot(p,5,6)
      betasq=1._dp-4._dp*mt**2/s56
      if (betasq >= 0._dp) then
        beta=sqrt(betasq)
        bp=0.5_dp*(1._dp+beta)
        bm=1._dp-bp
      else
        write(6,*) 'betasq < 0 in qqb_ttw.f, betasq=',betasq
        call flush(6)
        stop
      endif


c---create rodrigo momenta
      do j=1,6
      do nu=1,4
      if (j==5) then
      q(j,nu)=(bp/beta)*p(5,nu)-(bm/beta)*p(6,nu)
      elseif (j==6) then
      q(j,nu)=(bp/beta)*p(6,nu)-(bm/beta)*p(5,nu)
      else
      q(j,nu)=p(j,nu)
      endif
      enddo
      enddo

      s34=2._dp*dot(q,3,4)
      call spinoru(6,q,za,zb)
      call ttWprod(1,2,mqqbprod)
      call ttWprod(2,1,mqbqprod)

      do hb=1,2
      do hc=1,2
      mtotqqb(hb,hc)=czip
      mtotqbq(hb,hc)=czip
      mtotqqbv(hb,hc)=czip
      mtotqbqv(hb,hc)=czip
      do j1=1,2
      do j2=1,2
      mtotqqb(hb,hc)=mtotqqb(hb,hc)+
     & mtop(hb,j1)*mqqbprod(j1,j2)*manti(j2,hc)
      mtotqbq(hb,hc)=mtotqbq(hb,hc)+
     & mtop(hb,j1)*mqbqprod(j1,j2)*manti(j2,hc)
      mtotqqbv(hb,hc)=mtotqqbv(hb,hc)
     & +mtopv(hb,j1)*mqqbprod(j1,j2)*manti(j2,hc)
     & +mtop(hb,j1)*mqqbprod(j1,j2)*mantiv(j2,hc)
      mtotqbqv(hb,hc)=mtotqbqv(hb,hc)
     & +mtopv(hb,j1)*mqbqprod(j1,j2)*manti(j2,hc)
     & +mtop(hb,j1)*mqbqprod(j1,j2)*mantiv(j2,hc)
      enddo
      enddo
      enddo
      enddo

      prop=s34/cplx2(s34-wmass**2,wmass*wwidth)
     & /cplx2(zip,mt*twidth)**2
      fac=CF*ason2pi
      fac=V*gwsq**6*gsq**2*abs(prop)**2*fac

c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

      qqb=0._dp
      qbq=0._dp
      do hb=1,2
      do hc=1,2
      qqb=qqb+fac*aveqq*real(mtotqqb(hb,hc)*conjg(mtotqqbv(hb,hc)))
      qbq=qbq+fac*aveqq*real(mtotqbq(hb,hc)*conjg(mtotqbqv(hb,hc)))
      enddo
      enddo


c----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
          if ((j > 0) .and. (k < 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo
      return
      end
