!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine genVphoton(r,p,wt,*)
      implicit none
      include 'types.f'
c---generate two particle phase space and x1,x2 integration
c---p1+p2 --> V(p3+p4+p5)
c----
c---- with p5 generated using tau as a variable of integration,
c---- with minimum value taucut

      include 'constants.f'
      include 'mxpart.f'
      include 'limits.f'
      include 'vegas_common.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'debug.f'
      include 'kprocess.f'
      include 'nwz.f'
      include 'leptcuts.f'
      integer ipick
      real(dp):: r(mxdim),p(mxpart,4),Q(4),p3(4),p4(4),p5(4),p45(4)
      real(dp):: wbw,wt,wt345,wt45,Qsq,rtshat,Qsqmin,Qsqmax
      real(dp), parameter:: wt0=one/twopi

      p(:,:)=zip

      wt=zip

      Qsqmin=max(wsqmin,one) ! ensure minimum value of m(34)>1 GeV
      Qsqmax=min(wsqmax,sqrts**2*0.9999d0) ! ensure maximum value a bit below s
c--- generate invariant mass of Q=p3+p4+p5
      if (n3==0) then
         wbw=one
         call pick(2,Qsq,Qsqmin,Qsqmax,r(1),wbw)
      elseif (n3==1) then
         if (gammptmin < mass3) then
           call breitw(r(1),Qsqmin,Qsqmax,mass3,width3,Qsq,wbw)
         else ! artificially inflate width to make sampling easier
           call breitw(r(1),Qsqmin,Qsqmax,mass3,100._dp*width3,Qsq,wbw)
         endif
      endif

      rtshat=sqrt(Qsq)

c--- generate pa+pb -> Q
      call genQ(rtshat,r(2),1.e-8_dp,p,wt)

c      write(6,*) 'rtshat=',rtshat
c      write(6,*) 'p1',p(1,4),p(1,1),p(1,2),p(1,3)
c      write(6,*) 'p2',p(2,4),p(2,1),p(2,2),p(2,3)
c      write(6,*) 'p3',p(3,4),p(3,1),p(3,2),p(3,3)
c      write(6,*) 'p4',p(4,4),p(4,1),p(4,2),p(4,3)

c--- now decay Q into 3,4,5
      Q=p(1,:)+p(2,:)
      call phi1_2m_nobw(zip,r(3),r(4),r(5),zip,Q,p3,p45,wt345,*999)
      call phi3m0(r(6),r(7),p45,p4,p5,wt45,*999)
      wt=wt0*wbw*wt*wt345*wt45/twopi

c      write(6,*) 'ph1',ph(1,4),ph(1,1),ph(1,2),ph(1,3)
c      write(6,*) 'ph2',ph(2,4),ph(2,1),ph(2,2),ph(2,3)
c      write(6,*) 'ph3',ph(3,4),ph(3,1),ph(3,2),ph(3,3)
c      write(6,*) 'ph4',ph(4,4),ph(4,1),ph(4,2),ph(4,3)

c--- translate to momenta to be returned
      p(1,:)=-p(1,:)
      p(2,:)=-p(2,:)
      if ((kcase == kWgamma) .or. (kcase == kWgaj_a) .or. (kcase == kWga_ew)) then
        if (nwz == +1) then
          ipick=1 ! W+, charged lepton is p4
        else
          ipick=2 ! W-, charged lepton is p3
        endif
      endif
      if (kcase == kZgamma) then
c radiation from both leptons -- pick at random, noting that
c r(ndim+1) is a uniform random variable not adapted by VEGAS
        if (r(ndim+1) < 0.5_dp) then
          ipick=1
        else
          ipick=2
        endif
      endif
      if (ipick == 1) then
        p(3,:)=p3(:)
        p(4,:)=p4(:)
      else
        p(3,:)=p4(:)
        p(4,:)=p3(:)
      endif
      p(5,:)=p5(:)
      p(6,:)=zip

      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts

c--- the factor below will be added later, in lowint
      wt=wt*xx(1)*xx(2)*sqrts**2

c      call writeout(p)
c      pause

      if   ( (xx(1) > one) .or. (xx(2) > one)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in genVphoton',xx(1),xx(2)
        return 1
      endif

      return

  999 wt=zip
      return 1

      end
