!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function threemassbox(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq)
      implicit none
      include'types.f'
c     Three external offshellnesses
c     calculated for a limited order of arguments

c     $I_4^{D=4-2 \eps12lon}(0,\pd^2,\pt^2,\pq^2;s_{12},s_{23};0,0,0,0)$}

c----%\cite{Bern:1993kr}
c----\bibitem{Bern:1993kr}
c----  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c----  %``Dimensionally regulated pens23gon integrals,''
c----  Nucl.\ Phys.\ B {\bf 412}, 751 (1994)
c----  [arXiv:hep-ph/9306240].
c----  %%CIs23TION = HEP-PH 9306240;%%
c----  Eqs. (I.15)
c-----or from /hep-ph/0508308 v3 Eqn (A27)
c-----v3 corrects previous versions.
c                         [                   s12     p4sq ]
c                         [   0       0     - ---   - ---- ]
c                         [                    2       2   ]
c                         [                                ]
c                         [                   p2sq    s23  ]
c                         [   0       0     - ----  - ---  ]
c                         [                    2       2   ]
c                    y5 = [                                ]
c                         [   s12     p2sq            p3sq ]
c                         [ - ---   - ----    0     - ---- ]
c                         [    2       2               2   ]
c                         [                                ]
c                         [   p4sq    s23     p3sq         ]
c                         [ - ----  - ---   - ----    0    ]
c                         [    2       2       2           ]
      include 'constants.f'

      real(dp):: r,s12,s23,p1sq,p2sq,p3sq,p4sq,m1sq,m2sq,m3sq,m4sq
      complex(dp):: threemassbox,lnrat,fac,Li2(6),qlLi2omrat,qlLi2omx2
c      logical:: landau

c      s12=2._dp*Y(1,3)
c      s23=2._dp*Y(2,4)
c      p2sq=2._dp*Y(2,3)
c      p3sq=2._dp*Y(3,4)
c      p4sq=2._dp*Y(1,4)

       if (
     &    (p1sq  /=  0._dp) .or.
     &    (m1sq  /=  0._dp) .or.
     &    (m2sq  /=  0._dp) .or.
     &    (m3sq  /=  0._dp) .or.
     &    (m4sq  /=  0._dp) ) then
        write(6,*) 'Unimplemented configuration for threemassbox'
        write(6,*) p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq
        stop
      endif

      r=1._dp-p2sq*p4sq/(s12*s23)
c      write(6,*) 'threemassbox: r=',r

c--- May need to revisit this for special cases???

c     Use expansion only in cases where signs (s12,s23,p2sq,p4sq) are not
c     ++-- or --++
c      landau=((s12gn(1._dp,s12)  ==  s12gn(1._dp,s23))
c     & .and. (s12gn(1._dp,p2sq)  ==  s12gn(1._dp,p4sq))
c     & .and. (s12gn(1._dp,s12)  /=  s12gn(1._dp,p2sq)))
c      if ((abs(r)  <  1d-6) .and. (landau .eqv. .false.)) then
c---expanded case
c      Ires(-2)=czip
c      Ires(-1)=-cmplx((1._dp+0.5_dp*r)/(s12*s23),kind=dp)
c      Ires(0)=Ires(-1)*(lnrat(musq,s12)+lnrat(p3sq,s23)-cmplx(2._dp,kind=dp)
c     & -cmplx(1._dp+p4sq/s23,kind=dp)*qlL0(p4sq,s23))
c     & +cmplx(r/(s12*s23),kind=dp)*(qlL1(p4sq,s23)-qlL0(p4sq,s23)-cone)
c      else

c---General case
      fac=cmplx(1._dp/(s12*s23-p2sq*p4sq),kind=dp)

      Li2(1)=qlLi2omrat(-p2sq,-s12)
      Li2(2)=qlLi2omrat(-p4sq,-s23)
      Li2(3)=qlLi2omx2(-p2sq,-p4sq,-s12,-s23)

c--- the form below is the 6-d box that can be
c--- obtained by analogy with the usual hard box
      threemassbox=fac*(
     & -ctwo*(Li2(1)+Li2(2)-Li2(3))-lnrat(-s12,-s23)**2
     & +lnrat(-s12,-p3sq)*lnrat(-s12,-p4sq)
     & +lnrat(-s23,-p2sq)*lnrat(-s23,-p3sq))

c      endif
      return
      end

c Routines taken from QCDLoop (v1)
      double complex function qlLi2omrat(x,y)
      implicit none
      include 'types.f'
c     expression for dilog(1-(x-i*ep)/(y-i*ep)) for real x and y
c     Hence arguments are typically negative invariants
      real(dp):: x,y,omarg,arg,ddilog
      complex(dp):: lnrat,wlog
      real(dp), parameter:: pi=3.1415926535897932385_dp,pisq=pi*pi,pisqo6=pisq/6._dp
      omarg=x/y
      arg=1._dp-omarg
      if (arg  >  1._dp) then
      wlog=lnrat(x,y)
      qlLi2omrat=cmplx(pisqo6-ddilog(omarg),kind=dp)-log(arg)*wlog
      else
      qlLi2omrat=cmplx(ddilog(arg),kind=dp)
      endif
      return
      end

      double complex function qlLi2omx2(v,w,x,y)
      implicit none
      include 'types.f'
c     expression for dilog(1-(v-i*ep)*(w-i*ep)/(x-i*ep)/(y-i*ep))
c     for real v,w,x and y
      real(dp):: v,w,x,y,omarg,arg,ddilog
      complex(dp):: lnrat,lnarg,lnomarg,prod
      real(dp), parameter:: pi=3.1415926535897932385_dp,pisq=pi*pi,pisqo6=pisq/6._dp

      qlLi2omx2=(0._dp,0._dp)

      arg=(v*w)/(x*y)
      omarg=1._dp-arg

      if (arg  <=  1._dp) then
         if (arg  ==  0._dp .or. arg  ==  1._dp) then
            prod=0._dp
         else
            lnarg=lnrat(v,x)+lnrat(w,y)
            lnomarg=cmplx(log(omarg),kind=dp)
            prod=lnarg*lnomarg
         endif
         qlLi2omx2=cmplx(pisqo6-ddilog(arg),kind=dp)-prod
      elseif (arg  >  1._dp) then
         arg=(x*y)/(v*w)
         lnarg=-lnrat(v,x)-lnrat(w,y)
         lnomarg=cmplx(log(1._dp-arg),kind=dp)
         qlLi2omx2=-cmplx(pisqo6-ddilog(arg),kind=dp)+lnarg*lnomarg
     &   -0.5_dp*lnarg**2
      endif
      return
      end
