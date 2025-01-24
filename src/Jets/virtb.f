!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function virtb(s,t,u)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'masses.f'
c     Virtual corrections taken from Ellis and Sexton
c  %\cite{Ellis:1985er}
c  \bibitem{Ellis:1985er}
c   R.~K.~Ellis and J.~C.~Sexton,
c  %``QCD Radiative Corrections To Parton Parton Scattering,''
c  Nucl.\ Phys.\ B {\bf 269}, 445 (1986).
c  %%CITATION = NUPHA,B269,445;%%
c     In units of
c     as/2/pi*(4 \pi)^epsilon*Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)
c      \approx as/2/pi*(4 \pi)^epsilon /Gamma(1-ep)
      real(dp):: virtb,a4,b4,theta,x,s,t,u,Ws,Wt,Wu,Ws2,Wt2,Wu2,Wmu,q2,TfNf
      complex(dp):: A6texact
      a4(s,t,u) =2d0*v*(s**2+u**2)/t**2
      b4(s,t,u) =-8d0*cf*(s**2/u/t)
      theta(x)=half+sign(half,x)

      TfNf=TR*dfloat(nf)
      q2=musq
      Ws=log(abs(s/q2))
      Wt=log(abs(t/q2))
      Wu=log(abs(u/q2))
      Wmu=log(abs(musq/q2))
      Ws2=Ws**2-theta(s)*pisq
      Wt2=Wt**2-theta(t)*pisq
      Wu2=Wu**2-theta(u)*pisq


      virtb=
     & +(CF*(-4d0*EPINV*EPINV2-EPINV*(6d0+8d0*WS-8d0*WU-4d0*WT))
     & +xn*EPINV*(4d0*WS-2d0*WU-2d0*WT))*A4(S,T,U)
     & +(CF*(-16d0-2d0*WT2+WT*(6d0+8d0*WS-8d0*WU))
     & +xn*(85d0/9d0+PISQ+2d0*WT*(WT+WU-2d0*WS)
     & +11d0/3d0*(WMU-WT))
     & +TfNf*(4d0/3d0*(WT-WMU)-20d0/9d0))*A4(S,T,U)
     & -4d0*V*CF*(S**2-U**2)/T**2
     & *(2d0*PISQ+2d0*WT2+WS2+WU2-2d0*WS*WT-2d0*WT*WU)
     & +xn*V*((S**2-U**2)/T**2
     & *(3d0*PISQ+3d0*WT2+2d0*WS2+WU2-4d0*WS*WT-2d0*WT*WU))
     & +4d0*V*(CF*((WU-WS)-(U-S)/T*(2d0*WT-WS-WU))
     & +xn*(-0.5d0*S/T*(WT-WU)+U/T*(WT-WS)))

      virtb=virtb
     & +(CF*(-4d0*EPINV*EPINV2-EPINV*(6d0+8d0*WS-8d0*WT-4d0*WU))
     & +xn*EPINV*(4d0*WS-2d0*WT-2d0*WU))*A4(S,U,T)
     & +(CF*(-16d0-2d0*WU2+WU*(6d0+8d0*WS-8d0*WT))
     & +xn*(85d0/9d0+PISQ+2d0*WU*(WU+WT-2d0*WS)
     & +11d0/3d0*(WMU-WU))
     & +TfNf*(4d0/3d0*(WU-WMU)-20d0/9d0))*A4(S,U,T)
     & -4d0*V*CF*(S**2-T**2)/U**2
     & *(2d0*PISQ+2d0*WU2+WS2+WT2-2d0*WS*WU-2d0*WU*WT)
     & +xn*V*((S**2-T**2)/U**2
     & *(3d0*PISQ+3d0*WU2+2d0*WS2+WT2-4d0*WS*WU-2d0*WU*WT))
     & +4d0*V*(CF*((WT-WS)-(T-S)/U*(2d0*WU-WS-WT))
     & +xn*(-0.5d0*S/U*(WU-WT)+T/U*(WU-WS)))

      virtb=virtb
     & +CF*(-4d0*EPINV*EPINV2-EPINV*(6d0+4d0*WS-4d0*WT-4d0*WU))*B4(S,T,U)
     & +xn*(2d0*EPINV*(2d0*WS-WT-WU))*B4(S,T,U)
     & +CF*(-16d0-3d0/2d0*(WT2+WU2+2d0*WT*WU)
     & +2d0*WS*(WU+WT)+2d0*(WT+WU)-0.5d0*PISQ)*B4(S,T,U)
     & +xn*(85d0/9d0+5d0/4d0*(WT2+WU2+2d0*WT*WU)
     & -WT*WU-2d0*WS*(WT+WU)-4d0/3d0*(WT+WU)
     & +5d0/4d0*PISQ+11d0/3d0*WMU)*B4(S,T,U)
     & +TfNf*(-20d0/9d0+2d0/3d0*(WT+WU-2d0*WMU))*B4(S,T,U)
     & +((PISQ+WT2+WU2-2d0*WT*WU)*U*T/2d0/S**2
     & +U/S*WT+T/S*WU)*B4(S,T,U)/xn

       virtb=virtb-6d0*PISQ*theta(t)*A4(S,T,U)

c--- Contribution of a top-quark loop
      virtb=virtb+real(A6texact(t,mt**2),dp)*(A4(S,T,U)+half*B4(S,T,U))
     &           +real(A6texact(u,mt**2),dp)*(A4(S,U,T)+half*B4(S,U,T))

      return
      end
