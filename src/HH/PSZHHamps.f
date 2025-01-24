!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine PSZHHamps(pa,pb,pc,pd,amp)
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'scalarselect.f'
c--    Formula taken from Plehn, Spira and Zerwas
c--    Nucl. Phys. B479 (1996) 46

      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      real(dp):: pa(4),pb(4),pc(4),pd(4),pcsq,pdsq,
     & ss,tt,uu,S,T,U,T1,T2,U1,U2,rhoc,rhod,mQsq,lambdaHHH
      complex(dp):: Ftriangle,Fbox,Gbox,
     & Dabc,Dbac,Dacb,Cab,Cbc,Cac,Ccd,Cbd,Cad,
     & Cdelta,Cbox,amp(2)

      mQsq=mt**2
      ss=(pa(4)+pb(4))**2
     &  -(pa(1)+pb(1))**2-(pa(2)+pb(2))**2-(pa(3)+pb(3))**2
      tt=(pa(4)+pc(4))**2
     &  -(pa(1)+pc(1))**2-(pa(2)+pc(2))**2-(pa(3)+pc(3))**2
      uu=(pb(4)+pc(4))**2
     &  -(pb(1)+pc(1))**2-(pb(2)+pc(2))**2-(pb(3)+pc(3))**2
      S=ss/mQsq
      T=tt/mQsq
      U=uu/mQsq
      lambdaHHH=3._dp*(hmass/zmass)**2
      rhoc=hmass**2/mQsq
      rhod=hmass**2/mQsq
      T1=T-rhoc
      U1=U-rhoc
      T2=T-rhod
      U2=U-rhod
c      ss=(pa+pb)**2
c      tt=(pa+pc)**2
c      uu=(pb+pc)**2
      pcsq=pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2
      pdsq=pd(4)**2-pd(1)**2-pd(2)**2-pd(3)**2
      Cab=loopI3(0._dp,0._dp,ss,mQsq,mQsq,mQsq,musq,0)
      Cbc=loopI3(0._dp,pcsq,uu,mQsq,mQsq,mQsq,musq,0)
      Cac=loopI3(0._dp,pcsq,tt,mQsq,mQsq,mQsq,musq,0)
      Ccd=loopI3(pcsq,pdsq,ss,mQsq,mQsq,mQsq,musq,0)
      Cbd=loopI3(0._dp,pdsq,tt,mQsq,mQsq,mQsq,musq,0)
      Cad=loopI3(0._dp,pdsq,uu,mQsq,mQsq,mQsq,musq,0)
      Dabc=loopI4(0._dp,0._dp,pcsq,pdsq,ss,uu,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dbac=loopI4(0._dp,0._dp,pcsq,pdsq,ss,tt,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dacb=loopI4(0._dp,pcsq,0._dp,pdsq,tt,uu,mQsq,mQsq,mQsq,mQsq,musq,0)


      Ftriangle=2._dp/S*(2._dp+(4._dp-S)*mQsq*Cab)

      Fbox=(4._dp*S+8._dp*S*mQsq*Cab-2._dp*S*(S+rhoc+rhod-8._dp)
     & *mQsq**2*(Dabc+Dbac+Dacb)
     & +(rhoc+rhod-8._dp)*mQsq*(T1*Cac+U1*Cbc+U2*Cad+T2*Cbd
     & -(T*U-rhoc*rhod)*mQsq*Dacb))/S**2
      Gbox=((T**2+rhoc*rhod-8._dp*T)*mQsq
     & *(S*Cab+T1*Cac+T2*Cbd-S*T*mQsq*Dbac)
     & +(U**2+rhoc*rhod-8._dp*U)*mQsq
     & *(S*Cab+U1*Cbc+U2*Cad-S*U*mQsq*Dabc)
     & -(T**2+U**2-2._dp*rhoc*rhod)*(T+U-8._dp)*mQsq*Ccd
     & -2._dp*(T+U-8._dp)*(T*U-rhoc*rhod)*mQsq**2
     & *(Dabc+Dbac+Dacb))/(S*(T*U-rhoc*rhod))

      Cdelta=lambdaHHH*zmass**2/cmplx(ss-hmass**2,hmass*hwidth,kind=dp)
      Cbox=cone

      amp(1)=Cdelta*Ftriangle+Cbox*Fbox
      amp(2)=Cbox*Gbox
      return
      end
