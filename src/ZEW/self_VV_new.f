!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c -- formulism in the following subroutines and functions can be found in
c -- (B.1-4,6) [Fortshr. Phys. 38 (1990) 3 Hollik]
c*************************************************************************
c -- self energy
      subroutine znen(slf,alpha,sw2,psq,ep,flag)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'masses.f'
      include 'constants.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'scalarselect.f'
      character(len=2):: flag
      real(dp):: psq,mzsq,mwsq,mhsq,sw2,cw2,djitieD,alpha,
     &     Qf(9),gfp(9),gfm(9),fmsq(9),vf(9),af(9)
      complex(dp):: slf,wphi,fermi,cihfenF
      integer ep,j,k
      external djitieD,cihfenF
      include 'cplx.h'
c      complex(dp):: dfermi,fermi1,fermi2
c      common/em/aemmz
c      common/lepton/mme,mmmu,mmtau
c      common/quark/mup,mdwn,mchrm,mstrg,mbot,mtop

      mzsq = zmass**2
      mwsq = wmass**2
      mhsq = hmass**2

c      sw2 = xw
      cw2 = 1._dp - sw2

      Qf(1:3) = -1._dp
      Qf(4:8) = Q(1:5)*sqrt(xn)
      Qf(9) = Q(2)*sqrt(xn)

      gfp(1:3) = re
      gfm(1:3) = le
      gfp(4:8) = r(1:5)*sqrt(xn)
      gfm(4:8) = l(1:5)*sqrt(xn)
      gfp(9) = r(2)*sqrt(xn)
      gfm(9) = l(2)*sqrt(xn)

      vf = (gfm + gfp)/2._dp
      af = (gfm - gfp)/2._dp

c      print*, 'vf, d: ',
c     &     vf(4)/sqrt(xn), Qf(4)/sqrt(xn),xn,
c     &     (-0.5_dp - 2._dp*sw2*Qf(4)/sqrt(xn))/2._dp/sqrt(sw2*cw2)
c      print*, 'vf, u: ',
c     &     vf(5)/sqrt(xn), Qf(5)/sqrt(xn),
c     &     (0.5_dp - 2._dp*sw2*Qf(5)/sqrt(xn))/2._dp/sqrt(sw2*cw2)

c      print*, 'af, u: ',
c     &     af(5)/sqrt(xn), Qf(5)/sqrt(xn),
c     &     1._dp/4._dp/sqrt(sw2*cw2)

      fmsq(1) = mel**2
      fmsq(2) = mmu**2
      fmsq(3) = mtau**2
      fmsq(4) = md**2
      fmsq(5) = mu**2
      fmsq(6) = ms**2
      fmsq(7) = mc**2
      fmsq(8) = mb**2
      fmsq(9) = mt**2


      fermi = cplx2(0._dp, 0._dp)
c      fermi1 = cplx2(0._dp, 0._dp)
c      fermi2 = cplx2(0._dp, 0._dp)
c      dfermi = cplx2(0._dp, 0._dp)
      if(ep == 0) then
         if(flag == 'AA') then
            do j = 1, 9
               fermi = fermi + 4._dp/3._dp*Qf(j)**2*(
     &              + psq*djitieD(fmsq(j),ep)
     &              + (psq + 2._dp*fmsq(j))*cihfenF(psq,fmsq(j),fmsq(j))
     &              - psq/3._dp)
c               dfermi = dfermi + 4._dp/3._dp*Qf(j)**2*(
c     &              + djitieD(fmsq(j),ep)
c     &              + cihfenF(0._dp,fmsq(j),fmsq(j))
c     &              + (2._dp*fmsq(j))*ccdb0(0._dp,fmsq(j),fmsq(j))
c     &              - 1._dp/3._dp)
c               print*, cihfenF(psq,fmsq(j),fmsq(j)),
c     &              loopI2(psq,fmsq(j),fmsq(j),musq,ep)
c     &              -loopI2(0._dp,fmsq(j),fmsq(j),musq,ep)
c               pause
            end do
            slf = fermi - 3._dp*psq*djitieD(mwsq,ep)
     &           - (3._dp*psq + 4._dp*mwsq)*cihfenF(psq,mwsq,mwsq)
            slf = alpha/4._dp/pi*slf
c            dslf = dfermi - 3._dp*djitieD(mwsq,ep)
c     &           - 3._dp*cihfenF(0._dp,mwsq,mwsq)
c     &           - 4._dp*mwsq*ccdb0(0._dp,mwsq,mwsq)
c            dslf = alpha/4._dp/pi*dslf
c            ct = -psq*real(dslf,dp)
         else if(flag == 'AZ') then
            do j = 1, 9
               fermi = fermi - 4._dp/3._dp*Qf(j)*vf(j)*(
     &              + psq*djitieD(fmsq(j),ep)
     &              + (psq + 2._dp*fmsq(j))*cihfenF(psq,fmsq(j),fmsq(j))
     &              - psq/3._dp)
c               fermi1 = fermi1 - 4._dp/3._dp*Qf(j)*vf(j)*(
c     &              + mzsq*djitieD(fmsq(j),ep)
c     &              + 3._dp*mzsq*cihfenF(mzsq,fmsq(j),fmsq(j))
c     &              - mzsq/3._dp)
c               fermi2 = fermi2 - 4._dp/3._dp*Qf(j)*vf(j)*(
c     &              + 2._dp*fmsq(j)*cihfenF(0._dp,fmsq(j),fmsq(j)))
            end do
            wphi = djitieD(mwsq,ep)*((3._dp*cw2 + 1._dp/6._dp)*psq
     &           + 2._dp*mwsq)/sqrt(cw2*sw2)
     &           + cihfenF(psq,mwsq,mwsq)*((3._dp*cw2 + 1._dp/6._dp)*psq
     &           + (4._dp*cw2 + 4._dp/3._dp)*mwsq)/sqrt(cw2*sw2)
     &           + psq/9._dp/sqrt(cw2*sw2)

c            wphi1 = djitieD(mwsq,ep)*((3._dp*cw2 + 1._dp/6._dp)*mzsq
c     &           + 2._dp*mwsq)/sqrt(cw2*sw2)
c     &           + cihfenF(mzsq,mwsq,mwsq)*((3._dp*cw2
c     &           + 1._dp/6._dp)*mzsq
c     &           + (4._dp*cw2 + 4._dp/3._dp)*mwsq)/sqrt(cw2*sw2)
c     &           + mzsq/9._dp/sqrt(cw2*sw2)
c            wphi2 = djitieD(mwsq,ep)*2._dp*mwsq/sqrt(cw2*sw2)
c     &           + cihfenF(0._dp,mwsq,mwsq)*(4._dp*cw2
c     &           + 4._dp/3._dp)*mwsq/sqrt(cw2*sw2)

            slf = alpha/4._dp/pi*(fermi + wphi)
c            ct  = 0._dp !for now
         else if(flag == 'ZZ') then
            fermi = cplx2(0._dp, 0._dp)
            do j = 1, 9
               fermi = fermi + 4._dp/3._dp*(
     &              + (vf(j)**2 + af(j)**2)*(
     &              + psq*djitieD(fmsq(j),ep)
     &              + (psq + 2._dp*fmsq(j))*cihfenF(psq,fmsq(j),fmsq(j))
     &              - psq/3._dp)
     &              - 6._dp*af(j)**2*fmsq(j)*(djitieD(fmsq(j),ep)
     &              + cihfenF(psq,fmsq(j),fmsq(j)))
     &              )
            end do
c            fermi = cplx2(0._dp, 0._dp)
            do j = 1, 3
               fermi = fermi + 8._dp/3._dp*af(j)**2*psq*(
     &              djitieD(fmsq(j),ep) + 5._dp/3._dp
     &              - log(cplx2(-psq/fmsq(j), -1.e-12_dp)))
c               fermi = fermi + 8._dp/3._dp*af(j)**2*psq*(
c     &              log(musq/psq) + 5._dp/3._dp)
            end do
            slf = fermi + djitieD(mwsq,ep)*(
     &           + (3._dp - 19._dp/6._dp/sw2 + 1._dp/6._dp/cw2)*psq
     &           + (4._dp + 1._dp/cw2 - 1._dp/sw2)*mzsq)
     &           + 1._dp/12._dp/cw2/sw2*(
     &           + cihfenF(psq,mwsq,mwsq)*(
     &           - cw2**2*(40._dp*psq + 80._dp*mwsq)
     &           + (cw2 - sw2)**2*(8._dp*mwsq + psq) + 12._dp*mwsq)
     &           + cihfenF(psq,mhsq,mzsq)*(
     &           + 10._dp*mzsq - 2._dp*mhsq + psq
     &           + (mhsq - mzsq)**2/psq)
     &           - 2._dp*mhsq*log(mhsq/mwsq)
     &           - 2._dp*mzsq*log(mzsq/mwsq)
     &           + (10._dp*mzsq - 2._dp*mhsq + psq)*(1._dp
     &           - (mhsq + mzsq)/(mhsq - mzsq)*half*log(mhsq/mzsq)
     &           - half*log(mhsq*mzsq/mwsq**2))
     &           + 2._dp/3._dp*psq*(1._dp + (cw2 - sw2)**2 - 4._dp*cw2**2)
     &           )
c            slf = fermi
            slf = alpha/4._dp/pi*slf
c            print*, 'top: ', cihfenF(mzsq,fmsq(9),fmsq(9))
c            ct = 0._dp
         else if(flag == 'WW') then
            do j = 1, 3
               fermi = fermi + 1._dp/3._dp*(
     &              + (psq - 3._dp/2._dp*fmsq(j))*djitieD(fmsq(j),ep)
     &              + cihfenF(psq,0._dp,fmsq(j))*(psq - fmsq(j)/2._dp
     &              - fmsq(j)**2/2._dp/psq) +2._dp/3._dp*psq
     &              - fmsq(j)/2._dp)
            end do
c            write(6,*) 'sigww lepton only: ', alpha/4._dp/pi/sw2*fermi
c            fermi = cplx2(zip,zip)
c -- j is down-type quark and k is its up-type isospin partner
c            if(1==1) then
            do j = 4, 8, 2
               k = j + 1
               fermi = fermi + (
     &              + djitieD(fmsq(k),ep)/2._dp*(psq - 5._dp/2._dp*fmsq(k)
     &              - fmsq(j)/2._dp) + djitieD(fmsq(j),ep)/2._dp*(psq
     &              - 5._dp/2._dp*fmsq(j) - fmsq(k)/2._dp)
     &              + cihfenF(psq,fmsq(k),fmsq(j))*(psq - (fmsq(k)
     &              + fmsq(j))/2._dp - (fmsq(k) - fmsq(j))**2/2._dp/psq)
     &              - psq/3._dp)
               if(fmsq(j) /= fmsq(k)) then
                  fermi = fermi + (psq - (fmsq(k)
     &                 + fmsq(j))/2._dp)*(1._dp - (fmsq(k)
     &                 + fmsq(j))/(fmsq(k)
     &                 - fmsq(j))*log(fmsq(k)/fmsq(j))/2._dp)
               end if
            end do
c            write(6,*) 'sigww quark only: ', alpha/4._dp/pi/sw2*fermi
c            fermi = cplx2(zip,zip)
c         end if
            slf = fermi - djitieD(mwsq,ep)/3._dp*(19._dp/2._dp*psq
     &           + 3._dp*mwsq*(1._dp - sw2/cw2))
     &           + cihfenF(psq,mzsq,mwsq)*(sw2**2*mzsq - cw2/3._dp*(
     &           + 7._dp*mzsq + 7._dp*mwsq + 10._dp*psq
     &           - 2._dp*(mzsq - mwsq)**2/psq) - 1._dp/6._dp*(mwsq
     &           + mzsq - psq/2._dp - (mzsq - mwsq)**2/2._dp/psq))
     &           + cihfenF(psq,0._dp,mwsq)*sw2/3._dp*(-4._dp*mwsq
     &           - 10._dp*psq + 2._dp*mwsq**2/psq)
     &           + cihfenF(psq,mhsq,mwsq)/6._dp*(5._dp*mwsq - mhsq
     &           + psq/2._dp + (mhsq - mwsq)**2/2._dp/psq)
     &           + mzsq/(mzsq - mwsq)*log(mzsq/mwsq)*(
     &           + cw2/3._dp*(7._dp*mzsq + 7._dp*mwsq + 10._dp*psq
     &           - 4._dp*(mzsq - mwsq)) - sw2**2*mzsq + (2._dp*mwsq
     &           - psq/2._dp)/6._dp)
     &           - mhsq/(mhsq - mwsq)*log(mhsq/mwsq)*(
     &           + 2._dp/3._dp*mwsq + psq/12._dp) - cw2/3._dp*(7._dp*mzsq
     &           + 7._dp*mwsq + 32._dp/3._dp*psq) + sw2**2*mzsq
     &           + (5._dp/3._dp*psq + 4._dp*mwsq - mzsq - mhsq)/6._dp
     &           - sw2/3._dp*(4._dp*mwsq + 32._dp/3._dp*psq)
c            slf = fermi
            slf = alpha/4._dp/pi/sw2*slf
c            write(6,*) 'sigww boson only:', slf
c            print*, 'F(mwsq,mz,mw)', cihfenF(mwsq,mzsq,mwsq)
c            print*, 'F(mwsq,0,mw)', cihfenF(mwsq,0._dp,mwsq)
c            print*, 'F(mwsq,mh,mw)', cihfenF(mwsq,mhsq,mwsq)


         else
            write(6,*) "invalid flag: ", flag
            stop
         end if
      else
c -- to check the pole cancellation (ep = -1)
         if(flag == 'AA') then
            do j = 1, 9
               fermi = fermi
     &              + 4._dp/3._dp*Qf(j)**2*psq*djitieD(fmsq(j),ep)
            end do
            slf = fermi - 3._dp*psq*djitieD(mwsq,ep)
            slf = alpha/4._dp/pi*slf
         else if(flag == 'AZ') then
            do j = 1, 9
               fermi = fermi
     &              - 4._dp/3._dp*Qf(j)*vf(j)*psq*djitieD(fmsq(j),ep)
            end do
            slf = fermi + djitieD(mwsq,ep)/sqrt(cw2*sw2)*(
     &           + (3._dp*cw2 + 1._dp/6._dp)*psq + 2._dp*mwsq)
            slf = alpha/4._dp/pi*slf
         else if(flag == 'ZZ') then
            do j = 1, 9
               fermi = fermi +  4._dp/3._dp*(
     &              + (vf(j)**2 + af(j)**2)*psq*djitieD(fmsq(j),ep)
     &              - 6._dp*af(j)**2*fmsq(j)*djitieD(fmsq(j),ep)
     &              )
            end do
            do j = 1, 3
               fermi = fermi
     &              + 8._dp/3._dp*af(j)**2*psq*djitieD(fmsq(j),ep)
            end do
            slf = fermi + djitieD(mwsq,ep)*(
     &           + (3._dp - 19._dp/6._dp/sw2 + 1._dp/6._dp/cw2)*psq
     &           + (4._dp + 1._dp/cw2 - 1._dp/sw2)*mzsq)
            slf = alpha/4._dp/pi*slf
         else if(flag == 'WW') then
            do j = 1, 3
               fermi = fermi + djitieD(fmsq(j),ep)/3._dp*(psq
     &              - 3._dp/2._dp*fmsq(j))
            end do
            do j = 4, 8, 2
               k = j + 1
               fermi = fermi + 1._dp/2._dp*(
     &              + djitieD(fmsq(k),ep)*(psq - 5._dp/2._dp*fmsq(k)
     &              - fmsq(j)/2._dp) + djitieD(fmsq(j),ep)*(psq
     &              - 5._dp/2._dp*fmsq(j) - fmsq(k)/2._dp)
     &              )
            end do
            slf = fermi - djitieD(mwsq,ep)/3._dp*(
     &           19._dp/2._dp*psq + 3._dp*mwsq*(1._dp - sw2/cw2))
            slf = alpha/4._dp/pi/sw2*slf
         else
            write(6,*) "invalid flag: ", flag
            stop
         end if
      end if


c      print*, 'qcdloop1: ', loopI2(0._dp,fmsq(1),fmsq(1),musq,ep),
c     &     cihfenF(0._dp,fmsq(1),fmsq(1),musq,ep)
c      print*, 'qcdloop2: ', loopI2(0._dp,fmsq(2),fmsq(2),musq,ep),
c     &     cihfenF(0._dp,fmsq(2),fmsq(2),musq,ep)


      end subroutine znen



c -- sigW(psq = 0);
c -- since the expression of sigW in 'znen' is singular when psq = 0,
c -- instead of adding an if-condition, a separate subroutine is given here.
      subroutine sigw0(slf,alpha,sw2,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      real(dp):: alpha,sw2,cw2,djitieD,fmsq(9),mwsq,mzsq,mhsq
      complex(dp):: slf,fermi
      integer ep,j,k
c      common/lepton/mme,mmmu,mmtau
c      common/quark/mup,mdwn,mchrm,mstrg,mbot,mtop

      fmsq(1) = mel**2
      fmsq(2) = mmu**2
      fmsq(3) = mtau**2
      fmsq(4) = md**2
      fmsq(5) = mu**2
      fmsq(6) = ms**2
      fmsq(7) = mc**2
      fmsq(8) = mb**2
      fmsq(9) = mt**2

      mwsq = wmass**2
      mzsq = zmass**2
      mhsq = hmass**2
      cw2 = 1._dp - sw2

      fermi = (0._dp,0._dp)
      if(ep == 0) then
         do j = 1, 3
            fermi = fermi + 1._dp/3._dp*(
     &           - 3._dp/2._dp*fmsq(j)*djitieD(fmsq(j),ep)
     &           - fmsq(j)/2._dp - fmsq(j)/4._dp
     &           )
         end do
c         write(6,*) 'sigww lepton only: ', alpha/4._dp/pi/sw2*fermi
c         fermi = (zip,zip)
         do j = 4, 8, 2
            k = j + 1
            fermi = fermi + (
     &           + djitieD(fmsq(k),ep)/2._dp*(-5._dp/2._dp*fmsq(k)
     &           - fmsq(j)/2._dp) + djitieD(fmsq(j),ep)/2._dp*(
     &           -5._dp/2._dp*fmsq(j) - fmsq(k)/2._dp)
     &           )
            if(fmsq(j) /= fmsq(k)) then
               fermi = fermi - (fmsq(k)
     &              + fmsq(j))/2._dp*(1._dp - (fmsq(k)
     &              + fmsq(j))/(fmsq(k)
     &              - fmsq(j))*log(fmsq(k)/fmsq(j))/2._dp)
     &              - ((fmsq(k) + fmsq(j))/2._dp
     &              - fmsq(k)*fmsq(j)/(fmsq(k)
     &              - fmsq(j))*log(fmsq(k)/fmsq(j)))/2._dp
            end if
         end do
c         write(6,*) 'sigww quark only: ', alpha/4._dp/pi/sw2*fermi
c         fermi = (zip,zip)
         slf = fermi - djitieD(mwsq,ep)*mwsq*(1._dp - sw2/cw2)
     &        + mzsq/(mzsq - mwsq)*log(mzsq/mwsq)*(
     &        + cw2/3._dp*(7._dp*mzsq + 7._dp*mwsq
     &        - 4._dp*(mzsq - mwsq)) - sw2**2*mzsq + mwsq/3._dp)
     &        - mhsq/(mhsq - mwsq)*log(mhsq/mwsq)*2._dp/3._dp*mwsq
     &        - cw2/3._dp*(7._dp*mzsq + 7._dp*mwsq) + sw2**2*mzsq
     &        + (4._dp*mwsq - mzsq - mhsq)/6._dp
     &        - sw2/3._dp*(4._dp*mwsq)
     &        + (2._dp*cw2/3._dp + 1._dp/12._dp)*((mzsq + mwsq)/2._dp
     &        - mzsq*mwsq/(mzsq - mwsq)*log(mzsq/mwsq))
     &        + sw2*mwsq/3._dp + ((mhsq + mwsq)/2._dp - mhsq*mwsq/(mhsq
     &        - mwsq)*log(mhsq/mwsq))/12._dp
c         write(6,*) 'sigww boson only:', alpha/4._dp/pi/sw2*slf
      else
            do j = 1, 3
               fermi = fermi - djitieD(fmsq(j),ep)/2._dp*fmsq(j)
            end do
            do j = 4, 8, 2
               k = j + 1
               fermi = fermi + 1._dp/2._dp*(
     &              + djitieD(fmsq(k),ep)*(- 5._dp/2._dp*fmsq(k)
     &              - fmsq(j)/2._dp) + djitieD(fmsq(j),ep)*(
     &              - 5._dp/2._dp*fmsq(j) - fmsq(k)/2._dp)
     &              )
            end do
            slf = fermi - djitieD(mwsq,ep)*mwsq*(1._dp - sw2/cw2)
         end if

      slf = alpha/4._dp/pi/sw2*slf

      end subroutine sigw0


c -- derivative of gamma self energy
      subroutine dznenaa(dz,alpha,psq,ep)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'masses.f'
      include 'constants.f'
      include 'ewcharge.f'
      real(dp):: psq,djitieD,mwsq,alpha,Qf(9),fmsq(9)
      complex(dp) :: ccdb0
      complex(dp):: dz,cihfenF,dfermi,dslf
      integer ep,j
      external djitieD,cihfenF
c      common/em/aemmz
c      common/lepton/mme,mmmu,mmtau
c      common/quark/mup,mdwn,mchrm,mstrg,mbot,mtop

      mwsq = wmass**2
      Qf(1:3) = -1._dp
      Qf(4:8) = Q(1:5)*sqrt(xn)
      Qf(9) = Q(2)*sqrt(xn)

      fmsq(1) = mel**2
      fmsq(2) = mmu**2
      fmsq(3) = mtau**2
      fmsq(4) = md**2
      fmsq(5) = mu**2
      fmsq(6) = ms**2
      fmsq(7) = mc**2
      fmsq(8) = mb**2
      fmsq(9) = mt**2

      dfermi = (0._dp,0._dp)
      if(ep == 0) then
         do j = 1, 9
            dfermi = dfermi + 4._dp/3._dp*Qf(j)**2*(
     &           + djitieD(fmsq(j),ep)
     &           + cihfenF(psq,fmsq(j),fmsq(j))
     &           + (psq + 2._dp*fmsq(j))*ccdb0(psq,fmsq(j),fmsq(j))
     &           - 1._dp/3._dp)
         end do
         dslf = dfermi - 3._dp*djitieD(mwsq,ep)
     &        - 3._dp*cihfenF(psq,mwsq,mwsq)
     &        - (3._dp*psq + 4._dp*mwsq)*ccdb0(psq,mwsq,mwsq)
      else if(ep == -1) then
         do j = 1, 9
            dfermi = dfermi + 4._dp/3._dp*Qf(j)**2*(
     &           + djitieD(fmsq(j),ep))
         end do
         dslf = dfermi - 3._dp*djitieD(mwsq,ep)
      else
         write(6,*) 'only single pole exists, try ep = 0/-1'
         return
      end if

      dz = alpha/4._dp/pi*dslf

      return

      end subroutine dznenaa




c -- derivative of self-energy ZZ
      subroutine dznenzz(dz,alpha,sw2,psq,ep)
      include 'types.f'
      include 'nf.f'
      include 'masses.f'
      include 'constants.f'
      include 'zcouple.f'
      include 'cplx.h'
      real(dp):: psq,mzsq,mwsq,mhsq,sw2,cw2,djitieD,alpha,
     &     gfp(9),gfm(9),fmsq(9),vf(9),af(9)
      complex(dp):: dz,dslf,dfermi,cihfenF
      complex(dp) :: ccdb0
      integer ep,j
      external djitieD,cihfenF


      mzsq = zmass**2
      mwsq = wmass**2
      mhsq = hmass**2

c      sw2 = xw
      cw2 = 1._dp - sw2

      gfp(1:3) = re
      gfm(1:3) = le
      gfp(4:8) = r(1:5)*sqrt(xn)
      gfm(4:8) = l(1:5)*sqrt(xn)
      gfp(9) = r(2)*sqrt(xn)
      gfm(9) = l(2)*sqrt(xn)

      vf = (gfm + gfp)/2._dp
      af = (gfm - gfp)/2._dp


      fmsq(1) = mel**2
      fmsq(2) = mmu**2
      fmsq(3) = mtau**2
      fmsq(4) = md**2
      fmsq(5) = mu**2
      fmsq(6) = ms**2
      fmsq(7) = mc**2
      fmsq(8) = mb**2
      fmsq(9) = mt**2


      dfermi = cplx2(0._dp,0._dp)
      if(ep == 0) then
         do j = 1, 9
            dfermi = dfermi + 4._dp/3._dp*(
     &           + (vf(j)**2 + af(j)**2)*(
     &           + djitieD(fmsq(j),ep)
     &           + (psq + 2._dp*fmsq(j))*ccdb0(psq,fmsq(j),fmsq(j))
     &           + cihfenF(psq,fmsq(j),fmsq(j))
     &           - 1._dp/3._dp)
     &           - 6._dp*af(j)**2*fmsq(j)*ccdb0(psq,fmsq(j),fmsq(j))
     &           )
         end do

c         print*, 'dfermi --- 1: ', dfermi

         do j = 1, 3
            dfermi = dfermi + 8._dp/3._dp*af(j)**2*(
     &           djitieD(fmsq(j),ep) + 5._dp/3._dp
     &           - log(cplx2(-psq/fmsq(j), -1.e-12_dp))
     &           - 1._dp
     &           )
         end do

c         print*, 'dfermi --- 2: ', dfermi

         dslf = dfermi + djitieD(mwsq,ep)*(
     &        + (3._dp - 19._dp/6._dp/sw2 + 1._dp/6._dp/cw2))
     &        + 1._dp/12._dp/cw2/sw2*(
     &        + ccdb0(psq,mwsq,mwsq)*(
     &        - cw2**2*(40._dp*psq + 80._dp*mwsq)
     &        + (cw2 - sw2)**2*(8._dp*mwsq + psq) + 12._dp*mwsq)
     &        + cihfenF(psq,mwsq,mwsq)*(
     &        - cw2**2*40._dp + (cw2 - sw2)**2)
     &        + ccdb0(psq,mhsq,mzsq)*(
     &        + 10._dp*mzsq - 2._dp*mhsq + psq
     &        + (mhsq - mzsq)**2/psq)
     &        + cihfenF(psq,mhsq,mzsq)*(
     &        + 1._dp - (mhsq - mzsq)**2/psq**2)
     &        + (1._dp
     &        - (mhsq + mzsq)/(mhsq - mzsq)*half*log(mhsq/mzsq)
     &        - half*log(mhsq*mzsq/mwsq**2))
     &        + 2._dp/3._dp*(1._dp + (cw2 - sw2)**2 - 4._dp*cw2**2)
     &        )

      else if(ep == -1) then
         do j = 1, 9
            dfermi = dfermi + 4._dp/3._dp*(
     &           + (vf(j)**2 + af(j)**2)*(
     &           + djitieD(fmsq(j),ep))
     &           )
         end do

         do j = 1, 3
            dfermi = dfermi + 8._dp/3._dp*af(j)**2*(
     &           djitieD(fmsq(j),ep)
     &           )
         end do

         dslf = dfermi + djitieD(mwsq,ep)*(
     &        + (3._dp - 19._dp/6._dp/sw2 + 1._dp/6._dp/cw2))

      else
         write(6,*) 'only single pole exists, try ep = 0/-1'
         return
      end if

c      print*, dslf
      dz = alpha/4._dp/pi*dslf

      return

      end subroutine dznenzz





c -- delta(alpha)(Mz)
      subroutine damz(da,alpha,ep)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'masses.f'
      include 'constants.f'
      include 'ewcharge.f'
      real(dp):: da,djitieD,mzsq,mwsq,alpha,Qf(8),fmsq(8)
      complex(dp):: cihfenF,fermi,slf
      integer ep,j
      external djitieD,cihfenF
c      common/em/aemmz
c      common/lepton/mme,mmmu,mmtau
c      common/quark/mup,mdwn,mchrm,mstrg,mbot,mtop

      mzsq = zmass**2
      mwsq = wmass**2
      Qf(1:3) = -1._dp
      Qf(4:8) = Q(1:5)*sqrt(xn)

      fmsq(1) = mel**2
      fmsq(2) = mmu**2
      fmsq(3) = mtau**2
      fmsq(4) = md**2
      fmsq(5) = mu**2
      fmsq(6) = ms**2
      fmsq(7) = mc**2
      fmsq(8) = mb**2

      fermi = (0._dp,0._dp)
      if(ep == 0) then
         do j = 1, 8
            fermi = fermi + 4._dp/3._dp*Qf(j)**2*(
     &           + djitieD(fmsq(j),ep)*mzsq
     &           + (mzsq + 2._dp*fmsq(j))*cihfenF(mzsq,fmsq(j),fmsq(j))
     &           - mzsq/3._dp)
         end do
         slf = fermi - 3._dp*djitieD(mwsq,ep)*mzsq
     &        - (3._dp*mzsq + 4._dp*mwsq)*cihfenF(mzsq,mwsq,mwsq)
      else
         do j = 1, 8
            fermi = fermi + 4._dp/3._dp*Qf(j)**2*(
     &           + djitieD(fmsq(j),ep))*mzsq
         end do
         slf = fermi - 3._dp*djitieD(mwsq,ep)*mzsq
      end if

      da = -alpha/4._dp/pi*real(slf,dp)/mzsq

      end subroutine damz


c -- self energy counter term
      subroutine znenteghaon(tgh,alpha,sw2,psq,ep,flag)
      implicit none
      include 'types.f'
      include 'masses.f'
      character(len=2):: flag
      real(dp):: tgh,psq,sw2,cw2,mzsq,mwsq,dzg(2),dzz(2),dzgz(2),
     &     dww(2),rdmz,rdmw,alpha
      complex(dp):: sigz0,dsig0,dmz,dmw
      integer ep

      mzsq = zmass**2
      mwsq = wmass**2

c      sw2 = xw
      cw2 = 1._dp - sw2

      call znen(sigz0,alpha,sw2,0._dp,ep,'AZ')
      call znen(dmz,alpha,sw2,mzsq,ep,'ZZ')
      call znen(dmw,alpha,sw2,mwsq,ep,'WW')
      call dznenaa(dsig0,alpha,0._dp,ep)
c      print*, 'dsig0: ', dsig0, ep
c      print*, 'dmz', dmz
c      print*, 'dmw', dmw

      rdmz = real(dmz,dp)
      rdmw = real(dmw,dp)

      dzg(2) = -real(dsig0,dp)
      dzg(1) = dzg(2) - sqrt(sw2/cw2)*real(sigz0,dp)/mzsq
      dzz(2) = dzg(2)
     &     - 2._dp*(cw2 - sw2)/sqrt(sw2*cw2)*real(sigz0,dp)/mzsq
     &     + (cw2 - sw2)/sw2*(rdmz/mzsq - rdmw/mwsq)
      dzz(1) = dzg(2)
     &     - (3._dp*cw2 - 2._dp*sw2)/sqrt(sw2*cw2)*real(sigz0,dp)/mzsq
     &     + (cw2 - sw2)/sw2*(rdmz/mzsq - rdmw/mwsq)
      dzgz = sqrt(cw2*sw2)/(cw2 - sw2)*(dzz - dzg)
      dww = (sw2*dzg - cw2*dzz)/(sw2 - cw2)

c      print*, 'psq: ', psq
c      print*, 'tepi: ', rdmz/mzsq - rdmw/mwsq,
c     &     sqrt(sw2/cw2)*(3._dp*dzgz(2) - 2._dp*dzgz(1))

      if(flag == 'AA') then
         tgh = psq*dzg(2)
c         print*, 'tgh: ', tgh
      else if(flag == 'AZ') then
         tgh = -psq*dzgz(2) +(dzgz(1) - dzgz(2))*mzsq
      else if(flag == 'ZZ') then
         tgh = dzz(2)*(psq - mzsq) - rdmz
      else if(flag == 'WW') then
         tgh = dww(2)*(psq - mwsq) - rdmw
      else
         write(6,*) 'invalid flag: ', flag
         stop
      end if

      end subroutine znenteghaon



c -- integral F in (B.6)
      function cihfenF(psq,m1sq,m2sq)
      use loopI2_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'cplx.h'
      real(dp):: psq,m1sq,m2sq
      complex(dp):: mu2,som,cihfenF

      mu2 = cplx2(musq/sqrt(m1sq*m2sq),0._dp)
      if(psq == 0._dp) then
         cihfenF = cplx2(0._dp,0._dp)
      else if((m1sq /= 0._dp) .and. (m2sq /= 0._dp)) then
         if(m1sq == m2sq) then
            cihfenF =
     &           +  loopI2(psq,m1sq,m2sq,musq,0)
     &           - log(mu2)
c Series expansion for small m1sq/psq
            if (m1sq/psq < 1.e-8_dp) then
c              write(6,*) 'm1sq/psq',m1sq/psq,cihfenF/(2._dp-log(psq/m1sq)+im*pi+2._dp*m1sq/psq*(1._dp-im*pi+log(psq/m1sq)))
              cihfenF=2._dp-log(psq/m1sq)+im*pi+2._dp*m1sq/psq*(1._dp-im*pi+log(psq/m1sq))
            endif
         else
            cihfenF = cplx1(-1._dp
     &           + (m1sq + m2sq)/(m1sq - m2sq)*log(m1sq/m2sq)/2._dp)
            cihfenF = cihfenF + loopI2(psq,m1sq,m2sq,musq,0)
     &           - log(mu2)
         end if
      else if((m1sq == 0._dp) .and. (m2sq /= 0._dp)) then
         som = cplx2(1._dp - psq/m2sq, -1.e-12_dp)
         cihfenF = cone + (m2sq/psq - cone)*log(som)
      else if((m1sq /= 0._dp) .and. (m2sq == 0._dp)) then
         som = cplx2(1._dp - psq/m1sq, -1.e-12_dp)
         cihfenF = cone + (m1sq/psq - cone)*log(som)
      else
         cihfenF = cplx2(0._dp,0._dp)
      end if

      end function cihfenF



c -- Delta in (B.1)
      function djitieD(msq,ep)
      implicit none
      include 'types.f'
      include 'scale.f'
      real(dp):: msq,djitieD
      integer ep


      if(ep == -1) then
         djitieD = 1._dp
      else if(ep == 0) then
         djitieD = log(musq/msq)
      else
         djitieD = 0._dp
      end if

      end function djitieD





c      block data fermionmass
c      implicit none
c      include 'types.f'
c      real(dp):: mme,mmmu,mmtau,mup,mdwn,mchrm,mstrg,mbot,mtop
c      common/lepton/mme,mmmu,mmtau
c      common/quark/mup,mdwn,mchrm,mstrg,mbot,mtop
c      data mme,mmmu,mmtau/+0.51099892e-3_dp,+105.658369e-3_dp,+1.777_dp/
c      data mup,mdwn,mchrm,mstrg,mbot,mtop/+0.066_dp,+0.066_dp,+1.2_dp,
c     &     +0.15_dp,+4.6_dp,+173.2_dp/
c      end block data fermionmass
