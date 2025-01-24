!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function getet(E,px,py,pz)
      implicit none
      include 'types.f'
      real(dp):: getet
c--- given (E,px,py,pz) for a four-vector, calculates the corresponding
c--- Et or Pt, depending on the parameter that is set in mdata.f

      real(dp):: E,px,py,pz,etsq
      include 'useet.f'

      if (useEt) then
c--- this is the formula for Et
        etsq=px**2+py**2
        getet=sqrt(etsq)*E/sqrt(etsq+pz**2)
      else
c--- this is the formula for pt
        getet=sqrt(px**2+py**2)
      endif

      return
      end

      function pt(j,p)
      implicit none
      include 'types.f'
      real(dp):: pt
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt

      include 'mxpart.f'
      integer:: j
      real(dp):: p(mxpart,4),getet

      pt=getet(p(j,4),p(j,1),p(j,2),p(j,3))

      return
      end

      function ptpure(p)
        use types
        implicit none

        real(dp), intent(in) :: p(4)
        real(dp) :: ptpure, getet

        ptpure = getet(p(4),p(1),p(2),p(3))

      end function


      function pttwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: pttwo
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt

      include 'mxpart.f'
      integer:: j,k
      real(dp):: p(mxpart,4),getet

      pttwo=getet(p(j,4)+p(k,4),p(j,1)+p(k,1),
     &            p(j,2)+p(k,2),p(j,3)+p(k,3))

      return
      end

      function ptthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: ptthree
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt

      include 'mxpart.f'
      integer:: j,k,m
      real(dp):: p(mxpart,4),getet

      ptthree=getet(p(j,4)+p(k,4)+p(m,4),p(j,1)+p(k,1)+p(m,1),
     &              p(j,2)+p(k,2)+p(m,2),p(j,3)+p(k,3)+p(m,3))

      return
      end

      function ptfour(j,k,m,n,p)
      implicit none
      include 'types.f'
      real(dp):: ptfour
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt

      include 'mxpart.f'
      integer:: j,k,m,n
      real(dp):: p(mxpart,4),getet

      ptfour=getet(p(j,4)+p(k,4)+p(m,4)+p(n,4),
     &             p(j,1)+p(k,1)+p(m,1)+p(n,1),
     &             p(j,2)+p(k,2)+p(m,2)+p(n,2),
     &             p(j,3)+p(k,3)+p(m,3)+p(n,3))

      return
      end

      function ptsix(j1,k1,m1,j2,k2,m2,p)
      implicit none
      include 'types.f'
      real(dp):: ptsix
c--- This routine is now just a wrapper to getet, which
c--- computes either pt or et, depending on the value of useEt

      include 'mxpart.f'
      integer:: j1,k1,m1,j2,k2,m2
      real(dp):: p(mxpart,4),getet

      ptsix=getet(p(j1,4)+p(k1,4)+p(m1,4)+p(j2,4)+p(k2,4)+p(m2,4),
     &            p(j1,1)+p(k1,1)+p(m1,1)+p(j2,1)+p(k2,1)+p(m2,1),
     &            p(j1,2)+p(k1,2)+p(m1,2)+p(j2,2)+p(k2,2)+p(m2,2),
     &            p(j1,3)+p(k1,3)+p(m1,3)+p(j2,3)+p(k2,3)+p(m2,3))

      return
      end

      pure function etarappure(p)
        use types
        implicit none

        real(dp) :: etarappure
        real(dp), intent(in) :: p(4)

        etarappure = sqrt(p(1)**2+p(2)**2+p(3)**2)
        etarappure = (etarappure+p(3))/(etarappure-p(3))

        if (etarappure < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
          etarappure = 100._dp
        else
          etarappure = 0.5_dp*log(etarappure)
        endif
      end function

c---returns the value of the pseudorapidity
      function etarap(j,p)
        use types
        implicit none

        real(dp) :: etarap
        include 'mxpart.f'

        integer, intent(in) :: j
        real(dp), intent(in) :: p(mxpart,4)

        real(dp) :: etarappure

        etarap = etarappure(p(j,:))
      end

      pure function aetarap(j,p)
          use types
      implicit none
c---returns the absolute value of the pseudorapidity
      include 'mxpart.f'

      real(dp):: aetarap
      integer, intent(in) :: j
      real(dp), intent(in) :: p(mxpart,4)
      real(dp), parameter :: tiny = 1.e-9_dp

      aetarap=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
      aetarap=p(j,3)/aetarap
      if  ((1d0+aetarap < tiny) .or. (1d0-aetarap < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
          aetarap=100d0
      else
          aetarap=0.5_dp*abs(log((1._dp+aetarap)/(1._dp-aetarap)))
      endif

      return
      end

      pure function yrappure(p)
      implicit none
      include 'types.f'
      real(dp):: yrappure
c---returns the value of the rapidity
      real(dp), intent(in):: p(4)

      yrappure=(p(4)+p(3))/(p(4)-p(3))
      if ((p(4) < 1.e-13_dp) .or. (yrappure < 1.e-13_dp)) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrappure=100._dp
      else
      yrappure=0.5_dp*log(yrappure)
      endif
      return
      end

      function yrap(j,p)
      implicit none
      include 'types.f'
      real(dp):: yrap

c---returns the value of the rapidity
      include 'mxpart.f'
      integer:: j
      real(dp):: p(mxpart,4)
      yrap=(p(j,4)+p(j,3))/(p(j,4)-p(j,3))
      if ((p(j,4) < 1.e-13_dp) .or. (yrap < 1.e-13_dp)) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrap=100._dp
      else
      yrap=0.5_dp*log(yrap)
      endif
      return
      end

      function ayrap(j,p)
c---returns the absolute value of the rapidity
      implicit none
      include 'types.f'
      real(dp):: ayrap
      include 'mxpart.f'
      integer:: j
      real(dp):: p(mxpart,4)
      real(dp), parameter :: tiny = 1.e-9_dp

      ayrap=p(j,3)/p(j,4)
      if  ((1d0+ayrap < tiny) .or. (1d0-ayrap < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
          ayrap=100d0
      else
          ayrap=0.5_dp*abs(log((1._dp+ayrap)/(1._dp-ayrap)))
      endif

      return
      end

c--- this is the rapidity of pair j,k
      function yraptwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: yraptwo

      include 'mxpart.f'
      integer:: j,k
      real(dp):: p(mxpart,4)
      yraptwo=(p(j,4)+p(k,4)+p(j,3)+p(k,3))
     &       /(p(j,4)+p(k,4)-p(j,3)-p(k,3))
      if (yraptwo < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yraptwo=100._dp
      else
      yraptwo=0.5_dp*log(yraptwo)
      endif

      return
      end

c--- this is the pseudo-rapidity of pair j,k
      function etaraptwo(j,k,p)
      implicit none
      include 'types.f'
      real(dp):: etaraptwo

      include 'mxpart.f'
      integer:: j,k
      real(dp):: p(mxpart,4)

      etaraptwo=sqrt((p(j,1)+p(k,1))**2+(p(j,2)+p(k,2))**2
     &               +(p(j,3)+p(k,3))**2)
      if (abs(etaraptwo)-abs(p(j,3)+p(k,3)) < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etaraptwo=100._dp
      else
      etaraptwo=(etaraptwo+p(j,3)+p(k,3))
     &         /(etaraptwo-p(j,3)-p(k,3))
      etaraptwo=0.5_dp*log(etaraptwo)
      endif

      return
      end

      function yrapthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: yrapthree
c--- this is the rapidity of the combination j+k+m

      include 'mxpart.f'
      integer:: j,k,m
      real(dp):: p(mxpart,4)
      yrapthree=(p(j,4)+p(k,4)+p(m,4)+p(j,3)+p(k,3)+p(m,3))
     &         /(p(j,4)+p(k,4)+p(m,4)-p(j,3)-p(k,3)-p(m,3))
      if (yrapthree < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapthree=100._dp
      else
      yrapthree=0.5_dp*log(yrapthree)
      endif

      return
      end

      function etarapthree(j,k,m,p)
      implicit none
      include 'types.f'
      real(dp):: etarapthree
c--- this is the pseudo-rapidity of the combination j+k+m

      include 'mxpart.f'
      integer:: j,k,m
      real(dp):: p(mxpart,4)

      etarapthree=
     &    sqrt((p(j,1)+p(k,1)+p(m,1))**2+(p(j,2)+p(k,2)+p(m,2))**2
     &         +(p(j,3)+p(k,3)+p(m,3))**2)
      if (abs(etarapthree)-abs(p(j,3)+p(k,3)+p(m,3)) < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      etarapthree=100._dp
      else
      etarapthree=(etarapthree+p(j,3)+p(k,3)+p(m,3))
     &           /(etarapthree-p(j,3)-p(k,3)-p(m,3))
      etarapthree=0.5_dp*log(etarapthree)
      endif

      return
      end

      function yrapfour(j,k,m,n,p)
      implicit none
      include 'types.f'
      real(dp):: yrapfour
c--- this is the rapidity of the combination j+k+m+n

      include 'mxpart.f'
      integer:: j,k,m,n
      real(dp):: p(mxpart,4)
      yrapfour=(p(j,4)+p(k,4)+p(m,4)+p(n,4)+p(j,3)+p(k,3)+p(m,3)+p(n,3))
     &        /(p(j,4)+p(k,4)+p(m,4)+p(n,4)-p(j,3)-p(k,3)-p(m,3)-p(n,3))
      if (yrapfour < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapfour=100._dp
      else
      yrapfour=0.5_dp*log(yrapfour)
      endif

      return
      end

      function yrapsix(j1,k1,m1,j2,k2,m2,p)
      implicit none
      include 'types.f'
      real(dp):: yrapsix
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2

      include 'mxpart.f'
      integer:: j1,k1,m1,j2,k2,m2
      real(dp):: p(mxpart,4)
      yrapsix=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3))
     &       /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)
     &        +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3))
      if (yrapsix < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapsix=100._dp
      else
      yrapsix=0.5_dp*log(yrapsix)
      endif

      return
      end

      function yrapseven(j1,k1,m1,j2,k2,m2,n,p)
      implicit none
      include 'types.f'
      real(dp):: yrapseven
c--- this is the rapidity of the combination j1+k1+m1+j2+k2+m2+n

      include 'mxpart.f'
      integer:: j1,k1,m1,j2,k2,m2,n
      real(dp):: p(mxpart,4)
      yrapseven=(p(j1,4)+p(k1,4)+p(m1,4)+p(j1,3)+p(k1,3)+p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)+p(j2,3)+p(k2,3)+p(m2,3)+p(n,3))
     &         /(p(j1,4)+p(k1,4)+p(m1,4)-p(j1,3)-p(k1,3)-p(m1,3)+p(n,4)
     &          +p(j2,4)+p(k2,4)+p(m2,4)-p(j2,3)-p(k2,3)-p(m2,3)-p(n,3))
      if (yrapseven < 1.e-13_dp) then
c-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      yrapseven=100._dp
      else
      yrapseven=0.5_dp*log(yrapseven)
      endif
       return
       end
c -- RR: mass definitions

      function puremass(p)
          use types
          implicit none

          real(dp), intent(in) :: p(4)
          real(dp) :: puremass

       puremass=sqrt( p(4)**2-p(1)**2-p(2)**2-p(3)**2 )

      end function

       function onemass(j,p)
       implicit none
      include 'types.f'
      real(dp):: onemass

       include 'mxpart.f'
       integer:: j
       real(dp):: p(mxpart,4)
       onemass=p(j,4)**2-p(j,1)**2-p(j,2)**2-p(j,3)**2
       onemass=sqrt(onemass)
       return
       end

      function twomass(j,k1,p)
      implicit none
      include 'types.f'
      real(dp):: twomass

      include 'mxpart.f'
      integer:: j,k1
      real(dp):: p(mxpart,4),pjk(4)
      pjk(:)=p(j,:)+p(k1,:)
      twomass=pjk(4)**2-pjk(1)**2-pjk(2)**2-pjk(3)**2
      twomass=sqrt(twomass)
      return
      end

      function threemass(j,k1,k2,p)
      implicit none
      include 'types.f'
      real(dp):: threemass

      include 'mxpart.f'
      integer:: j,k1,k2
      real(dp):: p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)
      threemass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      threemass=sqrt(threemass)
      return
      end

      function fourmass(j,k1,k2,k3,p)
      implicit none
      include 'types.f'
      real(dp):: fourmass

      include 'mxpart.f'
      integer:: j,k1,k2,k3
      real(dp):: p(mxpart,4),ptot(4)
      ptot(:)=p(j,:)+p(k1,:)+p(k2,:)+p(k3,:)
      fourmass=ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2
      fourmass=sqrt(fourmass)
      return
      end






