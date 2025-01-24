!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dmmonojn_scal(j1,j2,j5,p,n1,qqbg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
c----- ME squared for q(j1)+qb(j2)+x(j3)+x~(j4)+g(j5) X n_nu
      integer:: j1,j2,j5,j3,j4
      real(dp):: n1(4),p(mxpart,4),q(mxpart,4)
      complex(dp)::amp(2,2,2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: qqbg(2)
c      real(dp):: z1jetn
c----- return as funciton of quark line helicity
      integer:: h1,h2,h3
      complex(dp):: amp_p(2),amp_dec(2,2)
      real(dp):: bp,beta,s34
      real(dp):: fac
      complex(dp):: s(mxpart,mxpart)

      fac=0.25d0

      j3=3
      j4=4

      call checkndotp(p,n1,j5)

      if(xmass>1d-8) then
c---------generate massless phase space
         call gen_masslessvecs(p,q,j3,j4)
c---------generate spinors
         call spinoru(5,q,za,zb)
         call spinork(5,q,zab,zba,n1)
      else
c--------massless dm can use usual spinoru
         call spinoru(5,p,za,zb)
         call spinork(5,p,zab,zba,n1)

      endif


      do h1=1,5
         do h2=1,5
            s(h1,h2)=za(h1,h2)*zb(h2,h1)
         enddo
      enddo

cbp = 1/2(1+beta)
c beta = sqrt(1-4xmass**2/s34)
      s34=Dble(za(j3,j4)*zb(j4,j3))
      beta=sqrt(1d0-4d0*xmass**2/s34)
      bp=0.5d0*(one+beta)

c      write(6,*) 'In gvec ',j3,j4,bp,s34,xmass
      call dm_scal_decay(j3,j4,za,zb,bp,amp_dec)
c------ bulid amplitudes
c------- helicity conserving

      amp_p(1)=(za(j1,j2)*zab(j1,j1))/s(j1,j5) -
     &  (za(j2,j5)*zab(j1,j5))/s(j1,j5) -
     &  (za(j1,j2)*zab(j2,j2))/s(j2,j5)
     &     - (za(j1,j5)*zab(j2,j5))/s(j2,j5)
      amp_p(2)=-((zb(j2,j1)*zba(j1,j1))/s(j1,j5)) +
     &  (zb(j5,j2)*zba(j1,j5))/s(j1,j5) +
     &  (zb(j2,j1)*zba(j2,j2))/s(j2,j5)
     &     + (zb(j5,j1)*zba(j2,j5))/s(j2,j5)



      do h1=1,2
         do h2=1,2
            do h3=1,2
               amp(h1,h2,h3)=amp_p(h1)*amp_dec(h2,h3)
            enddo
         enddo
      enddo


       qqbg(:)=0d0
       do h1=1,2
          do h2=1,2
             do h3=1,2
                 qqbg(h1)=qqbg(h1)+fac*abs(amp(h1,h2,h3))**2
              enddo
           enddo
        enddo

c        write(6,*) j1,j2,j5
c        write(6,*) one/16d0*abs(amp(2,2,1))**2/s(j3,j4)**2
c     &/z1jetn(j1,j2,j5,p,n1)
c        write(6,*) one/16d0*abs(amp(2,1,2))**2/s(j3,j4)**2
c     &/z1jetn(j2,j1,j5,p,n1)
c        pause

      return
      end
