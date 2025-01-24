!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine zzggamps(j1,j2,j3,j4,j5,j6,j7,j8,fAB,fBA)
c     Amplitudes for the process
c     0 --> qbar(j1) + q(j2)+e^-(j3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)+g(p8)
c     only s34,s56,propagators included.
c     f(2,2,2,2,2) \equiv f(polg1,polg2,polq,pol34,pol56)
c     only a-type diagrams
c     pol=1=LH, pol=2=RH
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'

      integer::j1,j2,j3,j4,j5,j6,j7,j8,polg1,polg2,p3,p4,p5,p6
      complex(dp):: a8amm(2,2),a8app(2,2),a8amp(2,2),a8apm(2,2)
      complex(dp):: fAB(2,2,2,2,2),fBA(2,2,2,2,2)
      integer,parameter:: minus=1,pplus=2


c---a8treea is a routine setup for WW where boson emission is ordered
c---for Z's we have to include both orders 3,4 <--> 5,6.
c   we will do the other ordering first
c   taking care to swap h34 and h56 i.e. a8amp gives "minus,pplus" etc.
c---Handle all cases with negative helicity quark line
c---Remember that the a8tree's are in DKS notation,
c-  so we swap (1<->2) and (5<-->6) to get to our standard notation
      p3=j5
      p4=j6
      p5=j3
      p6=j4
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,za,zb,a8amm)
      call a8treea(j2,j1,p4,p3,p6,p5,j7,j8,za,zb,a8apm)
      call a8treea(j2,j1,p3,p4,p5,p6,j7,j8,za,zb,a8amp)
      call a8treea(j2,j1,p4,p3,p5,p6,j7,j8,za,zb,a8app)
      do polg1=1,2
      do polg2=1,2
      fAB(polg1,polg2,minus,minus,minus)=+a8amm(polg1,polg2)
      fAB(polg1,polg2,minus,pplus,pplus)=+a8app(polg1,polg2)
      fAB(polg1,polg2,minus,minus,pplus)=+a8apm(polg1,polg2)
      fAB(polg1,polg2,minus,pplus,minus)=+a8amp(polg1,polg2)
      enddo
      enddo

c---Other color order swap (7<-->8)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,za,zb,a8amm)
      call a8treea(j2,j1,p4,p3,p6,p5,j8,j7,za,zb,a8apm)
      call a8treea(j2,j1,p3,p4,p5,p6,j8,j7,za,zb,a8amp)
      call a8treea(j2,j1,p4,p3,p5,p6,j8,j7,za,zb,a8app)
      do polg1=1,2
      do polg2=1,2
      fBA(polg1,polg2,minus,minus,minus)=+a8amm(polg2,polg1)
      fBA(polg1,polg2,minus,pplus,pplus)=+a8app(polg2,polg1)
      fBA(polg1,polg2,minus,minus,pplus)=+a8apm(polg2,polg1)
      fBA(polg1,polg2,minus,pplus,minus)=+a8amp(polg2,polg1)
      enddo
      enddo

c-- opposite quark helicity by taking complex conjugation (za<->zb)
c-- and swapping all helicities
c      call a8treeb(j2,j1,j4,j3,j5,p6,j7,j8,zb,za,a8b)
c      fAB(2,polg1,polg2,mplus)=+a8b(3-polg1,3-polg2)
      call a8treea(j2,j1,p4,p3,p5,p6,j7,j8,zb,za,a8amm)
      call a8treea(j2,j1,p3,p4,p5,p6,j7,j8,zb,za,a8apm)
      call a8treea(j2,j1,p4,p3,p6,p5,j7,j8,zb,za,a8amp)
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,zb,za,a8app)
      do polg1=1,2
      do polg2=1,2
      fAB(polg1,polg2,pplus,minus,minus)=+a8amm(3-polg1,3-polg2)
      fAB(polg1,polg2,pplus,pplus,pplus)=+a8app(3-polg1,3-polg2)
      fAB(polg1,polg2,pplus,minus,pplus)=+a8apm(3-polg1,3-polg2)
      fAB(polg1,polg2,pplus,pplus,minus)=+a8amp(3-polg1,3-polg2)
      enddo
      enddo

c---Other color order swap (7<-->8)
      call a8treea(j2,j1,p4,p3,p5,p6,j8,j7,zb,za,a8amm)
      call a8treea(j2,j1,p3,p4,p5,p6,j8,j7,zb,za,a8apm)
      call a8treea(j2,j1,p4,p3,p6,p5,j8,j7,zb,za,a8amp)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,zb,za,a8app)
      do polg1=1,2
      do polg2=1,2
      fBA(polg1,polg2,pplus,minus,minus)=+a8amm(3-polg2,3-polg1)
      fBA(polg1,polg2,pplus,pplus,pplus)=+a8app(3-polg2,3-polg1)
      fBA(polg1,polg2,pplus,minus,pplus)=+a8apm(3-polg2,3-polg1)
      fBA(polg1,polg2,pplus,pplus,minus)=+a8amp(3-polg2,3-polg1)
      enddo
      enddo

c     Now we have to add in the normalorder for boson
      p3=j3
      p4=j4
      p5=j5
      p6=j6
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,za,zb,a8amm)
      call a8treea(j2,j1,p4,p3,p6,p5,j7,j8,za,zb,a8apm)
      call a8treea(j2,j1,p3,p4,p5,p6,j7,j8,za,zb,a8amp)
      call a8treea(j2,j1,p4,p3,p5,p6,j7,j8,za,zb,a8app)
      do polg1=1,2
      do polg2=1,2
      fAB(polg1,polg2,minus,minus,minus)=
     &fAB(polg1,polg2,minus,minus,minus)+a8amm(polg1,polg2)

      fAB(polg1,polg2,minus,pplus,pplus)=
     &fAB(polg1,polg2,minus,pplus,pplus)+a8app(polg1,polg2)

      fAB(polg1,polg2,minus,minus,pplus)=
     &fAB(polg1,polg2,minus,minus,pplus)+a8amp(polg1,polg2)

      fAB(polg1,polg2,minus,pplus,minus)=
     &fAB(polg1,polg2,minus,pplus,minus)+a8apm(polg1,polg2)
      enddo
      enddo

c---Other color order swap (7<-->8)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,za,zb,a8amm)
      call a8treea(j2,j1,p4,p3,p6,p5,j8,j7,za,zb,a8apm)
      call a8treea(j2,j1,p3,p4,p5,p6,j8,j7,za,zb,a8amp)
      call a8treea(j2,j1,p4,p3,p5,p6,j8,j7,za,zb,a8app)
      do polg1=1,2
      do polg2=1,2
      fBA(polg1,polg2,minus,minus,minus)=
     &fBA(polg1,polg2,minus,minus,minus)+a8amm(polg2,polg1)

       fBA(polg1,polg2,minus,pplus,pplus)=
     &fBA(polg1,polg2,minus,pplus,pplus)+a8app(polg2,polg1)

      fBA(polg1,polg2,minus,minus,pplus)=
     &fBA(polg1,polg2,minus,minus,pplus)+a8amp(polg2,polg1)

      fBA(polg1,polg2,minus,pplus,minus)=
     &fBA(polg1,polg2,minus,pplus,minus)+a8apm(polg2,polg1)
      enddo
      enddo

c-- opposite quark helicity by taking complex conjugation (za<->zb)
c-- and swapping all helicities
c      call a8treeb(j2,j1,j4,j3,j5,j6,j7,j8,zb,za,a8b)
c      fAB(2,polg1,polg2,mplus)=+a8b(3-polg1,3-polg2)
      call a8treea(j2,j1,p4,p3,p5,p6,j7,j8,zb,za,a8amm)
      call a8treea(j2,j1,p3,p4,p5,p6,j7,j8,zb,za,a8apm)
      call a8treea(j2,j1,p4,p3,p6,p5,j7,j8,zb,za,a8amp)
      call a8treea(j2,j1,p3,p4,p6,p5,j7,j8,zb,za,a8app)
      do polg1=1,2
      do polg2=1,2
      fAB(polg1,polg2,pplus,minus,minus)=
     &fAB(polg1,polg2,pplus,minus,minus)+a8amm(3-polg1,3-polg2)

      fAB(polg1,polg2,pplus,pplus,pplus)=
     &fAB(polg1,polg2,pplus,pplus,pplus)+a8app(3-polg1,3-polg2)

      fAB(polg1,polg2,pplus,minus,pplus)=
     &fAB(polg1,polg2,pplus,minus,pplus)+a8amp(3-polg1,3-polg2)

      fAB(polg1,polg2,pplus,pplus,minus)=
     &fAB(polg1,polg2,pplus,pplus,minus)+a8apm(3-polg1,3-polg2)
      enddo
      enddo

c---Other color order swap (7<-->8)
      call a8treea(j2,j1,p4,p3,p5,p6,j8,j7,zb,za,a8amm)
      call a8treea(j2,j1,p3,p4,p5,p6,j8,j7,zb,za,a8apm)
      call a8treea(j2,j1,p4,p3,p6,p5,j8,j7,zb,za,a8amp)
      call a8treea(j2,j1,p3,p4,p6,p5,j8,j7,zb,za,a8app)
      do polg1=1,2
      do polg2=1,2
      fBA(polg1,polg2,pplus,minus,minus)=
     &fBA(polg1,polg2,pplus,minus,minus)+a8amm(3-polg2,3-polg1)

      fBA(polg1,polg2,pplus,pplus,pplus)=
     &fBA(polg1,polg2,pplus,pplus,pplus)+a8app(3-polg2,3-polg1)

      fBA(polg1,polg2,pplus,minus,pplus)=
     &fBA(polg1,polg2,pplus,minus,pplus)+a8amp(3-polg2,3-polg1)

      fBA(polg1,polg2,pplus,pplus,minus)=
     &fBA(polg1,polg2,pplus,pplus,minus)+a8apm(3-polg2,3-polg1)
      enddo
      enddo


      return
      end

