!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c--- LO amplitudes for the process 0 -> q(i1)+qb(i2)+g(i3)+g(i4)+gamma(i5)
      subroutine qqbgg_ga(i1,i2,i3,i4,i5,za,zb,qqb_gagg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'mmsq_cs_ga.f'

      integer:: i1,i2,i3,i4,i5,h1,h2,h3,h4
      real(dp):: qqb_gagg
      complex(dp) :: amp34(2,2,2,2),amp43(2,2,2,2),
     & ampQED(2,2,2,2)


      call qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp34)
      call qqbgg_ga_amp(i1,i2,i4,i3,i5,za,zb,amp43)

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  ampQED(h1,h2,h3,h4)=amp34(h1,h2,h3,h4)+
     &                 amp43(h1,h3,h2,h4)
               enddo
            enddo
         enddo
      enddo

c----- sum over polarizations
      mmsq_cs_ga(:)=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
       mmsq_cs_ga(1)= mmsq_cs_ga(1) +real(amp34(h1,h2,h3,h4)
     &                  *conjg(amp34(h1,h2,h3,h4)),dp)
       mmsq_cs_ga(2)= mmsq_cs_ga(2) +real(amp43(h1,h2,h3,h4)
     &                  *conjg(amp43(h1,h2,h3,h4)),dp)
       mmsq_cs_ga(0)= mmsq_cs_ga(0) +real(ampQED(h1,h2,h3,h4)
     &                  *conjg(ampQED(h1,h2,h3,h4)),dp)
      enddo
      enddo
      enddo
      enddo

      mmsq_cs_ga(0)=-1._dp/xn**2*mmsq_cs_ga(0)
      qqb_gagg=mmsq_cs_ga(1)+mmsq_cs_ga(2)+mmsq_cs_ga(0)
c      write(6,*) 'lord: amp 34 bit',mmsq_cs_ga(1)
c      write(6,*) 'lord: amp 43 bit',mmsq_cs_ga(2)
c      write(6,*) 'lord: QED bit',mmsq_cs_ga(0)
c      write(6,*) 'total ',qqb_gagg
      return
      end

      subroutine qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'

      integer:: i1,i2,i3,i4,i5,h2,h3,h4
      complex(dp)::amp(2,2,2,2),qqbgg_gaMHV,qqbgg_gMHVadj,qqbgg_gMHV,
     & sign

 !     qqb_gagg=zip

      amp(2,2,2,2)=czip
      amp(1,1,1,1)=czip
      amp(1,2,2,2)=czip
      amp(2,1,1,1)=czip

      amp(1,2,2,1)=qqbgg_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp(1,2,1,2)=qqbgg_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      amp(1,1,2,2)=qqbgg_gMHV(i1,i2,i3,i4,i5,za,zb)

      amp(1,1,1,2)=-qqbgg_gaMHV(i2,i1,i4,i3,i5,zb,za)
      amp(1,2,1,1)=-qqbgg_gMHVadj(i2,i1,i4,i3,i5,zb,za)
      amp(1,1,2,1)=-qqbgg_gMHV(i2,i1,i4,i3,i5,zb,za)

c      amp(2,1,1,2)=-qqbgg_gaMHV(i1,i2,i3,i4,i5,zb,za)
c      amp(2,1,2,1)=-qqbgg_gMHVadj(i1,i2,i3,i4,i5,zb,za)
c      amp(2,2,1,1)=-qqbgg_gMHV(i1,i2,i3,i4,i5,zb,za)

c      amp(2,2,2,1)=qqbgg_gaMHV(i2,i1,i4,i3,i5,za,zb)
c      amp(2,1,2,2)=qqbgg_gMHVadj(i2,i1,i4,i3,i5,za,zb)
c      amp(2,2,1,2)=qqbgg_gMHV(i2,i1,i4,i3,i5,za,zb)

c      fullamp(:,:,:,:)=amp(:,:,:,:)

c for speed, use trick to get remaining amplitudes; fix sign by explicit check of crossing
      if (((i1<3).and.(i2<3)) .or. ((i1>2).and.(i2>2))) then
        sign=+cone
      else
        sign=-cone
      endif

      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(2,h2,h3,h4)=sign*conjg(amp(1,3-h2,3-h3,3-h4))
      enddo
      enddo
      enddo

c      qqb_gagg=zip
c      do h1=1,2
c         do h2=1,2
c            do h3=1,2
c               do h4=1,2
c                  write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4),fullamp(h1,h2,h3,h4),amp(h1,h2,h3,h4)/fullamp(h1,h2,h3,h4)
c               enddo
c            enddo
c         enddo
c      enddo
c      pause

c      qqb_gagg=zip
c      do h1=1,2
c         do h2=1,2
c            do h3=1,2
c               do h4=1,2
c                  write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4)
c                  qqb_gagg=qqb_gagg
c     &                +real(amp(h1,h2,h3,h4)*conjg(amp(h1,h2,h3,h4)),dp)
c               enddo
c            enddo
c         enddo
c      enddo
c      pause
      return
      end

      function qqbgg_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) qqbgg_gaMHV
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gaMHV=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
      return
      end


      function qqbgg_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) qqbgg_gMHV
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gMHV=-(za(i1,i3)**3
     &     /(za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)))
      return
      end

      function qqbgg_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp)  qqbgg_gMHVadj
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gMHVadj=(za(i4,i1)**2*za(i4,i2))
     & /(za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4))
      return
      end
