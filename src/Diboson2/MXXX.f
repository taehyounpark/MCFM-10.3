!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine MXXX(order,qqb,qtype,p1,p2,p5,p6,p7,p8,qqbAj,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'nf.f'
      include 'nwz.f'
      include 'zcouple_cms.f'
      include 'srdiags.f'
      include 'masses.f'
      include 'ABCF.f'
      include 'coupfac_vv.f'
      include 'Qform.f'
      complex(dp):: amp(2,2,2),
     & ampA(2,2,2),ampB(2,2,2),ampC(2,2,2),ampF(2,2,2),ampS(2,2,2),
     & qqb5678(2,2,2),qqb7856(2,2,2),A6treeA,
     & Fa(5:8,2),Fb(5:8),qqbAj(0:2,4,10),
     & prop12,prop34,prop56
      integer p1,p2,p5,p6,p7,p8,h,h12,h78
      integer:: order,qtype
      integer,parameter:: dntype=1,uptype=2
      integer:: qtype2,sig
      logical:: qqb
c      complex(dp):: A6b_1,A6b_2,A6b_3,A6b_4
      complex(dp):: v2(2),cl1,cl2,en1,en2
      real(dp):: wwflag,q1

      q1=rq1

c      FAC=-two*(zesq)**2/zxw

      if     (nwz==-1) then
        cl1=cone
        cl2=czip
        en1=zle
        en2=zln
      elseif (nwz==+1) then
        cl1=czip
        cl2=cone
        en1=zln
        en2=zle
      endif
      wwflag=1._dp
      v2(1)=zl1
      v2(2)=zr1
c      cotw=sqrt((one-zxw)/zxw)


      if (kcase == kWZbbar) then
        if (qtype == uptype) then
            qtype2=dntype
            sig=-1
        elseif (qtype == dntype) then
            qtype2=uptype
            sig=+1
        endif

      else
      qtype2=qtype
      sig=+1
      endif

      call VVampfill(order,A,qtype,p1,p2,p5,p6,p7,p8,qqbAj,ampA)
      call VVampfill(order,B,qtype2,p1,p2,p5,p6,p7,p8,qqbAj,ampB)
      call VVampfill(order,C,qtype,p1,p2,p5,p6,p7,p8,qqbAj,ampC)
      call VVampfill(order,F,qtype,p1,p2,p5,p6,p7,p8,qqbAj,ampF)
      amp(:,:,:)=(ampA(:,:,:)+ampB(:,:,:)+ampC(:,:,:)+sig*ampF(:,:,:))

      if (srdiags) then
          amps(:,:,:)=czip
          if (kcase==kZZlept) then
          call qqbVVMXXXtreeS(.false.,p1,p2,p5,p6,p7,p8,qqb5678)
          call qqbVVMXXXtreeS(.true., p1,p2,p7,p8,p5,p6,qqb7856)

c     Multiply by sr diagrams by order-dependent factor
          ampS(:,:,:)=coupfac(5,qtype,:,:,:)*qqb5678(:,:,:)
     &               +coupfac(6,qtype,:,:,:)*qqb7856(:,:,:)

          elseif (kcase == kWWqqbr) then

          Fa(5,2)=A6treea(5,6,1,2,4,3,za,zb)
          Fa(6,2)=A6treea(3,4,5,6,2,1,za,zb)
          Fa(7,2)=A6treea(5,6,3,4,2,1,za,zb)
          Fa(8,2)=A6treea(3,4,1,2,6,5,za,zb)
          Fa(5,1)=A6treea(5,6,2,1,4,3,za,zb)
          Fa(6,1)=A6treea(3,4,5,6,1,2,za,zb)
          Fa(7,1)=A6treea(5,6,3,4,1,2,za,zb)
          Fa(8,1)=A6treea(3,4,2,1,6,5,za,zb)
          do h12=1,2
            if (qqb) then  !Calculating a qqb process, rather than qbq
               h=h12
            else
               h=3-h12
            endif
          ampS(h12,1,1)=coupfac(5,qtype,h12,1,1)*Fa(5,h)
     &                 +coupfac(6,qtype,h12,1,1)*Fa(6,h)
     &                 +coupfac(7,qtype,h12,1,1)*Fa(7,h)
     &                 +coupfac(8,qtype,h12,1,1)*Fa(8,h)
      enddo
      elseif (kcase == kWZbbar) then
          Fa(5,1)=A6treea(3,4,5,6,1,2,za,zb)
          Fa(5,2)=A6treea(3,4,6,5,1,2,za,zb)

          Fa(6,1)=A6treea(3,4,2,1,6,5,za,zb)
          Fa(6,2)=A6treea(3,4,2,1,5,6,za,zb)

          Fa(7,1)=A6treea(3,4,5,6,2,1,za,zb)
          Fa(7,2)=A6treea(3,4,6,5,2,1,za,zb)

          Fa(8,1)=A6treea(3,4,1,2,6,5,za,zb)
          Fa(8,2)=A6treea(3,4,1,2,5,6,za,zb)

          Fb(5)=A6treea(5,6,2,1,4,3,za,zb)
          Fb(6)=A6treea(5,6,3,4,1,2,za,zb)
          Fb(7)=A6treea(5,6,1,2,4,3,za,zb)
          Fb(8)=A6treea(5,6,3,4,2,1,za,zb)

          prop12=s(1,2)/cmplx(s(1,2)-wmass**2,wmass*wwidth,kind=dp)
          prop34=s(3,4)/cmplx(s(3,4)-wmass**2,wmass*wwidth,kind=dp)
          prop56=s(5,6)/cmplx(s(5,6)-zmass**2,zmass*zwidth,kind=dp)

          if (qqb) then
             do h78=1,2
             ampS(1,1,h78)=prop12*(
     &       (en1*Fa(5,h78)+en2*Fa(6,h78))*v2(h78)*prop56
     &       +q1*(-1._dp)*(cl1*Fa(5,h78)+cl2*Fa(6,h78))
     &       +(2-h78)*wwflag*0.5_dp/zxw*prop34*(cl1*Fb(5)+cl2*Fb(6)))
             enddo
          else
             do h78=1,2
             ampS(1,1,h78)=prop12*(
     &       (en1*Fa(7,h78)+en2*Fa(8,h78))*v2(h78)*prop56
     &       +q1*(-1._dp)*(cl1*Fa(7,h78)+cl2*Fa(8,h78))
     &       +(2-h78)*wwflag*0.5_dp/zxw*prop34*(cl1*Fb(7)+cl2*Fb(8)))
             enddo
          endif
      endif   ! End kcase if
      amp(:,:,:)=amp(:,:,:)+Qform(order)*ampS(:,:,:)
      endif       ! End srdiags if
      return
      end
