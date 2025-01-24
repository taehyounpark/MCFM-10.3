!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use pentbox_generic
      implicit none
      include 'Inc/Cred.f'
      integer,dimension(4,4,4)::I5to4
      integer::ppp,p1,p2,p3,p4
      real(dp)::mtsq

      I5to4(:,:,:)=0
      ppp=0
      do p1=1,4
      do p2=1,4
      if (p2==p1) cycle
      do p3=1,4
      if ((p3==p2).or.(p3==p1)) cycle
      ppp=ppp+1
      I5to4(p1,p2,p3)=ppp
      p4=10-p1-p2-p3
      call pentbox(p1,p2,p3,p4,mtsq,Cred(:,ppp))
      enddo
      enddo
      enddo
      return

