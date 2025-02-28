      subroutine pvCcheck(rank,q1,q2,m0s,m1s,m2s,
     & FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/TRydef.f'
      include 'lib/TensorReduction/Include/TRmetric.f'
      include 'lib/TensorReduction/Include/pvverbose.f'
c      include 'pvbadpoint.f'
      integer n2,n3,n4,n5,n6,ep,nu,rank,epmin
      real(dp):: q1(4),q2(4),p2(4),f1,f2,Cacc
      real(dp):: pvSDDDD,pvSDDPP,pvSDDPK,pvSPKKK,pvSPPKK
      real(dp):: q1Dq1,q2Dq2,q1Dq2,s12,m0s,m1s,m2s
      real(dp):: sing2(-2:0),sing3(-2:0),sing4(-2:0),
     & sing5(-2:0),sing6(-2:0)
      complex(dp):: FB01(-2:0),FB11(y1max,-2:0),FB21(y2max,-2:0),
     & FB31(y3max,-2:0),FB41(y4max,-2:0),FB51(y5max,-2:0),
     & FB61(y6max,-2:0)
      complex(dp):: FB02(-2:0),FB12(y1max,-2:0),FB22(y2max,-2:0),
     & FB32(y3max,-2:0),FB42(y4max,-2:0),FB52(y5max,-2:0),
     & FB62(y6max,-2:0)
      complex(dp):: FB03(-2:0),FB13(y1max,-2:0),FB23(y2max,-2:0),
     & FB33(y3max,-2:0),FB43(y4max,-2:0),FB53(y5max,-2:0),
     & FB63(y6max,-2:0)
      complex(dp):: FB13a(y1max,-2:0),FB23a(y2max,-2:0),
     & FB33a(y3max,-2:0),FB43a(y4max,-2:0),FB53a(y5max,-2:0)
c     & ,FB63a(y6max,-2:0)
      complex(dp):: FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0),FC4(y4max,-2:0),FC5(y5max,-2:0),FC6(y6max,-2:0),
     & trhs,tq
      logical::failed
      parameter(epmin=0) ! Only check finite pieces
      include 'lib/TensorReduction/Include/cplx.h'      
      failed=.false.
      
      Cacc=1.e-8_dp

      q1Dq1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2Dq2=q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2
      q1Dq2=q1(4)*q2(4)-q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)
      s12=q1Dq1+q2Dq2-two*q1Dq2

c      write(6,'(a35,5(e12.5,a2))') 
c     . '(p1sq, p2sq, m1sq, m2sq, m3sq) = ( ',
c     .  q1Dq1,', ',s12,', ',m0s,', ',m1s,', ',m2s,' )'
      
      do ep=-2,0
      sing2(ep)=zip
      sing3(ep)=zip
      sing4(ep)=zip
      sing5(ep)=zip
      sing6(ep)=zip
      enddo
      do nu=1,4
      p2(nu)=q2(nu)-q1(nu)
      enddo  

      call pvBtensor(q2,m0s,m2s,FB01,FB11,FB21,FB31,FB41,FB51,FB61)
      call pvBtensor(q1,m0s,m1s,FB02,FB12,FB22,FB32,FB42,FB52,FB62)
      call pvBtensor(p2,m1s,m2s,FB03,FB13,FB23,FB33,FB43,FB53,FB63)

      if ((rank .eq. 2) .or. (rank .eq. 3))
     &  call pvswitch1(q1,FB03,FB13,FB13a)
      if ((rank .eq. 3) .or. (rank .eq. 4))
     &  call pvswitch2(q1,FB03,FB13,FB23,FB23a)
      if ((rank .eq. 4) .or. (rank .eq. 5))
     &  call pvswitch3(q1,FB03,FB13,FB23,FB33,FB33a)
      if ((rank .eq. 5) .or. (rank .eq. 6))
     &  call pvswitch4(q1,FB03,FB13,FB23,FB33,FB43,FB43a)
      if (rank .eq. 6)
     &  call pvswitch5(q1,FB03,FB13,FB23,FB33,FB43,FB53,FB53a)
c      call pvswitch6(q1,FB03,FB13,FB23,FB33,FB43,FB53,FB63,FB63a)

      f1=m1s-m0s-q1Dq1
      f2=m2s-m0s-q2Dq2

c--- check rank 1
      if (rank .eq. 1) then
      if (pvverbose) write(6,*) 'q1.FC1'
      do ep=epmin,0
      tq=  q1(4)*FC1(4,ep)
     &    -q1(1)*FC1(1,ep)
     &    -q1(2)*FC1(2,ep)
     &    -q1(3)*FC1(3,ep)
      trhs=
     & -half*(FB01(ep)-FB03(ep)+f1*FC0(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo

      if (pvverbose) write(6,*) 'q2.FC1'
      do ep=epmin,0
      tq=  q2(4)*FC1(4,ep)
     &    -q2(1)*FC1(1,ep)
     &    -q2(2)*FC1(2,ep)
     &    -q2(3)*FC1(3,ep)
      trhs=
     & -half*(FB02(ep)-FB03(ep)+f2*FC0(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      
      endif
      

c--- check for rank 2
      if (rank .eq. 2) then
      if (pvverbose) write(6,*) 'q1.FC2'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      tq=  q1(4)*FC2(y2(4,n2),ep)
     &    -q1(1)*FC2(y2(1,n2),ep)
     &    -q1(2)*FC2(y2(2,n2),ep)
     &    -q1(3)*FC2(y2(3,n2),ep)
      trhs=
     & -half*(FB11(n2,ep)-FB13a(n2,ep)+f1*FC1(n2,ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q2.FC2'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      tq= q2(4)*FC2(y2(4,n2),ep)
     &   -q2(1)*FC2(y2(1,n2),ep)
     &   -q2(2)*FC2(y2(2,n2),ep)
     &   -q2(3)*FC2(y2(3,n2),ep)
      trhs=
     & -half*(FB12(n2,ep)-FB13a(n2,ep)+f2*FC1(n2,ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
c      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC2'
      sing2(0)=-half 
c      do ep=epmin,0 
      ep=0
      tq=FC2(y2(4,4),ep)
     & -FC2(y2(1,1),ep)
     & -FC2(y2(2,2),ep)
     & -FC2(y2(3,3),ep)
     & -m0s*FC0(ep)-FB03(ep)
      trhs=
     & +cplx1(sing2(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
c      enddo

      endif
      

      if (rank .eq. 3) then
      if (pvverbose) write(6,*) 'q1.FC3'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      tq=  +q1(4)*FC3(y3(4,n2,n3),ep)
     &     -q1(1)*FC3(y3(1,n2,n3),ep)
     &     -q1(2)*FC3(y3(2,n2,n3),ep)
     &     -q1(3)*FC3(y3(3,n2,n3),ep)
      trhs=
     &    -half*(FB21(y2(n2,n3),ep)
     &           -FB23a(y2(n2,n3),ep)+f1*FC2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      
      if (pvverbose) write(6,*) 'q2.FC3'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      tq=+q2(4)*FC3(y3(4,n2,n3),ep)
     &     -q2(1)*FC3(y3(1,n2,n3),ep)
     &     -q2(2)*FC3(y3(2,n2,n3),ep)
     &     -q2(3)*FC3(y3(3,n2,n3),ep)
      trhs=
     &    -half*(FB22(y2(n2,n3),ep)
     &           -FB23a(y2(n2,n3),ep)+f2*FC2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC3'
      do ep=epmin,0
      do n3=1,4
      sing3(0)=+one/six*(q1(n3)+q2(n3))
 
      tq=FC3(y3(4,4,n3),ep)
     & -FC3(y3(1,1,n3),ep)
     & -FC3(y3(2,2,n3),ep)
     & -FC3(y3(3,3,n3),ep)
     & -m0s*FC1(n3,ep)-FB13a(n3,ep)
      trhs=
     & +cplx1(sing3(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
  
      endif
      
      if (rank .eq. 4) then
      if (pvverbose) write(6,*) 'q1.FC4'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      tq=  q1(4)*FC4(y4(4,n2,n3,n4),ep)
     &    -q1(1)*FC4(y4(1,n2,n3,n4),ep)
     &    -q1(2)*FC4(y4(2,n2,n3,n4),ep)
     &    -q1(3)*FC4(y4(3,n2,n3,n4),ep)
      trhs=
     & -half*(FB31(y3(n2,n3,n4),ep)
     &       -FB33a(y3(n2,n3,n4),ep)+f1*FC3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FC4'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      tq=  q2(4)*FC4(y4(4,n2,n3,n4),ep)
     &    -q2(1)*FC4(y4(1,n2,n3,n4),ep)
     &    -q2(2)*FC4(y4(2,n2,n3,n4),ep)
     &    -q2(3)*FC4(y4(3,n2,n3,n4),ep)
      trhs=
     & -half*(FB32(y3(n2,n3,n4),ep)
     &       -FB33a(y3(n2,n3,n4),ep)+f2*FC3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC4'
c      s12=q1Dq1+q2Dq2-2*q1Dq2
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      sing4(0)=
     &  -one/12._dp*q1(n3)*q1(n4)
     &  -one/12._dp*q2(n3)*q2(n4)
     &  -one/24._dp*(q1(n3)*q2(n4)+q1(n4)*q2(n3))
     &  +(one/48._dp*(s12+q1Dq1+q2Dq2)
     &  -one/12._dp*(m0s+m1s+m2s))*g(n3,n4)
      tq= 
     &     FC4(y4(4,4,n3,n4),ep)
     &    -FC4(y4(1,1,n3,n4),ep)
     &    -FC4(y4(2,2,n3,n4),ep)
     &    -FC4(y4(3,3,n3,n4),ep)
     &    -m0s*FC2(y2(n3,n4),ep)
     &    -FB23a(y2(n3,n4),ep)
      trhs=+cplx1(sing4(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
         
      endif
      
      if (rank .eq. 5) then
      if (pvverbose) write(6,*) 'q1.FC5'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      tq=  q1(4)*FC5(y5(4,n2,n3,n4,n5),ep)
     &    -q1(1)*FC5(y5(1,n2,n3,n4,n5),ep)
     &    -q1(2)*FC5(y5(2,n2,n3,n4,n5),ep)
     &    -q1(3)*FC5(y5(3,n2,n3,n4,n5),ep)
      trhs=
     & -half*(FB41(y4(n2,n3,n4,n5),ep)
     &       -FB43a(y4(n2,n3,n4,n5),ep)+f1*FC4(y4(n2,n3,n4,n5),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FC5'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      tq=
     &   q2(4)*FC5(y5(4,n2,n3,n4,n5),ep)
     &  -q2(1)*FC5(y5(1,n2,n3,n4,n5),ep)
     &  -q2(2)*FC5(y5(2,n2,n3,n4,n5),ep)
     &  -q2(3)*FC5(y5(3,n2,n3,n4,n5),ep)
      trhs=
     & -half*(FB42(y4(n2,n3,n4,n5),ep)
     &       -FB43a(y4(n2,n3,n4,n5),ep)+f2*FC4(y4(n2,n3,n4,n5),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC5'
c      s12=q1Dq1+q2Dq2-2*q1Dq2
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      do n5=n4,4
      sing5(0)=
     &  -one/240._dp
     &  *(two*s12-five*m0s+two*(q1Dq1-five*m1s)+(q2Dq2-five*m2s))
     &  *(+g(n3,n4)*q1(n5)+g(n4,n5)*q1(n3)+g(n5,n3)*q1(n4))
     &  -one/240._dp
     &  *(two*s12-five*m0s+(q1Dq1-five*m1s)+two*(q2Dq2-five*m2s))
     &  *(+g(n3,n4)*q2(n5)+g(n4,n5)*q2(n3)+g(n5,n3)*q2(n4))
     &  + one/20._dp*q1(n3)*q1(n4)*q1(n5)
     &  + one/20._dp*q2(n3)*q2(n4)*q2(n5)
     & +one/60._dp
     & *(q1(n3)*q2(n4)*q2(n5)+q1(n4)*q2(n5)*q2(n3)+q1(n5)*q2(n3)*q2(n4))
     & +one/60._dp
     & *(q2(n3)*q1(n4)*q1(n5)+q2(n4)*q1(n5)*q1(n3)+q2(n5)*q1(n3)*q1(n4))

      tq=FC5(y5(4,4,n3,n4,n5),ep)
     &    -FC5(y5(1,1,n3,n4,n5),ep)
     &    -FC5(y5(2,2,n3,n4,n5),ep)
     &    -FC5(y5(3,3,n3,n4,n5),ep)
     &    -m0s*FC3(y3(n3,n4,n5),ep)
     &    -FB33a(y3(n3,n4,n5),ep)
      trhs=+cplx1(sing5(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo

      endif

      
      if (rank .eq. 6) then
      if (pvverbose) write(6,*) 'q1.FC6'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      tq=q1(4)*FC6(y6(4,n2,n3,n4,n5,n6),ep)
     &  -q1(1)*FC6(y6(1,n2,n3,n4,n5,n6),ep)
     &  -q1(2)*FC6(y6(2,n2,n3,n4,n5,n6),ep)
     &  -q1(3)*FC6(y6(3,n2,n3,n4,n5,n6),ep)

      trhs=
     &   -half*(FB51(y5(n2,n3,n4,n5,n6),ep)
     & -FB53a(y5(n2,n3,n4,n5,n6),ep)+f1*FC5(y5(n2,n3,n4,n5,n6),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      
      if (pvverbose) write(6,*) 'q2.FC6'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      tq =q2(4)*FC6(y6(4,n2,n3,n4,n5,n6),ep)
     &   -q2(1)*FC6(y6(1,n2,n3,n4,n5,n6),ep)
     &   -q2(2)*FC6(y6(2,n2,n3,n4,n5,n6),ep)
     &   -q2(3)*FC6(y6(3,n2,n3,n4,n5,n6),ep)
    
      trhs=
     & -half*(FB52(y5(n2,n3,n4,n5,n6),ep)
     & -FB53a(y5(n2,n3,n4,n5,n6),ep)+f2*FC5(y5(n2,n3,n4,n5,n6),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC6'
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4

      sing6(0)=
     & -pvSDDDD(n3,n4,n5,n6)/2880._DP
     & *(two*s12**2-six*s12*m0s+30._dp*m0s**2
     & +two*s12*(q1Dq1-six*m1s+q2Dq2-six*m2s)
     & -six*m0s*(two*q1Dq1-five*m1s+two*q2Dq2-five*m2s)
     & +two*(q1Dq1**2-six*q1Dq1*m1s+15._dp*m1s**2)
     & +two*(q2Dq2**2-six*q2Dq2*m2s+15._dp*m2s**2)
     & +(q1Dq1*q2Dq2-six*q1Dq1*m2s+15._dp*m1s*m2s)
     & +(q1Dq1*q2Dq2-six*q2Dq2*m1s+15._dp*m1s*m2s))

     & +pvSDDPP(n3,n4,n5,n6,q1)/720._dp
     & *(three*s12-six*m0s+three*(q1Dq1-six*m1s)+(q2Dq2-six*m2s))

     & +pvSDDPP(n3,n4,n5,n6,q2)/720._dp
     & *(three*s12-six*m0s+(q1Dq1-six*m1s)+three*(q2Dq2-six*m2s))

     & +pvSDDPK(n3,n4,n5,n6,q1,q2)/720._dp
     & *(two*s12-three*m0s+(q1Dq1-six*m1s)+(q2Dq2-six*m2s))

     & -q1(n3)*q1(n4)*q1(n5)*q1(n6)/30._dp
     & -q2(n3)*q2(n4)*q2(n5)*q2(n6)/30._dp
     & -pvSPKKK(n3,n4,n5,n6,q1,q2)/120._dp
     & -pvSPKKK(n3,n4,n5,n6,q2,q1)/120._dp
     & -pvSPPKK(n3,n4,n5,n6,q2,q1)/180._dp

      tq = 
     & +FC6(y6(4,4,n3,n4,n5,n6),ep)
     & -FC6(y6(1,1,n3,n4,n5,n6),ep)
     & -FC6(y6(2,2,n3,n4,n5,n6),ep)
     & -FC6(y6(3,3,n3,n4,n5,n6),ep)
     &  -m0s*FC4(y4(n3,n4,n5,n6),ep)


      trhs=
     & -FB43a(y4(n3,n4,n5,n6),ep)
     & +cplx1(sing6(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      
      endif

      return
      end
