      SUBROUTINE prop(ip,xp,p,R1,R2,xcos,npropin,nprt)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /mass/ xmass(3),g1,g2,zembed,zbisc
!      dimension p(npropin),jt(npropin+3),pv(npropin+3),hc(nprt)
!      dimension ft(npropin+3),ut(npropin+3),jb(nprt),pb(nprt)
      double precision, dimension (npropin)  :: p
      double precision, dimension (32)  :: pv
      double precision, dimension (npropin)  :: ft
      double precision, dimension (npropin+3)  :: ut
      double precision, dimension (nprt)  :: pb
      double precision, dimension (nprt)  :: hc    
      integer, dimension (32)  :: jt
      integer, dimension (nprt)  :: jb
      DATA RZ/0.9576257   /,RHO/75.48992362  /
      DATA TOANG/0.5291772/
      DATA X1/1.0D0/,X0/0.0D0/,TINY/9.0D-15/,X2/2.0D0/
      real*8  :: xp(1:npropin),xp_t(1:npropin),tempx,deltax
      integer :: ip(1:npropin)
      real*8  :: fitfactordeltax=0.0001d0
      save icall, jt, pv

!      print *,'I am prop'

         F1= X1/G1
         F2= X1/G2
         F12= X1 - F1*F2
         P1= R1*(X1-F1)/(G2*F12)
         P2= R2*(X1-F2)/(G1*F12)
         S1= R1-P1
         S2= R2-P2
         Q1= SQRT(P1*P1 + S2*S2 + X2*P1*S2*XCOS)/(X1-G1)
         Q2= SQRT(P2*P2 + S1*S1 + X2*P2*S1*XCOS)/(X1-G2)
         Q3= SQRT(P1*P1 + P2*P2 - X2*P1*P2*XCOS)
         COST = (Q1*Q1 + Q2*Q2 - Q3*Q3)/(X2*Q1*Q2)
         THETA = ACOS(COST)

      THETAeq = (180.d0 - RHO)*3.141592654/180.d0
      Req= RZ/TOANG
      Y1 = Q1 - Req
      Y2 = Cos(THETA) - Cos(THETAeq)
      Y3 = Q2 - Req

      S1 = (Y1 + Y3)/2.0d0
      S2 =  Y2
      S3 = (Y1 - Y3)/2.0d0 
      !
      xp_t = xp
      !
      att1=Q1*0.5291772d0
      att2=Q2*0.5291772d0
      !
      call poten(xp_t,v,att1,att2,THETA)
      !
      ncol=0
      !
      p = 0 
      !
      do  i=1,npropin
        if (ip(i) > 0) then
          !
          vleft  = v
          vright = v
          !
          ncol=ncol+1
          tempx=xp_t(i)
          deltax=fitfactordeltax*abs(tempx)
          if (deltax .le. 1e-15) deltax=1e-6
          !
          xp_t(i)=tempx+deltax
          !
          call poten(xp_t,vright,att1,att2,THETA)
          !
          xp_t(i)=tempx-deltax
          !
          call poten(xp_t,vleft,att1,att2,THETA)
          !
          p(i)=(vright-vleft)/(2.0d0*deltax)
          !
          xp_t(i)=tempx
          !
        endif
      enddo ! --- ncol
      !
      !p(:) = p(:)*219474.624d0
      !
      return
      !
      end subroutine prop