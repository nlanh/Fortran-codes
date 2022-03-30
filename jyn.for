      program yjn

      implicit real*8 (a-h,o-z)

      integer i,j,n,ni
      dimension x(10000),ayn(10000),dyn(10000),ajn(10000),djn(10000)
      
      open(unit=1, file='y.dat') 
      open(unit=2, file='yjn.inp')

!print*, 'Nhap n'
      read(2,*) n
!print*, 'Nhap xmax'
      read(2,*) xmax
!print*, 'Nhap step'
      read(2,*) step

      ni = int((xmax+0.00000001)/step)
      
      do i = 1,ni
        x(i) = i*step
        ajn(i) = sin(x(i))/x(i)
        ayn(i) = cos(x(i))/x(i)
      end do
      
      do i = 1,n
        call der1(ajn,djn,ni,step)
        call der1(ayn,dyn,ni,step)
        do j = 1,ni
            ajn(j) = djn(j)
            ayn(j) = dyn(j)
        enddo
      end do
      
      do i = 1,ni
        ajn(i)=(-1)**n*ajn(i)
        ayn(i)=(-1)**(n+1)*ayn(i)
        write(1,'(f10.2,2f20.4)') x(i), ajn(i), ayn(i)
      enddo

      endprogram
      
      SUBROUTINE DER1(Y,X,N,H)                                            
C NUMERICAL DERIVATION OF THE FUNCTION Y KNOWN AT N POINTS.             
C THIS SUBROUTINE REQUIRES AT LEAST 7 POINTS.                           
C AND RETURNS THE VALUE OF H*(D(Y)/DR) IN X.                            
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N)
      N3=N-3
      X(1)=(-147.D0*Y(1)+360.D0*Y(2)-450.D0*Y(3)+400.D0*Y(4)-225.D0*Y(5)
     1+72.D0*Y(6)-10.D0*Y(7))/60.D0                                     
      X(2)=(-10.D0*Y(1)-77.D0*Y(2)+150.D0*Y(3)-100.D0*Y(4)+50.D0*Y(5)-15
     1.D0*Y(6)+2.D0*Y(7))/60.D0                                         
      X(3)=(2.D0*Y(1)-24.D0*Y(2)-35.D0*Y(3)+80.D0*Y(4)-30.D0*Y(5)+8.D0*Y
     1(6)-Y(7))/60.D0                                                   
      DO I=4,N3                                                       
      X(I)=(45.D0*(Y(I+1)-Y(I-1))-9.D0*(Y(I+2)-Y(I-2))+Y(I+3)-Y(I-3))/60
     1.D0
      END DO                                                               
      X(N-2)=(Y(N-6)-8.D0*Y(N-5)+30.D0*Y(N-4)-80.D0*Y(N3)+35.D0*Y(N-2)+2
     14.D0*Y(N-1)-2.D0*Y(N))/60.D0                                      
      X(N-1)=(-2.D0*Y(N-6)+15.D0*Y(N-5)-50.D0*Y(N-4)+100.D0*Y(N3)-150.D0
     1*Y(N-2)+77.D0*Y(N-1)+10.D0*Y(N))/60.D0                            
      X(N)=(10.D0*Y(N-6)-72.D0*Y(N-5)+225.D0*Y(N-4)-400.D0*Y(N3)+450.D0*
     1Y(N-2)-360.D0*Y(N-1)+147.D0*Y(N))/60.D0
      DO I = 1, N
       X(I) = X(I)/H
      END DO
      RETURN
      END
