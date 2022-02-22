c	Runge Kutta cho he cac phuong trinh vi phan cap mot
      PROGRAM RK2
      IMPLICIT none
      
c     N: so luong phuong trinh, nsteps: so diem, tstep: buoc chia
c     y(1): vi tri dau, y(2): van toc dau
      REAL*8 t, tstep, y(5)
      INTEGER N,j,nsteps
      N = 2
      nsteps = 100
      tstep = 0.1
      y(1) = 1.0
      y(2) = 0.0
      
c     Mo file
      OPEN(6, FILE = 'RK4.dat')
      WRITE (6,'(2F15.5)') 0.d0, y(1)

c     Vong lap n buoc cua thuat toan Runge-Kutta
      DO 60 j = 1, nsteps
       t=j*tstep
       CALL rk4(t, y, tstep, N)
       WRITE (6,'(2F15.5)') t, y(1)
   60 CONTINUE
      CLOSE(6)
      STOP
      END
c--------------------Ket thuc chuong trinh chinh------------------------

c     Chuong trinh con Runge-Kutta 
      SUBROUTINE rk4(t, y, tstep, N)
      IMPLICIT none

      REAL*8 DERIV, h, t, tstep, y(5) 
      REAL*8 k1(5), k2(5),k3(5), k4(5), temp1(5), temp2(5), temp3(5)
      INTEGER i, N
      h = tstep/2.0

      DO 10 i = 1,N
       k1(i) = tstep * DERIV(t, y, i)
       temp1(i) = y(i) + 0.5*k1(i)
   10 CONTINUE

      DO 20 i = 1,N
       k2(i) = tstep * DERIV(t+h, temp1, i)
       temp2(i) = y(i) + 0.5*k2(i)
   20 CONTINUE

      DO 30 i = 1,N
       k3(i) = tstep * DERIV(t+h, temp2, i)
       temp3(i) = y(i) + k3(i)
   30 CONTINUE

      DO 40 i = 1,N
       k4(i) = tstep * DERIV(t+tstep, temp3, i)
       y(i) = y(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
   40 CONTINUE
 
      RETURN
      END

c     Ham so dao ham
      FUNCTION DERIV(t,temp,i)
      IMPLICIT none
      
      REAL*8 DERIV, t, temp(5)
      INTEGER i

      IF (i .EQ. 1) DERIV = temp(2)
      IF (i .EQ. 2) DERIV = -temp(1)

      RETURN
      END
