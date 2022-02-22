c	  Runge Kutta cho phuong trinh vi phan cap mot
      PROGRAM RK41
      IMPLICIT none
      
c     nsteps: so diem, tstep: buoc chia, y: gia tri y(0)
      REAL*8 t, y, tstep
      INTEGER j,nsteps
      nsteps = 100
      tstep = 0.1
      y = 1.0

c     Mo file
      OPEN(6, FILE='RK4.dat')
c
c 	  Vong lap nsteps buoc cua thuat toan Runga-Kutta
      DO 60 j = 0, nsteps
       t = j*tstep
       WRITE (6,'(2F15.5)') t,y
       CALL rk4(t, y, tstep)
   60 CONTINUE
   
      CLOSE(6)
      STOP
      END
c--------------------Ket thuc chuong trinh chinh------------------------

c     Chuong trinh con cho Runge-Kutta bac bon 
      SUBROUTINE rk4(t, y, tstep)
      IMPLICIT none
      REAL*8 DERIV, h, t, tstep, y 
      REAL*8 k1, k2, k3, k4
      h = tstep/2.0
      k1 = tstep * DERIV(t, y)
      k2 = tstep * DERIV(t + h, y + k1/2.0)
      k3 = tstep * DERIV(t + h, y + k2/2.0)
      k4 = tstep * DERIV(t + tstep, y + k3)
      y = y + (k1 + (2.*(k2 + k3)) + k4)/6.0
      RETURN
      END

c Ham so f(x,y) = y'
      FUNCTION DERIV(t, temp)
      IMPLICIT none

      REAL*8 DERIV, t, temp

      DERIV = t + temp
      
      RETURN
      END
