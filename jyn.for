      program yjn

      implicit real*8 (a-h,o-z)

      integer i,j,n,ni
      dimension x(10000),ayn(10000),dyn(10000),ajn(10000),djn(10000)
      
      open(unit=1, file='jyn.dat')

      print*, 'Enter n'
      read*, n
      print*, 'Enter xmax'
      read*, xmax
      print*, 'Enter step'
      read*, step

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
      
      subroutine der1(y,x,n,h)                                            
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n)
      n3=n-3
      x(1)=(-147.d0*y(1)+360.d0*y(2)-450.d0*y(3)+400.d0*y(4)-225.d0*y(5)
     1+72.d0*y(6)-10.d0*y(7))/60.d0                                     
      x(2)=(-10.d0*y(1)-77.d0*y(2)+150.d0*y(3)-100.d0*y(4)+50.d0*y(5)-15
     1.d0*y(6)+2.d0*y(7))/60.d0                                         
      x(3)=(2.d0*y(1)-24.d0*y(2)-35.d0*y(3)+80.d0*y(4)-30.d0*y(5)+8.d0*y
     1(6)-y(7))/60.d0                                                   
      do i=4,n3                                                       
      x(i)=(45.d0*(y(i+1)-y(i-1))-9.d0*(y(i+2)-y(i-2))+y(i+3)-y(i-3))/60
     1.d0
      end do                                                               
      x(n-2)=(y(n-6)-8.d0*y(n-5)+30.d0*y(n-4)-80.d0*y(n3)+35.d0*y(n-2)+2
     14.d0*y(n-1)-2.d0*y(n))/60.d0                                      
      x(n-1)=(-2.d0*y(n-6)+15.d0*y(n-5)-50.d0*y(n-4)+100.d0*y(n3)-150.d0
     1*y(n-2)+77.d0*y(n-1)+10.d0*y(n))/60.d0                            
      x(n)=(10.d0*y(n-6)-72.d0*y(n-5)+225.d0*y(n-4)-400.d0*y(n3)+450.d0*
     1y(n-2)-360.d0*y(n-1)+147.d0*y(n))/60.d0
      do i = 1, n
       x(i) = x(i)/h
      end do
      return
      end
