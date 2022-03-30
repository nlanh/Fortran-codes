      program Legendre

      implicit real*8 (a-h,o-z)
      dimension x(10000),pleg(10000),pleg1(10000),pleg2(10000)

      open(unit=1, file='pleg.dat')

      print*, 'Enter n'
      read*, n
      print*, 'Enter xmin'
      read*, xmin
      print*, 'Enter xmax'
      read*, xmax
      print*, 'Enter step'
      read*, step

      ni = int((xmax-xmin+0.00000001)/step+1)
      
      do i = 1, ni
        x(i)=xmin+(i-1)*step
        pleg2(i)=1.d0
        pleg1(i)=x(i)
        pleg(i)=pleg1(i)
        if(n.eq.0)pleg(i)=pleg2(i)
        if(n.eq.1)pleg(i)=pleg1(i)
      end do
      
      do i=2,n
        do j=1,ni
            pleg(j)=((2*i-1)*x(j)*pleg1(j)-(i-1)*pleg2(j))/i
            pleg2(j)=pleg1(j)
            pleg1(j)=pleg(j)
        end do
      end do
      
      do i=1,ni
        write(1,'((f8.2),(f20.4))') x(i), pleg(i)
      end do

      end program
