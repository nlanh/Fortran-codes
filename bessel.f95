program Bessel

implicit none

real*16 x,besj,xmax,step

integer*16 i,j,n,s,nmax,ni,fac1,fac2
open(unit=1, file='besj.dat')
open(unit=2, file='frho.dat')
open(unit=3, file='frho.dat')

print*, 'Enter n (order of bessel function, hint: 0)'
read*, n
print*, 'Enter nmax (infinity index for series, hint: 20)'
read*, nmax
print*, 'Enter xmax (maximum range of x for output, hint: 10)'
read*, xmax
print*, 'Enter step (step of x)'
read*, step

ni = int((xmax+0.00000001)/step)


do i=0,ni
	besj=0
	x=i*step
	do s=0,nmax
		fac1=1
		if (s.ne.0) then
		    fac1=1
			do j=1,s
				fac1=fac1*j
			end do
		end if
		fac2=1
		if ((n+s).ne.0) then
			fac2=1
			do j=1,n+s
				fac2=fac2*j
			end do
		end if
		besj = besj + (x/2)**(n+2*s)*(-1)**s/(fac1*fac2)
	end do
	write(1,'((f8.2),(f20.4))') x, besj
	
end do

end program
