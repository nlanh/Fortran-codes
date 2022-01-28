program Legendre

implicit none

real*16 x,pleg,xmin,xmax,step

integer*16 j,n,k,n2,fac1,fac2,fac3,fac4
integer ix,ni

open(unit=1, file='pleg.dat')

print*, 'Enter n (order of legendre polynomial, hint: 2)'
read*, n
print*, 'Enter xmin (maximum range of x for output, hint: -1)'
read*, xmin
print*, 'Enter xmax (maximum range of x for output, hint: 1)'
read*, xmax
print*, 'Enter step (step of x, hint: 0.01)'
read*, step

if(mod(n,2).eq.0)n2=n/2
if(mod(n,2).eq.1)n2=(n-1)/2

ni = int((xmax-xmin+0.00000001)/step)

do ix=0,ni
	pleg=0
	x=xmin+ix*step
	do k=0,n2
	    fac1=1
		do j=1,k
			fac1=fac1*j
		end do

		fac2=1
		if ((n-k).lt.0) then
			do j=n-k,-1
				fac2=fac2*j
			end do
		end if
		if ((n-k).gt.0) then
			do j=1,n-k
				fac2=fac2*j
			end do
		end if

		fac3=1
		if ((n-2*k).lt.0) then
			do j=n-2*k,-1
				fac3=fac3*j
			end do
		end if
		if ((n-2*k).gt.0) then
			do j=1,n-2*k
				fac3=fac3*j
			end do
		end if

		fac4=1
		if ((2*n-2*k).lt.0) then
			do j=2*n-2*k,-1
				fac4=fac4*j
			end do
		end if
		if ((2*n-2*k).gt.0) then
			do j=1,2*n-2*k
				fac4=fac4*j
			end do
		end if

		pleg = pleg + (-1)**k*fac4*x**(n-2*k)/(2**n*fac1*fac2*fac3)
	end do
	write(1,'((f8.2),(f20.4))') x, pleg
end do

end program
