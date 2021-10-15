program songuyento
implicit none
integer::l,m,n,i,j
	
print*,'Nhap n'
read*,n
print*,n,'bang tich cua cac thua so nguyen to sau nhan voi nhau'
	
m=int(n/2)
	do i=2,m
		if (mod(n,i).eq.0) then
			l=0
			do j=1,i
				if (mod(i,j).eq.0) then
					l = l+1
				end if
			end do
			if (l.eq.2) then
				do while (mod(n,i).eq.0)
					n=int(n/i)
					print*,i
				end do
			end if
		end if
	end do
end program
	
