subroutine cconsensus(clusters, n, K, ans, m)
implicit none

integer, intent(IN)					:: n, K, m		! dimensions
integer, dimension(K, n), intent(IN)	:: clusters		! data matrix
integer, dimension(m), intent(OUT)	:: ans			! output array
integer									:: i, ii, j, kk	! locals

do kk = 1, K
	j = 1
	do i = 1, n-1
		do ii = i+1, n
			if(clusters(kk, i) * clusters(kk, ii) > 0) then
				if(clusters(kk, i)  .eq. clusters(kk, ii)) then
					ans(j) = ans(j) + 1
				end if	
			end if	
			j = j +1
		end do
	end do
end do
return
stop
end 
