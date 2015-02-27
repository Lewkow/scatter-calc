
SUBROUTINE legendre(x, n)
	USE current_energy, ONLY : leg

	!----------------------------
	!	Inputs	
	!----------------------------

	REAL(KIND=8)								:: x

	INTEGER											:: n

	!----------------------------
	!	Internal	
	!----------------------------

	INTEGER											:: j
	INTEGER											:: L
	REAL(KIND=8)								:: nn

	L = n+1

	ALLOCATE(leg(L))	

	leg(1) = 1.0D0
	leg(2) = x

	DO j=2,(L-1)
		nn       = REAL(j-1)
		leg(j+1) = ( (2.0D0*nn+1.0D0)/(nn+1.0D0) )*x*leg(j) - (nn/(nn+1.0D0))*leg(j-1) 
	END DO	

END SUBROUTINE legendre
