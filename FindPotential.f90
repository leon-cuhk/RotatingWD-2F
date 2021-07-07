SUBROUTINE GRAVRHO
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Add the density !
IF(DM_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV	
			rho (i,j) = rho1(i,j) + rho2(i,j)
		END DO
	END DO
ELSE
	DO i = 1, KDIV
		DO j = 1, NDIV	
			rho (i,j) = rho2(i,j)
		END DO
	END DO
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDPOTENTIAL
USE DEFINITION 
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, m, l

! Initialization !
d1 = 0.0D0
d2 = 0.0D0
phi = 0.0D0

! Find d1 !
DO m = 0, lmax
	DO k = 1, NDIV
		DO i = 1, KDIV-2, 2
			d1(k,m) = d1(k,m) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))*(pn(i,2*m)*rho(i,k) + &
			4.0D0*pn(i+1,2*m)*rho(i+1,k) + pn(i+2,2*m)*rho(i+2,k))
		END DO
	END DO
END DO
	
! Find d2 !
DO m = 0, lmax
	DO j = 1, NDIV
		DO k = 1, NDIV-2, 2
			d2(m,j) = d2(m,j) + (1.0D0/6.0D0)*(r(k+2) - r(k))*(fp(k,j,2*m)*d1(k,m) + &
			4.0D0*fp(k+1,j,2*m)*d1(k+1,m) + fp(k+2,j,2*m)*d1(k+2,m))
		END DO
	END DO	
END DO

! Find potential !
DO i = 1, KDIV
	DO j = 1, NDIV
		DO m = 0, LMAX
			phi(i,j) = phi(i,j) + (-4.0D0*pi)*d2(m,j)*pn(i,2*m)
		END DO
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDFPOTENTIAL
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Initialize !
fhat2 = 0.0D0

! Assign for DM !
IF(dm_flag == 1) THEN
	fhat1 = 0.0D0
	DO i = 1, KDIV
		DO j = 1, NDIV
			fhat1(i,j) = -phi(i,j)
		END DO
	END DO
END IF

! For NM !
DO i = 1, KDIV
	DO j = 1, NDIV
		fhat2(i,j) = -phi(i,j) - h02new2*psi2(i,j)
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE POTENTIAL_FINAL
USE DEFINITION 
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, m, l

! For DM !
IF(dm_flag == 1) THEN

	! Initialization !
	d1 = 0.0D0
	d2 = 0.0D0
	phi1 = 0.0D0

	! Find d1 !
	DO m = 0, lmax
		DO k = 1, NDIV
			DO i = 1, KDIV-2, 2
				d1(k,m) = d1(k,m) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))*(pn(i,2*m)*rho1(i,k) + &
				4.0D0*pn(i+1,2*m)*rho1(i+1,k) + pn(i+2,2*m)*rho1(i+2,k))
			END DO
		END DO
	END DO
	
	! Find d2 !
	DO m = 0, lmax
		DO j = 1, NDIV
			DO k = 1, NDIV-2, 2
				d2(m,j) = d2(m,j) + (1.0D0/6.0D0)*(r(k+2) - r(k))*(fp(k,j,2*m)*d1(k,m) + &
				4.0D0*fp(k+1,j,2*m)*d1(k+1,m) + fp(k+2,j,2*m)*d1(k+2,m))
			END DO
		END DO	
	END DO

	! Find DM potential !
	DO i = 1, KDIV
		DO j = 1, NDIV
			DO m = 0, LMAX
				phi1(i,j) = phi1(i,j) + (-4.0D0*pi)*d2(m,j)*pn(i,2*m)
			END DO
		END DO
	END DO

END IF

! Initialization !
d1 = 0.0D0
d2 = 0.0D0
phi2 = 0.0D0

! Find d1 !
DO m = 0, lmax
	DO k = 1, NDIV
		DO i = 1, KDIV-2, 2
			d1(k,m) = d1(k,m) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))*(pn(i,2*m)*rho2(i,k) + &
			4.0D0*pn(i+1,2*m)*rho2(i+1,k) + pn(i+2,2*m)*rho2(i+2,k))
		END DO
	END DO
END DO
	
! Find d2 !
DO m = 0, lmax
	DO j = 1, NDIV
		DO k = 1, NDIV-2, 2
			d2(m,j) = d2(m,j) + (1.0D0/6.0D0)*(r(k+2) - r(k))*(fp(k,j,2*m)*d1(k,m) + &
			4.0D0*fp(k+1,j,2*m)*d1(k+1,m) + fp(k+2,j,2*m)*d1(k+2,m))
		END DO
	END DO	
END DO

! Find NM potential !
DO i = 1, KDIV
	DO j = 1, NDIV
		DO m = 0, LMAX
			phi2(i,j) = phi2(i,j) + (-4.0D0*pi)*d2(m,j)*pn(i,2*m)
		END DO
	END DO
END DO

END SUBROUTINE