SUBROUTINE EOSTABLE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Dummy variables !
REAL (DP) :: dummy

! Read the number of lines in the file !
nlines2 = 0 
OPEN (999, file = './Table/EOSTABLE.eos') 
DO 
    READ (999,*, END=10) 
    nlines2 = nlines2 + 1 
END DO 
10 CLOSE (999) 

! Allocate arrays !
ALLOCATE(rhotable2(0:nlines2-1))
ALLOCATE(htable2(0:nlines2-1))
ALLOCATE(xtable2(0:nlines2-1))

! Read !
OPEN(UNIT=999, FILE = './Table/EOSTABLE.eos', ACTION='READ')
DO i = 1, nlines2 
	IF(i == 1) THEN
		READ(999,*) ye2, dummy
	ELSE
		READ(999,*) rhotable2(i-2), xtable2(i-2), htable2(i-2)
	END IF
ENDDO
CLOSE(999)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE YETABLE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i

! Dummy variables !
REAL (DP) :: dummy

! Read the number of lines in the file !
yelines2 = 0 
OPEN (999, file = './Table/y_e_vs_rho_at_7ms_Tc5.0d09_from_VULCAN2D.dat') 
DO 
    READ (999,*, END=10) 
    yelines2 = yelines2 + 1 
END DO 
10 CLOSE (999) 

! Allocate arrays !
ALLOCATE(yetable2(0:yelines2-1))
ALLOCATE(rhoyetable2(0:yelines2-1))

! Read !
OPEN(UNIT=999, FILE = './Table/y_e_vs_rho_at_7ms_Tc5.0d09_from_VULCAN2D.dat', ACTION='READ')
DO i = 1, yelines2
	READ(999,*) rhoyetable2(yelines2-i), yetable2(yelines2-i)
	rhoyetable2(yelines2-i) = log10(rhoyetable2(yelines2-i))
ENDDO
CLOSE(999)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSRtoYE(rho_in, ye_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: rho_in
REAL (DP), INTENT(OUT) :: ye_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: rho_target

! Take log !
rho_target = log10(rho_in)

! Case by case !
IF(rho_target == rhoyetable2(0)) THEN

	! Table minimum !
	ye_out = yetable2(0)

ELSE

	! Binary search !
	left = 0
	right = yelines2 - 1
	DO
		IF(left > right) THEN
			IF(rhoyetable2(m) > rho_target) THEN
				m = m - 1
			ELSEIF(rhoyetable2(m) < rho_target) THEN
				m = m
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(rhoyetable2(m) < rho_target) THEN
			left = m + 1	
		ELSEIF(rhoyetable2(m) > rho_target) THEN
			right = m - 1
		ELSEIF(rhoyetable2(m) == rho_target) THEN
			EXIT
		END IF
	END DO
	IF(rhoyetable2(m) == rho_target) THEN
		ye_out = yetable2(m)
	ELSE
		!IF(m > 1 .AND. m < nlines2 - 3) THEN
		!	CALL AKIMA(rhoyetable2(m-2), rhoyetable2(m-1), rhoyetable2(m), rhoyetable2(m+1), rhoyetable2(m+2), rhoyetable2(m+3), & 
		!		yetable2(m-2), yetable2(m-1), yetable2(m), yetable2(m+1), yetable2(m+2), yetable2(m+3), rho_target, ye_out)
		!ELSE
			! Linear !
			CALL LINEAR(rhoyetable2(m), rhoyetable2(m+1), yetable2(m), yetable2(m+1), rho_target, ye_out)
		!END IF
	END IF

END IF

! Take power !
ye_out = ye_out

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSRtoH(rho_in, h_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: rho_in
REAL (DP), INTENT(OUT) :: h_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: rho_target

! Take log !
rho_target = log10(rho_in)

! Case by case !
IF(rho_target == rhotable2(0)) THEN

	! Table minimum !
	h_out = htable2(0)

ELSE

	! Binary search !
	left = 0
	right = nlines2 - 1
	DO
		IF(left > right) THEN
			IF(rhotable2(m) > rho_target) THEN
				m = m - 1
			ELSEIF(rhotable2(m) < rho_target) THEN
				m = m
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(rhotable2(m) < rho_target) THEN
			left = m + 1	
		ELSEIF(rhotable2(m) > rho_target) THEN
			right = m - 1
		ELSEIF(rhotable2(m) == rho_target) THEN
			EXIT
		END IF
	END DO
	IF(rhotable2(m) == rho_target) THEN
		h_out = htable2(m)
	ELSE
		!IF(m > 1 .AND. m < nlines2 - 3) THEN
		!	CALL AKIMA(rhotable2(m-2), rhotable2(m-1), rhotable2(m), rhotable2(m+1), rhotable2(m+2), rhotable2(m+3), & 
		!		htable2(m-2), htable2(m-1), htable2(m), htable2(m+1), htable2(m+2), htable2(m+3), rho_target, h_out)
		!ELSE
			! Linear !
			CALL LINEAR(rhotable2(m), rhotable2(m+1), htable2(m), htable2(m+1), rho_target, h_out)
		!END IF
	END IF

END IF

! Take power !
h_out = 10.0D0**(h_out)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSHtoR(h_in, rho_out)
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: h_in
REAL (DP), INTENT(OUT) :: rho_out

! Integer !
INTEGER :: left, right, m

! Target density !
REAL (DP) :: htarget

! Convert to CGS and take log !
htarget = log10(h_in*potential)

! Case by case !
IF(htarget == htable2(0)) THEN

	! Table minimum !
	rho_out = rhotable2(0)

ELSE

	! Binary search !
	left = 0
	right = nlines2 - 1
	DO
		IF(left > right) THEN
			IF(htable2(m) > htarget) THEN
				m = m - 1
			ELSEIF(htable2(m) < htarget) THEN
				m = m
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(htable2(m) < htarget) THEN
			left = m + 1	
		ELSEIF(htable2(m) > htarget) THEN
			right = m - 1
		ELSEIF(htable2(m) == htarget) THEN
			EXIT
		END IF
	END DO

	IF(htable2(m) == htarget) THEN
		rho_out = rhotable2(m)
	ELSE
		!IF(m > 1 .AND. m < nlines2 - 3 .AND. vul_flag /= 1) THEN
		!	CALL AKIMA(htable2(m-2), htable2(m-1), htable2(m), htable2(m+1), htable2(m+2), htable2(m+3), & 
		!		rhotable2(m-2), rhotable2(m-1), rhotable2(m), rhotable2(m+1), rhotable2(m+2), rhotable2(m+3), htarget, rho_out)
		!ELSE
			! Linear !
			CALL LINEAR(htable2(m), htable2(m+1), rhotable2(m), rhotable2(m+1), htarget, rho_out)
		!END IF
	END IF
	
END IF

! Convert to code unit !
rho_out = (10.0D0**(rho_out))!/(rhomax2)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EOSHtoX(h_in, x_out)
USE DEFINITION
IMPLICIT NONE

! Input and output !
REAL (DP), INTENT(IN) :: h_in
REAL (DP), INTENT(OUT) :: x_out

! Integer !
INTEGER :: left, right, m, n

! Target density !
REAL (DP) :: htarget

! Convert to number density !
htarget = log10(h_in*potential)

! Case by case !
IF(htarget == htable2(0)) THEN

	! Table minimum !
	x_out = xtable2(0)

ELSE

	! Binary search !
	left = 0
	right = nlines2 - 1
	DO
		IF(left > right) THEN
			IF(htable2(m) > htarget) THEN
				m = m - 1
			ELSEIF(htable2(m) < htarget) THEN
				m = m
			END IF
			EXIT
		END IF
		m = floor(REAL((right + left)/2))
		IF(htable2(m) < htarget) THEN
			left = m + 1	
		ELSEIF(htable2(m) > htarget) THEN
			right = m - 1
		ELSEIF(htable2(m) == htarget) THEN
			EXIT
		END IF
	END DO
	IF(htable2(m) == htarget) THEN
		x_out = xtable2(m)
	ELSE
		!IF(m > 1 .AND. m < nlines2 - 3) THEN
		!	CALL AKIMA(htable2(m-2), htable2(m-1), htable2(m), htable2(m+1), htable2(m+2), htable2(m+3), & 
		!		xtable2(m-2), xtable2(m-1), xtable2(m), xtable2(m+1), xtable2(m+2), xtable2(m+3), htarget, x_out)
		!ELSE
			! Linear !
			CALL LINEAR(htable2(m), htable2(m+1), xtable2(m), xtable2(m+1), htarget, x_out)
		!END IF
	END IF	

END IF

! Convert to cgs !
x_out = x_out

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing !
!datum and output the quantity that the user wished to interpolate 	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINEAR(x0, x1, y0, y1, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, y0, y1, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1
REAL (DP) :: de0, de1
REAL (DP) :: l0, l1

! Assign numerator !
nu0 = (x_in - x1)
nu1 = (x_in - x0)

! Assign denominator !
de0 = (x0 - x1)
de1 = (x1 - x0)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1

! Compute the output !
y_out = l0*y0 + l1*y1

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Akima spline interpolation. See Hiroshi Akima 1970 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AKIMA(xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: dm2, dm1, d0, dp1, dp2
REAL (DP) :: s0, s1

! Weights !
REAL (DP) :: w1, w2, w3, w4

! Coefficient of polynominal !
REAL (DP) :: p0, p1, p2, p3

! Temporal arrays !
REAL (DP) :: diff, temp

! Assign slopes !
dm2 = (ym1 - ym2)/(xm1 - xm2)
dm1 = (y0 - ym1)/(x0 - xm1)
d0 = (yp1 - y0)/(xp1 - x0)
dp1 = (yp2 - yp1)/(xp2 - xp1)
dp2 = (yp3 - yp2)/(xp3 - xp2)

! Assign weights (modified) !
w1 = abs(dp1 - d0) + 0.5D0*abs(dp1 + d0)
w2 = abs(dm1 - dm2) + 0.5D0*abs(dm1 + dm2)
w3 = abs(dp2 - dp1) + 0.5D0*abs(dp2 + dp1)
w4 = abs(d0 - dm1) + 0.5D0*abs(d0 + dm1)

! assign slopes !
IF(w1 == 0.0D0 .AND. w2 == 0.0D0) THEN
	s0 = 0.5D0*(dm1 + d0)
ELSE
	s0 = (w1*dm1 + w2*d0)/(w1 + w2)
END IF
IF(w3 == 0.0D0 .AND. w4 == 0.0D0) THEN
	s1 = 0.5D0*(d0 + dm1)
ELSE
	s1 = (w3*d0 + w4*dp1)/(w3 + w4)
END IF

! Assign temp !
diff = xp1 - x0
temp = x_in - x0

! assign coefficients !
p0 = y0
p1 = s0
p2 = (3.0D0*d0 - 2.0D0*s0 - s1)/diff
p3 = (s0 + s1 - 2.0D0*d0)/diff**2

! Output the interpolation !
y_out = p0 + p1*temp + p2*temp**2 + p3*temp**3

END SUBROUTINE