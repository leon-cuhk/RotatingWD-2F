SUBROUTINE ANALYSIS(sum_in, sum_out)
USE DEFINITION
IMPLICIT NONE

! Integer !
Integer :: i, j, k

! Output !
REAL (DP) :: sum_out

! Input !
REAL (DP), INTENT(IN), DIMENSION(KDIV,NDIV) :: sum_in

! Local arrays !
REAL (DP), DIMENSION(NDIV) :: Q_in

! Initialize !
sum_out = 0.0D0
Q_in = 0.0D0

DO j = 1, NDIV
	DO i = 1, KDIV-2, 2
		Q_in (j) = Q_in(j) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))* & 
		(sum_in(i,j) + 4.0D0*sum_in(i+1,j) + sum_in(i+2,j))
	END DO
END DO

DO j = 1, NDIV-2, 2
	sum_out = sum_out + (1.0D0/6.0D0)*(r(j+2) - r(j))* & 
	(r(j)**2*Q_in (j) + 4.0D0*r(j+1)**2*Q_in (j+1) + r(j+2)**2*Q_in (j+2))
END DO

sum_out = sum_out*4.0D0*pi

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDVOLUME
USE DEFINITION
IMPLICIT NONE

! Integer !
Integer :: i, j, k

! Local arrays !
REAL (DP), DIMENSION(NDIV) :: Q_in

! Local arrays !
REAL (DP), DIMENSION(KDIV,NDIV) :: factor

! Initialize !
volume2 = 0.0D0
Q_in = 0.0D0
factor = 0.0D0

! Do for NM !
DO i = 1, KDIV
	DO j = 1, NDIV
		If(rho2(i,j) > 0.0D0) THEN
			factor(i,j) = 1.0D0
		END IF
	END DO
END DO

DO j = 1, NDIV
	DO i = 1, KDIV-2, 2
		Q_in (j) = Q_in(j) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))* & 
		(factor(i,j) + 4.0D0*factor(i+1,j) + factor(i+2,j))
	END DO
END DO

DO j = 1, NDIV-2, 2
	volume2 = volume2 + (1.0D0/6.0D0)*(r(j+2) - r(j))* & 
	(r(j)**2*Q_in (j) + 4.0D0*r(j+1)**2*Q_in (j+1) + r(j+2)**2*Q_in (j+2))
END DO

volume2 = volume2*4.0D0*pi

IF(dm_flag == 1) THEN

	! Initialize !
	volume1 = 0.0D0
	Q_in = 0.0D0
	factor = 0.0D0

	! Do for DM !
	DO i = 1, KDIV
		DO j = 1, NDIV
			If(rho1(i,j) > 0.0D0) THEN
				factor(i,j) = 1.0D0	
			END IF
		END DO
	END DO

	DO j = 1, NDIV
		DO i = 1, KDIV-2, 2
			Q_in (j) = Q_in(j) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))* & 
			(factor(i,j) + 4.0D0*factor(i+1,j) + factor(i+2,j))
		END DO
	END DO

	DO j = 1, NDIV-2, 2
		volume1 = volume1 + (1.0D0/6.0D0)*(r(j+2) - r(j))* & 
		(r(j)**2*Q_in (j) + 4.0D0*r(j+1)**2*Q_in (j+1) + r(j+2)**2*Q_in (j+2))
	END DO

	volume1 = volume1*4.0D0*pi


	! Initialize !
	vtotal = 0.0D0
	Q_in = 0.0D0
	factor = 0.0D0

	! Do for total volume !
	DO i = 1, KDIV
		DO j = 1, NDIV
			If(rho1(i,j) > 0.0D0 .OR. rho2(i,j) > 0.0D0) THEN
				factor(i,j) = 1.0D0	
			END IF
		END DO
	END DO

	DO j = 1, NDIV
		DO i = 1, KDIV-2, 2
			Q_in (j) = Q_in(j) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))* & 
			(factor(i,j) + 4.0D0*factor(i+1,j) + factor(i+2,j))
		END DO
	END DO

	DO j = 1, NDIV-2, 2
		vtotal = vtotal + (1.0D0/6.0D0)*(r(j+2) - r(j))* & 
		(r(j)**2*Q_in (j) + 4.0D0*r(j+1)**2*Q_in (j+1) + r(j+2)**2*Q_in (j+2))
	END DO

	vtotal = vtotal*4.0D0*pi

ELSE
	
	vtotal = volume2

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDQUANTITY_1
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! For DM !
IF(dm_flag == 1) THEN

	! Initialize !
	ke1 = 0.0D0

END IF

! Initialize !
ke2 = 0.0D0

! Find NM quantity !
DO i = 1, KDIV
	DO j = 1, NDIV
		If(rho2(i,j) > 0.0D0) THEN	
			ke2(i,j) = 0.5D0*rho2(i,j)*velp2(i,j)**2
		END IF
	END DO
END DO

! Initialize !
grav = 0.0D0

! Find total gravitational energy density !
DO i = 1, KDIV
	DO j = 1, NDIV
		If(rho(i,j) > 0.0D0) THEN	
			grav(i,j) = rho(i,j)*0.5D0*phi(i,j)
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDQUANTITY_2
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

IF(dm_flag == 1) THEN

	! Initialize !
	i1 = 0.0D0
	q1 = 0.0D0
	grav1 = 0.0D0
	pphi1 = 0.0D0

	! Find DM quantity !
	DO i = 1, KDIV
		DO j = 1, NDIV
			If(rho1(i,j) > 0.0D0) THEN
				i1(i,j) = rho1(i,j)*r_polar(i,j)**2
				q1(i,j) = rho1(i,j)*r(j)**2*pn(i,2)
				grav1(i,j) = rho1(i,j)*(0.5D0*phi1(i,j) + phi2(i,j))
			END IF
		END DO
	END DO

END IF

! Initialize !
i2 = 0.0D0
q2 = 0.0D0
grav2 = 0.0D0
pphi2 = 0.0D0

! Find NM quantity !
DO i = 1, KDIV
	DO j = 1, NDIV
		If(rho2(i,j) > 0.0D0) THEN	
			i2(i,j) = rho2(i,j)*r_polar(i,j)**2
			q2(i,j) = rho2(i,j)*r(j)**2*pn(i,2)
			grav2(i,j) = rho2(i,j)*(0.5D0*phi2(i,j) + phi1(i,j))
			pphi2(i,j) = rho2(i,j)*velp2(i,j)*r_polar(i,j)
		END IF
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ANGULAR(array_in, array_out)
USE DEFINITION
IMPLICIT NONE

! Integer !
Integer :: i, j, k

! Input !
REAL (DP), INTENT(IN), DIMENSION(KDIV,NDIV) :: array_in

! Output !
REAL (DP), INTENT(OUT), DIMENSION(NDIV) :: array_out

! Initialize !
array_out = 0.0D0

! Simpson's Rule !
DO j = 1, NDIV
	DO i = 1, KDIV - 2, 2
		array_out (j) = array_out (j) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))* & 
			(array_in(i,j) + 4.0D0*array_in(i+1,j) + array_in(i+2,j))
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDRADIUS
USE DEFINITION
IMPLICIT NONE

! Equatorial radius !
r_equatorial = SQRT(hmax2/hhmaxnew2/gconst/rhomax2)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDSCALE
USE DEFINITION
IMPLICIT NONE

! Assign !
vol = r_equatorial**3
length = r_equatorial
mass = (r_equatorial**3*rhomax2)
persecond = SQRT(gconst*rhomax2)
energy = (gconst*r_equatorial**5*rhomax2**2)
potential = (gconst*r_equatorial**2*rhomax2)
clighthat = (clight/(r_equatorial*persecond))
pressure = (gconst*r_equatorial**2*rhomax2**2)
momentum = SQRT(gconst*r_equatorial**10*rhomax2**3)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDGLOBAL
USE DEFINITION
IMPLICIT NONE

! REAL !
REAL (DP) :: dummy

! DM axis-ratio !
IF(DM_flag == 1) THEN
	axratio1 = (r_axis1/r_equator1)
END IF

! Dimensionless variables !
IF(DM_flag == 1) THEN
	ibar = itotal1 + itotal2
	qbar = qtotal1 + qtotal2
	masst = mass1 + mass2
ELSE
	ibar = itotal2
	qbar = qtotal2
	masst = mass2
END IF
rbar = 1.0D0
ibar = log(ibar*clighthat**4/masst**3)
qbar = log(abs(qbar*masst*clighthat**2/j2**2))

! Stability parameter !
stab = (kine1 + kine2)/ABS(gravw)

! Compute quantity !
IF(rigid2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2)
ELSEIF(vconst2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2*r_equatorial)
ELSEIF(jconst2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2)*(r_equatorial**2)
ELSEIF(kepler2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2)*(r_equatorial**(1.5D0))
ELSEIF(yoon2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2)
ELSEIF(awd2 == 1) THEN
	h0out = SQRT(h02new2)*SQRT(gconst*rhomax2)
END IF

! Compute masses !
m1out = mass1*(mass/solar)
m2out = mass2*(mass/solar)
masst = m1out + m2out

! Radius !
reout = (r_equatorial/rscale)
re1 = r(ra1)*reout

! Find maximum DM density !
rhoscale1 = log10(maxval(rho1)*rhomax2)	

! Angular momentum !
j2out = j2*(momentum/jscale)

! Equatorial velocity !
vmax2 = maxval(velp2)*SQRT(gconst*rhomax2*r_equatorial)
wmax2 = maxval(omega2)*SQRT(gconst*rhomax2)

! Binding energies !
bind = (kine1 + kine2 + int1 + int2 + gravw)*(energy/escale)

! viral test !
viral = (ABS(2.0D0*(kine1 + kine2) + gravw + 3.0D0*(pint1 + pint2))/ABS(gravw))

! Compactness !
IF(DM_flag == 1) THEN
	dummy = (mass1 + mass2)/(max(r(ra1), r(ra2))*clighthat**2)
ELSE
	dummy = (mass2)/(r(ra2)*clighthat**2)
END IF

! Tidal love number !
IF(axratio2 == 1.0D0) THEN
	ktidal = log(ktidal*(2.0D0/3.0D0)/dummy**5)
ELSE
	ktidal = 0.0D0
END IF

END SUBROUTINE