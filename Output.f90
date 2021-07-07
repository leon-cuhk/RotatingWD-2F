SUBROUTINE OUTPUTLOG
USE DEFINITION
IMPLICIT NONE

! Output for NM !
WRITE (11, *) '---------------------------------------------'
WRITE (11, *) '|Print Out Result (in cgs, except for Mass) |'
WRITE (11, *) '---------------------------------------------'
WRITE (11, *) 'Axis-ratio', axratio2
WRITE (11, *) 'Log NM Central Density', rhoscale2
WRITE (11, *) 'Equatorial radius', r_equatorial
WRITE (11, *) 'NM mass', (mass2*mass/solar)
WRITE (11, *) 'NM volume', (volume2*vol)
WRITE (11, *) 'NM Angular momentum', (j2*momentum)
WRITE (11, *) 'NM Kinetic energy', (kine2*energy)
WRITE (11, *) 'NM Gravitational energy', (gravw2*energy)
WRITE (11, *) 'NM 3 times pressure integral', (3.0D0*pint2*energy)
WRITE (11, *) 'Pressure NM maximum', (maxval(p2)*pressure)
WRITE (11, *) 'Viral test', (ABS(2.0D0*(kine1 + kine2) + gravw + 3.0D0*(pint1 + pint2))/ABS(gravw))
WRITE (11, *) '------------'
WRITE (11, *) '|End Report|'
WRITE (11, *) '------------'

! Output for DM !
IF(dm_flag == 1) THEN
	WRITE (12, *) '---------------------------------------------'
	WRITE (12, *) '|Print Out Result (in cgs, except for Mass) |'
	WRITE (12, *) '---------------------------------------------'
	WRITE (12, *) 'Size-ratio', sizeratio
	WRITE (12, *) 'Log DM Central Density', rhoscale1
	WRITE (12, *) 'Equatorial radius', r_equatorial
	WRITE (12, *) 'DM mass', (mass1*mass/solar)
	WRITE (12, *) 'Absolute DM mass error', abserror 
	WRITE (12, *) 'DM volume', (volume1*vol)
	WRITE (12, *) 'DM Angular momentum', (j1*momentum)
	WRITE (12, *) 'DM Kinetic energy', (kine1*energy)
	WRITE (12, *) 'DM Gravitational energy', (gravw1*energy)
	WRITE (12, *) 'DM 3 times pressure integral', (3.0D0*pint1*energy)
	WRITE (12, *) 'Pressure DM maximum', (maxval(p1)*pressure)
	WRITE (12, *) 'Viral test', ABS(2.0D0*(kine1 + kine2) + gravw + 3.0D0*(pint1 + pint2))/ABS(gravw)
	WRITE (12, *) '------------'
	WRITE (12, *) '|End Report|'
	WRITE (12, *) '------------'
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUTPROFILE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j

! For NM, Density  profile
WRITE (41, *) KDIV, NDIV, rmax
DO j = 1, NDIV
	WRITE (41, 701) (rho2(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
END DO
WRITE (41, *)

! Entalphy profile
WRITE (42, *) KDIV, NDIV, rmax
DO j = 1, NDIV
	WRITE (42, 701) (hhatnew2(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
END DO
WRITE (42, *)

! Potential profile
WRITE (43, *) KDIV, NDIV, rmax
DO j = 1, NDIV
	WRITE (43, 701) (phi2(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
END DO
WRITE (43, *)

! For DM !
IF(dm_flag == 1) THEN

	! For DM, Density  profile
	WRITE (51, *) KDIV, NDIV, rmax
	DO j = 1, NDIV
		WRITE (51, 701) (rho1(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (51, *)

	! Entalphy profile
	WRITE (52, *) KDIV, NDIV, rmax
	DO j = 1, NDIV
		WRITE (52, 701) (hhatnew1(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (52, *)

	! Potential profile
	WRITE (53, *) KDIV, NDIV, rmax
	DO j = 1, NDIV
		WRITE (53, 701) (phi1(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (53, *)

	! For total, density profile
	WRITE (61, *) KDIV, NDIV, rmax
	DO j = 1, NDIV
		WRITE (61, 701) (rho(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (61, *)

	! Potential profile
	WRITE (62, *) KDIV, NDIV, rmax
	DO j = 1, NDIV
		WRITE (62, 701) (phi(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (62, *)

END IF

! Format !
701 FORMAT (1200ES16.8)

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUTGLOBAL
USE DEFINITION
IMPLICIT NONE

! For relation curves !
WRITE (21, 701) masst, m2out, reout, j2out, h0out, wmax2, vmax2, stab, bind, axratio2, ibar, qbar, rbar, log10(rhomax2), ktidal, viral

! For DM !
IF(dm_flag == 1) THEN
	WRITE (22, 701) masst, m1out, rhoscale1, re1, stab, axratio1
END IF

! Format !
701 FORMAT (1200ES16.8)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUTDMNM
USE DEFINITION
IMPLICIT NONE

! For relation curves !
WRITE (71, 701) m1out, m2out, masst, abserror, stab

! Format !
701 FORMAT (1200ES16.8)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OUTPUTFINAL
USE DEFINITION
IMPLICIT NONE

! Open !
OPEN (UNIT = 999, FILE = './Profile/Star_WENO_Parameter.dat', STATUS = 'REPLACE')

! Write !
WRITE (999,*) 'DMFlag', dm_flag
WRITE (999,*) 'DMMass', mass_dm
IF(rigid2 == 1) THEN
	WRITE (999,*) 'Rotation', ' rigid'
ELSEIF(vconst2 == 1) THEN
	WRITE (999,*) 'Rotation', ' vconst'
ELSEIF(jconst2 == 1) THEN
	WRITE (999,*) 'Rotation', ' jconst'
ELSEIF(kepler2 == 1) THEN
	WRITE (999,*) 'Rotation', ' kepler'
ELSEIF(yoon2 == 1) THEN
	WRITE (999,*) 'Rotation', ' yoon'
ELSEIF(awd2 == 1) THEN
	WRITE (999,*) 'Rotation', ' awd'
END IF
WRITE (999,*) 'rmax', rmax
WRITE (999,*) 'KDIV', KDIV
WRITE (999,*) 'NDIV', NDIV
WRITE (999,*) 'nrho', n_rho
WRITE (999,*) 'rhostart', rhostart
WRITE (999,*) 'rhoend', rhoend
WRITE (999,*) 'drho', drho
WRITE (999,*) 'n_axis', n_axis
WRITE (999,*) 'axstart', axstart
WRITE (999,*) 'axend', axend
WRITE (999,*) 'daxratio', daxratio

! Close !
CLOSE (999)

END SUBROUTINE