!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program generates rotating star according to the method of !
! self-consistent iterations by Hachisu 1986 apj, 61:479-507	  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM ROTATING
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE (*,*) '-----------------------'
WRITE (*,*) '|Rotating Star Program|'
WRITE (*,*) '-----------------------'
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read EOS table !
IF(fermi_flag == 1) THEN
	ye2 = 5.0E-1_DP
ELSE
	CALL EOSTABLE
	IF(yerho_flag == 1) THEN
		IF(vul_flag == 1) THEN
			CALL YETABLE
		ELSE
			CALL GETCONST
		END IF
	END IF
END IF

! Initialize aray !
CALL INITIAL
CALL FINDCONST
CALL ATMENTHALPY
CALL FINDLEGENDRE
CALL FINDEXPANSION

! Find rotational potential !
CALL FINDROTATION

! Get NM boundary !
CALL NMAXIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Start computing model !
IF(dm_flag == 1) THEN
	WRITE (*,*) 'Start Solving For Two-Fluid Star'
	WRITE (*,*)
	CALL SCF2F
ELSE
	WRITE (*,*) 'Start Solving For One-Fluid Star'
	WRITE (*,*)
	CALL SCF1F
END IF

! Output final !
CALL OUTPUTFINAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print out !
WRITE (*,*) 'Done!'
 
END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKMODEL
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Check equatorial radius !
DO j = 1, NDIV
	IF(rho2(1,j) > rhoa2 .AND. rho2(1,j+1) == rhoa2) THEN
		r_equator2 = r(j+1)
		EXIT 
	END IF
END DO

! Check axis radius !
DO j = 1, NDIV
	IF(rho2(KDIV,j) > rhoa2 .AND. rho2(KDIV,j+1) == rhoa2) THEN
		r_axis2 = r(j+1)
		EXIT 
	END IF
END DO

! For DM !
DO j = 1, NDIV
	IF(rho1(1,j) > rhoa1 .AND. rho1(1,j+1) == rhoa1) THEN
		r_equator1 = r(j+1)
		EXIT 
	END IF
END DO
DO j = 1, NDIV
	IF(rho1(KDIV,j) > rhoa1 .AND. rho1(KDIV,j+1) == rhoa1) THEN
		r_axis1 = r(j+1)
		EXIT 
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WRITE (*,*) 'DM Radius'
!WRITE (*,*) r_equator1, r(ra1)
!WRITE (*,*) 'NM Radius'
!WRITE (*,*) r_equator2, r(ra2)
!WRITE (*,*) r_axis2, r(rb2)
!WRITE (*,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Determine whether the model is successful !
IF(dm_flag == 1) THEN
	IF(rb2 == 1) THEN
		IF(r_equator1 == r(ra1) .AND. r_equator2 == r(ra2)) THEN
			success_flag = .true.
		ELSE
			success_flag = .false.
		END IF
	ELSE	
		IF(r_equator1 == r(ra1) .AND. r_equator2 == r(ra2) .AND. r_axis2 == r(rb2)) THEN
			success_flag = .true.
		ELSE
			success_flag = .false.
		END IF
	END IF
ELSE
	IF(rb2 == 1) THEN
		IF(r_equator2 == r(ra2)) THEN
			success_flag = .true.
		ELSE
			success_flag = .false.
		END IF
	ELSE	
		IF(r_equator2 == r(ra2) .AND. r_axis2 == r(rb2)) THEN
			success_flag = .true.
		ELSE
			success_flag = .false.
		END IF
	END IF
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECKENTHALPY
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Modify density !
DO i = 1, KDIV
	DO j = 1, NDIV
		IF(r(j) > rout) THEN
			hhatnew2(i,j) = 0.0D0 !hhata2
		END IF
	END DO
END DO

! For DM !
IF(DM_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(r(j) > rout*r(ra1)) THEN
				hhatnew1(i,j) = 0.0D0 !hhata1
			END IF
		END DO
	END DO
END IF

END SUBROUTINE