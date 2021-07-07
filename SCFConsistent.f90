SUBROUTINE SCF1F
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, l, m, n, o

! For iterations !
REAL (DP) :: chatold2, chatnew2, fhmax2

! Potentail at boundaries point !
REAL (DP) :: psia2, psib2, phia2, phib2
REAL (DP) :: fa2, fmax2

! Exit criteria !
REAL (DP) :: criteria1, criteria2, criteria3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the loop in generating a lot of models !
DO l = 1, n_rho
	
	! Assign NM maximum density !
	rhoscale2 = rhostart + (DBLE(l) - 1.0D0)*drho
	rhomax2 = 1.0D1**(rhoscale2)

	! Maximum Enthalpy !
	CALL MAXENTHALPY

	! Openfile !
	IF (output_results .eqv. .true.) THEN
		CALL OPENFILE_LOG
	END IF
	CALL OPENFILEGLOBAL

	! Do the loop in axis ratio !
	DO m = axstart, axend, daxratio

		! Get NM boundary !
		rb2 = m

		! Exit if the boundary is out of range !
		IF(rb2 < 1) THEN
			EXIT
		END IF

		! Find axis ratio !
		axratio2 = (r(rb2)/r(ra2))

		! Print out !
		WRITE (*,*) 'Start iterations for ...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Assign uniform density at the beginning of iteration 
		CALL INITIALRHO

		! Add the density !
		CALL GRAVRHO
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find Potential !
		CALL FINDPOTENTIAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Potential ar boundary point !
		phia2 = phi(1,ra2)
		phib2 = phi(KDIV,rb2)
		psia2 = psi2(1,ra2)
		psib2 = psi2(KDIV,rb2)

		! Find h02 for NM !
		IF(axratio2 == 1.0D0) THEN
			h02new2 = 0.0D0
		ELSE
			h02new2 = -(phia2 - phib2)/(psia2 - psib2)
		END IF

		! Find modified potential !
		CALL FINDFPOTENTIAL

		! Assign f-potential !
		fmax2 = maxval(fhat2)
		fa2 = fhat2(1,ra2)
		
		! Find chat !
		chatnew2 = (hconst2*fmax2 - fa2)/(1.0D0 - hconst2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Initialize !
		success_flag = .true.

		! Start the loop !
		DO n = 1, nmax

			! First backup !
			h02old2 = h02new2
			chatold2 = chatnew2
			hhmaxold2 = hhmaxnew2

			! Find enthalpy !
			DO i = 1, KDIV
				DO j = 1, NDIV
					hhatnew2(i,j) = chatnew2 - phi(i,j) - h02new2*psi2(i,j)
				END DO
			END DO
			
			! Assign hhmax !
			hhmaxnew2 = maxval(hhatnew2)
			
			! Find equatorial radius !
			CALL FINDRADIUS

			! Find all the scaling constant !
			CALL FINDSCALE

			! Atmosphere !
			CALL ATMOSPHERE
			
			! Check the enthalpy !
			CALL CHECKENTHALPY

			! Get new density !
			CALL FINDDENSITY
	
			! Override enthalpy !
			hhatnew2(1,ra2) = hhata2
			hhatnew2(KDIV,rb2) = hhata2
			
			! Override density !
			rho2(1,ra2) = rhoa2
			rho2(KDIV,rb2) = rhoa2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Override total density !
			CALL GRAVRHO

			! Find gravitational potential !
			CALL FINDPOTENTIAL

			! Potential ar boundary point !
			phia2 = phi(1,ra2)
			phib2 = phi(KDIV,rb2)
			
			! Find h02 !
			IF(axratio2 == 1.0D0) THEN
				h02new2 = 0.0D0
			ELSE
				h02new2 = -(phia2 - phib2)/(psia2 - psib2)
			END IF
			
			! Find modified potential !
			CALL FINDFPOTENTIAL
			
			! Assign f-potential !
			fmax2 = maxval(fhat2)
			fa2 = fhat2(1,ra2)

			! Find chat !
			chatnew2 = (hconst2*fmax2 - fa2)/(1.0D0 - hconst2)

			! Exit conditions !
			IF(n > 1) THEN
				criteria1 = ABS((h02new2 - h02old2)/h02old2)
				criteria2 = ABS((chatnew2 - chatold2)/chatold2)
				criteria3 = ABS((hhmaxnew2 - hhmaxold2)/hhmaxold2)
				IF(axratio2 == 1.0D0) THEN
					IF(criteria2 < tor .AND. criteria3 < tor) EXIT
				ELSE
					IF(criteria1 < tor .AND. criteria2 < tor .AND. criteria3 < tor) EXIT
				END IF
			END IF

			! Not sucessful if the model has NaN !
			IF(ieee_is_nan(h02new2) .OR. ieee_is_nan(chatnew2) .OR. ieee_is_nan(hhmaxnew2)) THEN
				WRITE (*,*) 'There is NaN'
				success_flag = .false.
				EXIT 
			END IF

			! Not sucessful if the model failed to converge !
			If(n == nmax) THEN
				WRITE (*,*) 'Maximum iteration reached'
				success_flag = .false.
				EXIT
			END IF

		END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Be sure to check density !
		IF(success_flag .eqv. .true.) THEN
			CALL CHECKMODEL
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Mass shedding occurs'
			END IF
		END IF

		! Exit if failed to converge !
		IF(critical_flag .eqv. .true.) THEN
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Critical rotation encountered'
				WRITE (*,*) 'Move on to the next density'
				WRITE (*,*) 
				EXIT
			END IF
		END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find rotational velocity and momentum !
		CALL FINDVELOCITY

		! Find physical quantity !
		CALL FINDQUANTITY_1

		! Find total kinetic energy !
		CALL ANALYSIS(ke2, kine2)

		! Find total gravitational energy !
		CALL ANALYSIS(grav, gravw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find the individual potential !
		CALL POTENTIAL_FINAL

		! Find volume !
		CALL FINDVOLUME	
		
		! Find pressure !
		CALL FINDPRESSURE

		! Find energy !
		CALL FINDEPSILON

		! Find physical quantity !
		CALL FINDQUANTITY_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find total mass !
		CALL ANALYSIS(rho2, mass2)

		! Find total gravitational energy !
		CALL ANALYSIS(grav2, gravw2)

		! Find total angular momentum !
		CALL ANALYSIS(pphi2, j2)

		! Find the pressure integral !
		CALL ANALYSIS(p2, pint2)

		! Internal energy !
		CALL ANALYSIS(eps2, int2)

		! Moment of inertia !
		CALL ANALYSIS(i2, itotal2)

		! Mass quadrupole moment !
		CALL ANALYSIS(q2, qtotal2)

		! Tidal Love Number !
		IF(axratio2 == 1.0D0) THEN
			CALL TIDAL
		END IF

		! Global quantities !
		CALL FINDGLOBAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Output log !
		IF(success_flag .eqv. .true.) THEN
			IF (output_results .eqv. .true.) THEN	
				CALL OUTPUTLOG
			END IF
		END IF

		! Openfile !
		IF(output_profile .eqv. .true.) THEN
			CALL OPENFILE_PROFILE
		END IF

		! Output profile !
		IF(output_profile .eqv. .true.) THEN
			CALL OUTPUTPROFILE
		END IF

		! Close profile !
		IF(output_profile .eqv. .true.) THEN
			CALL CLOSEFILE_PROFILE
		END IF

		! Output global !
		IF(success_flag .eqv. .true.) THEN	
			CALL OUTPUTGLOBAL
		END IF

		! Print out !
		WRITE (*,*) 'Done for...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END DO

	! Close log file !
	IF (output_results .eqv. .true.) THEN
		CALL CLOSEFILE_LOG
	END IF

	! Openfile !
	CALL CLOSEFILEGLOBAL

END DO
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCF2F
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, l, m, n, o

! Dummy integers !
INTEGER :: raold, dra

! Conditional integers !
INTEGER :: count1, count2

! For bisection method !
REAL (DP) :: checkm, checkmlast, masstemp
REAL (DP) :: err_new, err_min

! For iterations !
REAL (DP) :: chatold1, chatnew1
REAL (DP) :: chatold2, chatnew2

! Potentail at boundaries point !
REAL (DP) :: psia1, psib1, phia1, phib1
REAL (DP) :: psia2, psib2, phia2, phib2
REAL (DP) :: fa1, fmax1
REAL (DP) :: fa2, fmax2

! Exit criteria !
REAL (DP) :: criteria1, criteria2, criteria3
REAL (DP) :: criteria4, criteria5, criteria6

! Condition !
LOGICAL :: condition1, condition2

! Sign and temporal step size !
REAL (DP) :: signs, dtemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the loop in generating a lot of models !
DO l = 1, n_rho
	
	! Assign NM maximum density !
	rhoscale2 = rhostart + (DBLE(l) - 1.0D0)*drho
	rhomax2 = 1.0D1**(rhoscale2)

	! Maximum Enthalpy !
	CALL MAXENTHALPY

	! Openfile !
	IF (output_results .eqv. .true.) THEN
		CALL OPENFILE_LOG
	END IF
	CALL OPENFILEGLOBAL

	! Do the loop in axis ratio !
	DO m = axstart, axend, daxratio

		! Get NM boundary !
		rb2 = m

		! Exit if the boundary is out of range !
		IF(rb2 < 1) THEN
			EXIT
		END IF

		! Find axis ratio !
		axratio2 = (r(rb2)/r(ra2))

		! Print out !
		WRITE (*,*) 'Start iterations for ...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

		! Set the initial DM boundary !
		ra1 = ra2 !NDIV

		! Set initial step size !
		dsize = -NINT(dgrid*DBLE(NDIV))

		! Set an arbitary large error !
		err_new = 1.0D30
		err_min = 1.0D30

		! Initialize !
		masstemp = 0.0D0
		count1 = 0
		count2 = 0
		condition1 = .false.
		condition2 = .false.	

		! Bisection method
		DO o = 0, omax

			! Get DM boundary !
			CALL DMAXIS

			! Exit condition !
			IF(ra1 <= 0) THEN
				success_flag = .false.
				WRITE (*,*) 'No converged model'
				WRITE (*,*) 'The last DM Mass', masstemp
				WRITE (*,*)
				EXIT
			END IF

			! Print out !
			WRITE (*,*) 'Start bisections for ...'
			WRITE (*,*) 'Size ratio', sizeratio
			WRITE (*,*) 'Initial step size', dsize
			WRITE (*,*) 'DM boundary', ra1
			WRITE (*,*) 'Iteration step', o
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Assign uniform density at the beginning of iteration 
			CALL INITIALRHO

			! Add the density !
			CALL GRAVRHO
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Find Potential !
			CALL FINDPOTENTIAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! For DM !
			phia1 = phi(1,ra1)
			phib1 = phi(KDIV,rb1)
			psia1 = psi1(1,ra1)
			psib1 = psi1(KDIV,rb1)

			! Potential ar boundary point !
			phia2 = phi(1,ra2)
			phib2 = phi(KDIV,rb2)
			psia2 = psi2(1,ra2)
			psib2 = psi2(KDIV,rb2)

			! Find h02 for NM !
			IF(axratio2 == 1.0D0) THEN
				h02new2 = 0.0D0
			ELSE
				h02new2 = -(phia2 - phib2)/(psia2 - psib2)
			END IF

			! Modified potential !
			CALL FINDFPOTENTIAL

			! Assign f-potential !
			fmax2 = maxval(fhat2)
			fa2 = fhat2(1,ra2)

			! Find chat for NM !
			chatnew2 = (hconst2*fmax2 - fa2)/(1.0D0 - hconst2)

			! Find chat for DM !
			chatnew1 = hconst1*(fmax2 + chatnew2) + phia1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Initialize !
			success_flag = .true.

			! Start the loop !
			DO n = 1, nmax

				! First backup !
				h02old2 = h02new2
				chatold2 = chatnew2
				hhmaxold2 = hhmaxnew2

				! For DM !
				h02old1 = h02new1
				chatold1 = chatnew1
				hhmaxold1 = hhmaxnew1

				! Find enthalpy !
				DO i = 1, KDIV
					DO j = 1, NDIV
						hhatnew2(i,j) = chatnew2 - phi(i,j) - h02new2*psi2(i,j)
					END DO
				END DO

				! Find enthalpy for DM !
				DO i = 1, KDIV
					DO j = 1, NDIV
						hhatnew1(i,j) = chatnew1 - phi(i,j) !- h02new1*psi1(i,j)
					END DO
				END DO

				! Maximum enthalpy for DM !
				hhmaxnew1 = maxval(hhatnew1)

				! Assign hhmax !
				hhmaxnew2 = maxval(hhatnew2)
				
				! Find equatorial radius !
				CALL FINDRADIUS

				! Find all the scaling constant !
				CALL FINDSCALE

				! Atmosphere !
				CALL ATMOSPHERE
				
				! Check the enthalpy !
				CALL CHECKENTHALPY

				! Get new density !
				CALL FINDDENSITY

				! Override DM enthalpy !
				hhatnew1(1,ra1) = hhata1
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!hhatnew1(KDIV,rb2) = hhata1
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				! Override NM enthalpy 	!
				hhatnew2(1,ra2) = hhata2
				hhatnew2(KDIV,rb2) = hhata2

				! Override DM density !
				rho1(1,ra1) = rhoa1
				!!!!!!!!!!!!!!!!!!!!!!!!
				!rho1(KDIV,rb2) = rhoa1
				!!!!!!!!!!!!!!!!!!!!!!!!

				! Override NM density !
				rho2(1,ra2) = rhoa2
				rho2(KDIV,rb2) = rhoa2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				! Total density !
				CALL GRAVRHO		

				! Find gravitational potential !
				CALL FINDPOTENTIAL
	
				! For DM !
				phia1 = phi(1,ra1)
				phib1 = phi(KDIV,rb1)
				psia1 = psi1(1,ra1)
				psib1 = psi1(KDIV,rb1)

				! Potential ar boundary point !
				phia2 = phi(1,ra2)
				phib2 = phi(KDIV,rb2)
				psia2 = psi2(1,ra2)
				psib2 = psi2(KDIV,rb2)
			
				! Find h02 !
				IF(axratio2 == 1.0D0) THEN
					h02new2 = 0.0D0
				ELSE
					h02new2 = -(phia2 - phib2)/(psia2 - psib2)
				END IF

				! Find modified potential !
				CALL FINDFPOTENTIAL

				! Assign f-potential !
				fmax2 = maxval(fhat2)
				fa2 = fhat2(1,ra2)

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Find h02 for DM !
				!h02new1 = -(phia1 - phib1)/(psia1 - psib1)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				! Find chat for NM !
				chatnew2 = (hconst2*fmax2 - fa2)/(1.0D0 - hconst2)

				! Find chat for DM !
				chatnew1 = hconst1*(fmax2 + chatnew2) + phia1
			
				! Exit conditions !
				IF(n > 1) THEN
					criteria1 = ABS((h02new2 - h02old2)/h02old2)
					criteria2 = ABS((chatnew2 - chatold2)/chatold2)
					criteria3 = ABS((hhmaxnew2 - hhmaxold2)/hhmaxold2)
					criteria4 = ABS((chatnew1 - chatold1)/chatold1)
					criteria5 = ABS((hhmaxnew1 - hhmaxold1)/hhmaxold1)
					IF(axratio2 == 1.0D0) THEN
						IF(criteria2 < tor .AND. criteria3 < tor .AND. criteria4 < tor .AND. criteria5 < tor) EXIT
					ELSE
						IF(criteria1 < tor .AND. criteria2 < tor .AND. criteria3 < tor .AND. criteria4 < tor .AND. criteria5 < tor) EXIT
					END IF
				END IF

				! Not sucessful if the model has NaN !
				IF(ieee_is_nan(chatnew1) .OR. ieee_is_nan(hhmaxnew1)) THEN
					WRITE (*,*) 'There is NaN'
					success_flag = .false.
					EXIT 
				END IF

				! Not sucessful if the model has NaN !
				IF(ieee_is_nan(h02new2) .OR. ieee_is_nan(chatnew2) .OR. ieee_is_nan(hhmaxnew2)) THEN
					WRITE (*,*) 'There is NaN'
					success_flag = .false.
					EXIT 
				END IF

				! Not sucessful if the model failed to converge !	
				If(n == nmax) THEN
					WRITE (*,*) 'Maximum iteration reached'
					success_flag = .false.
					EXIT 
				END IF

			END DO

			! Find total mass !
			CALL ANALYSIS(rho1, mass1)

			! DM mass in solar mass !
			masstemp = mass1*(mass/solar)
			WRITE (*,*) 'DM Mass', masstemp
			WRITE (*,*)

			! Check the deviation from target mass !
			checkmlast = checkm
			checkm = mass_dm - masstemp
			WRITE (*,*) 'DM Mass', masstemp
			WRITE (*,*) 'Mass Error', ABS(checkm/mass_dm)

			! Assign boundary !
			dra = ra1 - raold
			WRITE (*,*) 'Difference in DM boundary', dra
			WRITE (*,*)

			! Assign condition and error !
			IF(dra == -1) THEN
				count1 = count1 + 1
				err_new = ABS(checkm/mass_dm)
				err_min = MIN(err_min, err_new)
			END IF
			IF(dra == 1) THEN
				count2 = count2 + 1
				err_new = ABS(checkm/mass_dm)
				err_min = MIN(err_min, err_new)
			END IF
			IF(count1 >= 2) THEN
				condition1 = .true.
			END IF
			IF(count2 >= 2) THEN
				condition2 = .true.
			END IF
			
			IF((condition1 .eqv. .true.) .AND. (condition2 .eqv. .true.)) THEN
				IF(err_new == err_min .AND. ABS(dra) == 1) THEN
					abserror = err_min
					WRITE (*,*) 'Sucess, DM Mass', masstemp
					WRITE (*,*) 'Sucess, Error', err_new
					WRITE (*,*) 'Sucess, Min Error', err_min
					WRITE (*,*)
					EXIT
				END IF
			END IF
	
			! Convergence failure !
			IF(o == omax) THEN
				success_flag = .false.
				EXIT 
			END IF

			! Be sure to go into right direction !	
			IF(o == 0) THEN
				If(checkm > 0.0D0) THEN
					dsize = -dsize
				END IF
			END IF

			! Bisection !
               		if(checkmlast * checkm < 0.0D0 .and. o /= 0) then
				signs = sign(1.0D0, dsize)
                  		dtemp = DBLE(dsize) * dratio
				IF(ABS(dtemp) < 1) THEN
					dsize = -NINT(signs)
				ELSE
					dsize = -NINT(dtemp)
				END IF
			END IF

			! Backup DM radius !
			raold = ra1

			! Update DM radius !
			ra1 = ra1 + dsize

			! Override DM radius !
			IF(raold /= 1 .AND. ra1 <= 1) THEN
				ra1 = 3
			ELSEIF(raold /= NDIV .AND. ra1 >=NDIV) THEN
				ra1 = NDIV - 2
			END IF

		END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Be sure to check density !
		IF(success_flag .eqv. .true.) THEN
			CALL CHECKMODEL
		END IF

		! Exit if failed to converge !
		IF(critical_flag .eqv. .true.) THEN
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Critical rotation encountered'
				WRITE (*,*) 'Move on to the next density'
				WRITE (*,*) 
				EXIT
			END IF
		END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find rotational velocity and momentum !
		CALL FINDVELOCITY

		! Find remaining physical quantity !
		CALL FINDQUANTITY_1

		! Find total kinetic energy !
		kine1 = 0.0D0
		CALL ANALYSIS(ke2, kine2)

		! Find total gravitational energy !
		CALL ANALYSIS(grav, gravw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find the individual potential !
		CALL POTENTIAL_FINAL	

		! Find volume !
		CALL FINDVOLUME

		! Find pressure !
		CALL FINDPRESSURE

		! Find energy !
		CALL FINDEPSILON

		! Find remaining physical quantity !
		CALL FINDQUANTITY_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Find total mass !
		CALL ANALYSIS(rho2, mass2)

		! Find total gravitational energy !
		CALL ANALYSIS(grav2, gravw2)
		CALL ANALYSIS(grav1, gravw1)

		! Find total angular momentum !
		j1 = 0.0D0
		CALL ANALYSIS(pphi2, j2)

		! Find the pressure integral !
		CALL ANALYSIS(p2, pint2)
		CALL ANALYSIS(p1, pint1)

		! Internal energy !
		CALL ANALYSIS(eps2, int2)
		CALL ANALYSIS(eps1, int1)

		! Moment of inertia !
		CALL ANALYSIS(i1, itotal2)
		CALL ANALYSIS(i2, itotal2)

		! Mass quadrupole moment !
		CALL ANALYSIS(q1, qtotal2)
		CALL ANALYSIS(q2, qtotal2)

		! Tidal Love Number !
		IF(axratio2 == 1.0D0) THEN
			CALL TIDAL
		END IF

		! Global quantities !
		CALL FINDGLOBAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Output log !
		IF(success_flag .eqv. .true.) THEN
			IF (output_results .eqv. .true.) THEN	
				CALL OUTPUTLOG
			END IF
		END IF

		! Openfile !
		IF(output_profile .eqv. .true.) THEN
			CALL OPENFILE_PROFILE
		END IF

		! Output profile !
		IF(output_profile .eqv. .true.) THEN
			CALL OUTPUTPROFILE
		END IF

		! Close profile !
		IF(output_profile .eqv. .true.) THEN
			CALL CLOSEFILE_PROFILE
		END IF

		! Output global !
		IF(success_flag .eqv. .true.) THEN	
			CALL OUTPUTGLOBAL
		END IF

		! Print out !
		WRITE (*,*) 'Done for...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END DO

	! Close log file !
	IF (output_results .eqv. .true.) THEN
		CALL CLOSEFILE_LOG
	END IF

	! Openfile !
	CALL CLOSEFILEGLOBAL

END DO
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
