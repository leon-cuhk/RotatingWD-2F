SUBROUTINE OPENFILE_LOG
USE DEFINITION
IMPLICIT NONE

! Open file !
IF(dm_flag == 1) THEN
	OPEN (UNIT = 11, FILE = './Parameter/Star_WENO_Results_NMCentralDensity_'//trim(str(rhoscale2))//'_DMMass_'//trim(str(mass_dm))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 12, FILE = './Parameter/Star_WENO_Results_NMCentralDensity_'//trim(str(rhoscale2))//'_DMMass_'//trim(str(mass_dm))//'_DM.dat', STATUS = 'REPLACE')
ELSE
	OPEN (UNIT = 11, FILE = './Parameter/Star_WENO_Results_NMCentralDensity_'//trim(str(rhoscale2))//'_NM.dat', STATUS = 'REPLACE')
END IF

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILE_LOG
USE DEFINITION
IMPLICIT NONE

! Close file !
IF(dm_flag == 1) THEN
	CLOSE (11)
	CLOSE (12)
ELSE
	CLOSE (11)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILEGLOBAL
USE DEFINITION 
IMPLICIT NONE

! Open file !
IF(dm_flag == 1) THEN
	OPEN (UNIT = 21, FILE = './Parameter/Star_WENO_Global_NMCentralDensity_'//trim(str(rhoscale2))//'_DMMass_'//trim(str(mass_dm))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 22, FILE = './Parameter/Star_WENO_Global_NMCentralDensity_'//trim(str(rhoscale2))//'_DMMass_'//trim(str(mass_dm))//'_DM.dat', STATUS = 'REPLACE')
ELSE
	OPEN (UNIT = 21, FILE = './Parameter/Star_WENO_Global_NMCentralDensity_'//trim(str(rhoscale2))//'_NM.dat', STATUS = 'REPLACE')
END IF

! Header !
WRITE (21, *) '------------------------------------------------------------------------------------------------------------'
WRITE (21, *) ' masst, m2out, reout, j2out, h0out, wa2, va2, stab, bind, axratio2, ibar, qbar, rbar, rhomax2, tidal, viral '
WRITE (21, *) '------------------------------------------------------------------------------------------------------------'

! For DM !
IF(dm_flag == 1) THEN
	WRITE (22, *) '----------------------------------------------'
	WRITE (22, *) ' masst, m1out, rhoscale1, re1, stab, axratio1 '
	WRITE (22, *) '----------------------------------------------'
END IF

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILEGLOBAL
USE DEFINITION 
IMPLICIT NONE

! Close file !
IF(dm_flag == 1) THEN
	CLOSE (21)
	CLOSE (22)
ELSE
	CLOSE (21)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILE_PROFILE
USE DEFINITION
IMPLICIT NONE

! Open file !
IF(dm_flag == 1) THEN
	OPEN (UNIT = 41, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 42, FILE = './Profile/Star_WENO_Enthalpy_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 43, FILE = './Profile/Star_WENO_Potential_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_NM.dat', STATUS = 'REPLACE')

	! For DM !
	OPEN (UNIT = 51, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_DM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 52, FILE = './Profile/Star_WENO_Enthalpy_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_DM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 53, FILE = './Profile/Star_WENO_Potential_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_DM.dat', STATUS = 'REPLACE')

	! For Total !
	OPEN (UNIT = 61, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_TOTAL.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 62, FILE = './Profile/Star_WENO_Potential_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DMMass_'//trim(str(mass_dm))//'_TOTAL.dat', STATUS = 'REPLACE')
ELSE
	OPEN (UNIT = 41, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 42, FILE = './Profile/Star_WENO_Enthalpy_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_NM.dat', STATUS = 'REPLACE')
	OPEN (UNIT = 43, FILE = './Profile/Star_WENO_Potential_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_NM.dat', STATUS = 'REPLACE')
END IF

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILE_PROFILE
USE DEFINITION
IMPLICIT NONE

! Close file !
IF(dm_flag == 1) THEN
	CLOSE (41)
	CLOSE (42)
	CLOSE (43)

	! Close file !
	CLOSE (51)
	CLOSE (52)
	CLOSE (53)

	! Close file !
	CLOSE (61)
	CLOSE (62)
ELSE
	CLOSE (41)
	CLOSE (42)
	CLOSE (43)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILEDMNM
USE DEFINITION 
IMPLICIT NONE

! Open file !
OPEN (UNIT = 71, FILE = './Parameter/Star_WENO_DMNM_DMMass_'//trim(str(targetdm))//'.dat', STATUS = 'REPLACE')

! Header !
WRITE (71, *) '-------------------------------------'
WRITE (71, *) ' m1out, m2out, masst, abserror, stab '
WRITE (71, *) '-------------------------------------'

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILEDMNM
USE DEFINITION
IMPLICIT NONE

! Close file !
CLOSE (71)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!