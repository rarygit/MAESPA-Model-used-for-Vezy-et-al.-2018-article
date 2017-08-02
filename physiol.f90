! PHYSIOL.FOR

!=======================================================================================
! Copyright 2015 Remko Duursma, Belinda Medlyn, Mathias Christina, Guerric le Maire
!---------------------------------------------------------------------------------------
! this file is part of MAESPA.
!
! MAESPA is free software: you can redistribute it and/or modify
! it under the terms of the gnu general public license as published by
! the free software foundation, either version 2 of the license, or
! (at your option) any later version.
!
! MAESPA is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the
! gnu general public license for more details.
!
! you should have received a copy of the gnu general public license
! along with MAESPA.  if not, see <http://www.gnu.org/licenses/>.
!=======================================================================================
    
!**********************************************************************
! This file contains all physiology calculations - photosynthesis,
! stomatal conductance, transpiration, respiration.
!
! The main subroutines (called externally) are:
! PSTRANSP - calls the other functions, does the leaf temp iteration calcn
! RESP - calculates maintenance respiration rate
! GRESP - calculates growth respiration rate
! CALCRMW - calculates stem maintenance respiration rate per unit mass
! CALCFBIOM - calculates foliar biomass from leaf area & SLA
! CALCWBIOM - calculates woody biomass from height & diameter
! For water balance: not finished
!   ETCAN - calculates canopy transpiration rate
!   PartitionPPT - partitions precip between drainage & canopy storage
!
! Subsidiary subroutines are
! 1. Photosynthesis
!   PHOTOSYN - calculates photosynthesis from the FvC model
!   GAMMAFN - calculates T dependence of Gamma*
!   KMFN - calculates T dependence of Km
!   ARRH - Arrhenius T dependence
!   JMAXTFN - calculates T dependence of Jmax
!   VCMAXTFN - calculates T dependence of Vcmax
! 2. Conductances  & transpiration
!   GBHFREE - calculates conductance to heat through free convection
!   GBHFORCED - calculates conductance to heat through forced convection
!   GRADIATION - calculates radiation conductance
!   GBCAN - calculates boundary layer conductance of canopy
!   PENMON - implements the Penman-Monteith equation
!**********************************************************************


!**********************************************************************
SUBROUTINE PSTRANSPIF(IDAY,IHOUR,RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25,   &
                    IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,K10F,RTEMP, &
                    DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,SMD1,SMD2,WC1,WC2,SOILDATA,    &
                    SWPEXP,FSOIL,GSMIN,GNIGHT,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,VPARA,VPARB,VPARC,     &
                    VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT,TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP,       &
                    WEIGHTEDSWP,KTOT,HMSHAPE,PSIL,ETEST,ETDEFICIT,CI,ISMAESPA,ISNIGHT,G02,G12,NEWTUZET,EV,  &
                    drycan,CANOPY_STORE_I) !glm canopy evap
!
! 'Interface' to PSTRANSP (new subroutine, Feb. 2011). 
! Calculates (numericall) the leaf water potential for the Tuzet model; 
! otherwise (at the moment), proceeds to call PSTRANSP.
!**********************************************************************
      
    USE maestcom
    IMPLICIT NONE

    INTEGER MODELGS,SOILDATA,WSOILMETHOD,ITER
    INTEGER IECO,ITERMAX,NSIDES,VFUN
    INTEGER IDAY,IHOUR
    REAL JMAX25,I0,LHV,MINROOTWP,KTOT,PSIL,K10F,GNIGHT
    REAL TLEAF,TAIR,DLEAF,VPD,VMLEAF,VMFD,RHLEAF,RH,CS,CA
    REAL SLOPE,GRADN,RDFIPT,TUIPT,TDIPT,GBHU,PRESS,WIND
    REAL WLEAF,PAR,TMOVE,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC
    REAL DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP
    REAL TBELOW,GSREF,GSMIN,D0,VK1,VK2,VPD1,VPD2,VMDF0
    REAL GSJA,GSJB,T0,TREF,TMAX,SOILMOISTURE,EMAXLEAF
    REAL SMD1,SMD2,WC1,WC2,SWPEXP,FSOIL,G0,D0L,GAMMA,G1
    REAL GSC,ALEAF,RD,WEIGHTEDSWP,GBHF,GBH,GH,VMFD0,GBV,GSV,GV
    REAL ET,RNET,GBC,TDIFF,TLEAF1,FHEAT,ETEST,SF,PSIV,HMSHAPE
    REAL PSILIN,PLANTK,TOTSOILRES,MINLEAFWP,CI
    REAL TMP,VPARA,VPARB,VPARC,VPDMIN,GK,ETDEFICIT
    REAL G02,G12
    REAL EV,drycan,CANOPY_STORE_I !glm canopy evap
    INTEGER NEWTUZET
    LOGICAL ISMAESPA,ISNIGHT


    ! Find leaf water potential that matches Tuzet model (tuzet gs = f(psi), and psi = f(gs)).
    ! The result is PSILIN, which is then used below to (re-)estimate all gas exchange variables.
    IF(MODELGS.EQ.6)THEN
    
        ! Do not find leaf water potential at night.
        IF(.NOT.ISNIGHT)THEN
            CALL PSILFIND(RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25,&
                    IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0, &
                    Q10F,K10F,RTEMP,DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,    &
                    SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,GSMIN,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,   &
                    VPARA,VPARB,VPARC,VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT,TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP,  &
                    WEIGHTEDSWP,HMSHAPE,PSILIN,ETDEFICIT,ETEST, IDAY, IHOUR,G02,G12,NEWTUZET, &
                    EV,drycan,CANOPY_STORE_I) !glm canopy evap
        ELSE
            PSILIN = WEIGHTEDSWP   ! Not entirely correct due to nighttime transpiration, but has no consequences.
        ENDIF
    
    ENDIF

    CALL PSTRANSP(IDAY,IHOUR,RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25,&
                    IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0, &
                    Q10F,K10F,RTEMP,DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,    &
                    SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,GSMIN,GNIGHT,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,   &
                    VPARA,VPARB,VPARC,VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT,TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP,   &
                    WEIGHTEDSWP,KTOT,HMSHAPE,PSILIN,PSIL,ETDEFICIT,ETEST,CI,ISMAESPA,ISNIGHT,G02,G12,NEWTUZET, &
                    EV,drycan,CANOPY_STORE_I) !glm canopy evap
    

END SUBROUTINE PSTRANSPIF

!**********************************************************************
SUBROUTINE PSTRANSP(iday,ihour,RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25,&
                    IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0, &
                    Q10F,K10F,RTEMP,DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,    &
                    SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,GSMIN,GNIGHT,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,   &
                    VPARA,VPARB,VPARC,VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT, &
                    TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP,  &
                    WEIGHTEDSWP,KTOT,HMSHAPE,PSILIN,PSIL,ETDEFICIT,ETEST,CI,ISMAESPA,ISNIGHT, &
                    G02,G12,NEWTUZET,EV,drycan,CANOPY_STORE_I) !glm canopy evap
! This subroutine calculates leaf photosynthesis and transpiration.
! These may be calculated by
! (1) assuming leaf temperature = air temperature, Cs = Ca and Ds = Da
! (2) using iterative scheme of Leuning et al (1995) (PC&E 18:1183-1200)
! to calculate leaf temp, Cs & Ca.
! Setting ITERMAX = 0 gives (1); ITERMAX > 0 (suggest 100) gives (2).
!**********************************************************************
      
    USE maestcom
    IMPLICIT NONE

    INTEGER MODELGS,SOILDATA,WSOILMETHOD,ITER
    INTEGER IECO,ITERMAX,NSIDES,VFUN
    INTEGER IDAY,IHOUR
    REAL JMAX25,I0,LHV,MINROOTWP,KTOT,PSIL,K10F
    REAL PLANTK,TOTSOILRES,MINLEAFWP,GNIGHT
    REAL TLEAF,TAIR,DLEAF,VPD,VMLEAF,VMFD,RHLEAF,RH,CS,CA
    REAL SLOPE,GRADN,RDFIPT,TUIPT,TDIPT,GBHU,PRESS,WIND
    REAL WLEAF,PAR,TMOVE,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC
    REAL DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP
    REAL TBELOW,GSREF,GSMIN,D0,VK1,VK2,VPD1,VPD2,VMDF0
    REAL GSJA,GSJB,T0,TREF,TMAX,SOILMOISTURE,EMAXLEAF
    REAL SMD1,SMD2,WC1,WC2,SWPEXP,FSOIL,G0,D0L,GAMMA,G1
    REAL GSC,ALEAF,RD,WEIGHTEDSWP,GBHF,GBH,GH,VMFD0,GBV,GSV,GV
    REAL ET,RNET,GBC,TDIFF,TLEAF1,FHEAT,ETEST,SF,PSIV,HMSHAPE
    REAL PSILIN,CI,VPARA,VPARB,VPARC,VPDMIN,GK,RD0ACC,ETDEFICIT
    REAL EV, drycan,DRYCAN_I,GVEV,CANOPY_STORE_I  !glm canopy evap
    REAL TLEAF2,DELTATLEAFIN,DELTATLEAF,DELTATLEAF2 !glm secant method of convergence
    
    LOGICAL ISMAESPA,ISNIGHT,FAILCONV
    INTEGER NEWTUZET
    REAL G02,G12 
    CHARACTER*70 errormessage
    
    REAL, EXTERNAL :: SATUR
    REAL, EXTERNAL :: GRADIATION
    REAL, EXTERNAL :: GBHFORCED
    REAL, EXTERNAL :: GBHFREE
    REAL, EXTERNAL :: PENMON
    REAL, EXTERNAL :: RESP

    FAILCONV = .FALSE.
    
    ! Hydraulic conductance of the entire soil-to-leaf pathway
    IF(ISMAESPA)THEN
        KTOT = 1./(TOTSOILRES + 1./PLANTK)
    ENDIF

    ! Set initial values of leaf temp and surface CO2 & VPD
    TLEAF = TAIR
    DLEAF = VPD
    VMLEAF = VMFD
    RHLEAF = RH
    CS = CA

    ! Following calculations do not depend on TLEAF
    ! Latent heat of water vapour at air temperature (J mol-1)
    LHV = (H2OLV0 - 2.365E3 * TAIR) * H2OMW
    ! Const s in Penman-Monteith equation  (Pa K-1)
    SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
    ! Radiation conductance (mol m-2 s-1)
    GRADN = GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
    
    ! Boundary layer conductance for heat - single sided, forced convection
    GBHU = GBHFORCED(TAIR,PRESS,WIND,WLEAF)

    !**********************************************************************
    ITER = 0  ! Counter for iterations - finding leaf temperature
100 CONTINUE  ! Return point for iterations

    
!    TLEAF = TAIR
!    TLEAF2 = TAIR+1 
!    DELTATLEAF = 0.01     ! normally should run the code once with TLEAF initial values before
!    DELTATLEAF2 = -0.1    ! normally should run the code once with TLEAF2 initial values before
!    DELTATLEAFIN = 1.0
!    
!    ITERTLEAF = 0  ! Counter for iterations 
!    DO WHILE (((ABS(DELTATLEAFIN).GT.(0.02)).OR.(ABS(DELTATLEAF).GT.(0.02))).AND.(ITERTLEAF.LE.ITERMAX))   
!        
!    !Secant method for convergence of TLEAF
!   IF (DELTATLEAF.NE.DELTATLEAF2) THEN
!    DELTATLEAFIN = (TLEAF2-TLEAF)*DELTATLEAF/(DELTATLEAF-DELTATLEAF2) !glm secant method        
!    TLEAF2=TLEAF          ! xl=rtsec
!    TLEAF=TLEAF+DELTATLEAFIN     ! rtsec = rtsec+dx
!    DELTATLEAF2=DELTATLEAF    ! fl=f
!   ELSE
!        DELTATLEAFIN=0
!        DELTATLEAF=0
!        DELTATLEAF2=0
!    ENDIF
 
    
    
    
    IF(.NOT.ISNIGHT)THEN
    CALL PHOTOSYN(PAR,TLEAF,TMOVE,CS,RHLEAF,DLEAF,VMLEAF,JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,&
                    EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,K10F,RTEMP,DAYRESP,TBELOW,&
                    MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0,GSJA,GSJB,T0,TREF,TMAX,&
                    WSOILMETHOD,SOILMOISTURE,EMAXLEAF,SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,&
                    G0,D0L,GAMMA,VPDMIN,G1,GK,GSC,ALEAF,RD,MINLEAFWP,KTOT,WEIGHTEDSWP, & 
                    VPARA,VPARB,VPARC,VFUN,SF,PSIV,HMSHAPE,PSILIN,PSIL,CI,ISMAESPA,G02,G12,NEWTUZET)
    ELSE
        ! Nighttime calculations
        GSC = GNIGHT
        
        ! Maximum transpiration rate
        EMAXLEAF = KTOT * (WEIGHTEDSWP - MINLEAFWP)
        IF(EMAXLEAF.LT.0.0)EMAXLEAF = 0.0
        
        ! Hydraulic constraints on night-time conductance
        ETEST = 1000 * (VPD/PATM) * GSC * GSVGSC
        
        IF(ETEST.GT.EMAXLEAF)THEN
            GSV = 1E-03 * EMAXLEAF / (VPD/PATM)
            GSC = GSV / GSVGSC
            
            ! Cuticular conductance
            IF(GSC.LT.GSMIN)THEN
                GSC = GSMIN
                GSV = GSVGSC*GSC 
            ENDIF
        ENDIF
        
        CI = CS
        RD = RESP(RD0,RD0ACC,TLEAF,TMOVE,Q10F,K10F,RTEMP,DAYRESP,TBELOW)
        ALEAF = -RD
    ENDIF
    
    ! Boundary layer conductance for heat - single sided, free convection
    GBHF = GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
    
    ! Total boundary layer conductance for heat
    GBH = GBHU + GBHF
    ! Total conductance for heat - two-sided
    GH = 2.*(GBH + GRADN)
    ! Total conductance for water vapour
    GBV = GBVGBH*GBH
    
    ! Stomatal conductance to water vapour
    GSV = GSVGSC*GSC
    
    ! Boundary + stomatal conductance to water vapour
    GV = (GBV*GSV)/(GBV+GSV)

    ! Call Penman-Monteith equation !(if leaf 100% dry)  !glm canopy evap
    ET = PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GV)
    
    ! calcul 2 !glm
    ! GAMMA = CPAIR*AIRMA*PRESS/LHV
    ! IF (GV.GT.0.0) THEN
    !     ET = (CPAIR * AIRMA / GAMMA) * (SATUR(TLEAF) - (SATUR(TAIR) - VPD))/((1/GH)+(1/GV))
    !     ET = ET/LHV
    ! ELSE
    !     ET = 0.0
    ! ENDIF
    
    ! Call Penman-Monteith equation for evaporation of leaf surface water (if 100% wet leaf)  !glm canopy evap
    GVEV = 1./(1./(1E09) + 1./GBV) !stomatal infinite conductance + boundary  !glm canopy evap
    EV = PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GVEV)  !glm canopy evap

    !calcul 2 !glm
    !  GAMMA = CPAIR*AIRMA*PRESS/LHV
    !  IF (GV.GT.0.0) THEN
    !     EV = (CPAIR * AIRMA / GAMMA) * (SATUR(TLEAF) - (SATUR(TAIR) - VPD))/((1/GH)+(1/GVEV))
    !     EV = EV/LHV
    !     !print*,ET,ET2        
    ! ELSE
    !     EV = 0.0
    ! ENDIF

    !modify EV to be <= CANOPY_STORE_I; modification through change in drycan parameter, to allow a compensation with ET !glm canopy evap
       
    ! IF (((1-drycan)*EV*SPERHR *18 * 1E-03).gt.CANOPY_STORE_I) THEN !to kg m-2 t-1       !glm canopy evap
    !     drycan = MIN(1.0,MAX(0.0,1-(CANOPY_STORE_I/(SPERHR *18 * 1E-03))/EV))                              !glm canopy evap
    ! ENDIF

    ! glm canopy evap
    ! RV: DRYCAN_I is now used here instead of drycan because it can be modified for the considered voxel only.
    ! EV in mol.m-2.s-1, convert to mm. CANOPY_STORE_I in mm
    IF (((1-drycan)*EV*SPERHR * 18 * 1E-03).gt.CANOPY_STORE_I) THEN       
        DRYCAN_I = MIN(1.0,MAX(0.0,1-(CANOPY_STORE_I/(SPERHR *18 * 1E-03))/EV))
    ELSE 
        DRYCAN_I = drycan
    ENDIF

    ! End of subroutine if no iterations wanted.
    IF (ITERMAX.EQ.0) GOTO 200
    

    ! Otherwise, calculate new TLEAF, DLEAF, RHLEAF & CS
    GBC = GBH/GBHGBC
    IF(.NOT.ISNIGHT)THEN
        CS = CA - ALEAF/GBC
    ENDIF
    !TDIFF = (RNET - ET*LHV) / (CPAIR * AIRMA * GH)  
    ! changed to weighted average between evap and transp in function of drycan ratio. Drycan computed at canopy scale (this is a limit) !glm canopy evap
    TDIFF = (RNET - DRYCAN_I*ET*LHV - (1-DRYCAN_I)*EV*LHV) / (CPAIR * AIRMA * GH)  !glm canopy evap
    
    TLEAF1 = TAIR + TDIFF/4 ! divide by 4 to slow down convergence and avoid big changes
!    DELTATLEAF = TLEAF-(TAIR + TDIFF)
    
    DLEAF = ET * PRESS / GV
    RHLEAF = 1. - DLEAF/SATUR(TLEAF1)
    VMLEAF = DLEAF/PRESS*1E-3
    
    ! Check to see whether convergence achieved or failed
    IF (ABS(TLEAF - TLEAF1).LT.TOL/4) GOTO 200

    IF (ITER.GT.ITERMAX) THEN
        !WRITE(ERRORMESSAGE, '(I4,A,I2,A)') IDAY,'  ', IHOUR, ' FAILED CONVERGENCE IN PSTRANSP'
        !CALL SUBERROR(ERRORMESSAGE,IWARN,0)
        !print*,TLEAF,TLEAF1,CS,PSIL
        FAILCONV = .TRUE.
	    GOTO 200
    END IF

    ! Update temperature & do another iteration
    TLEAF = TLEAF1
    ITER = ITER + 1
    GOTO 100
    
!    ITERTLEAF = ITERTLEAF + 1 !glm secant method
!    END DO !TLEAF WHILE LOOP !glm secant method
    
    
    ! Sensible heat flux
!200 FHEAT = RNET - LHV*ET
200 FHEAT = RNET - DRYCAN_I*LHV*ET - (1-DRYCAN_I)*LHV*EV ! glm canopy evap
    
    ! Transpiration minus supply by soil/plant (EMAX) must be drawn from plant reserve
    ETDEFICIT = (VPD/PATM) * GSV *1E03 - EMAXLEAF
    IF(ETDEFICIT.LT.0.0) ETDEFICIT = 0.0
    
    ! Return ET,EV in umol m-2 s-1
    ET = ET*1E6
    EV = EV*1E6 !glm canopy evap
    

    ! Calculate leaf water potential
    ! supposed to be the one of the dry part of the leaf !glm canopy evap!
    IF(ISMAESPA)THEN
       ! Used to use ET without boundary layer; not sure why.
       ! ETEST = 1E06 * (VPD/PATM) * GSV
       ! PSIL = WEIGHTEDSWP - (ETEST/1000)/KTOT
       ! PSIL = min(0.,WEIGHTEDSWP - (ET/1000)/KTOT) ! RV 04-2017. 
        PSIL = WEIGHTEDSWP - (ET/1000)/KTOT ! RV 04-2017. 
       ! PSIL can be positive when ET is negative (e.g. dew)
    ELSE
        PSIL = 0.0
    ENDIF
    
    ET = DRYCAN_I*ET ! (output = actual reduced Transp)!glm canopy evap 
    EV = (1-DRYCAN_I)*EV !(output = actual reduced evaporation)!glm canopy evap 

    
    RETURN
    END SUBROUTINE PSTRANSP


!**********************************************************************
SUBROUTINE PHOTOSYN(PAR,TLEAF,TMOVE,CS,RH,VPD,VMFD, &
        JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        THETA,AJQ,RD0,Q10F,K10F,RTEMP,DAYRESP,TBELOW, &
        MODELGS,GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0, &
        GSJA,GSJB,T0,TREF,TMAX,WSOILMETHOD,SOILMOISTURE, &
        EMAXLEAF,SMD1,SMD2,WC1,WC2, &
        SOILDATA,SWPEXP, &
        FSOIL,G0,D0L,GAMMA,VPDMIN,G1,GK,GS,ALEAF,RD,MINLEAFWP,KTOT,WEIGHTEDSWP, &
        VPARA,VPARB,VPARC,VFUN,SF,PSIV,HMSHAPE,PSILIN,PSIL,CI,ISMAESPA, &
        G02,G12,NEWTUZET)
! This subroutine calculates photosynthesis according to the ECOCRAFT
! agreed formulation of the Farquhar & von Caemmerer (1982) equations.
! Stomatal conductance may be calculated according to the Jarvis,
! Ball-Berry or BB-Leuning models.
! NB ALEAF is NET leaf photosynthesis.
! FSOIL is now output (RAD 2008).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER MODELGS,IQERROR,IECO,SOILDATA, WSOILMETHOD,IQERROR1
    INTEGER VFUN
    REAL PAR,TLEAF,CS,RH,VPD,VMFD,PSIL,KTOT
    REAL JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN
    REAL THETA,AJQ,RD0,Q10F,K10F,RTEMP,DAYRESP,TBELOW
    REAL GSREF,GSMIN,I0,D0,VK1,VK2,VPD1,VPD2,VMFD0
    REAL GSJA,GSJB,T0,TREF,TMAX
    REAL G0,D0L,GAMMA,G1,MINLEAFWP
    REAL GS,ALEAF,RD,MINROOTWP,RD0ACC
    REAL TMOVE,RESP,FSOIL,SOILMOISTURE,SMD1,SMD2,WC1,WC2,SWPEXP
    REAL VPDG,ETEST,WEIGHTEDSWP,EMAXLEAF,GSV
    REAL SF,PSIV,HMSHAPE,PSILIN,VJMOD
    REAL GAMMASTAR,KM,JMAX,VCMAX,J,VJ
    REAL A,B,C,AC,AJ,GSDIVA,CIC,CIJ,FPSIF
    REAL KMFN,JMAXTFN,CHK,CI,tmp,VPDMIN
    REAL VPARA,VPARB,VPARC,GK
    REAL G02, G12, GS2,ALEAF2, GSDIVA2
    INTEGER NEWTUZET
    LOGICAL ISMAESPA
    REAL, EXTERNAL :: GAMMAFN
    REAL, EXTERNAL :: VCMAXTFN
    REAL, EXTERNAL :: QUADM
    REAL, EXTERNAL :: CALCFSOIL
    REAL, EXTERNAL :: QUADP
    REAL, EXTERNAL :: FPSIL
    REAL, EXTERNAL :: VJMAXWFN


    ! Calculate photosynthetic parameters from leaf temperature.
    GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensati
    KM = KMFN(TLEAF,IECO)                             ! Michaelis-Ment
    JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Po
    VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Ma
    RD = RESP(RD0,RD0ACC,TLEAF,TMOVE,Q10F,K10F,RTEMP,DAYRESP,TBELOW)
    
    ! Effect of soil water stress on Vcmax and Jmax.
    IF(ISMAESPA)THEN
        VJMOD = VJMAXWFN(WEIGHTEDSWP, VPARA, VPARB, VPARC, VFUN)
        VCMAX = VCMAX * VJMOD
        JMAX = JMAX * VJMOD
    ENDIF
    
    ! Actual electron transport rate
    J = QUADM(THETA,-(AJQ*PAR+JMAX),AJQ*PAR*JMAX,IQERROR)
    
    ! RuBP regeneration rate
    VJ = J/4.0

    ! Deal with extreme cases
    IF ((JMAX.LE.0.0).OR.(VCMAX.LE.0.0)) THEN
        ALEAF = -RD        
        GS = G0
        ALEAF2 = -RD        
        GS2 = G0

        RETURN
    END IF

    ! Calculate soil moisture modifying factor
    ! MAESPA uses Emax approach, or Tuzet. MAESTRA calculates it various ways based on input data.
    IF(.NOT.ISMAESPA.AND.WSOILMETHOD.NE.0)THEN
        FSOIL = CALCFSOIL(WSOILMETHOD,SOILMOISTURE,SOILDATA,SMD1,SMD2,WC1,WC2,SWPEXP)
    ELSE
        FSOIL = 1.0
    ENDIF
    
        
        ! Note that MODELGS=5 is not implemented (but it is in Maestra).
        IF (MODELGS.EQ.2) THEN
            ! Ball-Berry model
            GSDIVA = G1 * RH / (CS - GAMMA) * FSOIL
            IF (NEWTUZET.EQ.1) THEN
                GSDIVA2 = G12 * RH / (CS - GAMMA) * FSOIL
            END IF
        ELSE IF (MODELGS.EQ.3) THEN
            ! Ball-Berry-Leuning model
            GSDIVA = G1 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
            IF (NEWTUZET.EQ.1) THEN
                GSDIVA2 = G12 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
            END IF
        ELSE IF (MODELGS.EQ.4) THEN
            IF(VPD.LT.VPDMIN)THEN
                VPDG = VPDMIN/1000.0
            ELSE
                VPDG =VPD/1000.0
            ENDIF
            ! Old, linearized version of BBOpti model
            !GSDIVA = G1 / (CS - GAMMA) / VPDG**(1-GK)
            
            ! Full version
            ! NOTE: 1.6 (from corrigendum to Medlyn et al 2011) is missing here,
            ! because we are calculating conductance to CO2!
            GSDIVA = (1.0 + G1 / VPDG**(1-GK)) / CS
            IF (NEWTUZET.EQ.1) THEN
                GSDIVA2 = (1.0 + G12 / VPDG**(1-GK)) / CS
            END IF
            
        ELSE IF (MODELGS.EQ.6) THEN
            IF(VPD.LT.VPDMIN)THEN
                VPDG = VPDMIN/1000.0
            ELSE
                VPDG =VPD/1000.0
            ENDIF
            FPSIF = FPSIL(PSILIN,SF,PSIV)
            GSDIVA = (G1 / (CS - GAMMA)) * FPSIF
            IF (NEWTUZET.EQ.1) THEN
                GSDIVA2 = (G12/ (CS -GAMMA)) * FPSIF
            END IF
        END IF
           ! print*,'GSDIVA',GSDIVA,CS,GAMMA

        ! Following calculations are used for both BB & BBL models.
        !-----------
        !glm WITH G1
        !-----------
        
        ! Solution when Rubisco activity is limiting
        A = G0 + GSDIVA * (VCMAX - RD)
        B = (1. - CS*GSDIVA) * (VCMAX - RD) + G0 * (KM - CS)- GSDIVA * (VCMAX*GAMMASTAR + KM*RD)
        C = -(1. - CS*GSDIVA) * (VCMAX*GAMMASTAR + KM*RD) - G0*KM*CS

        CIC = QUADP(A,B,C,IQERROR)

        IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
            AC = 0.0
        ELSE
            AC = VCMAX * (CIC - GAMMASTAR) / (CIC + KM)
        END IF
 
        ! Solution when electron transport rate is limiting
        A = G0 + GSDIVA * (VJ - RD)
        B = (1. - CS*GSDIVA) * (VJ - RD) + G0 * (2.*GAMMASTAR - CS) &
            - GSDIVA * (VJ*GAMMASTAR + 2.*GAMMASTAR*RD)
        C = -(1. - CS*GSDIVA) * GAMMASTAR * (VJ + 2.*RD) &
            - G0*2.*GAMMASTAR*CS
        CIJ = QUADP(A,B,C,IQERROR)

        AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        IF (AJ-RD.LT.1E-6) THEN        ! Below light compensation point
            CIJ = CS
            AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        END IF

        ALEAF = AMIN1(AC,AJ) - RD  ! Solution for Ball-Berry model
        GS = G0 + GSDIVA*ALEAF
        !print*,'ALEAF',ALEAF,AC,AJ,RD
        !-----------
        !glm WITH G2
        !-----------
        GS2 = G02 + GSDIVA2*ALEAF !test to speed up
        
        IF (GS2.GT.GS) THEN 
        ! Solution when Rubisco activity is limiting
        A = G02 + GSDIVA2 * (VCMAX - RD)
        B = (1. - CS*GSDIVA2) * (VCMAX - RD) + G02 * (KM - CS)- &
              GSDIVA2 * (VCMAX*GAMMASTAR + KM*RD)
        C = -(1. - CS*GSDIVA2) * (VCMAX*GAMMASTAR + KM*RD) - G02*KM*CS

        CIC = QUADP(A,B,C,IQERROR)

        IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
                    AC = 0.0
        ELSE
                    AC = VCMAX * (CIC - GAMMASTAR) / (CIC + KM)
        END IF
 
        ! Solution when electron transport rate is limiting
        A = G02 + GSDIVA2 * (VJ - RD)
        B = (1. - CS*GSDIVA2) * (VJ - RD) + G02 * (2.*GAMMASTAR - CS) &
                    - GSDIVA2 * (VJ*GAMMASTAR + 2.*GAMMASTAR*RD)
        C = -(1. - CS*GSDIVA2) * GAMMASTAR * (VJ + 2.*RD) &
                    - G02*2.*GAMMASTAR*CS
        CIJ = QUADP(A,B,C,IQERROR)

        AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        IF (AJ-RD.LT.1E-6) THEN        ! Below light compensation point
                    CIJ = CS
                    AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        END IF

        ALEAF2 = AMIN1(AC,AJ) - RD  ! Solution for Ball-Berry model
        GS2 = G02 + GSDIVA2*ALEAF2
        
        !-----------
        !glm take the max
        !-----------
        IF (GS2.GT.GS) THEN                        
                GS = GS2
                ALEAF = ALEAF2
        ENDIF  ! new model gs
        !-------------
        
        END IF !fin test
                    

        ! Set nearly zero conductance (for numerical reasons).
        IF (GS.LT.GSMIN) GS = GSMIN

        ! Maximum transpiration rate
        EMAXLEAF = KTOT * (WEIGHTEDSWP - MINLEAFWP)
        IF(EMAXLEAF.LT.0.0)EMAXLEAF = 0.0

        ! If E > Emax, set E to Emax, and corresponding gs and A.
        ! Do not use this routine when the Tuzet model of gs (6) is used.
        IF(ISMAESPA)THEN
            IF(WSOILMETHOD.EQ.1.AND.MODELGS.NE.6)THEN

                ! Leaf transpiration in mmol m-2 s-1  -  ignoring boundary layer effects!
                ETEST = 1000 * (VPD/PATM) * GS * GSVGSC

                ! Leaf water potential
                PSIL = WEIGHTEDSWP - ETEST/KTOT
    
                IF(ETEST > EMAXLEAF)THEN

                    ! Just for output:
                    FSOIL = EMAXLEAF / ETEST

                    ! Gs in mol m-2 s-1
                    GSV = 1E-03 * EMAXLEAF / (VPD/PATM)
                    GS = GSV / GSVGSC

                    ! GS cannot be lower than minimum (cuticular conductance)
                    IF(GS.LT.GSMIN)THEN
                        GS = GSMIN
                    ENDIF

                    ! Now that GS is known, solve for CI and A as in the Jarvis model.
                    ! Photosynthesis when Rubisco is limiting
                    A = 1./GS
                    B = (RD - VCMAX)/GS - CS - KM
                    C = VCMAX * (CS - GAMMASTAR) - RD * (CS + KM)

                    A = 1./GS
                    B = (0.0 - VCMAX)/GS - CS - KM
                    C = VCMAX * (CS - GAMMASTAR)
                    
                    AC = QUADM(A,B,C,IQERROR1)
      
                    ! Photosynthesis when electron transport is limiting
                    A = 1./GS
                    B = (RD - VJ)/GS - CS - 2*GAMMASTAR
                    C = VJ * (CS - GAMMASTAR) - RD * (CS + 2*GAMMASTAR)
                    AJ = QUADM(A,B,C,IQERROR)

                    ALEAF = AMIN1(AC,AJ)       ! Emax model solution.
                    
                ENDIF ! if (E>EMAX)
            ENDIF
        ENDIF
    
    ! Return CI.
    IF(GS.GT.0.AND.ALEAF.GT.0)THEN
        CI = CS - ALEAF/GS
    ELSE
        CI = CS
    ENDIF
    
    
    RETURN
END SUBROUTINE PHOTOSYN

!**********************************************************************
REAL FUNCTION GAMMAFN(TLEAF, IECO)
! This subroutine calculates Gamma(star), or the CO2 compensation point
! in the absence of non-photorespiratory respiration.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IECO
    REAL TLEAF
    REAL, EXTERNAL :: ARRH
    
    IF (IECO.EQ.1) THEN
        ! Ecocraft fomulation; based on Brooks & Farquhar and von Caemmerer et al.
        ! If TLEAF < -1.0 then calculate Gamma for T = -1 (quadratic not applicable)
        IF (TLEAF.LT.-1.0) THEN
            GAMMAFN = 36.9 + 1.88*(-26.0) + 0.036*(-26.0)*(-26.0)
        ELSE
            GAMMAFN = 36.9 + 1.88*(TLEAF-25) + 0.036*(TLEAF-25)*(TLEAF-25)
            
        END IF
    ELSE      ! Bernacchi et al 2001 PCE 24: 253-260
            GAMMAFN = ARRH(42.75,37830.0,TLEAF,25.0)
    ENDIF
    RETURN
END FUNCTION GAMMAFN

!**********************************************************************
REAL FUNCTION KMFN(TLEAF,IECO)
! This subroutine calculates Km, or the effective Michaelis-Menten
! coefficient of Rubisco activity.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER IECO
    REAL OI,KC25,KO25,KCEA,KOEA,KC,KO, TLEAF
    REAL, EXTERNAL :: ARRH
    
    OI = 205000         ! Oxygen partial pressure (umol mol-1)
    IF (IECO.EQ.1) THEN
    ! Physiological constants - values agreed by Ecocraft - Badger & Collatz values
        KC25 = 404          ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 248000       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 59400        ! Temp. response of Kc (J mol-1)
        KOEA = 36000        ! Temp. response of Ko (J mol-1)
    ELSE !  Bernacchi et al 2001 PCE 24: 253-260
        KC25 = 404.9        ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 278400       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 79430        ! Temp. response of Kc (J mol-1)
        KOEA = 36380        ! Temp. response of Ko (J mol-1)
    END IF

    
      ! Km <- exp(38.05-79430/(8.314*(Tleaf+273)))*(1+210/exp(20.3-36380/(8.314*(Tleaf+273))))
      ! ARRH = KT*EXP(EA*(T-TREF)/(RCONST*(T-ABSZERO)*(TREF-ABSZERO)))
      ! ARRH(KT,EA,T,TREF)
    
    ! This function is well-behaved for TLEAF < 0.0
    KC = ARRH(KC25,KCEA,TLEAF,25.0)
    KO = ARRH(KO25,KOEA,TLEAF,25.0)
    KMFN = KC * (1. + OI/KO)

    RETURN
END FUNCTION KMFN


!**********************************************************************
REAL FUNCTION JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)
! This subroutine calculates the potential electron transport rate
! (Jmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL JMAX25,TLEAF,EAVJ,EDVJ,DELSJ, TVJUP, TVJDN
    REAL, EXTERNAL :: TK
    
    ! This function is well-behaved for TLEAF < 0.0
    JMAXTFN = JMAX25 * EXP((TLEAF-25)*EAVJ/(RCONST*TK(TLEAF)*TK(25.)))  &
                    * (1.+EXP((DELSJ*TK(25.)-EDVJ)/(RCONST*TK(25.))))   &
                    / (1.+EXP((DELSJ*TK(TLEAF)-EDVJ)/(RCONST*TK(TLEAF))))

    ! Function allowing Vcmax to be forced linearly to zero at low T -
    ! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
        JMAXTFN = 0.0
    ELSE IF (TLEAF.LT.TVJUP) THEN
        JMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*JMAXTFN
    END IF

    RETURN
END FUNCTION JMAXTFN


!**********************************************************************
REAL FUNCTION VJMAXWFN(SWP, VPARA, VPARB, VPARC, VFUN)
! Dependance of Vcmax and Jmax on soil water potential (SWP).
! Modifier function (returns value between 0 and 1).
!**********************************************************************
USE maestcom
IMPLICIT NONE
REAL SWP, VPARA, VPARB, VPARC
INTEGER VFUN

IF(VFUN.EQ.0)THEN
    VJMAXWFN = 1
ENDIF

! Simple linear dependance. 
! VPARA is SWP where Vcmax is zero, VPARB is SWP where Vcmax is one.
IF(VFUN.EQ.1)THEN

    IF(SWP.LT.VPARA)THEN
        VJMAXWFN = 0
    ELSEIF(SWP.GT.VPARB) THEN
        VJMAXWFN = 1
    ELSE 
        VJMAXWFN =  (SWP - VPARA) / (VPARB - VPARA)
    ENDIF
ENDIF

! As in Zhou, et al. Agricultural and Forest Meteorology. 2013.
IF(VFUN.EQ.2)THEN
    
    VJMAXWFN = (1 + EXP(VPARA*VPARB) )/ (1 + EXP(VPARA*(VPARB-SWP)))
    
ENDIF


RETURN
END FUNCTION VJMAXWFN





!**********************************************************************
REAL FUNCTION VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)
! This subroutine calculates the maximum Rubisco activity
! (Vcmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL VCMAX25,TLEAF,EAVC,EDVC,DELSC
    REAL TVJDN,TVJUP
    REAL, EXTERNAL :: TK

    ! There is still disagreement as to whether this function has an
    ! optimum or not. Both forms are available here. If EDVC <= 0 (default 0)
    ! then the no-optimum form is used.
    ! Both functions are well-behaved for TLEAF < 0.0

    IF (EDVC.LE.0.0) THEN
        VCMAXTFN = VCMAX25 * EXP(EAVC*(TLEAF - 25)/ (TK(25.)*RCONST*TK(TLEAF)))
    ELSE
        VCMAXTFN = VCMAX25 * EXP((TLEAF-25.)*EAVC/(RCONST*TK(TLEAF)*TK(25.))) &
                    * (1.+EXP((DELSC*TK(25.)-EDVC)/(RCONST*TK(25.)))) &
                    / (1.+EXP((DELSC*TK(TLEAF)-EDVC)/(RCONST*TK(TLEAF))))
    END IF

    ! Function allowing Vcmax to be forced linearly to zero at low T -
    ! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
        VCMAXTFN = 0.0
    ELSE IF (TLEAF.LT.TVJUP) THEN
        VCMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*VCMAXTFN
    END IF

    RETURN
END FUNCTION VCMAXTFN


!**********************************************************************
REAL FUNCTION RESP(RD0,RD0ACC,TLEAF,TMOVE,Q10F,K10F,RTEMP,DAYRESP,TBELOW)
! This function calculates respiration from temperature
! using a Q10 (exponential) formulation.
!**********************************************************************

    IMPLICIT NONE      
    REAL RD0,TLEAF,Q10F,RTEMP,DAYRESP,K10F,RD0ACC,TMOVE,TBELOW

    IF (TLEAF.GE.TBELOW) THEN
        RD0ACC = RD0 * EXP(K10F*(TMOVE-RTEMP))
        RESP =  RD0ACC * EXP(Q10F * (TLEAF-RTEMP)) * DAYRESP
    ELSE
        RESP = 0.0
    END IF

    RETURN
END FUNCTION RESP


!**********************************************************************
REAL FUNCTION ARRH(KT,EA,T,TREF)
! The Arrhenius function.
! KT is the value at Tref deg C; Ea the activation energy (J mol-1) and T the temp (deg C).
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL KT, EA, T, TREF
    
    ARRH = KT*EXP(EA*(T-TREF)/(RCONST*(T-ABSZERO)*(TREF-ABSZERO)))
    RETURN
END FUNCTION ARRH


!**********************************************************************
REAL FUNCTION ETCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,RNET,VPD,GSCAN,STOCKING,TREEH,TOTLAI)
! Calculate transpiration by applying Penman-Monteith to whole canopy.
! Returns umol m-2 s-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL LHV,WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,RNET,VPD,GSCAN,STOCKING
    REAL GB,GSV,RNETM2,SLOPE,GH,GV,TREEH, TOTLAI, GBCANMS1,GBCANMS2,CMOLAR
    REAL, EXTERNAL :: GBCAN
    REAL, EXTERNAL :: HEATEVAP
    REAL, EXTERNAL :: SATUR
    REAL, EXTERNAL :: PENMON
    REAL, EXTERNAL :: TK

    ! Get boundary layer conductance, GBCANMS give conductance in m s-1
!    GB = GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
    CALL GBCANMS(WIND,ZHT,Z0HT,ZPD,TREEH,TOTLAI,GBCANMS1,GBCANMS2)
    CMOLAR = PRESS / (RCONST * TK(TAIR))
    GB = GBCANMS1*CMOLAR    ! in mol m-2 s-1
        
    ! Convert mol CO2/tree/s to mol H2O/m2/s
    GSV = GSCAN*GSVGSC*STOCKING
    RNETM2 = RNET*STOCKING

    IF (GB*GSV.GT.0.0) THEN
        ! Latent heat of water vapour at air temperature (J mol-1)
        LHV = HEATEVAP(TAIR) * H2OMW

        ! Const s in Penman-Monteith equation  (Pa K-1)
        SLOPE = (SATUR(TAIR + 0.1) - SATUR(TAIR)) / 0.1
        
        ! Call Penman-Monteith
        GH = GB
        GV = 1./(1./GSV + 1./GB)
        ETCAN = PENMON(PRESS,SLOPE,LHV,RNETM2,VPD,GH,GV)*1E6
    ELSE
        ETCAN = 0.0
    END IF

    RETURN
END FUNCTION ETCAN

!**********************************************************************
REAL FUNCTION PENMON(PRESS,SLOPE,LHV,RNET,VPD,GH,GV)
! This subroutine calculates evapotranspiration by leaves using the Penman-Monteith equation.
! Inputs:      PRESS atmospheric pressure, Pa
!            SLOPE slope of VPD/T curve, Pa K-1
!            LHV latent heat of water at air T, J mol-1
!            RNET net radiation, J m-2 s-1
!            VPD vapour pressure deficit of air, Pa
!            GH boundary layer conductance to heat (free & forced & radiative components), mol m-2 s-1
!            GV conductance to water vapour (stomatal & bdry layer components), mol m-2 s-1
! Result in mol H2O m-2 s-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL PRESS,SLOPE,LHV,RNET,VPD,GH,GV,GAMMA
    REAL ET
    
    GAMMA = CPAIR*AIRMA*PRESS/LHV

    IF (GV.GT.0.0) THEN
        ET = (SLOPE * RNET + VPD * GH * CPAIR * AIRMA) / (SLOPE + GAMMA * GH/GV)
    ELSE
        ET = 0.0
    END IF
    PENMON = ET / LHV
    
    IF (PENMON.LT.0.0) PENMON = 0.0            ! BM 12/05 Should not be negative
    
      RETURN
END FUNCTION PENMON


!**********************************************************************
REAL FUNCTION GRADIATION(TAIR,RDFIPT,TUIPT,TDIPT)
! Returns the 'radiation conductance' at given temperature.
! Formula from Ying-Ping's version of Maestro.
! See also Jones (1992) p. 108.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL TAIR,RDFIPT,TDIPT,TUIPT
    REAL, EXTERNAL :: TK

    GRADIATION = 4.*SIGMA*(TK(TAIR)**3.) * RDFIPT/TDIPT * EMLEAF * (TDIPT + TUIPT) &
                / (CPAIR * AIRMA)

    RETURN
END FUNCTION GRADIATION


!**********************************************************************
REAL FUNCTION GBHFORCED(TAIR,PRESS,WIND,WLEAF)
! Boundary layer conductance for heat - single sided, forced convection
! in mol m-2 s-1
! See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL TAIR,PRESS,WIND,WLEAF,CMOLAR
    REAL, EXTERNAL :: TK

    CMOLAR = PRESS / (RCONST * TK(TAIR))
    GBHFORCED = 0.003 * SQRT(WIND/WLEAF) * CMOLAR

    RETURN
END FUNCTION GBHFORCED


!**********************************************************************
REAL FUNCTION GBHFREE(TAIR,TLEAF,PRESS,WLEAF)
! Boundary layer conductance for heat - single sided, free convection
! in mol m-2 s-1
! See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL CMOLAR,PRESS,TAIR,TLEAF,GRASHOF,WLEAF
    REAL, EXTERNAL :: TK

    CMOLAR = PRESS / (RCONST * TK(TAIR))
    IF ((TLEAF-TAIR).NE.0.0) THEN
        GRASHOF = 1.6E8 * ABS(TLEAF-TAIR) * (WLEAF**3.) ! Grashof number
        GBHFREE = 0.5 * DHEAT * (GRASHOF**0.25) / WLEAF * CMOLAR
    ELSE
        GBHFREE = 0.0
    END IF

    RETURN
END FUNCTION GBHFREE


!**********************************************************************
REAL FUNCTION GBCAN(WIND,ZHT,Z0HT,ZPD,PRESS,TAIR)
! Canopy boundary layer conductance (from Jones 1992 p 68)
! in mol m-2 s-1
! ZHT =
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL WIND,ZHT,Z0HT,ZPD,PRESS,TAIR,CMOLAR
    REAL, EXTERNAL :: TK
    
    IF (Z0HT.GT.0.0) THEN

        ! Formula from Jones 1992 p 68
        GBCAN = WIND*(VONKARMAN**2)/(LOG((ZHT - ZPD)/Z0HT))**2
        
        ! Convert from m s-1 to mol m-2 s-1
        CMOLAR = PRESS / (RCONST * TK(TAIR))
        GBCAN = GBCAN * CMOLAR
    ELSE
        GBCAN = 0.0
    END IF

    RETURN
END FUNCTION GBCAN


!**********************************************************************
SUBROUTINE GBCANMS(WIND,ZHT,Z0HT,ZPD, TREEH, TOTLAI, GBCANMS1, GBCANMS2)
! Canopy boundary layer conductance (from Jones 1992 p 68)
! in m s-1
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL WIND,ZHT,ZSTAR,Z0HT,ZPD,WIND2
    REAL Cd, X, TOTLAI, ZPD2, TREEH, Z0, KH, ALPHA, Z0HT2
    REAL GBCANMS1, GBCANMS2, GBCANMSINI, GBCANMSROU
    REAL ALPHA1, WINDSTAR, ZW, COAT, GBCANMS3,Z0H,GBCANMS1MIN
    
    ! In this model, we assumed 2 aerodynamic conductances in series
    ! 1) from the atmosphere to the canopy, based on Van de Griend 1989
    !       this is actually 2 conductances, one in the inertial sublayer and one in the roughness sublayer
    ! 2) Within the canopy to the soil, based on CHoudhury & Monteith 1988
    
    ! Aerodynamic conductance from the canopy to the atmosphere
    ! Formula from Jones 1992 p 68, aerodynamic conductance air-canopy - air, adapted from the CASTANEA model (Dufrene et al., 2005)
    ! OLIVER & MAYHEAD (1974) :
    ZPD2 = 0.75 * TREEH
    Z0 = 0.1 *  TREEH
    ! ZPD2 = 0.67 * TREEH
    ! Z0 = 0.046 *  TREEH
    ! RV: Z0 and ZPD2 are unknown... Carefull, model is sensible to z0.
    ! NB: Lettau (1969) proposed an other way: Z0= 0.5 . h* . s / S where h* is canopy height, s is
    ! the average silhouette area (projected area of the tree on a vertical plane) and S the specific 
    ! area, with S= A/n, where A is the total plot area and n the number of trees.
    ! For ZPD2, Verhoef et al. (1997) said d= 0.67 . h is a good proxy without any prior knowledge. 
    
    ZSTAR=ZHT
    WIND2= WIND
    ! ZSTAR can be corrected if lower than maximum tree height (TREEH) because GBCANMS have to be computed
    ! for the whole canopy.
    
    ! RV: Observations indicate that z* lies 1–2 times the height of the canopy above the canopy (Garratt 1992, Harman&Finnigan 2007).
    ! In their paper z* is the height of the roughness sublayer.
    
    
    ! RV & GLM 05/2017 : If ZSTAR (=ZHT) is below TREEH (e.g. wind measurements are made under canopy)
    ! wind is increased to obtain a wind above the canopy for GBH (Wind in the roughness layer). 
    IF (ZSTAR.LE.TREEH) THEN
        WIND2 = WIND * EXP(0.13155 * (TREEH/ZSTAR -1))
        ZSTAR= TREEH*2.0
    ENDIF
    
    ! Aerodynamic conductance between the atmosphere and the canopy
    ! Reference Wind used in the conductance calculation
    ! (this is ustar, the friction velocity)
    WINDSTAR = WIND2 * VONKARMAN / log((ZSTAR-ZPD2)/Z0)
    
    ! We supposed 2 aerodynamic conductances, one in the inertial layer (from ZSTAR to ZW)
    ! and another one in the roughness layer from ZSTAR to TREEH
    ! According to Van de Griend 1989, we can assumed that :
    ALPHA1 = 1.5
    ZW = ZPD2 + ALPHA1 * (TREEH-ZPD2)
    
    ! Aerodynamic conductance in the inertial sublayer (Van de Griend 1989)
    GBCANMSINI = WINDSTAR*VONKARMAN /(LOG((ZSTAR - ZPD2)/(ZW - ZPD2)))
    ! Aerodynamic conductance in the roughness layer
    ! The roughness layer is located between TREEH and a height ZW, according to 
    ! GBCANMSROU = WINDSTAR*VONKARMAN * ((ZW - TREEH)/(ZW - ZPD2))
    GBCANMSROU = WINDSTAR*VONKARMAN / ((ZW - TREEH)/(ZW - ZPD2)) !glm 03/2016
    
    ! Total aerodynamic conductance between the canopy and the atmosphere
    GBCANMS1MIN = 0.0123 ! RV, GBCANMS1 for WIND= 0.035 and CANOPY HEIGHT at 25 m
    GBCANMS1 = 1/ (1/GBCANMSINI + 1/GBCANMSROU)
    IF (GBCANMS1.LE.GBCANMS1MIN) THEN
        GBCANMS1 = GBCANMS1MIN
    ENDIF
           
    ! Aerodynamic conductance between the soil surface to the the canopy, 2nd conductance term from choudhury et al. 1988   
    ! based on an exponential decrease of wind speed with height
    ALPHA = 2
    Z0HT2 = 0.01

    ! Assuming uniform vegetation, the aerodynamic conductivity at the top of the canopy KH,
    ! and following Van de Griend 1989
    KH = ALPHA1 * VONKARMAN * WINDSTAR * (TREEH - ZPD2) 

    ! Aerodynamic conductance soil-air below canopy according to Chourdhury et al., 1988
    GBCANMS2 = ALPHA * KH / ( TREEH * exp(ALPHA) * (exp(-ALPHA * Z0HT2/TREEH)  -  exp(-ALPHA * (ZPD2+Z0) / TREEH) ) )

     ! 2nd alternative to GBCANMS1
      ! Atenuation coefficient, from Lafleur & Roux 1989
!      COAT = 2.6 * TOTLAI**0.36
!        IF (COAT.LT.1) COAT = 1
!        IF (COAT.GT.3) COAT = 3

      ! Intermediate values
!        WINDSTAR = WIND * VONKARMAN / log((ZSTAR-ZPD2)/Z0HT2)
!        Z0H = Z0*exp(-6.27*VONKARMAN*(WINDSTAR**(1/3)))
!      
!        GBCANMS2 = 1/  ( log((ZSTAR-ZPD2)/Z0)/(WIND*VONKARMAN**2) * (log((ZSTAR-ZPD2)/(TREEH-ZPD2)) +  &
!               (TREEH/(COAT*(TREEH-ZPD2)))*   (exp(COAT*(1-(ZPD2+Z0H)/TREEH))- 1)))     
!    GBCANMS1=0.02
!    GBCANMS2=0.01

    RETURN

END SUBROUTINE GBCANMS


!**********************************************************************
SUBROUTINE CALCWBIOM(IDAY,HT,DIAM,COEFFT,EXPONT,WINTERC,WBIOM,WBINC)
! Calculate the woody biomass (kg DW) on the given day from the height
! (m) and diameter (m). Also calculate the increment in woody biomass
! since previous day (g DW). Needed to calculate woody respiration.
!**********************************************************************
    IMPLICIT NONE
    INTEGER IDAY
    REAL PREVWBIOM,HT,DIAM,COEFFT,EXPONT,WINTERC,WBIOM,WBINC
    
    PREVWBIOM = WBIOM
    WBIOM = COEFFT*HT*(DIAM**EXPONT) + WINTERC
    IF (IDAY.EQ.0) PREVWBIOM = WBIOM
    WBINC = (WBIOM - PREVWBIOM)*1E3

    RETURN
END SUBROUTINE CALCWBIOM

!**********************************************************************
SUBROUTINE CALCFBIOM(IDAY,NOSLADATES,FOLLAY,SLA,PROP,NOLAY,NOAGEP, &
                        FBIOM,FBINC)
! Calculate foliage biomass from SLA and leaf area - done in layers.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER NOSLADATES,IDAY,NOLAY,NOAGEP,I,J
    REAL FOLLAY(MAXLAY)
    REAL SLA(MAXLAY,MAXC)
    REAL PROP(MAXC)
    REAL FBIOM,FBINC,PREVFBIOM

    IF (NOSLADATES.GT.0) THEN
        PREVFBIOM = FBIOM
        FBIOM = 0.0
        DO I = 1,NOLAY
            DO J = 1,NOAGEP
                FBIOM = FBIOM + FOLLAY(I)*PROP(J)/SLA(I,J)
            END DO
        END DO    
        IF (IDAY.EQ.0) PREVFBIOM = FBIOM
        FBINC = (FBIOM - PREVFBIOM)*1E3
    ELSE
        FBIOM = 0.0
        FBINC = 0.0
    END IF

    RETURN
END SUBROUTINE CALCFBIOM


!**********************************************************************
REAL FUNCTION CALCRMW(MODELRW,COLLA,COLLK,STEMSDW,DIAM,HT,STEMFORM,RMWAREA,WBIOM,RMW)
! Calculate stem respiration rate per unit biomass if necessary
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    INTEGER MODELRW
    REAL COLLA,COLLK,STEMSDW,DIAM,HT,STEMFORM,RMWAREA,WBIOM,RMW
    REAL STEMAREA

    IF (MODELRW.EQ.1) THEN
        RMW = COLLA*EXP(COLLK*DIAM*100.0)*STEMSDW
    ELSE IF (MODELRW.EQ.2) THEN
        STEMAREA = STEMFORM*PI*(DIAM**2)*HT
        RMW = RMWAREA*STEMAREA/WBIOM
    END IF

    CALCRMW = RMW
    RETURN
END FUNCTION CALCRMW


!**********************************************************************
REAL FUNCTION GRESP(BINC,EFFY)
! Calculate the growth respiration. Use increment in biomass
! (g DW tree-1 d-1) and the efficiency of growth (g g-1 C).
! Returns a value in mol CO2 tree-1 d-1.
!**********************************************************************

    USE maestcom
    IMPLICIT NONE
    REAL BINC,EFFY

    IF (BINC.GT.0.0) THEN
        GRESP = BINC * EFFY / GCPERMOL * CPERDW
    ELSE
        GRESP = 0.0
    ENDIF

    RETURN
END FUNCTION GRESP

!**********************************************************************
REAL FUNCTION CALCFSOIL(WSOILMETHOD,SOILMOISTURE,SOILDATA,SMD1,SMD2,WC1,WC2,SWPEXP)
! Calculate the effect of soil water content on stomatal conductance
! Two alternative forms:
! Granier & Loustau 1994 Fs = 1-a exp(b SMD) where SMD is soil moisture deficit, dimnless
! Otherwise negative exponential function of soil water potential
!**********************************************************************

    USE maestcom
    USE metcom
    IMPLICIT NONE
    INTEGER SOILDATA, WSOILMETHOD, SETFSOIL, IOERROR
    REAL SMD1,SMD2,WC1,WC2,SWPEXP,SOILMOISTURE,SOILMOIST

    CALCFSOIL = 1.0
    SETFSOIL = 0

    ! Using Emax method, set FSOIL to unity, re-calculate in PHOTOSYN subroutine.
    IF(WSOILMETHOD .EQ. 1)THEN
        CALCFSOIL = 1.0
        SETFSOIL = 1
    ENDIF

    ! Use volumetric water content, modelled or measured.
    ! Soil water is volumetric content (RAD 2008), use threshold response.
    IF(WSOILMETHOD .EQ. 2)THEN
        IF(SOILDATA .EQ. CONTENT .OR. SOILDATA .EQ. SIMULATED) THEN
            IF((WC1+WC2) .LT. 1E-09)THEN
                CALL SUBERROR('Error: Need to set WC1 and WC2.', IFATAL,IOERROR)
            ENDIF
            IF(WC1 .GT. WC2)THEN
                CALL SUBERROR('Error: WC1 needs to be smaller than WC2.',IFATAL,IOERROR)
            END IF
            CALCFSOIL = -WC1/(WC2-WC1) + SOILMOISTURE/(WC2-WC1)       
            IF(CALCFSOIL.GT.1.)CALCFSOIL = 1.                         
            IF(CALCFSOIL.LT.0.)CALCFSOIL = 0.                         
            SETFSOIL = 1
        END IF
    END IF

    ! Exponential relationship with potential: parameter = SWPEXP
    IF (SOILDATA.EQ.POTENTIAL .AND. SWPEXP.GT.0.0) THEN
        CALCFSOIL = EXP(SWPEXP*SOILMOISTURE)
        SETFSOIL = 1
    END IF

    ! Exponential relationship with deficit: params SMD1,SMD2
    IF(WSOILMETHOD .EQ. 4 .OR. SOILDATA .EQ. DEFICIT) THEN
    ! See ASSIGNSOILWATER, if this method used, soil water is
    ! already converted to deficit there.
    ! Data is already deficit:
        IF (SMD1.GT.0.0) THEN
            CALCFSOIL = 1.0 - SMD1*EXP(SMD2*SOILMOISTURE)
            ! Linear decline with increasing deficit: pUT SMD1 < 0, param SMD2
            SETFSOIL = 1
        ELSE IF (SMD2.GT.0.0) THEN
            IF ((1.0-SOILMOIST).LT.SMD2) THEN
                CALCFSOIL = (1.0-SOILMOIST)/SMD2
                SETFSOIL = 1
            END IF
        END IF
    ENDIF
    IF (CALCFSOIL.LT.0.0) CALCFSOIL = 0.0
    IF(SETFSOIL.EQ.0)THEN
        CALL SUBERROR('Error calculating FSOIL in CALCFSOIL. Check parameters.',IFATAL,IOERROR)
    END IF

    RETURN
END FUNCTION CALCFSOIL



!**********************************************************************
REAL FUNCTION FPSIL(PSIL,SF,PSIV)

! Tuzet et al. 2003 leaf water potential function.
!**********************************************************************

IMPLICIT NONE
REAL PSIL,SF,PSIV

FPSIL = (1 + EXP(SF*PSIV) )/ (1 + EXP(SF*(PSIV-PSIL)))

END FUNCTION FPSIL



!**********************************************************************

SUBROUTINE PSILFIND(RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25,&
                    IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0, &
                    Q10F,K10F,RTEMP,DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,    &
                    SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,GSMIN,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,   &
                    VPARA,VPARB,VPARC,VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT,TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP,  &
                    WEIGHTEDSWP,HMSHAPE,PSILIN,ETDEFICIT,ETEST,iday,ihour,G02,G12,NEWTUZET,EV,drycan,CANOPY_STORE_I) !glm canopy evap
                    
!**********************************************************************
        USE maestcom
        IMPLICIT NONE

        INTEGER MODELGS,SOILDATA,WSOILMETHOD,ITER
        INTEGER IECO,ITERMAX,NSIDES,VFUN,iday,ihour
        REAL JMAX25,I0,LHV,PSIL,K10F
        REAL VPARA,VPARB,VPARC
        REAL MINLEAFWP,TOTSOILRES,PLANTK, ETDEFICIT
        REAL TLEAF,TAIR,DLEAF,VPD,VMLEAF,VMFD,RHLEAF,RH,CS,CA
        REAL SLOPE,GRADN,RDFIPT,TUIPT,TDIPT,GBHU,PRESS,WIND
        REAL WLEAF,PAR,TMOVE,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC
        REAL DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP
        REAL TBELOW,GSREF,D0,VK1,VK2,VPD1,VPD2,VMDF0
        REAL GSJA,GSJB,T0,TREF,TMAX,SOILMOISTURE,EMAXLEAF
        REAL SMD1,SMD2,WC1,WC2,SWPEXP,FSOIL,G0,D0L,GAMMA,G1
        REAL GSC,ALEAF,RD,WEIGHTEDSWP,GBHF,GBH,GH,VMFD0,GBV,GSV,GV
        REAL ET,RNET,GBC,TDIFF,TLEAF1,FHEAT,ETEST,SF,PSIV,HMSHAPE
        REAL PSILIN,T1,T2,XACC,GK,TMP,VPDMIN
        REAL EXTRAPARS(EXTRAPARDIM)
        REAL G02,G12,GSMIN
        REAL EV,drycan,CANOPY_STORE_I !glm canopy evap
        INTEGER NEWTUZET
        INTEGER EXTRAINT(10)
        REAL, EXTERNAL :: ZBRENT
        REAL, EXTERNAL :: PSILOBJFUN
        
        
        ! INIT
        PSILIN = -0.1
        PSIL = -0.1
          
        EXTRAINT(1) = IECO
        EXTRAINT(2) = MODELGS
        EXTRAINT(3) = WSOILMETHOD
        EXTRAINT(4) = NSIDES
        EXTRAINT(5) = ITERMAX
        EXTRAINT(6) = VFUN
        EXTRAINT(7) = IDAY  !modification
        EXTRAINT(8) = IHOUR !modification        

        EXTRAPARS(1) = RDFIPT
        EXTRAPARS(2) = TUIPT
        EXTRAPARS(3) = TDIPT
        EXTRAPARS(4) = RNET
        EXTRAPARS(5) = WIND
        EXTRAPARS(6) = PAR
        EXTRAPARS(7) = TAIR
        EXTRAPARS(8) = TMOVE
        EXTRAPARS(9) = CA
        EXTRAPARS(10) = RH
        EXTRAPARS(11) = VPD
        EXTRAPARS(12) = VMFD        
        EXTRAPARS(13) = PRESS
        EXTRAPARS(14) = JMAX25
        EXTRAPARS(15) = EAVJ
        EXTRAPARS(16) = EDVJ
        EXTRAPARS(17) = DELSJ       
        EXTRAPARS(18) = VCMAX25
        EXTRAPARS(19) = EAVC
        EXTRAPARS(20) = EDVC
        EXTRAPARS(21) = DELSC
        EXTRAPARS(22) = TVJUP
        EXTRAPARS(23) = TVJDN
        EXTRAPARS(24) = THETA
        EXTRAPARS(25) = AJQ
        EXTRAPARS(26) = RD0
        EXTRAPARS(27) = Q10F
        EXTRAPARS(28) = K10F
        EXTRAPARS(29) = RTEMP
        EXTRAPARS(30) = DAYRESP
        EXTRAPARS(31) = TBELOW
        EXTRAPARS(32) = TOTSOILRES
        EXTRAPARS(33) = SOILMOISTURE
        EXTRAPARS(34) = SMD1
        EXTRAPARS(35) = SMD2       
        EXTRAPARS(36) = WC1        
        EXTRAPARS(37) = WC2        
        EXTRAPARS(38) = SOILDATA        
        EXTRAPARS(39) = SWPEXP        
        EXTRAPARS(40) = FSOIL        
        EXTRAPARS(41) = G0        
        EXTRAPARS(42) = D0L        
        EXTRAPARS(43) = GAMMA        
        EXTRAPARS(44) = G1        
        EXTRAPARS(45) = WLEAF        
        EXTRAPARS(46) = SF        
        EXTRAPARS(47) = PSIV        
        EXTRAPARS(48) = PLANTK        
        EXTRAPARS(49) = WEIGHTEDSWP        
        EXTRAPARS(50) = HMSHAPE        
        EXTRAPARS(51) = VPARA
        EXTRAPARS(52) = VPARB
        EXTRAPARS(53) = VPARC
        EXTRAPARS(54) = VPDMIN
        EXTRAPARS(55) = GK
        EXTRAPARS(56) = NEWTUZET
        EXTRAPARS(57) = G02
        EXTRAPARS(58) = G12
        EXTRAPARS(59) = GSMIN
        EXTRAPARS(60) = ETDEFICIT
        EXTRAPARS(61) = drycan
        EXTRAPARS(62) = CANOPY_STORE_I
        
        ! Set bounds for root-finding
        T1 = -100.0
        T2 = 0.0

        ! Error tolerance.
        XACC = 1E-03
        
        PSILIN = ZBRENT(PSILOBJFUN,T1,T2,XACC,EXTRAPARS,EXTRAINT)

END


!**********************************************************************
REAL FUNCTION PSILOBJFUN(PSILIN, EXTRAPARS, EXTRAINT)

! New PSIL finder:
! Make wrapper function that takes PSIL as first argument, EXTRAPARS an array
! second argument, and returns the squared difference in PSILIN and PSIL
! Then : PSILOPT = PSIFLFINDER(PSIL, EXTRAPARS)
! Load and unload parameters from EXTRAPARS...
!**********************************************************************
        USE maestcom
        IMPLICIT NONE

        INTEGER MODELGS,SOILDATA,WSOILMETHOD,ITER
        INTEGER IECO,ITERMAX,NSIDES,VFUN
        REAL JMAX25,I0,LHV,MINROOTWP,KTOT,PSIL,K10F
        REAL TLEAF,TAIR,DLEAF,VPD,VMLEAF,VMFD,RHLEAF,RH,CS,CA
        REAL SLOPE,GRADN,RDFIPT,TUIPT,TDIPT,GBHU,PRESS,WIND
        REAL WLEAF,PAR,TMOVE,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC
        REAL DELSC,TVJUP,TVJDN,THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP
        REAL TBELOW,GSREF,GSMIN,D0,VK1,VK2,VPD1,VPD2,VMDF0
        REAL GSJA,GSJB,T0,TREF,TMAX,SOILMOISTURE,EMAXLEAF
        REAL SMD1,SMD2,WC1,WC2,SWPEXP,FSOIL,G0,D0L,GAMMA,G1
        REAL GSC,ALEAF,RD,WEIGHTEDSWP,GBHF,GBH,GH,VMFD0,GBV,GSV,GV
        REAL ET,RNET,GBC,TDIFF,TLEAF1,FHEAT,ETEST,SF,PSIV,HMSHAPE
        REAL PSILIN,TOTSOILRES,PLANTK,MINLEAFWP,CI,GK,ETDEFICIT
        REAL VPARA,VPARB,VPARC,VPDMIN,GNIGHT
        REAL EV,drycan,CANOPY_STORE_I !glm canopy evap
        LOGICAL ISMAESPA,ISNIGHT
        integer iday,ihour, NEWTUZET
        REAL G02,G12
        REAL EXTRAPARS(EXTRAPARDIM)
        INTEGER EXTRAINT(10)
        LOGICAL EXTRALOGIC(10)

          IECO = EXTRAINT(1)
          MODELGS = EXTRAINT(2)
          WSOILMETHOD = EXTRAINT(3)
          NSIDES = EXTRAINT(4)
          ITERMAX = EXTRAINT(5)
          VFUN = EXTRAINT(6)
          IDAY = EXTRAINT(7)
          IHOUR = EXTRAINT(8)
          
          RDFIPT =EXTRAPARS(1)
          TUIPT =EXTRAPARS(2)
          TDIPT =EXTRAPARS(3)
          RNET =EXTRAPARS(4)
          WIND =EXTRAPARS(5)
          PAR =EXTRAPARS(6)
          TAIR =EXTRAPARS(7)
          TMOVE =EXTRAPARS(8)
          CA =EXTRAPARS(9)
          RH =EXTRAPARS(10)
          VPD =EXTRAPARS(11)
          VMFD = EXTRAPARS(12)
          PRESS =EXTRAPARS(13)
          JMAX25 =EXTRAPARS(14)
          EAVJ =EXTRAPARS(15)
          EDVJ =EXTRAPARS(16)
          DELSJ =     EXTRAPARS(17) 
          VCMAX25 =EXTRAPARS(18)
          EAVC =EXTRAPARS(19) 
          EDVC =EXTRAPARS(20)
          DELSC =EXTRAPARS(21)
          TVJUP =EXTRAPARS(22)
          TVJDN =EXTRAPARS(23)
          THETA =EXTRAPARS(24)
          AJQ =EXTRAPARS(25)
          RD0 =EXTRAPARS(26)
          Q10F =EXTRAPARS(27)
          K10F =EXTRAPARS(28)
          RTEMP =EXTRAPARS(29)
          DAYRESP =EXTRAPARS(30)
          TBELOW =EXTRAPARS(31)
          TOTSOILRES =EXTRAPARS(32)
          SOILMOISTURE =EXTRAPARS(33)
          SMD1 =EXTRAPARS(34)
          SMD2 =   EXTRAPARS(35)   
          WC1 = EXTRAPARS(36)       
          WC2 =EXTRAPARS(37)       
          SOILDATA =EXTRAPARS(38)
          SWPEXP = EXTRAPARS(39) 
          FSOIL =EXTRAPARS(40)     
          G0 = EXTRAPARS(41)  
          D0L = EXTRAPARS(42)      
          GAMMA = EXTRAPARS(43)       
          G1 = EXTRAPARS(44)   
          WLEAF =  EXTRAPARS(45)      
          SF = EXTRAPARS(46)     
          PSIV = EXTRAPARS(47)      
          PLANTK = EXTRAPARS(48)      
          WEIGHTEDSWP =EXTRAPARS(49)
          HMSHAPE = EXTRAPARS(50) 
          VPARA = EXTRAPARS(51)
          VPARB = EXTRAPARS(52)
          VPARC = EXTRAPARS(53)
          VPDMIN = EXTRAPARS(54)
          GK = EXTRAPARS(55)
          NEWTUZET = EXTRAPARS(56)
          G02 = EXTRAPARS(57)
          G12 = EXTRAPARS(58)
          GSMIN = EXTRAPARS(59)
          ETDEFICIT = EXTRAPARS(60)
          drycan  = EXTRAPARS(61)
          CANOPY_STORE_I  = EXTRAPARS(62)
          
          ISMAESPA = .TRUE.
          ISNIGHT = .FALSE.

        MINLEAFWP = 0  ! Not used in tuzet.
        CI = 0
        KTOT = 0  ! output, not needed.
        GNIGHT = 0.01  ! not needed as input because psilfind not executed during the night
        
        CALL PSTRANSP(iday,ihour,RDFIPT,TUIPT,TDIPT,RNET,WIND,PAR,TAIR,TMOVE,CA,RH,VPD,VMFD,PRESS,JMAX25, &
             IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN,THETA,AJQ,RD0, &
             Q10F,K10F,RTEMP,DAYRESP,TBELOW,MODELGS,WSOILMETHOD,EMAXLEAF,SOILMOISTURE,    &
             SMD1,SMD2,WC1,WC2,SOILDATA,SWPEXP,FSOIL,GSMIN,GNIGHT,G0,D0L,GAMMA,VPDMIN,G1,GK,WLEAF,NSIDES,   &
             VPARA,VPARB,VPARC,VFUN,SF,PSIV,ITERMAX,GSC,ALEAF,RD,ET,FHEAT, &
             TLEAF,GBH,PLANTK,TOTSOILRES,MINLEAFWP, WEIGHTEDSWP,KTOT,HMSHAPE,PSILIN,PSIL,ETDEFICIT,ETEST,CI,ISMAESPA,ISNIGHT,&
             G02,G12,NEWTUZET,EV,drycan,CANOPY_STORE_I) !glm canopy evap

        PSILOBJFUN = (PSILIN - PSIL)**2
        !print*,'PSILOBJFUN',ihour,PSILIN,PSIL,PSILOBJFUN

END
