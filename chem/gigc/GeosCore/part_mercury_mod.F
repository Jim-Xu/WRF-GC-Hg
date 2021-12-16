!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mercury_mod.F
!
! !DESCRIPTION: Module MERCURY\_MOD contains variables and routines for the 
!  GEOS-CHEM mercury simulation.  Many choices of reaction mechanism and 
!  model processes can be selected with logical switches located in 
!  INIT\_MERCURY.
!\\
!\\
! !INTERFACE:
!
      MODULE PART_MERCURY_MOD
!
! !USES:
!
      USE PhysConstants           ! Physical constants
      USE PRECISION_MOD           ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: PARTITIONHG
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: partitionhg
!
! !DESCRIPTION: Subroutine PARTITIONHG2 splits Hg(II) into gas and aerosol 
!  portions  according to the thermodynamic equilibrium determined by 
!  temperature and aerosol surface area.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PARTITIONHG( am_I_Root, Input_Opt, State_Chm, 
     &                         State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : SAFE_DIV  , GEOS_CHEM_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ind_ 
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  26-Jun-2019 - Y. Zhang    _ Modify to GCHP
!  12-Jul-2010 - H. Amos     - Add option to partition Hg2 according to Fg/Fp
!  02-Jan-2011 - H. Amos     - Add/clean up comments
!  04 Jan 2012 - H. Amos     - modify algorithms to reflect the fact that 
!                              anthropogenic Hg(p) is now emitted as Hg(II) 
!                              (i.e. it's no longer considered refractory). 
!  28-Mar-2013 - C. Holmes   - Since we now assume that all HgP and Hg2 
!                              actively partition between gas and aerosol, 
!                              we no longer need to 'reverse' partition. 
!  03 Jun 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  18 Nov 2014 - M. Yannetti - Took out OCEAN_MERCURY_MOD calling Fp
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER          :: I, J, L, N, IND
      REAL(fp)         :: FGas, Hg2TOT
      REAL(fp)         :: PM25, R, Hg2, HgP
      INTEGER          :: IND_SO4, IND_NH4, IND_NIT, IND_BCPI, IND_BCPO
      INTEGER          :: IND_OCPI, IND_OCPO, IND_DST1, IND_Hg2, IND_HgP

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! PARTITIONHG2 begins here!
      !=================================================================

      ! Assume success
      RC  =  GC_SUCCESS 

      ! Point to the chemical species array [kg/kg dry]
      Spc => State_Chm%Species

      ! Get index of tracers
      IND_SO4  = IND_('SO4')
      IND_NH4  = IND_('NH4')
      IND_NIT  = IND_('NIT')
      IND_BCPI = IND_('BCPI')
      IND_BCPO = IND_('BCPO')
      IND_OCPI = IND_('OCPI')
      IND_OCPO = IND_('OCPO')
      IND_DST1 = IND_('DST1')
      IND_Hg2  = IND_('Hg2')
      IND_HgP  = IND_('HgP')

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, IND, PM25, R, FGas, Hg2TOT )
!$OMP+PRIVATE( Hg2, HgP                               )
      DO L=1, LLPAR
      DO J=1, JJPAR
      DO I=1, IIPAR
         
         ! Get PM2.5 concentrations, ug/m3
         ! Initialize
         PM25 =  Spc(I,J,L,IND_SO4) +   ! kg/kg dry
     &           Spc(I,J,L,IND_NH4) + 
     &           Spc(I,J,L,IND_NIT) + 
     &           Spc(I,J,L,IND_BCPI) + 
     &           Spc(I,J,L,IND_BCPO) + 
     &           Spc(I,J,L,IND_OCPI) + 
     &           Spc(I,J,L,IND_OCPO) + 
     &           Spc(I,J,L,IND_DST1)

         ! transfer from kg/kg dry to ug/m3
         IF ( State_Chm%Spc_Units == 'kg/kg dry') THEN
            PM25 = 1.0e9_fp * PM25 * State_Met%AIRDEN(I,J,L)
         ELSEIF ( State_Chm%Spc_Units == 'kg') THEN
            PM25 = 1.0e9_fp * PM25 / State_Met%AIRVOL(I,J,L)
         ELSE
            write(*,*) State_Chm%Spc_Units
         ENDIF

         ! Get Hg2 and HgP
         Hg2 = Spc(I,J,L,IND_Hg2)
         HgP = Spc(I,J,L,IND_HgP)

         ! Calculate R = HgP_ads/RGM (i.e. the ratio of Hg2
         ! adsorbed onto aerosol to Hg2 in the  gas phase) 
         R = PM25*
     &     ( 10e+0_fpp**( ( 2.5e+3_fpp / State_Met%T(I,J,L)) 
     &     - 10e+0_fpp ))

         ! Fraction of Hg(II) in the gas phase (unitless)
         FGas = 1e+0_fpp / (R + 1e+0_fpp)   

         ! Total Hg(II) (gas +aerosol)
         HG2TOT = Hg2 + HgP

         ! Gas portion
         Spc(I,J,L,IND_Hg2) = Hg2TOT * Fgas 
         Spc(I,J,L,IND_HgP) = Hg2TOT * (1e+0_fp - Fgas)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Free pointer memory
      Spc => NULL()

      ! Return to calling program
      END SUBROUTINE PARTITIONHG
!EOC
      END MODULE PART_MERCURY_MOD
