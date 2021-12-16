!WRF:MODEL_LAYER:INITIALIZATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   WRF-GCHP
!   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
!
!   WRF & GCHP are (c) their original authors.
!   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
!
!   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng, 2018-01
!   Peking University, School of Physics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
!
!   This Chemical Interface (chem/) is written after comprehensive study of
!   the original chem_driver.f & chemics_init.f from WRF-Chem v3.6.1
!   which is (c) their respective authors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   MAJOR CHANGES:
!   - Jan 2018. H.P. Lin - Initial Version
!   - Apr 2018. H.P. Lin - Restructure chemics_init to use WRF memorder XZY & get IJK from module_domain,
!                          as there seems to be stack incompatibility between chemics_init & start_em.
!   - May 2018. H.P. Lin - Allow ICs to be passed in from WRF, while passing b/g values from GC if no ICs.
!   - Jun 2018. H.P. Lin - Now use GIGC_Stateful_Mod
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chem_init(id, chem, emis_ant, scalar, dt, bioemdt, photdt, chemdt, stepbioe, &
                    stepphot, stepchem, stepfirepl, plumerisefire_frq, z_at_w, xlat, xlong, &
                    g, aerwrf, config_flags, grid, alt, t, p, CONVFAC, ttday, tcosz, julday, gmt, &
                    tauaer1, tauaer2, tauaer3, tauaer4, &
                    gaer1, gaer2, gaer3, gaer4, &
                    waer1, waer2, waer3, waer4, &
                    l2aer, l3aer, l4aer, l5aer, l6aer, l7aer, &
                    extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
                    extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
                    extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
                    extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
                    tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
                    tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
                    tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
                    tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
                    dgnum4d, dgnumwet4d, dgnum_a1, dgnum_a2, dgnum_a3, &
                    dgnumwet_a1, dgnumwet_a2, dgnumwet_a3, &
                    pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
                    last_chem_time_year, last_chem_time_month, &
                    last_chem_time_day, last_chem_time_hour, &
                    last_chem_time_minute, last_chem_time_second, &
                    chem_in_opt, kemit, num_vert_mix, &
                    oldids, oldide, oldjds, oldjde, oldkds, oldkde, &
                    oldims, oldime, oldjms, oldjme, oldkms, oldkme, &
                    oldits, oldite, oldjts, oldjte, oldkts, oldkte)

      ! WRF Core Dependencies
      use module_domain
      use module_configure
      use module_state_description
      use module_model_constants, only: DEGRAD
      use module_dm

      ! WRF Physics Dependencies
      ! Just maintaining old code to not break CAM Physics Support.
      use module_cam_support, only: gas_pcnst_modal_aero, gas_pcnst_modal_aero_pos

      ! WRF Chemistry Dependencies
      use module_tropopause, only: tropopause_init
      use module_input_chem_data, only: get_last_gas, last_chem_time, setup_gasprofile_maps

      !! GIGC: Grid Independent GEOS-Chem
      ! GIGC_Chunk_Mod provides the GEOS-Chem Column Model Code
      use GIGC_Chunk_Mod

      ! HEMCO: Harvard-NASA Emissions Component
      use HCO_TYPES_MOD, only: ConfigObj

      ! Input_Opt_Mod provides the Input_Opt% type (OptInput) struct
      use Input_Opt_Mod, only: OptInput
      use State_Chm_Mod, only: ChmState
      use State_Met_Mod, only: MetState
      use State_Diag_Mod, only: DgnState

      ! Use our conversion module for GEOS-Chem code
      use GIGC_Convert_State_Mod, only: GIGC_Set_WRF, GIGC_IdxSetup

      ! Use our stateful module for Grid-Independent GEOS-Chem States
      use GIGC_Stateful_Mod, only: GIGC_State_Get_Status

      implicit none

      ! provides information on WRF DMPAR parallel (MPI) process information
      ! namely whether this is the master process, using wrf_dm_on_monitor()
      logical, external :: wrf_dm_on_monitor

      real, intent(in) :: bioemdt, photdt, chemdt, dt, gmt
      integer, intent(in) :: plumerisefire_frq
      integer, intent(in) :: chem_in_opt
      integer, intent(inout) :: num_vert_mix
      integer, intent(in) :: id, julday, kemit, &
                             last_chem_time_year, &
                             last_chem_time_month, &
                             last_chem_time_day, &
                             last_chem_time_hour, &
                             last_chem_time_minute, &
                             last_chem_time_second, &
                             oldids, oldide, oldjds, oldjde, oldkds, oldkde, &
                             oldims, oldime, oldjms, oldjme, oldkms, oldkme, &
                             oldits, oldite, oldjts, oldjte, oldkts, oldkte

      integer :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte
                 
      integer, intent(out) :: stepbioe, stepphot, stepchem, stepfirepl
      type(grid_config_rec_type), intent(in) :: config_flags
      type(domain), intent(inout) :: grid

      ! MEMORDER XZY: %s-em31 is I, %s-em32 is K, %s-em33 is J
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), &
         intent(INOUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
         tauaer1, tauaer2, tauaer3, tauaer4, &
         extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
         extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
         extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
         extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
         tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
         tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
         tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
         tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
         gaer1, gaer2, gaer3, gaer4, &
         waer1, waer2, waer3, waer4

      !-- Arrays needed to output dgnum when restarting
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 3), intent(INOUT) :: dgnum4d, dgnumwet4d
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(INOUT) :: dgnum_a1, dgnum_a2, dgnum_a3, dgnumwet_a1, dgnumwet_a2, dgnumwet_a3
      !-- end dgnum restart arrays

      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 1:4), &
         intent(INOUT) :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(IN) :: z_at_w, t, p, alt, convfac
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(INOUT) :: chem
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, 1) :: moist_dummy
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_emis_ant), intent(INOUT) :: emis_ant
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(INOUT) :: scalar
      real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(INOUT) :: aerwrf
      real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33), intent(INOUT) :: ttday, tcosz, xlat, xlong
      real, intent(IN) :: g

      ! Local Variables
      character*256 :: mminlu_loc
      character*256 :: message_txt
      type(WRFU_TimeInterval) :: tmpTimeInterval
      integer :: i, j, k, l, numgas, n, kk, nv, gigc_ptr

      !! GEOS-Chem Integration Variables
      ! GEOS_CHEM_RC (RC in GEOS-Chem speak) is the return result code for all GIGC Subroutines.
      ! 0 = Success, > 0 is error
      integer :: GEOS_CHEM_RC
      logical :: IsInit

      type(OptInput)       :: Input_Opt      ! Input Options (same for all domains)
      type(MetState)       :: State_Met      ! Meteorology state
      type(ChmState)       :: State_Chm      ! Chemistry state
      type(DgnState)       :: State_Diag     ! Diagnostics State
      type(ConfigObj), pointer :: HcoConfig


      real(4), allocatable :: lonCtr(:, :)  ! Lon Centers (rad)
      real(4), allocatable :: latCtr(:, :)  ! Lat Centers (rad)
      real(4), allocatable :: lonEdge(:, :) ! Lon Edges (rad)
      real(4), allocatable :: latEdge(:, :) ! Lat Edges (rad)
      integer              :: IM, JM
      integer              :: II, JJ       ! Loop variable
      integer              :: nymdB, nhmsB, nymdE, nhmsE

      ! From module_aerosols_sorgam / dep_simple, hplin 2-25-2018
      ! num_chem is from module_state_description
      logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)

      ! From module ctrans_grell, hplin 2-25-2018
      real                 :: epsilc

      ! WRF DM (MPI) Parallel Information - is master process?
      logical              :: Am_I_Root
      integer              :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm

      !!!!!!!!!!!!! %% EXECUTABLE SECTION %% !!!!!!!!!!!!!!
      call get_ijk_from_grid(grid, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte)

      ! Do not touch ide, jde
      ite = min(ite, ide - 1)
      jte = min(jte, jde - 1)

      ! Domain sizing
      IM = ite - its + 1
      JM = jte - jts + 1

      ! Prep GEOS_CHEM ReturnCode (RC)
      GEOS_CHEM_RC = 0

      if(wrf_dm_on_monitor()) then
         Am_I_Root = .true.
      else
         Am_I_Root = .false.
      endif

      call wrf_get_nproc(WRF_DM_NProc)
      call wrf_get_myproc(WRF_DM_MyProc)
      call wrf_get_dm_communicator(WRF_DM_Comm)

      ! Pass some HPC Information to Input_Opt...
      Input_Opt%HPC     = .true.
      Input_Opt%RootCPU = Am_I_Root
      Input_Opt%myCPU   = WRF_DM_MyProc
      Input_Opt%NPES    = WRF_DM_NProc
      Input_Opt%MPICOMM = WRF_DM_Comm

      Input_Opt%LPRT    = .true. ! Print debug information? (LD70)

      ! Have we ran this before...?
      call GIGC_State_Get_Status(Am_I_Root, grid%id, IsInit)
      if(IsInit) then
        call wrf_message("GIGC chemics_init: Domain already initialized, skipping.")
        return
      endif

      ! For parallelism (dirty) in chemics_init let's just pass dummies
      ! to lon/latCtr for now. lonCtr and latCtrs' correct values will be
      ! passed onto GIGC_CHUNK_MOD at a later time when doing Chunk_Run.
      !
      ! A value that is enough for any chunk-size should be specified here.
      ! TODO hplin 4/27/2018
      allocate(lonCtr (IM,   JM))
      allocate(latCtr (IM,   JM))
      allocate(lonEdge(IM+1, JM))
      allocate(latEdge(IM,   JM+1))
      ! write(6, *) "its, ite, jts, jte", its, ite, jts, jte
      ! write(6, *) "ims, ime, jms, jme", ims, ime, jms, jme
      ! write(6, *) "im, jm", (ite-its+1), (jte-jts+1)

      ! For every grid box populate with relevant info...
      ! FIXME this breaks multiple domain configuration - beware! (hplin 5/11/18)
      do j = jts, jte
         do i = its, ite
            ! Loop indexes in GC
            II = i - its + 1
            JJ = j - jts + 1

            ! Prepare lat/lon grid box center information in [rad]
            ! Grab lat lon directly from grid%XLAT, grid%XLONG [deg] (DEGRAD is piconst/180. in model_constants)
            lonCtr(II, JJ) = (grid%XLONG(i, j)) * DEGRAD
            latCtr(II, JJ) = (grid%XLAT (i, j)) * DEGRAD

            ! Grab grid%XLONG_U, grid%XLAT_V for lonEdge, latEdge
            lonEdge(II, JJ) = grid%XLONG_U(i, j) * DEGRAD
            if(i .eq. ite) then
              lonEdge(II+1, JJ) = grid%XLONG_U(i+1, j) * DEGRAD
            endif

            latEdge(II, JJ) = grid%XLAT_V(i, j) * DEGRAD
            if(j .eq. jte) then
              latEdge(II, JJ+1) = grid%XLAT_V(i, j+1) * DEGRAD
            endif
         enddo
      enddo

      ! Get time... YYYYMMDD
      nymdB = config_flags%start_year * 10000 + config_flags%start_month * 100 + config_flags%start_day
      nhmsB = config_flags%start_hour * 10000 + config_flags%start_minute * 100 + config_flags%start_second

      nymdE = config_flags%end_year * 10000 + config_flags%end_month * 100 + config_flags%end_day
      nhmsE = config_flags%end_hour * 10000 + config_flags%end_minute * 100 + config_flags%end_second


#ifdef CHEM_DBG_I
      call print_chem_species_index(config_flags%chem_opt)
#endif

      write(message_txt, *) "   This is domain ", grid%id
      call wrf_message("*****************************************************************")
      call wrf_message("  __          _______  ______       _____  _____ _    _ _____    ")
      call wrf_message("  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   ")
      call wrf_message("   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  ")
      call wrf_message("    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   ")
      call wrf_message("     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       ")
      call wrf_message("      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       ")
      call wrf_message("*****************************************************************")
      call wrf_message("         WRF-GC MODEL, VERSION 1811.11 - WRF BINDINGS            ")
      call wrf_message(" For errors, suggestions & feedback email linhaipeng@pku.edu.cn  ")
      call wrf_message("*****************************************************************")
      call wrf_message(message_txt)


      numgas = get_last_gas(config_flags%chem_opt)

      ! WRF-Chem used to perform configuration flags validation here.
      ! This can be done too in the future to prevent "gotchas"

      IF (config_flags%chem_opt == 0 .AND. config_flags%aer_ra_feedback .NE. 0) THEN
         call wrf_error_fatal("chemics_init: If no chemistry, aer_ra_feedback must be 0")
      ENDIF

      CALL nl_get_mminlu(1, mminlu_loc)

      IF (trim(mminlu_loc) /= 'USGS' .and. trim(mminlu_loc) /= 'MODIFIED_IGBP_MODIS_NOAH') THEN
         print *, mminlu_loc
         call wrf_error_fatal("chemics_init: Chemistry routines require USGS or MODIS_NOAH land use maps.")
      ELSE
         IF (trim(mminlu_loc) == 'USGS' .and. grid%num_land_cat <= 23) THEN
            call wrf_error_fatal("chemics_init: USGS land use map should have 24 or more categories.")
         ELSEIF (trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' .and. grid%num_land_cat <= 19) THEN
            call wrf_error_fatal("chemics_init: MODIS_NOAH land use map should have 20 or more categories.")
         ENDIF
      ENDIF

      ! Catch some unsupported options.
      ! In early development.
      if (config_flags%progn > 0) then
         call wrf_error_fatal("Pumpkin chemics_init: progn > 0 option is not supported by this Chemistry abstraction layer.")
      endif

      if (config_flags%have_bcs_upper) then
         call wrf_error_fatal("Pumpkin chemics_init: have_bcs_upper = .true. is not supported by this Chemistry abstraction layer.")
      endif

      ! Don't run if there isn't anything to do...
      if(config_flags%chem_opt == 0) then
        CALL wrf_debug(15, 'Pumpkin chemics_init: **** NO CHEMISTRY IS USED IN THIS WRF RUN ****')
        return
      endif

      ! Set the num_vert_mix variable if using ACM
      num_vert_mix = 0
      IF (config_flags%bl_pbl_physics == ACMPBLSCHEME) THEN
         mix_select: select case(config_flags%chem_opt)
            case (radm2sorg_aq, radm2sorg_aqchem, racmsorg_aq)
               ! if you have aqueous chemistry, then the num_vert_mix needs to be # of gas species & not # of chemistry_species.
               num_vert_mix = numgas
            case default
               num_vert_mix = num_chem
         end select mix_select

         if (num_vert_mix .gt. config_flags%ndepvel) then
            write (message_txt, '(A30,2(I8,2x))') 'chem_init: num_vert_mix and ndepvel ', num_vert_mix, config_flags%ndepvel
            call wrf_message(trim(message_txt))
            call wrf_error_fatal("Pumpkin chemics_init: FATAL - num_vert_mix > ndepvel ")
         endif
      endif

      ! Calculate "timesteps" for out intent.
      stepbioe = nint(bioemdt*60./dt)
      stepphot = nint(photdt*60./dt)
      stepchem = nint(chemdt*60./dt)
      stepfirepl = nint(plumerisefire_frq*60/dt)
      stepbioe = max(stepbioe, 1)
      stepphot = max(stepphot, 1)
      stepchem = max(stepchem, 1)
      stepfirepl = max(stepfirepl, 1)
      call wrf_debug(15, 'Pumpkin chemics_init: Starting chemistry init')

      ! Setup the cross reference mappings between the default profiles and
      ! the gas mechanism species (wig, 2-May-2007)
      ! Sets up the cross reference mapping indices and fractional
      ! apportionment of the default species profiles for use with
      ! ICs and BCs. (module_input_chem_data)
      call setup_gasprofile_maps(config_flags%chem_opt, numgas)

      ! Aerosol Information (is_aerosol array) allocation.
      if (.not. allocated(is_aerosol)) then
         allocate (is_aerosol(num_chem))
      else
         if (size(is_aerosol) /= num_chem) &
            call wrf_error_fatal("The number of chemistry species has changed between nests. Use same chem_opt for all domains.")
      end if

      ! The top level in WRF does not contain data and should NOT be passed to GEOS-Chem
      ! nor the stateful conversion routines (8/7/18, hplin)
      kte = kte - 1

      write(message_txt, *) "Pumpkin chemics_init: its, ite, jts, jte, kts, kte: ", its, " ", ite, " ", jts, " ", jte, " ", kts, " ", kte, "domain", grid%id
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chemics_init: ims, ime, jms, jme, kms, kme: ", ims, " ", ime, " ", jms, " ", jme, " ", kms, " ", kme, "domain", grid%id
      call wrf_debug(15, message_txt)
     
      write(message_txt, *) "Pumpkin chemics_init: ids, ide, jds, jde, kds, kde: ", ids, " ", ide, " ", jds, " ", jde, " ", kds, " ", kde, "domain", grid%id
      call wrf_debug(15, message_txt)

      ! Prepare GIGC
      ! Call the GEOS-Chem Column Code for Initialization 
      ! GIGC_Chunk_Init will call GIGC_State_Init to initialize the respective stateful variable values.
      ! Also call GIGC_State_Boot to read the configuration variables as necessary ...
      call GIGC_Chunk_Init(am_I_Root = Am_I_Root,           & ! Are we on the root PET?   
                           I_LO      = 1,                   & ! Min lon index on this PET
                           J_LO      = 1,                   & ! Min lat index on this PET
                           I_HI      = ite-its+1,           & ! Max lon index on this PET - for IIPAR
                           J_HI      = jte-jts+1,           & ! Max lat index on this PET - for JJPAR
                           IM        = ite-its+1,           & ! # lons   on this PET
                           JM        = jte-jts+1,           & ! # lats   on this PET
                           LM        = kte-kts+1,           & ! # levels on this PET - do not pass top level (see deduction above)
                           ID        = grid%id,             & ! Domain identifier, this PET
                           IM_WORLD  = ite-its+1,           & ! # lons   in global grid
                           JM_WORLD  = jte-jts+1,           & ! # lats   in global grid
                           LM_WORLD  = kte-kts+1,           & ! # levels in global grid - do not pass top level (see deduction above)
                           nymdB     = nymdB,               & ! YYYYMMDD @ start of run
                           nhmsB     = nhmsB,               & ! hhmmss   @ start of run
                           nymdE     = nymdE,               & ! YYYYMMDD @ end of run
                           nhmsE     = nhmsE,               & ! hhmmss   @ end of run
                           tsChem    = chemdt*60.,          & ! Chemical timestep [s], REAL
                           tsDyn     = chemdt*60.,          & ! Dynamic  timestep [s], REAL
                           lonCtr    = lonCtr,              & ! Lon centers [radians]
                           latCtr    = latCtr,              & ! Lat centers [radians]
                           lonEdge   = lonEdge,             & ! Lon edges [radians]
                           latEdge   = latEdge,             & ! Lat edges [radians]
                           myPET     = WRF_DM_MyProc,       & ! Local PET
                           Input_Opt = Input_Opt,           & ! Input Options obj
                           State_Chm = State_Chm,           & ! Chemistry State obj
                           State_Diag= State_Diag,          & ! Diagnostics State obj
                           State_Met = State_Met,           & ! Meteorology State obj
                           HcoConfig = HcoConfig,           & ! HEMCO Configuration Object
                           RC        = GEOS_CHEM_RC,        & ! GC Status Code (Success?)
                           MPI_COMM  = WRF_DM_Comm)           ! MPI Communicator

      write(message_txt, *) "GIGC chemics_init: Return from GIGC_Chunk_Init return code ", GEOS_CHEM_RC
      call wrf_debug(1, message_txt)

      write(message_txt, *) "GIGC chemics_init: GIGC_Chunk_Init completed for ID #", grid%id
      call wrf_debug(1, message_txt)

      ! Debug chemistry state and species information on PET (1, 1, 1) lowest lev.
      ! WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP AFTER CHUNK-INIT/ CHMICS %%%%"
      ! do N = 1, MIN(State_Chm(grid%id)%nSpecies, 8)
      !     WRITE(6, *) "N:", N, "Name:", State_Chm(grid%id)%SpcData(N)%Info%Name
      !     ! WRITE(6, *) "Full Name:", State_Chm%SpcData(N)%Info%FullName
      !     WRITE(6, *) "Value at PET (1,1,1):", State_Chm(grid%id)%Species(1, 1, 1, N)
      !     WRITE(6, *) "Value at PET (7,1,1):", State_Chm(grid%id)%Species(7, 1, 1, N)
      ! enddo

      ! Output some debug information regarding Input_Opt object.
      write(6, *) "GIGC chemics_init: Following information regarding Input_Opt% (for all doms)"
      write(6, *) "After domain", grid%id
      write(6, *) "HPC: ", Input_Opt%HPC
      write(6, *) "LPRT: ", Input_Opt%LPRT
      write(6, *) "DoConv: ", Input_Opt%LCONV
      write(6, *) "DoDryDep: ", Input_Opt%LDRYD
      write(6, *) "DoEmis: ", Input_Opt%LEMIS
      write(6, *) "DoChem: ", Input_Opt%LCHEM
      write(6, *) "DoTurb: ", Input_Opt%LTURB
      write(6, *) "DoWetDep: ", Input_Opt%LWETD

      ! For debug only -- this will print out information for chemical species & ID#s
      ! useful if we are updating GEOS-Chem and need to reseed the WRF registry...
      ! WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP - DEBUG %%%%"
      ! do N = 1, State_Chm%nSpecies
      !     WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
      ! enddo
      ! End debug only

      ! Setup the chemistry state indexes in GIGC_Convert_State_Mod
      call GIGC_IdxSetup(Am_I_Root)

      ! Run the conversion routines from GEOS-Chem column chemistry code
      ! to store background values
      !! WARNING: scalar, num_scalar is passed into moist
      if(config_flags%have_bcs_chem .eq. .false.) then
        call GIGC_Set_WRF(Am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, 1, moist_dummy, & 
                          its, ite, jts, jte, &
                          ide, jde, &
                          kts, kte, &
                          Input_Opt, State_Met, State_Chm, State_Diag)
      endif

      ! Fill top layer of chem% (K is not paralleled so kde or kte will all work fine)
      do nv = 1, num_chem
         do j = jts, jte
            do i = its, ite
               chem(i, kde, j, nv) = chem(i, kde - 1, j, nv)
            enddo
         enddo
      enddo

      ! Setup time info for chemistry
      if (config_flags%restart) then
         call wrf_debug(15, "Setting last_chem_time from restart file")
         call WRFU_TimeSet(last_chem_time(id), &
                           YY=last_chem_time_year, &
                           MM=last_chem_time_month, &
                           DD=last_chem_time_day, &
                           H=last_chem_time_hour, &
                           M=last_chem_time_minute, &
                           S=last_chem_time_second)
      else
         call wrf_debug(15, "Setting last_chem_time to model start time-dt")
         call WRFU_TimeIntervalSet(tmpTimeInterval, s_=real(dt, 8))
         last_chem_time(id) = domain_get_current_time(grid) - tmpTimeInterval
      end if

      ! Setup Tropopause Diagnosis Data (module_tropopause)
      call tropopause_init( id, xlat, xlong, config_flags,   &
                            ids, ide, jds,jde, kds,kde,      &
                            ims, ime, jms,jme, kms,kme,      &
                            its, ite, jts,jte, kts,kte       )

      ! Keep the following for compatibility with WRF CAM/MAM Physics Scheme.
      gas_pcnst_modal_aero_pos = max(1, gas_pcnst_modal_aero) ! BSINGH - Fix for non-MAM simulations
END SUBROUTINE chem_init

!-----------------------------------------------------------------------
!  CHEMISTRY DEBUG: print_chem_species_index(chem_opt)
!  Outputs Index Numbers for Chemistry species
!-----------------------------------------------------------------------

#ifdef CHEM_DBG_I
   subroutine print_chem_species_index(chem_opt)
      use module_state_description
      implicit none

      integer, intent(in) :: chem_opt

      print *, 'Index numbers for chemistry species:'
!
! Gas species...
!
      select case (chem_opt)
      case (233)
         print *, p_no, "no"
         print *, p_o3, "o3"
      end select
   end subroutine print_chem_species_index
#endif
