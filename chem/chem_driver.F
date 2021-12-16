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
!   the original chem_driver.f from WRF-Chem v3.6.1
!   which is (c) their respective authors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chem_driver.F
!
! !DESCRIPTION: Subroutine chem_driver is the entry point for WRF to call chem
!  routines and adapted to enter (GI)GC entry point.
!\\
!\\
! !INTERFACE:
subroutine chem_driver(grid, config_flags &
#include "dummy_new_args.inc"
)

!
! !USES:
!
   !! WRF Interfaces
   use module_domain, only: domain
   use module_configure
   use module_driver_constants
   use module_machine
   use module_tiles
   use module_dm
   use module_model_constants
   use module_state_description

   !! Chemistry stub interfaces
   ! module_chem_utilities: Provides chem_prep, which passes grid meteo information to chemistry
   ! to local variables below.
   use module_chem_utilities

   ! module_input_chem_data: get chemistry information, namely last_chem_time (last runtime of chemistry),
   ! get_last_gas, last gas species index #, and a chem_dbg function
   use module_input_chem_data, only: last_chem_time, &
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
                                    chem_dbg, &
#endif
      get_last_gas

   ! upper_bc_driver: provides upper boundary conditions
   use module_upper_bc_driver, only: upper_bc_driver

   ! module_tropopause: diagnoses tropopause location
   use module_tropopause, only: tropopause_driver

   !! GIGC: Grid Independent GEOS-Chem
   ! GIGC_Chunk_Mod provides the GEOS-Chem Column Model Code
   use GIGC_Chunk_Mod
   use PRECISION_MOD ! Variable Precision code (fp) parameter


   ! Input_Opt_Mod provides the Input_Opt% type (OptInput) struct
   use Input_Opt_Mod, only: OptInput
   use State_Chm_Mod, only: ChmState
   use State_Met_Mod, only: MetState
   use State_Diag_Mod, only: DgnState

   ! Use our conversion module for GEOS-Chem code
   use GIGC_Convert_State_Mod

   ! Use our stateful module for GEOS-Chem code
   use GIGC_Stateful_Mod

   ! End of module use declarations
   implicit none

!
! !INPUT/OUTPUT PARAMETERS:
!
   ! provides information on WRF DMPAR parallel (MPI) process information
   ! namely whether this is the master process, using wrf_dm_on_monitor()
   logical, external :: wrf_dm_on_monitor

   !! Input data.
   type(domain), target :: grid

   !! Definitions of dummy arguments to solve
# include <dummy_new_decl.inc>
# define NO_I1_OLD

   ! namestep configuration flags (config_flags%chem_opt, for example)
   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
!
! !REMARKS:
!  This code is based upon the original interface for WRF-Chem and adapted to become
!  an abstraction layer for generic chemistry in WRF. See the "Pumpkin" project.
!
! !REVISION HISTORY:
!  30 May 2018 - H.P. Lin   - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   !-- ids           start index for i in domain
   !-- ide           end index for i in domain
   !-- jds           start index for j in domain
   !-- jde           end index for j in domain
   !-- kds           start index for k in domain
   !-- kde           end index for k in domain

   !-- ims           start index for i in memory
   !-- ime           end index for i in memory
   !-- jms           start index for j in memory
   !-- jme           end index for j in memory
   !-- kms           start index for k in memory
   !-- kme           end index for k in memory

   !-- its           start index for i in tile
   !-- ite           end index for i in tile
   !-- jts           start index for j in tile
   !-- jte           end index for j in tile
   !-- kts           start index for k in tile
   !-- kte           end index for k in tile

   !-- ips           start index for i in patch
   !-- ipe           end index for i in patch
   !-- jps           start index for j in patch
   !-- jpe           end index for j in patch
   integer :: ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              ips, ipe, jps, jpe, kps, kpe, &
              its, ite, jts, jte, kts, kte

   !! GEOS-Chem Integration Variables
   ! GEOS_CHEM_RC (RC in GEOS-Chem speak) is the return result code for all GIGC Subroutines.
   ! 0 = Success, > 0 is error
   integer :: GEOS_CHEM_RC

   ! Input Options
   TYPE(OptInput) :: Input_Opt      ! Input Options
   TYPE(MetState) :: State_Met      ! Meteorology state
   TYPE(ChmState) :: State_Chm      ! Chemistry state
   TYPE(DgnState) :: State_Diag     ! Diagnostics state

   ! Operator Options
   TYPE(GIGC_Chunk_Operators) :: GIGC_Ops

   ! For passing grid location (latlon) information
   real(4), allocatable :: lonCtr(:, :)  ! Lon Centers (rad)
   real(4), allocatable :: latCtr(:, :)  ! Lat Centers (rad)
   real(4), allocatable :: lonEdge(:, :) ! Lon Edges (rad)
   real(4), allocatable :: latEdge(:, :) ! Lat Edges (rad)

   ! For tiling information
   ! # of lons, lats, levs (I, J, L) in-PET
   integer :: IM, JM, LM

   ! Loop variables
   integer :: II, JJ

   ! Local Scalars
   integer :: i, j, k, l, numgas, nv, n, nr, ktau
   integer :: ijulian, nymd, nhms

   ! Check if first time running - if yes, skip because PBL Height is not computed in the first WRF timestep
   logical, save :: FIRST = .TRUE.

   ! Memorder XZY: %s-em31 is I, %s-em32 is K, %s-em33 is J
   ! Local Variables computed for chemistry usage.
   !
   ! To pass physics (meteo) variable to CHEMISTRY, edit this to create new local variables
   ! (sized as ims:ime, kms:kme, jms:jme)
   !
   ! Modify the call to module_chem_utilities::chem_prep,
   ! and modify the calls to external_driver as necessary.
   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      p_phy, u_phy, v_phy, t_phy, dz8w, t8w, p8w, rho, z_at_w, vvel, zmid, rh, pi_phy

   ! PBL Height, dimension is (I x J) obviously...
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: pbl_h

   ! Variables for calculating adaptive timesteps.
   TYPE(WRFU_TimeInterval) :: tmpTimeInterval
   real(KIND=8) :: curr_secs
   real(KIND=8) :: real_time_r8 ! ext. function in adapt_timestep_em.F
   logical      :: do_chemstep

   integer :: debug_level
   integer :: ij
   integer :: num_3d_m, ic, num_3d_c, num_3d_s

   ! From module_ctrans_grell - used for setting a minimum range on chemistry species.
   real :: epsilc
   parameter(epsilc=1.e-30)

   real :: chem_minval

   character(len=256) :: current_date_char
   integer :: WRF_dateY, WRF_dateM, WRF_dateD, WRF_dateH, WRF_dateI, WRF_dateS
   ! Y-M-D H:I:S

   real     :: WRF_minutesElapsed
   real(f4) :: WRF_hoursElapsed
   
   ! For debug
   character*256 :: message_txt
   character*256 :: debug_format

   ! For timing and scalability tests
   ! Note: Cannot use config_flags%max_dom in a SAVE statement, so we just save the outermost 8
   !       domains' information.
   real(KIND=8), dimension(8), save ::        &
                         WRFGC_Overhead_Time, &          ! WRF-GC Coupling Overhead Time, Generic [s]
                         WRFGC_GC_Time,       &          ! GEOS-Chem Column Code Time [s]
                         WRFGC_Phys_Time,     &          ! WRF-GC Coupling, Physics/Parameterization [s]
                         WRFGC_Diag_Time,     &          ! WRF-GC Coupling, Diagnostics [s]
                         WRF_Total_Time                  ! Total WRF Time (Estimate from MPI) [s]
   real(KIND=8), save :: WRF_Time_Last_Call
   real(KIND=8)       :: WRFGC_Time_Temp_Start, WRFGC_Time_Temp_End

   intrinsic max, min

   ! WRF DM (MPI) Parallel Information - is master process?
   logical :: am_I_Root
   integer :: WRF_DM_MyProc, WRF_DM_NProc, WRF_DM_Comm

   ! ---------------- Code begins here ---------------- !

   ! Measure time WRF between chem_driver calls
   if(FIRST) then
      ! Do not measure first timestep (initialization time is not taken into account)
   else
      ! The time spent in WRF is essentially the time between chem_driver calls.
      ! Though we cannot distinguish between domain times...
      !
      ! ...wrf does stuff...
      ! in chem_driver
      ! out chem_driver --- measure time 1
      ! ...wrf does stuf...
      ! in chem_driver  --- measure time 2 => time in WRF is 2 - 1
      ! out chem_driver --- measure time 1
      WRF_Total_Time(grid%id) = WRF_Total_Time(grid%id) + (MPI_Wtime() - WRF_Time_Last_Call)
   endif

   if(wrf_dm_on_monitor()) then
      am_I_Root = .true.
   else
      am_I_Root = .false.
   endif

   call wrf_get_nproc(WRF_DM_NProc)
   call wrf_get_myproc(WRF_DM_MyProc)
   call wrf_get_dm_communicator(WRF_DM_Comm)

   ! Get debug level
   call nl_get_debug_level(1, debug_level)
   call set_wrf_debug_level(debug_level)

   ! Options to run operators in GIGC_Chunk_Run.
   ! It now accepts configuration in the WRF namelist for each individual process,
   ! so you can turn off turbulence if it causes spurious results
   ! (hplin, 8/15/18)
   GIGC_Ops%Conv   = config_flags%gc_do_convection
   GIGC_Ops%Emis   = config_flags%gc_do_hemco
   GIGC_Ops%Tend   = .true.
   GIGC_Ops%Turb   = config_flags%gc_do_pblmix
   GIGC_Ops%Chem   = config_flags%gc_do_chemistry
   GIGC_Ops%DryDep = config_flags%gc_do_drydep
   GIGC_Ops%WetDep = config_flags%gc_do_wetdep

   GIGC_Ops%GCDiagn = config_flags%gc_do_gcdiagn

   ! Setup the adaptive timestep for the chem routines. Most of this follows
   ! what is in solve_em, except for the call to adjust time_step.
   !
   ! The necessary variables exist for the EM core and using the adaptive
   ! techniques will work even with a constant time step. In fact, they
   ! prevent issues with restarts and changed time steps. So, we will
   ! always use them with the EM core.
   ktau = grid%itimestep
   tmpTimeInterval = domain_get_time_since_sim_start(grid)
   curr_secs = real_time_r8(tmpTimeInterval)
   ijulian = ifix(grid%julian)

   ! initializing diagnostics and macros
   if (ktau == 1) then
      grid%conv_ct(:, :, :, :) = 0.
      grid%chem_ct(:, :, :, :) = 0.
      grid%vmix_ct(:, :, :, :) = 0.
   endif

   ! Check whether to do chemistry steps (do_chemstep logical)
   do_chemstep = .false.
   if (ktau == 1) then
      do_chemstep = .true.
      grid%ktauc = 1
   else
      if ((grid%chemdt <= 0) .or. &
          (curr_secs + real(grid%dt, 8) + 0.01 >= &
           (int(curr_secs/real(grid%chemdt*60., 8) + 1, 8)*real(grid%chemdt*60., 8))) &
          ) then
         do_chemstep = .true.
         grid%ktauc = grid%ktauc + 1
         last_chem_time(grid%id) = domain_get_current_time(grid)

         ! Don't be fooled, this is for ESMF, it's to Get data into WRFU's Ancient ESMF bindings,
         ! not get from grid! hplin 4/25/18
         call WRFU_TimeGet(last_chem_time(grid%id), &
                           YY=grid%last_chem_time_year, &
                           MM=grid%last_chem_time_month, &
                           DD=grid%last_chem_time_day, &
                           H=grid%last_chem_time_hour, &
                           M=grid%last_chem_time_minute, &
                           S=grid%last_chem_time_second)
      endif
   endif

   call get_ijk_from_grid(grid, &
                          ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          ips, ipe, jps, jpe, kps, kpe)

   call domain_clock_get(grid, current_timestr=current_date_char)

   ! Note this gets MINUTES first then convert it to a seconds number...
   call domain_clock_get(grid, minutesSinceSimulationStart=WRF_minutesElapsed)
   WRF_hoursElapsed = WRF_minutesElapsed / 60.0_f4

   ! write(6, *) "Pumpkin chem_driver: current_date_char = ", current_date_char
   ! 2018-02-22_00:00:00

   read (current_date_char(1:4), FMT='(I4)') WRF_dateY
   read (current_date_char(6:7), FMT='(I2)') WRF_dateM
   read (current_date_char(9:10), FMT='(I2)') WRF_dateD
   read (current_date_char(12:13), FMT='(I2)') WRF_dateH
   read (current_date_char(15:16), FMT='(I2)') WRF_dateI
   read (current_date_char(18:19), FMT='(I2)') WRF_dateS

   ! Compute nymd (YYYYMMDD) & nhms (hhmmss) for GEOS-Chem Column Code
   nymd = WRF_dateY * 10000 + WRF_dateM * 100 + WRF_dateD
   nhms = WRF_dateH * 10000 + WRF_dateI * 100 + WRF_dateS

   grid%raincv_b(:, :) = grid%raincv(:, :)

   ! Number of moist variables
   num_3d_m = num_moist
   num_3d_c = num_chem
   num_3d_s = num_scalar
   numgas = get_last_gas(config_flags%chem_opt)

   ! Compute these starting and stopping locations for each tile and number of tiles.
   call set_tiles(grid, ids, ide, jds, jde, ips, ipe, jps, jpe)

   chem_minval = epsilc ! chem_minval can be case dependant and set below...
   chem_select: select case(config_flags%chem_opt)
      case (233)
         call wrf_debug(15, 'GIGC chem_driver: chem_opt = 233. Running GEOS-Chem HP chemistry option.')
      case default
         call wrf_error_fatal("Pumpkin chem_driver: Unrecognized chem_opt. WRF-GCHP is chem_opt=233.")
   end select chem_select


   ! Some leftovers from original WRF chem_select. Don't touch unless you know what you are doing (I don't)
   ! $ hplin 20180327 $
   do j = jps, min(jde - 1, jpe)
      do k = kps, kpe
         do i = ips, min(ide - 1, ipe)
            vvel(i, k, j) = grid%w_2(i, k, j)
            zmid(i, k, j) = grid%z(i, k, j)
         enddo
      enddo
   enddo
   do j = jps, min(jde - 1, jpe)
      do i = ips, min(ide - 1, ipe)
         pbl_h(i, j) = grid%pblh(i, j)
      enddo
   enddo

!------------------------------------------------------------------------
! Main chemistry tile loop
!------------------------------------------------------------------------
   chem_tile_loop_1: DO ij = 1, grid%num_tiles
      ! hplin 5/7/2018: I have no idea why it loops only until ide-1, but
      ! I have removed this so allocs to XMID (GC) in chemics_init are not
      ! broken afterwards.
      !
      ! hplin 5/15/2018: Learn everyday. WRF does everything for a reason.
      ! Do NOT touch IDE, JDE -- ALWAYS allocate ONE less when you encounter
      ! the end boundaries. You have been warned.

      its = grid%i_start(ij)
      ite = min(grid%i_end(ij), ide - 1)
      ! ite = grid%i_end(ij)
      jts = grid%j_start(ij)
      ! jte = grid%j_end(ij)
      jte = min(grid%j_end(ij), jde - 1)

      kts = kps
      kte = min(kpe, kde - 1)

      ! Compute IM, JM, LM (lon, lat, lev) on this particular run instance
      IM = ite - its + 1
      JM = jte - jts + 1
      LM = kte - kts + 1

      ! Start WRF-GC Coupling Physics Timing (Timer-ID: chem_driver_3)
      ! Timer: WRFGC_Phys_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      CALL wrf_debug(15, 'Pumpkin chem_driver: Calling module_chem_utilities::chem_prep')

      ! chem_prep @ module_chem_utilities.f
      ! The chem_prep routine reads in the grid & configuration flags (grid%...)
      ! And does conversions onto memorder XZY local variables, so they can be used by chemistry.
      ! This is considered for timer purposes part of the physics/parameterization process.
      CALL chem_prep(config_flags, &
                     grid%u_2, grid%v_2, grid%p, grid%pb, &
                     grid%alt, grid%ph_2, grid%phb, grid%t_2, &
                     moist, num_3d_m, rho, &
                     p_phy, pi_phy, u_phy, v_phy, &
                     p8w, t_phy, t8w, grid%z, z_at_w, &
                     dz8w, rh, grid%fnm, grid%fnp, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, &
                     kps, kpe)

      ! End WRF-GC Coupling Physics Timing (Timer-ID: chem_driver_3)
      ! Timer: WRFGC_Phys_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Phys_Time(grid%id) = WRFGC_Phys_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_3

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (kps <= CHEM_DBG_K .and. kpe >= CHEM_DBG_K)) then
         call wrf_debug(15, 'Pumpkin chem_driver: Calling chem_dbg before doing chemistry')
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit)
      end if
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CHEMISTRY DRIVER/ MECHANISM (do_timestep)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (do_chemstep) then
      if(am_I_Root) then
         call wrf_debug(15, 'Pumpkin chem_driver: THIS IS THE MASTER PROCESS - CALLING CHEMISTRY')
      else
         write(message_txt, *) "Pumpkin chem_driver: Calling chemistry on subprocess ", WRF_DM_MyProc, "(+1) of ", WRF_DM_NProc
         call wrf_debug(15, message_txt)
      endif

      write(message_txt, *) 'GIGC chem_driver: Before State_Met conversion IM = ', IM, ' JM = ', JM
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chem_driver: its, ite, jts, jte, kts, kte: ", its, " ", ite, " ", jts, " ", jte, " ", kts, " ", kte
      call wrf_debug(15, message_txt)

      write(message_txt, *) "Pumpkin chem_driver: ims, ime, jms, jme, kms, kme: ", ims, " ", ime, " ", jms, " ", jme, " ", kms, " ", kme
      call wrf_debug(15, message_txt)
     
      write(message_txt, *) "Pumpkin chem_driver: ids, ide, jds, jde, kds, kde: ", ids, " ", ide, " ", jds, " ", jde, " ", kds, " ", kde
      call wrf_debug(15, message_txt)

      ! Start WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_0)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      ! Allocate array for lon/lat centers for passing lat/lon geo. information to GIGC. (hplin, 4/24/18)
      ! Added grid box edge information (xfeng, 11/11/18)
      allocate(lonCtr (IM,   JM))
      allocate(latCtr (IM,   JM))
      allocate(lonEdge(IM+1, JM))
      allocate(latEdge(IM,   JM+1))

      ! For every grid box populate with relevant info...
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

      ! Debugging only -- output longitude, latitudes for grid boxes, first in WRF format (so edges)
      ! then in GEOS-Chem format as diagnosed by WRF-GC coupling chem_driver
      ! write(6, *) "=============================================================="
      ! write(6, *) "D I A G N O S T I C   W R F   G R I D   C O O R D I N A T E S"
      ! write(6, *) "=============================================================="

      ! write(6, *) "its, ite, jts, jte", its, ite, jts, jte
      ! write(6, *) "IM, JM", IM, JM

      ! do i = its, ite+1
      !    do j = jts, jte+1
      !       write(6, *) "(", grid%XLONG(i, j), ",", grid%XLAT(i, j), ")"
      !    enddo
      !    write(6, *) " "
      ! enddo

      ! write(6, *) "===  G E O S - C H E M    L O N / L A T    C E N T E R S  ==="

      ! do i = 1, IM+1
      !    do j = 1, JM+1
      !       write(6, *) "(", lonCtr(i, j), ",", latCtr(i, j), ")"
      !    enddo
      !    write(6, *) " "
      ! enddo

      ! write(6, *) "GIGC chem_driver lonCtr debug:", lonCtr
      ! write(6, *) "GIGC chem_driver latCtr debug:", latCtr
      ! End debugging

      ! Retrieve stateful information from GIGC_Stateful_Mod (passing in grid#)
      call GIGC_State_Get_Opt(am_I_Root, Input_Opt)
      call GIGC_State_Get_Met(am_I_Root, grid%id, State_Met, GEOS_CHEM_RC)
      call GIGC_State_Get_Chm(am_I_Root, grid%id, State_Chm, GEOS_CHEM_RC)
      call GIGC_State_Get_Diag(am_I_Root, grid%id, State_Diag, GEOS_CHEM_RC)

      ! Pass some HPC Information to Input_Opt...
      Input_Opt%HPC = .true.
      Input_Opt%RootCPU = am_I_Root
      Input_Opt%myCPU = WRF_DM_MyProc
      Input_Opt%NPES = WRF_DM_NProc
      Input_Opt%haveImpRst = .true.

      ! write(6, *) 'GIGC chem_driver: NYMD, NHMS debug: nymd =', nymd, 'nhms =', nhms, 'ijulian =', ijulian
      if(WRF_dateM .eq. 2 .and. WRF_dateD .eq. 15 .and. WRF_dateH .eq. 3 .and. WRF_dateI .eq. 10) then
            call wrf_debug(1, "        ___iiiii___            ")
            call wrf_debug(1, "       |           |           ")
            call wrf_debug(1, "     __|___________|__         ")
            call wrf_debug(1, "    |^^^^^^^^^^^^^^^^^|        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    |                 |        ")
            call wrf_debug(1, "    ~~~~~~~~~~~~~~~~~~~        ")
            call wrf_debug(1, " hplin 02/15/1997 wrf|gchp     ")
      endif

      ! Diagnose the tropopause location using WRF-Chem's tropopause_driver.
      call wrf_debug(15, 'Pumpkin chem_driver: Diagnosing tropopause location')
      call tropopause_driver(grid%id, grid%dt, current_date_char,           &
                               t_phy, p_phy, p8w, zmid, z_at_w,             &
                               grid%tropo_lev, grid%tropo_p,  grid%tropo_z, &
                               ids, ide, jds, jde, kds, kde,                &
                               ims, ime, jms, jme, kms, kme,                &
                               its, min(grid%i_end(ij), ide - 1), jts, min(grid%j_end(ij), jde - 1), kts, kte)
      call wrf_debug(15, 'Pumpkin chem_driver: After tropopause_driver')

      ! WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP BEFORE CHUNK-CNV %%%%"
      ! do N = 1, MIN(State_Chm%nSpecies, 8)
      !     WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
      !     ! WRITE(6, *) "Full Name:", State_Chm%SpcData(N)%Info%FullName
      !     WRITE(6, *) "Value at PET (1,1,1):", State_Chm%Species(1, 1, 1, N)
      !     WRITE(6, *) "Value at PET (2,1,1):", State_Chm%Species(2, 1, 1, N)
      !     WRITE(6, *) "Value at PET (7,1,1):", State_Chm%Species(7, 1, 1, N)
      ! enddo

      ! Run the conversion routines for GEOS-Chem column chemistry code
      ! Now passes full tile information for conversion (hplin, 5/7/2018)
      call GIGC_Get_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        dz8w, p8w, t8w, pi_phy, &
                        f_qc, f_qi, f_qs,   &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        ids, jds, kds, &
                        curr_secs, &
                        Input_Opt, State_Met, State_Chm, &
                        WRFGC_Phys_Time(grid%id))

      call wrf_debug(15, 'GIGC chem_driver: After GIGC_Get_WRF conversion')

      ! End WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_0)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_0

      ! Start WRF-GC GC Column Timing (Timer-ID: chem_driver_1)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      ! If first timestep, some meteorology data isn't quite ready and needs some computing
      ! before running GEOS-Chem. Hence skip first timestep.
      if(.not. FIRST) then

         call wrf_debug(15, 'GIGC chem_driver: Before GIGC_Chunk_Run')
         ! Run the GEOS-Chem column chemistry code for the given phase
         ! $ hplin 4/3/2018 1419 $
         call GIGC_Chunk_Run( am_I_Root  = am_I_Root,                      & ! Is this the root PET?
                              IM         = IM,                             & ! # of lons on this PET
                              JM         = JM,                             & ! # of lats on this PET
                              LM         = LM,                             & ! # of levs on this PET
                              ID         = grid%id,                        & ! Domain identifier, this PET
                              nymd       = nymd,                           & ! Current YYYYMMDD date - INTEGER
                              nhms       = nhms,                           & ! Current hhmmss time - INTEGER
                              year       = WRF_dateY,                      & ! Current year - INTEGER
                              month      = WRF_dateM,                      & ! Current month - INTEGER
                              day        = WRF_dateD,                      & ! Current day - INTEGER
                              dayOfYr    = ijulian,                        & ! Current day of year - INTEGER (grid%julian)
                              hour       = WRF_dateH,                      & ! Current hour - INTEGER
                              minute     = WRF_dateI,                      & ! Current minute - INTEGER
                              second     = WRF_dateS,                      & ! Current second - INTEGER
                              hElapsed   = WRF_hoursElapsed,               & ! Runs hours elapsed [r4]
                              lonCtr     = lonCtr,                         & ! Column Lon Centers [rad]
                              latCtr     = latCtr,                         & ! Column Lat Centers [rad]
                              lonEdge    = lonEdge,                        & ! Column Lon Edges [rad]
                              latEdge    = latEdge,                        & ! Column Lat Edges [rad]
                              Input_Opt  = Input_Opt,                      & ! Input Options
                              State_Chm  = State_Chm,                      & ! Chemistry State
                              State_Met  = State_Met,                      & ! Meteorology State
                              State_Diag = State_Diag,                     & ! Diagnostics State
                              Operators  = GIGC_Ops,                       & ! Run operator options
                              IsChemTime = .true.,                         & ! Is it time for chem?
                              RC         = GEOS_CHEM_RC)                     ! Success or failure?
         call wrf_debug(15, 'GIGC chem_driver: After GIGC_Chunk_Run')

      endif

      ! End WRF-GC GC Column Timing (Timer-ID: chem_driver_1)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_GC_Time(grid%id) = WRFGC_GC_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_1

      ! Start WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_2)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_Start = MPI_Wtime()

      ! Run the conversion routines from GEOS-Chem column chemistry code
      call GIGC_Set_WRF(am_I_Root, config_flags, grid, num_chem, chem, num_scalar, scalar, num_moist, moist, &
                        its, ite, jts, jte, &
                        ide, jde, &
                        kts, kte, &
                        Input_Opt, State_Met, State_Chm, State_Diag)

      call wrf_debug(15, 'GIGC chem_driver: After GIGC_Set_WRF conversion')

      ! Prevent memory leak and deallocate GIGC arrays for next move
      call GIGC_Chunk_Cleanup(config_flags, kds, kde, Input_Opt, State_Met, State_Chm)

      ! Save state information into GIGC_Stateful_Mod (hplin, 6/12/18)
      call GIGC_State_Set_Met(am_I_Root, grid%id, State_Met, GEOS_CHEM_RC)
      call GIGC_State_Set_Chm(am_I_Root, grid%id, State_Chm, GEOS_CHEM_RC)
      call GIGC_State_Set_Diag(am_I_Root, grid%id, State_Diag, GEOS_CHEM_RC)

      ! End WRF-GC Coupling Overhead Timing (Timer-ID: chem_driver_2)
      ! Timer: WRFGC_Overhead_Time
      WRFGC_Time_Temp_End = MPI_Wtime()
      WRFGC_Overhead_Time(grid%id) = WRFGC_Overhead_Time(grid%id) + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
      ! End Timer chem_driver_2


      ! Output diagnostic debug output for timers (do_chemstep)
      ! Debug output for timers --
      if(debug_level .ge. 1) then
         ! Output timer information
         debug_format = "(2x, a30, 4x, f14.3)"
         write(6, *) "================================================="
         write(6, *) "W R F - G C   T I M E R S"
         write(6, *) "   Domain", grid%id, "(Units: s)"
         write(6, *) "================================================="
         write(6, debug_format) "=> Coupler Component          ", WRFGC_Overhead_Time(grid%id)
         write(6, debug_format) " -- Physics & Parameterization", WRFGC_Phys_Time(grid%id)
         write(6, debug_format) " -- Diagnostics Out           ", WRFGC_Diag_Time(grid%id)
         write(6, debug_format) "=> GEOS-Chem Column           ", WRFGC_GC_Time(grid%id)
         write(6, debug_format) "===> Total Chemistry          ", WRFGC_Overhead_Time(grid%id)+WRFGC_GC_Time(grid%id)
         write(6, debug_format) "=> WRF                        ", WRF_Total_Time(grid%id)
         write(6, debug_format) "===> Total WRF-GC             ", WRFGC_Overhead_Time(grid%id)+WRFGC_GC_Time(grid%id)+WRF_Total_Time(grid%id)
      endif
   end if ! Chemistry time step check

   ! Fill top level to prevent spurious interpolation results (no extrapolation)
   do nv = 1, num_chem
      do j = jts, jte
         do i = its, ite
            chem(i, kpe, j, nv) = chem(i, kte, j, nv)
         enddo
      enddo
   enddo

!-- set upper boundary condition
   if (config_flags%have_bcs_upper) then
      call wrf_debug(15, 'Pumpkin chem_driver: set upper boundary condition')
      call upper_bc_driver(grid%id, grid%dt, current_date_char, &
                           chem, p_phy, p8w, grid%tropo_lev, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)
   endif

   call wrf_debug(15, 'Pumpkin chem_driver: Exiting chem_driver (chem_tile_loop_1)')
end do chem_tile_loop_1

! No longer first call
   FIRST = .FALSE.

! For WRF time measurement
   WRF_Time_Last_Call = MPI_Wtime()

end subroutine chem_driver

!
! While I am waiting for unit tests to complete, I would like to bury a message here at
! this lonely corner of the code. (hplin, 8/11/18)
!
! I would like to take this opportunity to thank those along the long and tedious,
! no doubt frustrating journey of the WRF-GC (formerly "WRF-GCHP") project.
!
! This project was envisioned by my group's leader, Prof. Tzung-May Fu, and assigned to
! me and my fellow collaborator Feng Xu, in the winter of 2017. Without Prof. Fu's
! tremendous support and trust that we could make it, this project would never come into
! fruition as it has today. Thank you.
!
! I'd also like to thank Feng Xu, who has provided the day-to-day scientific support
! that was crucial as I looked at WRF & GEOS-Chem from a programmer's perspective only.
! Thanks for bearing with me as I broke, built, and broke things again and again,
! gradually building the core architecture that makes WRF-GC today.
!
! Thanks to the Harvard Atmospheric Chemistry Modeling Group & the GEOS-Chem support team,
! without whose help, especially integrating the changes (or disruptions) WRF-GC 
! has "constributed" to the GEOS-Chem code base, WRF-GC would not be in this position today.
! Special thanks to Prof. Daniel Jacob, who provided me with an amazing opportunity to
! be in Harvard in Summer '18 to work closer with the GEOS-Chem community.
!
! I am deeply indebted to my family, who although hear nothing specific about my work,
! have no doubt had more than enough of my hair-pulling mumbling while I worked late
! weekend nights. Their unconditional support of my endeavors throughout my life have
! been vital to make me who I am today.
!
! Last but not least, thanks to my fellow friend and companion "Angela" Huang, whom
! I had met mid-April at the peak of my frustrations with the WRF-to-GCHP interface.
! Long walks and rants at night in the PKU campus have helped provide the emotional
! support necessary to push this project further.
!
! Thank you everyone.
!                                                             - Haipeng Lin, 2018.8.11
! 
!
