!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   WRF-GCHP
!   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
!
!   WRF & GCHP are (c) their original authors.
!   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
!
!   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng <fengx7@pku.edu.cn>, 2018-04
!   Peking University, School of Physics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_convert_state_mod.F
!
! !DESCRIPTION: Module GIGC\_Convert\_State\_Mod performs conversion of WRF
!  data structures to/from Grid-Independent GEOS-Chem (GIGC) State structures.
!  (hplin, 4/3/18)
!\\
!\\fyxz
! !INTERFACE:
module GIGC_Convert_State_Mod

!
! !USES:
!
    ! WRF
    use module_domain
    use module_dm
    use module_model_constants
    use module_state_description
    use module_chem_utilities, only: calc_zenithb, getpar, get_cloud_optical_depth
    use module_convection_prep

    ! GEOS-Chem
    use GIGC_Chunk_Mod
    use PRECISION_MOD
    use Species_Mod,    only: MISSING_INT
    use ERROR_MOD,      only: IT_IS_NAN
    use Input_Opt_Mod,  only: OptInput
    use State_Chm_Mod,  only: ChmState, Ind_
    use State_Met_Mod,  only: MetState
    use State_Diag_Mod, only: DgnState

    implicit none
    private

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: GIGC_IdxSetup
    public :: GIGC_Get_WRF
    public :: GIGC_Set_WRF
    public :: GIGC_Chunk_Cleanup
!
! !PRIVATE VARIABLES:
!
    integer :: gi_no, gi_o3, gi_pan, gi_co, gi_alk4, gi_isop, gi_hno3, gi_h2o2, gi_acet, gi_mek, gi_ald2, gi_rcho, gi_mvk, gi_macr, gi_npmn, &
               gi_ppn, gi_r4n2, gi_prpe, gi_c3h8, gi_ch2o, gi_c2h6, gi_n2o5, gi_hno4, gi_mp, gi_dms, gi_so2, gi_so4, gi_so4s, gi_msa, &
               gi_nh3, gi_nh4, gi_nit, gi_nits, gi_bcpi, gi_ocpi, gi_bcpo, gi_ocpo, gi_dst1, gi_dst2, gi_dst3, gi_dst4, gi_sala, gi_salc
    integer :: gi_br2, gi_br, gi_bro, gi_hobr, gi_hbr, gi_brno2, gi_brno3, gi_chbr3, gi_ch2br2, gi_ch3br, gi_mpn, gi_isopnd, gi_isopnb, gi_moba, &
               gi_propnn, gi_hac, gi_glyc, gi_mvkn, gi_macrn, gi_map, gi_no2, gi_no3, gi_hno2, gi_benz, gi_tolu, gi_xyle, gi_mtpa, gi_limo
    integer :: gi_mtpo, gi_soap, gi_soas, gi_eoh, gi_mgly, gi_glyx, gi_acta, gi_hpald, gi_dhdn, gi_ethln, gi_hcooh, gi_iepoxa, gi_iepoxb, &
               gi_iepoxd, gi_isn1, gi_ripa, gi_ripb, gi_ripd, gi_imae, gi_soaie, gi_soame, gi_soagx, gi_soamg, gi_lvoc, gi_lvocoa, gi_isn1og
    integer :: gi_isn1oa, gi_monits, gi_monitu, gi_honit, gi_ionita, gi_monita, gi_indiol, gi_ipmn, gi_hc187, gi_n2o, gi_ocs, gi_ch4, gi_brcl, &
               gi_hcl, gi_ccl4, gi_ch3cl, gi_ch3ccl3, gi_cfc113, gi_cfc114, gi_cfc115, gi_hcfc123, gi_hcfc141b, gi_hcfc142b, gi_cfc11, gi_cfc12, &
               gi_hcfc22, gi_h1211, gi_h1301, gi_h2402, gi_cl, gi_clo, gi_hocl, gi_clno3, gi_clno2, gi_cloo, gi_oclo, gi_cl2, gi_cl2o2, gi_h2o
    integer :: gi_brsala, gi_brsalc, gi_chcl3, gi_ch2cl2, gi_ch3i, gi_ch2i2, gi_ch2icl, gi_ch2ibr, gi_hoi, gi_i2, gi_ibr, gi_icl, gi_i, gi_io, &
               gi_hi, gi_oio, gi_ino, gi_iono, gi_iono2, gi_i2o2, gi_i2o3, gi_i2o4, gi_isala, gi_isalc, gi_aeri, gi_pfe, gi_co2, gi_lbro2h, &
               gi_lbro2n, gi_lisopoh, gi_lisopno3, gi_ltro2h, gi_ltro2n, gi_lxro2h, gi_lxro2n, gi_pyac, gi_so4h1, gi_so4h2, gi_pox, gi_lox
    integer :: gi_pco, gi_lco, gi_pso4, gi_lch4, gi_ph2o2, gi_dhdc, gi_pmnn, gi_tro2, gi_bro2, gi_n, gi_xro2, gi_imao3, gi_etp, gi_ra3p, gi_rb3p, &
               gi_rp, gi_pp, gi_prpn, gi_pip, gi_r4p, gi_dhpcarp, gi_iap, gi_hpc52o2, gi_vrp, gi_atooh, gi_dhmob, gi_isnp, gi_maop, gi_mrp, gi_inpn
    integer :: gi_mobaoo, gi_diboo, gi_h, gi_macrno2, gi_roh, gi_isnohoo, gi_mvkoo, gi_isnoob, gi_gaoo, gi_ch3choo, gi_iepoxoo, gi_mglyoo, gi_prn1, &
               gi_mgloo, gi_a3o2, gi_man2, gi_isnooa, gi_po2, gi_isopndo2, gi_b3o2, gi_macroo, gi_r4n1, gi_maopo2, gi_ato2, gi_rco3, gi_limo2, &
               gi_olnd, gi_olnn, gi_ko2, gi_isopnbo2, gi_pio2, gi_hc5oo, gi_nmao3, gi_vro2, gi_ch2oo, gi_hc5, gi_rio2, gi_eto2, gi_ino2, gi_r4o2
    integer :: gi_mro2, gi_mco3, gi_mo2, gi_o1d, gi_o, gi_ho2, gi_oh, gi_h2, gi_moh, gi_n2, gi_o2, gi_rcooh
!   add mercury, yxzhang, 06/20/2019
    integer :: gi_hg0, gi_hg2, gi_hgp
!
! !REMARKS:
!  This module provides Get\_WRF and Set\_WRF routines that will "Get" and "Set"
!  variables from/to WRF, respectively, to/from the State_Met, State_Chm vars.
!
! !REVISION HISTORY:
!  24 Apr 2018 - H.P. Lin  - First crack.
!  25 Apr 2018 - H.P. Lin  - Refactored to support WRF (i,j) bounds, dbg options,
!                            and WRF scalar, chem, moist arrays from _new_args.inc
!  7 May 2018  - H.P. Lin  - Allocate arrays when necessary only & support tiling instead of 1x1 col
!  17 May 2018 - H.P. Lin  - Added ProTeX headers and some variables from xfeng, 5/18.
!  12 Oct 2018 - H.P. Lin  - Add State_Diag for GEOS-Chem NetCDF Diagnostics
!  14 Nov 2018 - X. Feng   - Add Cumulus parameterization for ntiedtke
!  21 Nov 2018 - H.P. Lin  - Massive performance optimizations
!------------------------------------------------------------------------------
!BOC
contains
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_IdxSetup
!
! !DESCRIPTION: Subroutine GIGC\_IdxSetup saves the State_Chm IND_ results into
!  a local scalar variable set so that we do not check the map every time from G-C,
!  when converting between internal states. (hplin, 11/21/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_IdxSetup(am_I_Root)
!
! !INPUT PARAMETERS:
!
        logical, intent(in) :: am_I_Root

        ! Code begins here!
        gi_no = IND_('NO')
        gi_o3 = IND_('O3')
        gi_pan = IND_('PAN')
        gi_co = IND_('CO')
        gi_alk4 = IND_('ALK4')
        gi_isop = IND_('ISOP')
        gi_hno3 = IND_('HNO3')
        gi_h2o2 = IND_('H2O2')
        gi_acet = IND_('ACET')
        gi_mek = IND_('MEK')
        gi_ald2 = IND_('ALD2')
        gi_rcho = IND_('RCHO')
        gi_mvk = IND_('MVK')
        gi_macr = IND_('MACR')
        gi_npmn = IND_('NPMN')
        gi_ppn = IND_('PPN')
        gi_r4n2 = IND_('R4N2')
        gi_prpe = IND_('PRPE')
        gi_c3h8 = IND_('C3H8')
        gi_ch2o = IND_('CH2O')
        gi_c2h6 = IND_('C2H6')
        gi_n2o5 = IND_('N2O5')
        gi_hno4 = IND_('HNO4')
        gi_mp = IND_('MP')
        gi_dms = IND_('DMS')
        gi_so2 = IND_('SO2')
        gi_so4 = IND_('SO4')
        gi_so4s = IND_('SO4S')
        gi_msa = IND_('MSA')
        gi_nh3 = IND_('NH3')
        gi_nh4 = IND_('NH4')
        gi_nit = IND_('NIT')
        gi_nits = IND_('NITS')
        gi_bcpi = IND_('BCPI')
        gi_ocpi = IND_('OCPI')
        gi_bcpo = IND_('BCPO')
        gi_ocpo = IND_('OCPO')
        gi_dst1 = IND_('DST1')
        gi_dst2 = IND_('DST2')
        gi_dst3 = IND_('DST3')
        gi_dst4 = IND_('DST4')
        gi_sala = IND_('SALA')
        gi_salc = IND_('SALC')
        gi_br2 = IND_('BR2')
        gi_br = IND_('BR')
        gi_bro = IND_('BRO')
        gi_hobr = IND_('HOBR')
        gi_hbr = IND_('HBR')
        gi_brno2 = IND_('BRNO2')
        gi_brno3 = IND_('BRNO3')
        gi_chbr3 = IND_('CHBR3')
        gi_ch2br2 = IND_('CH2BR2')
        gi_ch3br = IND_('CH3BR')
        gi_mpn = IND_('MPN')
        gi_isopnd = IND_('ISOPND')
        gi_isopnb = IND_('ISOPNB')
        gi_moba = IND_('MOBA')
        gi_propnn = IND_('PROPNN')
        gi_hac = IND_('HAC')
        gi_glyc = IND_('GLYC')
        gi_mvkn = IND_('MVKN')
        gi_macrn = IND_('MACRN')
        gi_map = IND_('MAP')
        gi_no2 = IND_('NO2')
        gi_no3 = IND_('NO3')
        gi_hno2 = IND_('HNO2')
        gi_benz = IND_('BENZ')
        gi_tolu = IND_('TOLU')
        gi_xyle = IND_('XYLE')
        gi_mtpa = IND_('MTPA')
        gi_limo = IND_('LIMO')
        gi_mtpo = IND_('MTPO')
        gi_soap = IND_('SOAP')
        gi_soas = IND_('SOAS')
        gi_eoh = IND_('EOH')
        gi_mgly = IND_('MGLY')
        gi_glyx = IND_('GLYX')
        gi_acta = IND_('ACTA')
        gi_hpald = IND_('HPALD')
        gi_dhdn = IND_('DHDN')
        gi_ethln = IND_('ETHLN')
        gi_hcooh = IND_('HCOOH')
        gi_iepoxa = IND_('IEPOXA')
        gi_iepoxb = IND_('IEPOXB')
        gi_iepoxd = IND_('IEPOXD')
        gi_isn1 = IND_('ISN1')
        gi_ripa = IND_('RIPA')
        gi_ripb = IND_('RIPB')
        gi_ripd = IND_('RIPD')
        gi_imae = IND_('IMAE')
        gi_soaie = IND_('SOAIE')
        gi_soame = IND_('SOAME')
        gi_soagx = IND_('SOAGX')
        gi_soamg = IND_('SOAMG')
        gi_lvoc = IND_('LVOC')
        gi_lvocoa = IND_('LVOCOA')
        gi_isn1og = IND_('ISN1OG')
        gi_isn1oa = IND_('ISN1OA')
        gi_monits = IND_('MONITS')
        gi_monitu = IND_('MONITU')
        gi_honit = IND_('HONIT')
        gi_ionita = IND_('IONITA')
        gi_monita = IND_('MONITA')
        gi_indiol = IND_('INDIOL')
        gi_ipmn = IND_('IPMN')
        gi_hc187 = IND_('HC187')
        gi_n2o = IND_('N2O')
        gi_ocs = IND_('OCS')
        gi_ch4 = IND_('CH4')
        gi_brcl = IND_('BRCL')
        gi_hcl = IND_('HCL')
        gi_ccl4 = IND_('CCL4')
        gi_ch3cl = IND_('CH3CL')
        gi_ch3ccl3 = IND_('CH3CCL3')
        gi_cfc113 = IND_('CFC113')
        gi_cfc114 = IND_('CFC114')
        gi_cfc115 = IND_('CFC115')
        gi_hcfc123 = IND_('HCFC123')
        gi_hcfc141b = IND_('HCFC141B')
        gi_hcfc142b = IND_('HCFC142B')
        gi_cfc11 = IND_('CFC11')
        gi_cfc12 = IND_('CFC12')
        gi_hcfc22 = IND_('HCFC22')
        gi_h1211 = IND_('H1211')
        gi_h1301 = IND_('H1301')
        gi_h2402 = IND_('H2402')
        gi_cl = IND_('CL')
        gi_clo = IND_('CLO')
        gi_hocl = IND_('HOCL')
        gi_clno3 = IND_('CLNO3')
        gi_clno2 = IND_('CLNO2')
        gi_cloo = IND_('CLOO')
        gi_oclo = IND_('OCLO')
        gi_cl2 = IND_('CL2')
        gi_cl2o2 = IND_('CL2O2')
        gi_h2o = IND_('H2O')
        gi_brsala = IND_('BRSALA')
        gi_brsalc = IND_('BRSALC')
        gi_chcl3 = IND_('CHCL3')
        gi_ch2cl2 = IND_('CH2CL2')
        gi_ch3i = IND_('CH3I')
        gi_ch2i2 = IND_('CH2I2')
        gi_ch2icl = IND_('CH2ICL')
        gi_ch2ibr = IND_('CH2IBR')
        gi_hoi = IND_('HOI')
        gi_i2 = IND_('I2')
        gi_ibr = IND_('IBR')
        gi_icl = IND_('ICL')
        gi_i = IND_('I')
        gi_io = IND_('IO')
        gi_hi = IND_('HI')
        gi_oio = IND_('OIO')
        gi_ino = IND_('INO')
        gi_iono = IND_('IONO')
        gi_iono2 = IND_('IONO2')
        gi_i2o2 = IND_('I2O2')
        gi_i2o3 = IND_('I2O3')
        gi_i2o4 = IND_('I2O4')
        gi_isala = IND_('ISALA')
        gi_isalc = IND_('ISALC')
        gi_aeri = IND_('AERI')
        gi_pfe = IND_('PFE')
        gi_co2 = IND_('CO2')
        gi_lbro2h = IND_('LBRO2H')
        gi_lbro2n = IND_('LBRO2N')
        gi_lisopoh = IND_('LISOPOH')
        gi_lisopno3 = IND_('LISOPNO3')
        gi_ltro2h = IND_('LTRO2H')
        gi_ltro2n = IND_('LTRO2N')
        gi_lxro2h = IND_('LXRO2H')
        gi_lxro2n = IND_('LXRO2N')
        gi_pyac = IND_('PYAC')
        gi_so4h1 = IND_('SO4H1')
        gi_so4h2 = IND_('SO4H2')
        gi_pox = IND_('POX')
        gi_lox = IND_('LOX')
        gi_pco = IND_('PCO')
        gi_lco = IND_('LCO')
        gi_pso4 = IND_('PSO4')
        gi_lch4 = IND_('LCH4')
        gi_ph2o2 = IND_('PH2O2')
        gi_dhdc = IND_('DHDC')
        gi_pmnn = IND_('PMNN')
        gi_tro2 = IND_('TRO2')
        gi_bro2 = IND_('BRO2')
        gi_n = IND_('N')
        gi_xro2 = IND_('XRO2')
        gi_imao3 = IND_('IMAO3')
        gi_etp = IND_('ETP')
        gi_ra3p = IND_('RA3P')
        gi_rb3p = IND_('RB3P')
        gi_rp = IND_('RP')
        gi_pp = IND_('PP')
        gi_prpn = IND_('PRPN')
        gi_pip = IND_('PIP')
        gi_r4p = IND_('R4P')
        gi_dhpcarp = IND_('DHPCARP')
        gi_iap = IND_('IAP')
        gi_hpc52o2 = IND_('HPC52O2')
        gi_vrp = IND_('VRP')
        gi_atooh = IND_('ATOOH')
        gi_dhmob = IND_('DHMOB')
        gi_isnp = IND_('ISNP')
        gi_maop = IND_('MAOP')
        gi_mrp = IND_('MRP')
        gi_inpn = IND_('INPN')
        gi_mobaoo = IND_('MOBAOO')
        gi_diboo = IND_('DIBOO')
        gi_h = IND_('H')
        gi_macrno2 = IND_('MACRNO2')
        gi_roh = IND_('ROH')
        gi_isnohoo = IND_('ISNOHOO')
        gi_mvkoo = IND_('MVKOO')
        gi_isnoob = IND_('ISNOOB')
        gi_gaoo = IND_('GAOO')
        gi_ch3choo = IND_('CH3CHOO')
        gi_iepoxoo = IND_('IEPOXOO')
        gi_mglyoo = IND_('MGLYOO')
        gi_prn1 = IND_('PRN1')
        gi_mgloo = IND_('MGLOO')
        gi_a3o2 = IND_('A3O2')
        gi_man2 = IND_('MAN2')
        gi_isnooa = IND_('ISNOOA')
        gi_po2 = IND_('PO2')
        gi_isopndo2 = IND_('ISOPNDO2')
        gi_b3o2 = IND_('B3O2')
        gi_macroo = IND_('MACROO')
        gi_r4n1 = IND_('R4N1')
        gi_maopo2 = IND_('MAOPO2')
        gi_ato2 = IND_('ATO2')
        gi_rco3 = IND_('RCO3')
        gi_limo2 = IND_('LIMO2')
        gi_olnd = IND_('OLND')
        gi_olnn = IND_('OLNN')
        gi_ko2 = IND_('KO2')
        gi_isopnbo2 = IND_('ISOPNBO2')
        gi_pio2 = IND_('PIO2')
        gi_hc5oo = IND_('HC5OO')
        gi_nmao3 = IND_('NMAO3')
        gi_vro2 = IND_('VRO2')
        gi_ch2oo = IND_('CH2OO')
        gi_hc5 = IND_('HC5')
        gi_rio2 = IND_('RIO2')
        gi_eto2 = IND_('ETO2')
        gi_ino2 = IND_('INO2')
        gi_r4o2 = IND_('R4O2')
        gi_mro2 = IND_('MRO2')
        gi_mco3 = IND_('MCO3')
        gi_mo2 = IND_('MO2')
        gi_o1d = IND_('O1D')
        gi_o = IND_('O')
        gi_ho2 = IND_('HO2')
        gi_oh = IND_('OH')
        gi_h2 = IND_('H2')
        gi_moh = IND_('MOH')
        gi_n2 = IND_('N2')
        gi_o2 = IND_('O2')
        gi_rcooh = IND_('RCOOH')
!       Add mercury, yxzhang, 06/20/2019
        gi_hg0 = IND_('HG0')
        gi_hg2 = IND_('HG2')
        gi_hgp = IND_('HGP')      
    end subroutine GIGC_IdxSetup
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_Get_WRF
!
! !DESCRIPTION: Subroutine GIGC\_Get\_WRF converts WRF data structures (grid, chem,
!  moist, scalar) according to WRF configuration (config_flags) to GIGC data structures
!  to call the GIGC_Chunk_Mod routines (State_Met, State_Chm). (hplin, 4/24/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_Get_WRF(am_I_Root, &
        config_flags, grid, &
        num_chem, chem, num_scalar, scalar, num_moist, moist, &
        dz8w, p8w, t8w, pi_phy, &
        f_qc, f_qi, f_qs, &
        its, ite, jts, jte, &
        ide, jde, &
        kts, kte, &
        ids, jds, kds, &
        curr_secs, &
        Input_Opt, State_Met, State_Chm, &
        WRFGC_Phys_Time)
        use PRESSURE_MOD,  only: ACCEPT_EXTERNAL_APBP

        implicit none

!
! !INPUT PARAMETERS:
!
        logical, intent(in) :: am_I_Root
        type(grid_config_rec_type), intent(in) :: config_flags

        type(domain), target :: grid
        integer, intent(in)  :: num_chem, num_scalar, num_moist ! # of WRF chem, moist variables. Check "dummy_new_args.inc" & "actual_new_args.inc" in WRF compile.
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(in)   :: chem
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(in) :: scalar
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_moist), intent(in)  :: moist
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33), intent(in)  :: dz8w, p8w, pi_phy, t8w
        integer, intent(in) :: its, ite, jts, jte, ide, jde, kts, kte, ids, jds, kds

        real(KIND=8) :: curr_secs
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm

        real(KIND=8), intent(inout)   :: WRFGC_Phys_Time  ! Cumulative time spent in physics/parameterization calculation time

!
! !REMARKS:
!  This code was originally created for a 1x1 (truly 1-D) column conversion, but later
!  changed to perform computation on the WRF parallelization tile (i,j,kts~e). These
!  are mapped to the GIGC_Chunk arrays II, JJ, k, dimension IM, JM, LM.
!
!  Please do take note that in WRF, memorder is XZY so arrays are (i, k, j), whereas in
!  GEOS-Chem the memorder is XYZ.
!
!  Also, the State_Chm conversion routines are auto-generated. Please use
!  https://github.com/jimmielin/wrf-gchp-dev-tools to generate them.
!  (hplin, 5/17/18)
!
! !REVISION HISTORY:
!  17 May 2018 - H.P. Lin   - Added ProTeX headers
!  29 May 2018 - X. Feng    - Fixed original extrapolated dummy %PEDGE
!  17 Aug 2018 - H.P. Lin   - Diagnose ApBp from WRF hybrid_grid and use those if available
!  14 Nov 2018 - X. Feng    - Add ntiedtke cumulus parameterization
!  14 Nov 2018 - H.P. Lin   - Add timers for physics parameterization & diagnostic output
!  01 Dec 2018 - H.P. Lin   - Add land type use code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: debug_level      ! WRF Debug Level
        integer :: GEOS_CHEM_RC     ! GEOS-Chem Return Code (0 = Success)

        integer :: IM, II           ! # of lons on this PET
        integer :: JM, JJ           ! # of lats on this PET
        integer :: LM               ! # of levs on this PET
        integer :: i, j, k          ! Loop variable over xts, xte (x = i, j, k)
        integer :: N

        real :: zen, pardb, pardif  ! For computing %SUNCOSmid, %PARDB, %PARDF
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: cmfmc, pmflxrain, pmflxsnow    ! For computing %CMFMC, %PFLCU, %PFICU
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: dqrcu  ! For computing %DQRCU
        integer, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: cloud_bot ! For computing %DQRCU
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: dtrainu, dtraind, dtrain ! For computing %DTRAIN
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: reevapcn ! For computing %REEVAPCN
        logical, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: cu_act_flag ! For cu_ntiedtke
        real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: raincv_diag, pratec_diag ! For cu_ntiedtke
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: rqccuten_diag, rqicuten_diag, rqvcuten_diag, &
                                                                                          rthcuten_diag, rucuten_diag, rvcuten_diag ! For preventing re-writing the values of grid%
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: taucldc, taucldi, optd !For computing %TAUCLW, %TAUCLI, %OPTD
        logical, optional, intent(in) :: f_qc, f_qi, f_qs ! For computing %TAUCLW, %TAUCLI, %OPTD        

        ! For hybrid-grid in WRF
        real(fp), allocatable :: Ap(:)
        real(fp), allocatable :: Bp(:)

        ! For land use type in WRF
        character*256 :: mminlu_loc

        ! For temporary timers...
        real(KIND=8)          :: WRFGC_Time_Temp_Start, WRFGC_Time_Temp_End

        logical, save :: FIRST = .true.      ! Is first run?

        ! Code begins below...
        call nl_get_debug_level(1, debug_level)
        call set_wrf_debug_level(debug_level)

        ! Get land use map type
        ! USGS or MODIFIED_IGBP_MODIS_NOAH, as validated by chemics_init.F
        call nl_get_mminlu(1, mminlu_loc) 

        IM = ite - its + 1
        JM = jte - jts + 1
        LM = kte - kts + 1

        ! write(6, *) "GIGC_Convert_State_Mod ide =", ide, "jde =", jde
        ! write(6, *) "GIGC_Convert_State_Mod its =", its, "ite =", ite, " IM =", IM
        ! write(6, *) "GIGC_Convert_State_Mod jts =", jts, "jte =", jte, " JM =", JM
        ! write(6, *) "GIGC_Convert_State_Mod kts =", kts, "kte =", kte, " LM =", LM

        ! Start WRF-GC Coupling Physics Timing (Timer-ID: gigc_get_0)
        ! Timer: WRFGC_Phys_Time
        WRFGC_Time_Temp_Start = MPI_Wtime()

        ! Compute required variables from WRF cumulus parameterization
        if(config_flags%cu_physics .eq. 16) then
        
            ! Compute required variables from WRF cumulus parameterization
            ! Preventing rewriting the values of the output from physical cumulus
            ! parameterization
            raincv_diag   = grid%raincv
            pratec_diag   = grid%pratec
            rqccuten_diag = grid%rqccuten
            rqicuten_diag = grid%rqicuten
            rqvcuten_diag = grid%rqvcuten
            rthcuten_diag = grid%rthcuten
            rucuten_diag  = grid%rucuten
            rvcuten_diag  = grid%rvcuten
        
            call cu_ntiedtke(                                         &
                 grid%dt, grid%itimestep, grid%stepcu                 &
                ,raincv_diag, pratec_diag, grid%qfx, grid%hfx         &
                ,grid%u_phy, grid%v_phy, grid%w_2, grid%t_phy         &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qv)            &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qc)            &
                ,moist(grid%sm31,grid%sm32,grid%sm33,p_qi)            &
                ,pi_phy, grid%rho                                     &
                ,grid%rqvften, grid%rthften                           &
                ,dz8w, grid%p_hyd, grid%p_hyd_w, grid%xland           & 
                ,cu_act_flag, grid%dx                                 &
                ,ids, ide, jds, jde, kds, kte+1                       &
                ,grid%sm31, grid%em31, grid%sm33, grid%em33           &
                ,grid%sm32, grid%em32                                 &
                ,its, ite, jts, jte, kts, kte                         &
                ,rthcuten_diag, rqvcuten_diag, rqccuten_diag          &
                ,rqicuten_diag, rucuten_diag, rvcuten_diag            &
                ,pmflxrain, pmflxsnow, cmfmc, dqrcu, cloud_bot        &
                ,dtrainu, dtraind, dtrain, reevapcn                   )
        endif
        
        if(config_flags%ra_lw_physics .eq. 4) then
           call get_cloud_optical_depth(t8w, dz8w                     & 
               ,grid%t_phy, grid%p_hyd_w                              &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qc)             &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qi)             &
               ,moist(grid%sm31,grid%sm32,grid%sm33,p_qs)             &
               ,taucldc, taucldi, optd                                &
               ,f_qc, f_qi, f_qs                                      &
               ,grid%sm31, grid%em31, grid%sm33, grid%em33            &
               ,grid%sm32, grid%em32                                  &
               ,its, ite, jts, jte, kts, kte                          )
        endif
        
        ! See if are using WRF hybrid-grid and diagnose Ap, Bp (hplin, 8/17/18)
        ! FIXME: This calls a lower-level G-C function available in -DMODEL_WRF
        ! Note that Ap is in [hPa]s and you should see a reasonable set of values in pressure_mod
        if(config_flags%hybrid_opt .eq. 2) then
            allocate(Ap(LM+1))
            allocate(Bp(LM+1))

            do k = kts, kte
                Ap(k) = (grid%c4f(k) + (1 - grid%c3f(k)) * grid%p_top) * .01_fp
                Bp(k) = grid%c3f(k)
            enddo

            Ap(LM+1) = grid%p_top / p1000mb * .01_fp
            Bp(LM+1) = 0

            ! Debug the vertical grid (copied from pressure_mod.F, G-C)
            ! WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
            ! WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
            ! WRITE( 6, '(a,/)' ) 'D I A G N O S E D   F R O M   W R F   M O D E L'
            ! WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) Ap(:)
            ! WRITE( 6, '(a)'   )
            ! WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) Bp(:)
            ! WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

            call ACCEPT_EXTERNAL_APBP(am_I_Root, Ap, Bp, GEOS_CHEM_RC)

            deallocate(Ap)
            deallocate(Bp)
        else
            write(6, *) "====================================================="
            write(6, *) "|               W A R N I N G (WRF-GC)              |"
            write(6, *) "====================================================="
            write(6, *) "  GEOS-CHEM USES A VERTICAL HYBRID-SIGMA GRID.       "
            write(6, *) "  WRF MUST BE CONFIGURED TO USE THIS GRID using the  "
            write(6, *) "  namelist option &dynamics: hybrid_opt = 2, and     "
            write(6, *) "  RECOMPILING with ./configure -hyb.                 "
            write(6, *) "                                                     "
            write(6, *) "  WE WERE *NOT* ABLE TO DETECT THIS IN YOUR WRF-GC   "
            write(6, *) "  CONFIGURATION, WHICH MEANS THE VERTICAL LEVELS MAY "
            write(6, *) "  BE INACCURATE AND OUTRIGHT WRONG. PLEASE CHECK.    "
            write(6, *) "====================================================="
        endif

        ! End WRF-GC Coupling Physics Timing (Timer-ID: gigc_get_0)
        ! Timer: WRFGC_Phys_Time
        WRFGC_Time_Temp_End = MPI_Wtime()
        WRFGC_Phys_Time = WRFGC_Phys_Time + (WRFGC_Time_Temp_End - WRFGC_Time_Temp_Start)
        ! End Timer gigc_get_0

        call wrf_debug(100, "GIGC_Convert_State_Mod before State_Met% Conversion")

        !--------------------------------------------------------------
        ! 2-D Variables (Surface Fields)
        !--------------------------------------------------------------
        ! Loop over WRF (IIPAR, JJPAR would be OK too but we are leading with WRF vars on this one)
        do j = jts, jte
        do i = its, ite
            II = i - its + 1
            JJ = j - jts + 1

            ! ALBD              Visible surface albedo [1]
            State_Met%ALBD(II, JJ) = grid%albedo(i, j)

            ! CLDFRC            Column cloud fraction [1]
            State_Met%CLDFRC(II, JJ) = grid%CLDT(i, j)

            ! CLDTOPS           Max cloud top height [levels]
            ! Used for G-C bpch diagnostics, not used for GIGC
            
            ! EFLUX             Latent heat flux [W/m^2]
            State_Met%EFLUX(II, JJ) = grid%LH(i, j)

            ! FRCLND            Olson land fraction [1]
            State_Met%FRCLND(II, JJ) = State_Met%FRLAND(II, JJ)

            ! FRLAND            Fraction of land [1]
            ! FROCEAN           Fraction of ocean [1]
            ! FRLAKE            Fraction of lake [1]
            ! FRSNO             Surface snow fraction [1]
            ! FRLANDIC          Fraction of land ice [1]
            if(grid%lu_mask(i, j) .eq. 0 .and. grid%lakemask(i, j) .eq. 0 .and. grid%snowh(i, j) .eq. 0) then
                 State_Met%FRLAND(II, JJ) = 1.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%lu_mask(i, j) .eq. 1) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 1.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%lakemask(i, j) .eq. 1) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 0.0_fp
                 State_Met%FRLAKE(II, JJ) = 1.0_fp
                 State_Met%FRLANDIC(II, JJ) = 0.0_fp
            else if(grid%snowh(i, j) .ge. 0) then
                 State_Met%FRLAND(II, JJ) = 0.0_fp
                 State_Met%FROCEAN(II, JJ) = 0.0_fp
                 State_Met%FRSNO(II, JJ) = 1.0_fp
                 State_Met%FRLAKE(II, JJ) = 0.0_fp
                 State_Met%FRLANDIC(II, JJ) = 1.0_fp
            endif

            ! FRSEAICE          Surface sea ice fraction [1] | Sea ice fraction [1]
            State_Met%FRSEAICE(II, JJ) = grid%icefrac_gc(i, j)

            ! LandTypeFrac      Olson landmap fraction per type (I,J,TYPE) [1]
            !
            ! XLAI              MODIS LAI per land type for this month, dims I,J,N
            ! Note that the third dimension is actually land type# - 1, as dry dep indexes
            ! start with 0.
            !
            ! IMPORTANT: State_Met%XLAI is average LAI for this surface type
            !  DIVIDED BY THE FRACTIONAL COVERAGE OF THIS SURFACE TYPE. Which means that
            !  it is calculated as
            !    S = 1...%IREG(I,J) (land use type minus 1 for drydep)
            !     State_Met%XLAI(I,J,S) = State_Met%XLAI_NATIVE(I,J,S+1) / landFrac
            !  Where landFrac is State_Met%LandTypeFrac(I,J,S+1)
            !  This means that it is scaled to be larger or equal to the original LAI
            !
            ! For now (hplin, xfeng with discussion 12/1/18) we decide to assign 1.0 fractions to land types
            ! based on grid%lu_index, so we don't have to scale XLAI from LAI for now.
            if(mminlu_loc .eq. 'USGS') then
                select case(int(grid%lu_index(i, j)))
                case (1, 31, 32, 33) ! Urban, residential
                    State_Met%LandTypeFrac(II, JJ, 2 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 2 ) = grid%lai(i, j)
                case (2, 3, 4, 5, 6) ! Dryland Cropland, Cr/Grassland Mosaic, Cr/Woodland Mosaic
                    State_Met%LandTypeFrac(II, JJ, 39) = 1.0_fp ! Cool Irrigated Cropland
                    State_Met%XLAI        (II, JJ, 39) = grid%lai(i, j)
                case (7, 8, 9) ! Grassland, Shrubland, Mixed Sh/Gr
                    State_Met%LandTypeFrac(II, JJ, 43) = 1.0_fp ! Cold Grassland
                    State_Met%XLAI        (II, JJ, 43) = grid%lai(i, j)
                case (10) ! Savanna
                    State_Met%LandTypeFrac(II, JJ, 44) = 1.0_fp ! Woody Savanna
                    State_Met%XLAI        (II, JJ, 44) = grid%lai(i, j)
                case (11) ! Deciduous Broadleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 6 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 6 ) = grid%lai(i, j)
                case (12) ! Deciduous Needleleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 5 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 5 ) = grid%lai(i, j)
                case (13) ! Evergreen Broadleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 7 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 7 ) = grid%lai(i, j)
                case (14) ! Evergreen Needleleaf Forest
                    State_Met%LandTypeFrac(II, JJ, 4 ) = 1.0_fp ! Coniferous Forest
                    State_Met%XLAI        (II, JJ, 4 ) = grid%lai(i, j)
                case (15) ! Mixed forest
                    State_Met%LandTypeFrac(II, JJ, 25) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 25) = grid%lai(i, j)
                case (16) ! Water bodies
                    State_Met%LandTypeFrac(II, JJ, 1 ) = 1.0_fp ! Inland/Sea Water
                    State_Met%XLAI        (II, JJ, 1 ) = grid%lai(i, j)
                case (18) ! Wooded Wetland
                    State_Met%LandTypeFrac(II, JJ, 14) = 1.0_fp ! Wooded Wet Swamp
                    State_Met%XLAI        (II, JJ, 14) = grid%lai(i, j)
                case (19) ! Barren or Sparsely Vegetated
                    State_Met%LandTypeFrac(II, JJ, 3 ) = 1.0_fp ! Low Sparse Grassland
                    State_Met%XLAI        (II, JJ, 3 ) = grid%lai(i, j)
                case (21) ! Wooded Tundra
                    State_Met%LandTypeFrac(II, JJ, 64) = 1.0_fp ! Wooden Tundra
                    State_Met%XLAI        (II, JJ, 64) = grid%lai(i, j)
                case (23) ! Bare Ground Tundra
                    State_Met%LandTypeFrac(II, JJ, 54) = 1.0_fp ! Barren Tundra
                    State_Met%XLAI        (II, JJ, 54) = grid%lai(i, j)
                case (24) ! Snow or Ice
                    State_Met%LandTypeFrac(II, JJ, 13) = 1.0_fp ! Glacier Ice
                    State_Met%XLAI        (II, JJ, 13) = grid%lai(i, j)
                case (25) ! Playa
                    State_Met%LandTypeFrac(II, JJ, 8 ) = 1.0_fp ! Bare Desert
                    State_Met%XLAI        (II, JJ, 8 ) = grid%lai(i, j)
                case default
                    ! default: Urban
                    State_Met%LandTypeFrac(II, JJ, 2 ) = 1.0_fp
                    State_Met%XLAI        (II, JJ, 2 ) = grid%lai(i, j)
                end select
            else if(mminlu_loc .eq. 'MODIFIED_IGBP_MODIS_NOAH') then

            endif

            ! LandTypeFrac      Olson landmap fraction per type (I,J,TYPE) [1]
            ! We translate each of the WRF land use maps to the Olson landmaps used by G-C below
            ! ... in the future, TBD
            ! if(mminlu_loc .eq. 'USGS') then
            !     ! 1 Inland/sea water
            !     !  = 16 Water Bodies
            !     State_Met%LandTypeFrac(II, JJ, 1) = grid%landusef(i, 16, j)

            !     ! 2 Urban
            !     !  = 1 Urban and Built-up Land
            !     State_Met%LandTypeFrac(II, JJ, 2) = grid%landusef(i, 1, j)

            !     ! 3 Low Sparse Grassland

            !     ! 4 Coniferous Forest

            !     ! 5 Deciduous Conifer Forest
            !     !  = 12 Deciduous Needleleaf Forest
            !     State_Met%LandTypeFrac(II, JJ, 5) = grid%landusef(i, 12, j)

            !     ! 6 Deciduous Broadleaf Forest
            !     !  = 11 Deciduous Broadleaf Forest
            !     State_Met%LandTypeFrac(II, JJ, 6) = grid%landusef(i, 11, j)

            !     ! 7 Evergreen Broadleaf Forests
            !     !  = 13 Evergreen Broadleaf Forests
            !     State_Met%LandTypeFrac(II, JJ, 7) = grid%landusef(i, 13, j)

            !     ! 8 Tall Grasses and Shrubs

            !     ! 9 Bare Desert

            !     ! 10 Upland Tundra

            !     ! 11 Irrigated Grassland

            !     ! 12 Semi Desert

            !     ! 13 Glacier Ice

            !     ! 14 Wooded Wet Swamp
            ! else if(mminlu_loc .eq. 'MODIFIED_IGBP_MODIS_NOAH') then
            ! endif
            
            ! GWETROOT          Root soil wetness [1]
            State_Met%GWETROOT(II, JJ) = grid%sm100200(i, j)
            
            ! GWETTOP           Top soil moisture [1]
            State_Met%GWETTOP(II, JJ) = grid%sm000010(i, j)
            
            ! HFLUX             Sensible heat flux [W/m^2]
            State_Met%HFLUX(II, JJ) = grid%hfx(i, j)

            ! LAI               Leaf area index [m^2/m^2] (online)
            State_Met%LAI(II, JJ) = grid%lai(i, j)

            ! MODISLAI          Daily LAI computed from monthly offline MODIS [m2/m2]
            State_Met%MODISLAI(II, JJ) = grid%lai(i, j)

            ! ITY
            ! LWI               Land/water indices [1]
            ! If LWI = 0 and ALBEDO < 0.695, is WATER box.
            !    LWI = 1 .. is LAND box.
            !    LWI = 2 .. is ICE box. (dao_mod.F)
            ! Currently, we are using LWI = 1 as a kludge. Fixme this from WPS data.
            
            if(grid%lu_mask(i, j) .eq. 0) then
                State_Met%LWI(II, JJ) = 1
            else if(grid%lu_mask(i, j) .eq. 1) then
                State_Met%LWI(II, JJ) = 0
            endif
            
            ! PBLH              PBL height [m]
            State_Met%PBLH(II, JJ) = grid%pblh(i, j)

            ! PBL_TOP_L         PBL layer k-index [1] | WRF: kpbl
            State_Met%PBL_TOP_L(II, JJ) = grid%kpbl(i, j)

            ! PHIS              Surface geopotential height [m^2/s^2] | geopotential height at max pressure [m]
            ! State_Met%PHIS(II, JJ) = 9.80665*(grid%ght_max_p(i, j))
            State_Met%PHIS(II, JJ) = grid%ph_2(i, kts, j) + grid%phb(i, kts, j)
            
            ! PRECANV           Anvil precipitation @ ground [kg/m^2/s]
            State_Met%PRECANV(II, JJ) = (grid%snowncv(i, j) + grid%graupelncv(i, j) + grid%hailncv(i, j))/grid%dt
            
            ! PRECCON           Convective precipitation @ ground [kg/m^2/s] | PRECIP rate from cumulus scheme [mm/s]*water density [kg/m^3]
            State_Met%PRECCON(II, JJ) = grid%pratec(i, j)
            
            ! PRECTOT           Total precipitation @ ground [kg/m^2/s] | rain+ice+snow+graupel precipitation rate at all levels [kg/m^2/s]
            State_Met%PRECTOT(II, JJ) = (grid%rainncv(i, j) + grid%snowncv(i, j) + grid%graupelncv(i, j) + grid%hailncv(i, j))/grid%dt + grid%pratec(i, j)
            
            ! PRECLSC           Non-anvil large-scale precipitation @ ground [kg/m^2/s]
            State_Met%PRECLSC(II, JJ) = grid%rainncv(i, j)/grid%dt

            ! PS1_WET           Wet sfc press at dt start [hPa] | surface pressure [Pa]
            ! PS2_WET           Wet sfc press at dt end [hPa]
            ! PSC2_WET          Wet interpolated sfc press [hPa]
            State_Met%PS1_WET(II, JJ) = .01*grid%psfc(i, j)
            State_Met%PS2_WET(II, JJ) = .01*(grid%psfc(i, j) + grid%dpsdt(i, j) * grid%dt)
            State_Met%PSC2_WET(II, JJ) = .01*grid%psfc(i, j)

            ! PS1_DRY           Dry sfc press at dt start [hPa] | Surface pressure - water vapor partial pressure
            ! PS2_DRY           Dry sfc press at dt end [hPa]
            ! PSC2_DRY          Dry interpolated sfc press [hPa]
            !
            ! from http://wiki.seas.harvard.edu/geos-chem/index.php/Moisture_Updates_in_v11-01 (Legacy implementation)
            ! GeosCore/dao_mod.F (New implementation)
            !
            ! FIXME (hplin): As this is online coupling, we do not have access to next dt info.
            ! Possibly port a diagnostics Pa/s tendency from WRFV3/
            ! module_diag_misc? Anyway, if grid%dpsdt (sfc prs tendency) is present, add it to PS2_DRY
            
            ! State_Met%PS1_DRY(II, JJ) = .01*grid%psfc(i, j) - max(moist(i, kts, j, p_qv) * (1._fp/(1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) * 1.6072e-3_fp * grid%psfc(i, j)
            State_Met%PS1_DRY(II, JJ) = .01*grid%psfc(i, j) - 28.97e-3_fp*max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) /&
                                    (18.01e+0_fp + 10.96e-3_fp * max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp)) * grid%psfc(i, j)
            State_Met%PS2_DRY(II, JJ) = State_Met%PS1_DRY(II, JJ) + .01 * grid%dpsdt(i, j) * grid%dt
            ! State_Met%PSC2_DRY(II, JJ) = .01*grid%psfc(i, j) - max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) * 1.6072e-3_fp * grid%psfc(i, j)
            State_Met%PSC2_DRY(II, JJ) = .01*grid%psfc(i, j) - 28.97e-3_fp*max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp) /&
                                    (18.01e+0_fp + 10.96e-3_fp * max(moist(i, kts, j, p_qv) * (1._fp / (1._fp + moist(i, kts, j, p_qv))), 1.0e-30_fp)) * grid%psfc(i, j)
            
            ! RADLWG            Net LW radiation @ ground [W/m^2] | Downward l/w at gnd sfc (WRF)
            ! State_Met%RADLWG(II, JJ) = grid%glw(i, j)
            ! Not used in G-C, comment out from 12.0.0

            ! RADLWG            Net LW radiation @ ground [W/m^2] | Instantaneous downwelling l/w at bottom - upwelling l/w at bottom (ra_lw_physics == 3/4/24)
            !                                                     | Surface l/w downwelling flux - Surface l/w upwelling flux (ra_lw_physics == 5) 
            ! if(config_flags%ra_lw_physics .eq. 3 .or. config_flags%ra_lw_physics .eq. 4 .or. config_flags%ra_lw_physics .eq. 24) then
            !     State_Met%RADLWG(II, JJ) = grid%lwdnb(i, j) - grid%lwupb(i, j)
            ! else if(config_flags%ra_lw_physics .eq. 5) then
            !     State_Met%RADLWG(II, JJ) = grid%slwdn(i, j) - grid%slwup(i, j)
            ! endif
            ! Not used in G-C, comment out from 12.0.0

            ! RADSWG            Solar radiation @ ground [W/m^2] | Downward s/w at gnd sfc (WRF)
            ! State_Met%RADSWG(II, JJ) = grid%swdown(i, j)
            ! Not used in G-C, comment out from 12.0.0

            ! SEAICE00/.../90   Sea Ice Coverage
            ! FIXME (hplin): Set as 1.0 for 00 for now, used only for specialty simulations
            State_Met%SEAICE00(II, JJ) = 1.0_fp

            ! SLP               Sea level pressure [hPa]
            State_Met%SLP(II, JJ) = grid%pslv_gc(i, j) / 100.

            ! SNICE             Fraction of snow/ice [1]

            ! SNODP             Snow depth [m] | Physical snow depth [m]
            State_Met%SNODP(II, JJ) = grid%snowh(i, j)

            ! SNOMAS            Snow mass [kg/m^2] | Snow water equivalent [kg/m^2]
            State_Met%SNOMAS(II, JJ) = grid%acsnow(i, j)
            
            ! SUNCOS            COS(Solar Zenith Angle) current time [1] | COS of solar zenith angle [1]
            State_Met%SUNCOS(II, JJ) = grid%coszen(i, j)
            
            ! SUNCOSmid
            ! FIXME: SUNCOSmid is at midpt of chem timestep, potentially fix by not using %SUNCOS
            call calc_zenithb(grid%xlat(i, j), grid%xlong(i, j)*(-1.), grid%julday, grid%dt, grid%gmt, curr_secs, zen)
            State_Met%SUNCOSmid(II, JJ) = cos(zen)

            ! PARDR             Surface downwelling PAR direct flux [W/m^2] | from WRF module_bioemi_beis314.F:subroutine getpar
            ! PARDF             Surface downwelling PAR diffusive flux [W/m^2] | from WRF module_bioemi_beis314.F:subroutine getpar
            ! Requires zen from above
            call getpar(grid%swdown(i, j), 0.01*(grid%p(i, kts, j)+grid%pb(i, kts, j)), zen, pardb, pardif)
            State_Met%PARDR(II, JJ) = pardb
            State_Met%PARDF(II, JJ) = pardif

            ! SWGDN             Incident radiation @ ground [W/m^2] (RADSWG)
            State_Met%SWGDN(II, JJ) = grid%swdown(i, j)
            
            ! TO3               Total overhead O3 column [DU]
            if (.not. FIRST) then
                State_Met%TO3(II, JJ) = 0.0_fp
                do k = kts, kte
                   State_Met%TO3(II, JJ) =  State_Met%TO3(II, JJ) + chem(i, k, j, p_o3) * State_Met%AD(II, JJ, k) / (grid%dx * grid%dy / grid%msftx(i, j) / grid%msfty(i, j)) * 22.4e+0_fp / 28.97e+0_fp * 0.1_fp
                enddo
            endif

            ! TROPP             Tropopause pressure [hPa] | tropo_p [Pa]
            State_Met%TROPP(II, JJ) = grid%tropo_p(i, j) / 100.

            ! TS                Surface temperature [K]
            State_Met%TS(II, JJ) = grid%t2(i, j)

            ! TSKIN             Surface skin temperature [K]
            State_Met%TSKIN(II, JJ) = grid%tsk(i, j)
            
            ! U10M              E/W wind speed @ 10m height [m/s] | u10
            State_Met%U10M(II, JJ) = grid%u10(i, j)

            ! USTAR             Fraction velocity [m/s] | ust
            State_Met%USTAR(II, JJ) = grid%ust(i, j)
            
            ! UVALBEDO          UV surface albedo [1] | albedo
            State_Met%UVALBEDO(II, JJ) = grid%albedo(i, j)
            
            ! V10M              N/S wind speed @ 10m height [m/s] | v10
            State_Met%V10M(II, JJ) = grid%v10(i, j)

            ! Z0                Surface roughness length [m] | Time varying roughness length [m]
            State_Met%Z0(II, JJ) = grid%znt(i, j)
            
            ! CNV_FRC           Convective fraction, used by ESMF_ mode [1]

            !----------------------------------------------------------
            ! 3-D variables filled on bottom/top levels
            !----------------------------------------------------------
            ! Fill top levels
            State_Met%PFLCU(II, JJ, LM+1) = 0._fp
            State_Met%PFICU(II, JJ, LM+1) = 0._fp
            State_Met%PFLLSAN(II, JJ, LM+1) = 0._fp
            State_Met%PFILSAN(II, JJ, LM+1) = 0._fp

            State_Met%PEDGE(II, JJ, kts) = 0.01_fp*grid%psfc(i, j)
!            State_Met%PEDGE(II, JJ, 1) = .01*((grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 3, j) - grid%phb(i, 3, j))/                                         &
!                                               (grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))*(grid%p(i, 1, j) + grid%pb(i, 1, j)) +  &
!                                               (grid%ph_2(i, 3, j) + grid%phb(i, 3, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))/                                        &
!                                               (grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))*(grid%p(i, 2, j) + grid%pb(i, 2, j)))

            State_Met%PEDGE(II, JJ, LM + 1) = .01*(exp((grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 3, j) - grid%phb(i, 3, j))/                                &
                                                      (grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))*                                 &
                                                   log(grid%p(i, kte, j) + grid%pb(i, kte, j)) +                                                                         &
                                                      (grid%ph_2(i, 3, j) + grid%phb(i, 3, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))/                                 &
                                                      (grid%ph_2(i, 1, j) + grid%phb(i, 1, j) - grid%ph_2(i, 2, j) - grid%phb(i, 2, j))*                                 &
                                                   log(grid%p(i, kte-1, j) + grid%pb(i, kte-1, j))))
        enddo
        enddo

        call wrf_debug(100, "GIGC_Convert_State_Mod after 2-D State_Met% Conversion")

        ! Loop over WRF tile but in GEOS-Chem order as we are writing to GEOS-Chem
        ! Remember, for 3-D looping over (I, J, K) loops should be written outermost to innermost
        ! to make the most efficient use of the column-first order in Fortran
        !
        ! Top levels and lowermost levels (edge conditions) are now filled in the 2-D loop that
        ! does not loop over k
        do k = kts, kte
        do j = jts, jte
        do i = its, ite
            II = i - its + 1
            JJ = j - jts + 1

            !--------------------------------------------------------------
            ! 3-D Variables
            !
            ! For 3-D Variables, there is one unified loop that performs
            ! the conversions over vertical levs.
            !--------------------------------------------------------------

            ! AREA_M2           Grid box surface area [m^2]
            State_Met%AREA_M2(II, JJ, 1) = grid%dx * grid%dy / grid%msftx(i, j) / grid%msfty(i, j)

            ! PEDGE             Wet air pressure at (vertial) level edges [hPa]
            if(k .ne. kts) then ! Start from kts + 1
                State_Met%PEDGE(II, JJ, k) = .01*(grid%fnm(k)*(grid%p(i, k, j) + grid%pb(i, k, j)) + grid%fnp(k)*(grid%p(i, k-1, j) + grid%pb(i, k-1, j)))
            endif

            ! CLDF              3-D cloud fraction [1]
            State_Met%CLDF(II, JJ, k) = grid%cldfra(i, k, j)

            ! CMFMC                               Upward moist convective mass flux [kg/m^2/s] (max:LM+1)
            if(config_flags%cu_physics .eq. 10) then
                State_Met%CMFMC(II, JJ, k) = grid%mfup_cup(i, k, j)
            else if(config_flags%cu_physics .eq. 7) then
                State_Met%CMFMC(II, JJ, k) = grid%cmfmcdzm(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%CMFMC(II, JJ, k) = cmfmc(i, k, j)
            endif

            ! T                 Temperature [K]
            ! These calculations are linked, do not change their order and place them first
            State_Met%T(II, JJ, k) = (grid%t_2(i, k, j) + t0)*((grid%p(i, k, j) + grid%pb(i, k, j))/p1000mb)**rcp

            ! DQRCU             Convective precipitation production rate [kg/kg(dry air)/s] | dq/dt due to deep convective rainout [kg/kg/s]
            if(config_flags%cu_physics .eq. 7) then
                State_Met%DQRCU(II, JJ, k) = grid%zmntprpd(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%DQRCU(II, JJ, k) = dqrcu(i, k, j)
            endif
            
            ! DQRLSAN           Large scale precipitation production rate [kg/kg(dry air)/s] | Rate of conversion of condensate to precipitation [kg/kg/s]
            if(config_flags%mp_physics .eq. 10 .or. config_flags%mp_physics .eq. 8) then
                State_Met%DQRLSAN(II, JJ, k) = grid%rainprod(i, k, j)
            else if(config_flags%mp_physics .eq. 11) then
                State_Met%DQRLSAN(II, JJ, k) = grid%prain3d(i, k, j)
            endif
            
            ! DTRAIN            Detrainment flux [kg/m^2/s]
            if(config_flags%cu_physics .eq. 7) then
                State_Met%DTRAIN(II, JJ, k) = grid%du3d(i, k, j)*grid%dp3d(i, k, j)*100._fp/9.8_fp
                if(State_Met%DTRAIN(II, JJ, k) > 0.5) then
                    State_Met%DTRAIN(II, JJ, k) = 0.5
                endif
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%DTRAIN(II, JJ, k) = dtrain(i, k, j)
            endif

            ! OMEGA                  Updraft velocity [Pa/s]
            if(State_Met%AIRVOL(II, JJ, k) .ne. 0) then
                State_Met%OMEGA(II, JJ, k) = -9.8_fp*State_Met%AD(II, JJ, k)/State_Met%AIRVOL(II, JJ, k)*grid%w_2(i, k, j)
            endif

            ! QI                Ice mixing ratio [kg/kg dry air] || cloud ice mixing ratio [kg/kg]
            ! QL                Water mixing ratio [kg/kg dry air] || cloud water mixing ratio [kg/kg]
            ! QI, QL (WRF registry ikjftb on moist)
            State_Met%QI(II, JJ, k) = moist(i, k, j, p_qi)
            State_Met%QL(II, JJ, k) = moist(i, k, j, p_qc)

            ! REEVAPCN
            if(config_flags%cu_physics .eq. 7) then
                State_Met%REEVAPCN(II, JJ, k) = abs(grid%evapcdp3d(i, k, j))
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%REEVAPCN(II, JJ, k) = abs(reevapcn(i, k, j)*State_Met%AD(II, JJ, k)/State_Met%AREA_M2(II, JJ, k))
            endif
            
            ! REEVAPLS
            if(config_flags%mp_physics .eq. 10) then
                State_Met%REEVAPLS(II, JJ, k) = abs(grid%evapprod(i, k, j))
            else if(config_flags%mp_physics .eq. 11) then
                State_Met%REEVAPLS(II, JJ, k) = abs(grid%nevapr3d(i, k, j))
            endif

            ! RH                Relative humidity [%]
            State_Met%RH(II, JJ, k) = max(.1, min(.95, moist(i, k, j, p_qv)/ &
                                   (3.80*exp(17.27*(State_Met%T(II, JJ, k) - 273.)/ &
                                             (State_Met%T(II, JJ, k) - 36.))/(.01*(grid%p(i, k, j) + grid%pb(i, k, j)))))) * 100

            ! SPHU              Specific humidity [g H2O/kg tot air]
            ! SPHU1             ... at ts start [g/kg]
            ! SPHU2             ... at ts end [g/kg]
            ! Compute specific humidities based on moisture data
            ! from module_cam_mam_gas_wetdep_driver:168
            State_Met%SPHU(II, JJ, k) = max(moist(i, k, j, p_qv) * (1._fp / (1._fp + moist(i, k, j, p_qv))), 1.0e-30_fp)
            State_Met%SPHU1(II, JJ, k) = State_Met%SPHU(II, JJ, k)
            if(config_flags%cu_physics .eq. 7) then
                State_Met%SPHU2(II, JJ, k) = State_Met%SPHU(II, JJ, k) + grid%zmdq(i, k, j)*1.0e3_fp*grid%dt
            else ! Unsupported physics option
                State_Met%SPHU2(II, JJ, k) = State_Met%SPHU(II, JJ, k)
            endif

            ! U                 E/W component of wind [m/s] (from WRF U - X cmp)
            ! V                 N/S component of wind [m/s] (from WRF V - X cmp)
            if(i .ne. ide) then
                State_Met%U(II, JJ, k) = 0.5*(grid%u_2(i, k, j) + grid%u_2(i + 1, k, j))
            else
                State_Met%U(II, JJ, k) = grid%u_2(i, k, j)
            endif

            if(j .ne. jde) then
                State_Met%V(II, JJ, k) = 0.5*(grid%v_2(i, k, j) + grid%v_2(i, k, j + 1))
            else
                State_Met%V(II, JJ, k) = grid%v_2(i, k, j)
            endif

            ! TAUCLI                 Optical depth of ice clouds [1] | Cloud optical thickness for ice [1]
            ! TAUCLW                 Optical depth of H20 clouds [1] | Cloud optical thickness for water [1]
            ! OPTD                   Total in-cloud optical depth (visible band) 
            if(config_flags%ra_lw_physics .eq. 3 .or. config_flags%ra_lw_physics .eq. 5) then
                State_Met%TAUCLI(II, JJ, k) = grid%taucldi(i, k, j)
                State_Met%TAUCLW(II, JJ, k) = grid%taucldc(i, k, j)
                State_Met%OPTD(II, JJ, k) = State_Met%TAUCLI(II, JJ, k) + State_Met%TAUCLW(II, JJ, k)
            else if(config_flags%ra_lw_physics .eq. 4) then
                State_Met%TAUCLI(II, JJ, k) = taucldi(i, k, j)
                State_Met%TAUCLW(II, JJ, k) = taucldc(i, k, j)
                State_Met%OPTD(II, JJ, k) = optd(i, k, j)
            endif

            ! PFICU/PFLCU                         Downward flux of ice/liquid precipitation (convective) [kg/m^2/s] 
            if(config_flags%cu_physics .eq. 7) then
                State_Met%PFLCU(II, JJ, k) = grid%zmflxprc(i, k, j)
                State_Met%PFICU(II, JJ, k) = grid%zmflxsnw(i, k, j)
            else if(config_flags%cu_physics .eq. 16) then
                State_Met%PFLCU(II, JJ, k) = pmflxrain(i, k, j)
                State_Met%PFICU(II, JJ, k) = pmflxsnow(i, k, j)
            endif

            ! PFILSAN/PFLLSAN                    Downward flux of ice/liquid precipitation (Large-scale & Anvil) [kg/m^2/s]
            State_Met%PFLLSAN(II, JJ, k) = grid%precr(i, k, j)
            State_Met%PFILSAN(II, JJ, k) = grid%preci(i, k, j) + grid%precs(i, k, j) + grid%precg(i, k, j)

            ! TMPU1                              Temperature at start of timestep [K]
            ! TMPU2                              Temperature at end of timestep [K]
            ! Not used for main G-C, only in Met Fields by FlexGrid, but for safety set to %T
            State_Met%TMPU1(II, JJ, k) = State_Met%T(II, JJ, k)
            State_Met%TMPU2(II, JJ, k) = State_Met%T(II, JJ, k)

            ! UPDVVEL                            Updraft vertical velocity [hPa/s]
            if(State_Met%AIRVOL(II, JJ, k) .ne. 0) then
                State_Met%UPDVVEL(II, JJ, k) = -9.8_fp*State_Met%AD(II, JJ, k)/State_Met%AIRVOL(II, JJ, k)*grid%w_2(i, k, j)
            endif

            ! The following code simply updates WRF diagnostics data for testing...
            grid%diagikj0(i, k, j) = State_Met%CMFMC(II, JJ, k)
            grid%diagikj1(i, k, j) = State_Met%DTRAIN(II, JJ, k)
            grid%diagikj2(i, k, j) = State_Met%PFLCU(II, JJ, k)
            grid%diagikj3(i, k, j) = State_Met%PFLLSAN(II, JJ, k)
            grid%diagikj4(i, k, j) = State_Met%DQRCU(II, JJ, k)
            grid%diagikj5(i, k, j) = State_Met%AD(II, JJ, k)
            grid%diagikj6(i, k, j) = State_Met%AIRVOL(II, JJ, k)
            grid%diagikj7(i, k, j) = State_Met%UPDVVEL(II, JJ, k)
            grid%diagij0(i, j) = State_Met%PRECTOT(II, JJ)
            grid%diagij1(i, j) = State_Met%PRECCON(II, JJ)
            grid%diagij2(i, j) = State_Met%PRECLSC(II, JJ)
            grid%diagij3(i, j) = State_Met%TO3(II, JJ)

            ! This code ONLY runs after all levels for given grid box are filled (k = LM)
            if(k .eq. kte) then
                ! A kludge that stops GEOS-CHEM CGM ("ChemGrid_Mod") GET_TPAUSE_LEVEL from breaking,
                ! when the T-PAUSE could not be found.
                ! WRF has a weaker stratosphere routine, so some grid boxes might NOT have a T-PAUSE at
                ! all. In this case, force CGM to recognize LLPAR (top level, LM) as T-PAUSE.
                ! Is T_PAUSE cond: State_Met%TROPP(I, J) >= State_Met%PEDGE(I, J, L_CURRENT)
                if(minval(State_Met%PEDGE(II, JJ, 1:LM)) .ge. State_Met%TROPP(II, JJ)) then
                    State_Met%TROPP(II, JJ) = minval(State_Met%PEDGE(II, JJ, 1:LM))
                    write(6, *) "GIGC_Convert_State_Mod fixed T-PAUSE for GEOS-CHEM-CGM, set T-PAUSE to LLPAR"
                endif
            endif

            ! Update chemistry concentrations in GEOS-Chem from WRF IC values
            if(.not. FIRST) then
                State_Chm%Species(II, JJ, k, gi_no) = chem(i, k, j, p_no) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_o3) = chem(i, k, j, p_o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pan) = chem(i, k, j, p_pan) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_co) = chem(i, k, j, p_co) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_alk4) = chem(i, k, j, p_alk4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isop) = chem(i, k, j, p_isoprene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno3) = chem(i, k, j, p_hno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o2) = chem(i, k, j, p_h2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acet) = chem(i, k, j, p_acet) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mek) = chem(i, k, j, p_mek) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ald2) = chem(i, k, j, p_ald2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_rcho) = chem(i, k, j, p_rcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvk) = chem(i, k, j, p_mvk) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macr) = chem(i, k, j, p_macr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_npmn) = chem(i, k, j, p_npmn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ppn) = chem(i, k, j, p_ppn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_r4n2) = chem(i, k, j, p_r4n2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_prpe) = chem(i, k, j, p_prpe) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c3h8) = chem(i, k, j, p_c3h8) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2o) = chem(i, k, j, p_hcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c2h6) = chem(i, k, j, p_c2h6) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o5) = chem(i, k, j, p_n2o5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno4) = chem(i, k, j, p_hno4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mp) = chem(i, k, j, p_ch3ooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dms) = chem(i, k, j, p_dms) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so2) = chem(i, k, j, p_so2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4) = chem(i, k, j, p_so4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4s) = chem(i, k, j, p_so4s) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_msa) = chem(i, k, j, p_msa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh3) = chem(i, k, j, p_nh3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh4) = chem(i, k, j, p_nh4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nit) = chem(i, k, j, p_nit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nits) = chem(i, k, j, p_nits) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpi) = chem(i, k, j, p_bcpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpi) = chem(i, k, j, p_ocpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpo) = chem(i, k, j, p_bcpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpo) = chem(i, k, j, p_ocpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst1) = chem(i, k, j, p_dst1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst2) = chem(i, k, j, p_dst2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst3) = chem(i, k, j, p_dst3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst4) = chem(i, k, j, p_dst4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_sala) = chem(i, k, j, p_sala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_salc) = chem(i, k, j, p_salc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_br2) = chem(i, k, j, p_br2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_br) = chem(i, k, j, p_br) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bro) = chem(i, k, j, p_bro) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hobr) = chem(i, k, j, p_hobr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hbr) = chem(i, k, j, p_hbr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brno2) = chem(i, k, j, p_brno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brno3) = chem(i, k, j, p_brno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_chbr3) = chem(i, k, j, p_chbr3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2br2) = chem(i, k, j, p_ch2br2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3br) = chem(i, k, j, p_ch3br) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mpn) = chem(i, k, j, p_mpn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isopnd) = chem(i, k, j, p_isopnd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isopnb) = chem(i, k, j, p_isopnb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_moba) = chem(i, k, j, p_moba) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_propnn) = chem(i, k, j, p_propnn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hac) = chem(i, k, j, p_hac) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_glyc) = chem(i, k, j, p_glyc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvkn) = chem(i, k, j, p_mvkn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macrn) = chem(i, k, j, p_macrn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_map) = chem(i, k, j, p_map) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no2) = chem(i, k, j, p_no2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no3) = chem(i, k, j, p_no3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno2) = chem(i, k, j, p_hno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_benz) = chem(i, k, j, p_benzene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_tolu) = chem(i, k, j, p_toluene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_xyle) = chem(i, k, j, p_xylenes) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mtpa) = chem(i, k, j, p_mtpa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_limo) = chem(i, k, j, p_limon) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mtpo) = chem(i, k, j, p_mtpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soap) = chem(i, k, j, p_soap) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soas) = chem(i, k, j, p_soas) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_eoh) = chem(i, k, j, p_eoh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mgly) = chem(i, k, j, p_mgly) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_glyx) = chem(i, k, j, p_glyx) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acta) = chem(i, k, j, p_acta) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hpald) = chem(i, k, j, p_hpald) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dhdn) = chem(i, k, j, p_dhdn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ethln) = chem(i, k, j, p_ethln) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcooh) = chem(i, k, j, p_hcooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxa) = chem(i, k, j, p_iepoxa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxb) = chem(i, k, j, p_iepoxb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iepoxd) = chem(i, k, j, p_iepoxd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1) = chem(i, k, j, p_isn1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripa) = chem(i, k, j, p_ripa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripb) = chem(i, k, j, p_ripb) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ripd) = chem(i, k, j, p_ripd) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_imae) = chem(i, k, j, p_imae) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soaie) = chem(i, k, j, p_soaie) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soame) = chem(i, k, j, p_soame) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soagx) = chem(i, k, j, p_soagx) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_soamg) = chem(i, k, j, p_soamg) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_lvoc) = chem(i, k, j, p_lvoc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_lvocoa) = chem(i, k, j, p_lvocoa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1og) = chem(i, k, j, p_isn1og) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isn1oa) = chem(i, k, j, p_isn1oa) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monits) = chem(i, k, j, p_monits) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monitu) = chem(i, k, j, p_monitu) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_honit) = chem(i, k, j, p_honit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ionita) = chem(i, k, j, p_ionita) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_monita) = chem(i, k, j, p_monita) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_indiol) = chem(i, k, j, p_indiol) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ipmn) = chem(i, k, j, p_ipmn) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hc187) = chem(i, k, j, p_hc187) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o) = chem(i, k, j, p_n2o) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocs) = chem(i, k, j, p_ocs) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch4) = chem(i, k, j, p_ch4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brcl) = chem(i, k, j, p_brcl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcl) = chem(i, k, j, p_hcl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ccl4) = chem(i, k, j, p_ccl4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3cl) = chem(i, k, j, p_ch3cl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3ccl3) = chem(i, k, j, p_ch3ccl3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc113) = chem(i, k, j, p_cfc113) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc114) = chem(i, k, j, p_cfc114) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc115) = chem(i, k, j, p_cfc115) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc123) = chem(i, k, j, p_hcfc123) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc141b) = chem(i, k, j, p_hcfc141b) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc142b) = chem(i, k, j, p_hcfc142b) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc11) = chem(i, k, j, p_cfc11) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cfc12) = chem(i, k, j, p_cfc12) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hcfc22) = chem(i, k, j, p_hcfc22) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h1211) = chem(i, k, j, p_h1211) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h1301) = chem(i, k, j, p_h1301) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2402) = chem(i, k, j, p_h2402) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl) = chem(i, k, j, p_cl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clo) = chem(i, k, j, p_clo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hocl) = chem(i, k, j, p_hocl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clno3) = chem(i, k, j, p_clno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_clno2) = chem(i, k, j, p_clno2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cloo) = chem(i, k, j, p_cloo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oclo) = chem(i, k, j, p_oclo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl2) = chem(i, k, j, p_cl2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_cl2o2) = chem(i, k, j, p_cl2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o) = chem(i, k, j, p_h2o) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brsala) = chem(i, k, j, p_brsala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_brsalc) = chem(i, k, j, p_brsalc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_chcl3) = chem(i, k, j, p_chcl3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2cl2) = chem(i, k, j, p_ch2cl2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch3i) = chem(i, k, j, p_ch3i) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2i2) = chem(i, k, j, p_ch2i2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2icl) = chem(i, k, j, p_ch2icl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2ibr) = chem(i, k, j, p_ch2ibr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hoi) = chem(i, k, j, p_hoi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2) = chem(i, k, j, p_i2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ibr) = chem(i, k, j, p_ibr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_icl) = chem(i, k, j, p_icl) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i) = chem(i, k, j, p_i) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_io) = chem(i, k, j, p_io) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hi) = chem(i, k, j, p_hi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oio) = chem(i, k, j, p_oio) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ino) = chem(i, k, j, p_ino) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iono) = chem(i, k, j, p_iono) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_iono2) = chem(i, k, j, p_iono2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o2) = chem(i, k, j, p_i2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o3) = chem(i, k, j, p_i2o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_i2o4) = chem(i, k, j, p_i2o4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isala) = chem(i, k, j, p_isala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isalc) = chem(i, k, j, p_isalc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_aeri) = chem(i, k, j, p_aeri) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pfe) = chem(i, k, j, p_pfe) * 1.0e-6_fp

                ! Add mercury, yxzhang, 06/20/2019
                State_Chm%Species(II, JJ, k, gi_hg0) = chem(i, k, j, p_hg0) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hg2) = chem(i, k, j, p_hg2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hgp) = chem(i, k, j, p_hgp) * 1.0e-6_fp

                ! A lot of non-advected species were removed here
                ! as we do NOT want them to be advected by WRF, so their concentrations
                ! should NOT be updated from WRF values. (hplin, 8/23/18)
            endif
        enddo ! Loop over i
        enddo ! Loop over j
        enddo ! Loop over k

        call wrf_debug(100, "GIGC_Convert_State_Mod after 3-D State_Met% Conversion")

        ! Unless we have existing initial / boundary conditions passed in from wrfinput, then update
        ! GEOS-Chem concentrations with existing ones.
        ! Note: Do not merge this with the below as species are limited here.
        if(config_flags%have_bcs_chem .and. FIRST) then
            call wrf_debug(1, "GIGC_Convert_State_Mod updating chemistry concentrations in GEOS-Chem from WRF IC values")
            do k = kts, kte
            do j = jts, jte
            do i = its, ite
                II = i - its + 1
                JJ = j - jts + 1

                State_Chm%Species(II, JJ, k, gi_no) = chem(i, k, j, p_no) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_o3) = chem(i, k, j, p_o3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_pan) = chem(i, k, j, p_pan) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_co) = chem(i, k, j, p_co) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_alk4) = chem(i, k, j, p_alk4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_isop) = chem(i, k, j, p_isoprene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno3) = chem(i, k, j, p_hno3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_h2o2) = chem(i, k, j, p_h2o2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acet) = chem(i, k, j, p_acet) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mek) = chem(i, k, j, p_mek) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ald2) = chem(i, k, j, p_ald2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mvk) = chem(i, k, j, p_mvk) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_macr) = chem(i, k, j, p_macr) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c3h8) = chem(i, k, j, p_c3h8) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch2o) = chem(i, k, j, p_hcho) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_c2h6) = chem(i, k, j, p_c2h6) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_n2o5) = chem(i, k, j, p_n2o5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hno4) = chem(i, k, j, p_hno4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mp) = chem(i, k, j, p_ch3ooh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dms) = chem(i, k, j, p_dms) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so2) = chem(i, k, j, p_so2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_so4) = chem(i, k, j, p_so4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh3) = chem(i, k, j, p_nh3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nh4) = chem(i, k, j, p_nh4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no2) = chem(i, k, j, p_no2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_no3) = chem(i, k, j, p_no3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_tolu) = chem(i, k, j, p_toluene) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_limo) = chem(i, k, j, p_limon) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_mgly) = chem(i, k, j, p_mgly) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ch4) = chem(i, k, j, p_ch4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hc5) = chem(i, k, j, p_hc5) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ho2) = chem(i, k, j, p_ho2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_oh) = chem(i, k, j, p_oh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_prpe) = chem(i, k, j, p_prpe) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_nit) = chem(i, k, j, p_nit) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpi) = chem(i, k, j, p_bcpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpi) = chem(i, k, j, p_ocpi) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_bcpo) = chem(i, k, j, p_bcpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_ocpo) = chem(i, k, j, p_ocpo) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst1) = chem(i, k, j, p_dst1) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst2) = chem(i, k, j, p_dst2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst3) = chem(i, k, j, p_dst3) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_dst4) = chem(i, k, j, p_dst4) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_sala) = chem(i, k, j, p_sala) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_salc) = chem(i, k, j, p_salc) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hac) = chem(i, k, j, p_hac) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_map) = chem(i, k, j, p_map) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_eoh) = chem(i, k, j, p_eoh) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_acta) = chem(i, k, j, p_acta) * 1.0e-6_fp
                ! add mercury, yxzhang, 06/20/2019
                State_Chm%Species(II, JJ, k, gi_hg0) = chem(i, k, j, p_hg0) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hg2) = chem(i, k, j, p_hg2) * 1.0e-6_fp
                State_Chm%Species(II, JJ, k, gi_hgp) = chem(i, k, j, p_hgp) * 1.0e-6_fp
            enddo
            enddo
            enddo
        endif

        call wrf_debug(100, "GIGC_Convert_State_Mod after State Conversion")

        if(debug_level .ge. 999 .and. FIRST) then
            ! Debug State_Chm
            write(6, *) "%%%%%% GIGC_Convert_State_Mod Chemistry State Export %%%%%%"
            write(6, *) "Units: ", State_Chm%Spc_Units
            do N = 1, State_Chm%nSpecies
                WRITE(6, *) "N:", N, State_Chm%SpcData(N)%Info%Name, &
                            "Is_Adv:", State_Chm%SpcData(N)%Info%Is_Advected
            enddo
        endif

        FIRST = .false.

    end subroutine GIGC_Get_WRF
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_Set_WRF
!
! !DESCRIPTION: Subroutine GIGC\_Set\_WRF converts GIGC data structures to
!  WRF data structures (grid, chem, moist, scalar) so that they are saved into
!  wrfout outputs and advected by the WRF dynamics core. (hplin, 4/24/18)
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_Set_WRF(am_I_Root, &
        config_flags, grid, &
        num_chem, chem, num_scalar, scalar, num_moist, moist, &
        its, ite, jts, jte, &
        ide, jde, &
        kts, kte, &
        Input_Opt, State_Met, State_Chm, State_Diag)
        use aerosol_mod, only:   PM25
        implicit none

!
! !INPUT/OUTPUT PARAMETERS:
!
        logical, intent(in) :: am_I_Root
        type(grid_config_rec_type), intent(in) :: config_flags
        type(domain), target :: grid

        integer, intent(in) :: num_chem, num_scalar, num_moist ! # of WRF chem, moist variables. Check "dummy_new_args.inc" & "actual_new_args.inc" in WRF compile.
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem), intent(inout)   :: chem
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_scalar), intent(inout) :: scalar
        real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_moist), intent(inout)  :: moist

        integer, intent(in) :: its, ite, jts, jte, ide, jde, kts, kte

        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm
        type(DgnState), intent(inout) :: State_Diag
!
! !REMARKS:
!  Please do take note that in WRF, memorder is XZY so arrays are (i, k, j), whereas in
!  GEOS-Chem the memorder is XYZ.
!
!  Also, the State_Chm conversion routines are auto-generated. Please use
!  https://github.com/jimmielin/wrf-gchp-dev-tools to generate them.
!  (hplin, 5/17/18)
!
! !REVISION HISTORY:
!  17 May 2018 - H.P. Lin   - Added ProTeX headers
!  12 Oct 2018 - H.P. Lin   - Add State_Diag for G-C netCDF diagnostics into WRF-GC output (wrfout)
!  13 Nov 2018 - H.P. Lin   - Add a temporarily kludge for MPI striping; merge in aerosol pm25
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: debug_level      ! WRF Debug Level
        integer :: N

        integer :: IM, II           ! # of lons on this PET
        integer :: JM, JJ           ! # of lats on this PET
        integer :: LM               ! # of levs on this PET
        integer :: i, j, k          ! Loop variable over xts, xte (x = i, j, k)

        ! Code below
        call nl_get_debug_level(1, debug_level)
        call set_wrf_debug_level(debug_level)

        ! # of lon, lat, lev in this PET
        IM = ite - its + 1
        JM = jte - jts + 1
        LM = kte - kts + 1

        ! write(6, *) "GIGC_Convert_State_Mod::Set grid%id", grid%id
        ! write(6, *) "GIGC_Convert_State_Mod::Set ide =", ide, "jde =", jde
        ! write(6, *) "GIGC_Convert_State_Mod::Set its =", its, "ite =", ite, " IM =", IM
        ! write(6, *) "GIGC_Convert_State_Mod::Set jts =", jts, "jte =", jte, " JM =", JM
        ! write(6, *) "GIGC_Convert_State_Mod::Set kts =", kts, "kte =", kte, " LM =", LM

        ! write(6, *) "Sanity check: State_Chm% species sizing"
        ! write(6, *) "IM =", size(State_Chm%Species, 1)
        ! write(6, *) "JM =", size(State_Chm%Species, 2)
        ! write(6, *) "LM =", size(State_Chm%Species, 3)
        ! write(6, *) "SPC =", size(State_Chm%Species, 4)

        ! write(6, *) "Sanity check: State_Met%PEDGE sizing"
        ! write(6, *) "IM =", size(State_Met%PEDGE, 1)
        ! write(6, *) "JM =", size(State_Met%PEDGE, 2)
        ! write(6, *) "LM+1 =", size(State_Met%PEDGE, 3)

        ! When doing State_Chm% conversions, you can use State_Chm_Mod::IND_
        ! to retrieve a species ID for N value, e.g. id_O3 = IND_('O3'), so
        ! State_Chm%Species(:, :, :, id_O3) holds the O3 information.

        ! Fix negatives
        where(State_Chm%Species < 0.0e0)
            State_Chm%Species = 1.0e-36
        endwhere

        call wrf_debug(1, "GIGC_Convert_State_Mod updating chemistry concs/diags in WRF from GEOS-Chem values")

        ! For now, "WRF-GC" style diagnostics from State_Diag are ALWAYS done.
        ! If you use gc_do_gcdiagn, it only turns on/off the G-C/C GCHP style HISTORY.rc diagns,
        ! it does not affect this.
        ! Loop over WRF tile
        do j = jts, jte
        do k = kts, kte
        do i = its, ite
            ! write(6, *) "Debug GIGC_Convert_State_Mod::Set loop i, its", i, its
            ! write(6, *) "Debug GIGC_Convert_State_Mod::Set loop j, jts", j, jts

            II = i - its + 1
            JJ = j - jts + 1

            ! Loop over WRF levels (same as GEOS-Chem levels since there is no global indexing here)
            ! write(6, *) "Debug GIGC_Convert_State_Mod::Set in k loop IIJJk", II, JJ, k
            ! Some species are particularly problematic, for some reason, in wet-deposition RESUSPENSION
            ! in middle levs. We patch these values out for known species so WetDep does not
            ! complain about them...
            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 7))) then
                State_Chm%Species(II, JJ, k, 8) = 1.0e-36_fp
            endif

            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 8))) then
                State_Chm%Species(II, JJ, k, 8) = 1.0e-36_fp
            endif

            if(IT_IS_NAN(State_Chm%Species(II, JJ, k, 26))) then
                State_Chm%Species(II, JJ, k, 26) = 1.0e-36_fp
            endif

            chem(i, k, j, p_no) = State_Chm%Species(II, JJ, k, gi_no) * 1.0e+6_fp
            chem(i, k, j, p_o3) = State_Chm%Species(II, JJ, k, gi_o3) * 1.0e+6_fp
            chem(i, k, j, p_pan) = State_Chm%Species(II, JJ, k, gi_pan) * 1.0e+6_fp
            chem(i, k, j, p_co) = State_Chm%Species(II, JJ, k, gi_co) * 1.0e+6_fp
            chem(i, k, j, p_alk4) = State_Chm%Species(II, JJ, k, gi_alk4) * 1.0e+6_fp
            chem(i, k, j, p_isoprene) = State_Chm%Species(II, JJ, k, gi_isop) * 1.0e+6_fp
            chem(i, k, j, p_hno3) = State_Chm%Species(II, JJ, k, gi_hno3) * 1.0e+6_fp
            chem(i, k, j, p_h2o2) = State_Chm%Species(II, JJ, k, gi_h2o2) * 1.0e+6_fp
            chem(i, k, j, p_acet) = State_Chm%Species(II, JJ, k, gi_acet) * 1.0e+6_fp
            chem(i, k, j, p_mek) = State_Chm%Species(II, JJ, k, gi_mek) * 1.0e+6_fp
            chem(i, k, j, p_ald2) = State_Chm%Species(II, JJ, k, gi_ald2) * 1.0e+6_fp
            chem(i, k, j, p_rcho) = State_Chm%Species(II, JJ, k, gi_rcho) * 1.0e+6_fp
            chem(i, k, j, p_mvk) = State_Chm%Species(II, JJ, k, gi_mvk) * 1.0e+6_fp
            chem(i, k, j, p_macr) = State_Chm%Species(II, JJ, k, gi_macr) * 1.0e+6_fp
            chem(i, k, j, p_npmn) = State_Chm%Species(II, JJ, k, gi_npmn) * 1.0e+6_fp
            chem(i, k, j, p_ppn) = State_Chm%Species(II, JJ, k, gi_ppn) * 1.0e+6_fp
            chem(i, k, j, p_r4n2) = State_Chm%Species(II, JJ, k, gi_r4n2) * 1.0e+6_fp
            chem(i, k, j, p_prpe) = State_Chm%Species(II, JJ, k, gi_prpe) * 1.0e+6_fp
            chem(i, k, j, p_c3h8) = State_Chm%Species(II, JJ, k, gi_c3h8) * 1.0e+6_fp
            chem(i, k, j, p_hcho) = State_Chm%Species(II, JJ, k, gi_ch2o) * 1.0e+6_fp
            chem(i, k, j, p_c2h6) = State_Chm%Species(II, JJ, k, gi_c2h6) * 1.0e+6_fp
            chem(i, k, j, p_n2o5) = State_Chm%Species(II, JJ, k, gi_n2o5) * 1.0e+6_fp
            chem(i, k, j, p_hno4) = State_Chm%Species(II, JJ, k, gi_hno4) * 1.0e+6_fp
            chem(i, k, j, p_ch3ooh) = State_Chm%Species(II, JJ, k, gi_mp) * 1.0e+6_fp
            chem(i, k, j, p_dms) = State_Chm%Species(II, JJ, k, gi_dms) * 1.0e+6_fp
            chem(i, k, j, p_so2) = State_Chm%Species(II, JJ, k, gi_so2) * 1.0e+6_fp
            chem(i, k, j, p_so4) = State_Chm%Species(II, JJ, k, gi_so4) * 1.0e+6_fp
            chem(i, k, j, p_so4s) = State_Chm%Species(II, JJ, k, gi_so4s) * 1.0e+6_fp
            chem(i, k, j, p_msa) = State_Chm%Species(II, JJ, k, gi_msa) * 1.0e+6_fp
            chem(i, k, j, p_nh3) = State_Chm%Species(II, JJ, k, gi_nh3) * 1.0e+6_fp
            chem(i, k, j, p_nh4) = State_Chm%Species(II, JJ, k, gi_nh4) * 1.0e+6_fp
            chem(i, k, j, p_nit) = State_Chm%Species(II, JJ, k, gi_nit) * 1.0e+6_fp
            chem(i, k, j, p_nits) = State_Chm%Species(II, JJ, k, gi_nits) * 1.0e+6_fp
            chem(i, k, j, p_bcpi) = State_Chm%Species(II, JJ, k, gi_bcpi) * 1.0e+6_fp
            chem(i, k, j, p_ocpi) = State_Chm%Species(II, JJ, k, gi_ocpi) * 1.0e+6_fp
            chem(i, k, j, p_bcpo) = State_Chm%Species(II, JJ, k, gi_bcpo) * 1.0e+6_fp
            chem(i, k, j, p_ocpo) = State_Chm%Species(II, JJ, k, gi_ocpo) * 1.0e+6_fp
            chem(i, k, j, p_dst1) = State_Chm%Species(II, JJ, k, gi_dst1) * 1.0e+6_fp
            chem(i, k, j, p_dst2) = State_Chm%Species(II, JJ, k, gi_dst2) * 1.0e+6_fp
            chem(i, k, j, p_dst3) = State_Chm%Species(II, JJ, k, gi_dst3) * 1.0e+6_fp
            chem(i, k, j, p_dst4) = State_Chm%Species(II, JJ, k, gi_dst4) * 1.0e+6_fp
            chem(i, k, j, p_sala) = State_Chm%Species(II, JJ, k, gi_sala) * 1.0e+6_fp
            chem(i, k, j, p_salc) = State_Chm%Species(II, JJ, k, gi_salc) * 1.0e+6_fp
            chem(i, k, j, p_br2) = State_Chm%Species(II, JJ, k, gi_br2) * 1.0e+6_fp
            chem(i, k, j, p_br) = State_Chm%Species(II, JJ, k, gi_br) * 1.0e+6_fp
            chem(i, k, j, p_bro) = State_Chm%Species(II, JJ, k, gi_bro) * 1.0e+6_fp
            chem(i, k, j, p_hobr) = State_Chm%Species(II, JJ, k, gi_hobr) * 1.0e+6_fp
            chem(i, k, j, p_hbr) = State_Chm%Species(II, JJ, k, gi_hbr) * 1.0e+6_fp
            chem(i, k, j, p_brno2) = State_Chm%Species(II, JJ, k, gi_brno2) * 1.0e+6_fp
            chem(i, k, j, p_brno3) = State_Chm%Species(II, JJ, k, gi_brno3) * 1.0e+6_fp
            chem(i, k, j, p_chbr3) = State_Chm%Species(II, JJ, k, gi_chbr3) * 1.0e+6_fp
            chem(i, k, j, p_ch2br2) = State_Chm%Species(II, JJ, k, gi_ch2br2) * 1.0e+6_fp
            chem(i, k, j, p_ch3br) = State_Chm%Species(II, JJ, k, gi_ch3br) * 1.0e+6_fp
            chem(i, k, j, p_mpn) = State_Chm%Species(II, JJ, k, gi_mpn) * 1.0e+6_fp
            chem(i, k, j, p_isopnd) = State_Chm%Species(II, JJ, k, gi_isopnd) * 1.0e+6_fp
            chem(i, k, j, p_isopnb) = State_Chm%Species(II, JJ, k, gi_isopnb) * 1.0e+6_fp
            chem(i, k, j, p_moba) = State_Chm%Species(II, JJ, k, gi_moba) * 1.0e+6_fp
            chem(i, k, j, p_propnn) = State_Chm%Species(II, JJ, k, gi_propnn) * 1.0e+6_fp
            chem(i, k, j, p_hac) = State_Chm%Species(II, JJ, k, gi_hac) * 1.0e+6_fp
            chem(i, k, j, p_glyc) = State_Chm%Species(II, JJ, k, gi_glyc) * 1.0e+6_fp
            chem(i, k, j, p_mvkn) = State_Chm%Species(II, JJ, k, gi_mvkn) * 1.0e+6_fp
            chem(i, k, j, p_macrn) = State_Chm%Species(II, JJ, k, gi_macrn) * 1.0e+6_fp
            chem(i, k, j, p_map) = State_Chm%Species(II, JJ, k, gi_map) * 1.0e+6_fp
            chem(i, k, j, p_no2) = State_Chm%Species(II, JJ, k, gi_no2) * 1.0e+6_fp
            chem(i, k, j, p_no3) = State_Chm%Species(II, JJ, k, gi_no3) * 1.0e+6_fp
            chem(i, k, j, p_hno2) = State_Chm%Species(II, JJ, k, gi_hno2) * 1.0e+6_fp
            chem(i, k, j, p_benzene) = State_Chm%Species(II, JJ, k, gi_benz) * 1.0e+6_fp
            chem(i, k, j, p_toluene) = State_Chm%Species(II, JJ, k, gi_tolu) * 1.0e+6_fp
            chem(i, k, j, p_xylenes) = State_Chm%Species(II, JJ, k, gi_xyle) * 1.0e+6_fp
            chem(i, k, j, p_mtpa) = State_Chm%Species(II, JJ, k, gi_mtpa) * 1.0e+6_fp
            chem(i, k, j, p_limon) = State_Chm%Species(II, JJ, k, gi_limo) * 1.0e+6_fp
            chem(i, k, j, p_mtpo) = State_Chm%Species(II, JJ, k, gi_mtpo) * 1.0e+6_fp
            chem(i, k, j, p_soap) = State_Chm%Species(II, JJ, k, gi_soap) * 1.0e+6_fp
            chem(i, k, j, p_soas) = State_Chm%Species(II, JJ, k, gi_soas) * 1.0e+6_fp
            chem(i, k, j, p_eoh) = State_Chm%Species(II, JJ, k, gi_eoh) * 1.0e+6_fp
            chem(i, k, j, p_mgly) = State_Chm%Species(II, JJ, k, gi_mgly) * 1.0e+6_fp
            chem(i, k, j, p_glyx) = State_Chm%Species(II, JJ, k, gi_glyx) * 1.0e+6_fp
            chem(i, k, j, p_acta) = State_Chm%Species(II, JJ, k, gi_acta) * 1.0e+6_fp
            chem(i, k, j, p_hpald) = State_Chm%Species(II, JJ, k, gi_hpald) * 1.0e+6_fp
            chem(i, k, j, p_dhdn) = State_Chm%Species(II, JJ, k, gi_dhdn) * 1.0e+6_fp
            chem(i, k, j, p_ethln) = State_Chm%Species(II, JJ, k, gi_ethln) * 1.0e+6_fp
            chem(i, k, j, p_hcooh) = State_Chm%Species(II, JJ, k, gi_hcooh) * 1.0e+6_fp
            chem(i, k, j, p_iepoxa) = State_Chm%Species(II, JJ, k, gi_iepoxa) * 1.0e+6_fp
            chem(i, k, j, p_iepoxb) = State_Chm%Species(II, JJ, k, gi_iepoxb) * 1.0e+6_fp
            chem(i, k, j, p_iepoxd) = State_Chm%Species(II, JJ, k, gi_iepoxd) * 1.0e+6_fp
            chem(i, k, j, p_isn1) = State_Chm%Species(II, JJ, k, gi_isn1) * 1.0e+6_fp
            chem(i, k, j, p_ripa) = State_Chm%Species(II, JJ, k, gi_ripa) * 1.0e+6_fp
            chem(i, k, j, p_ripb) = State_Chm%Species(II, JJ, k, gi_ripb) * 1.0e+6_fp
            chem(i, k, j, p_ripd) = State_Chm%Species(II, JJ, k, gi_ripd) * 1.0e+6_fp
            chem(i, k, j, p_imae) = State_Chm%Species(II, JJ, k, gi_imae) * 1.0e+6_fp
            chem(i, k, j, p_soaie) = State_Chm%Species(II, JJ, k, gi_soaie) * 1.0e+6_fp
            chem(i, k, j, p_soame) = State_Chm%Species(II, JJ, k, gi_soame) * 1.0e+6_fp
            chem(i, k, j, p_soagx) = State_Chm%Species(II, JJ, k, gi_soagx) * 1.0e+6_fp
            chem(i, k, j, p_soamg) = State_Chm%Species(II, JJ, k, gi_soamg) * 1.0e+6_fp
            chem(i, k, j, p_lvoc) = State_Chm%Species(II, JJ, k, gi_lvoc) * 1.0e+6_fp
            chem(i, k, j, p_lvocoa) = State_Chm%Species(II, JJ, k, gi_lvocoa) * 1.0e+6_fp
            chem(i, k, j, p_isn1og) = State_Chm%Species(II, JJ, k, gi_isn1og) * 1.0e+6_fp
            chem(i, k, j, p_isn1oa) = State_Chm%Species(II, JJ, k, gi_isn1oa) * 1.0e+6_fp
            chem(i, k, j, p_monits) = State_Chm%Species(II, JJ, k, gi_monits) * 1.0e+6_fp
            chem(i, k, j, p_monitu) = State_Chm%Species(II, JJ, k, gi_monitu) * 1.0e+6_fp
            chem(i, k, j, p_honit) = State_Chm%Species(II, JJ, k, gi_honit) * 1.0e+6_fp
            chem(i, k, j, p_ionita) = State_Chm%Species(II, JJ, k, gi_ionita) * 1.0e+6_fp
            chem(i, k, j, p_monita) = State_Chm%Species(II, JJ, k, gi_monita) * 1.0e+6_fp
            chem(i, k, j, p_indiol) = State_Chm%Species(II, JJ, k, gi_indiol) * 1.0e+6_fp
            chem(i, k, j, p_ipmn) = State_Chm%Species(II, JJ, k, gi_ipmn) * 1.0e+6_fp
            chem(i, k, j, p_hc187) = State_Chm%Species(II, JJ, k, gi_hc187) * 1.0e+6_fp
            chem(i, k, j, p_n2o) = State_Chm%Species(II, JJ, k, gi_n2o) * 1.0e+6_fp
            chem(i, k, j, p_ocs) = State_Chm%Species(II, JJ, k, gi_ocs) * 1.0e+6_fp
            chem(i, k, j, p_ch4) = State_Chm%Species(II, JJ, k, gi_ch4) * 1.0e+6_fp
            chem(i, k, j, p_brcl) = State_Chm%Species(II, JJ, k, gi_brcl) * 1.0e+6_fp
            chem(i, k, j, p_hcl) = State_Chm%Species(II, JJ, k, gi_hcl) * 1.0e+6_fp
            chem(i, k, j, p_ccl4) = State_Chm%Species(II, JJ, k, gi_ccl4) * 1.0e+6_fp
            chem(i, k, j, p_ch3cl) = State_Chm%Species(II, JJ, k, gi_ch3cl) * 1.0e+6_fp
            chem(i, k, j, p_ch3ccl3) = State_Chm%Species(II, JJ, k, gi_ch3ccl3) * 1.0e+6_fp
            chem(i, k, j, p_cfc113) = State_Chm%Species(II, JJ, k, gi_cfc113) * 1.0e+6_fp
            chem(i, k, j, p_cfc114) = State_Chm%Species(II, JJ, k, gi_cfc114) * 1.0e+6_fp
            chem(i, k, j, p_cfc115) = State_Chm%Species(II, JJ, k, gi_cfc115) * 1.0e+6_fp
            chem(i, k, j, p_hcfc123) = State_Chm%Species(II, JJ, k, gi_hcfc123) * 1.0e+6_fp
            chem(i, k, j, p_hcfc141b) = State_Chm%Species(II, JJ, k, gi_hcfc141b) * 1.0e+6_fp
            chem(i, k, j, p_hcfc142b) = State_Chm%Species(II, JJ, k, gi_hcfc142b) * 1.0e+6_fp
            chem(i, k, j, p_cfc11) = State_Chm%Species(II, JJ, k, gi_cfc11) * 1.0e+6_fp
            chem(i, k, j, p_cfc12) = State_Chm%Species(II, JJ, k, gi_cfc12) * 1.0e+6_fp
            chem(i, k, j, p_hcfc22) = State_Chm%Species(II, JJ, k, gi_hcfc22) * 1.0e+6_fp
            chem(i, k, j, p_h1211) = State_Chm%Species(II, JJ, k, gi_h1211) * 1.0e+6_fp
            chem(i, k, j, p_h1301) = State_Chm%Species(II, JJ, k, gi_h1301) * 1.0e+6_fp
            chem(i, k, j, p_h2402) = State_Chm%Species(II, JJ, k, gi_h2402) * 1.0e+6_fp
            chem(i, k, j, p_cl) = State_Chm%Species(II, JJ, k, gi_cl) * 1.0e+6_fp
            chem(i, k, j, p_clo) = State_Chm%Species(II, JJ, k, gi_clo) * 1.0e+6_fp
            chem(i, k, j, p_hocl) = State_Chm%Species(II, JJ, k, gi_hocl) * 1.0e+6_fp
            chem(i, k, j, p_clno3) = State_Chm%Species(II, JJ, k, gi_clno3) * 1.0e+6_fp
            chem(i, k, j, p_clno2) = State_Chm%Species(II, JJ, k, gi_clno2) * 1.0e+6_fp
            chem(i, k, j, p_cloo) = State_Chm%Species(II, JJ, k, gi_cloo) * 1.0e+6_fp
            chem(i, k, j, p_oclo) = State_Chm%Species(II, JJ, k, gi_oclo) * 1.0e+6_fp
            chem(i, k, j, p_cl2) = State_Chm%Species(II, JJ, k, gi_cl2) * 1.0e+6_fp
            chem(i, k, j, p_cl2o2) = State_Chm%Species(II, JJ, k, gi_cl2o2) * 1.0e+6_fp
            chem(i, k, j, p_h2o) = State_Chm%Species(II, JJ, k, gi_h2o) * 1.0e+6_fp
            chem(i, k, j, p_brsala) = State_Chm%Species(II, JJ, k, gi_brsala) * 1.0e+6_fp
            chem(i, k, j, p_brsalc) = State_Chm%Species(II, JJ, k, gi_brsalc) * 1.0e+6_fp
            chem(i, k, j, p_chcl3) = State_Chm%Species(II, JJ, k, gi_chcl3) * 1.0e+6_fp
            chem(i, k, j, p_ch2cl2) = State_Chm%Species(II, JJ, k, gi_ch2cl2) * 1.0e+6_fp
            chem(i, k, j, p_ch3i) = State_Chm%Species(II, JJ, k, gi_ch3i) * 1.0e+6_fp
            chem(i, k, j, p_ch2i2) = State_Chm%Species(II, JJ, k, gi_ch2i2) * 1.0e+6_fp
            chem(i, k, j, p_ch2icl) = State_Chm%Species(II, JJ, k, gi_ch2icl) * 1.0e+6_fp
            chem(i, k, j, p_ch2ibr) = State_Chm%Species(II, JJ, k, gi_ch2ibr) * 1.0e+6_fp
            chem(i, k, j, p_hoi) = State_Chm%Species(II, JJ, k, gi_hoi) * 1.0e+6_fp
            chem(i, k, j, p_i2) = State_Chm%Species(II, JJ, k, gi_i2) * 1.0e+6_fp
            chem(i, k, j, p_ibr) = State_Chm%Species(II, JJ, k, gi_ibr) * 1.0e+6_fp
            chem(i, k, j, p_icl) = State_Chm%Species(II, JJ, k, gi_icl) * 1.0e+6_fp
            chem(i, k, j, p_i) = State_Chm%Species(II, JJ, k, gi_i) * 1.0e+6_fp
            chem(i, k, j, p_io) = State_Chm%Species(II, JJ, k, gi_io) * 1.0e+6_fp
            chem(i, k, j, p_hi) = State_Chm%Species(II, JJ, k, gi_hi) * 1.0e+6_fp
            chem(i, k, j, p_oio) = State_Chm%Species(II, JJ, k, gi_oio) * 1.0e+6_fp
            chem(i, k, j, p_ino) = State_Chm%Species(II, JJ, k, gi_ino) * 1.0e+6_fp
            chem(i, k, j, p_iono) = State_Chm%Species(II, JJ, k, gi_iono) * 1.0e+6_fp
            chem(i, k, j, p_iono2) = State_Chm%Species(II, JJ, k, gi_iono2) * 1.0e+6_fp
            chem(i, k, j, p_i2o2) = State_Chm%Species(II, JJ, k, gi_i2o2) * 1.0e+6_fp
            chem(i, k, j, p_i2o3) = State_Chm%Species(II, JJ, k, gi_i2o3) * 1.0e+6_fp
            chem(i, k, j, p_i2o4) = State_Chm%Species(II, JJ, k, gi_i2o4) * 1.0e+6_fp
            chem(i, k, j, p_isala) = State_Chm%Species(II, JJ, k, gi_isala) * 1.0e+6_fp
            chem(i, k, j, p_isalc) = State_Chm%Species(II, JJ, k, gi_isalc) * 1.0e+6_fp
            chem(i, k, j, p_aeri) = State_Chm%Species(II, JJ, k, gi_aeri) * 1.0e+6_fp
            chem(i, k, j, p_pfe) = State_Chm%Species(II, JJ, k, gi_pfe) * 1.0e+6_fp
            chem(i, k, j, p_co2) = State_Chm%Species(II, JJ, k, gi_co2) * 1.0e+6_fp
            chem(i, k, j, p_lbro2h) = State_Chm%Species(II, JJ, k, gi_lbro2h) * 1.0e+6_fp
            chem(i, k, j, p_lbro2n) = State_Chm%Species(II, JJ, k, gi_lbro2n) * 1.0e+6_fp
            chem(i, k, j, p_lisopoh) = State_Chm%Species(II, JJ, k, gi_lisopoh) * 1.0e+6_fp
            chem(i, k, j, p_lisopno3) = State_Chm%Species(II, JJ, k, gi_lisopno3) * 1.0e+6_fp
            chem(i, k, j, p_ltro2h) = State_Chm%Species(II, JJ, k, gi_ltro2h) * 1.0e+6_fp
            chem(i, k, j, p_ltro2n) = State_Chm%Species(II, JJ, k, gi_ltro2n) * 1.0e+6_fp
            chem(i, k, j, p_lxro2h) = State_Chm%Species(II, JJ, k, gi_lxro2h) * 1.0e+6_fp
            chem(i, k, j, p_lxro2n) = State_Chm%Species(II, JJ, k, gi_lxro2n) * 1.0e+6_fp
            chem(i, k, j, p_pyac) = State_Chm%Species(II, JJ, k, gi_pyac) * 1.0e+6_fp
            chem(i, k, j, p_so4h1) = State_Chm%Species(II, JJ, k, gi_so4h1) * 1.0e+6_fp
            chem(i, k, j, p_so4h2) = State_Chm%Species(II, JJ, k, gi_so4h2) * 1.0e+6_fp
            chem(i, k, j, p_pox) = State_Chm%Species(II, JJ, k, gi_pox) * 1.0e+6_fp
            chem(i, k, j, p_lox) = State_Chm%Species(II, JJ, k, gi_lox) * 1.0e+6_fp
            chem(i, k, j, p_pco) = State_Chm%Species(II, JJ, k, gi_pco) * 1.0e+6_fp
            chem(i, k, j, p_lco) = State_Chm%Species(II, JJ, k, gi_lco) * 1.0e+6_fp
            chem(i, k, j, p_pso4) = State_Chm%Species(II, JJ, k, gi_pso4) * 1.0e+6_fp
            chem(i, k, j, p_lch4) = State_Chm%Species(II, JJ, k, gi_lch4) * 1.0e+6_fp
            chem(i, k, j, p_ph2o2) = State_Chm%Species(II, JJ, k, gi_ph2o2) * 1.0e+6_fp
            chem(i, k, j, p_dhdc) = State_Chm%Species(II, JJ, k, gi_dhdc) * 1.0e+6_fp
            chem(i, k, j, p_pmnn) = State_Chm%Species(II, JJ, k, gi_pmnn) * 1.0e+6_fp
            chem(i, k, j, p_tro2) = State_Chm%Species(II, JJ, k, gi_tro2) * 1.0e+6_fp
            chem(i, k, j, p_bro2) = State_Chm%Species(II, JJ, k, gi_bro2) * 1.0e+6_fp
            chem(i, k, j, p_n) = State_Chm%Species(II, JJ, k, gi_n) * 1.0e+6_fp
            chem(i, k, j, p_xro2) = State_Chm%Species(II, JJ, k, gi_xro2) * 1.0e+6_fp
            chem(i, k, j, p_imao3) = State_Chm%Species(II, JJ, k, gi_imao3) * 1.0e+6_fp
            chem(i, k, j, p_etp) = State_Chm%Species(II, JJ, k, gi_etp) * 1.0e+6_fp
            chem(i, k, j, p_ra3p) = State_Chm%Species(II, JJ, k, gi_ra3p) * 1.0e+6_fp
            chem(i, k, j, p_rb3p) = State_Chm%Species(II, JJ, k, gi_rb3p) * 1.0e+6_fp
            chem(i, k, j, p_rp) = State_Chm%Species(II, JJ, k, gi_rp) * 1.0e+6_fp
            chem(i, k, j, p_pp) = State_Chm%Species(II, JJ, k, gi_pp) * 1.0e+6_fp
            chem(i, k, j, p_prpn) = State_Chm%Species(II, JJ, k, gi_prpn) * 1.0e+6_fp
            chem(i, k, j, p_pip) = State_Chm%Species(II, JJ, k, gi_pip) * 1.0e+6_fp
            chem(i, k, j, p_r4p) = State_Chm%Species(II, JJ, k, gi_r4p) * 1.0e+6_fp
            chem(i, k, j, p_dhpcarp) = State_Chm%Species(II, JJ, k, gi_dhpcarp) * 1.0e+6_fp
            chem(i, k, j, p_iap) = State_Chm%Species(II, JJ, k, gi_iap) * 1.0e+6_fp
            chem(i, k, j, p_hpc52o2) = State_Chm%Species(II, JJ, k, gi_hpc52o2) * 1.0e+6_fp
            chem(i, k, j, p_vrp) = State_Chm%Species(II, JJ, k, gi_vrp) * 1.0e+6_fp
            chem(i, k, j, p_atooh) = State_Chm%Species(II, JJ, k, gi_atooh) * 1.0e+6_fp
            chem(i, k, j, p_dhmob) = State_Chm%Species(II, JJ, k, gi_dhmob) * 1.0e+6_fp
            chem(i, k, j, p_isnp) = State_Chm%Species(II, JJ, k, gi_isnp) * 1.0e+6_fp
            chem(i, k, j, p_maop) = State_Chm%Species(II, JJ, k, gi_maop) * 1.0e+6_fp
            chem(i, k, j, p_mrp) = State_Chm%Species(II, JJ, k, gi_mrp) * 1.0e+6_fp
            chem(i, k, j, p_inpn) = State_Chm%Species(II, JJ, k, gi_inpn) * 1.0e+6_fp
            chem(i, k, j, p_mobaoo) = State_Chm%Species(II, JJ, k, gi_mobaoo) * 1.0e+6_fp
            chem(i, k, j, p_diboo) = State_Chm%Species(II, JJ, k, gi_diboo) * 1.0e+6_fp
            chem(i, k, j, p_h) = State_Chm%Species(II, JJ, k, gi_h) * 1.0e+6_fp
            chem(i, k, j, p_macrno2) = State_Chm%Species(II, JJ, k, gi_macrno2) * 1.0e+6_fp
            chem(i, k, j, p_roh) = State_Chm%Species(II, JJ, k, gi_roh) * 1.0e+6_fp
            chem(i, k, j, p_isnohoo) = State_Chm%Species(II, JJ, k, gi_isnohoo) * 1.0e+6_fp
            chem(i, k, j, p_mvkoo) = State_Chm%Species(II, JJ, k, gi_mvkoo) * 1.0e+6_fp
            chem(i, k, j, p_isnoob) = State_Chm%Species(II, JJ, k, gi_isnoob) * 1.0e+6_fp
            chem(i, k, j, p_gaoo) = State_Chm%Species(II, JJ, k, gi_gaoo) * 1.0e+6_fp
            chem(i, k, j, p_ch3choo) = State_Chm%Species(II, JJ, k, gi_ch3choo) * 1.0e+6_fp
            chem(i, k, j, p_iepoxoo) = State_Chm%Species(II, JJ, k, gi_iepoxoo) * 1.0e+6_fp
            chem(i, k, j, p_mglyoo) = State_Chm%Species(II, JJ, k, gi_mglyoo) * 1.0e+6_fp
            chem(i, k, j, p_prn1) = State_Chm%Species(II, JJ, k, gi_prn1) * 1.0e+6_fp
            chem(i, k, j, p_mgloo) = State_Chm%Species(II, JJ, k, gi_mgloo) * 1.0e+6_fp
            chem(i, k, j, p_a3o2) = State_Chm%Species(II, JJ, k, gi_a3o2) * 1.0e+6_fp
            chem(i, k, j, p_man2) = State_Chm%Species(II, JJ, k, gi_man2) * 1.0e+6_fp
            chem(i, k, j, p_isnooa) = State_Chm%Species(II, JJ, k, gi_isnooa) * 1.0e+6_fp
            chem(i, k, j, p_po2) = State_Chm%Species(II, JJ, k, gi_po2) * 1.0e+6_fp
            chem(i, k, j, p_isopndo2) = State_Chm%Species(II, JJ, k, gi_isopndo2) * 1.0e+6_fp
            chem(i, k, j, p_b3o2) = State_Chm%Species(II, JJ, k, gi_b3o2) * 1.0e+6_fp
            chem(i, k, j, p_macroo) = State_Chm%Species(II, JJ, k, gi_macroo) * 1.0e+6_fp
            chem(i, k, j, p_r4n1) = State_Chm%Species(II, JJ, k, gi_r4n1) * 1.0e+6_fp
            chem(i, k, j, p_maopo2) = State_Chm%Species(II, JJ, k, gi_maopo2) * 1.0e+6_fp
            chem(i, k, j, p_ato2) = State_Chm%Species(II, JJ, k, gi_ato2) * 1.0e+6_fp
            chem(i, k, j, p_rco3) = State_Chm%Species(II, JJ, k, gi_rco3) * 1.0e+6_fp
            chem(i, k, j, p_limo2) = State_Chm%Species(II, JJ, k, gi_limo2) * 1.0e+6_fp
            chem(i, k, j, p_olnd) = State_Chm%Species(II, JJ, k, gi_olnd) * 1.0e+6_fp
            chem(i, k, j, p_olnn) = State_Chm%Species(II, JJ, k, gi_olnn) * 1.0e+6_fp
            chem(i, k, j, p_ko2) = State_Chm%Species(II, JJ, k, gi_ko2) * 1.0e+6_fp
            chem(i, k, j, p_isopnbo2) = State_Chm%Species(II, JJ, k, gi_isopnbo2) * 1.0e+6_fp
            chem(i, k, j, p_pio2) = State_Chm%Species(II, JJ, k, gi_pio2) * 1.0e+6_fp
            chem(i, k, j, p_hc5oo) = State_Chm%Species(II, JJ, k, gi_hc5oo) * 1.0e+6_fp
            chem(i, k, j, p_nmao3) = State_Chm%Species(II, JJ, k, gi_nmao3) * 1.0e+6_fp
            chem(i, k, j, p_vro2) = State_Chm%Species(II, JJ, k, gi_vro2) * 1.0e+6_fp
            chem(i, k, j, p_ch2oo) = State_Chm%Species(II, JJ, k, gi_ch2oo) * 1.0e+6_fp
            chem(i, k, j, p_hc5) = State_Chm%Species(II, JJ, k, gi_hc5) * 1.0e+6_fp
            chem(i, k, j, p_rio2) = State_Chm%Species(II, JJ, k, gi_rio2) * 1.0e+6_fp
            chem(i, k, j, p_eto2) = State_Chm%Species(II, JJ, k, gi_eto2) * 1.0e+6_fp
            chem(i, k, j, p_ino2) = State_Chm%Species(II, JJ, k, gi_ino2) * 1.0e+6_fp
            chem(i, k, j, p_r4o2) = State_Chm%Species(II, JJ, k, gi_r4o2) * 1.0e+6_fp
            chem(i, k, j, p_mro2) = State_Chm%Species(II, JJ, k, gi_mro2) * 1.0e+6_fp
            chem(i, k, j, p_mco3) = State_Chm%Species(II, JJ, k, gi_mco3) * 1.0e+6_fp
            chem(i, k, j, p_mo2) = State_Chm%Species(II, JJ, k, gi_mo2) * 1.0e+6_fp
            chem(i, k, j, p_o1d) = State_Chm%Species(II, JJ, k, gi_o1d) * 1.0e+6_fp
            chem(i, k, j, p_o) = State_Chm%Species(II, JJ, k, gi_o) * 1.0e+6_fp
            chem(i, k, j, p_ho2) = State_Chm%Species(II, JJ, k, gi_ho2) * 1.0e+6_fp
            chem(i, k, j, p_oh) = State_Chm%Species(II, JJ, k, gi_oh) * 1.0e+6_fp
            chem(i, k, j, p_h2) = State_Chm%Species(II, JJ, k, gi_h2) * 1.0e+6_fp
            chem(i, k, j, p_moh) = State_Chm%Species(II, JJ, k, gi_moh) * 1.0e+6_fp
            chem(i, k, j, p_n2) = State_Chm%Species(II, JJ, k, gi_n2) * 1.0e+6_fp
            chem(i, k, j, p_o2) = State_Chm%Species(II, JJ, k, gi_o2) * 1.0e+6_fp
            chem(i, k, j, p_rcooh) = State_Chm%Species(II, JJ, k, gi_rcooh) * 1.0e+6_fp

            ! add mercury, yxzhang, 06/20/2019
            chem(i, k, j, p_hg0) = State_Chm%Species(II, JJ, k, gi_hg0) * 1.0e+6_fp
            chem(i, k, j, p_hg2) = State_Chm%Species(II, JJ, k, gi_hg2) * 1.0e+6_fp
            chem(i, k, j, p_hgp) = State_Chm%Species(II, JJ, k, gi_hgp) * 1.0e+6_fp

            ! Diagnose Particulate Matter (PM2.5 Concentration in ug/m^3)
            ! Per GeosCore/aerosol_mod.F & http://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem
            !
            ! PM25 = 1.33 (NH4 + NIT  + SO4) + BCPI + BCPO +
            !        2.10 (OCPO + 1.16 OCPI) + 1.16 SOA    +
            !        DST1 + 0.38 DST2 + 1.86 SALA + SOAGX + SOAMG
            !
            if(allocated(PM25)) then
                grid%pm2_5_dry(i, k, j) = PM25(II, JJ, k) * 1.0e+9_fp ! kg/m^3 -> \mug/m^3
            endif

            !-----------------------------------------------------------------
            !         *****  N E T C D F   D I A G N O S T I C S  *****
            !    Needs to be done inside loop so we do not iterate over grid
            !         more than once. i, k, j in WRF --> II, JJ, k in G-C
            !-----------------------------------------------------------------

            !
            ! 3-dimensional values (G-C I, J, K    --> WRF I, K, J)
            !

            ! Dry deposition fluxes in chemistry, mixing
            ! Uses nDryDep map
            ! grid%drydepchm(i, k, j) = State_Diag%DryDepChm(II, JJ, nDryDepSpecies)
            ! grid%drydepmix(i, k, j) = State_Diag%DryDepMix(II, JJ, nDryDepSpecies)

            ! Dry deposition velocity
            ! grid%drydepvel(i, k, j) = State_Diag%DryDepVel(II, JJ, nDryDepSpecies)

            ! PRLMixFrac: Fraction of BL occupied by level (mixing)
            ! grid%pblmixfrc(i, k, j) = State_Diag%PBLMixFrac(II, JJ, k)

            ! PrecipFracLS: Frac of box in LS precip (wetdep)
            ! grid%precipfracls(i, k, j) = State_Diag%PrecipFracLS(II, JJ, k)

            !
            ! 4-dimensional values (G-C I, J, K, N --> WRF I, K, J)
            ! OK so this is a bit tricky. NetCDF diagnostics for GEOS-Chem, when in 4-D values
            ! cannot be passed as arrays to WRF directly, so you can choose UP to 4 species
            ! by gc_diagn_spc_n0...3 which will correspond to GEOS-Chem species IDs
            ! and these species will be diagnosed by WRF-GC.
            !
            ! This is a STOPGAP solution until we figure out output by G-C correctly.
            ! hplin, 10/13/18
            !

            ! 2-dimensional values (G-C I, J, N --> WRF I, J)
            ! Are for budget diagnostics summed over the Full column, PBL & Troposphere
            ! for respective G-C operator processes --
            !   Emissions + Dry Deposition
            !   Mixing
            !   Convection
            !   Chemistry
            !   Wet Deposition
            if(config_flags%gc_diagn_spc_n0 .ne. 0 .and. config_flags%gc_diagn_spc_n0 .le. State_Chm%nSpecies) then
                !grid%pblflux_n0(i, k, j) = State_Diag%PBLFlux(II, JJ, k, config_flags%gc_diagn_spc_n0)
                grid%cldcnvflx_n0(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n0)
                !grid%rainfrcls_n0(i, k, j) = State_Diag%RainFracLS(II, JJ, k, config_flags%gc_diagn_spc_n0)
                !grid%washfrcls_n0(i, k, j) = State_Diag%WashFracLS(II, JJ, k, config_flags%gc_diagn_spc_n0)

                ! Budget Diagnostics

                ! EmisDryDep uses nAdvect map
                grid%gcemisdrydep_full_n0(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcemisdrydep_trop_n0(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n0(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                ! Mixing uses nAdvect map
                grid%gcmixing_full_n0(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcmixing_trop_n0(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcmixing_pbl_n0(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                ! Convection uses nAdvect map
                grid%gcconv_full_n0(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcconv_trop_n0(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcconv_pbl_n0(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                ! Chemistry uses nAdvect map
                grid%gcchem_full_n0(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcchem_trop_n0(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)
                grid%gcchem_pbl_n0(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%AdvectID)

                ! Wet Deposition uses nWetDep map
                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n0(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%gcwetdep_trop_n0(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%gcwetdep_pbl_n0(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%wetlscnvfrc_n0(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                    grid%wetlscnv_n0(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)

                    grid%wetlossls_n0(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n1 .ne. 0 .and. config_flags%gc_diagn_spc_n1 .le. State_Chm%nSpecies) then
                !grid%pblflux_n1(i, k, j) = State_Diag%PBLFlux(II, JJ, k, config_flags%gc_diagn_spc_n1)
                grid%cldcnvflx_n1(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n1)
                !grid%rainfrcls_n1(i, k, j) = State_Diag%RainFracLS(II, JJ, k, config_flags%gc_diagn_spc_n1)
                !grid%washfrcls_n1(i, k, j) = State_Diag%WashFracLS(II, JJ, k, config_flags%gc_diagn_spc_n1)

                ! Budget Diagnostics
                grid%gcemisdrydep_full_n1(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcemisdrydep_trop_n1(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n1(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcmixing_full_n1(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcmixing_trop_n1(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcmixing_pbl_n1(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcconv_full_n1(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcconv_trop_n1(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcconv_pbl_n1(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)

                grid%gcchem_full_n1(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcchem_trop_n1(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)
                grid%gcchem_pbl_n1(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%AdvectID)


                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n1(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%gcwetdep_trop_n1(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%gcwetdep_pbl_n1(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)

                    grid%wetlscnvfrc_n1(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                    grid%wetlscnv_n1(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)

                    grid%wetlossls_n1(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n2 .ne. 0 .and. config_flags%gc_diagn_spc_n2 .le. State_Chm%nSpecies) then
                !grid%pblflux_n2(i, k, j) = State_Diag%PBLFlux(II, JJ, k, config_flags%gc_diagn_spc_n2)
                grid%cldcnvflx_n2(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n2)
                !grid%rainfrcls_n2(i, k, j) = State_Diag%RainFracLS(II, JJ, k, config_flags%gc_diagn_spc_n2)
                !grid%washfrcls_n2(i, k, j) = State_Diag%WashFracLS(II, JJ, k, config_flags%gc_diagn_spc_n2)

                ! Budget Diagnostics
                grid%gcemisdrydep_full_n2(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcemisdrydep_trop_n2(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n2(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcmixing_full_n2(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcmixing_trop_n2(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcmixing_pbl_n2(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcconv_full_n2(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcconv_trop_n2(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcconv_pbl_n2(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                grid%gcchem_full_n2(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcchem_trop_n2(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)
                grid%gcchem_pbl_n2(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%AdvectID)

                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n2(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%gcwetdep_trop_n2(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%gcwetdep_pbl_n2(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)

                    grid%wetlscnvfrc_n2(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                    grid%wetlscnv_n2(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)

                    grid%wetlossls_n2(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%WetDepId)
                endif
            endif

            if(config_flags%gc_diagn_spc_n3 .ne. 0 .and. config_flags%gc_diagn_spc_n3 .le. State_Chm%nSpecies) then
                !grid%pblflux_n3(i, k, j) = State_Diag%PBLFlux(II, JJ, k, config_flags%gc_diagn_spc_n3)
                grid%cldcnvflx_n3(i, k, j) = State_Diag%CloudConvFlux(II, JJ, k, config_flags%gc_diagn_spc_n3)
                !grid%rainfrcls_n3(i, k, j) = State_Diag%RainFracLS(II, JJ, k, config_flags%gc_diagn_spc_n3)
                !grid%washfrcls_n3(i, k, j) = State_Diag%WashFracLS(II, JJ, k, config_flags%gc_diagn_spc_n3)

                ! Budget Diagnostics
                grid%gcemisdrydep_full_n3(i, j) = State_Diag%BudgetEmisDryDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcemisdrydep_trop_n3(i, j) = State_Diag%BudgetEmisDryDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcemisdrydep_pbl_n3(i, j)  = State_Diag%BudgetEmisDryDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcmixing_full_n3(i, j) = State_Diag%BudgetMixingFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcmixing_trop_n3(i, j) = State_Diag%BudgetMixingTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcmixing_pbl_n3(i, j)  = State_Diag%BudgetMixingPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcconv_full_n3(i, j) = State_Diag%BudgetConvectionFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcconv_trop_n3(i, j) = State_Diag%BudgetConvectionTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcconv_pbl_n3(i, j)  = State_Diag%BudgetConvectionPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                grid%gcchem_full_n3(i, j) = State_Diag%BudgetChemistryFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcchem_trop_n3(i, j) = State_Diag%BudgetChemistryTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)
                grid%gcchem_pbl_n3(i, j)  = State_Diag%BudgetChemistryPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%AdvectID)

                if(State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId .ne. MISSING_INT) then
                    grid%gcwetdep_full_n3(i, j) = State_Diag%BudgetWetDepFull(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%gcwetdep_trop_n3(i, j) = State_Diag%BudgetWetDepTrop(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%gcwetdep_pbl_n3(i, j)  = State_Diag%BudgetWetDepPBL(II, JJ, &
                                                  State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)

                    grid%wetlscnvfrc_n3(i, k, j) = State_Diag%WetLossConvFrac(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                    grid%wetlscnv_n3(i, k, j) = State_Diag%WetLossConv(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)

                    grid%wetlossls_n3(i, k, j) = State_Diag%WetLossLS(II, JJ, k, &
                                                   State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%WetDepId)
                endif
            endif
        enddo
        enddo
        enddo

        !-----------------------------------------------------------------
        !         *****  D E B U G G I N G   R O U T I N E S  *****
        !-----------------------------------------------------------------

        ! After conversion if debug level is appropriately high
        if(debug_level .ge. 100) then
            ! Debug chemistry state and species information on PET (1, 1, 1) lowest lev.
            WRITE(6, *) "%%%% REGISTERED STATE_CHM SPECIES DUMP %%%%"
            do N = 1, MIN(State_Chm%nSpecies, 15)
                WRITE(6, *) "N:", N, "Name:", State_Chm%SpcData(N)%Info%Name
                ! WRITE(6, *) "Full Name:", State_Chm%SpcData(N)%Info%FullName
                WRITE(6, *) "Value at PET (1,1,1):", State_Chm%Species(1, 1, 1, N)
            enddo
        endif

        ! For hplin debug only...
        if(debug_level .ge. 1 .and. debug_level .lt. 100) then
            WRITE(6, *) "%%%% State_Chm IMPORTANT SPECIES DUMP - DEBUG WRF|GCHP %%%%"
            WRITE(6, *) "At I = 1, J = 1, L = 1, units: ", State_Chm%Spc_Units
            WRITE(6, *) "O3:", State_Chm%Species(1, 1, 1, IND_('O3'))
            ! WRITE(6, *) "O2:", State_Chm%Species(1, 1, 1, IND_('O2'))
            WRITE(6, *) "CO:", State_Chm%Species(1, 1, 1, IND_('CO'))
            WRITE(6, *) "NH3:", State_Chm%Species(1, 1, 1, IND_('NH3'))
            WRITE(6, *) "NO2:", State_Chm%Species(1, 1, 1, IND_('NO2'))

            WRITE(6, *) "Sums at L = 1"
            WRITE(6, *) "O3:", SUM(State_Chm%Species(:, :, 1, IND_('O3')))
            WRITE(6, *) "CO:", SUM(State_Chm%Species(:, :, 1, IND_('CO')))
            WRITE(6, *) "NH3:", SUM(State_Chm%Species(:, :, 1, IND_('NH3')))

            WRITE(6, *) "%%%% State_Diag DUMP - DEBUG WRF|GC %%%%"
            WRITE(6, *) "At I = 1, J = 1 (L = 1) - diagnostic states"
            WRITE(6, *) "Species#", config_flags%gc_diagn_spc_n0, " ",        &
                                    config_flags%gc_diagn_spc_n1, " ",        &
                                    config_flags%gc_diagn_spc_n2, " ",        &
                                    config_flags%gc_diagn_spc_n3
            WRITE(6, *) "Species Names -"
            IF(config_flags%gc_diagn_spc_n0 .ne. 0 .and. config_flags%gc_diagn_spc_n0 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "0: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n0)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n1 .ne. 0 .and. config_flags%gc_diagn_spc_n1 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "1: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n1)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n2 .ne. 0 .and. config_flags%gc_diagn_spc_n2 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "2: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n2)%Info%Name
            ENDIF
            IF(config_flags%gc_diagn_spc_n3 .ne. 0 .and. config_flags%gc_diagn_spc_n3 .le. State_Chm%nSpecies) THEN
                WRITE(6, *) "3: ", State_Chm%SpcData(config_flags%gc_diagn_spc_n3)%Info%Name
            ENDIF
            WRITE(6, *) "Budget Diags:", State_Diag%Archive_Budget,           &
                                         State_Diag%Archive_BudgetEmisDryDep, &
                                         State_Diag%Archive_BudgetMixing,     &
                                         State_Diag%Archive_BudgetConvection, &
                                         State_Diag%Archive_BudgetChemistry,  &
                                         State_Diag%Archive_BudgetWetDep
            WRITE(6, *) "(Glo|EmisDryDep|Mix|Conv|Chem|WetDep)"
            WRITE(6, *) "Values - FULL for N0", config_flags%gc_diagn_spc_n0
            WRITE(6, *) State_Diag%BudgetEmisDryDepFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetMixingFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetConvectionFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetChemistryFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) State_Diag%BudgetWetDepFull(1, 1, config_flags%gc_diagn_spc_n0)
            WRITE(6, *) "-- end of wrf|gc temporary diagnostics debug --"
        endif

        ! if(IM .ge. 12 .and. JM .ge. 5) then
        !     WRITE(6, *) "At I = 12, J = 5, L = :"
        !     WRITE(6, *) "H2O2:", State_Chm%Species(12, 5, :, IND_('H2O2'))
        !     WRITE(6, *) "SO2:", State_Chm%Species(12, 5, :, IND_('SO2'))
        ! endif
    end subroutine GIGC_Set_WRF
!EOC
!------------------------------------------------------------------------------
!                             The WRF-GCHP Project                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GIGC_Chunk_Cleanup
! 
! !DESCRIPTION: Subroutine GIGC\_Chunk\_Cleanup deallocates stateless information
!  in the GEOS-Chem Column Code, freeing up information for use in the next
!  timestep. This is crucial to the correct operation of the nested domain system.
!\\
!\\
! !INTERFACE:
!
    subroutine GIGC_Chunk_Cleanup(config_flags, kts, kte, Input_Opt, State_Met, State_Chm)
        implicit none

!
! !INPUT PARAMETERS:
!
        type(grid_config_rec_type), intent(in) :: config_flags
        integer, intent(in) :: kts, kte
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(OptInput), intent(inout) :: Input_Opt
        type(MetState), intent(inout) :: State_Met
        type(ChmState), intent(inout) :: State_Chm
!
! !REMARKS:
!  See notes below in code why some DEALLOCATE-ions are necessary.
!
! !REVISION HISTORY:
!  01 May 2018 - H.P. Lin   - First crack.
!EOP
!------------------------------------------------------------------------------
!BOC
        ! Deallocate all pointers.
        !
        ! Why are we not using State_Met_Mod::Cleanup_State_Met, you may ask?
        ! In ESMF_ mode (GEOS-Chem H.P.) State_Met_Mod cleanup performs NULLIFY on the variables,
        ! as opposed to GEOS-Chem Classic DEALLOCATE. This is because in GCHP, the POINTERs
        ! actually point to TARGETs, thus NULLIFY is the correct operation. But the way
        ! MET fields are set into GIGC_Convert_State_Mod in WRF-GCHP are more GCC and we directly
        ! ALLOCATE to them, and thus the right operation is not the one specified in Cleanup_State_Met.
        !
        ! It is an unfortunate side effect of different and mismatching design decisions, but whatever.

        ! -- 2-D fields -- !

        ! -- 3-D fields -- !

    end subroutine GIGC_Chunk_Cleanup
!EOC

end module GIGC_Convert_State_Mod
