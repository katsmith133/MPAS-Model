!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_cvmix
!
!> \brief MPAS ocean vertical mixing interface to new CVMix that runs on GPU
!> \author Luke Van Roekel
!> \date   16 January 2019
!> \details
!>  This module contains the routines for calls into CVMix-Kokkos
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_cvmix2_kokkos

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log

   use ocn_constants

   use,intrinsic :: iso_c_binding

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_coefs_cvmix2_kokkos_build, &
             ocn_vmix_cvmix2_kokkos_init,        &
             ocn_vmix_cvmix2_kokkos_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type, bind(c) :: CVMix_mixingOpts
     real (C_DOUBLE) :: bl1, bl2, bl3, bl4, prandtl
     real (C_DOUBLE) :: PP_nu_zero, PP_alpha, PP_exp, PP_nu_b, PP_kappa_b, KPP_nu_zero, KPP_Ri_zero, KPP_exp
     real (c_double) :: convect_diff, convect_visc, BVsqr_convect
     logical (C_BOOL) :: bryanLewis, lBruntVaisala, lShearKPP, lShearPP, &
                         lEkman, lMonOb, lKPP, computelangmuir
     logical (C_BOOL) :: LANGMUIR_LWF16, LANGMUIR_lifk17, parabolicNonLocal
     real (c_double) :: Ri_crit, surf_layer_ext, c_s, minVtsqr, rho_o, minOSBLUnderIce
     real (c_double) :: a_s, c_m, a_m, non_local_coeff, numRiSmoothPasses
   end type CVMix_mixingOpts

   INTERFACE
     SUBROUTINE f_cvmix2_dummy(nCells,nVertLevels,layerThickness,mixingOpts) &
                     BIND(c, NAME='c_compute_cvmix2_mixing_dummy')

             USE mpas_kind_types
             USE, INTRINSIC :: ISO_C_BINDING
             import :: CVMix_mixingOpts
             IMPLICIT none
             integer ( c_int ), value :: nCells,nVertLevels
             real(kind=c_double) :: layerThickness(nVertLevels,nCells)
             type (CVMix_mixingOpts) :: mixingOpts
     END SUBROUTINE f_cvmix2_dummy
   END INTERFACE


   INTERFACE
     SUBROUTINE f_kokkos_initialize(nVertLevels,nCols) &
       BIND(c, NAME='c_kokkos_initialize')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT none
       integer(c_int),value :: nCols,nVertLevels
     END SUBROUTINE f_kokkos_initialize
   END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE
     SUBROUTINE f_kokkos_finalize() &
       BIND(c, NAME='c_kokkos_finalize')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
     END SUBROUTINE f_kokkos_finalize
   END INTERFACE

   interface
     subroutine f_cvmix2_compute_mixing(nCells,nVertLevels,visc,diff, layerThickness, U, V, &
       NLT, N2, RI_grad, rho, fCor, sfcBuoyancyForcing, sfcFrictionVelocity, maxLevel, &
       OSBL, efactor, lasl, wind10, mixingOpts)                                             &
       BIND(c, NAME='c_compute_cvmix2_mixing') 

        use mpas_kind_types
        use,intrinsic :: iso_c_binding
        import :: CVMix_mixingOpts
        implicit none
        integer ( c_int ), value :: nCells,nVertLevels
        real(kind=RKIND) :: visc(nVertLevels+1,nCells), diff(nVertLevels+1,nCells),      &   
           layerThickness(nVertLevels,nCells), NLT(nVertLevels+1,nCells),       &
           N2(nVertLevels+1,nCells),                 &
           Ri_grad(nVertLevels+1,nCells), sfcBuoyancyForcing(nCells),           &
           sfcFrictionVelocity(nCells), U(nVertLevels, nCells), V(nVertLevels,nCells),  &
           rho(nVertLevels,nCells), fCor(nCells),                  &
           OSBL(nCells), efactor(nCells),lasl(nCells),   &   
           wind10(nCells)
        integer (c_int) :: maxLevel(nCells)
        type (CVMix_mixingOpts) :: mixingOpts
     end subroutine f_cvmix2_compute_mixing 
   end interface

   logical :: cvmixOn, cvmixBackgroundOn, cvmixConvectionOn, cvmixKPPOn
   real (kind=RKIND) :: backgroundVisc, backgroundDiff

   type (CVMix_mixingOpts) :: mixingOpts

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_cmvix_build
!
!> \brief   Computes mixing coefficients using CVMix
!> \author  Todd Ringler
!> \date    04 February 2013
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers by calling CVMix routines.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_cvmix2_kokkos_build(meshPool, statePool, forcingPool, diagnosticsPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

! Things needed
! nCells, nVertLevs, viscosity array, diff array, layer thickness
! velocityZonal, velocityMeridional, nonlocal T, nonLocal S, BruntVaisala, gradient rich, p_rho, 
! fcor, sfcBuoy, sfcFric, maxLevelCell, OSBLdepth, langmuir efactor, lasl, 10m wind, mixingOpts
! efactor and lasl sent in doesn't matter what's inn them.

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool   !< Input/Output: diagnostic information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer :: &
        maxLevelCell, nEdgesOnCell

      real (kind=RKIND), dimension(:), pointer :: &
        surfaceBuoyancyForcing, surfaceFrictionVelocity, fCell, &
        boundaryLayerDepth, iceFraction, windSpeed10m, lasl

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        zMid, zTop, potentialDensity, bulkRichardsonNumber, RiTopOfCell, BruntVaisalaFreqTop,  &
        velocityMeridional, velocityZonal

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer, pointer :: index_vertNonLocalFluxTemp, config_cvmix_num_ri_smooth_loops

      logical, pointer :: config_use_cvmix_shear, config_use_cvmix_convection, config_use_cvmix_kpp
      logical, pointer :: config_use_cvmix_fixed_boundary_layer, config_cvmix_use_BLD_smoothing
      real (kind=RKIND), pointer :: config_cvmix_kpp_stop_OBL_search, config_cvmix_kpp_criticalBulkRichardsonNumber
      real (kind=RKIND), pointer :: config_cvmix_kpp_boundary_layer_depth, config_cvmix_kpp_surface_layer_extent
      real (kind=RKIND), pointer :: configure_cvmix_kpp_minimum_OBL_under_sea_ice
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme, config_cvmix_kpp_matching

      integer :: k, i, iCell, timeLevel, nCells
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray

      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. cvmixOn ) return

      !
      ! set parameters
      !
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_criticalBulkRichardsonNumber', &
                                config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_num_ri_smooth_loops', config_cvmix_num_ri_smooth_loops)
      call mpas_pool_get_config(ocnConfigs, 'configure_cvmix2_kpp_minimum_OBL_under_sea_ice', &
                                configure_cvmix_kpp_minimum_OBL_under_sea_ice)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional) 

      !
      ! set pointers for fields related to position on sphere
      !
      call mpas_pool_get_array(meshPool, 'fCell', fCell)

      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      !
      ! set pointers for fields related ocean state
      !
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)

      !
      ! set pointers for fields related to ocean forcing state
      !
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_array(forcingPool, 'windSpeed10m', windSpeed10m)
      call mpas_pool_get_array(diagnosticsPool, 'LaSL', lasl)

      !
      ! set pointers for fields related forcing at ocean surface
      !
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)

      !
      ! set pointers for viscosity/diffusivity and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)


      !
      ! set pointers for nonlocal flux and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_vertNonLocalFluxTemp', index_vertNonLocalFluxTemp)

      nCells = nCellsArray( size(nCellsArray) )

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertViscTopOfCell(:, iCell) = 0.0_RKIND
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
         vertNonLocalFlux(:, :, iCell) = 0.0_RKIND
         boundaryLayerDepth(iCell) = layerThickness(1,iCell)
      end do
      !$omp end do

      call mpas_timer_start('cvmix2.0 compute coefficients', .false.)
     
      call f_cvmix2_compute_mixing(nCells, nVertLevels, vertViscTopOfCell,      &
                  vertDiffTopOfCell, layerThickness, velocityZonal,       &
                  velocityMeridional, vertNonLocalFlux,            &
                  BruntVaisalaFreqTop,          &
                  RiTopOfCell, potentialDensity, fCell,                   &
                  surfaceBuoyancyForcing, surfaceFrictionVelocity, &
                  maxLevelCell, boundaryLayerDepth, iceFraction, lasl,    &
                  windSpeed10m, mixingOpts)

      call mpas_timer_stop('cvmix2.0 compute coefficients')

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_cvmix2_kokkos_build!}}}

   subroutine ocn_vmix_cvmix2_kokkos_finalize()!{{{

           call f_kokkos_finalize()

   end subroutine ocn_vmix_cvmix2_kokkos_finalize!}}}

!***********************************************************************
!
!  routine ocn_vmix_cvmix2_kokkos_init
!
!> \brief   Initializes ocean vertical mixing quantities for cvmix 2.0
!> \author  Luke Van Roekel 
!> \date    19 February 2019
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean for cvmix2
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_cvmix2_kokkos_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (mpas_pool_type), pointer :: meshPool

      integer, intent(out) :: err !< Output: error flag

      integer :: nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, pointer :: nVertLevels, config_cvmix_num_ri_smooth_loops
      type (block_type), pointer :: block

      ! CVMix
      logical, pointer :: config_use_cvmix

      ! background
      logical, pointer :: config_use_cvmix_background
      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion
      real (kind=RKIND), pointer :: config_cvmix_prandtl_number

      ! Shear configs
      logical, pointer :: config_use_cvmix_shear
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme,   &
                                          config_cvmix2_background_mixing_scheme
      real (kind=RKIND), pointer :: config_cvmix_shear_PP_nu_zero, config_cvmix_shear_PP_alpha, config_cvmix_shear_PP_exp, &
                                    config_cvmix_shear_KPP_nu_zero, config_cvmix_shear_KPP_Ri_zero, config_cvmix_shear_KPP_exp

      ! Convection configs
      logical, pointer ::  config_use_cvmix_convection
      real (kind=RKIND), pointer :: config_cvmix_convective_diffusion, config_cvmix_convective_viscosity, &
                                    config_cvmix_convective_triggerBVF
      logical, pointer :: config_cvmix_convective_basedOnBVF, config_cvmix_kpp_use_enhanced_diff

      ! KPP configs
      logical, pointer :: config_use_cvmix_kpp, config_cvmix2_computeLangmuir, &
                                    config_cvmix2_addLFK17_langmuir, config_cvmix2_useParabolicNonlocal
      logical, pointer :: config_cvmix_kpp_EkmanOBL, config_cvmix_kpp_MonObOBL
      real (kind=RKIND), pointer :: config_cvmix_kpp_criticalBulkRichardsonNumber,      &
                                    config_cvmix_kpp_surface_layer_extent,              &
                                    config_cvmix2_minVTsqr, config_density0,            &
                                    config_cvmix2_as_value, config_cvmix2_cs_value,     &
                                    config_cvmix2_cm_value, config_cvmix2_nonLocalTermStrength, &
                                    configure_cvmix_kpp_minimum_OBL_under_sea_ice,      &
                                    config_cvmix2_am_value, config_cvmix2_bl_val1,      &
                                    config_cvmix2_bl_val2,  config_cvmix2_bl_val3,      &
                                    config_cvmix2_bl_val4
      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      ! Options shared with original CVMix driver
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_prandtl_number', config_cvmix_prandtl_number)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_criticalBulkRichardsonNumber', &
                                config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_num_ri_smooth_loops', config_cvmix_num_ri_smooth_loops)
      call mpas_pool_get_config(ocnConfigs, 'configure_cvmix2_kpp_minimum_OBL_under_sea_ice', &
                                configure_cvmix_kpp_minimum_OBL_under_sea_ice)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2', config_use_cvmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix2_constant', config_use_cvmix_background)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_PP_nu_zero', config_cvmix_shear_PP_nu_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_PP_alpha', config_cvmix_shear_PP_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_PP_exp', config_cvmix_shear_PP_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_KPP_nu_zero', config_cvmix_shear_KPP_nu_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_KPP_Ri_zero', config_cvmix_shear_KPP_Ri_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_shear_KPP_exp', config_cvmix_shear_KPP_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_convective_basedOnBVF', config_cvmix_convective_basedOnBVF)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_convective_triggerBVF', config_cvmix_convective_triggerBVF)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_EkmanOBL', config_cvmix_kpp_EkmanOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_MonObOBL', config_cvmix_kpp_MonObOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_convective_diffusion', config_cvmix_convective_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_convective_viscosity', config_cvmix_convective_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_use_enhanced_diff', config_cvmix_kpp_use_enhanced_diff)

      ! New options for CVMix 2.0
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_background_mixing_scheme',   &
                                config_cvmix2_background_mixing_scheme)
      !computeLanmuir (Logical), addLFK17(logical), PNL (log), (reals follow) epsilon, cs, minVTsqr, rho_o,
      !(reals) as,cm,am,non_local_coef
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_surface_layer_extent', &
                                config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_kpp_minimumVt2', config_cvmix2_minVTsqr)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_computeLangmuir', config_cvmix2_computeLangmuir)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_addLFK17_langmuir', config_cvmix2_addLFK17_langmuir)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_useParabolicNonlocal', config_cvmix2_useParabolicNonlocal)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_surface_layer_extent',    &
                                config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_cs_value', config_cvmix2_cs_value)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_minVTsqr', config_cvmix2_minVTsqr)
      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_as_value', config_cvmix2_as_value)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_cm_value', config_cvmix2_cm_value)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_am_value', config_cvmix2_am_value)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_nonLocalTermStrength', config_cvmix2_nonLocalTermStrength)

      cvmixOn = config_use_cvmix
      if (.not.config_use_cvmix) return

      ! init backround mixing
      backgroundDiff = 0.0_RKIND
      backgroundVisc = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         nCells = nCellsArray( size(nCellsArray) )
         call f_kokkos_initialize(nVertLevels, nCells)
         block => block % next
      end do

      if(config_cvmix2_background_mixing_scheme == 'bryanLewis') then
         call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_bl_val1', config_cvmix2_bl_val1)
         call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_bl_val2', config_cvmix2_bl_val2)
         call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_bl_val3', config_cvmix2_bl_val3)
         call mpas_pool_get_config(ocnConfigs, 'config_cvmix2_bl_val4', config_cvmix2_bl_val4)
      ! set Bryan and Lewis background profile specific values
         mixingOpts%bl1 = config_cvmix2_bl_val1
         mixingOpts%bl2 = config_cvmix2_bl_val2
         mixingOpts%bl3 = config_cvmix2_bl_val3
         mixingOpts%bl4 = config_cvmix2_bl_val4
         mixingOpts%bryanLewis = .true.
         mixingOpts%prandtl = config_cvmix_prandtl_number
      elseif(config_cvmix2_background_mixing_scheme == 'constant') then
         backgroundVisc = config_cvmix_background_viscosity
         backgroundDiff = config_cvmix_background_diffusion
      else
         call mpas_log_write("unknown option for config_cvmix2_background_mixing_scheme, allowed values are"  &
                   // " 'bryanLewis' or 'constant'", MPAS_LOG_CRIT)
      endif
      cvmixConvectionOn = config_use_cvmix_convection
      cvmixKPPOn = config_use_cvmix_kpp

      
      !
      ! When CVMix is turned on, all other vertical mixing schemes should be off
      ! Test to make sure this is the case.
      !
      ! test here, err=1 if a problem

      !
      ! pull nVertLevels out of the mesh structure
      !
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      !
      ! initialize shear-based mixing
      !
      if (config_use_cvmix_shear) then
        if (.not. config_use_cvmix_background .and. trim(config_cvmix_shear_mixing_scheme) == 'PP') then
            call mpas_log_write("config_use_cvmix_shear cannot be used with with config_cvmix_shear_mixing_scheme = 'PP'" &
                 // "       without config_use_cvmix_background = .true.", MPAS_LOG_CRIT)
            err = 1
            return
        end if
            mixingOpts%lShearKPP = .false.
            mixingOpts%lShearPP = .false.
            if( trim(config_cvmix_shear_mixing_scheme) == 'KPP' ) then
                mixingOpts%lShearKPP = .true.
                mixingOpts%KPP_nu_zero = config_cvmix_shear_KPP_nu_zero
                mixingOpts%KPP_Ri_zero = config_cvmix_shear_KPP_Ri_zero
                mixingOpts%KPP_exp = config_cvmix_shear_KPP_exp
            elseif ( trim(config_cvmix_shear_mixing_scheme) == 'PP' ) then
                mixingOpts%lShearPP = .true.
                mixingopts%PP_alpha = config_cvmix_shear_PP_alpha
                mixingOpts%PP_exp = config_cvmix_shear_PP_exp
                mixingOpts%PP_nu_zero = config_cvmix_shear_PP_nu_zero
                mixingOpts%PP_nu_b = config_cvmix_background_viscosity 
                mixingOpts%PP_kappa_b = config_cvmix_background_diffusion
            else 
                call mpas_log_write("unknown option for config_cvmix2_shear_mixing_scheme, allowed values are"  &
                        // " 'KPP' or 'PP'", MPAS_LOG_CRIT)
            endif

            mixingOpts%numRiSmoothPasses = config_cvmix_num_ri_smooth_loops
      endif

      !
      ! initialize convective mixing
      !
      if (config_use_cvmix_convection) then

        ! config_cvmix_convective_basedOnBVF is not supported at this time
        if (.not.config_cvmix_convective_basedOnBVF) then
            call mpas_log_write("config_cvmix_convective_basedOnBVF = .false. is not supported. Change to true.", MPAS_LOG_CRIT)
            err = 1
            return
        endif

        mixingOpts%lBruntVaisala = config_use_cvmix_convection
        mixingOpts%BVsqr_convect = config_cvmix_convective_triggerBVF
        mixingOpts%convect_diff = config_cvmix_convective_diffusion
        mixingOpts%convect_visc = config_cvmix_convective_viscosity

      endif

      !
      ! FUTURE: initialize tidal mixing
      !

      !
      ! FUTURE: initialize double diffusion
      !

      !
      ! FUTURE: initialize KPP boundary layer scheme
      !


      if (config_use_cvmix_kpp) then
         mixingOpts%Ri_crit = config_cvmix_kpp_criticalBulkRichardsonNumber
         mixingOpts%lEkman = config_cvmix_kpp_EkmanOBL
         mixingOpts%lMonOb = config_cvmix_kpp_MonObOBL
         mixingOpts%lKPP = config_use_cvmix_kpp
         mixingOpts%surf_layer_ext = config_cvmix_kpp_surface_layer_extent
         mixingOpts%c_s = config_cvmix2_cs_value
         mixingOpts%minVtsqr = config_cvmix2_minVTsqr
         mixingOpts%rho_o = config_density0
         mixingOpts%minOSBLUnderIce = configure_cvmix_kpp_minimum_OBL_under_sea_ice
         mixingOpts%a_s = config_cvmix2_as_value
         mixingOpts%c_m = config_cvmix2_cm_value 
         mixingOpts%a_m = config_cvmix2_am_value
         mixingOpts%non_local_coeff = config_cvmix2_nonLocalTermStrength
         mixingOpts%parabolicNonLocal = config_cvmix2_useParabolicNonlocal

        !langmuir stuff
         mixingOpts%computeLangmuir = config_cvmix2_computeLangmuir
         if(mixingOpts%computeLangmuir) then
            mixingOpts%LANGMUIR_LWF16 = .true.     
            if( config_cvmix2_addLFK17_langmuir ) then
               mixingOpts%LANGMUIR_lifk17 = config_cvmix2_addLFK17_langmuir
            endif
         else
            mixingOpts%LANGMUIR_LWF16 = .false.
         endif

      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix2_kokkos_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix2_kokkos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
