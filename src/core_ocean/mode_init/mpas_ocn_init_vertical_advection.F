! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_vertical_advection
!
!> \brief MPAS ocean initialize case -- CVMix Unit Test
!> WSwSBF means Wind Stress with Surface Buoyancy Forcing
!> \author Todd Ringler
!> \date   04/23/2015
!> \details
!>  This module contains the routines for initializing the
!>  the cvmix WSwSBF unit test configuration. This in a
!>  single column configuration
!
!-----------------------------------------------------------------------

module ocn_init_vertical_advection

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_init_cell_markers
   use ocn_init_vertical_grids

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

   public :: ocn_init_setup_vertical_advection, &
             ocn_init_validate_vertical_advection

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_vertical_advection
!
!> \brief   Setup for cvmix WSwSBF unit test configuration
!> \author  Todd Ringler
!> \date    04/23/2015
!> \details
!>  This routine sets up the initial conditions for the cvmix WSwSBF
!unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_vertical_advection(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool,statePool
      type (mpas_pool_type), pointer :: diagnosticsPool

      type (mpas_pool_type), pointer :: tracersPool

      integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve,nEdgesSolve, nEdge
      integer, pointer :: index_temperature, index_salinity,index_tracer1

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth,refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), pointer :: salinityRestore,bottomDepth, angleEdge
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness,restingThickness, normalVelocity, vertAleTransportTop
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers,debugTracers

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      integer :: iCell, iEdge, iVertex, k, kML

      real (kind=RKIND) :: gaussMean, BLdepth, r_pt

      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_vertical_advection_vertical_grid

      integer, pointer ::                 config_vertical_advection_vert_levels

      real (kind=RKIND), pointer ::config_vertical_advection_temperature, &
                                   config_vertical_advection_salinity,&
                                          config_vertical_advection_bottom_depth,&
                                          config_vertical_advection_layer_thickness,          &
                                          config_vertical_advection_vert_velocity
      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs,'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('vertical_advection')) return

      ! build the vertical grid
      ! intent(out) is interfaceLocations. An array ranging from 0 to 1
      call mpas_pool_get_config(domain % configs,'config_vertical_advection_vertical_grid', config_vertical_advection_vertical_grid)
      call mpas_pool_get_dimension(domain % blocklist % dimensions,'nVertLevelsP1', nVertLevelsP1)
      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid(config_vertical_advection_vertical_grid,interfaceLocations)

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs,'config_vertical_advection_temperature', &
                                config_vertical_advection_temperature)
      call mpas_pool_get_config(domain % configs,'config_vertical_advection_salinity', config_vertical_advection_salinity)
      call mpas_pool_get_config(domain % configs,'config_vertical_advection_bottom_depth', config_vertical_advection_bottom_depth)
      call mpas_pool_get_config(domain % configs,'config_vertical_advection_layer_thickness',config_vertical_advection_layer_thickness)
      call mpas_pool_get_config(domain % configs, 'config_vertical_advection_vert_velocity', config_vertical_advection_vert_velocity )

      ! load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh',meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh',verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state',statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics',diagnosticsPool)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels',nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve',nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdge', nEdge)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve',nEdgesSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature',index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity',index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1',index_tracer1)

        call mpas_pool_get_array(meshPool, 'refBottomDepth',refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights',vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
        call mpas_pool_get_array(statePool, 'normalVelocity',normalVelocity)
        
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness',restingThickness)

        call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers',debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness',layerThickness, 1)

        call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_vertical_advection_bottom_depth *interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * config_vertical_advection_bottom_depth *(interfaceLocations(k) + interfaceLocations(k+1))
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND
        debugTracers(:,:,:) = 0.0_RKIND
        vertAleTransportTop(:,:) = 0.0_RKIND

        do iCell = 1, nCellsSolve

          do k=1,nVertLevels+1
            vertAleTransportTop(k, iCell) = config_vertical_advection_vert_velocity 
          enddo

          if(associated(activeTracers) ) then

            do k=1,nVertLevels
              activeTracers(index_temperature, k, iCell) = config_vertical_advection_temperature
              activeTracers(index_salinity, k, iCell) = config_vertical_advection_salinity
            enddo

          endif ! if (associated(activeTracer))

          ! Set bottomDepth
          bottomDepth(iCell) = config_vertical_advection_bottom_depth

          ! Set layerThickness
          do k = 1, nVertLevels
            layerThickness(k, iCell) = config_vertical_advection_bottom_depth * (interfaceLocations(k+1) - interfaceLocations(k))
            restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          ! Set tracer based on gaussian distribution
          do k = 1, nVertLevels
            gaussMean = config_vertical_advection_bottom_depth -100.0_RKIND- config_vertical_advection_layer_thickness
            debugTracers(:,k,iCell) = exp(-((refZMid(k) + gaussMean) / (config_vertical_advection_layer_thickness))**2)
          enddo

          ! Set bottomDepth
          bottomDepth(iCell) = config_vertical_advection_bottom_depth

          ! Set maxLevelCell
          maxLevelCell(iCell) = nVertLevels

     end do  ! do iCell
     block_ptr => block_ptr % next
   end do

   deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_vertical_advection!}}}

!***********************************************************************
!
!  routine ocn_init_validate_vertical_advection
!
!> \brief   Validation for CVMix WSwSBF mixing unit test case
!> \author  Doug Jacobsen
!> \date    04/01/2015
!> \details
!>  This routine validates the configuration options for the CVMix
!WSwSBF mixing unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_vertical_advection(configPool, packagePool,iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels,config_vertical_advection_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration',config_init_configuration)

      if(config_init_configuration .ne. trim('vertical_advection')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels',config_vert_levels)
      call mpas_pool_get_config(configPool,'config_vertical_advection_vert_levels', config_vertical_advection_vert_levels)

      if(config_vert_levels <= 0 .and. config_vertical_advection_vert_levels > 0) then
         config_vert_levels = config_vertical_advection_vert_levels
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for CVMix WSwSBF unittest case. Not given a usable value for ' &
                          // 'vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_vertical_advection!}}}

!***********************************************************************

end module ocn_init_vertical_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
