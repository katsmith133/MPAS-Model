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
!  ocn_init_horizontal_advection
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

module ocn_init_horizontal_advection

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

   public :: ocn_init_setup_horizontal_advection, &
             ocn_init_validate_horizontal_advection

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_horizontal_advection
!
!> \brief   Setup for cvmix WSwSBF unit test configuration
!> \author  Todd Ringler
!> \date    04/23/2015
!> \details
!>  This routine sets up the initial conditions for the cvmix WSwSBF
!unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_horizontal_advection(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: temperature, salinity

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool,statePool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool

      type (mpas_pool_type), pointer :: tracersPool

      integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve,nEdgesSolve, nEdges, nVerticesSolve
      integer, pointer :: index_temperature, index_salinity,index_tracer1

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth,refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), pointer :: bottomDepth, angleEdge
      real (kind=RKIND), dimension(:), pointer :: fCell, fEdge, fVertex,xCell, yCell
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness,restingThickness, normalVelocity
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracers,debugTracers

      integer :: iCell, iEdge, iVertex, k, kML

      real (kind=RKIND) :: xDist, yDist, r_pt

      character (len=StrKIND), pointer :: config_init_configuration


      real (kind=RKIND), pointer ::config_horizontal_advection_temperature, &
                                   config_horizontal_advection_salinity,&
                                   config_horizontal_advection_zonalC, &
                                   config_horizontal_advection_meridionalC,&
                                   config_horizontal_advection_x_cent, &
                                   config_horizontal_advection_y_cent, &
                                   config_horizontal_advection_tracer_radius, &
                                   config_horizontal_advection_gaussian_width
      ! assume no error
      iErr = 0

      ! get and test if this is the configuration specified
      call mpas_pool_get_config(domain % configs,'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('horizontal_advection')) return

      ! load the remaining configuration parameters
      call mpas_pool_get_config(domain % configs,'config_horizontal_advection_temperature', &
                                config_horizontal_advection_temperature)
      call mpas_pool_get_config(domain % configs,'config_horizontal_advection_salinity', config_horizontal_advection_salinity)
      call mpas_pool_get_config(domain % configs,'config_horizontal_advection_zonalC', config_horizontal_advection_zonalC)
      call mpas_pool_get_config(domain % configs,'config_horizontal_advection_meridionalC', config_horizontal_advection_meridionalC)
      call mpas_pool_get_config(domain % configs,'config_horizontal_advection_tracer_radius', config_horizontal_advection_tracer_radius)
      call mpas_pool_get_config(domain % configs, 'config_horizontal_advection_y_cent', config_horizontal_advection_y_cent)
      call mpas_pool_get_config(domain % configs, 'config_horizontal_advection_x_cent', config_horizontal_advection_x_cent)
      call mpas_pool_get_config(domain % configs, 'config_horizontal_advection_gaussian_width', &
                                config_horizontal_advection_gaussian_width)

      ! load data that required to initialize the ocean simulation
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh',meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh',verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state',statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics',diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing',forcingPool)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels',nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve',nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve',nEdgesSolve)
        call mpas_pool_get_dimension(meshPool, 'nVerticesSolve',nVerticesSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature',index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity',index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1',index_tracer1)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
        call mpas_pool_get_array(statePool, 'normalVelocity',normalVelocity,1)

        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness',restingThickness)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

        call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers',debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness',layerThickness, 1)

        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = 300.0_RKIND 
        end do

        refZMid(1) = -50.0_RKIND
        do k = 2, nVertLevels
          refZMid(k) = refZMid(k-1) - 100.0_RKIND
        enddo

        do iCell = 1, nCellsSolve
           if(associated(activeTracers) ) then
              activeTracers(index_temperature,:,iCell) = config_horizontal_advection_temperature
              activeTracers(index_salinity,:,iCell) = config_horizontal_advection_salinity
           endif ! if (associated(activeTracer))

        ! as a place holder, have some debug tracer in the top few
        ! layers and zero below

        if(associated(debugTracers)) then
          if(config_horizontal_advection_gaussian_width > 0) then
            !use a gaussian hill for the tracer
            xDist = xCell(iCell) - config_horizontal_advection_x_cent
            yDist = yCell(iCell) - config_horizontal_advection_y_cent
            debugTracers(:,:,iCell) = exp(-config_horizontal_advection_gaussian_width*(xDist**2.0 + yDist**2.0))
          else
            !Use the bar tracer distribution
            if(config_horizontal_advection_zonalC .ne. 0) then
              r_pt = abs(xCell(iCell) - config_horizontal_advection_x_cent)
            elseif (config_horizontal_advection_meridionalC .ne. 0) then 
              r_pt = abs(yCell(iCell) - config_horizontal_advection_y_cent)
            else
              print *, 'error: either meridC or zonalC must be zero.'
            endif
            debugTracers(:,:,iCell) = 0.0_RKIND

            if(r_pt < config_horizontal_advection_tracer_radius) then
              debugTracers(:,:,iCell) = 1.0_RKIND
            endif
          endif
        endif

        ! Set layerThickness
        do k = 1, nVertLevels
           layerThickness(k, iCell) = 100.0_RKIND 
           restingThickness(k, iCell) = layerThickness(k, iCell)
        end do
        ! Set Coriolis parameter
        fCell(iCell) = 0.0_RKIND 

        ! Set bottomDepth
        bottomDepth(iCell) = 300.0_RKIND 

        ! Set maxLevelCell
        maxLevelCell(iCell) = nVertLevels

     end do  ! do iCell

     do iEdge=1, nEdgesSolve
        normalVelocity(:, iEdge) = ( cos(angleEdge(iEdge)) * config_horizontal_advection_zonalC) + &
                                   ( sin(angleEdge(iEdge)) * config_horizontal_advection_meridionalC )
     enddo
     block_ptr => block_ptr % next
   end do

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_horizontal_advection!}}}

!***********************************************************************
!
!  routine ocn_init_validate_horizontal_advection
!
!> \brief   Validation for CVMix WSwSBF mixing unit test case
!> \author  Doug Jacobsen
!> \date    04/01/2015
!> \details
!>  This routine validates the configuration options for the CVMix
!WSwSBF mixing unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_horizontal_advection(configPool, packagePool,iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels,config_horizontal_advection_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration',config_init_configuration)

      if(config_init_configuration .ne. trim('horizontal_advection')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels',config_vert_levels)

      config_vert_levels = 3

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_horizontal_advection!}}}

!***********************************************************************

end module ocn_init_horizontal_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
