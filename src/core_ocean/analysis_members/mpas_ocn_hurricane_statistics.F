! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_hurricane_statistics
!
!> \brief MPAS ocean analysis mode member: hurricane statistics 
!> \author Luke Van Roekel & Karthik Balaguru
!> \date   July 2019
!> \details
!>  MPAS ocean analysis mode member: hurricane statistics
!>
!  computes dynamic temperature and tropical cyclone heat potential 
!     
!-----------------------------------------------------------------------

module ocn_hurricane_statistics

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_init_hurricane_statistics, &
             ocn_compute_hurricane_statistics, &
             ocn_restart_hurricane_statistics, &
             ocn_finalize_hurricane_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_hurricane_statistics
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Luke Van Roekel
!> \date    August 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_hurricane_statistics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_init_hurricane_statistics!}}}

!***********************************************************************
!
!  routine ocn_compute_hurricane_statistics
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Luke Van Roekel & Karthik Balaguru
!> \date    July 2019 
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_hurricane_statistics(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      type (mpas_pool_type), pointer :: mixedLayerDepthsAMPool,hurricaneStatisticsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: hurricaneStatisticsAM
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve
      integer :: k, iCell, i, refIndex, refLevel(1)
      integer, pointer :: index_temperature
      integer, dimension(:), pointer :: maxLevelCell

      logical :: found_t26

      real (kind=RKIND), dimension(:), pointer :: tThreshMLD, t26Depth
      real (kind=RKIND), dimension(:), pointer :: dThreshMLD, dynamicTemperature
      real (kind=RKIND), dimension(:), pointer :: tropicalCycloneHeatPotential
      real (kind=RKIND), dimension(:), pointer :: surfaceFrictionVelocity
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:), pointer :: zTop, zMid, BruntVaisalaFreqTop, layerThickness
      integer, pointer :: interp_type
      integer :: interp_local
      real (kind=RKIND) :: mldTemp,dTempThres, dDenThres, dTempGrad, dDenGrad
      real (kind=RKIND) :: dz,temp_ref_lev, den_ref_lev, dV, dVm1, dVp1, localVals(6)
      real (kind=RKIND) :: tcPotentialSum,t26Temp, depthDiff, obukhovLength
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_config(ocnConfigs, 'config_AM_mixedLayerDepths_enable', config_AM_mixedLayerDepths_enable)
      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
      call mpas_get_timeInterval(timeStepESMF, dt=dtSim)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'hurricaneStatisticsAM', hurricaneStatisticsAMPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

      call mpas_pool_get_config(domain % configs, 'config_AM_mixedLayerDepths_Tthreshold', tThresholdFlag)
      call mpas_pool_get_config(domain % configs, 'config_AM_hurricaneStatistics_interp_method', interp_type)
      call mpas_pool_get_config(domain % configs, 'config_AM_hurricaneStatistics_crit_friction_velocity', &
                                                   config_AM_hurricaneStatistics_crit_friction_velocity)
      call mpas_pool_get_config(domain % configs. 'config_AM_hurricaneStatistics_Tdy_averaging_depth', &
                                                   config_AM_hurricaneStatistics_Tdy_averaging_depth)

      if (interp_type == 1) interp_local = 1
      if (interp_type == 2) interp_local = 2
      if (interp_type == 3) interp_local = 3

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mixedLayerDepthsAM', mixedLayerDepthsAMPool)
         call mpas_pool_get_subpool(block % structs, 'hurricaneStatisticsAM', hurricaneStatisticsAMPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(tracersPool, 'activeTracers', tracers, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', &
                        potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
         call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dThreshMLD', dThreshMLD)
         call mpas_pool_get_array(mixedLayerDepthsAMPool, 'tThreshMLD', tThreshMLD)
         call mpas_pool_get_array(hurricaneStatisticsAMPool, 'tropicalCycloneHeatPotential', &
                                                  tropicalCycloneHeatPotential)
         call mpas_pool_get_array(hurricaneStatisticsAMPool, 'dynamicTemperature', &
                                                  dynamicTemperature)
         call mpas_pool_get_array(hurricaneStatisticsAMPool, 't26depth', t26depth) 
        !find depth of the 26C isotherm first

         do iCell = 1,nCellsSolve

               !Initialize RefIndex for cases of very shallow columns
               refIndex = maxLevelCell(iCell)

               if(tracers(index_temperature,1,iCell) < 26.0_RKIND) THEN
                 tropicalCycloneHeatPotential(iCell) = -9999.0_RKIND
               else ! compute value as surface temperature is warm enough
               found_t26 = .false.

               do k=2, maxLevelCell(iCell)-1
                  if(tracers(index_temperature,k,iCell) < 26.0_RKIND) then
                     localvals(2:3) = zMid(k-1:k,iCell)
                     dV = tracers(index_temperature,k-1,iCell)
                     dVp1 = tracers(index_temperature,k,iCell)
                     call interp_bw_levels(localVals(2),localVals(3), dV, dVp1, 26.0_RKIND,   &
                                           interp_local, t26Temp)
                     t26Temp = max(t26Temp, zMid(k,iCell))
                     t26Depth(iCell) = abs(min(t26Temp,zMid(k-1,iCell)))
                     found_t26=.true.
                     refIndex=k
                     exit
                   endif
               enddo
               !if 26 isotherm not found set the potential to -9999
               !integrate from surface to 26C isotherm

               tcPotentialSum = 0.0_RKIND
               do k=1,refIndex-1
                  tcPotentialSum = tcPotentialSum + layerThickness(k,iCell)*tracers(index_temperature,k,iCell)
               enddo
               ! add the remaining bit to the 26C depth
               depthDiff = t26Depth(iCell) - abs(zTop(refIndex,iCell))
               tcPotentialSum = tcPotentialSum + depthDiff*tracers(index_temperature,refIndex-1,iCell)

               tropicalCycloneHeatPotential(iCell) = tcPotentialSum / t26Depth(iCell)
               enddo !iCell

          endif !end TC heat potential code 

          !next compute the Tdy where L = D_MLD + ((2*u*^3*dt)/(kappa*N2))**(1/3), where N2 is X m below the MLD
          !steps
          !1) average bruntvaisala from MLD (ref k-1) to MLD+50 or bottom whichever is smaller
          !2) compute L
          !3) integrate temp from surface to L -- DONE!

          do iCell = 1,nCellsSolve

            do k=1, maxLevelCell(iCell)-1
              if(abs(zTop(k,iCell)) > dThreshMLD(iCell)) then
                startIndex=k
                stopAverage = min(dThreshMLD(iCell) + 50.0_RKIND, abs(zMid(maxLevelCell(iCell),iCell)))
                exit
              endif
            enddo

! average starts with a bit of k-1
            BVF_sum = (abs(zTop(refIndex,iCell)) - dThreshMLD(iCell)) * BruntVaisalaFreqTop(refIndex-1,iCell)
            do k=refIndex,maxLevelCell(iCell)-1
              if(abs(zTop(k,iCell)) > stopAverage) then
                BVF_sum = BVF_sum + (stopAverage - abs(zTop(k-1,iCell)))*BruntVaisalaFreqTop(k,iCell)
                exit
              else
                BVF_sum = BVF_sum + layerThickness(k,iCell) * BruntVaisalaFreqTop(k,iCell)
              endif
            enddo

            BVF_sum = BVF_sum / (stopAverage - dThreshMLD)
            obukhovLength = dThreshMLD(iCell) + ((2.0_RKIND*surfaceFrictionVelocity(iCell)**3.0_RKIND*dt) / (0.4_RKIND*BVF_sum))**(1.0/3.0)

            tcPotentialSum = 0.0_RKIND
            do k=1,maxLevelCell(iCell)-1
              if(abs(zTop(k+1,iCell)) > obukhovLength) then
                 depthDiff = obukhovLength - abs(zTop(k,iCell))
                 tcPotentialSum = tcPotentialSum + depthDiff*tracers(index_temperature,k,iCell)
                 exit
              else
                 tcPotentialSum = tcPotentialSum + layerThickness(k,iCell)*tracers(index_temperature,k,iCell)
              endif
              dynamicTemperature(iCell) = tcPotentialSum / obukhovLength
            enddo
        block => block % next
      end do

   end subroutine ocn_compute_hurricane_statistics!}}}

!***********************************************************************
!
!  routine interp_bw_levels
!
!> \brief   Interpolates between model layers
!> \author  Luke Van Roekel
!> \date    September 2015
!> \details
!>  This routine conducts computations to compute various field values
!>     between model levels (in pressure or depth) or could interpolate
!>      between temperature/salinity/density values.  Interpolations are
!>      of the form
!>       y = coeffs(1)*x^3 + coeffs(2)*x^2 + coeffs(3)*x + coeffs(4)
!
!-----------------------------------------------------------------------

   subroutine interp_bw_levels(y0,y1,x0,x1,xT,interp_f,yT,xm1,ym1)!{{{

   integer,intent(in) :: interp_f  ! linear, quadratic, or spline
   real(kind=RKIND),intent(in)  :: y0,y1,x0,x1,xT
   real(kind=RKIND),intent(inout) :: yT
   real(kind=RKIND),optional,intent(in) :: xm1,ym1
                ! these values are to match the slope at a given point

!------------------------------------------------------------------------
!
!  Local variables for the interpolations
!
!------------------------------------------------------------------------

   real(kind=RKIND) :: coeffs(4)   ! stores the coefficients for the interp
   real(kind=RKIND) :: Minv(4,4)   ! holds values for computing quad and spline
   real(kind=RKIND) :: det
   real(kind=RKIND) :: rhs(4)
   integer :: k,k2

   coeffs(:) = 0.0_RKIND
   Minv(:,:) = 0.0_RKIND
   rhs(:)    = 0.0_RKIND


   select case (interp_f)

       case (1)

          coeffs(2) = (y1-y0)/(x1-x0)
          coeffs(1) = y0 - coeffs(2)*x0
       case (2)

          det = -(x1-x0)**2
          rhs(1) = y0
          rhs(2) = y1

          if(present(xm1) .and. present(ym1)) then
             rhs(3) = (y0-ym1)/(x0-xm1)
          else
             rhs(3) = 0.0_RKIND
          endif

          Minv(1,1) = -1.0_RKIND/det
          Minv(1,2) = 1.0_RKIND/det
          Minv(1,3) = -1.0_RKIND/(x1-x0)
          Minv(2,1) = 2.0_RKIND*x0/det
          Minv(2,2) = -2.0_RKIND*x0/det
          Minv(2,3) = (x1+x0)/(x1-x0)
          Minv(3,1) = -(x0**2)/det
          Minv(3,2) = x1*(2.0_RKIND*x0-x1)/det
          Minv(3,3) = -x1*x0/(x1-x0)

          do k=1,3
             do k2=1,3
                coeffs(k2) = coeffs(k2) + Minv(4-k2,k)*rhs(k)
             enddo
          enddo

        case (3)
          det = -(x1-x0)**3
          rhs(1) = y1
          rhs(2) = y0
          if(present(xm1) .and. present(ym1)) then
             rhs(3) = (y0-ym1)/(x0-xm1)
          else
             rhs(3) = 0.0_RKIND
          endif

          rhs(4) = (y1-y0)/(x1-x0)

          Minv(1,1) = 2.0_RKIND/det
          Minv(1,2) = -2.0_RKIND/det
          Minv(1,3) = (x0-x1)/det
          Minv(1,4) = (x0-x1)/det
          Minv(2,1) = -3.0_RKIND * (x1+x0)/det
          Minv(2,2) = 3.0_RKIND*(x1+x0)/det
          Minv(2,3) = (x1-x0)*(2.0_RKIND*x1+x0)/det
          Minv(2,4) = (x1-x0)*(2.0_RKIND*x0+x1)/det
          Minv(3,1) = 6.0_RKIND*x1*x0/det
          Minv(3,2) = -6.0_RKIND*x1*x0/det
          Minv(3,3) = -x1*(x1-x0)*(2.0_RKIND*x0+x1)/det
          Minv(3,4) = -x0*(x1-x0)*(2.0_RKIND*x1+x0)/det
          Minv(4,1) = -(x0**2)*(3.0_RKIND*x1-x0)/det
          Minv(4,2) = -(x1**2)*(-3.0_RKIND*x0+x1)/det
          Minv(4,3) = x0*(x1**2)*(x1-x0)/det
          Minv(4,4) = x1*(x0**2)*(x1-x0)/det

          do k=1,4
             do k2=1,4
                coeffs(k2) = coeffs(k2) + Minv(5-k2,k)*rhs(k)
             enddo
          enddo

     end select

     yT = coeffs(4)*xT**3 + coeffs(3)*xT**2 + coeffs(2)*xT + coeffs(1)
   end subroutine interp_bw_levels!}}}

!***********************************************************************
!
!  routine ocn_restart_hurricane_statistics
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Luke Van Roekel & Karthik Balaguru
!> \date    July 2019
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_hurricane_statistics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_restart_hurricane_statistics!}}}

!***********************************************************************
!
!  routine ocn_finalize_hurricane_statistics
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Luke Van Roekel & Karthik Balaguru
!> \date    July 2019 
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_hurricane_statistics(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_finalize_hurricane_statistics!}}}

end module ocn_hurricane_statistics

! vim: foldmethod=marker
