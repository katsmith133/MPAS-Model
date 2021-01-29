!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_compute_nonhydrostatic_pressure
!
!> \brief MPAS ocean module to solve for the nonhydrostatic pressure correction
!> \author Luke Van Roekel
!> \date January 2021
!> \details
!    This module will derive the nonhydrostatic pressure correction
!     uses petsc with the ILU preconditioner for solving
!
!-----------------------------------------------------------------------

module ocn_compute_nonhydrostatic_pressure

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_constants
  use mpas_log

  use ocn_constants
  use ocn_config

  !PETSC includes for the matrix solve
  use petscsys
  use petscvec
  use petscmat
  use petscpc
  use petscksp

  implicit none
  private
  save

  !--------------------------------------------------------------------
  !
  ! Public member functions
  !
  !--------------------------------------------------------------------

  public :: ocn_nonhydrostatic_pressure_tend


  !Things to do
  ! add PETSC paths here
  ! make sure a simple pass through builds with PETSC stuff
  ! construct PETSC RHS vector
  ! create the mapping of index to a matrix that is nCells+nVertLevels in each direction
  ! create the matrix
  ! separate the preconditioner creation into a separate subroutine to see if it needs to be done each timestep
  ! do the solve, update normal velocity and recompute vertical velocity to be divergence free
  ! need: big routine to come in and then call out to -- preconditioner -- matrix solve -- update u/w

  ! init routine needed to set up index function and init a counter for the pre-condition solve
  contains

  subroutine ocn_nonhydrostatic_pressure_tend(normalVelocityNew, normalVelocityCur, &
      verticalVelocityNew, verticalVelocityCur, layerThickness, layerThicknessEdge, &
      nonHydrostaticPressure)

    real,dimension(:,:), intent(in) :: normalVelocityCur, verticalVelocityCur, &
      layerThickness, layerThicknessEdge

    real,dimension(:,:), intent(inout) :: nonHydrostaticPressure, normalVelocityNew, &
      verticalVelocityNew

    ! calculate pre conditioner -- add some counter to do every so often
    ! build RHS vector and matrix on LHS
    ! do the actual solve
    ! compute new normal velocity from nonhydrostatic pressure
    ! compute new vertical velocity from div of new normal V

  end subroutine ocn_nonhydrostatic_pressure_tend

end module ocn_compute_nonhydrostatic_pressure

