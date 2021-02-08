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
  use ocn_mesh

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

  public :: ocn_nonhydrostatic_pressure_tend, &
            ocn_nonhydrostatic_solver_init,   &
            ocn_nonhydrostatic_pressure_build_matvec

  !--------------------------------------------------------------------
  !
  ! Private module variables
  !
  !--------------------------------------------------------------------

  double precision :: norm ! norm of solution error
  PetscInt ::  i, j, II, JJ, m, n, its
  PetscInt ::  Istart, Iend, ione, globalM
  PetscErrorCode :: petsc_err
  PetscMPIInt :: petsc_rank, petsc_size
  PetscBool :: flag
  PetscScalar :: v, one, neg_one
  Vec :: x, b, u ! approximate solution, right hand side vector, exact solution vector
  Mat :: A ! Matrix that defines the system
  KSP :: ksp ! krylov subspace method context
  PetscRandom :: rctx ! random number generator
  PC :: pc ! preconditioner context
  PCType :: ptype ! flag to set type of preconditioner MAYBE make this a NL option?
  PetscReal :: tol 

  !Things to do
  ! construct PETSC RHS vector
  ! create the mapping of index to a matrix that is nCells+nVertLevels in each direction
  ! create the matrix
  ! separate the preconditioner creation into a separate subroutine to see if it needs to be done each timestep
  ! do the solve, update normal velocity and recompute vertical velocity to be divergence free
  ! need: big routine to come in and then call out to -- preconditioner -- matrix solve -- update u/w

  ! init routine needed to set up index function and init a counter for the pre-condition solve
  contains

    subroutine ocn_nonhydrostatic_pressure_tend(normalvelocitynew, normalvelocitycur, &
      verticalvelocitynew, verticalvelocitycur, layerthickness, layerthicknessedge, &
      nonHydrostaticPressure, mpi_comm)

    real,dimension(:,:), intent(in) :: normalVelocityCur, verticalVelocityCur, &
      layerThickness, layerThicknessEdge

    real,dimension(:,:), intent(inout) :: nonHydrostaticPressure, normalVelocityNew, &
      verticalVelocityNew

    integer, intent(in) :: mpi_comm

    ! calculate pre conditioner -- add some counter to do every so often
    ! build RHS vector and matrix on LHS
    ! do the actual solve
    ! compute new normal velocity from nonhydrostatic pressure
    ! compute new vertical velocity from div of new normal V

  end subroutine ocn_nonhydrostatic_pressure_tend

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_pressure_build_matvec 
  !
  !> \brief   fills matrix and vector for pressure solve 
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  Fills matrix on LHS for pressure solve and the RHS vector
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_pressure_build_matvec(normalVelocityNew,  &
      normalVelocityCur, verticalVelocityNew, verticalVelocityCur,        &
      layerThickness, layerthicknessedge, mpi_comm)

    real,dimension(:,:), intent(in) :: normalVelocityCur, verticalVelocityCur, &
      layerThickness, layerThicknessEdge

    real,dimension(:,:), intent(inout) :: normalVelocityNew, verticalVelocityNew

    integer, intent(in) :: mpi_comm



  end subroutine ocn_nonhydrostatic_pressure_build_matvec

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_solver_init
  !
  !> \brief   Fills matrices and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>   Fills the matrix for the nonhydrostatic solve
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_fill_matVec(normalVelocityNew,  &
      normalVelocityCur, verticalVelocityNew, verticalVelocityCur,        &
      layerThickness, layerthicknessedge, mpi_comm)

    real,dimension(:,:), intent(in) :: normalVelocityCur, verticalVelocityCur, &
      layerThickness, layerThicknessEdge

    real,dimension(:,:), intent(inout) :: normalVelocityNew, verticalVelocityNew

    integer, intent(in) :: mpi_comm

  !--------------------------------------------------------------------
  !
  ! subroutine variables
  !
  !--------------------------------------------------------------------

   integer :: iEdge, iCell, k, matIndex

   do iCell = 1,nCells
     !add the points with iCell and k, k+1, k-1i
     !at the top of the domain and bottom (maxLevelCell+1, fill things in with zeros below
     ! should only contain an icell entry to get derivatives right.
      do iEdge=1,nEdgesOnCell(iCell)
   
    !need to loop over cells and then edges on cells with mapping of
    ! iSpot = iCell*(nVertLevels - 1) + nVertLevels
    !Loop over cells and then edges within each to get the row done
    !

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_solver_init
  !
  !> \brief   Initializes matrices and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  Initializese the matrix for the nonhydrostatic solve, automatically
  !>  determines what is owned by processor.  Does not fill values.
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_init()

    call PetscInitialize(PETSC_NULL_CHARACTER,petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc initialize failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    m = nCellsOwned + nVertLevels
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-m',m,flag,petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc OptionsGetInt failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatCreate(PETSC_COMM_WORLD, A, petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatCreate failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatSetSizes(A, m, m, PETSC_DETERMINE, PETSC_DETERMINE, petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetSizes failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatSetFromOptions(A,petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetFromOptions failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatSetUp(A,petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetUp failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatGetOwnershipRange(A,Istart,Iend,petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetOwnershipRange failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call MatGetSize(A, globalM, globalM, petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetSize failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/))
    endif

    call VecCreate(PETSC_COMM_WORLD, b, petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/)
    endif

    call VecSetSizes(b, m, globalM, petsc_error)
    if(petsc_error .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_error/)
    endif

  end subroutine ocn_nonhydrostatic_solver_init
end module ocn_compute_nonhydrostatic_pressure

