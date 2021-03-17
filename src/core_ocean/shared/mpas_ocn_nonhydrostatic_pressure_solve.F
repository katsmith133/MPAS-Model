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

#define PETSC_USE_FORTRAN_MODULES 1 
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscmat.h"
!#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petsc.h"

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
  use petsc

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
            ocn_nonhydrostatic_solver_fill_matvec, &
            ocn_nonhydrostatic_solver_finalize, &
            ocn_nonhydrostatic_pressure_solve, &
            ocn_nonhydrostatic_pressue_update_velocity

  !--------------------------------------------------------------------
  !
  ! Private module variables
  !
  !--------------------------------------------------------------------

  double precision :: norm ! norm of solution error
  PetscInt ::  i, j, II, JJ, m, n, its, Istart, Iend
  PetscInt ::  row, col, ione, globalM, first, last
  PetscErrorCode :: petsc_err, ierr2
  PetscMPIInt :: petsc_rank, petsc_size
  PetscBool :: flag
  PetscScalar :: v, one, neg_one
  Vec :: x, b, u ! approximate solution, right hand side vector, exact solution vector
  Mat :: Amat ! Matrix that defines the system
  Mat :: Pmat ! Preconditioner Matrix
  KSP :: ksp ! krylov subspace method context
  PetscRandom :: rctx ! random number generator
  PC :: pc ! preconditioner context
  PCType :: ptype ! flag to set type of preconditioner MAYBE make this a NL option?
  PetscReal :: tol 

  contains

    subroutine ocn_nonhydrostatic_pressure_tend(normalVelocityNew, &
      verticalVelocityNew, layerThickness, layerThicknessEdge, &
      nonhydrostaticPressure, mpi_comm, dt)

    real,dimension(:,:), intent(in) :: normalVelocityNew, verticalVelocityNew, &
      layerThickness, layerThicknessEdge

    real(kind=RKIND), intent(in) :: dt

    real,dimension(:,:), intent(inout) :: nonhydrostaticPressure

    integer, intent(in) :: mpi_comm

    logical :: reuse_preconditioner 

    reuse_preconditioner = .false.

    call ocn_nonhydrostatic_solver_fill_matVec(normalVelocityNew, &
      verticalVelocityNew, layerThickness, layerThicknessEdge, mpi_comm, dt)

    call ocn_nonhydrostatic_pressure_solve(mpi_comm, reuse_preconditioner, nonhydrostaticPressure)

  end subroutine ocn_nonhydrostatic_pressure_tend

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_pressure_update_velocity
  !
  !> \brief  Update the normal and vertical velocity 
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  creates KSP solver context, gets preconditioner and solves
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_pressue_update_velocity(normalVelocityNew, &
      verticalVelocityNew, layerThickness, layerThicknessEdge,           &
      nonhydrostaticPressure, dt)

    real(KIND=RKIND),intent(in) :: dt

    real(KIND=RKIND),dimension(:,:), intent(inout) ::  &
      normalVelocityNew,    &
      verticalVelocityNew

    real(KIND=RKIND),dimension(:,:),intent(in) :: &
      layerThicknessEdge, &
      layerThickness,     &
      nonhydrostaticPressure

    integer :: cell1, cell2, nEdges, iCell, iEdge, nCells, i, k

    real(kind=RKIND), dimension(nVertLevels) :: div_hu

    real(kind=RKIND) :: invAreaCell, r_tmp

    nEdges = nEdgesHalo(2)

    ! loop over edges to update normalVelocity
    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(k, cell1, cell2)
    do iEdge=1,nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       do k = 1,maxLevelEdgeTop(iEdge)
          normalVelocityNew(k,iEdge) = normalVelocityNew(k,iEdge) - dt* &
            (nonhydrostaticPressure(k,cell2) - nonhydrostaticPressure(k,cell1)) / &
            dcEdge(iEdge)*edgeMask(k,iEdge)
       end do
    end do
    !$omp end do
    !$omp end parallel

    nCells = nCellsHalo( 1 )

    !$omp parallel
    !$omp do schedule(runtime) &
    !$omp private(i, r_tmp, iEdge, invAreaCell, k, div_hu)
    do iCell = 1, nCells
       div_hu(:) = 0.0_RKIND
       invAreaCell = 1.0_RKIND / areaCell(iCell)
       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          do k = 1,maxLevelCell(iCell)
             r_tmp = dvEdge(iEdge)*normalVelocityNew(k,iEdge)*invAreaCell
             div_hu(k) = div_hu(k) - layerThicknessEdge(k,iEdge)*edgeSignOnCell(i,iCell)*r_tmp
          end do
       end do

       verticalVelocityNew(maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
       do k=maxLevelCell(iCell),1,-1
          verticalVelocityNew(k,iCell) = verticalVelocityNew(k+1,iCell) - div_hu(k)
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine ocn_nonhydrostatic_pressue_update_velocity

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_pressure_solve_mat
  !
  !> \brief   
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  creates KSP solver context, gets preconditioner and solves
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_pressure_solve(mpi_comm, reuse_preconditioner, nonhydrostaticPressure)

    integer, intent(in) :: mpi_comm
    logical, intent(in) :: reuse_preconditioner
    real(kind=RKIND),dimension(:,:), intent(out) :: nonhydrostaticPressure
    PetscScalar :: x_array(m)
    PetscOffset :: i_x
    Vec :: uu

    integer :: vecSpot, iCell, k, reason


    !Get version number as some function calls are dependent on that
    call KSPCreate(mpi_comm, ksp, petsc_err)
    if(petsc_err .ne. 0) then
       call mpas_log_write("error: petsc create solver context failed, error code = $i", mpas_log_crit, &
             intArgs=(/petsc_err/))
    end if

    if(reuse_preconditioner) then
       call KSPSetReusePreconditioner(ksp, PETSC_TRUE, petsc_err)
    else
       call KSPSetReusePreconditioner(ksp, PETSC_FALSE, petsc_err)
    end if

    call KSPSetOperators(ksp,Amat,Amat,petsc_err)
    if(petsc_err .ne. 0) then
       call mpas_log_write("error: petsc set operaters failed, error code = $i", mpas_log_crit, &
             intArgs=(/petsc_err/))
    end if

    call KSPSetType(ksp, config_nonhydrostatic_solver_type, petsc_err)
    if(petsc_err .ne. 0) then
       call mpas_log_write("error: petsc set type failed, error code = $i", mpas_log_crit, &
             intArgs=(/petsc_err/))
    end if

    if(.not. reuse_preconditioner) then
       call KSPGetPC(ksp,pc,petsc_err)
       if(petsc_err .ne. 0) then
          call mpas_log_write("error: petsc getPC failed, error code = $i", mpas_log_crit, &
             intArgs=(/petsc_err/))
       end if

       call PCSetType(pc,config_nonhydrostatic_preconditioner,petsc_err)
       if(petsc_err .ne. 0) then
          call mpas_log_write("error: petsc PCsetType failed, error code = $i", mpas_log_crit, &
             intArgs=(/petsc_err/))
       end if
    end if

    !now that things are set up, we solve the matrix equation
    call KSPSolve(ksp, b, x, petsc_err)
    if(petsc_err .ne. 0) then
       call mpas_log_write("error: petsc solve failed, error code = $i", mpas_log_crit, &
          intArgs=(/petsc_err/))
    end if

    !check convergence of solution
    call KSPGetConvergedReason(ksp,reason,petsc_err)
    if(reason < 0) then
       call mpas_log_write("error: petsc failed to converge", mpas_log_crit)
    end if

    !Set next initial guess at solution equal to current solution
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,petsc_err)
    if(petsc_err < 0) then
       call mpas_log_write("error: petsc set next initial guess failed, error code = $i", mpas_log_crit, &
          intArgs=(/petsc_err/))
    end if

    !now update the nonhydrostaticPressure array unpack into k,iCell
!    call VecGetOwnershipRange(x,first,last,petsc_err) 
    !$omp parallel
    !$omp do schedule(runtime) private(k,vecSpot)
    do iCell=1,nCellsOwned
       do k=1,maxLevelCell(iCell)
         vecSpot = (mpasToPetscMap(iCell)-1)*nVertLevels + k - 1
          call VecGetValues(x,1,vecSpot,nonhydrostaticPressure(k,iCell),petsc_err)
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine ocn_nonhydrostatic_pressure_solve

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_pressure_fill_matvec 
  !
  !> \brief   fills matrix and vector for pressure solve 
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  Fills matrix on LHS for pressure solve and the RHS vector
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_fill_matVec(normalVelocityNew, &
      verticalVelocityNew, layerThickness, layerThicknessEdge, mpi_comm, dt)

    real(kind=RKIND),dimension(:,:), intent(in) :: normalVelocityNew, verticalVelocityNew
    
    real(kind=RKIND),dimension(:,:), intent(in) :: layerThickness, layerThicknessEdge

    real(kind=RKIND), intent(in) :: dt

    integer, intent(in) :: mpi_comm

  !--------------------------------------------------------------------
  !
  ! subroutine variables
  !
  !--------------------------------------------------------------------

   integer :: jj,cell1, cell2, neighborInd, j, iEdge, iCell, k, matIndex
   PetscScalar :: insertVal
   PetscInt :: i
   real :: termSign, edge_tmp

   jj=maxval(cellsOnEdge(:,1:nCellsOwned))
   !Make sure all entries are zero so there is no unintended overwrite
   call MatZeroEntries(Amat, petsc_err)
   do iCell = 1,nCellsOwned
     !add the points with iCell and k, k+1, k-1i
     !at the top of the domain and bottom (maxLevelCell+1, fill things in with zeros below
     ! should only contain an icell entry to get derivatives right. better BCs now
     ! the conversion to matrix space is i = nVertLevels*(iCell-1) + nVertLevelsi
     do j=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(j, iCell)
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        if(cell1 == iCell) then
          termSign = -1
          neighborInd = cell2
        else
          termSign = 1
          neighborInd = cell1
        end if

        edge_tmp = dvEdge(iEdge)*edgeSignOnCell(j,iCell)*dt/(rho_sw*dcEdge(iEdge)*areaCell(iCell))

        do k=1,maxLevelEdgeTop(iEdge)
          row = (mpasToPetscMap(iCell)-1)*nVertLevels + k - 1 
          col = (mpasToPetscMap(cell2)-1)*nVertLevels + k - 1
          insertVal = edge_tmp*layerThicknessEdge(k,iEdge)
          call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

          col = (mpasToPetscMap(cell1)-1)*nVertLevels + k - 1
          insertVal = -edge_tmp*layerThicknessEdge(k,iEdge)
          call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)
        end do
     end do !End loop over edges

     ! So now do the iCell, k (and k-1 k + 1)
     k = 1
     col = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
     row = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
     ! k value
     insertVal = -dt*(1.0_RKIND / (0.5_RKIND*rho_sw*(layerThickness(k,iCell) + layerThickness(k+1,iCell) + 1.0E-15_RKIND)))
     call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

     ! k + 1 value
     col = (mpasToPetscMap(iCell) - 1)*nVertLevels + (k+1) - 1
     insertVal = dt / (0.5_RKIND*rho_sw*(layerThickness(k,iCell) + layerThickness(k+1,iCell) + 1.0E-15_RKIND))
     call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

     do k=2,maxLevelCell(iCell)-1
        col = (mpasToPetscMap(iCell) - 1)*nVertLevels + (k-1) - 1
        row = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
        insertVal = dt / (0.5_RKIND*rho_sw*(layerThickness(k-1,iCell) + layerThickness(k,iCell)) + 1.0E-15_RKIND)
        call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

        ! k value
        col = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
        insertVal = -dt*(1.0/(0.5_RKIND*rho_sw*(layerThickness(k-1,iCell) + layerThickness(k,iCell) + 1.0E-15_RKIND)) + &
                         1.0/(0.5_RKIND*rho_sw*(layerThickness(k,iCell) + layerThickness(k+1,iCell) + 1.0E-15_RKIND)))
        call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

        ! k+1 value
        col = (mpasToPetscMap(iCell) - 1)*nVertLevels + (k+1) - 1
        insertVal = dt/(0.5_RKIND*rho_sw*(layerThickness(k,iCell) + layerThickness(k+1,iCell) + 1.0E-15_RKIND))
        call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)
     end do !k loop

     k = maxLevelCell(iCell)
     !k-1 value
     col = (mpasToPetscMap(iCell) - 1)*nVertLevels + (k-1) - 1
     row = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
     insertVal = dt / (0.5_RKIND*rho_sw*(layerThickness(k-1,iCell) + layerThickness(k,iCell)) + 1.0E-15_RKIND)
     call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)

     !k value
     col = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1 
     insertVal = -dt*(1.0/(0.5_RKIND*rho_sw*(layerThickness(k-1,iCell) + layerThickness(k,iCell) + 1.0E-15_RKIND)))
     call MatSetValues(Amat, 1, row, 1, col, insertVal, ADD_VALUES, petsc_err)
  end do !iCell loop

  call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,petsc_err)
  call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,petsc_err)

  ! next fill RHS vector 'b'
  call vecZeroEntries(b,petsc_err)

  do iCell=1,nCellsOwned
     do k=1,maxLevelCell(iCell)
        insertVal = verticalVelocityNew(k,iCell) - verticalVelocityNew(k+1,iCell)
        row = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
        call VecSetValues(b, 1, row, insertVal, ADD_VALUES, petsc_err)
     enddo
  enddo

  do iCell=1,nCellsOwned
     do j=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(j,iCell)
        do k=1,maxLevelEdgeTop(iEdge)
           insertVal = layerThicknessEdge(k,iEdge)*dvEdge(iEdge)*edgeSignOnCell(j,iCell)* &
                          normalVelocityNew(k,iEdge)/areaCell(iCell)
           row = (mpasToPetscMap(iCell) - 1)*nVertLevels + k - 1
           call VecSetValues(b, 1, row, insertVal, ADD_VALUES, petsc_err)
        end do
     end do
  end do

  call VecAssemblyBegin(b,petsc_err)
  call VecAssemblyEnd(b,petsc_err)

  end subroutine ocn_nonhydrostatic_solver_fill_matVec

  !--------------------------------------------------------------------
  ! IIIIII
  ! routine ocn_nonhydrostatic_solver_init
  !
  !> \brief   Initializes matricesIIIIII and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>  Initializese the matrix for the nonhydrostatic solve, automatically
  !>  determines what is owned by processor.  Does not fill values.
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_init(mpi_comm, ierr)

    integer, intent(in) :: mpi_comm
    integer, intent(inout) :: ierr
    PetscInt :: major, minor, subminor, release, locRow, locCol, locals(nCellsOwned)
    integer :: commsize, globalCells
    real(kind=RKIND) :: decimals, version
    ISLocalToGlobalMapping :: mapping

    call MPI_comm_size(mpi_comm, commsize, ierr)

    call PetscInitialize(PETSC_NULL_CHARACTER,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("error: petsc initialize failed, error code = $i", mpas_log_crit, &
        intargs=(/petsc_err/))
      ierr = petsc_err
    endif

    call PetscGetVersionNumber(major, minor, subminor, release, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc GetVersionNumber failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    end if

    !Solver routines require version to be greater than 3.5
    decimals = floor(log10(float(minor))) + 1
    version = float(major) + float(minor) / 10.0_RKIND**decimals
    if(version < 3.05) then
      call mpas_log_write("ERROR: Petsc Version required to be greater than 3.5, you are using $r", &
        MPAS_LOG_CRIT, realArgs=(/version/))
      ierr = petsc_err
    end if

    m = nCellsOwned * nVertLevels
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-m',m,flag,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc OptionsGetInt failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatCreate(mpi_comm, Amat, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatCreate failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatSetSizes(Amat, m, m, PETSC_DETERMINE, PETSC_DETERMINE, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetSizes failed, error code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatSetType( Amat, MATAIJ, petsc_err )
    if( commsize == 1) then
      call MatSetType( Amat, MATAIJ, petsc_err )
    else
      call MatSetType( Amat, MATMPIAIJ, petsc_err )
    end if

    call MatSetFromOptions(Amat,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetFromOptions failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatSetUp(Amat,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatSetUp failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatGetOwnershipRange(Amat,Istart,Iend,petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetOwnershipRange failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call MatGetSize(Amat, globalM, globalM, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc MatGetSize failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif


    call MatGetLocalSize(Amat,locRow,locCol,petsc_err)
    call VecCreate(mpi_comm, b, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif


    call VecSetSizes(b, m, globalM, petsc_err)
    if ( commsize == 1 ) then
       call VecSetType(b, VECSEQ, petsc_err)
    else
       call VecSetType(b, VECMPI, petsc_err)
    end if

    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecDuplicate(b, x, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecCreate Duplicate, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecSetSizes(x, m, globalM, petsc_err)
    if(petsc_err .ne. 0) then
      call mpas_log_write("ERROR: Petsc VecSetSizes failed, error_code = $i", MPAS_LOG_CRIT, &
        intArgs=(/petsc_err/))
      ierr = petsc_err
    endif

    call VecGetOwnershipRange(x,first,last,petsc_err)

    end subroutine ocn_nonhydrostatic_solver_init

  !--------------------------------------------------------------------
  ! 
  ! routine ocn_nonhydrostatic_solver_finalize
  !
  !> \brief   Destroys matrices and vectors for PETSC solve
  !> \author  Luke Van Roekel
  !> \date    February 2021
  !> \details
  !>   Destroys matrices and frees memory, finalizes Petsc 
  !
  !--------------------------------------------------------------------

  subroutine ocn_nonhydrostatic_solver_finalize()

    call MatDestroy(Amat,petsc_err)
    call VecDestroy(b,petsc_err)
    call VecDestroy(x,petsc_err)
    call KSPDestroy(ksp,ierr2)
    call PetscFinalize(ierr2)

  end subroutine ocn_nonhydrostatic_solver_finalize  

end module ocn_compute_nonhydrostatic_pressure

