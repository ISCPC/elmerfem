MODULE VESolver

#include "../config.h"

#ifdef HAVE_VESOLVER
#define VES_MODE_GATHER_ON_VH  0
#define VES_MODE_GATHER_ON_VE  1
#define VES_MODE_SYMMETRIC     2

#define VESOLVER_HS            0
#define VESOLVER_BICGSTAB2   100
#define VESOLVER_DUMMY      9999

#define VES_MATINFO_TAG     1000
#define VES_MATDATA_TAG     1001
#define VES_MATGATHER_TAG   1002

    USE Types
    USE Lists
    IMPLICIT NONE

    INTERFACE
        SUBROUTINE VESolver_Deactivate() bind (c,name="vesolver_deactivate")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE

CONTAINS

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Activate(comm, nprocs, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: comm, nprocs
        INTEGER, intent(out) :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Activate_C(comm, nprocs) bind (c,name="vesolver_activate_f")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: comm, nprocs
            END FUNCTION
        END INTERFACE

        err = VESolver_Activate_C(comm, nprocs)
!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Activate
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Solve(solver, neq, Values, Rows, Cols, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: solver, neq
        REAL(8), DIMENSION(:) :: Values, b, x
        INTEGER, DIMENSION(:) :: Rows, Cols
        !REAL(8), DIMENSION(:), TARGET CONTIG :: Values, b, x
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTEGER :: i, nnz
        INTEGER :: ves_comm
        INTEGER :: tag=VES_MATDATA_TAG

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Solve_C(solver, mtype, neq, &
                Rows, Cols, Values, b, x, res) bind (c,name="vesolver_solve")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: solver, mtype, neq
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        nnz = Rows(neq+1)-1
#if 0
        write(*,*) 'INFO:Matrix_A:info:', 1, neq, Rows(neq+1)-1
        DO i=1,neq+1
        !DO i=1,10
          write(*,*) 'INFO:Matrix_A:Rows: ', i, i, Rows(i)
        END DO
        DO i=1,nnz
        !DO i=1,10
          write(*,*) 'INFO:Matrix_A:Cols_Value: ', i, Cols(i), Values(i)
        END DO

        DO i=1,neq
        !DO i=1,10
            write(*,*) 'INFO:Vector_b: ', i, i, b(i)
        END DO

        DO i=1,neq
            write(*,*) 'INFO:Vector_x: ', i, i, x(i)
        END DO
#endif
        err = VESolver_Solve_C(solver, 17, neq, Rows, Cols, Values, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Solve
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_PSolve(mode, solver, neq, nRows, Values, Rows, Cols, Rorder, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: mode, solver, neq, nRows
        REAL(8), DIMENSION(:), intent(in) :: Values, b, x
        INTEGER, DIMENSION(:), intent(in) :: Rows, Cols, Rorder
        !REAL(8), DIMENSION(:), TARGET CONTIG :: Values, b, x
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTEGER :: i, nnz
        INTEGER :: ves_comm
        INTEGER :: tag=VES_MATDATA_TAG
        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_PSolve_C(mode, solver, mtype, neq, &
                nRows, Rows, Cols, Values, Rorder, b, x, res) bind (c,name="vesolver_psolve")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: mode, solver, mtype, neq, nRows
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols, Rorder
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        nnz = Rows(nRows+1)-1
#if 0
        write(*,*) 'INFO:Matrix_A:info:', 1, neq, nRows, nnz
        DO i=1,nRows+1
        !DO i=1,10
          write(*,*) 'INFO:Matrix_A:Rows: ', i, Rows(i)
        END DO
        DO i=1,nnz
        !DO i=1,10
          write(*,*) 'INFO:Matrix_A:Cols_Value: ', i, Cols(i), Values(i)
        END DO

        DO i=1,nRows
        !DO i=1,10
            write(*,*) 'INFO:Vector_b: ', i, b(i)
        END DO

        DO i=1,neq
            write(*,*) 'INFO:Rorder: ', i, Rorder(i)
        END DO
#endif
        err = VESolver_PSolve_C(mode, solver, 17, neq, nRows, &
            Rows, Cols, Values, Rorder, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_PSolve
!--------------------------------------------------------------------------------

!================================================================================
!
! API for elmer
!
!================================================================================
!--------------------------------------------------------------------------------
    SUBROUTINE VEDirectSolver( A, x, b, Solver )
!--------------------------------------------------------------------------------
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
        !REAL(KIND=dp), TARGET :: x(*), b(*)
        TYPE(Matrix_t), TARGET :: A
        INTEGER :: n, nz, err
        REAL(8) :: res

        !res = 1.0e-10
        res = 1.0e-8
        n = A % Numberofrows

        WRITE( Message, * ) 'INFO: VEDirectSolver called.'
        CALL Info( 'VEDirectSolver', Message, Level=5 )
        CALL VESolver_Activate(A % Comm, 1, err)
        CALL VESolver_Deactivate()
        CALL VESolver_Solve(VESOLVER_HS, n, A % Values, A % Rows, A % Cols, b, x, res, err)

!-------------------------------------------------------------------------------
    END SUBROUTINE VEDirectSolver
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VEIterSolver( A, x, b, Solver )
!--------------------------------------------------------------------------------
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
        !REAL(KIND=dp),  TARGET :: x(*), b(*)
        TYPE(Matrix_t), TARGET :: A
        INTEGER :: n, nz, err
        REAL(8) :: res

        !res = 1.0e-10
        res = 1.0e-8
        n = A % Numberofrows

        WRITE( Message, * ) 'INFO: VEIterSolver called.'
        CALL Info( 'VEIterSolver', Message, Level=5 )
        CALL VESolver_Activate(A % Comm, 1, err)
        CALL VESolver_Solve(VESOLVER_BICGSTAB2, n, A % Values, A % Rows, A % Cols, b, x, res, err)
        CALL VESolver_Deactivate()

!-------------------------------------------------------------------------------
    END SUBROUTINE VEIterSolver
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!    SUBROUTINE VEParIterSolver( SourceMatrix, ParallelInfo, DOFs, XVec, &
!              RHSVec, Solver, SParMatrixDesc )
!-------------------------------------------------------------------------------
!       TYPE (Matrix_t) :: SourceMatrix
!       TYPE (ParallelInfo_t) :: ParallelInfo
!       INTEGER :: DOFs
!       REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec
!       TYPE (Solver_t) :: Solver
!       TYPE (SParIterSolverGlobalD_t), POINTER :: SParMatrixDesc

!-------------------------------------------------------------------------------
    SUBROUTINE VEParIterSolver( A, x, b, Solver )
!-------------------------------------------------------------------------------
        USE SParIterSolve
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
        !REAL(KIND=dp),  TARGET :: x(*), b(*)
        TYPE(Matrix_t), TARGET :: A

        INTEGER :: allocstat, nOwned
        INTEGER, DIMENSION(:), POINTER CONTIG :: Cols, Rorder, Owner
        REAL(KIND=dp), DIMENSION(:), POINTER CONTIG :: x1
        INTEGER :: n, neq, err, i, mode
        INTEGER :: nRows, nnz
        REAL(8) :: res=1.0e-8

        LOGICAL :: GotIt
        CHARACTER(LEN=MAX_NAME_LEN) :: Method

        ! Set variable
        nRows = A % Numberofrows
        nnz = A % Rows(nRows+1)-1

        ! Check parameter
        Method = ListGetString( Solver % Values, 'Linear System Parallelize Mode',GotIt )
        IF ( .NOT. GotIt ) Method = 'GatherOnVH'

        CALL Info('VEParIterSolver','Using mode: '//TRIM(Method),Level=5)
        SELECT CASE( Method )
          CASE( 'gatheronve' )
            mode = VES_MODE_GATHER_ON_VE

          CASE( 'symmetric' )
            mode = VES_MODE_SYMMETRIC

          CASE DEFAULT
            mode = VES_MODE_GATHER_ON_VH
        END SELECT

        ! Set up continuous numbering for the whole computation domain
        n = SIZE(A % ParallelInfo % GlobalDOFs)
        ALLOCATE(A % Gorder(n), Owner(n), STAT=allocstat)
        IF (allocstat /= 0) THEN
             CALL Fatal('VEParIterSolver', &
                        'Memory allocation for VEParSolver global numbering failed')
        END IF
        CALL ContinuousNumbering(A % ParallelInfo, A % Perm, A % Gorder, Owner, nOwn=nOwned)

        ! Compute the number of global dofs
        CALL MPI_ALLREDUCE(nOwned, neq, 1, MPI_INTEGER, MPI_SUM, A % Comm, err)
        DEALLOCATE(Owner)

        ! vesolver_send_matrix_data_distributed
        ALLOCATE(Cols(nnz), Rorder(neq), x1(neq), STAT=allocstat)
        IF (allocstat /= 0) THEN
             CALL Fatal('VEParIterSolver', &
                        'Memory allocation for VEParSolver Cols, Rorder failed')
        END IF

        DO i=1,nnz
            Cols(i) = A % Gorder(A % Cols(i))
        END DO

        Rorder = 0
        DO i=1,nRows
            Rorder(A % Gorder(i)) = i
        END DO

        !
        ! Call common solver function
        !
        CALL VESolver_Activate(A % comm, 1, err)
        CALL VESolver_PSolve(mode, VESOLVER_BICGSTAB2, &
            neq, nRows, A % Values, A % Rows, Cols, Rorder, &
            b, x1, res, err)
        CALL VESolver_Deactivate()

        ! Distribute solution
        DO i=1,nRows
            x(i)=x1(A % Gorder(i))
        END DO

        ! free buffers
        DEALLOCATE(Cols)
        DEALLOCATE(Rorder)
        DEALLOCATE(x1)
!-------------------------------------------------------------------------------
    END SUBROUTINE VEParIterSolver
!-------------------------------------------------------------------------------

#endif

END MODULE VESolver
