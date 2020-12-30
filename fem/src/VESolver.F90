MODULE VESolver
    USE Types
    IMPLICIT NONE
    public VESolver_Solve, VESolver_Init

    INCLUDE "mpif.h"

    INTERFACE
        integer(C_INT32_T) FUNCTION VESolver_Get_Local_comm() bind (c,name="vesolver_get_local_comm_f")
            USE iso_c_binding
        END FUNCTION
    END INTERFACE

    INTERFACE
        SUBROUTINE VESolver_Fini() bind (c,name="vesolver_fini")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE

CONTAINS

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Solve(solver, mtype, neq, nnz, Rows, Cols, Values, b, x, res, comm, err)
!--------------------------------------------------------------------------------
        INTEGER :: solver, mtype, neq, nnz, err
        !INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
        !REAL(8), POINTER CONTIG :: Values(:), b(:), x(:)
        INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: Values, x, b
        REAL(KIND=dp) :: res
        INTEGER :: comm

        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
        INTEGER :: ve_master, vl_comm, myrank
        INTEGER :: tag=1, i

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Get_VEMaster_C() bind (c,name="vesolver_get_vemaster")
                USE iso_c_binding
            END FUNCTION
        END INTERFACE

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Get_VL_comm() bind (c,name="vesolver_get_vl_comm_f")
                USE iso_c_binding
            END FUNCTION
        END INTERFACE

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Send_Request_C(solver, mtype, nrow, ncol, nval, res) &
                bind (c,name="vesolver_send_request")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: solver, nrow, ncol, nval
                INTEGER(c_int), VALUE :: mtype
                REAL(c_double), VALUE :: res
            END FUNCTION
        END INTERFACE
!
!        write(*,*) 'INFO:Matrix_A:info:', 1, neq, Rows(neq+1)-1
!        !DO i=1,neq+1
!        DO i=1,10
!          write(*,*) 'INFO:Matrix_A:Rows: ', i, i, Rows(i)
!        END DO
!        !DO i=1,Rows(n+1)-1
!        DO i=1,10
!          write(*,*) 'INFO:Matrix_A:Cols_Value: ', i, Cols(i), Values(i)
!        END DO
!
!        !DO i=1,neq
!        DO i=1,10
!            write(*,*) 'INFO:Vector_b: ', i, i, b(i)
!        END DO
!
!        DO i=neq-10,neq
!            write(*,*) 'INFO:Vector_x: ', i, i, x(i)
!        END DO
!
        ve_master = VESolver_Get_VEmaster_C();
        vl_comm = VESolver_Get_VL_comm();
        CALL MPI_COMM_RANK(vl_comm, myrank, err)

        WRITE( Message, * ) 'INFO: Send Request to VE (rank, solver, mtype, res) = (', &
            ve_master, solver, mtype, res, ').'
        CALL Info( 'VESolver_Solve', Message, Level=5 )
        err = VESolver_Send_Request_C(solver, mtype, neq, neq, nnz, res)

        IF (myrank.eq.0) THEN
            ! Send Matrix A
            CALL MPI_SEND(Values, nnz, MPI_DOUBLE, ve_master, tag, vl_comm, err);
            CALL MPI_SEND(Cols, nnz, MPI_INTEGER, ve_master, tag, vl_comm, err);
            CALL MPI_SEND(Rows, neq+1, MPI_INTEGER, ve_master, tag, vl_comm, err);

            ! Send Matrix b
            CALL MPI_Send(b, neq, MPI_DOUBLE, ve_master, tag, vl_comm, err);

            ! Recive Matrix x
            CALL MPI_RECV(x, neq, MPI_DOUBLE, ve_master, tag, vl_comm, status, err);
        END IF
        CALL MPI_BARRIER(comm)
!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Solve
!--------------------------------------------------------------------------------

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
        nz = A % Rows(n+1)-1

        WRITE( Message, * ) 'INFO: VEDirectSolver called.'
        CALL Info( 'VEDirectSolver', Message, Level=5 )
        CALL VESolver_Solve(0, 17, n, nz, A % Rows, A % Cols, A % Values, b, x, res, A % comm, err)

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
        nz = A % Rows(n+1)-1

        WRITE( Message, * ) 'INFO: VEIterSolver called.'
        CALL Info( 'VEIterSolver', Message, Level=5 )
        CALL VESolver_Solve(100, 17, n, nz, A % Rows, A % Cols, A % Values, b, x, res, A % comm, err)

!-------------------------------------------------------------------------------
    END SUBROUTINE VEIterSolver
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE VEParIterSolver( SourceMatrix, ParallelInfo, DOFs, XVec, &
              RHSVec, Solver, SParMatrixDesc )
!-------------------------------------------------------------------------------
       TYPE (Matrix_t) :: SourceMatrix
       TYPE (ParallelInfo_t) :: ParallelInfo
       INTEGER :: DOFs
       REAL(KIND=dp), DIMENSION(:) :: XVec, RHSVec
       TYPE (Solver_t) :: Solver
       TYPE (SParIterSolverGlobalD_t), POINTER :: SParMatrixDesc
                                                                                                                               
#ifdef PARALLEL_FOR_REAL
       CALL SParIterSolver( SourceMatrix, ParallelInfo, XVec, &
                 RHSVec, Solver, SParMatrixDesc )
#endif
!-------------------------------------------------------------------------------
    END SUBROUTINE VEParIterSolver
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!   VESolver_Init
!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Init(comm, err)
        INTEGER :: comm, err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Init_C(comm) bind (c,name="vesolver_init_f")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: comm
            END FUNCTION
        END INTERFACE

        err = VESolver_Init_C(comm)
    END SUBROUTINE VESolver_Init

END MODULE VESolver
