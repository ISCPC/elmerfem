MODULE VESolver

#include "../config.h"

#ifdef HAVE_VESOLVER
    USE Types
    !USE SParIterSolve
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
    SUBROUTINE VESolver_Solve(solver, mtype, A, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER :: solver, mtype
        TYPE(Matrix_t), TARGET :: A
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        !REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: Values, x, b
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x, b
        REAL(KIND=dp) :: res
        INTEGER :: comm

        INTEGER :: neq, nnz, err
        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
        INTEGER :: ve_master, vl_comm
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
            integer(C_INT32_T) FUNCTION VESolver_Send_Request_C(solver, mtype, nprocs, res) &
                bind (c,name="vesolver_send_request")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: solver, mtype, nprocs
                REAL(c_double), VALUE :: res
            END FUNCTION
        END INTERFACE

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Send_Matrix_Info_C(mtype, nrow, ncol, nval, neq) &
                bind (c,name="vesolver_send_matrix_info")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: mtype
                INTEGER(c_int), VALUE :: nrow, ncol, nval, neq
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

        WRITE( Message, * ) 'INFO: Send Request to VE (rank, solver, mtype, res) = (', &
            ve_master, solver, mtype, res, ').'
        CALL Info( 'VESolver_Solve', Message, Level=5 )

        neq = A % Numberofrows
        nnz = A % Rows(neq+1)-1
        comm = A % comm

        ! vesolver_send_request
        err = VESolver_Send_Request_C(solver, mtype, 1, res)

        ! vesolver_send_matrix_info
        err = VESolver_Send_Matrix_Info_C(mtype, neq, neq, nnz, neq)

        ! vesolver_send_matrix_data
        ! Send Matrix A
        CALL MPI_SEND(A % Values, nnz, MPI_DOUBLE, ve_master, tag, vl_comm, err);
        CALL MPI_SEND(A % Cols, nnz, MPI_INTEGER, ve_master, tag, vl_comm, err);
        CALL MPI_SEND(A % Rows, neq+1, MPI_INTEGER, ve_master, tag, vl_comm, err);

        ! Send Matrix b
        CALL MPI_Send(b, neq, MPI_DOUBLE, ve_master, tag, vl_comm, err);

        ! Recive Matrix x
        CALL MPI_RECV(x, neq, MPI_DOUBLE, ve_master, tag, vl_comm, status, err);
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
        !n = A % Numberofrows
        !nz = A % Rows(n+1)-1

        WRITE( Message, * ) 'INFO: VEDirectSolver called.'
        CALL Info( 'VEDirectSolver', Message, Level=5 )
        !CALL VESolver_Solve(0, 17, n, nz, A % Rows, A % Cols, A % Values, b, x, res, A % comm, err)
        CALL VESolver_Solve(0, 17, A, b, x, res, err)

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
        !n = A % Numberofrows
        !nz = A % Rows(n+1)-1

        WRITE( Message, * ) 'INFO: VEIterSolver called.'
        CALL Info( 'VEIterSolver', Message, Level=5 )
        !CALL VESolver_Solve(100, 17, n, nz, A % Rows, A % Cols, A % Values, b, x, res, A % comm, err)
        CALL VESolver_Solve(100, 17, A, b, x, res, err)

!-------------------------------------------------------------------------------
    END SUBROUTINE VEIterSolver
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_PSolve(solver, mtype, A, b, x, res, err)
!--------------------------------------------------------------------------------
        USE SParIterSolve
        INTEGER :: solver, mtype
        TYPE(Matrix_t), TARGET :: A
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x, b
        REAL(KIND=dp) :: res
        INTEGER :: comm

        INTEGER :: allocstat, nOwned
        INTEGER, DIMENSION(:), POINTER CONTIG :: Cols, Rorder, Owner
        REAL(KIND=dp), DIMENSION(:), POINTER CONTIG :: x1
        INTEGER :: n, neq, nnz, err, nRows
        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
        INTEGER :: ve_master, vl_comm, nprocs
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
            integer(C_INT32_T) FUNCTION VESolver_Send_Request_C(solver, mtype, nprocs, res) &
                bind (c,name="vesolver_send_request")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: solver, mtype, nprocs
                REAL(c_double), VALUE :: res
            END FUNCTION
        END INTERFACE

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Send_Matrix_Info_C(mtype, nrow, ncol, nval, neq) &
                bind (c,name="vesolver_send_matrix_info")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: mtype
                INTEGER(c_int), VALUE :: nrow, ncol, nval, neq
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

        nRows = A % Numberofrows
        nnz = A % Rows(nRows+1)-1

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
 
        WRITE( *, * ) 'INFO: VESolver_PSolve', nRows, neq
        WRITE( Message, * ) 'INFO: Send Request to VE (rank, solver, mtype, res) = (', &
            ve_master, solver, mtype, res, nRows, nnz, neq, ').'
        CALL Info( 'VESolver_PSolve', Message, Level=5 )

        ! vesolver_send_request
        CALL MPI_COMM_SIZE(A % comm, nprocs, err)
        err = VESolver_Send_Request_C(solver, mtype, nprocs, res)

        ! vesolver_send_matrix_info
        err = VESolver_Send_Matrix_Info_C(mtype, nRows , nRows, nnz, neq)

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

        ! Send Matrix A
        CALL MPI_SEND(A % Values, nnz, MPI_DOUBLE, ve_master, tag, vl_comm, err);
        CALL MPI_SEND(Cols, nnz, MPI_INTEGER, ve_master, tag, vl_comm, err);
        CALL MPI_SEND(A % Rows, nRows+1, MPI_INTEGER, ve_master, tag, vl_comm, err);
        CALL MPI_SEND(Rorder, neq, MPI_INTEGER, ve_master, tag, vl_comm, err);

        ! Send Matrix b
        CALL MPI_Send(b, nRows, MPI_DOUBLE, ve_master, tag, vl_comm, err);
        DEALLOCATE(Cols)
        DEALLOCATE(Rorder)

        ! Recive Matrix x
        CALL MPI_RECV(x1, neq, MPI_DOUBLE, ve_master, tag, vl_comm, status, err);
        CALL MPI_BARRIER(A % Comm)

        ! Distribute solution
        DO i=1,nRows
            x(i)=x1(A % Gorder(i))
        END DO
        DEALLOCATE(x1)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_PSolve
!--------------------------------------------------------------------------------

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
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
        !REAL(KIND=dp),  TARGET :: x(*), b(*)
        TYPE(Matrix_t), TARGET :: A
        INTEGER :: n, nz, err
        REAL(8) :: res

        res = 1.0e-8

        WRITE( Message, * ) 'INFO: VEParIterSolver called.'
        CALL Info( 'VEParIterSolver', Message, Level=5 )
        CALL VESolver_PSolve(100, 4096+17, A, b, x, res, err)

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
#endif

END MODULE VESolver
