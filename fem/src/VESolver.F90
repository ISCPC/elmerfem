MODULE VESolver

#include "../config.h"

#ifdef HAVE_VESOLVER
#define VES_MODE_GATHER_ON_VH  0
#define VES_MODE_GATHER_ON_VE  1
#define VES_MODE_SYMMETRIC     2

#define VESOLVER_HS            0
#define VESOLVER_BICGSTAB2   100
#define VESOLVER_DUMMY      9999

    USE Types
    USE Lists
    USE VESolverAPI
    IMPLICIT NONE

CONTAINS

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
        CALL VESolver_Solve(VESOLVER_HS, n, A % Values, A % Rows, A % Cols, b, x, res, err)
        CALL VESolver_Deactivate()

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
    SUBROUTINE VEParSolver( A, x, b, Solver )
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
        INTEGER :: solverId

        LOGICAL :: GotIt
        CHARACTER(LEN=MAX_NAME_LEN) :: Method

        ! Set variable
        nRows = A % Numberofrows
        nnz = A % Rows(nRows+1)-1

        ! Check parameter
        Method = ListGetString( Solver % Values, 'Linear System Solver',GotIt )
        IF ( .NOT. GotIt ) Method = 'veiterative'
        CALL Info('VEParIterSolver','Solver type: '//TRIM(Method),Level=5)
        SELECT CASE(Method)
          CASE('veiterative')
            solverId = VESOLVER_BICGSTAB2
          CASE('vedirect')
            solverId = VESOLVER_HS
          CASE DEFAULT
            solverId = VESOLVER_DUMMY
        END SELECT

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
        CALL VESolver_PSolve(mode, solverId, &
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
    END SUBROUTINE VEParSolver
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE VEParDirectSolver( A, x, b, Solver )
!------------------------------------------------------------------------------
    USE SParIterSolve
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), TARGET :: x(*), b(*)

    INTEGER mtype, n, neq, nrhs, ierror
    INTEGER i, j, k, nz
    LOGICAL :: Found, matsym, matpd

    INTEGER, ALLOCATABLE :: Owner(:)
    INTEGER allocstat, nOwned, nl, nt, ip
    INTEGER :: rind, lrow, rptr, rsize, lind, tind
    INTEGER :: nzutd, nhalo

    LOGICAL :: Factorize, FreeFactorize
    CHARACTER(LEN=MAX_NAME_LEN) :: threads, mat_type

    REAL(KIND=dp), ALLOCATABLE :: aa(:)
    INTEGER, ALLOCATABLE  :: ia(:), ja(:)
    REAL*8 :: res

    INTEGER, DIMENSION(:), POINTER CONTIG :: iperm, Order
    REAL(KIND=dp), ALLOCATABLE :: dbuf(:), RHS(:)

    REAL(KIND=dp), DIMENSION(:), POINTER CONTIG :: b1, x1
    INTEGER :: nprocs

    CALL Info( 'VEParDirectSolver:', 'Called.', Level=5 )

    ! Set matrix type for Pardiso
    !mtype = HS_UNSYMMETRIC

    ! Set up continuous numbering for the whole computation domain
    n = SIZE(A % ParallelInfo % GlobalDOFs)
    ALLOCATE(A % Gorder(n), Owner(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('VEParDirectSolver', &
                    'Memory allocation for global numbering failed')
    END IF
    CALL ContinuousNumbering(A % ParallelInfo, A % Perm, A % Gorder, Owner, nOwn=nOwned)

    ! Compute the number of global dofs
    CALL MPI_ALLREDUCE(nOwned, neq, &
                       1, MPI_INTEGER, MPI_SUM, A % Comm, ierror)
    DEALLOCATE(Owner)

    ! Find bounds of domain
    nl = A % Gorder(1)
    nt = A % Gorder(1)
    DO i=2,n
        ! NOTE: Matrix is structurally symmetric
        rind = A % Gorder(i)
        nl = MIN(rind, nl)
        nt = MAX(rind, nt)
    END DO

    ! Allocate temp storage for global numbering
    ALLOCATE(Order(n), iperm(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('VEParDirectSolver', &
                    'Memory allocation for global numbering failed')
    END IF

    ! Sort global numbering to build matrix
    Order(1:n) = A % Gorder(1:n)
    DO i=1,n
        iperm(i)=i
    END DO
    CALL SortI(n, Order, iperm)

    ! Allocate storage for CPardiso matrix
    nhalo = (nt-nl+1)-n
    nz = A % Rows(A % NumberOfRows+1)-1
!    IF (mtype.eq.HS_UNSYMMETRIC) THEN
        ALLOCATE(ia(nt-nl+2), &
                 ja(nz+nhalo), &
                 aa(nz+nhalo), &
                STAT=allocstat)
!    ELSE
!        nzutd = ((nz-n)/2)+1 + n
!        ALLOCATE(ia(nt-nl+2), &
!                 ja(nzutd+nhalo), &
!                 aa(nzutd+nhalo), &
!                 STAT=allocstat)
!    END IF
    IF (allocstat /= 0) THEN
        CALL Fatal('VEParDirectSolver', &
                   'Memory allocation for CPardiso matrix failed')
    END IF

    ! Build distributed CRS matrix
    ia(1) = 1
    lrow = 1      ! Next row to add
    rptr = 1      ! Pointer to next row to add, equals ia(lrow)
    lind = Order(1)-1 ! Row pointer for the first round

    ! Add rows of matrix 
    DO i=1,n
      ! Skip empty rows
      tind = Order(i)
      rsize = (tind-lind)-1

      DO j=1,rsize
        ia(lrow+j)=rptr
      END DO
      lrow = lrow + rsize

      ! Add next row
      rind = iperm(i)
      lind = A % rows(rind)
      tind = A % rows(rind+1)
      rsize = tind-lind
      DO j=lind, tind-1
        ja(rptr+(j-lind))=A % Gorder(A % Cols(j))
        aa(rptr+(j-lind))=A % values(j)
      END DO

      ! Sort column indices
      CALL SortF(rsize, ja(rptr:rptr+rsize), aa(rptr:rptr+rsize))

      ! Set up row pointers
      rptr = rptr + rsize
      lrow = lrow + 1
      ia(lrow) = rptr

      lind = Order(i) ! Store row index for next round
    END DO

    ! Deallocate temp storage
    DEALLOCATE(Order, iperm)

#if 0
    CALL MPI_Barrier(A % Comm)
    write(*,*) 'INFO: nl, nt, ni, nj = ', nl, nt, nt-nl+2, nz+nhalo
    write(*,*) 'INFO:Number of Rows, nonzero = ', A % NumberOfRows, &
        SIZE(A % ParallelInfo % GlobalDOFs), A % Rows(n+1)-1, neq, ia(i+1)-1
    write(*,*) 'INFO:Location(Col,Row) = ', A % Cols(1), A % Rows(1)
    write(*,*) 'INFO:Location(Col,Row) = ', A % Cols(2), A % Rows(2)
    write(*,*) 'INFO:Values = ', A % Values(1), A % Values(2)
    DO i=1,nt-nl+1
      DO j=ia(i),ia(i+1)-1
        write(*,*) 'INFO: (i,j,v) = ', i, j, ia(i), ja(j), aa(j)
      END DO
    END DO
    write(*,*) 'Done'
    CALL flush(6)
    CALL MPI_Barrier(A % Comm)
    CALL Fatal('VEParDirectSolver','Stop due to matrix output')
    RETURN
#endif

#if 1
    ALLOCATE( b1(nt-nl+1) )
    DO i=1,A % NumberOfRows
      ip = A % Gorder(i) - nl + 1
      b1(ip) = b(i)
    END DO
    !
    ! Call common solver function
    !
    ALLOCATE( x1(neq) )
    CALL MPI_COMM_SIZE(A % Comm, nprocs, ierror)
    CALL VESolver_Activate(A % comm, nprocs, ierror)
    CALL VESolver_PSolve_dcsr(VES_MODE_SYMMETRIC, VESOLVER_HS, neq, aa, ia, ja, nl, nt, &
        b1, x1, res, ierror)
    CALL VESolver_Deactivate()

    ! Distribute solution
    DO i=1,A % NumberOfRows
        x(i)=x1(A % Gorder(i))
    END DO

   ! free buffers
    DEALLOCATE(b1)
    DEALLOCATE(x1)
#else
    ! Compute factorization if necessary
    ! initialize handle
    !mtype = HS_UNSYMMETRIC
    CALL PHS_init_handle(A % HShandle, neq, neq, mtype, HS_DCSR, HS_DIST_A, nl, nt, A % Comm, ierror)

    ! Perform preprocess
    CALL PHS_preprocess_rd(A % HShandle, ia, ja, aa, ierror);

    ! Perform factorization
    CALL PHS_factorize_rd(A % HShandle, ia, ja, aa, ierror)

    ! ------------------------------------------
    ALLOCATE( dbuf(neq) )
    DO i=1,A % NumberOfRows
      ip = A % Gorder(i)
      dbuf(ip) = b(i)
    END DO
    ALLOCATE( RHS(neq) )
    CALL MPI_ALLREDUCE( dbuf, RHS, neq, MPI_DOUBLE_PRECISION, MPI_SUM, A % Comm, ierror )

    ! Perform solve
    nrhs      = 1 ! Use only one RHS
    res       = 1.0e-8
    CALL PHS_solve_rd(A % HShandle, ia, ja, aa, nrhs, RHS, RHS, res, ierror);

    ! Distribute the solution to all:
    ! -------------------------------
    CALL MPI_BCAST( RHS, neq, MPI_DOUBLE_PRECISION, 0, A % Comm, ierror )

    ! Select the values which belong to us:
    ! -------------------------------------
    DO i=1,A % NumberOfRows
      ip = A % Gorder(i)
      x(i) = RHS(ip)
    END DO

    DEALLOCATE(RHS)
    DEALLOCATE(dbuf)
#endif

!------------------------------------------------------------------------------
  END SUBROUTINE VEParDirectSolver
!------------------------------------------------------------------------------

#endif /* HAVE_VESOLVER */

END MODULE VESolver
