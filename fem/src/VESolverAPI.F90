!/*****************************************************************************/
! *
! *  VESolver API Fortran module
! *
! *  Copyright 2021 - , Shunji Uno <shunji_uno@iscpc.jp>
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
MODULE VESolverAPI

    IMPLICIT NONE
    public VESolver_Init, VESolver_Fini

    INTERFACE
        SUBROUTINE VESolver_Fini() bind (c,name="vesolver_fini")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE

    INTERFACE
        SUBROUTINE VESolver_Deactivate() bind (c,name="vesolver_deactivate")
            USE iso_c_binding
        END SUBROUTINE
    END INTERFACE


CONTAINS

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_Init(comm, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: comm
        INTEGER, intent(out) :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_Init_C(comm) bind (c,name="vesolver_init_f")
                USE iso_c_binding
                INTEGER(c_int), VALUE :: comm
            END FUNCTION
        END INTERFACE

        err = VESolver_Init_C(comm)
!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Init
!--------------------------------------------------------------------------------

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

        !WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        err = VESolver_Solve_C(solver, 17, neq, Rows, Cols, Values, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_Solve
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
    SUBROUTINE VESolver_PSolve(mode, solver, neq, nRows, Values, Rows, Cols, order, b, x, res, err)
!--------------------------------------------------------------------------------
        INTEGER, intent(in) :: mode, solver, neq, nRows
        REAL(8), DIMENSION(:), intent(in) :: Values, b, x
        INTEGER, DIMENSION(:), intent(in) :: Rows, Cols, order
        !REAL(8), DIMENSION(:), TARGET CONTIG :: Values, b, x
        !INTEGER, DIMENSION(:), TARGET CONTIG :: Rows, Cols
        REAL(8) :: res
        INTEGER :: err

        INTERFACE
            integer(C_INT32_T) FUNCTION VESolver_PSolve_C(mode, solver, neq, &
                nRows, Rows, Cols, Values, order, b, x, res) bind (c,name="vesolver_psolve")
                USE iso_c_binding
                INTEGER(c_int), VALUE, intent(in) :: mode, solver, neq, nRows
                REAL(c_double), DIMENSION(*), intent(in) :: Values, b, x
                INTEGER(c_int), DIMENSION(*), intent(in) :: Rows, Cols, order
                REAL(c_double), VALUE, intent(in) :: res
            END FUNCTION
        END INTERFACE

        !WRITE(*,*) 'INFO: Send Request to VE (solver, res) = (', solver, res, ').'
        err = VESolver_PSolve_C(mode, solver, neq, nRows, &
            Rows, Cols, Values, order, b, x, res)

!--------------------------------------------------------------------------------
    END SUBROUTINE VESolver_PSolve
!--------------------------------------------------------------------------------

END MODULE VESolverAPI
