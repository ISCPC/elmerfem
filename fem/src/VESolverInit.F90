MODULE VESolverInit

#include "../config.h"

#ifdef HAVE_VESOLVER
    IMPLICIT NONE
    public VESolver_Init, VESolver_Fini

    INTERFACE
        SUBROUTINE VESolver_Fini() bind (c,name="vesolver_fini")
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

#endif

END MODULE VESolverInit
