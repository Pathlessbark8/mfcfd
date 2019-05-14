C The following to set the AMPI_FORTRANCOMPATIBLE
C  if it was placed in the configure command.
#include "ampi/userIF/libConfig.h"

#ifdef AMPI_FORTRANCOMPATIBLE
      subroutine ADTOOL_AMPI_fortranSetupTypes(adouble, areal)
#include "ampi/ampif.h"
        INTEGER adouble
        INTEGER areal
C for now just set this to be scalars of the matching type 
        AMPI_ADOUBLE_PRECISION=MPI_DOUBLE_PRECISION
        adouble=AMPI_ADOUBLE_PRECISION
        AMPI_AREAL=MPI_REAL
        areal=AMPI_AREAL
      end subroutine

      subroutine ADTOOL_AMPI_fortranCleanupTypes(adouble, areal)
#include "ampi/ampif.h"
        INTEGER adouble
        INTEGER areal
C for now do nothing
      end subroutine
#endif
