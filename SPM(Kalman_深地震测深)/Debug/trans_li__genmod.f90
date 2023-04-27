        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANS_LI__genmod
          INTERFACE 
            SUBROUTINE TRANS_LI(NELEMENT,NDATA,NUNK,IROW,IRTN,ICOL0,ICOL&
     &,ICTN)
              INTEGER(KIND=4), INTENT(IN) :: NUNK
              INTEGER(KIND=4), INTENT(IN) :: NDATA
              INTEGER(KIND=4), INTENT(IN) :: NELEMENT
              INTEGER(KIND=4), INTENT(IN) :: IROW(NELEMENT)
              INTEGER(KIND=4), INTENT(IN) :: IRTN(NDATA)
              INTEGER(KIND=4), INTENT(OUT) :: ICOL0(NELEMENT)
              INTEGER(KIND=4), INTENT(OUT) :: ICOL(NELEMENT)
              INTEGER(KIND=4), INTENT(OUT) :: ICTN(NUNK)
            END SUBROUTINE TRANS_LI
          END INTERFACE 
        END MODULE TRANS_LI__genmod
