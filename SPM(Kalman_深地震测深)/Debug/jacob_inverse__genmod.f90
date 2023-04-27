        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:07 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE JACOB_INVERSE__genmod
          INTERFACE 
            SUBROUTINE JACOB_INVERSE(IDREC,ORDER,IDSRC,V)
              USE GLOBALP
              INTEGER(KIND=4) :: IDREC
              INTEGER(KIND=4) :: ORDER
              INTEGER(KIND=4) :: IDSRC
              REAL(KIND=4) :: V(NNX*NNZ)
            END SUBROUTINE JACOB_INVERSE
          END INTERFACE 
        END MODULE JACOB_INVERSE__genmod
