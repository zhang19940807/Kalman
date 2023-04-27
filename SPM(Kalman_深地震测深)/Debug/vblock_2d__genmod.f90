        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VBLOCK_2D__genmod
          INTERFACE 
            SUBROUTINE VBLOCK_2D(NOB,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,V,X0,Z0, &
     &BV,NB)
              INTEGER(KIND=4) :: NOB
              REAL(KIND=4) :: XMIN
              REAL(KIND=4) :: XMAX
              REAL(KIND=4) :: ZMIN
              REAL(KIND=4) :: ZMAX
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: DZ
              REAL(KIND=4) :: V(*)
              REAL(KIND=4) :: X0
              REAL(KIND=4) :: Z0
              REAL(KIND=4) :: BV(*)
              INTEGER(KIND=4) :: NB(*)
            END SUBROUTINE VBLOCK_2D
          END INTERFACE 
        END MODULE VBLOCK_2D__genmod
