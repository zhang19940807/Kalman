        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RAY2D_REF__genmod
          INTERFACE 
            SUBROUTINE RAY2D_REF(XMIN,XMAX,ZMIN,ZMAX,DX,DZ,NSIDEX,NSIDEY&
     &,NS,NS0,NR0,XS,ZS,NF,NNF,XR,ZR,V,VS,NNFF1,NNUM1,LU1,LU2,LU3,NDC)
              INTEGER(KIND=4) :: NNUM1
              INTEGER(KIND=4) :: NNFF1
              INTEGER(KIND=4) :: NNF
              INTEGER(KIND=4) :: NF
              REAL(KIND=4) :: XMIN
              REAL(KIND=4) :: XMAX
              REAL(KIND=4) :: ZMIN
              REAL(KIND=4) :: ZMAX
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: DZ
              INTEGER(KIND=4) :: NSIDEX
              INTEGER(KIND=4) :: NSIDEY
              INTEGER(KIND=4) :: NS
              INTEGER(KIND=4) :: NS0
              INTEGER(KIND=4) :: NR0
              REAL(KIND=4) :: XS(*)
              REAL(KIND=4) :: ZS(*)
              REAL(KIND=4) :: XR(1:NNF,1:NF)
              REAL(KIND=4) :: ZR(1:NNF,1:NF)
              REAL(KIND=4) :: V(*)
              REAL(KIND=4) :: VS(*)
              INTEGER(KIND=4) :: LU1(NNFF1,NNUM1)
              INTEGER(KIND=4) :: LU2(NNFF1)
              INTEGER(KIND=4) :: LU3(NNFF1,NNUM1)
              INTEGER(KIND=4) :: NDC(*)
            END SUBROUTINE RAY2D_REF
          END INTERFACE 
        END MODULE RAY2D_REF__genmod
