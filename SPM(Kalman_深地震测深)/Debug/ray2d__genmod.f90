        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RAY2D__genmod
          INTERFACE 
            SUBROUTINE RAY2D(XMIN,XMAX,ZMIN,ZMAX,DX,DZ,NSIDE,NS,XS,ZS,NR&
     &,XR,ZR,V,MSR)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: NS
              REAL(KIND=4) :: XMIN
              REAL(KIND=4) :: XMAX
              REAL(KIND=4) :: ZMIN
              REAL(KIND=4) :: ZMAX
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: DZ
              INTEGER(KIND=4) :: NSIDE
              REAL(KIND=4) :: XS(*)
              REAL(KIND=4) :: ZS(*)
              REAL(KIND=4) :: XR(*)
              REAL(KIND=4) :: ZR(*)
              REAL(KIND=4) :: V(*)
              INTEGER(KIND=4) :: MSR(1:NS,1:NR)
            END SUBROUTINE RAY2D
          END INTERFACE 
        END MODULE RAY2D__genmod
