        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 23 19:40:06 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TIAOGEZI__genmod
          INTERFACE 
            SUBROUTINE TIAOGEZI(XMIN,XMAX,ZMIN,ZMAX,LL,NFT,X,Z,NGR,NF,  &
     &NNX,NNZ,FF,DX,DZ,MJ,MJ1,MJ2,MJ3)
              INTEGER(KIND=4) :: NNZ
              INTEGER(KIND=4) :: NNX
              INTEGER(KIND=4) :: NF
              INTEGER(KIND=4) :: NFT
              INTEGER(KIND=4) :: LL
              REAL(KIND=4) :: XMIN
              REAL(KIND=4) :: XMAX
              REAL(KIND=4) :: ZMIN
              REAL(KIND=4) :: ZMAX
              REAL(KIND=4) :: X(NFT)
              REAL(KIND=4) :: Z(NFT)
              INTEGER(KIND=4) :: NGR(LL,NF)
              INTEGER(KIND=4) :: FF(1:LL+1,1:(NNX-1)*(NNZ-1))
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: DZ
              INTEGER(KIND=4) :: MJ(1:LL+1)
              INTEGER(KIND=4) :: MJ1(1:LL,1:10*NNX)
              INTEGER(KIND=4) :: MJ2(1:LL)
              INTEGER(KIND=4) :: MJ3(NNX-1,LL)
            END SUBROUTINE TIAOGEZI
          END INTERFACE 
        END MODULE TIAOGEZI__genmod
