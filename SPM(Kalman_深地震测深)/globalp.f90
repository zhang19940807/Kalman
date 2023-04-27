MODULE globalp
    IMPLICIT NONE
    REAL,ALLOCATABLE  :: switch(:)
    INTEGER,PARAMETER :: IntPara =  1000
    INTEGER :: NUNK,NELEA,nnx,nnz,inv,ini
    REAL :: sigma,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,DMX,DMZ
    REAL,ALLOCATABLE :: TTO(:,:,:),TTcal(:,:,:),cov(:,:),dm(:),T0(:)
    REAL,ALLOCATABLE :: AAA(:),XR1(:,:),ZR1(:,:),ZRI(:)
    INTEGER,ALLOCATABLE :: noderay(:,:,:),IROW(:)
    Contains
    SUBROUTINE bilinear2D(nv,dsx,dsz,biv)
    IMPLICIT NONE
    INTEGER :: i,j,k
    REAL :: RealTemp1,RealTemp2,RealTemp3,dsx,dsz,biv
    REAL, DIMENSION(2,2) :: nv
    REAL :: produ
    biv=0.0
    DO i=1,2
        DO k=1,2
            RealTemp1 = 1.0-ABS(((i-1)*dx-dsx)/dx)
            RealTemp3 = 1.0-ABS(((k-1)*dz-dsz)/dz)
            produ = RealTemp1 * RealTemp3
            biv=biv+nv(i,k)*produ
        ENDDO
    ENDDO
    END SUBROUTINE bilinear2D
    END MODULE globalp
    