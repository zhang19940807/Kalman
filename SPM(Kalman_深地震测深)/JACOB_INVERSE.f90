SUBROUTINE JACOB_inverse(IDrec,order,IDsrc,V)
    use globalp
    IMPLICIT NONE
    INTEGER :: i,j,k,i1,IDrec,order,IDsrc,line,M1,IBLOCK,IB,ix,iz,np,nn,Intemp1,Intemp2
    INTEGER :: lr(10),i2,Unit
    REAL :: Dpos(3),A1(3),A2(3),MEB(4),NB(4),BV(2,2),EPS = 1E-6,V(nnx*nnz)
    REAL :: RealTemp1,RealTemp2,RealTemp3,vsrc,SN,DT(4),gradt_in_perp,gradt_out_perp
    REAL :: XO,YO,ZO,XI,YI,ZI,XX0,YY0,ZZ0,Z1,Z2,Y1,Y2,X1,X2,X0,Y0,Z0,VO,VI
    REAL :: bw(2),np_x(1:3),np_z(1:3),gradt_in(1:2),gradt_out(1:2),normal(1:2),maxv,maxi
    REAL ::vel_out,vel_in,norm,w,geo_factor,drx,drz,XOO,ZOO,XII,ZII,XRAY(IntPara),ZRAY(IntPara) 
    AAA = 0.0
    NELEA = 0 
    IROW = 0
    maxv = 0.0
    maxi = 0.0
    IF(inv ==1) THEN
        IF(order ==1)  Unit = 22
        IF(order ==2 .or. order ==3) Unit=23
        read(Unit,*) Intemp1,Intemp2
        Read(Unit,*) noderay(order,IDsrc,IDrec)
        Do line=1,noderay(order,IDsrc,IDrec)
            READ(Unit,*) XRAY(line),ZRAY(line)
        ENDDO
        DO line=1,noderay(order,IDsrc,IDrec)-1
            XO = XRAY(line)
            ZO = ZRAY(line)
            XI = XRAY(line+1)
            ZI = ZRAY(line+1)
            XX0 = 0.5*(XO+XI)
            ZZ0 = 0.5*(ZO+ZI) 
            M1=0
            DO  I=1,NNZ-1
                Z1=ZMIN+FLOAT(I-1)*DZ
                Z2=Z1+DZ
                IF((ZZ0.LT.Z1).OR.(ZZ0.GT.Z2)) Cycle
                DO K=1,NNX-1
                    X1=XMIN+FLOAT(K-1)*DX
                    X2=X1+DX
                    IF((XX0.LT.X1).OR.(XX0.GT.X2)) Cycle
                    M1=M1+1
                    IBLOCK=(I-1)*(NNX-1) + K
                    MEB(M1)=IBLOCK
                ENDDO
            ENDDO
            IF(M1.EQ.0)THEN
                WRITE(*,*)'WRONG RAY-PATH !'
                PAUSE
                STOP
            ENDIF
            SN = SQRT((XO-XI)*(XO-XI)+(ZO-ZI)*(ZO-ZI))
            SN = -0.5 * SN
            DO IB=1,M1
                IBLOCK=MEB(IB)
                IF(mod(IBLOCK,( NNX - 1))  == 0) Then
                    ix = nnx - 1
                    iz = IBLOCK/(nnx-1)
                ELSE
                    ix = mod(IBLOCK,(nnx-1))
                    iz = IBLOCK/(nnx-1) +1 
                ENDIF
                nb(1) = ix+(iz-1)*nnx  
                nb(2) = nb(1) + 1
                nb(3) = nb(1) + nnx
                nb(4) = nb(3) + 1
                DO i1=1,2
                    Do i2=1,2
                        bv(i1,i2)=v(nb((i1-1)*2 + i2)) 
                    ENDDO
                ENDDO
                X0 = (ix -1) * dx
                Z0 = (iz -1) * dz
                XOO = XO - X0
                ZOO = ZO - Z0
                XII = XI - X0
                ZII = ZI - Z0
                Call bilinear2D(BV,XOO,ZOO,VO)
                Call bilinear2D(BV,XII,ZII,VI)
                XOO = XOO/dx
                ZOO = ZOO/dz
                XII = XII/dx
                ZII = ZII/dz
                sn = sn /((VO+VI)/2.0)**2  ! -L/(v*v)
                DT(1)=SN*((1.-XOO)*(1.-ZOO)+(1.-XII)*(1.-ZII))
                DT(2)=SN*(XOO*(1.-ZOO)+XII*(1.-ZII))
                DT(3)=SN*(ZOO*(1.-XOO)+ZII*(1.-XII))
                DT(4)=SN*(XOO*ZOO+XII*ZII)
                DO i=1,4
                    AAA(nb(i))=AAA(nb(i))+dt(i)
                ENDDO
            ENDDO  
        ENDDO
    ENDIF
    !--计算走时关于界面深度的偏导数---------------------------------------------------
    IF(order /=1)  THEN
        IF(ini==1) THEN
            read(21,*) np,(lr(i),i=1,np)        
            Do i=1,np  !交点个数循
                DO j=1,3 
                    read(21,*) np_x(j),np_z(j)
                ENDDO 
            !计算入射的单位向量
                gradt_in(1)=np_z(2)-np_z(1)
                gradt_in(2)=np_x(2)-np_x(1)
            !-----------------------------------------------------
                x0=0.5*(np_x(2)+np_x(1))
                z0=0.5*(np_z(2)+np_z(1))

                ix=int((x0-xmin)/dx + eps/2.0)+1
                if(ix == nnx) ix = ix-1
                iz=int((z0-zmin)/dz + eps/2.0)+1
                if(iz == nnz) iz = iz-1
                nb(1) = ix + (iz-1) * nnx
                nb(2) = nb(1)+1
                nb(3) = nb(1)+nnx
                nb(4) = nb(3)+1
                do i1=1,2
                    Do i2 =1,2  
                        bv(i1,i2)=v(nb((i1-1)*2 +i2))
                    ENDDO 
                enddo
                drx = x0 - (ix -1) * dx
                drz = z0 - (iz -1) * dz
                call bilinear2D(bv,drx,drz,vel_in)
            !-----------------------------------------------------
                norm = sqrt(sum(gradt_in**2))
                gradt_in=gradt_in/(vel_in*norm)
                !计算射线与界面交点处的法向量
                call interface_normal(np_x(2),lr(i),normal(1),normal(2))

            gradt_in_perp = dot_product(gradt_in,normal)
            gradt_out(1)=np_z(3)-np_z(2)
            gradt_out(2)=np_x(3)-np_x(2)
            !-----------------------------------------------------
            x0=0.5*(np_x(3)+np_x(2))
            z0=0.5*(np_z(3)+np_z(2))

            ix=int((x0-xmin)/dx+eps/2.0)+1
            if(ix==nnx) ix=ix-1
            iz=int((z0-zmin)/dz+eps/2.0)+1
            if(iz==nnz) iz=iz-1

            nb(1)=ix+nnx*(iz-1)
            nb(2)=nb(1)+1
            nb(3)=nb(1)+nnx
            nb(4)=nb(3)+1
            do i1=1,2
                Do i2 =1,2 
                    bv(i1,i2)=v(nb((i1-1)*2 + i2))
                ENDDO
            ENDDO
            drx = x0 -(ix - 1)*dx - xmin
            drz = z0 -(iz - 1)*dz - zmin
            call bilinear2D(bv,drx,drz,vel_out)
            norm=sqrt(sum(gradt_out**2))
            gradt_out=gradt_out/(vel_out*norm)
            gradt_out_perp = dot_product(gradt_out,normal)
            geo_factor=(gradt_in_perp-gradt_out_perp)*normal(1)   !normal(1)=dot_product(normal,Wz),Wz=(1,0,0)-z方向的单位向量
            ix=int((np_x(2)-xmin)/dx)+1           
            if(ix==nnx)ix=ix-1
            w=(np_x(2)-((ix-1)*dx - xmin))/dx
            bw(1)=1-w
            bw(2)=w
            do k=1,2  
                nn = inv * nnx * nnz + nnx * (lr(i)-1) + ix + k-1
                AAA(nn) = AAA(nn) + geo_factor *bw(k)
            enddo
            ENDDO
        ENDIF
    ENDIF   
    DO k=1,NUNK
        IF( ABS(AAA(k)) > eps ) THEN
            NELEA = NELEA + 1
            IROW(NELEA) = k
        ENDIF
    ENDDO
    Call KF_FOMULAS(order,IDrec,IDrec,IDsrc)
    RETURN
    END SUBROUTINE JACOB_inverse
    !----------------------KF_FOMULAS----------------------------
    Subroutine KF_FOMULAS(order,IDrec,AnothSta,IDsrc)
    use globalp
    IMPLICIT NONE
    Integer :: i,j,k,ii,jj,kk,order,IDrec,AnothSta,IDsrc
    Real :: dt,GcovG,Gdm,realTemp
    Real,allocatable :: Gcov(:),covG(:),covGG(:,:),covGGcov(:,:),Kagain(:)
    allocate(Gcov(NUNK))
    allocate(covG(NUNK),covGG(NUNK,NUNK),covGGcov(NUNK,NUNK))
    allocate(Kagain(NUNK))
    Kagain = 0.0
    Gcov = 0.0
    covG = 0.0
    covGG = 0.0
    covGGcov = 0.0
    dt = TTO(order,IDrec,IDsrc) - TTcal(order,IDrec,IDsrc) - T0(IDsrc)
    Do i=1,NUNK
        Do j=1,NELEA
            Gcov(i) = Gcov(i) + AAA(IROW(j))*cov(IROW(j),i)
            !                   AAA(j) * cov(j,i)
            Kagain(i) = Kagain(i) + cov(i,IROW(j)) * AAA(IROW(j))
        ENDDO
    ENDDO
    GcovG = 0.0
    Gdm = 0.0
    Do i=1,NELEA
        GcovG = GcovG + Gcov(IROW(i)) * AAA(IROW(i))
        Gdm = Gdm + AAA(IROW(i))*dm(IROW(i))
    ENDDO
    GcovG = GcovG + sigma * sigma
    Kagain = Kagain/GcovG
    !dm(:,k+1)= dm(:,k) +  K*(dt(k)  - g(k,:)*dm(:,k)) ;
    realTemp = dt - Gdm
    dm = dm + Kagain * realTemp
    Do i=1,NUNK
        Do j=1,NELEA
            covG(i) = covG(i) + cov(i,IROW(j)) * AAA(IROW(j))
        ENDDO
    ENDDO
    DO i=1,NUNK 
        Do j=1,NELEA
            covGG(i,IROW(j)) = covG(i) * AAA(IROW(j))
        ENDDO
    ENDDO
    DO i=1,NUNK 
        DO j=1,NUNK
            DO ii=1,NELEA
                covGGcov(i,j) = covGGcov(i,j) + covGG(i,IROW(ii)) * cov(IROW(ii),j)
            ENDDO
        ENDDO
    ENDDO
    cov = cov - covGGcov/GcovG
    DEALLOCATE(Gcov,covG,covGG,covGGcov,Kagain)
    END Subroutine KF_FOMULAS

!*****************************************************
! this subroutine returns the upward normal of an interface at a horizontal position
! using 1D bicubic spline interpolation
!*****************************************************
      subroutine interface_normal(x,iface,norm_z,norm_x)
      use globalp
      implicit none
      integer  i,j,ix,iface,ib

       INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

       real x,norm_z,norm_x
       real  u,v,bv(2),h,dhdv,bpv(2),dhdx,norm


       ix = floor((x-xmin)/dx)+1	 
       if(ix==nnx) ix = ix-1
       v = (x-((ix-1)*dx - xmin))/dx
       ib=(iface-1)*nnx+ix

       h=0.0_dp
       dhdv=0.0_dp
	   
        ! bu(1)=1-u
	    ! bu(2)=u
        bv(1)=1-v
	    bv(2)=v

	    ! bpu(1)=-1
	    ! bpu(2)=1
	    bpv(1)=-1
	    bpv(2)=1
          
	    do j=1,2 !x
            !do i=1,2 !y
               h   =h   + bv(j)* zrI(ib + j-1)
               dhdv=dhdv+ bpv(j)*zrI(ib + j-1)
	      !end do
          end do
       dhdx=dhdv/dx

       norm=sqrt(1.0_dp+dhdx**2)
       
       norm_z = 1.0_dp/norm
       norm_x = -dhdx/norm

      end subroutine interface_normal