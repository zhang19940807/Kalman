!c-------------------------------------------------------------------c
!c     This program is doing 2D inversion                            c
!c     with the first P (or S) arrival times by a MSPM               c
!c                                                                   c
!c     by Chao-ying Bai,May 07,2006                                  c
!c     Chang'an University,Xi'an, China                              c                    
!c-------------------------------------------------------------------c
      PROGRAM TOMOGRAPHY
      use globalp
	IMPLICIT NONE

	CHARACTER*80 LAB
	REAL ALPH,ALMDA,RMSDA,RMSMM,DEP,VUP,VDOWN,A0,RealTemp,drx
      INTEGER NFIL,NVA,ITSMAX,IWRONG,NSIDEX,NSIDEY,NS,&
      &NS0,NR0,ITS,IV,NDATA,NP,IE,IB,I,J,K,IX,IS,NNFF,NNUM,NNF,IY,NF,NGG,NNDATA,NPP, itype,wave
      REAL RMST(50),RMSM(50)    
      REAL,ALLOCATABLE::XM(:),ZM(:),VS(:),V(:),V0(:,:),XS(:), uncertain_v(:),&
     &      ZS(:),DZZ(:,:),V2(:,:),VF(:),V1(:),DDZZ(:),uncertain_interface(:)
     Real,allocatable :: XS0(:),ZS0(:),XR0(:),ZR0(:),ZRNODE(:,:),XRNODE(:,:)
     INTEGER,ALLOCATABLE::MSR(:,:),LU1(:,:),LU2(:),LU3(:,:),NDC(:),Trriger(:,:)
      
      OPEN(1,FILE='input/2d.inp')
      OPEN(3,FILE='2dray2.out')
!c---------------------------------------------------
!c     input constrained parameters                 
!c---------------------------------------------------
      WRITE(*,*)'       INPUT:'
      WRITE(*,*)'            routine definition '
      READ(1,200)LAB
      READ(1,*)INV,ini,NFIL,NVA,ALPH,ALMDA,ITSMAX,itype,wave,sigma
      IWRONG=0
      IF(NVA.LT.2)IWRONG=1
      IF((ALPH.LT.0.0).OR.(ALPH.GT.0.5))IWRONG=1
      IF(ALMDA.LT.0.0)IWRONG=1
      IF(ITSMAX.LT.0) IWRONG=1
      IF(IWRONG.EQ.1)THEN
      WRITE(*,*)'    wrong with routine definition !'
      ENDIF    
      
!c---------------------------------------------------------
!c     input the ray grid parameters       
!c---------------------------------------------------------
      WRITE(*,*)'            ray grid parameters'
      READ(1,200)LAB
      READ(1,*)XMIN,XMAX,ZMIN,ZMAX,DX,DZ,NSIDEX
      NSIDEY = NSIDEX
      IF(XMAX.LE.XMIN)IWRONG=1
      IF(ZMAX.LE.ZMIN)IWRONG=1
      IF(DX.LE.0.0)IWRONG=1
      IF(DZ.LE.0.0)IWRONG=1
      IF(NSIDEX.LE.0)IWRONG=1
	IF(NSIDEY.LE.0)IWRONG=1
      IF(IWRONG.EQ.1)THEN
      WRITE(*,*)'        wrong with ray grid input !'
      ENDIF  
      NNZ=INT((ZMAX-ZMIN)/DZ+0.5)+1
      NNX=INT((XMAX-XMIN)/DX+0.5)+1
      NP=NNZ*NNX
      NPP=(NNZ-1)*(NNX-1)
      DMX=DX/FLOAT(NSIDEX+1)
      DMZ=DZ/FLOAT(NSIDEY+1)

!c---------------------------------------------------
!c     input background velocity         
!c---------------------------------------------------
      ALLOCATE(V1(1:NNZ),VF(1:NNX),V2(1:NNZ,1:NNX))
!C      WRITE(*,*)'            background velocity'
!C      READ(1,200)LAB      
!C      READ(1,*)VUP,VDOWN,(V1(I),I=1,NNZ)
!C      DO I=1,NNZ
!C       IF(V1(I).LE.0.0)IWRONG=1
!C      ENDDO
!C      IF(IWRONG.EQ.1)THEN
!C      WRITE(*,*)'wrong with background velocity !'
!C      ENDIF
      
!c---------------------------------------------------------
!c     input source positions
!c---------------------------------------------------------
      WRITE(*,*)'            source positions'
      READ(1,200)LAB
      READ(1,*)NS0
	ALLOCATE(NDC(NS0)); NDC = 1
      allocate(XS0(NS0),ZS0(NS0),T0(NS0))
      DO I=1,NS0
      READ(1,*) XS0(I),ZS0(I)
      IF((XS0(I).LT.XMIN).OR.(XS0(I).GT.XMAX))IWRONG=1
      IF((ZS0(I).LT.ZMIN).OR.(ZS0(I).GT.ZMAX))IWRONG=1
      ENDDO
      T0 = 0.0
      IF(IWRONG.EQ.1)THEN
      WRITE(*,*)'       wrong source position !'
      ENDIF      
      
!c---------------------------------------------------------
!c     input receiver positions
!c---------------------------------------------------------
      WRITE(*,*)'            receiver positions'
      READ(1,200)LAB
      READ(1,*)NR0
      allocate(XR0(NR0),ZR0(NR0))
      allocate(Trriger(NS0,NR0))
      DO I=1,NR0
      READ(1,*)XR0(I),ZR0(I)
      IF((XR0(I).LT.XMIN).OR.(XR0(I).GT.XMAX))IWRONG=1
      IF((ZR0(I).LT.ZMIN).OR.(ZR0(I).GT.ZMAX))IWRONG=1
      ENDDO
      IF(IWRONG.EQ.1)THEN
      WRITE(*,*)'       wrong source position !'
      ENDIF 
!c----------------------------------------------------------
!c     sort out the actual total sources in reflected 
!c     or transformed arrivals
!c----------------------------------------------------------
      NS=NS0+NR0
	ALLOCATE(XS(1:NS),ZS(1:NS))
      xs(1:ns) = [xs0(1:ns0),xr0(1:nr0)]
      zs(1:ns) = [zs0(1:ns0),zr0(1:nr0)]
      
!c-----------------------------------------------------------
!c     source-receiver pairs
!c-----------------------------------------------------------     
      ALLOCATE(MSR(1:NS0,1:NR0)); MSR = 1
!C      !WRITE(*,*)'            source-geophone pair'
!C      !READ(1,200)LAB
!C      !DO I=1,NS0
!C      !  READ(1,"(2000i1)")(MSR(I,J),J=1,NR0)
!C      !ENDDO
      NNDATA=count(MSR/=0)
	NDATA=NNDATA
      
!c----------------------------------------------------------
!c     input starting P velocity model
!c----------------------------------------------------------
      ALLOCATE(V0(1:NNZ,1:NNX),VS(1:NNX*NNZ),V(1:NNX*NNZ))
	ALLOCATE(ZM(1:NNZ),XM(1:NNX))
      WRITE(*,*)'            P velocity model'
      open(2,file="input/veln.txt")
      DO K=1,NNX
            DO I =1,NNZ 
                  READ(2,*)XM(K),ZM(I),V0(I,K)
            ENDDO
      ENDDO
      close(2)
      IV=0
      DO 12 I=1,NNZ
      DO 12 K=1,NNX
      IV=IV+1
      V(IV)=V0(I,K)
      IF(V(IV).LT.0.0) then
	  IWRONG=1
	end if
  12  CONTINUE
      
      IF(IWRONG.EQ.1)THEN
      WRITE(*,*)'         wrong value with veolcity model !'
      STOP
      ENDIF
      vs = v / sqrt(3.0)
!c----------------------------------------------------------
!c     input starting S velocity model
!c----------------------------------------------------------


!c----------------------------------------------------------
!c     input Reflected wave type
!c----------------------------------------------------------
      
	READ(1,200)LAB
	WRITE(*,*) "            reflected wave type"
      READ(1,*)NNFF,NNUM
	ALLOCATE(LU1(NNFF,NNUM),LU2(NNFF),LU3(NNFF,NNUM))
	DO I=1,NNFF
	READ(1,*)IX,LU2(I),(LU1(I,K),K=1,NNUM)
      END DO	
	READ(1,200)LAB
	READ(1,*)NNFF,NNUM
	DO I=1,NNFF
	READ(1,*)IX,IY,(LU3(I,K),K=1,NNUM)
	END DO	
      if(itype==2) then
	  NDATA=NDATA*NNFF
      else if(itype==3) then
        NDATA=NDATA*NNFF+NDATA
      end if
       close(1)
!c----------------------------------------------
!c     input starting reflected line
!c----------------------------------------------
      WRITE(*,*)'            reflected line position'
      OPen(1,file="input/Interace.txt")
      READ(1,*)NNF,NF
      ALLOCATE (XR1(1:NNF,1:NF))
      ALLOCATE (ZR1(1:NNF,1:NF))
      ALLOCATE (ZRI(NNF*NF))
	ALLOCATE (DDZZ(1:NF)) 
	ALLOCATE (DZZ(1:NNF,1:NF))
	DO J=1,NNF
	READ(1,*)
	DO K=1,NF
	READ(1,*)XR1(J,K),ZR1(J,K)
	END DO
	END DO 
	ZR1=ABS(ZR1)
      close(1)
      IV =0 
      DO i=1,NNF
            DO j=1,nnx
                  IV = IV +1
                  ZRI(IV) = ZR1(i,j)
            ENDDO
      ENDDO
!c----------------------------------------------------
!c     input observed time data
!c      if do inverse routine         
!c-----------------------------------------------------

	NGG=INT(FLOAT(NP*NDATA)*0.3)
!c-----------------------------------------------------     
!c     end of the data input and begin loop here
!c-----------------------------------------------------
      OPEN(111,file="output/traveltime_difference")
      IF(inv==1) THEN
            allocate(noderay(wave,NS0,NR0))
      ENDIF
      ITS=0
  20  ITS=ITS+1
      WRITE(*,"(a4,i2,',')")'ITS=',ITS
!C-------------------------------------------------------
!C     DIRECT WAVE TRACING
!C-------------------------------------------------------
       if(itype/=2) CALL RAY2D(XMIN,XMAX,ZMIN,ZMAX,DX,DZ,&
     &           NSIDEX,NS0,XS0,ZS0,NR0,XR0,ZR0,V,MSR)
!c---------------------------------------------------------
!c---------------------------------------------------------
!c     do the ray trace by Moser method
!c---------------------------------------------------------
       if(itype/=1) CALL RAY2D_REF(XMIN,XMAX,ZMIN,ZMAX,DX,&
     &    DZ,NSIDEX,NSIDEY,NS,NS0,NR0,XS,ZS,NF,NNF,XR1,ZR1,V,VS,&
     &    NNFF,NNUM,LU1,LU2,LU3,NDC)
!c---------------------------------------------------
!c     do inversion
!c---------------------------------------------------   
      IF(INV + INI < 1) GOTO 600
                IF(ITS == 1) THEN
            NUNK = inv * nnx * nnz + ini * NNF * nnx 
            ALLOCATE(TTO(wave,NR0,NS0),TTcal(wave,NR0,NS0))
            OPEN(1,file="input/observe_time.txt")
            DO k=1,wave
                  DO i=1,NS0
                        READ(1,*) (TTO(k,j,i),j=1,NR0)
                  ENDDO
            ENDDO
            close(1)
        ENDIF
        open(1,file="output\2dtime-d.txt") 
        DO i=1,NS0
            READ(1,*) (TTcal(1,j,i),j=1,NR0)
        ENDDO
        close(1)
        IF(wave >1 .or. ini /= 0) THEN
            OPEN(1,file = "output\2dtime-d2.dat")
            DO i=2,wave
                  DO k=1,NS0
                        Read(1,*) (TTcal(i,j,k),j=1,NR0)
                  enddo
            ENDDO
            close(1)
        ENDIF
        RMSDA = sum(ABS(TTO - TTcal))/(wave * NS0 * NR0)
        WRITE(*,"(a,g0,a,g0)") ' RMST= ',RMSDA
        WRITE(111,*) ITS,RMSDA
        IF(ITS == 1) THEN
            allocate(ZRNODE(NNF,nnx),XRNODE(NNF,nnx))
            allocate(AAA(NUNK))
            allocate(uncertain_v(NP),cov(NUNK,NUNK))
            allocate(uncertain_interface(NNF*NF))
            allocate(dm(NUNK))
            allocate(IROW(NUNK))
            cov = 0.0
            IF(inv ==1) THEN
                  open(1,file="input/uncertain_v(X_Z).txt")
                  IV = 0
                  DO i=1,nnz
                        DO j=1,nnx
                              IV = IV +1 
                              Read(1,*) RealTemp,RealTemp, uncertain_v(IV)
                        ENDDO
                  ENDDO
                  DO i=1,NP
                        cov(i,i) = uncertain_v(i) ** 2
                  ENDDO
                  deallocate(uncertain_v)
                  close(1)
            ENDIF
            IF(ini ==1) THEN
                  open(1,file="input/uncertain_interface.txt")
                  DO i=1,NNF*NNX
                        Read(1,*) uncertain_interface(i)
                  ENDDO
                  DO i=1,NNF*NNX
                        cov(NUNK - NNF*NNX +i,NUNK - NNF*NNX +i) = uncertain_interface(i)**2
                  ENDDO
                  deallocate(uncertain_interface)
                  close(1)
            ENDIF
        ENDIF
        dm = 0.0
!c--------------------------------------------------
!c     calculate the Jacobian matrix and inverse       
!c--------------------------------------------------
      OPEN(21,file="temp/intersection.txt")    
      OPEN(22,file="output/2dpath-d.txt") 
      OPEN(23,file="output/2dpath-d2.dat") 
      DO k=1,wave
            DO i=1,NS0
                  DO j=1,NR0
                        Call JACOB_inverse(j,k,i,V)
                  ENDDO
            ENDDO
      ENDDO
      close(21)
      close(22)
      close(23)
      IF(inv ==1) V = V + dm(1:NP)
      IF(ini ==1) THEN
            Do i=1,NNF
                  Do j=1,NNX
                        ZRNODE(i,j) = ZR1(i,(j-1)*(NSIDEY+1) +1) + dm( inv * NP +(i-1) * NNX +j)
                  ENDDO
            ENDDO
            DO i=1,NNF
                  DO j=1,(NNX-1) * (NSIDEX +1) +1
                        ix = INT((j-1)*DMX/dx) +1
                        IF(ix == nnx) ix = ix -1
                        drx = (j-1)*DMX - (ix-1) * dx
                        drx = drx /dx
                        ZR1(i,j) = ZRNODE(i,ix) *(1.0 - drx) + ZRNODE(i,ix+1) * drx
                  enddo
            ENDDO
      ENDIF
!c--------------------------------------------------
!c     output velocity model/reflect line after each inversion   
!c--------------------------------------------------     


!c---------------------------------------------------   
!     inverse loop constrained here
!c---------------------------------------------------
      IF(ITS.LT.ITSMAX)GO TO 20
      close(111)
      IF(inv==1) THEN
            OPEN(1,file="output/new_veln_uncertain.txt")
            IV = 0
            Do i=1,nnz
                  Do j=1,nnx
                        IV = IV +1
                        WRITE(1,*) (j-1)*dx,(i-1)*dz,sqrt(cov(IV,IV))
                  ENDDO
            ENDDO
            close(1)
      ENDIF
      IF(ini==1) THEN
            OPEN(1,file="output/new_interface_uncertain.txt")
            IV = 0
            Do i=1,NNF
                  WRITE(1,*) "interface:",I
                  Do j=1,nnx
                        IV = IV +1
                        WRITE(1,*) (j-1)*dx,sqrt(cov(IV,IV))
                  ENDDO
            ENDDO
            close(1)
      ENDIF


!c-----------------------------------------------------
!c     filter the data by 2D Alpha filter
!c-----------------------------------------------------
       IF(NFIL.EQ.1) THEN
            IF(inv ==1) THEN
                  DO I=1,NNZ
                        VF(:)=V((I-1)*nnx +1 : I*nnx)
                        CALL FILT2(NNX,NVA,ALPH,VF)
                        V((I-1)*nnx +1 : I*nnx) = VF(:)
                  ENDDO
            ENDIF
            IF(ini ==1) THEN
                  DO i=1,NNF
                        CAll FILT2(NF,NVA,ALPH,ZR1(i,:))
                  ENDDO
            ENDIF
       END IF
       IF(inv ==1) THEN
            OPEN(1,file ="output/newveln.txt")
            IV =0
            DO i=1,nnz
                  DO j=1,nnx
                        IV = IV +1
                        v0(i,j) = V(IV)
                  ENDDO
            ENDDO
            DO i=1,nnx
                  DO j=1,nnz
                        Write(1,*)   (i-1)*dx ,(j-1)*dz,v0(j,i)
                  ENDDO
            ENDDO
            close(1)
      ENDIF
      IF(ini ==1) THEN
            OPEN(1,file="output/newinterface.bln")
            DO i=1,NNF
                  WRITE(1,*) NF,1
                  Do j=1,NF
                        WRITE(1,*) XR1(i,j),ZR1(i,j)
                  ENDDO
            ENDDO
            close(1) 
      ENDIF
!c---------------------------------------------------  
!c     output the convergence curve
!c---------------------------------------------------          
  30  FORMAT('--- ITS:',I5,' RMS(T):',F12.6,' -----')
  32  FORMAT(100(F7.2))
  34  FORMAT('-------------- CONVERGENCE -----------')
  36  FORMAT(1X,I3,1X,F16.8,1X,F16.8) 
      
	         
!C  40  DEALLOCATE(V,VS,TTO,MSR,ZZR,XR1,ZR1,ZZ,T0,XS,ZS,
!C     &	       TT1,XXR,LU1,LU2,DZZ,LU3,NDC,DDZZ)
     
!c-----------------------------------------------------
!c     output file notice
!c-----------------------------------------------------    
600   WRITE(*,*)
      WRITE(*,*)'   OUTPUT                        '
      WRITE(*,*)'   velocity structure:'
      WRITE(*,*)'                           2dray2.out'
      WRITE(*,*)'   final ray-path:'
      WRITE(*,*)'                           2dpath-d2.txt'
      WRITE(*,*)'   final travel time:' 
      WRITE(*,*)'                           2dtime-d2.txt'
      

	call cpu_time(alph)
	write(*,*) 'cpu_time:',alph
	write(3,*)
	write(3,*) 'cpu_time:',alph
      CLOSE(3)      
200   FORMAT(80A)  
      pause
      END PROGRAM TOMOGRAPHY

! c-------------------------------------------------------------------c
! c     This subroutine performs 2D ray tracing for any               c  
! c     volocity model. The irregular network method was              c
! c     used to compute the primary reflected traveltimes             c 
! c     and related raypaths of a grided model.                       c
! c                                                                   c
! c     Entries: (1) INP.................P or S wave;                 c
! c              (2) XMIN,XMAX,ZMIN,ZMAX............... range;        c
! c              (3) DX,DZ..................cell size of grid;        c
! c              (4) NSIDE........number of secondary nodes;          c
! c              (5) NS0,NR0.......actual source and receiver number; c
! c              (6) NS,XS,ZS.....acting source number and location;  c
! c              (7) NR,XR,ZR.........geophone number and locations;  c
! c              (8) NF,XR,ZR.....reflected point and location;       c
! c              (9) IFN,DF...reflected line number and sample length;c
! c              (8) V(*).................volocity parameters.        c
! c                                                                   c
! c     Results: (1) 2D traveltime  file:     2Dtime-pp2.dat;         c
! c                                        or 2dtime-ss2.dat;         c
! c              (2) 2D  raypath coordinates: 2Dpath-pp2.dat;         c
! c                                        or 2dpath-ss2.dat;         c
! c-------------------------------------------------------------------c 
      SUBROUTINE RAY2D_REF(XMIN,XMAX,ZMIN,ZMAX,DX,DZ,NSIDEX, &
     &    NSIDEY,NS,NS0,NR0,XS,ZS,NF,NNF,XR,ZR,V,VS,NNFF1,NNUM1,&
     &    LU1,LU2,LU3,NDC)

     use globalp,only : inv,noderay,ini,DMX,DMZ
      IMPLICIT NONE  
	    
      CHARACTER*80 OPP4
	INTEGER NSIDEX,NSIDEY,NS,NS0,NR0,NF,NNF,NNFF1,NNUM1,NNX,NNZ, &
     &        NP1,NP2,NPP,NP,NODE,NODE3,NODE1,I,J,K,KK,JJ,NBLOCK,IBP, &
     &        IS1,K1,KJ,IR1,NFT,NNCE,KJU,L,II,KGH,ILJ,LJ,JN,KLK,KKI,KKK, &
     &        JP,KL,KB,JOI,JJ1,NN222,IS,KJH,NNF1,IB1,IKJ,LIH,IS0, &
     &        IEER,M1,IB,NU,KN,NUU,IBLOCK,I1,NETP,J1,NOB,NOB1,IP2,IP1, &
     &        IDONE,IY,IT,IT1,IT2,IT3,IR,LL,N2,II1,IKT,IS01,IRE, &
     &        K0,I4,KY,M2,IRAY,I2,IU,ITR,ID
	REAL VFU,XR0,ZR0
	REAL    XMIN,XMAX,ZMIN,ZMAX,DX,DZ,E1,ZI,XK,XJ,ZJJ,XKK, &
     &        TR11,TR12,E,X0,Z0,A1,A3,D,Z1,Z2,X1,X2,XI,X11,X12,ZM1,ZM2, &
     &        XM1,XM2,TBMIN,TD,TTMIN,XP1,ZP1,V1,XP2,ZP2,V2,XX,ZZ,DD,V12, &
     &        T,T12,T21,XX1,XX2,ZZ1,ZZ2,XX0,ZZ0,TT_MIN,IR2,IR3
      REAL(KIND=8) TIMEST,TIMEEN
      REAL  XS(*),ZS(*),BV(4),V(*),VS(*),XR(1:NNF,1:NF),ZR(1:NNF,1:NF)

	INTEGER LU2(NNFF1),NDC(*),NB(1:4)
	INTEGER LU1(NNFF1,NNUM1),LU3(NNFF1,NNUM1)
      REAL,ALLOCATABLE:: X(:),Z(:),TT(:),xtest(:),ytest(:),TB(:),TT1(:)
	INTEGER,ALLOCATABLE::GKK(:,:),MJ2(:),MJ1(:,:),NBNZ(:,:), &
     &MJ3(:,:),ITANG(:),MEB6(:,:,:),MJ(:),NCELL_Z(:,:,:),&
     &ISII(:),NCE(:),NCE1(:),ISII1(:,:),ISII2(:,:),M23(:),M111(:), &
     &ITANG1(:),NGS(:),ICD(:),MEB0(:),MEB1(:),MEB2(:), &
     &MEB4(:,:),MEB5(:,:),ICD0(:,:),NGR(:,:),ICD1(:,:)
	TIMEST=0.0    
	CALL CPU_TIME(TIMEST) 

      NNX=INT((XMAX-XMIN)/DX+0.5)+1
      NNZ=INT((ZMAX-ZMIN)/DZ+0.5)+1
      
      NP1=((NNX-1)*(NSIDEX+1)+1)*NNZ
      NP2=(NNZ-1)*NSIDEY*NNX
      NPP=NP1+NP2+NS+NF
      NP=(NNZ-1)*(NNX-1)
!C      WRITE(*,*)'total cells=',NP,'total nodes=',NPP
 

	E1=MAX(DMX,DMZ)
	E1=E1/2.0
      
      ALLOCATE (X(NPP+NNF*NF))
      ALLOCATE (Z(NPP+NNF*NF))
      NODE=0
!c-----Creat primary node------------------   
      DO 222 I=1,NNZ
       ZI=ZMIN+FLOAT(I-1)*DZ
       IF(ZI.GT.ZMAX) GOTO 222
       DO 188 K=1,NNX
        XK=XMIN+FLOAT(K-1)*DX
        IF(XK.GT.XMAX) GOTO 188
      NODE=NODE+1
      X(NODE)=XK
      Z(NODE)=ZI
 188  CONTINUE
 222  CONTINUE
	NODE3=NODE       !NODE3 primary node
! !C	WRITE(*,*)'CHUSHIHUA'
! c-----Insert the secondary nodes-------------         
      DO 22 I=1,NNZ
       ZI=ZMIN+FLOAT(I-1)*DZ
       DO 20 J=1,NNX
        XJ=XMIN+FLOAT(J-1)*DX
        DO 8 JJ=1,NSIDEY+1
         IF(JJ.EQ.1) THEN
          DO 6 KK=2,NSIDEX+1
           ZJJ=ZI
           XKK=XJ+FLOAT(KK-1)*DMX
           IF(XKK.GT.XMAX)GO TO 6
           NODE=NODE+1
           X(NODE)=XKK
           Z(NODE)=ZJJ
   6     CONTINUE
         ELSE
          XKK=XJ
          ZJJ=ZI+FLOAT(JJ-1)*DMZ          
          IF(ZJJ.GT.ZMAX)GO TO 8
          NODE=NODE+1
          X(NODE)=XKK
          Z(NODE)=ZJJ
         ENDIF 
   8    CONTINUE 
  20  CONTINUE  
  22  CONTINUE
 !     OPEN(789,FILE='YT.TXT')
 !     DO I=1,NNF
	!  DO J=1,NF
	!	WRITE(789,*)XR(I,J),-ZR(I,J)
	!END DO
	!END DO
	!CLOSE(789)

! c---------------------------------------------------c
! c     loop over for the sources                     c
! c---------------------------------------------------c
      NBLOCK=(NNZ-1)*(NNX-1)
      IBP=2*NSIDEX+2*NSIDEY+4 !ibpÿ�������ı��ϵ��ܽڵ���
      
      ALLOCATE (NGS(NS))
      ALLOCATE (NGR(NNF,NF))

! C      ALLOCATE (NCELL(NBLOCK,IBP+500))
! C      ALLOCATE (NBN(NBLOCK))
      ALLOCATE (ICD(NPP+NNF*NF))
      ALLOCATE (ICD0(NS0,NR0))
      ALLOCATE (ICD1(10,NF))
	ALLOCATE (ISII1(1:10,1:NF),M23(10),M111(1:10),ISII2(1:10,1:NF)) 
! c      ALLOCATE (MEB3(NS,NF,NBLOCK))
!	ALLOCATE (MEB4(NNF,NF,NBLOCK+200),MEB5(NR,NBLOCK+200))
! C	ALLOCATE (TR11(NF),TR12(NS0,NNF,NF))
      ALLOCATE (TT(NPP+NNF*NF),TT1(NF))
      
      ALLOCATE (ISII(1:NR0)) 
! C      ALLOCATE (TRR(NS0,NNF,NR))
	ALLOCATE(GKK(1:NNF+1,1:NBLOCK),MJ(1:NNF+1),MJ2(1:NNF),MJ1(1:NNF,1:10*NNX),MJ3(NNX-1,NNF))
	GKK=0
	ngr=0
	ngs=0
! C	ncell=0
! C	nbn=0
	icd=0
	icd0=0
	icd1=0
! C	meb0=0
! C	meb1=0 
! C	meb2=0
	!meb4=0
	!meb5=0 
	tr11=0
	tr12=0
! C	tb=0
	tt=0
! C	trr=0 
	TT1=0

! c----------------------------------------------------c
! c     (2) find the grid No. for source & geophones   c
! c----------------------------------------------------c
      E=MIN(DMX,DMZ)/100.0
      DO 30 I=1,NS
      X0=XS(I)
      Z0=ZS(I)
      NGS(I)=0
      DO 25 J=1,NODE
      A1=(X0-X(J))
      A3=(Z0-Z(J))
      D=SQRT(A1*A1+A3*A3)
      IF(D.GT.E)GOTO 25
      NGS(I)=J
      GO TO 30
  25  CONTINUE
  30  CONTINUE
      K=0
      DO 32 I=1,NS
      IS1=NGS(I)
      IF(IS1.NE.0)GOTO 32
	
      K=K+1
      K1=NODE+K
      NGS(I)=K1
      X(K1)=XS(I)
      Z(K1)=ZS(I)
  32  CONTINUE
      NODE1=NODE+K
      DO KJ=1,NNF     !9000

      DO 355 I=1,NF
      X0=XR(KJ,I)
      Z0=ZR(KJ,I)
      NGR(KJ,I)=0
      DO 27 J=1,NODE
      A1=(X0-X(J))
      A3=(Z0-Z(J))
! c	write(*,*)'hhhhhhhhh',a1,a3
      D=SQRT(A1*A1+A3*A3)
      IF(D.GT.E)GO TO 27
      NGR(KJ,I)=J
! C	WRITE(*,*)J
      GO TO 355
  27  CONTINUE
  355  CONTINUE
   
      K=0
      DO 33 I=1,NF
      IR1=NGR(KJ,I)
      IF(IR1.NE.0)GO TO 33
      K=K+1
      K1=NODE1+K
      NGR(KJ,I)=K1
! C	WRITE(*,*)'TTTTTTTTTT',K1
      X(K1)=XR(KJ,I)
      Z(K1)=ZR(KJ,I)
  33  CONTINUE
	NODE1=NODE1+K
      END DO                !900

	NODE=NODE1

	NFT=0
	NFT=NPP+NNF*NF
!C	WRITE(*,*)'FENQU--COMPUTE'
	CALL TIAOGEZI(XMIN,XMAX,ZMIN,ZMAX,NNF,NFT,X,Z,NGR,NF,NNX,NNZ,GKK,&
     &DX,DZ,MJ,MJ1,MJ2,MJ3)             !�ϸ񽫸��ӷ���
!	DO I=1,NNF+1
!	  WRITE(*,*) MJ(I)
!	END DO
!	PAUSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC		
	NNCE=-1E5
	DO I=1,NNF
	IF(NNCE.LT.MJ(I))THEN
	NNCE=MJ(I)
	END IF
	END DO
      ALLOCATE(NCELL_Z(NNF,NNCE,5*NSIDEX+2*NSIDEY+10),NBNZ(1:NNF,NNCE))
	ALLOCATE(NCE(6*NSIDEX+10),NCE1(6*NSIDEX+10))
      ALLOCATE (MEB4(NF,NNCE*4),MEB5(NR0,NNCE*4),MEB6(NNUM1-1,NF,1000))
      ALLOCATE (MEB0(NNCE))
      ALLOCATE (MEB1(NNCE))
      ALLOCATE (MEB2(NNCE))
	ALLOCATE (TB(NNCE))
	ALLOCATE (xtest(NSIDEX+2),ytest(NSIDEX+2),ITANG(1:20),ITANG1(1:20))
! ccccccccccccccccccccccccccccccccccccCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C	NNCE=JJ*(NNX-1)
! !C	WRITE(*,*)'ERR',NNCE
! !C	WRITE(*,*)'CHA-DIAN---COMPUTE'
	DO I=1,NNF 
!	WRITE(*,*)I
	DO KJU=1,MJ(I)
!C	WRITE(*,*)GKK(I,KJU),MJ(I),NNX
	JJ=GKK(I,KJU)/(NNX-1)
	JJ1=GKK(I,KJU)-JJ*(NNX-1)
	IF(JJ1==0)THEN
	JJ1=NNX-1
	JJ=JJ-1
	END IF
!C	WRITE(*,*)JJ,JJ1
	Z1=ZMIN+JJ*DZ
	Z2=Z1+DZ
	X1=XMIN+(JJ1-1)*DX
	X2=X1+DX
!C	WRITE(*,*)X1,X2,Z1,Z2,GKK(I,KJU)
	L=0
	DO 371 II=1,NODE
      XI=X(II)
      ZI=Z(II)
      IF((XI.LT.(X1-E)).OR.(XI.GT.(X2+E)))GO TO 371
      IF((ZI.LT.(Z1-E)).OR.(ZI.GT.(Z2+E)))GO TO 371
      L=L+1
      NCELL_Z(I,KJU,L)=II
  371  CONTINUE
	NBNZ(I,KJU)=L
!C	WRITE(*,*)'LLLLL',L
	END DO
	END DO
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	DO I=1,NNF     !1
	KGH=0
	DO J=1,NNX-1    !2
	X11=(J-1)*DX
	X12=J*DX
	KK=0
	DO K=1,NF       !3
	IF(X(ngr(I,K))-X11.ge.-e1.AND.X(ngr(I,K))-X12.LE.0.0)THEN
	KK=KK+1
	xtest(kk)=X(ngr(I,K))
	YTEST(KK)=Z(ngr(I,K))
!c	WRITE(*,*)'ytrr',XTEST(KK),YTEST(KK)
	END IF
	END DO          !E3
	DO ILJ=1,MJ3(J,I)  !E10
	KGH=KGH+1
      LJ=MJ1(I,KGH)
	DO KJU=1,MJ(I)    !4
	JJ=GKK(I,KJU)
	IF(JJ==LJ)THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	JN=NBNZ(I,KJU)
	KLK=0
	KKI=0
	NCE=0
	NCE1=0
	DO KKK=1,JN         !DO5
	JP=NCELL_Z(I,KJU,KKK)
	DO KL=1,KK           !6
		IF(ABS(XTEST(KL)-X(JP)).LE.E1)THEN
	     IF(YTEST(KL).GE.Z(JP)-E1)THEN
	        KLK=KLK+1
			NCE(KLK)=JP
	     END IF
	    IF(Z(JP)+E1.GE.YTEST(KL))THEN
	          KKI=KKI+1
	          NCE1(KKI)=JP
	     END IF
	END IF
	END DO              !E6
	END DO              !E5
!C	WRITE(*,*)'KLK',KLK,KKI
	DO KKK=1,KLK
	NCELL_Z(I,KJU,KKK)=0
	NCELL_Z(I,KJU,KKK)=NCE(KKK)
	END DO
	NBNZ(I,KJU)=KLK

	IF(I.LT.NNF)THEN
	DO KKK=1,MJ(I+1)
	IF(JJ==GKK(I+1,KKK))THEN
	KB=KKK
	GO TO 13
! C	WRITE(*,*)KB,LJ,JJ
! C	PAUSE
	END IF
	END DO

13	DO KKK=1,KKI
	NCELL_Z(I+1,KB,KKK)=0
	NCELL_Z(I+1,KB,KKK)=NCE1(KKK)
	END DO
	NBNZ(I+1,KB)=KKI
	END IF
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	END IF
	END DO            !E4
	GO  TO 14
14	END DO             !E10
	END DO             !E2
	END DO             !E1
	DEALLOCATE(MJ1,MJ2,MJ3,XTEST,YTEST)

	 OPEN(13,FILE='output\2dpath-d2.dat')
       open(333,FILE='output\2dpath-d2.bln')
       OPEN(16,FILE='output\2dtime-d2.dat')
       OPEN(19,FILE='temp\intersection.txt')
	DO JOI=1,NNFF1       

	NN222=LU2(JOI)
! C	WRITE(*,*)NN222
! C	PAUSE
	DO I=1,NNUM1
	ITANG(I)=LU1(JOI,I) 
	ITANG1(I)=LU3(JOI,I)
! C		WRITE(*,*)'ITANG-',ITANG(I),NNUM1
	END DO

	 IS=0 
 211   IS=IS+1 
	KJH=0
	NNF1=IS
      IB1=0
	IKJ=NDC(IS)        
	ISII=0
	LIH=1

!c---- loop over for sources-------------------------------c                
	TT=1.E+10   
      IS0=NGS(IS)
	IF(ITANG(1).GE.IKJ)THEN
	IEER=1     !���в�
	ELSE
	IEER=2     !���в�
	END IF
      TT(IS0)=0.0
!c-----find the source blocks------------------------------c
      X0=X(IS0)
      Z0=Z(IS0)
!C	write(*,*)'source',x0,z0
!cccccccccccccccccccccccccccccccccccccccccccc   
!C      write(*,*)'Ѱ���ڵ�'    
872   M1=0
      IF(IKJ==0)THEN
	GO TO 9000
	END IF
!c	VT=V(IKJ)

      DO 54 I=1,NNZ-1
      Z1=ZMIN+FLOAT(I-1)*DZ
      Z2=Z1+DZ
      IF((Z0.LT.(Z1-E)).OR.(Z0.GT.(Z2+E))) GO TO 54
      DO 50 K=1,NNX-1
      X1=XMIN+FLOAT(K-1)*DX
      X2=X1+DX
      IF((X0.LT.(X1-E)).OR.(X0.GT.(X2+E))) GO TO 50     
      IB=(I-1)*(NNX-1)+K       !�ڵ����ڵĸ��Ӻ� 
!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      NU=0   
      DO KN=1,MJ(IKJ)
        IF(IB==GKK(IKJ,KN))THEN
			NU=1
          NUU=KN
          GO TO  11
        END IF
      END DO
	if(nu==0)go to 50
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
11    M1=M1+1
      MEB2(M1)=NUU
      MEB0(M1)=0
  50  CONTINUE
  54  CONTINUE
      IF(M1.EQ.0)THEN
      WRITE(*,*)'     Something wrong with the source position !'
      STOP
      ENDIF
!c-----loop over MEB2(M1)-----------------------------------c
      IBLOCK=0
      ZM1= 1.E+10
      ZM2=-1.E+10
      XM1= 1.E+10
      XM2=-1.E+10
!c-----determine the mini.times of blocks-------------------c   
  46  DO 70 I=1,M1
       DO 53 I1=1,M1
        IB=MEB2(I1)
	  NETP=NBNZ(IKJ,IB)
        TBMIN=1.E+10
        DO J1=1,M1
         IF(MEB0(J1).EQ.IB)GO TO 520
        ENDDO
        DO 51 J1=1,NETP
         II=NCELL_Z(IKJ,IB,J1)
         TD=TT(II)
         IF(TD.GE.TBMIN)GO TO 51
         TBMIN=TD
  51    CONTINUE				   
 520    TB(I1)=TBMIN
  53   CONTINUE
!C      WRITE(*,*)'Ѱ�ҵ�Ԫ������С��ʱ'
!c-----look for the mini.time block------------------------c
       II=I        
       TTMIN=1.E+10
       DO 540 I1=1,M1
        TD=TB(I1)
        IF(TD.GE.TTMIN)GO TO 540
        TTMIN=TD
        II=I1
 540   CONTINUE
       TB(II)=1.E+10
       MEB0(I)=MEB2(II)
!c-----start at the minimum-time block--------------------c	
       IBLOCK=IBLOCK+1
       NOB=MEB2(II)
       MEB1(IBLOCK)=NOB
       NETP=NBNZ(IKJ,NOB)
	 NOB1=GKK(IKJ,NOB)
	IF(ITANG1(LIH)==1)THEN
	 CALL VBLOCK_2D(NOB1,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,V,X0,Z0,BV,NB)
	ELSE
	 CALL VBLOCK_2D(NOB1,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,VS,X0,Z0,BV,NB)
	END IF
       DO 58 J=1,NETP
        IP1=NCELL_Z(IKJ,NOB,J)
        XP1=(X(IP1)-X0)/DX
        ZP1=(Z(IP1)-Z0)/DZ
        V1=VFU(XP1,ZP1,BV)
     
        IF(X(IP1).LT.XM1)XM1=X(IP1)
        IF(X(IP1).GT.XM2)XM2=X(IP1)
        IF(Z(IP1).LT.ZM1)ZM1=Z(IP1)
        IF(Z(IP1).GT.ZM2)ZM2=Z(IP1)        
        DO 56 K=J+1,NETP
         IP2=NCELL_Z(IKJ,NOB,K) 
         XP2=(X(IP2)-X0)/DX
         ZP2=(Z(IP2)-Z0)/DZ
         V2=VFU(XP2,ZP2,BV)
         XX=X(IP2)-X(IP1)
         ZZ=Z(IP2)-Z(IP1)         
         DD=SQRT(XX*XX+ZZ*ZZ)
        V12=0.5*(V1+V2)
!C	   V12=VT
         T=DD/V12
	   T12=TT(IP1)+T
         IF(T12.GE.TT(IP2))GO TO 55
         TT(IP2)=T12   !��Ӧ����ʱ
         ICD(IP2)=IP1  !����·��
  55     T21=TT(IP2)+T
         IF(T21.GE.TT(IP1))GO TO 56
         TT(IP1)=T21
         ICD(IP1)=IP2 
  56    CONTINUE
  58   CONTINUE
  70  CONTINUE
!C     WRITE(*,*)'��ʱѰ�����'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IBLOCK.GE.MJ(IKJ))GO TO 87         
!c-----look for new blocks by expanding the range----------c
      XX1=XM1-DX
      IF(XX1.LT.XMIN)THEN
      XM1=XMIN
      ELSE
      XM1=XX1
      ENDIF
      XX2=XM2+DX
      IF(XX2.GT.XMAX)THEN
      XM2=XMAX
      ELSE
      XM2=XX2
      ENDIF

      ZZ1=ZM1-DZ
      IF(ZZ1.LT.ZMIN)THEN
      ZM1=ZMIN
      ELSE
      ZM1=ZZ1
      ENDIF
      ZZ2=ZM2+DZ
      IF(ZZ2.GT.ZMAX)THEN
      ZM2=ZMAX
      ELSE
      ZM2=ZZ2
      ENDIF

!C	write(*,*)'NEW������'

      M1=0
      DO 88 I=1,NNZ-1
       Z1=ZMIN+FLOAT(I-1)*DZ
       Z2=Z1+DZ
       ZZ0=0.5*(Z1+Z2)
       IF((ZZ0.LT.ZM1).OR.(ZZ0.GT.ZM2)) GO TO 88

        DO 84 K=1,NNX-1
         X1=XMIN+FLOAT(K-1)*DX
         X2=X1+DX
         XX0=0.5*(X1+X2)
         IF((XX0.LT.XM1).OR.(XX0.GT.XM2)) GO TO 84
         IB=(I-1)*(NNX-1)+K
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
		NU=0        
	    DO KN=1,MJ(IKJ)
	      IF(IB==GKK(IKJ,KN))THEN
			NU=1
	        IB1=KN
	        GO TO  12
	      END IF
	   END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF(NU==0)GO TO 84
12        DO L=1,IBLOCK
          IDONE=MEB1(L)
          IF(IB1.EQ.IDONE)GO TO 84
         ENDDO
         M1=M1+1
         MEB2(M1)=IB1
         MEB0(M1)=0
  84    CONTINUE
  88  CONTINUE

      IF(M1.EQ.0)THEN
      WRITE(*,*)'Something wrong with the source position-subsource'
      stop
	END IF

      GO TO 46
!c---- pick the traveltime for each receiver-----------------------c
87	tt_min=1000000.0
      do it=1,NF
	IF(IKJ.GT.1.AND.IEER==2)THEN !��Դ�ڵ�����Ϊ���в�
	iy=ikj-1          
	it1=ngr(IKJ-1,it) !��������ϵĵ�
	ELSE !���в�
	iy=ikj
	it1=ngr(ikj,it)
	END IF

	if(tt(it1).lt.tt_min)then
	tt_min=tt(it1)               
	it2=it        
	it3=it1        
	end if
	end do


	x0=x(it3)
	z0=z(it3)    !�����ϵ������

	IF(IEER==1)THEN  !���в�
	IF(IKJ.EQ.ITANG(LIH))THEN 
	GO TO 871
	ELSE
!cccccccccccccccccccccccccccccccccccccc
	IKJ=IKJ+1
	GO TO 872
	ENDIF

	ELSE      !���в�

	IF(IY.EQ.ITANG(LIH))THEN
	GO TO 871
	ELSE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IKJ=IKJ-1
	IF(IKJ==0)THEN
	GO TO 871
	END IF
	GO TO 872
	ENDIF 
	END IF
! cccccccccccccccccccccccccccccccccccccccccccccc
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
 871   ID=0   
	IF(LIH==1)THEN           !90001  !LIH=1׷����Դ�����һ������
	DO I=1,NF
       IR=NGR(IY,I)  
	 MEB4(i,1)=ir !MEB4�洢���ߴ����Ľڵ�
       LL=1
  89   N2=ICD(IR)
       LL=LL+1
	 meb4(i,ll)=n2
        IF(IR.EQ.IS0) THEN
	   MEB4(I,LL)=IR
        GO TO 90
        ELSE
        ENDIF

        IF(N2.EQ.IS0) GO TO 90
        IF(LL.GE.NODE) THEN
         WRITE(*,*)'YOUR IPRE(I) OR NGS(I) EXISTS ERROR!'
         STOP
         ENDIF

       IR=N2
      GO TO 89

  90    II1=LL-1
        ICD1(1,I)=LL
      END DO
	END IF
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(LIH==NN222)THEN        !׷���������һ�μ�����������Ƕ�
	ID=0
	tr11=0.0
	meb5=0.0
	DO I=NS0+1,NS
       IR=NGS(I)
	 MEB5(i-NS0,1)=ir
       LL=1
  891   N2=ICD(IR)
       LL=LL+1
	 meb5(i-NS0,ll)=n2
	DO IKT=1,NF
	IS01=NGR(itang(LIH-1),IKT)
       IF(IR.EQ.IS01) THEN
	  ISII(I-NS0)=IKT
	  MEB5(i-NS0,LL)=IR
        GO TO 901	
       ELSE
       ENDIF
	END DO

!C       IF(N2.EQ.IS0) GO TO 901
        IF(LL.GE.NODE) THEN
         WRITE(*,*)'YOUR IPRE(I) OR NGS(I) EXISTS ERROR!'
         STOP
         ENDIF
        IR=N2
        GO TO 891
  901    II1=LL-1
        ICD1(NN222,i-NS0)=LL
 100  CONTINUE

      END DO
	END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	IF(LIH.LT.NN222.AND.LIH.GT.1)THEN !׷���м����Щ���߶�
	LL=0
	KJH=KJH+1
	ID=0
	DO I=1,NF
       IR=NGR(IY,I)
	 MEB6(LIH,i,1)=ir
       LL=1
  894   N2=ICD(IR)
       LL=LL+1
	 meb6(LIH,i,ll)=n2
	DO IKT=1,NF
	  IS01=NGR(ITANG(LIH-1),IKT)
	  ISII2(LIH,I)=I
       IF(N2.EQ.IS01) THEN
	  ISII1(LIH,I)=IKT
	  MEB6(LIH,I,LL)=N2
       GO TO 904
       ELSE
       ENDIF
	END DO
!C      IF(N2.EQ.IS0) GO TO 904
        IF(LL.GE.NODE) THEN
         WRITE(*,*)'YOUR IPRE(I) OR NGS(I) EXISTS ERROR!'
         STOP
         ENDIF
        IR=N2
        GO TO 894
  904    II1=LL-1
        ICD1(KJH+1,I)=LL
      END DO
	END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	IF(LIH==NN222)THEN
	      M1=0
	      DO 1000 IRE=1,NR0
	             K0=0
	            K0=ISII(IRE)
                  M23(NN222)=K0
                  M111(NN222)=ICD1(nn222,ire)
	            DO I4=NN222-1,2,-1	     
                        KY=ISII2(I4,K0)
                        M111(I4)=ICD1(I4,KY)
                        K0=ISII1(I4,Ky)
                        M23(I4)=K0
	            END DO 
                  M111(1)=ICD1(1,K0)
                  M23(1)=K0
                  WRITE(13,99)IRE,IS 
                  M2=0
                  DO i4=1,NN222
                        M2=M2+M111(I4)
	            END DO
       !"xxxxxxxxxxxxxxxxx"
                  m2=m2-nn222+1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  WRITE(13,94)M2-1
                  WRITE(333,*)M2-1,1
                  IF(inv > 0) noderay(JOI+1,IS,IRE) = M2-1
                  WRITE(19,*) NN222-1,(LU1(JOI,I1),I1=1,NN222-1) 
                  DO I1=M111(1),1,-1
                        IRAY=MEB4(K0,I1)
                        WRITE(13,93)X(IRAY),Z(IRAY)
                        WRITE(333,93)X(IRAY),Z(IRAY)
                  ENDDO
        !-----------write the intersection of ray path and interface----c
                  IR2=MEB4(K0,1)
                  XR0=X(IR2)
	            ZR0=Z(IR2)
                  DO I1=2,M111(1)
	                  IR1=MEB4(K0,I1)
		            IF(ABS(X(IR1)-XR0).GT.1.0D-4.OR.ABS(Z(IR1)-ZR0).GT.1.0D-4)THEN			    
			            write(19,93) X(IR1),Z(IR1) ! 反射点的上一个点的坐标
	                        GO TO 31
	                  ENDIF
	            ENDDO
31		      write(19,93) X(IR2),Z(IR2) ! 反射点坐标 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	            DO  I1=2,NN222-1	    
	                  DO I2=M111(I1)-1,1,-1
	                        IRAY=meb6(I1,M23(I1+1),I2)
                              WRITE(13,93)X(IRAY),Z(IRAY)
                              WRITE(333,93)X(IRAY),Z(IRAY)
	                  END DO
                        DO I2=M111(I1)-1,2,-1
		                  IR1=MEB6(I1,M23(I1+1),I2)
		                  IF(ABS(X(IR1)-XR0).GT.1.0D-4.OR.ABS(Z(IR1)-ZR0).GT.1.0D-4)THEN			    
		                        write(19,93) X(IR1),Z(IR1) ! 反射点的下一个点的坐标
	                              GO TO 35
	                        ENDIF
	                  ENDDO
35                      IR3 = MEB6(I1,M23(I1+1),1)
                        XR0 = X(IR3)
		            ZR0 = Z(IR3)
                        DO I2=2,M111(I1)-1
			            IR2=MEB6(I1,M23(I1+1),I2)
		                  IF(ABS(X(IR2)-XR0).GT.1.0D-4.OR.ABS(Z(IR2)-ZR0).GT.1.0D-4)THEN			    
			                  write(19,93) X(IR2),Z(IR2) ! 反射点的上一个点的坐标
	                              GO TO 39
	                        ENDIF
	                  ENDDO
39                      write(19,93) X(IR3),Z(IR3)
	            ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  DO I2=M111(nn222)-2,1,-1
                        IRAY=MEB5(IRE,I2)
                        WRITE(13,93)X(IRAY),Z(IRAY)
                        WRITE(333,93)X(IRAY),Z(IRAY)
                  ENDDO
                  DO I2=M111(NN222)-2,1,-1
		            IR1=MEB5(IRE,I2)
		            IF(ABS(X(IR1)-XR0).GT.1.0D-4.OR.ABS(Z(IR1)-ZR0).GT.1.0D-4)THEN			    
			            write(19,93) X(IR1),Z(IR1) ! 反射点的下一个点的坐标
	                        GO TO 1000
	                  END IF
                  ENDDO
1000	      ENDDO

      WRITE(16,"(10(g0,','))") TT(NGS(1+NS0:NS0+NR0))
      OPEN(101,file="traveltime_reflect.dat")
       write(101,"(3(g0,','))") (x(i),z(i),TT(i),i=1,nnx*nnz)
      CLOSE(101)
      
ENDIF
!!C      write(*,*)'����ʱ���'

	DO IU=1,NF
	ITR=NGR(iy,IU)
	TT1(IU)=TT(ITR)
	END DO

	TT=1E10

	DO IU=1,NF
	ITR=NGR(iy,IU)
	TT(ITR)=TT1(IU)
	END DO
!C      WRITE(*,*)'׼�����ڼ���'
	IKJ=IKJ
	IF(IKJ==0)THEN
	 GO TO 9000
	END IF

	IF(LIH.LE.NN222)THEN
	IF(IEER==1)THEN
	LIH=LIH+1
	IF(ITANG(LIH).GT.ITANG(LIH-1))THEN
	IEER=1
	ELSE
	IEER=2
	END IF
	GO TO 872
	ELSE
	LIH=LIH+1
	IF(ITANG(LIH).GT.ITANG(LIH-1))THEN
	IEER=1
	ELSE
	IEER=2
	END IF
	go to 872
	END IF
	ELSE
	GO TO 9000
	END IF
	CLOSE(67)

9000  IF(IS.LT.NS0)GO TO 211
      END DO

      CLOSE(13)
      CLOSE(16)
      close(19)
      close(333)
	DO I=1,NNF
	WRITE(OPP4,*)I
	OPP4=TRIM(ADJUSTL(OPP4))//'T.TXT'
!c	WRITE(*,*)OPP4
	OPEN(87,FILE=OPP4)
	CLOSE(87,STATUS='DELETE')
	END DO
      close(133)
  93  FORMAT(1X,F10.4,1X,F10.4)
  94  FORMAT(1X,I10)
  99  FORMAT(1X,I10,1X,I10)
 101  FORMAT(1X,'SOURCE NUMBER=',I10)
 
      DEALLOCATE(X,Z,TT,ICD,MEB0,MEB1,MEB2,NGS,NGR,TB,ICD0,MEB4,MEB5,ICD1) 
	OPEN(1006,FILE='TIME.TXT')
	CALL CPU_TIME(TIMEEN)
	WRITE(1006,*)TIMEEN-TIMEST
	CLOSE(1006)
      RETURN
      END Subroutine RAY2D_REF
! C-------------------------------------------------------------------------------
!C-------------------------------------------------------------------------------
	
      SUBROUTINE TIAOGEZI(XMIN,XMAX,ZMIN,ZMAX,LL,NFT,X,Z,NGR,NF,NNX,NNZ,&
 &  FF,DX,DZ,MJ,MJ1,MJ2,MJ3)
	          
	IMPLICIT NONE
	INTEGER NNX,NNZ,NF,LL,I,J,KJ,K,JI,KH,NFT
	REAL XMIN,XMAX,ZMAX,ZMIN,XX(8),YY(8),DX,DZ,YGEMAX,YGEMIN
	REAL YGEMIN1,XK,X2,ZI,Z2,X(NFT),Z(NFT)
!C	REAL XR(LL,NF),YR(LL,NF)
	INTEGER NGR(LL,NF),FF(1:LL+1,1:(NNX-1)*(NNZ-1))!FF(1:LL+1,1:(NNX-1)*NNX)
	INTEGER MJ(1:LL+1),MJ1(1:LL,1:10*NNX),&!MJ1(1:LL,1:10*NNX)
     &    	MJ2(1:LL),MJ3(NNX-1,LL)
	MJ=0
	MJ1=0
	MJ2=0
      YGEMIN1=0

	mj=0
	KH=0
	MJ3=0
	DO I=1,NNX-1
	   XK=XMIN+FLOAT(I-1)*DX
	   X2=XK+DX
	           
       DO  kj=1,ll+1        !������
			YGEMAX=0.0
	        YGEMIN=100000.0
	   IF(KJ.LT.LL+1)THEN
	     DO K=1,NF
	          IF(X(NGR(kj,K)).GE.XK.AND.X(NGR(kj,K)).LE.X2)THEN
	            YGEMAX=MAX(Z(NGR(kj,K)),YGEMAX)    !xĳ���䷴����z�������Сֵ
	            YGEMIN=MIN(Z(NGR(kj,K)),YGEMIN)
	          END IF
	     END DO
! C		WRITE(*,*)YGEMAX,YGEMIN
! C	    PAUSE
		YY(KJ)=YGEMIN
	    XX(KJ)=YGEMAX
	ELSE
	YY(KJ)=ZMAX
	XX(KJ)=ZMAX
	END IF
! C		WRITE(*,*)YGEMIN,YGEMAX,XK,X2
! C		PAUSE		
	  if(kj==1)then
	   ygemin1=Zmin
	  ELSE                 !1���ڵ�ļ�ֵ
		ygemin1=YY(KJ-1)
	 end if

	ygemin1=ygemin1-DZ  

	IF(ygemin1.LT.ZMIN)THEN
		YGEMIN1=ZMIN
	END IF

	   
      DO J=1,NNZ-1
	   ZI=ZMIN+FLOAT(J-1)*DZ
	   Z2=ZI+DZ
	 IF(zi.lE.XX(KJ).and.zi.gE.ygemin1)THEN 
		
         
	    IF(KJ.GT.1.AND.zi.gT.ygemin1.AND.zi.lE.XX(max(1,KJ-1)))THEN 
	    MJ2(KJ-1)=MJ2(KJ-1)+1
	    MJ1(KJ-1,MJ2(KJ-1))=I+(NNX-1)*(J-1)	
		MJ3(I,KJ-1)=MJ3(I,KJ-1)+1	  
	    END IF

	do ji=1,ll+1
	   if(ji==kj)then
	      mj(ji)=mj(ji)+1  
	      ff(kj,mj(ji))=I+(NNX-1)*(J-1)
	   end if
	end do	
      END IF

	END DO
	end do
	END DO
	end SUBROUTINE TIAOGEZI
! c------------------------------------------------------------c
! c                                                            c
! c     Find the local velocity for a given number of block    c
! c                                                            c
! c     Entries: (1) NOB............specified number of cell;  c
! c              (2) XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX.........    c
! c                  range of x-,y- and z-axies;               c
! c              (3) DX,DY,DZ.....................cell size;   c
! c              (4) V(*)...........input velocity model;      c
! c                                                            c
! c     Returns: (1) X0,Y0,Z0.......coordinates of the cell;   c
! c              (2) BV(8)....velocity samples of the cell.    c 
! c                                                            c
! c------------------------------------------------------------c
      SUBROUTINE VBLOCK_2D(NOB,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,V,X0,Z0,BV,NB) 
    	IMPLICIT NONE
	INTEGER NOB,NX,NZ,IM,I,K,IE,I1,I2,I3,I4
	REAL XMIN,XMAX,ZMIN,ZMAX,DX,DZ,X0,Z0
      REAL V(*),BV(*)
	INTEGER NB(*)

!c-----prepare some parameters---------------------c      
      NX=INT((XMAX-XMIN)/DX+0.5)
      NZ=INT((ZMAX-ZMIN)/DZ+0.5)
    
      IM=0
      DO 9 I=1,NZ
        DO 5 K=1,NX
         IM=IM+1
         IE=(I-1)*(NX+1)+K 
         IF(IM.NE.NOB)GO TO 5
         Z0=ZMIN+FLOAT(I-1)*DZ
         X0=XMIN+FLOAT(K-1)*DX
         I1=IE
         I2=IE+1
         I3=I1+(NX+1)
         I4=I3+1
	   NB(1)=I1
         NB(2)=I2
         NB(3)=I3
         NB(4)=I4
         BV(1)=V(I1)
         BV(2)=V(I2)
         BV(3)=V(I3)
         BV(4)=V(I4)
    5    CONTINUE
   9  CONTINUE
      GO TO 80
  80  RETURN
      END Subroutine VBLOCK_2D
! c-------------------------------------------------------------c
! c                                                             c
! c     Velocity function in a specified cell                   c
! c                                                             c
! c     Entries:                                                c
! c              (1) (x,y,z).....the point to be calculated;    c
! c              (2) BV(8).....velocity samples of the cell.    c
! c                                                             c
! c     Return:  VFU....velocity at a given point (x,y,z).      c
! c                                                             c
! c-------------------------------------------------------------c
      REAL FUNCTION VFU(X,Z,BV)
	IMPLICIT NONE
	REAL X,Z,VN1,VN2,VN3,VN4
	REAL BV(1:4)
      VN1=(1.-X)*(1.-Z)
      VN2=X*(1.-Z)
      VN3=(1.-X)*Z
      VN4=X*Z
      VFU=VN1*BV(1)+VN2*BV(2)+VN3*BV(3)+VN4*BV(4)
     
  40  RETURN
      END  FUNCTION VFU
! c-------------------------------------------------c
! c       THIS SUBROUTINE IS A 2-D FILTERING ONE    c
! c-------------------------------------------------c
      SUBROUTINE FILT2(N,NVA,ALPH,C)
	IMPLICIT NONE
	INTEGER N,NVA,IDX,I,K,K1,II1
	REAL ALPH,C(1:N),SS,AA,S
      REAL,ALLOCATABLE::A1(:)
      ALLOCATE(A1(1:N))
       IDX = FLOAT(NVA) / 2.0
       DO 20 I=1,N
       SS = C(I)
       IF (SS.EQ.1.) GOTO 20
       K = 0
       DO 28 K1=1,NVA
       II1 = I + (K1 - (IDX + 1))
       IF((II1.LE.0).OR.(II1.GT.N))GO TO 28
       AA = C(II1)
       IF (AA.EQ.1.) GOTO 28
       K = K + 1
       A1(K) = AA
  18   CONTINUE
  28   CONTINUE
       CALL FILT1(K, A1, ALPH, S)
       C(I) = S
  20   CONTINUE
       DEALLOCATE(A1)
       RETURN
       END
! c------------------------------------------------c
! c      THIS SUBROUTINE IS A 1-D FILTERING ONE    c
! c------------------------------------------------c
      SUBROUTINE FILT1(N, A, ALPH, S)
	IMPLICIT NONE
	INTEGER N,N1,N2,NN,I1,I2,K,I
	REAL ALPH,S,A1,XMIN
      REAL A(1:N)
       DO 4 K=1,N
       XMIN = 500000.0
       DO 3 I=K,N
       IF (A(I).GE.XMIN) GOTO 3
       XMIN = A(I)
       A1 = A(K)
       A(K) = XMIN
       A(I) = A1
  3    CONTINUE
  4    CONTINUE
  5    N1=N/2
       N2 = N1 * 2
       IF (N2.EQ.N) THEN
       NN = N - 1
       ELSE
       NN = N
       END IF
       N1 = ALPH * NN
       I1 = N1 + 1
       I2 = N - N1
       A1 = 0
       DO 7 I=I1,I2
  7    A1=A1+A(I)/(N-2*N1)
       S = A1
       RETURN
       END
! c------------------------------------------------------------------c
! c                                                                  c
! c     This subroutine performs 2D ray tracing for any              c  
! c     volocity model. The network method was used  to              c
! c     compute the minimum traveltimes and raypaths of              c
! c     a grided model.                                              c
! c                                                                  c
! c     Entries: (1) XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX......range;       c
! c              (2) DX,DY,DZ...............cell size of grid;       c
! c              (3) NSIDE..............density of ray points;       c
! c              (4) XS(NS),YS(NS),ZS(NR).....source location;       c
! c              (5) XR(NR),YR(NR),ZR(NR)..geophone locations;       c
! c              (6) MSR......matrix for source-geophone pair;       c
! c              (6) V(*).................volocity parameters.       c
! c                                                                  c
! c     Results: (1) 3D traveltime  file:  3Dtravel-time.txt;        c
! c              (2) 3D  raypath coordinates: 3Dray-path.txt;        c
! c                                                                  c
! c------------------------------------------------------------------c 
      SUBROUTINE RAY2D(XMIN,XMAX,ZMIN,ZMAX,DX,DZ,NSIDE,NS,XS,ZS,NR,XR,ZR,V,MSR)
      use globalp,only : inv,noderay
      IMPLICIT NONE
	INTEGER NSIDE,NS,NR,NNX,NNZ,NP1,NP2,NPP,NP,NODE,NBLOCK, &
     &        NOB,NODE1,NETP,IBLOCK,IDONE,IB,ID,IBP,IS,IS0,IS1,IR1, &
     &        IP1,IP2,IRR,IR,I,K,J,L,II,JJ,KK,K1,I1,J1,LL,N2, &
     &        II1,IRAY
      REAL XMIN,XMAX,ZMIN,ZMAX,DX,DZ,DMX,DMZ,XK,ZI,XJ,X0,Z0,ZJJ,&
     &     XKK,E,A1,A3,D,Z1,Z2,X1,X2,XI,M1,ZM1,ZM2,XM1,XM2,TD, &
     &     TBMIN,TTMIN,XP1,ZP1,V1,VFU,XP2,ZP2,V2,XX,ZZ,DD,V12,T, &
     &     T12,T21,XX1,XX2,ZZ1,ZZ2,ZZ0,XX0
      INTEGER NB(1:4),MSR(1:NS,1:NR)
	REAL BV(1:4),XS(*),ZS(*),XR(*),ZR(*),V(*)
	INTEGER,ALLOCATABLE::NCELL(:,:),NGS(:),NGR(:),NBN(:),ICD(:),MEB0(:),MEB1(:),MEB2(:)
      REAL,ALLOCATABLE:: X(:),Z(:),TT(:),TB(:),TR(:)

      
      OPEN(11,FILE='output/2dpath-d.txt')
      OPEN(12,FILE='output/2dtime-d.txt')
      open(13,file="output/2dpath.bln")


      NNX=INT((XMAX-XMIN)/DX+0.5)+1
      NNZ=INT((ZMAX-ZMIN)/DZ+0.5)+1
      
      NP1=((NNX-1)*(NSIDE+1)+1)*NNZ
      NP2=(NNZ-1)*NSIDE*NNX
      NPP=NP1+NP2+1000
      NP=(NNZ-1)*(NNX-1)

      DMX=DX/FLOAT(NSIDE+1)
      DMZ=DZ/FLOAT(NSIDE+1)
      
      ALLOCATE(X(1:NPP),Z(1:NPP))

      NODE=0
!c-----Creat primary node------------------   
      DO 222 I=1,NNZ
       ZI=ZMIN+FLOAT(I-1)*DZ
       IF(ZI.GT.ZMAX) GOTO 222
       DO 188 K=1,NNX
        XK=XMIN+FLOAT(K-1)*DX
        IF(XK.GT.XMAX) GOTO 188
      NODE=NODE+1
      X(NODE)=XK
      Z(NODE)=ZI
 188  CONTINUE
 222  CONTINUE
!c-----Insert the secondary nodes-------------         
      DO 22 I=1,NNZ
       ZI=ZMIN+FLOAT(I-1)*DZ
       DO 20 J=1,NNX
        XJ=XMIN+FLOAT(J-1)*DX
        DO 8 JJ=1,NSIDE+1
         IF(JJ.EQ.1) THEN
          DO 6 KK=2,NSIDE+1
           ZJJ=ZI
           XKK=XJ+FLOAT(KK-1)*DMX
           IF(XKK.GT.XMAX)GO TO 6
           NODE=NODE+1
           X(NODE)=XKK
           Z(NODE)=ZJJ
   6     CONTINUE
         ELSE
          XKK=XJ
          ZJJ=ZI+FLOAT(JJ-1)*DMZ          
          IF(ZJJ.GT.ZMAX)GO TO 8
          NODE=NODE+1
          X(NODE)=XKK
          Z(NODE)=ZJJ
         ENDIF 
   8    CONTINUE 
  20  CONTINUE  
  22  CONTINUE
! c---------------------------------------------------c
! c     loop over for the sources                     c
! c---------------------------------------------------c
      NBLOCK=(NNZ-1)*(NNX-1)
      IBP=4*NSIDE+4
      ALLOCATE(NGS(1:NS))
      ALLOCATE(NGR(1:NR))

      ALLOCATE(NCELL(1:NBLOCK,1:(IBP+100)))
      ALLOCATE(NBN(1:NBLOCK))
      ALLOCATE(ICD(1:NPP))
      ALLOCATE(MEB0(1:NBLOCK))
      ALLOCATE(MEB1(1:NBLOCK))
      ALLOCATE(MEB2(1:NBLOCK))
      ALLOCATE(TB(1:NBLOCK))
      ALLOCATE(TT(1:NPP))
      ALLOCATE(TR(1:NR))   
         
       IS=0
 211   IS=IS+1             
!  !      WRITE(*,*)'SOURCE NO =',IS,'LOCATION=',XS(IS),ZS(IS)
! c----------------------------------------------------c
! c     (2) find the grid No. for source & geophones   c
! c----------------------------------------------------c

      E=MIN(DMX,DMZ)/100.0
      DO 30 I=1,NS
      X0=XS(I)
      Z0=ZS(I)
      NGS(I)=0
      DO 25 J=1,NODE
      A1=(X0-X(J))
      A3=(Z0-Z(J))
      D=SQRT(A1*A1+A3*A3)
      IF(D.GT.E)GOTO 25
      NGS(I)=J
      GO TO 30
  25  CONTINUE
  30  CONTINUE
      K=0
      DO 32 I=1,NS
      IS1=NGS(I)
      IF(IS1.NE.0)GOTO 32
      K=K+1
      K1=NODE+K
      NGS(I)=K1
      X(K1)=XS(I)
      Z(K1)=ZS(I)
  32  CONTINUE
      NODE1=NODE+K
       
      ngr=0
      DO 35 I=1,NR
      X0=XR(I)
      Z0=ZR(I)
      DO 27 J=1,NODE
      A1=(X0-X(J))
      A3=(Z0-Z(J))
      D=SQRT(A1*A1+A3*A3)
      IF(D.GT.E)GO TO 27
      NGR(I)=J
      GO TO 35
  27  CONTINUE
  35  CONTINUE
      K=0
      DO 33 I=1,NR
      IR1=NGR(I)
      IF(IR1.NE.0)GO TO 33
      K=K+1
      K1=NODE1+K
      NGR(I)=K1
      X(K1)=XR(I)
      Z(K1)=ZR(I)
  33  CONTINUE
      NODE=NODE1+K
! c---------------------------------------------------------c
! c     (3)find the node-number for each block: NCELL(*,*)  c
! c---------------------------------------------------------c
      IBLOCK=0
      DO 38 I=1,NNZ-1
      Z1=ZMIN+FLOAT(I-1)*DZ
      Z2=Z1+DZ

      DO 34 K=1,NNX-1
      X1=XMIN+FLOAT(K-1)*DX
      X2=X1+DX
      IBLOCK=IBLOCK+1
      
      L=0
      DO 37 II=1,NODE
      XI=X(II)
      ZI=Z(II)
      IF((XI.LT.(X1-E)).OR.(XI.GT.(X2+E)))GO TO 37
      IF((ZI.LT.(Z1-E)).OR.(ZI.GT.(Z2+E)))GO TO 37
      L=L+1
      NCELL(IBLOCK,L)=II

      NBN(IBLOCK)=L
  37  CONTINUE
  34  CONTINUE
  38  CONTINUE
! c---- loop over for sources-------------------------------c                
! c-----initialization--------------------------------------c      
      DO I=1,NODE
      TT(I)=1.E+10
      ENDDO

      IS0=NGS(IS)
      TT(IS0)=0.0
!c-----find the source blocks------------------------------c
      X0=X(IS0)
      Z0=Z(IS0)
            
      M1=0
      DO 54 I=1,NNZ-1
      Z1=ZMIN+FLOAT(I-1)*DZ
      Z2=Z1+DZ
      IF((Z0.LT.(Z1-E)).OR.(Z0.GT.(Z2+E))) GO TO 54

      DO 50 K=1,NNX-1
      X1=XMIN+FLOAT(K-1)*DX
      X2=X1+DX
      IF((X0.LT.(X1-E)).OR.(X0.GT.(X2+E))) GO TO 50
      M1=M1+1
      IB=(I-1)*(NNX-1)+K
      MEB2(M1)=IB
      MEB0(M1)=0
  50  CONTINUE
  54  CONTINUE

      IF(M1.EQ.0)THEN
      WRITE(*,*)'     Something wrong with the source position !'
      STOP
      ENDIF
      
!c-----loop over MEB2(M1)-----------------------------------c
      IBLOCK=0
      ZM1= 1.E+10
      ZM2=-1.E+10
      XM1= 1.E+10
      XM2=-1.E+10
!c-----determine the mini.times of blocks-------------------c  
  46  DO 70 I=1,M1
       DO 53 I1=1,M1
        IB=MEB2(I1)
        NETP=NBN(IB)
        TBMIN=1.E+10
        DO J1=1,M1
         IF(MEB0(J1).EQ.IB)GO TO 520
        ENDDO

        DO 51 J1=1,NETP
         II=NCELL(IB,J1)
         TD=TT(II)
         IF(TD.GE.TBMIN)GO TO 51
         TBMIN=TD
  51    CONTINUE	
 520    TB(I1)=TBMIN
  53   CONTINUE
!c-----look for the mini.time block------------------------c
       II=I        
       TTMIN=1.E+10
       DO 540 I1=1,M1
        TD=TB(I1)
        IF(TD.GE.TTMIN)GO TO 540
        TTMIN=TD
        II=I1
 540   CONTINUE
       TB(II)=1.E+10
       MEB0(I)=MEB2(II)
!c-----start at the minimum-time block--------------------c	
       IBLOCK=IBLOCK+1
       NOB=MEB2(II)
       MEB1(IBLOCK)=NOB
       NETP=NBN(NOB)
       CALL VBLOCK_2D(NOB,XMIN,XMAX,ZMIN,ZMAX,DX,DZ,V,X0,Z0,BV,NB)   
           
       DO 58 J=1,NETP
        IP1=NCELL(NOB,J)
        XP1=(X(IP1)-X0)/DX
        ZP1=(Z(IP1)-Z0)/DZ
        V1=VFU(XP1,ZP1,BV)
     
        IF(X(IP1).LT.XM1)XM1=X(IP1)
        IF(X(IP1).GT.XM2)XM2=X(IP1)
        IF(Z(IP1).LT.ZM1)ZM1=Z(IP1)
        IF(Z(IP1).GT.ZM2)ZM2=Z(IP1)        
        DO 56 K=J+1,NETP
         IP2=NCELL(NOB,K) 
         XP2=(X(IP2)-X0)/DX
         ZP2=(Z(IP2)-Z0)/DZ
         V2=VFU(XP2,ZP2,BV)
         XX=X(IP2)-X(IP1)
         ZZ=Z(IP2)-Z(IP1)         
         DD=SQRT(XX*XX+ZZ*ZZ)
         V12=0.5*(V1+V2)
         T=DD/V12
         T12=TT(IP1)+T
         IF(T12.GE.TT(IP2))GO TO 55
         TT(IP2)=T12
         ICD(IP2)=IP1
  55     T21=TT(IP2)+T
         IF(T21.GE.TT(IP1))GO TO 56
         TT(IP1)=T21
         ICD(IP1)=IP2 
  56    CONTINUE
  58   CONTINUE
  70  CONTINUE
      IF(IBLOCK.EQ.NBLOCK)GO TO 87         
!c-----look for new blocks by expanding the range----------------c

      XX1=XM1-DX
      IF(XX1.LT.XMIN)THEN
      XM1=XMIN
      ELSE
      XM1=XX1
      ENDIF
      XX2=XM2+DX
      IF(XX2.GT.XMAX)THEN
      XM2=XMAX
      ELSE
      XM2=XX2
      ENDIF

      ZZ1=ZM1-DZ
      IF(ZZ1.LT.ZMIN)THEN
      ZM1=ZMIN
      ELSE
      ZM1=ZZ1
      ENDIF
      ZZ2=ZM2+DZ
      IF(ZZ2.GT.ZMAX)THEN
      ZM2=ZMAX
      ELSE
      ZM2=ZZ2
      ENDIF

      M1=0
      DO 88 I=1,NNZ-1
       Z1=ZMIN+FLOAT(I-1)*DZ
       Z2=Z1+DZ
       ZZ0=0.5*(Z1+Z2)
       IF((ZZ0.LT.ZM1).OR.(ZZ0.GT.ZM2)) GO TO 88

        DO 84 K=1,NNX-1
         X1=XMIN+FLOAT(K-1)*DX
         X2=X1+DX
         XX0=0.5*(X1+X2)
         IF((XX0.LT.XM1).OR.(XX0.GT.XM2)) GO TO 84
         IB=(I-1)*(NNX-1)+K
         DO L=1,IBLOCK
          IDONE=MEB1(L)
          IF(IB.EQ.IDONE)GO TO 84
         ENDDO
         M1=M1+1
         MEB2(M1)=IB
         MEB0(M1)=0
  84    CONTINUE
  88  CONTINUE
      GO TO 46
! !
  87   ID=0
       DO 100 I=1,NR
        IRR=MSR(IS,I)
        IF(IRR.EQ.0) GOTO 100 
!	WRITE(*,*) NGR(I),NGR(6)
        IR=NGR(I)
        WRITE(11,99)IS,I
        ID=ID+1
        TR(ID)=TT(IR)
        MEB1(1)=IR
        LL=1
  89    N2=ICD(IR)
        LL=LL+1
        MEB1(LL)=N2
        IF(IR.EQ.IS0) THEN
         MEB1(LL)=IR
         GO TO 90
        ELSE
        ENDIF
       IF(N2.EQ.IS0) GO TO 90
        IF(LL.GE.NODE) THEN
         WRITE(*,*)'YOUR IPRE(I) OR NGS(I) EXISTS ERROR!'
         STOP
         ENDIF
        IR=N2
        GO TO 89
  90    II1=LL-1
       WRITE(11,94)LL
       WRITE(13,*) LL ,1
       DO 95 I1=LL,1,-1
        IRAY=MEB1(I1)
        WRITE(11,93)X(IRAY),Z(IRAY)
        WRITE(13,*) X(IRAY),Z(IRAY)
  95   CONTINUE
 100  CONTINUE
      WRITE(12,92)(TR(I1),I1=1,ID)
      IF(IS.LT.NS)GO TO 211
      CLOSE(11)
      CLOSE(12)
  92  FORMAT(10(F8.4,1X))
  93  FORMAT(1X,F10.4,1X,F10.4)
  94  FORMAT(1X,I10)
  99  FORMAT(1X,I10,1X,I10)
 101  FORMAT(1X,'SOURCE NUMBER=',I10)
      open(101,file='TraveltimeField.txt')
      write(101,"(3(g0,','))") (x(i),z(i),TT(i),i=1,nnx*nnz)
      close(101)
      DEALLOCATE(X,Z,TT,TR,ICD,MEB0,MEB1,MEB2,NCELL,NGS,NGR,TB,NBN) 
      RETURN
      END Subroutine RAY2D



      subroutine trans_li(nelement,ndata,nunk,irow,irtn,icol0,icol,ictn)
      implicit none
      integer,intent(in):: nelement, ndata, nunk
      integer,intent(in):: irow(nelement),irtn(ndata)
      integer,intent(out)::icol0(nelement),icol(nelement),ictn(nunk)
      integer,allocatable::iloc(:), iray(:)
      integer i, m, n

     
      ictn = 0; icol = 0; icol0 = 0

      do i = 1, nelement
	      m=irow(i)  
	      ictn(m) = ictn(m) + 1
      end do

      allocate( iloc(nunk) )
      m = 0; iloc = 1
      do i = 2, nunk
	      m = m + ictn(i-1)
	      iloc(i) = m + 1
      end do
      allocate( iray(nelement) )
      n = 0
      do i = 1, ndata
	      m = irtn(i)
	      n = n + m
	      iray(n-m+1:n) = i
      end do

    
      do i = 1, nelement
      m=irow(i)
      n = iloc(m) 
      icol(n) = i
      icol0(n) = iray(i)
      iloc(m) = iloc(m) + 1
      end do

      deallocate( iloc, iray )

      end subroutine trans_li