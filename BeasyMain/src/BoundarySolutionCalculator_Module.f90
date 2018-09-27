module BoundarySolutionCalculator_Module
    use, intrinsic :: iso_fortran_env;
    use DataStructure_Module;
    use Utility_Module;                                        !==> TO DELETE?
    implicit none;
    private;
    public :: DDDMS_BSProc;
    public :: DDDMS_BoundarySolutionCalculator;
    public :: DDDMS_SMBoundarySolutionCalculator; 
    
    
    ! DDDMS_BoundarySolutionCalculator definition    
    type :: DDDMS_BSProc
    contains
        procedure                       :: calcBSProc;
    end type DDDMS_BSProc
    
    ! Abstract DDDMS_BoundarySolutionCalculator user-defined data type definition    
    type, abstract :: DDDMS_BoundarySolutionCalculator        
    contains 
        procedure(CalculateABTRACT), deferred :: calculateSM;
        procedure(CalculateABTRACT), deferred :: calculateDM;        
    end type DDDMS_BoundarySolutionCalculator

    
    ! DDDMS_BoundarySolutionCalculator abstract user-defined data type procedure definition
    
    abstract interface
        subroutine CalculateABTRACT(self, inputData, zonesBounSol, inputParamsData)
            import;
            class(DDDMS_BoundarySolutionCalculator) :: self;
            class(DDDMS_DataContainer)              :: inputData;
            integer(INT32), allocatable             :: zonesBounSol(:);            
            class(DDDMS_InputParams)                :: inputParamsData;
        end subroutine CalculateABTRACT
    end interface

    
    !concrete DDDMS_SMBoundarySolutionCalculator user-defined data type definition for shared memory architecture    
    type, extends(DDDMS_BoundarySolutionCalculator) :: DDDMS_SMBoundarySolutionCalculator
        
    contains 
        procedure :: calculateSM => calculateOnSM;
        procedure :: calculateDM => calculateOnDM;
        
    end type DDDMS_SMBoundarySolutionCalculator

    
    ! DDDMS_SMBoundarySolutionCalculator user-defined data type constructor definition    
    interface DDDMS_SMBoundarySolutionCalculator
        module procedure NewDDDMS_SMBoundarySolutionCalculator; ! ADD CONSTRUCTOR TO DDDMS_SMBoundarySolutionCalculator GENERIC INTERFACE
    end interface DDDMS_SMBoundarySolutionCalculator

contains
    
    ! DDDMS_SMBoundarySolutionCalculator user-defined data type constructor implementation    
    type(DDDMS_SMBoundarySolutionCalculator) function NewDDDMS_SMBoundarySolutionCalculator(self)
    
        implicit none;
        
        
        ! Declarative zone        
        class(DDDMS_SMBoundarySolutionCalculator) :: self;
        
        
        ! Body of the program       
        
    end function NewDDDMS_SMBoundarySolutionCalculator

    
    ! Calculate procedure implementation    
    subroutine calculateOnSM(self, inputData, zonesBounSol, inputParamsData)
    
        use, intrinsic :: iso_fortran_env;
        use omp_lib;        
        implicit none;
    
        
        !Declarative zone        
        class(DDDMS_SMBoundarySolutionCalculator) :: self;
        class(DDDMS_DataContainer)                :: inputData;
        class(DDDMS_InputParams)                  :: inputParamsData;
        integer(INT32), allocatable               :: zonesBounSol(:);                                        
        integer(INT32)                            :: i, j, k;
        integer(INT32)                            :: iZone;
        integer(INT32)                            :: status;                
        integer(INT32)                            :: IndepInterfNodesCount;
        integer(INT32)                            :: countStep;        
        real(REAL64)                              :: ALPHA;
        real(REAL64)                              :: BETA;
        real(REAL64), allocatable                 :: UL(:);        
        real(REAL64), allocatable                 :: TempVEC_(:);          
! APB   
        integer(INT32)                            :: intnod,intfnod,Bcnd;
        integer(INT64)                            :: idim;
        logical                                   :: infnodstat=.false.        
! APB       
        
        
        ! Body of the program        
!-------------------------------
! Initialisation
!-------------------------------        
        if (inputParamsData%analysType .eq. 1) then     ! Problem dimension
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif

        
        ! Loop over zones in the model for which a boundary solution was required
        
    !$OMP PARALLEL PRIVATE(i, iZone, j, k, status, IndepInterfNodesCount, countStep, Bcnd, intfnod, ALPHA, BETA, UL, TempVEC_)
    !$OMP DO                
        do i = 1, size( zonesBounSol )
            iZone              = zonesBounSol(i);
            IndepInterfNodesCount = count( inputData%IndependenInterfaceNodesArray(:, iZone) /= 0 ); ! number of independent nodes lying on the interface for the current zone            
           
            
            ! Allocate temporary array "UL" to collect interface solution concerning the current zone            
            allocate( UL( idim*IndepInterfNodesCount ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of UL array!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                uL(:) = 0.0;
            end if  
            
            
            ! Collect interface solution for the current zone from vector U_InterfaceSolution
            countStep = 1;
            do j = 1, size( inputData%IndependenInterfaceNodesArray(:, iZone), 1)                 
                if( inputData%IndependenInterfaceNodesArray(j, iZone) /= 0 ) then                                 
                    if (idim .eq. 1) then
                        UL(countStep) = inputData%U_InterfaceSolution(j);                    
                    else
                        UL((idim*countStep - 2):(idim*countStep)) = inputData%U_InterfaceSolution( (idim*j - 2):(idim*j) );
                    endif                     
                    countStep = countStep + 1;                    
                end if                
            end do            
                    
            inputData%ZonesData(iZone)%ULI(:)=UL(:);            
            
            
            !(1)--DGEMV: Computes matrix-vector product: CB <= (CB - DL*UL)            
            ALPHA = -1.0; 
            BETA  = +1.0;
            
            call DGEMV('N',                                         & !matrices a not be transposed or conjugate transposed before multiplication.
                       size( inputData%ZonesData( iZone )%DL, 1 ),  & !specifies the number of rows of the matrix a.
                       size( inputData%ZonesData( iZone )%DL, 2 ),  & !specifies the number of columns of the matrix a.
                       ALPHA,                                       & !specifies the scalar alpha.
                       inputData%ZonesData( iZone )%DL,             & !matrix a, size (lda, n).
                       size( inputData%ZonesData( iZone )%DL, 1 ),  & !leading dimension of array a.
                       UL,                                          & !vector x.
                       1,                                           & !specifies the increment for the elements of x.
                       BETA,                                        & !specifies the scalar beta. 
                       inputData%ZonesData( iZone )%CB,             & !vector y.
                       1 );                                           !specifies the increment for the elements of y.            
                        
            ! Rewrite interface results for the nodal directions where physical continuity BC's are applied
            ! Interface traction results(CB) are used instead of displacements
            countStep=0
            infnodstat=.false.
            do j = 1, size( inputData%IndependenInterfaceNodesArray(:, iZone), 1)
                intfnod=inputData%IndependenInterfaceNodesArray(j, iZone)
                infnodstat=inputParamsData%isNodMin(inputData,intfnod);
                if( intfnod /= 0) then
                    countStep=countStep+1                    
                    if (infnodstat) then    ! check if the node considered equals the lower node ID of the node pair in interface definition file                                                                  
                            do k=1,size(inputData%ZonesData( iZone )%CoordsNodesInZone(:,1))     ! Get zonal BC type index corresponding to global node ID
                                if (intfnod.eq.int(inputData%ZonesData( iZone )%CoordsNodesInZone(k,1))) then
                                    Bcnd=k                                    
                                    exit
                                endif
                            enddo                            
                          
                            ! check if the node is assigned physical continuity BC's
                            if( inputData%ZonesData( iZone )%BCsType(Bcnd, 1) == 39 ) then   
                                if (idim .eq. 1) then
                                    inputData%ZonesData(iZone)%ULI(countStep)=inputData%ZonesData(iZone)%CB(countStep)
                                else
                                    inputData%ZonesData(iZone)%ULI(idim*countStep - 2)=inputData%ZonesData(iZone)%CB(idim*countStep - 2)
                                endif
                            endif
                            if (idim .ne. 1) then
                                if( inputData%ZonesData(iZone)%BCsType(Bcnd, 2) == 40 ) then
                                    inputData%ZonesData(iZone)%ULI(idim*countStep - 1)=inputData%ZonesData(iZone)%CB(idim*countStep - 1)
                                endif                        
                                if( inputData%ZonesData(iZone)%BCsType(Bcnd, 3) == 41 ) then
                                    inputData%ZonesData(iZone)%ULI(idim*countStep)=inputData%ZonesData(iZone)%CB(idim*countStep)         
                                endif
                            endif
                    endif
                end if
            enddo           
            
            
            !(2)--DGEMV: computes matrix-vector product: XB <= B0L*(CB - DL*UL)            
            ALPHA = +1.0; 
            BETA  = +0.0;
            
            call DGEMV('N',                                         & !matrices a not be transposed or conjugate transposed before multiplication.
                       size( inputData%ZonesData( iZone )%B0L, 1 ), & !specifies the number of rows of the matrix a.
                       size( inputData%ZonesData( iZone )%B0L, 2 ), & !specifies the number of columns of the matrix a.
                       ALPHA,                                       & !specifies the scalar alpha.
                       inputData%ZonesData( iZone )%B0L,            & !matrix A, SIZE (LDA, N).
                       size( inputData%ZonesData( iZone )%B0L, 1 ), & !leading dimension of array a.
                       inputData%ZonesData( iZone )%CB,             & !vector X.
                       1,                                           & !specifies the increment for the elements of X.
                       BETA,                                        & !specifies the scalar beta. 
                       inputData%ZonesData( iZone )%XB,             & !vector Y.
                       1 );                                           !specifies the increment for the elements of Y.
            
            
            
            !(3)--DGEMV: computes matrix-vector product: XB <= XB - A0L*UL            
            ALPHA = -1.0; 
            BETA  = +1.0;
            
            call DGEMV('N',                                         & !matrices a not be transposed or conjugate transposed before multiplication.
                       size( inputData%ZonesData( iZone )%A0L, 1 ), & !specifies the number of rows of the matrix a.
                       size( inputData%ZonesData( iZone )%A0L, 2 ), & !specifies the number of columns of the matrix a.
                       ALPHA,                                       & !specifies the scalar alpha.
                       inputData%ZonesData( iZone )%A0L,            & !matrix a, size (LDA, N).
                       size( inputData%ZonesData( iZone )%A0L, 1 ), & !leading dimension of array a.
                       UL,                                          & !vector X.
                       1,                                           & !specifies the increment for the elements of X.
                       BETA,                                        & !specifies the scalar beta. 
                       inputData%ZonesData( iZone )%XB,             & !vector Y.
                       1 );                                           !specifies the increment for the elements of Y.
            
            
            ! (4)--DGEMV: Computes matrix-vector product: XB <= XB + B00*YB            
            ALPHA = +1.0; 
            BETA  = +1.0;
            
            call DGEMV('N',                                         & !matrices a not be transposed or conjugate transposed before multiplication.
                       size( inputData%ZonesData( iZone )%B00, 1 ), & !specifies the number of rows of the matrix A.
                       size( inputData%ZonesData( iZone )%B00, 2 ), & !specifies the number of columns of the matrix A.
                       ALPHA,                                       & !specifies the scalar ALPHA.
                       inputData%ZonesData( iZone )%B00,            & !matrix A, SIZE (LDA, N).
                       size( inputData%ZonesData( iZone )%B00, 1 ), & !leading dimension of array a.
                       inputData%ZonesData( iZone )%Y,              & !vector X.
                       1,                                           & !specifies the increment for the elements of X.
                       BETA,                                        & !specifies the scalar BETA. 
                       inputData%ZonesData( iZone )%XB,             & !vector Y.
                       1 );                                           !specifies the increment for the elements of Y.
            
            
            
            ! (5)--DGEMV: Computes matrix-vector product: XB <= (A00)^(-1)*XB            
            
            ! Allocate temporary support array tempvec_            
            allocate( TempVEC_( size( inputData%ZonesData( iZone )%XB, 1 ) ), STAT = status );
            if(status /= 0) then
                print*, "Allocation of array TempVEC_ failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            else
                TempVEC_(:) = 0.0;
            end if
            
            ALPHA = +1.0; 
            BETA  = +0.0;
            
            write(*,*)'DGEMV: Calculation of boundary solution for Zone: ',iZone
            
            call DGEMV('N',                                         & !matrices a not be transposed or conjugate transposed before multiplication.
                       size( inputData%ZonesData( iZone )%A00, 1 ), & !specifies the number of rows of the matrix A.
                       size( inputData%ZonesData( iZone )%A00, 2 ), & !specifies the number of columns of the matrix A.
                       ALPHA,                                       & !specifies the scalar ALPHA.
                       inputData%ZonesData( iZone )%A00,            & !matrix A, SIZE (LDA, N).
                       size( inputData%ZonesData( iZone )%A00, 1 ), & !leading dimension of array A.
                       inputData%ZonesData( iZone )%XB,             & !vector X.
                       1,                                           & !specifies the increment for the elements of X.
                       BETA,                                        & !specifies the scalar BETA. 
                       TempVEC_,                                    & !vector Y.
                       1 );                                           !specifies the increment for the elements of Y.
            
            
            ! Fill the sub-matrix XB using new values (A00)^(-1)*XB            
            inputData%ZonesData( iZone )%XB = TempVEC_;                     
            
            ! Deallocate temporary support array TEMPMAT_
            if( allocated( TempVEC_ ) ) then
                deallocate( TempVEC_ );
                if(status /= 0) then
                    print*, "Deallocation of array TempVEC_ failed!.";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if            
          
            
            ! Deallocate temporary array "UL" to collect interface solution concerning the current zone            
            if( allocated( UL ) ) then
                deallocate( UL );
                if(status /= 0) then
                    print*, "Deallocation of array UL failed!.";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
            
        end do  ! End of loop statement over zones in the model for which a boundary solution was required          
        
    !$OMP ENDDO
    !$OMP END PARALLEL    
        
    end subroutine calculateOnSM    
!
!**********************************************************************************************    
    subroutine calculateOnDM(self, inputData, zonesBounSol, inputParamsData)    
!**********************************************************************************************
!            
        
        use, intrinsic :: iso_fortran_env;
        use omp_lib;        
        implicit none;    
        
        include 'mpif.h'            
    
        ! Declarative zone
        class(DDDMS_SMBoundarySolutionCalculator) :: self;
        class(DDDMS_DataContainer)                :: inputData;
        class(DDDMS_InputParams)                  :: inputParamsData;
        class(DDDMS_BSProc),allocatable           :: BSProc;
        integer(INT32), allocatable               :: zonesBounSol(:);                                        
        integer(INT32)                            :: i, j, k;                
        integer(INT32)                            :: IndepInterfNodesCount;
        integer(INT32)                            :: iZone;
        real(REAL64)                              :: ALPHA;
        real(REAL64)                              :: BETA;
        real(REAL64), allocatable                 :: U(:);                
        integer(INT32)                            :: intnod,intfnod,Bcnd,stat;
        integer(INT64)                            :: idim;
        logical                                   :: infnodstat=.false.     
        integer                                   :: nprocs,rank,numtasks,tag,dest,source,ierr;
        integer                                   :: iszXB,iszU, iszULI
        integer                                   :: status(MPI_STATUS_SIZE);
        
        common /MPILST/ rank,numtasks        
        
!-------------------------------
! Initialisation
!-------------------------------                
        nprocs=size(zonesBounSol);
        dest=0;
        source=0;
        iszXB=0;
        iszU=0;
        iszULI=0;
        dest=0;
        tag=1                
        
        if (rank.eq.0) then      
            if (numtasks.ne.nprocs+1) then  ! Check if the chosen # of processors is enough to address all the zones
                write(*,*)'Execute the job with',nprocs+1, '# of processors'
                stop
            endif             
            do i=1,nprocs
                iZone=zonesBounSol(i);
                dest=dest+1                               
                iszXB=size(inputData%ZonesData(iZone)%XB)
                iszU=size(inputData%U_InterfaceSolution) 
                call MPI_Send(iszU,1,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
                call MPI_Send(inputData%U_InterfaceSolution,iszU,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)                
                
!                 Receiving computed LS components for each zone from worker processors  
                call MPI_Recv(inputData%ZonesData(iZone)%XB,iszXB,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr)
                iszULI=size(inputData%ZonesData(iZone)%ULI);                
                call MPI_Recv(inputData%ZonesData(iZone)%ULI,iszULI,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr)
            enddo                        
        elseif(rank.ge.0) then                            
                write(*,*)'From Processor#',rank                        
                iZone=zonesBounSol(rank);
                iszXB=size(inputData%ZonesData(iZone)%XB);
                call MPI_Recv(iszU,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,ierr)
                allocate( U(iszU), STAT = stat);
                call MPI_Recv(U,iszU,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,ierr)            
                !            
                ! Compute LSC for each of the zone in the worker processors
                allocate( BSProc, STAT = stat, SOURCE = DDDMS_BSProc());

                iZone=zonesBounSol(rank);
                iszXB=size(inputData%ZonesData(iZone)%XB);
                call BSProc%calcBSProc(inputData%ZonesData(iZone)%A00,inputData%ZonesData(iZone)%A0L,      &
                                        inputData%ZonesData(iZone)%B00,inputData%ZonesData(iZone)%B0L,     &
                                        inputData%ZonesData(iZone)%DL,inputData%ZonesData(iZone)%CB,       &
                                        inputData%ZonesData(iZone)%Y, inputData%ZonesData(iZone)%XB,       &
                                        U, inputData, inputParamsData, iZone, rank                         &             
                                        )
                deallocate(U);
                                    
                ! Sending LS components to the master processor            
                call MPI_Send(inputData%ZonesData(iZone)%XB,iszXB,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,ierr)              
                iszULI=size(inputData%ZonesData(iZone)%ULI);            
                call MPI_Send(inputData%ZonesData(iZone)%ULI,iszULI,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,ierr)                    
        endif        
        
    end subroutine calculateOnDM
!
!**********************************************************************************************************    
    subroutine calcBSProc(self, subMA00, subMA0L, subMB00, subMB0L, DL, CB, YB, XB,       &
                                U, inputData, inputParamsData, iZone, rank                &
                         )
!**********************************************************************************************************    
!    
        class(DDDMS_BSProc)            :: self;      
        class(DDDMS_DataContainer)     :: inputData;
        class(DDDMS_InputParams)       :: inputParamsData;
        real(REAL64)                   :: subMA00(:,:),subMA0L(:,:),subMB00(:,:),subMB0L(:,:);
        real(REAL64)                   :: DL(:,:), CB(:), YB(:), XB(:);
        real(REAL64), allocatable      :: UL(:),TempVEC_(:), U(:);
        integer(INT32)                 :: iZone,iszDL(2),iszB0L(2),iszA0L(2),iszB00(2),iszA00(2);
        integer(INT32)                 :: IndepInterfNodesCount,status,countStep,j,k,idim;
        real(REAL64)                   :: ALPHA,BETA        
        integer(INT32)                 :: intnod,intfnod,Bcnd;
        logical                        :: infnodstat=.false.
        integer                        :: rank
    
! Initialization        
        IndepInterfNodesCount = count( inputData%IndependenInterfaceNodesArray(:, iZone) /= 0 ); ! Number of independent nodes lying on the interface for the current zone            
        if (inputParamsData%analysType .eq. 1) then     ! Problem dimension
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif
           
! Allocate temporary array "UL" to collect interface solution concerning the current zone            
        allocate( UL( idim*IndepInterfNodesCount ), STAT = status );
        if(status /= 0) then
            print*, "Failed allocation of UL array!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            UL(:) = 0.0;
        end if              
        
! Collect interface solution for the current zone from vector U_InterfaceSolution        
        countStep = 1;        
        do j = 1, size( inputData%IndependenInterfaceNodesArray(:, iZone), 1)                             
            if( inputData%IndependenInterfaceNodesArray(j, iZone) /= 0 ) then                                 
                if (idim .eq. 1) then
                    UL(countStep) = U(j);                    
                else
                    UL((idim*countStep - 2):(idim*countStep)) = U( (idim*j - 2):(idim*j) );
                endif              
                countStep = countStep + 1;                    
            end if                
        end do                    
        inputData%ZonesData(iZone)%ULI(:)=UL(:);
        
!   DGEMV: computes matrix-vector product: CB = (CB - DL*UL)        
            ALPHA = -1.0; 
            BETA  = +1.0;
            iszDL=shape(DL)            
            call DGEMV('N', iszDL(1),iszDL(2),ALPHA,DL,iszDL(1),UL,1,BETA,CB,1)
                        
            ! Rewrite interface results for the nodal directions where physical continuity BC's are applied
            ! Interface traction results(CB) are used instead of displacements
            countStep=0
            infnodstat=.false.
            do j = 1, size( inputData%IndependenInterfaceNodesArray(:, iZone), 1)
                intfnod=inputData%IndependenInterfaceNodesArray(j, iZone)
                infnodstat=inputParamsData%isNodMin(inputData,intfnod);
                if( intfnod /= 0) then
                    countStep=countStep+1                    
                    if (infnodstat) then    ! check if the node considered equals the lower node ID of the node pair in interface definition file                                                                  
                            do k=1,size(inputData%ZonesData(iZone)%CoordsNodesInZone(:,1))     ! Get zonal BC type index corresponding to global node ID
                                if (intfnod.eq.int(inputData%ZonesData(iZone)%CoordsNodesInZone(k,1))) then
                                    Bcnd=k                                    
                                    exit
                                endif
                            enddo                                                      
                            ! check if the node is assigned physical continuity BC's
                            if( inputData%ZonesData( iZone )%BCsType(Bcnd, 1) == 39 ) then   
                                if (idim .eq. 1) then
                                    inputData%ZonesData(iZone)%ULI(countStep)=inputData%ZonesData(iZone)%CB(countStep)
                                else
                                    inputData%ZonesData(iZone)%ULI(idim*countStep - 2)=inputData%ZonesData(iZone)%CB(idim*countStep - 2)
                                endif
                            endif
                            if (idim .ne. 1) then
                                if( inputData%ZonesData(iZone)%BCsType(Bcnd, 2) == 40 ) then
                                    inputData%ZonesData(iZone)%ULI(idim*countStep - 1)=inputData%ZonesData(iZone)%CB(idim*countStep - 1)
                                endif                        
                                if( inputData%ZonesData(iZone)%BCsType(Bcnd, 3) == 41 ) then
                                    inputData%ZonesData(iZone)%ULI(idim*countStep)=inputData%ZonesData(iZone)%CB(idim*countStep)         
                                endif
                            endif
                    endif
                end if
            enddo            
            
! DGEMV: computes matrix-vector product: XB = B0L*(CB - DL*UL)            
            ALPHA = +1.0; 
            BETA  = +0.0;
            iszB0L=shape(subMB0L)            
            call DGEMV('N',iszB0L(1),iszB0L(2),ALPHA,subMB0L,iszB0L(1),CB,1,BETA,XB,1)
            
            
! DGEMV: computes matrix-vector product: XB = XB - A0L*UL            
            ALPHA = -1.0; 
            BETA  = +1.0;
            iszA0L=shape(subMA0L)            
            call DGEMV('N',iszA0L(1),iszA0L(2),ALPHA,subMA0L,iszA0L(1),UL,1,BETA,XB,1)                        

! DGEMV: computes matrix-vector product: XB = XB + B00*YB
            ALPHA = +1.0; 
            BETA  = +1.0;
            iszB00=shape(subMB00)
            call DGEMV('N',iszB00(1),iszB00(2),ALPHA,subMB00,iszB00(1),YB,1,BETA,XB,1)            

! DGEMV: computes matrix-vector product: XB = (A00)^(-1)*XB
            
            ! Allocate temporary support array TEMPVEC_            
            allocate( TempVEC_(size(XB)), STAT = status );
            if(status /= 0) then
                print*, "Allocation of array TempVEC_ failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            else
                TempVEC_(:) = 0.0;
            end if
            
            write(*,*)'DGEMV: Calculation of boundary solution for Zone: ',iZone
            ALPHA = +1.0; 
            BETA  = +0.0;
            iszA00=shape(subMA00)            
            call DGEMV('N',iszA00(1),iszA00(2),ALPHA,subMA00,iszA00(1),XB,1,BETA,TempVEC_,1)            
            
            ! Fill the sub-matrix XB using new values (A00)^(-1)*XB            
            XB=TempVEC_;                        
            
            ! Deallocate temporary support array TEMPMAT_            
            if( allocated( TempVEC_ ) ) then
                deallocate( TempVEC_ );
                if(status /= 0) then
                    print*, "Deallocation of array TempVEC_ failed!.";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
            
            ! Deallocate temporary array "UL" to collect interface solution concerning the current zone            
            if( allocated( UL ) ) then
                deallocate( UL );
                if(status /= 0) then
                    print*, "Deallocation of array UL failed!.";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if            
    end subroutine calcBSProc
    
end module BoundarySolutionCalculator_Module
    