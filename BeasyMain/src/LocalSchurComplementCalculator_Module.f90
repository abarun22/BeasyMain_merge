    module LocalSchurComplementCalculator_Module
    use, intrinsic :: iso_fortran_env;
    use DataStructure_Module;
    use Utility_Module;                                        !==> TO DELETE?
    use omp_lib;
    implicit none;
    private;
    public :: DDDMS_LSCProc;
    public :: DDDMS_LocalSchurComplementCalculator;
    public :: DDDMS_SMLocalSchurComplementCalculator;    

    
    ! DDDMS_LocalSchurComplementCalculator definition
    ! Abstract DDDMS_LocalSchurComplementCalculator user-defined data type definition    
    type :: DDDMS_LSCProc        
    contains
        procedure                       :: calcLSCProc;        
    end type DDDMS_LSCProc
    
    type, abstract :: DDDMS_LocalSchurComplementCalculator
        
    contains 
        procedure(CalculateABTRACT), deferred :: calculateSM;
        procedure(CalculateABTRACT), deferred :: calculateDM;        
    end type DDDMS_LocalSchurComplementCalculator
    
    
    ! DDDMS_LocalSchurComplementCalculator abstract user-defined data type procedure definition
    abstract interface
        subroutine CalculateABTRACT(self, dataContainer, inputParamsData)
            import;
            class(DDDMS_LocalSchurComplementCalculator) :: self;
            class(DDDMS_DataContainer)                  :: dataContainer;
            class(DDDMS_InputParams)                    :: inputParamsData;
        end subroutine CalculateABTRACT
    end interface

    
    ! Concrete DDDMS_LocalSchurComplementCalculator user-defined data type definition for shared memory architecture    
    type, extends(DDDMS_LocalSchurComplementCalculator) :: DDDMS_SMLocalSchurComplementCalculator
        
    contains 
        procedure :: calculateSM => calculateOnSM_SH;
        procedure :: calculateDM => calculateOnSM_DM;
        
    end type DDDMS_SMLocalSchurComplementCalculator

    
    ! DDDMS_SMLocalSchurComplementCalculator user-defined data type constructor definition    
    interface DDDMS_SMLocalSchurComplementCalculator
        module procedure NewDDDMS_SMLocalSchurComplementCalculator; ! add constructor to DDDMS_SMLocalSchurComplementCalculator generic interface
    end interface DDDMS_SMLocalSchurComplementCalculator   
    
contains

    
    ! DDDMS_SMLocalSchurComplementCalculator user-defined data type constructor implementation    
    type(DDDMS_SMLocalSchurComplementCalculator) function NewDDDMS_SMLocalSchurComplementCalculator(self)
    
    implicit none;
        
    !********************************************
    !DECLARATIVE ZONE
    !********************************************
    class(DDDMS_SMLocalSchurComplementCalculator) :: self;
        
    !********************************************
    !BODY OF THE PROGRAM
    !********************************************        
        
    end function NewDDDMS_SMLocalSchurComplementCalculator

    
    ! Calculate procedure implementation    
    subroutine calculateOnSM_SH(self, dataContainer, inputParamsData)
    
        implicit none;
    
        
        ! Declarative zone        
        class(DDDMS_SMLocalSchurComplementCalculator) :: self;
        class(DDDMS_DataContainer)                    :: dataContainer;
        class(DDDMS_InputParams)                      :: inputParamsData;
        class(DDDMS_LSCProc),allocatable              :: LSCProc;
        integer(INT64)                                :: N;               !==>VARIABLE TO DELETE?
	    integer(INT64)                                :: i;               !==>VARIABLE TO DELETE?
        integer(INT64)                                :: j;               !==>VARIABLE TO DELETE?
        integer(INT64)                                :: k;               !==>VARIABLE TO DELETE?
        real(REAL64), allocatable                     :: TempMAT_(:, :);  !==>VARIABLE TO DELETE?
        real(REAL64), allocatable                     :: TempVEC_(:);     !==>VARIABLE TO DELETE?
        real(REAL64)                                  :: temp;            !==>VARIABLE TO DELETE?
        real(REAL64)                                  :: ALPHA;
        real(REAL64)                                  :: BETA;
        real(REAL64)                                  :: GAMMA;
        integer(INT32)                                :: numOfZones,stat;
        integer(INT32)                                :: fid;
        integer(INT64), allocatable                   :: pivot(:);
        real(REAL64),   allocatable                   :: work(:);
        integer(INT32)                                :: info;
        integer(INT32)                                :: status;
        character(len = 4)                            :: seqstring;            !==> VARIABLE TO DELETE.
        integer(INT64)                                :: thread_number,gnumthrd
        integer			                      		  :: iszB0L(2),iszA00(2),iszAL0(2),iszTMAT(2),iszALL(2);
	    integer		   	   		   					  :: iszBLL(2),iszDL(2),iszB00(2),iszA0L(2),iszBL0(2);
    
        
        ! Body of the program
        numOfZones = size( dataContainer%ZonesData, 1 );        
       
        !$OMP PARALLEL PRIVATE(thread_number, N, pivot, work, TempMAT_, info, status, ALPHA, BETA, seqstring, TempVEC_)                               
        !$OMP DO
            do thread_number = 1, numOfZones                                    
                allocate( LSCProc, STAT = stat, SOURCE = DDDMS_LSCProc());
                call LSCProc%calcLSCProc(dataContainer%ZonesData(thread_number)%A00, dataContainer%ZonesData(thread_number)%AL0,   &
                                    dataContainer%ZonesData(thread_number)%A0L,  dataContainer%ZonesData(thread_number)%ALL,   &
                                    dataContainer%ZonesData(thread_number)%Y,    inputParamsData%gMatScaleFact,       &
                                    dataContainer%ZonesData(thread_number)%B00,  dataContainer%ZonesData(thread_number)%BL0,   &
                                    dataContainer%ZonesData(thread_number)%B0L,  dataContainer%ZonesData(thread_number)%BLL,   &
                                    dataContainer%ZonesData(thread_number)%DL,   dataContainer%ZonesData(thread_number)%CB     &
                                    )                                                            
            end do !END OF DO-WHILE LOOP STATEMENT INTERNAL TO PARALLEL ZONE.            
        !$OMP ENDDO                
        !$OMP END PARALLEL
       
    end subroutine calculateOnSM_SH
!
!**********************************************************************************************        
    subroutine calculateOnSM_DM(self, dataContainer, inputParamsData)
!**********************************************************************************************        
!   Objective: Message passing among the processors to handle the LSC computation concurrently        
!**********************************************************************************************        
!         
        include 'mpif.h'            
        
! Declarative zone        
        class(DDDMS_SMLocalSchurComplementCalculator) :: self;
        class(DDDMS_DataContainer)                    :: dataContainer;
        class(DDDMS_InputParams)                      :: inputParamsData;
        class(DDDMS_LSCProc),allocatable              :: LSCProc;
        integer(INT64)                                :: N;              
	    integer(INT64)                                :: i;              
        integer(INT64)                                :: j;              
        integer(INT64)                                :: k;                              
        integer(INT32)                                :: iworkproc,stat;
        integer                                       :: st_count,ierr,ierror
        integer                                       :: iszAL0,iszA00,iszALL,iszY,iszBL0,iszB00,iszBLL,iszDL,iszCB
        integer                                       :: dest,orig,status(MPI_STATUS_SIZE)    
        integer                                       :: st_source,st_tag,source,leng,nprocs        
        integer                                       :: rank, numtasks
        logical                                       :: flag;
        real(REAL64), allocatable                     :: A00(:,:),A0L(:,:),B00(:,:),B0L(:,:),Y(:);
        integer                                       :: szA00(2),szA0L(2),szB00(2),szB0L(2)
        
        common /MPILST/ rank,numtasks        		
        
! Initialization
        nprocs=inputParamsData%numberOfZones
        dest=0
        source=0
        iworkproc=0        

!        
!-----------------------------------------------------------------------------------------        
! If rank=0: Master process handles the computation
!            Send the sub-arrays needed for LSC computation to worker processors
! Else:     The worker processors receives the sub-arrays
!           Does the LSC computation independently in each of the worker processors
!           Each of them Send their LS components to the master processor
!-----------------------------------------------------------------------------------------
!        
        if (rank.eq.0) then      
            if (numtasks.ne.nprocs+1) then  ! Check if the chosen # of processors is enough to address all the zones
                write(*,*)'Execute the job with',nprocs+1, '# of processors'
!                call MPI_Abort(MPI_COMM_WORLD,ierr)
                stop
            endif            
            do i=1,nprocs
                dest=dest+1
                iszA00=size(dataContainer%ZonesData(i)%A00)                
                iszY=size(dataContainer%ZonesData(i)%Y)                
                iszDL=size(dataContainer%ZonesData(i)%DL)
                iszCB=size(dataContainer%ZonesData(i)%CB)

                ! Receiving computed LS components for each zone from worker processors
                iworkproc=iworkproc+1                
                call MPI_Recv(dataContainer%ZonesData(i)%DL,iszDL,MPI_DOUBLE_PRECISION,dest,1,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(dataContainer%ZonesData(i)%CB,iszCB,MPI_DOUBLE_PRECISION,dest,1,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(dataContainer%ZonesData(i)%A00,iszA00,MPI_DOUBLE_PRECISION,dest,1,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(dataContainer%ZonesData(i)%Y,iszY,MPI_DOUBLE_PRECISION,dest,1,MPI_COMM_WORLD,status,ierr)
                
                st_source=status(MPI_SOURCE)                
            enddo                        
        elseif(rank.ge.0) then                            
            write(*,*)'From Processor#',rank                        
            iszA00= size(dataContainer%ZonesData(rank)%A00)                        
            iszY  = size(dataContainer%ZonesData(rank)%Y)            
            iszDL = size(dataContainer%ZonesData(rank)%DL)
            iszCB = size(dataContainer%ZonesData(rank)%CB)
            
            ! Compute LSC for each of the zone in the worker processors
            allocate( LSCProc, STAT = stat, SOURCE = DDDMS_LSCProc());
            call LSCProc%calcLSCProc(dataContainer%ZonesData(rank)%A00, dataContainer%ZonesData(rank)%AL0,   &
                                    dataContainer%ZonesData(rank)%A0L,  dataContainer%ZonesData(rank)%ALL,   &
                                    dataContainer%ZonesData(rank)%Y,    inputParamsData%gMatScaleFact,       &
                                    dataContainer%ZonesData(rank)%B00,  dataContainer%ZonesData(rank)%BL0,   &
                                    dataContainer%ZonesData(rank)%B0L,  dataContainer%ZonesData(rank)%BLL,   &
                                    dataContainer%ZonesData(rank)%DL,   dataContainer%ZonesData(rank)%CB     &
                                    )                        
                                    
            ! Sending LS components to the master processor            
            call MPI_Send(dataContainer%ZonesData(rank)%DL,iszDL,MPI_DOUBLE_PRECISION,source,1,MPI_COMM_WORLD,ierr)
            call MPI_Send(dataContainer%ZonesData(rank)%CB,iszCB,MPI_DOUBLE_PRECISION,source,1,MPI_COMM_WORLD,ierr)
            call MPI_Send(dataContainer%ZonesData(rank)%A00,iszA00,MPI_DOUBLE_PRECISION,source,1,MPI_COMM_WORLD,ierr)
            call MPI_Send(dataContainer%ZonesData(rank)%Y,iszY,MPI_DOUBLE_PRECISION,source,1,MPI_COMM_WORLD,ierr)
        endif                
       
    end subroutine calculateOnSM_DM
!
!*****************************************************************************************************************
    subroutine calcLSCProc(self, subMA00, subMAL0, subMA0L, subMALL, YB, gMScf, &
                                 subMB00, subMBL0, subMB0L, subMBLL, DL, CB     &
                          )
!*****************************************************************************************************************
!                                 
        
        implicit none;
        
        class(DDDMS_LSCProc)            :: self;         
        integer(INT64)                  :: N
        integer(INT32)                  :: status,info
        integer(INT64)                  :: iszAL0(2),iszA0L(2),iszA00(2),iszALL(2)
        integer(INT64)                  :: iszB0L(2),iszBLL(2),iszB00(2),iszBL0(2)
        integer(INT64), allocatable     :: pivot(:);
        real(REAL64)                    :: subMA00(:,:),subMAL0(:,:),subMA0L(:,:),subMALL(:,:);
        real(REAL64)                    :: subMB00(:,:),subMBL0(:,:),subMB0L(:,:),subMBLL(:,:);
        real(REAL64)                    :: YB(:), DL(:,:), CB(:);
        real(REAL64)                    :: gMScf,ALPHA,BETA,GAMMA
        real(REAL64),   allocatable     :: work(:),TempMAT_(:, :),TempVEC_(:);

!   Get Row size required for MKL routines
        N=size(subMA00,1);
        
!  Scale right hand side vector by the factor defined in (*.key.txt) input files
        YB=YB/gMScf;
        
! Perform allocation for the pivot vector         
        allocate(pivot(N), STAT=status);
        if(status /= 0) then
            print*, "Allocation of array pivot failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            pivot(:) = 0.0;
        end if

! Allocate work vector
        allocate( work(N), STAT=status);
        if(status /= 0) then
            print*, "Allocation of array work failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            work(:) = 0.0;
        end if
        
! Perform matrix inversion (A00)^(-1) in two steps        
        call DGETRF(N, N, subMA00, N, pivot, info);   ! Step 1: LU factorization
        if(info /= 0) then
            pause;
            stop 'Matrix A00 is numerically singular!'
        end if
                
        call DGETRI(N, subMA00, N, pivot, work, N, info);  ! Step 2: Matrix Inversion
        if(info /= 0) then
            pause;
            stop 'Matrix inversion failed!'
        end if

        iszAL0=shape(subMAL0)
        iszA00=shape(subMA00)
        allocate(TempMAT_(iszAL0(1), iszAL0(2)), STAT=status)
        if(status /= 0) then
            print*, "Allocation of array _TempMAT failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            TempMAT_(:, :) = 0.0;
        end if
        
        ALPHA = 1.0; 
        BETA  = 0.0;
        GAMMA = 1.0;

! Perform AL0~ = AL0*(A00)^(-1) with the format [C = Alpha*A*B + Beta*C]        
        call DGEMM('N','N',iszAL0(1),iszA00(2),iszA00(1),               &
                    ALPHA,subMAL0,iszAL0(1),subMA00,iszA00(1),BETA,     &
                    TempMAT_, size(TempMAT_,1)                          & 
                  )
                    
        subMAL0=TempMAT_; ! Fill in the sub-matrix AL0 using new values AL0*(A00)^(-1)
        

! Deallocate temporary support array TEMPMAT_
        if( allocated( TempMAT_ ) ) then
            deallocate( TempMAT_ );
            if(status /= 0) then
                print*, "Deallocation of array TempMAT_ failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if

! Perform [BLL = AL0~*B0L - BLL] with the format [C = Alpha*A*B + Beta*C]
        ALPHA = +1.0; 
        BETA  = -1.0;
        iszB0L=shape(subMB0L)
        iszBLL=shape(subMBLL)
        call DGEMM('N','N',iszAL0(1),iszB0L(2),iszB0L(1),               &
                    ALPHA,subMAL0,iszAL0(1),subMB0L,iszB0L(1),BETA,     &
                    subMBLL, iszBLL(1)                                  & 
                  )

! Perform [ALL = ALL - AL0~*A0L] with the format [C = Alpha*A*B + Beta*C]
        ALPHA = -1.0; 
        BETA  = +1.0;        
        iszA0L=shape(subMA0L)
        iszALL=shape(subMALL)
        call DGEMM('N','N',iszAL0(1),iszA0L(2),iszA0L(1),               &
                    ALPHA,subMAL0,iszAL0(1),subMA0L,iszA0L(1),BETA,     &
                    subMALL, iszALL(1)                                  & 
                  )
                    
! Allocate temporary support array TEMPMAT_                    
        allocate( TempMAT_(iszB0L(2),iszB0L(1)), STAT = status );
        if(status /= 0) then
            print*, "Allocation of array _TempMAT failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            TempMAT_(:, :) = 0.0;
        end if

! Perform [BL0 = BL0 - AL0~*B00] with the format [C = Alpha*A*B + Beta*C] 
        ALPHA = -1.0; 
        BETA  = +1.0;
        iszB00=shape(subMB00)
        iszBL0=shape(subMBL0)
        call DGEMM('N','N',iszAL0(1),iszB00(2),iszB00(1),               &
                    ALPHA,subMAL0,iszAL0(1),subMB00,iszB00(1),BETA,     &
                    subMBL0, iszBL0(1)                                  & 
                  )

! Deallocate temporary support array TEMPMAT_
        if( allocated( TempMAT_ ) ) then
            deallocate( TempMAT_ );
            if(status /= 0) then
                print*, "Deallocation of array TempMAT_ failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if

! Deallocate temporary support array pivot
        if( allocated( pivot ) ) then
            deallocate( pivot, STAT = status );
            if(status /= 0) then
                print*, "Deallocation of array pivot failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if
            
        
! Deallocate temporary support array work        
        if( allocated( work ) ) then
            deallocate( work, STAT = status );
            if(status /= 0) then
                print*, "Deallocation of array work failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if

! Perform inversion of the matrix BLL and replace it with its inverse                
        N = iszBLL(1);                            

! Allocate the integer vector "pivot" containing the pivot indices                
        allocate(pivot(N), STAT=status);
        if(status /= 0) then
            print*, "Allocation of array pivot failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            pivot(:) = 0.0;
        end if            
                
! Allocate the real vector "work"                
        allocate(work(N), STAT=status);
        if(status /= 0) then
            print*, "Allocation of array work failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            work(:) = 0.0;
        end if
                
        call DGETRF(N, N, subMBLL, N, pivot, info);  ! Step 1: LU factorization
        if(info /= 0) then
            pause;
            stop 'Matrix BLL is numerically singular!'
        end if                                  
        
        call DGETRI(N, subMBLL, N, pivot, work, N, info);  ! Step 2: Matrix Inversion [BLL^(-1)]

! Deallocate temporary support array pivot
        if( allocated( pivot ) ) then
            deallocate( pivot, STAT = status );
            if(status /= 0) then
                print*, "Deallocation of array pivot failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if            
        
! Deallocate temporary support array work        
        if( allocated( work ) ) then
            deallocate( work, STAT = status );
            if(status /= 0) then
                print*, "Deallocation of array work failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if
        
        write(*,*)'DGEMV: Calculation of LHS LSC component DL'
        ALPHA = +1.0; 
        BETA  = +0.0;        
! Perform [DL = BLL^(-1)*ALL] with the format [C = Alpha*A*B + Beta*C]
        call DGEMM('N','N',iszBLL(1),iszALL(2),iszALL(1),               &
                    ALPHA,subMBLL,iszBLL(1),subMALL,iszALL(1),BETA,     &
                    DL, size(DL,1)                                      & 
                  )
                    
! Allocate temporary support array TEMPVEC_
        allocate(TempVEC_(iszAL0(1)), STAT=status);
        if(status /= 0) then
            print*, "Allocation of array TempVEC_ failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            TempVEC_(:) = 0.0;
        end if

! Perform [TempVEC_ = BL0*YB] with the format [C = Alpha*A*x+ Beta*y]
        ALPHA = +1.0; 
        BETA  = +0.0;
        call DGEMV('N',iszBL0(1),iszBL0(2),ALPHA,subMBL0,iszBL0(1),     &
                    YB,1,BETA,TempVEC_,1                                & 
                  )

! Perform [CB = BLL^(-1)*TempVEC_] with the format [C = Alpha*A*x+ Beta*y]                    
        write(*,*)'DGEMV: Calculation of LHS LSC component CB'
        ALPHA = +1.0; 
        BETA  = +0.0;
        call DGEMV('N',iszBLL(1),iszBLL(2),ALPHA,subMBLL,iszBLL(1),     &
                    TempVEC_,1,BETA,CB,1                                & 
                  )
            
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
                    
    end subroutine calcLSCProc


end module LocalSchurComplementCalculator_Module