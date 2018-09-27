!  BEASYMain.f90 
!
!  FUNCTIONS:
!  BEASYMain - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: BEASYMain
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program DirectDDMMain        
        
        use iso_fortran_env;
        use Utility_Module;     
        use DataStructure_Module;
        use SolverController_Module;
        use DataLoader_Module;
        use LocalSchurComplementCalculator_Module;
        use GlobalSchurComplementAssembler_Module;
        use InterfaceSolutionCalculator_Module;
        use BoundarySolutionCalculator_Module;
        use getResultsVector_Module;
        
        include 'mpif.h'    
        
!        implicit none;        
        
        ! Declarative zone
        
        character(len = 200)                                     :: inputSolverPathName;      !==>VARIABLE TO DELETE?
        character(len = 200)                                     :: inputModelFileName;       !==>VARIABLE TO DELETE?
        character(len = 200)                                     :: inputPathFileName;        !==>VARIABLE TO DELETE?
        character(len = 200)                                     :: inputBlockFileName;       !==>VARIABLE TO DELETE?
        character(len = 200)                                     :: currentLine;              !==>VARIABLE TO DELETE?
        character(len = 100)                                     :: currCommandArg;           !==>VARIABLE TO DELETE?
        character(len = 200)                                     :: commandLine;              !==>VARIABLE TO DELETE?
        character(len = 100),                        allocatable :: batchFileCommands(:);
        integer(INT32)                                           :: HPC_Arch;
        integer(INT32)                                           :: i;
        integer(INT32)                                           :: j;
        integer(INT32)                                           :: k;
        integer(INT32)                                           :: status;
        integer(INT64)                                           :: count;
        integer(INT32)                                           :: unit
        integer(INT64)                                           :: file_size;
        integer(INT64)                                           :: file_sizeMod;
        integer(INT64)                                           :: file_sizePath;
        integer(INT64)                                           :: file_sizeBlock;
        integer(INT32)                                           :: file_sizeSolver;
        integer(INT32)                                           :: indexPath;
        integer(INT32),                              allocatable :: zonesBounSol(:);
        integer(INT32)                                           :: currBatchCommand;
        logical                                                  :: file_exists;
        logical                                                  :: file_existsMod;
        logical                                                  :: file_existsPath;
        logical                                                  :: file_existsBlock;
        logical                                                  :: file_existsSolver;
        logical,                                     allocatable :: chechBatchCommandsSatus(:);
        class(DDDMS_SolverController),               allocatable :: controller;
        class(DDDMS_DataLoaderFromFile),             allocatable :: loader;
        type(DDDMS_DataContainer),                   allocatable :: dataContainer;
        class(DDDMS_InputParams),                    allocatable :: inputParamsData;    
        class(DDDMS_LocalSchurComplementCalculator), allocatable :: LocSchurComp_Calculator;
        class(DDDMS_GlobalSchurComplementAssembler), allocatable :: GloSchurComp_Assembler;
        class(DDDMS_InterfaceSolutionCalculator),    allocatable :: InterSol_Calculator;
        class(DDDMS_BoundarySolutionCalculator),     allocatable :: BoundSol_Calculator;
#ifdef __WIN64
	    type(DDDMS_getResultsVector),               allocatable :: ResVec;
#elif __LINUX__
		type(DDDMS_getResultsVector)		                 :: ResVec;
#endif    
        real(REAL64)                                             :: start_time, stop_time;
        integer(INT64)                                           :: nzones;
! APB
        integer(INT64)                                           :: ic1,ic2,elsc
        integer                                                  :: ierr,ierror
        integer                                                  :: rank, numtasks
        logical                                                  :: flag;
        common /MPILST/ rank,numtasks
! APB        
        
! MPI initialization                
        rank=0                        
        
! call utility procedure to collect data from command line arguments
        
        call Utility_CommandParameters(inputSolverPathName, &
                                       inputModelFileName,  &
                                       inputPathFileName,   &
                                       inputBlockFileName,  &
                                       batchFileCommands,   &
                                       HPC_Arch);     
        if (HPC_Arch.eq.2) then
            call MPI_INIT(ierr)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
        elseif (HPC_Arch.eq.1) then
            rank=0
        endif
                                       
! Query for the existence of input files and their sizes
        inquire( FILE = inputModelFileName , EXIST = file_existsMod );     !CHECK IF THE INPUT MODEL FILE INPUTMODELFILENAME EXISTS.
        inquire( FILE = inputModelFileName , SIZE = file_sizeMod );        !CHECK IF THE INPUT MODEL FILE INPUTMODELFILENAME IS EMPTY.
            
        inquire( FILE = inputPathFileName , EXIST = file_existsPath );     !CHECK IF THE INPUT FILE INPUTPATHFILENAME EXISTS.
        inquire( FILE = inputPathFileName , SIZE = file_sizePath );        !CHECK IF THE INPUT FILE INPUTPATHFILENAME IS EMPTY.
            
        inquire( FILE = inputBlockFileName , EXIST = file_existsBlock );   !CHECK IF THE DUMP SIZE INPUT FILE INPUTBLOCKFILENAME EXISTS.
        inquire( FILE = inputBlockFileName , SIZE = file_sizeBlock );      !CHECK IF THE DUMP SIZE INPUT FILE INPUTBLOCKFILENAME IS EMPTY.
            
        inquire( FILE = inputSolverPathName , EXIST = file_existsSolver ); !CHECK IF THE DUMP SIZE INPUT FILE INPUTBLOCKFILENAME EXISTS.
        inquire( FILE = inputSolverPathName , SIZE = file_sizeSolver );    !CHECK IF THE DUMP SIZE INPUT FILE INPUTBLOCKFILENAME IS EMPTY.                    
        

! Process batch file and execute the modules input by the user                
        if( (( HPC_Arch == 1 ) .or. ( HPC_Arch == 2 ))       .AND. &   !HPC ARCHITECTURE CAN BE SHARED MEMORY OR DISTRIBUTED MEMORY ONLY
            ( file_existsMod .AND. (file_sizeMod /= 0) )     .AND. &   !THE INPUT MODEL FILE INPUTMODELFILENAME EXISTS
            ( file_existsPath .AND. (file_sizePath /= 0) )   .AND. &   !THE INPUT FILE INPUTMODELFILENAME EXISTS
            ( file_existsBlock .AND. (file_sizeBlock /= 0) ) ) THEN    !THE INPUT FILE INPUTBLOCKFILENAME EXISTS
            
            ! Create solvercontroller            
            allocate( controller, STAT = status, SOURCE = DDDMS_SolverController( HPC_Arch ) );
            if (status /= 0) then
                print*, "Failed allocation of DDM controller!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if  

            
            ! Obtain data loader (load from file) from solvercontroller object            
            allocate( loader, STAT = status, SOURCE = controller%getDataLoader( inputPathFileName ) );
            if (status /= 0) then
                print*, "Failed allocation of DDM loader!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if 
        
        do k = 1, size( batchFileCommands, 1)
                
                currBatchCommand = 0;
                
                call To_upper( batchFileCommands(k) ); 

                currBatchCommand = currBatchCommand + merge(1, 0, index( batchFileCommands(k), "BEASY"      ) /= 0);
                currBatchCommand = currBatchCommand + merge(2, 0, index( batchFileCommands(k), "CALC_LSC"   ) /= 0);
                currBatchCommand = currBatchCommand + merge(3, 0, index( batchFileCommands(k), "ASSEMB_GSC" ) /= 0);
                currBatchCommand = currBatchCommand + merge(4, 0, index( batchFileCommands(k), "CALC_IS"    ) /= 0);
                currBatchCommand = currBatchCommand + merge(5, 0, index( batchFileCommands(k), "CALC_BS"    ) /= 0);
                currBatchCommand = currBatchCommand + merge(6, 0, index( batchFileCommands(k), "RES_VEC"    ) /= 0);
                
                commandLine      = trim( '"'//trim( inputSolverPathName(:) )//'"'//" -rootname "//'"'//trim( inputModelFileName( 1:index(inputModelFileName, '.dat') - 1) )//'"' ); 

! The Beasy command line here enables the creation of zonal A and B matrices. 
! To be revoked if the matrices are to be created during the execution of code                
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                
!                commandLine      = trim( '"'//trim( inputSolverPathName(:) )//'"'//" -rootname "//'"'//trim( inputModelFileName( 1:index(inputModelFileName, '.dat') - 1) )//'"'//" -DDMtask createZoneAB" ); 
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                
!              
               
                select case ( currBatchCommand )
                    case ( 1 )                        
                        if (rank.eq.0) then
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            write(OUTPUT_UNIT, '("BEASY Solver is being called............................: ")');
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            flush OUTPUT_UNIT;                                                     
                            call system_clock(count=ic1)  
    !                        call execute_command_line ( commandLine, wait = .true. ); !SYNCHRONOUS CALL TO BEASY SOLVER.
                            call system_clock(count=ic2)
                            elsc=ic2 - ic1                        
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';!                        
                            write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                            write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") (elsc);
                            write(OUTPUT_UNIT, '("-[s]) ")');                                             
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';
                            write(OUTPUT_UNIT, *);                        
                        endif
                    
                    case( 2 )
                        if (rank.eq.0) then
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            write(OUTPUT_UNIT, '("Loading input data required for the computation......: ")');
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            flush OUTPUT_UNIT;
                        endif                        
                        
                        ! Check data in folder pathname passed in input using check procedure of the loader object                
                        if( loader%CheckData() ) then                       
                            
                            ! Create input parameters data container                            
                            allocate( inputParamsData, STAT = status, SOURCE = DDDMS_InputParams() );
                            if (status /= 0) then
                                print*, "Failed allocation of inputParamsData data structure!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if                         
                            
                            ! Load inputparams                            
                            call inputParamsData%loadInputParams( inputPathFileName, inputBlockFileName );                        
                            
                            ! Create data container                            
#ifdef __WIN64
                            allocate( dataContainer, STAT = status, SOURCE = DDDMS_DataContainer() );
#elif __LINUX__
							allocate( dataContainer, STAT = status);
#endif
                            if (status /= 0) then
                                print*, "Failed allocation of DDM controller!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if                                                      
                            
                            ! Load data from a specific folder using loaddata procedure of the loader object                                                        
                            call system_clock(count=ic1)  
                            call loader%LoadData( dataContainer, inputParamsData );     ! APB: Polymorphism with type bound procedure
                            call system_clock(count=ic2)
                            elsc=ic2 - ic1              
                            if (rank.eq.0) then
                                write(OUTPUT_UNIT,*)'---------------------------------------------------------';
                                write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                                write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") elsc;
                                write(OUTPUT_UNIT, '("-[s]) ")');                                                        
                                write(OUTPUT_UNIT,*)'---------------------------------------------------------';
                                write(OUTPUT_UNIT, *);
                            endif
                            if (rank.eq.0) then
                                write(OUTPUT_UNIT,'("************************************************************************")');
                                write(OUTPUT_UNIT, '("Local Schur Complement calculation is being performed...: ")');
                                write(OUTPUT_UNIT,'("************************************************************************")');
                            endif

                            ! Obtain local schur complement calculator object (shared memory architecture) from solvercontroller object                            
#ifdef __WIN64
                            allocate( LocSchurComp_Calculator, STAT = status, SOURCE = controller%getLocalSchurComplementCalculator( HPC_Arch ) );
#elif __LINUX__
							allocate( LocSchurComp_Calculator, STAT = status, SOURCE = DDDMS_SMLocalSchurComplementCalculator());
#endif

                            if (status /= 0) then
                                print*, "Failed allocation of DDM loader!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if
                            
                            ! Calculate local schur complement for all zones involved in the model                            
                            call cpu_time(start_time);
                            call system_clock(count=ic1)  
                            if (HPC_ARCH.eq.1 .and. rank.eq.0) then                            
                                write(*,*)'Computation occurs in parallel mode in a shared environment'
                                call LocSchurComp_Calculator%calculateSM( dataContainer, inputParamsData );
                            elseif (HPC_ARCH.eq.2) then                            
                                if (rank.eq.0) write(*,*)'Computation occurs in parallel mode in a distributed environment'                                
                                call LocSchurComp_Calculator%calculateDM( dataContainer, inputParamsData );
                            endif
                            
                            
                            call system_clock(count=ic2)
                            call cpu_time(stop_time)
                            elsc=ic2 - ic1                       
                            if (rank.eq.0) then
                                write(OUTPUT_UNIT,*)'-------------------------------------------------';                            
                                write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                                write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") elsc;
                                write(OUTPUT_UNIT, '("-[s]) ")');                                                        
                                write(OUTPUT_UNIT,*)'-------------------------------------------------';
                                write(OUTPUT_UNIT, *);
                            endif
                        else
                            print*, "Check function returned an error.";
                            print*, "Please, check input files path name folder ";
                            pause;
                            stop;
                        endif ! End of if statement checking data
                        
                    case( 3 )
                        if (rank.eq.0) then
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            write(OUTPUT_UNIT, '("Global Schur Complement assembling is being performed...: ")');
                            write(OUTPUT_UNIT,'("************************************************************************")');                        
                            flush OUTPUT_UNIT;                                                
                            
                            ! Obtain global schur complement assembler from solvercontroller object                        
#ifdef __WIN64
                            allocate( GloSchurComp_Assembler, STAT = status, SOURCE = controller%getGlobalSchurComplementAssembler() );
#elif __LINUX__
                            allocate( GloSchurComp_Assembler, STAT = status, SOURCE = DDDMS_GlobalSchurComplementAssembler());
#endif

                            if (status /= 0) then
                                print*, "Failed allocation of DDM assembler!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if                          
                            
                            ! Call subroutine "assembly" to assemble all contributes provided by all zones in the model                            
                            call cpu_time(start_time);
                            call system_clock(count=ic1)                                                      
                            call GloSchurComp_Assembler%Assembly( dataContainer,  inputParamsData);
                            call system_clock(count=ic2)
                            call cpu_time(stop_time)
                            elsc=ic2 - ic1                        
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';                                                        
                            write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                            write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") elsc;
                            write(OUTPUT_UNIT, '("-[s]) ")');                                                    
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';
                            write(OUTPUT_UNIT, *);
                        endif
!                        
                    case( 4 )    
                        if (rank.eq.0) then
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            write(OUTPUT_UNIT, '("Interface solution is being calculated..................: ")');
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            flush OUTPUT_UNIT;
    !                                                    
                            ! Obtain interface solution calculator from solvercontroller object                            
#ifdef __WIN64
                            allocate( InterSol_Calculator, STAT = status, SOURCE = controller%getInterfaceSolutionCalculator(HPC_Arch) );
#elif __LINUX__
                            allocate( InterSol_Calculator, STAT = status, SOURCE = DDDMS_SMInterfaceSolutionCalculator());
#endif
                            if (status /= 0) then
                                print*, "Failed allocation of DDM interface solution calculator!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if      
                           
                            ! Call subroutine "calculate" to perform calculation of interface solution on interface                            
                            write(*,*)'Solve algebraic equations SGlob*UL=gGLOB for inteface solution'
                            call cpu_time(start_time);
                            call system_clock(count=ic1)                                                               
                            call InterSol_Calculator%Dense_InterfaceSolution( dataContainer,inputParamsData);  ! Solve dense systems using MKL blas/lapack routines
    !                        call InterSol_Calculator%Sparse_InterfaceSolution( dataContainer,inputParamsData); ! Solve sparse systems using MKL pardiso routines
                            call system_clock(count=ic2)
                            call cpu_time(stop_time)
                            elsc=ic2 - ic1                        
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';        !                            
                            write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                            write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") elsc;
                            write(OUTPUT_UNIT, '("-[s]) ")');                                                        
                            write(OUTPUT_UNIT,*)'-------------------------------------------------';
                            write(OUTPUT_UNIT, *);                        
                        endif
                        
                    case( 5 )    
                            if (rank.eq.0) then
                                write(OUTPUT_UNIT,'("************************************************************************")');
                                write(OUTPUT_UNIT, '("Boundary solution is being calculated...................: ")');
                                write(OUTPUT_UNIT,'("************************************************************************")');
                                flush OUTPUT_UNIT;                            
                            endif
                            
                            ! Obtain boundary solution calculator from solvercontroller object                            
#ifdef __WIN64
                            allocate( BoundSol_Calculator, STAT = status, SOURCE = controller%getBoundarySolutionCalculator(HPC_Arch) );                            
#elif __LINUX__
                            allocate( BoundSol_Calculator, STAT = status, SOURCE = DDDMS_SMBoundarySolutionCalculator());
#endif
                            if (status /= 0) then
                                print*, "Failed allocation of DDM boundary solution calculator!";
                                print*, "Errore code: ", status;
                                pause;
                                stop;
                            end if                              
                            
                            ! Perform calculation of boundary solution on zones involved in the model                            
                            allocate(zonesBounSol, SOURCE = Utility_ZoneBoundarySolution( batchFileCommands(k), inputParamsData ) );                                                                         
                            
                            call system_clock(count=ic1)                           
                            
                            if (HPC_ARCH.eq.1 .and. rank.eq.0) then                            
                                write(*,*)'Computation occurs in parallel mode in a shared environment'
                                call BoundSol_Calculator%calculateSM( dataContainer, zonesBounSol, inputParamsData );
                            elseif (HPC_ARCH.eq.2) then                            
                                if (rank.eq.0) write(*,*)'Computation occurs in parallel mode in a distributed environment'                                
                                call BoundSol_Calculator%calculateDM( dataContainer, zonesBounSol, inputParamsData );
                            endif                            
                            
                            call system_clock(count=ic2)                                                    
                            elsc=ic2 - ic1                        
                            if (rank.eq.0) then
                                write(OUTPUT_UNIT,*)'-------------------------------------------------';                            
                                write(OUTPUT_UNIT, '("(Elapsed System Time: ")', ADVANCE = "NO");
                                write(OUTPUT_UNIT, '(i20)', ADVANCE = "NO") elsc;
                                write(OUTPUT_UNIT, '("-[s]) ")');                                                    
                                write(OUTPUT_UNIT,*)'-------------------------------------------------';
                                write(OUTPUT_UNIT, *);
                            endif
                    case( 6 )    
                        if (rank.eq.0) then
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            write(OUTPUT_UNIT,*)'Creating results vector......';
                            write(OUTPUT_UNIT,'("************************************************************************")');
                            flush OUTPUT_UNIT;
                            nzones=size(dataContainer%ZonesData);
#ifdef __WIN64
                            allocate( ResVec, STAT = status, SOURCE = DDDMS_getResultsVector() );
#endif
                            allocate( ResVec%ZonalResults(nzones), STAT = status);
                            if (status.gt.0) then
                                write(*,*)'Allocation failed for results vector object!....'
                                stop;
                            endif
                            call ResVec%getResultsVector(dataContainer, inputParamsData);
                        endif
                        
                    case default
                        print *, "Entered command was not recognised.";
                        print *, "Check input command.";
                        
                end select
            end do ! End of loop over input batch file commands         
            
        else
            print*, "Check HPC hardware architecture selected.";
            pause;
            stop;
        endif ! close if statement on checking HPC hardware architecture.
        if (rank.eq.0)then
            write(*,*)'------------------------'
            write(*,*)'Execution halted'     
            write(*,*)'------------------------'
        endif
        if (HPC_Arch.eq.2) call MPI_FINALIZE(ierr)            
!        
    end program DirectDDMMain

