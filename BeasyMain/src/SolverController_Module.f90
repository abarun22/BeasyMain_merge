module SolverController_Module
    use, intrinsic :: iso_fortran_env;
    use DataLoader_Module;
    use LocalSchurComplementCalculator_Module;
    use GlobalSchurComplementAssembler_Module;
    use InterfaceSolutionCalculator_Module;
    use BoundarySolutionCalculator_Module;
    implicit none;
    private;
    public DDDMS_SolverController;
    
    ! DDDMS_SolverController        
    ! DDDMS_SolverController user-defined data type definition    
    type :: DDDMS_SolverController
        integer :: solverType;
    contains
        procedure :: getDataLoader;
        procedure :: getLocalSchurComplementCalculator;
        procedure :: getGlobalSchurComplementAssembler;
        procedure :: getInterfaceSolutionCalculator;
        procedure :: getBoundarySolutionCalculator;        
    end type DDDMS_SolverController
    
    ! DDDMS_SolverController user-defined data type constructor definition    
    interface DDDMS_SolverController
        module procedure NewDDDMS_SolverController;       ! Add constructor to DDDMS_SolverController generic interface
    end interface DDDMS_SolverController

contains
    
    ! DDDMS_SolverController user-defined data type constructor implementation    
    type(DDDMS_SolverController) function NewDDDMS_SolverController(self, solverTypeInput)
    
        implicit none;
            
        ! Declarative zone        
        class(DDDMS_SolverController) :: self;
        integer                       :: solverTypeInput;
        
        ! Body of the program        
        NewDDDMS_SolverController%solverType = solverTypeInput;
        print*, "Faccia di merda!";
        
    end function NewDDDMS_SolverController
   
    ! getDataLoader procedure implementation    
     function getDataLoader(self, pathFileName) result(dataLoader)
    
        implicit none;        
        
        ! Declarative zone        
        class(DDDMS_SolverController) :: self;
        class(DDDMS_DataLoaderFromFile), allocatable         :: dataLoader;
        integer                       :: status;
        character(len = *)            :: pathFileName;        
        
        ! Body of the program
        allocate( dataLoader, STAT = status, SOURCE = DDDMS_DataLoaderFromFile( pathFileName ) );        
        if (status /= 0) then
            print*, "Failed allocation of DDM loader!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if  
        
    end function getDataLoader

    ! getLocalSchurComplementCalculator procedure implementation    
    function getLocalSchurComplementCalculator(self, HPC_Arch) result(LSC_Calculator)
    
        implicit none;        
        
        ! Declarative zone        
        class(DDDMS_SolverController) :: self;
        class(*), allocatable         :: LSC_Calculator;
        integer(INT32)                :: status;
        integer(INT32)                :: HPC_Arch;
        
        ! Body of the program                      
        ! Shared memory hardware architecture        
!        if( HPC_Arch == 1 )then
            allocate( LSC_Calculator, STAT = status, SOURCE = DDDMS_SMLocalSchurComplementCalculator( ) );
            if (status /= 0) then
                print*, "Failed allocation of DDM LSC_Calculator!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if
                
        ! Distributed memory hardware architecture        
!        elseif( HPC_Arch == 2 )then
            
!        end if
                
    end function getLocalSchurComplementCalculator

    ! getGlobalSchurComplementAssembler procedure implementation    
    function getGlobalSchurComplementAssembler(self) result(GSC_Assembler)
    
        implicit none;
        
        !********************************************
        !DECLARATIVE ZONE
        !********************************************
        class(DDDMS_SolverController)                            :: self;
        class(DDDMS_GlobalSchurComplementAssembler), allocatable :: GSC_Assembler;
        integer(INT32)                                           :: status;
        
        !********************************************
        !BODY OF THE PROGRAM
        !********************************************  
        allocate( GSC_Assembler, STAT = status, SOURCE = DDDMS_GlobalSchurComplementAssembler( ) );
        if (status /= 0) then
            print*, "Failed allocation of DDM GSC_Assembler!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
                
    end function getGlobalSchurComplementAssembler
    
    
    ! getInterfaceSolutionCalculator PROCEDURE IMPLEMENTATION    
    function getInterfaceSolutionCalculator(self, HPC_Arch) result(IS_Calculator)
    
        implicit none;        
        
        ! Declarative zone        
        class(DDDMS_SolverController) :: self;
        class(*), allocatable         :: IS_Calculator;
        integer(INT32)                :: status;
        integer(INT32)                :: HPC_Arch;
        
        
        ! Body of the program            
        
        !Shared memory hardware architecture
        if( HPC_Arch == 1 .or. HPC_Arch == 2)then
            allocate( IS_Calculator, STAT = status, SOURCE = DDDMS_SMInterfaceSolutionCalculator( ) );
            if (status /= 0) then
                print*, "Failed allocation of DDM Interface Calculator!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if
            
        
        ! Distributed memory hardware architecture        
!        elseif( HPC_Arch == 2 )then
            
        end if
                
    end function getInterfaceSolutionCalculator
    
    ! getBoundarySolutionCalculator procedure implementation    
    function getBoundarySolutionCalculator(self, HPC_Arch) result(BS_Calculator)
    
        implicit none;
        
        
        ! Declarative zone        
        class(DDDMS_SolverController) :: self;
        class(*), allocatable         :: BS_Calculator;
        integer(INT32)                :: status;
        integer(INT32)                :: HPC_Arch;
        
        
        ! Body of the program               
        ! shared memory hardware architecture        
        if( HPC_Arch == 1 .or. HPC_Arch == 2 )then
            allocate( BS_Calculator, STAT = status, SOURCE = DDDMS_SMBoundarySolutionCalculator( ) );
            if (status /= 0) then
                print*, "Failed allocation of DDM Boundary Calculator!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if
            
        ! Distributed memory hardware architecture        
!        elseif( HPC_Arch == 2 )then
            
        end if
                
    end function getBoundarySolutionCalculator

end module SolverController_Module