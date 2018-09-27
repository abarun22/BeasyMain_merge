module Utility_Module
    use iso_fortran_env;
    use DataStructure_Module;
    implicit none;

contains  
    
    ! Procedure to print out data contained in a matrix of real coefficients    
    subroutine Utility_OutputDataMatrix(fileName, dataToPrint)

        implicit none;    
        
        ! Declarative zone        
        character(len = *), intent(in) :: fileName;
        character(len = 4)             :: seqstring; 
        class(*), intent(in)           :: dataToPrint(:, :);
        integer(INT32)                 :: fid;
        integer(INT64)                 :: k;
        integer(INT64)                 :: j;
        integer(INT64)                 :: i;
        
        ! body of the program        
        select type(dataToPrint)        ! APB:  Polymorphism with type bound Procedures
            
        type is ( real(REAL64) )            
            open( newunit = fid, ACTION = 'write', FILE = fileName, status = 'unknown' );       
            do j = 1, size( dataToPrint, 1)
                do k = 1, size( dataToPrint, 2)
                    write(fid, '(E12.6,X)', ADVANCE = 'NO') dataToPrint(j, k);
                end do
                write(fid, *);
            end do                        
            close( fid );
            
        type is ( integer(INT32) )            
            open( newunit = fid, ACTION = 'write', FILE = fileName, status = 'unknown' );       
            do j = 1, size( dataToPrint, 1)
                do k = 1, size( dataToPrint, 2)
                    write(fid, '(I6,X)', ADVANCE = 'NO') dataToPrint(j, k);
                end do
                write(fid, *);
            end do                        
            close( fid );
            
        type is ( integer(INT64) )            
            open( newunit = fid, ACTION = 'write', FILE = fileName, status = 'unknown' );       
            do j = 1, size( dataToPrint, 1)
                do k = 1, size( dataToPrint, 2)
                    write(fid, '(I6,X)', ADVANCE = 'NO') dataToPrint(j, k);
                end do
                write(fid, *);
            end do                        
            close( fid );

            
        end select
                
    end subroutine Utility_OutputDataMatrix
    
    ! Procedure to print out data contained in a vector of real coefficients    
    subroutine Utility_OutputDataVector(fileName, dataToPrint)

        implicit none;
        !********************************************
        ! Declarative zone
        !********************************************
        character(len = *), intent(in) :: fileName;
        character(len = 4)             :: seqstring; 
        real(REAL64), intent(in)       :: dataToPrint(:);
        integer(INT32)                 :: fid;
        integer(INT64)                 :: k;
        integer(INT64)                 :: j;
        integer(INT64)                 :: i;
        
        !********************************************
        ! Body of the program
        !********************************************    
        open( newunit = fid, ACTION = 'write', FILE = fileName, status = 'unknown' );       
        do j = 1, size( dataToPrint, 1)
            write(fid, "(E12.6,X)") dataToPrint(j);
        end do                        
        close( fid );
        
    end subroutine Utility_OutputDataVector

    ! Procedure to collect data from command line arguments    
    subroutine Utility_CommandParameters(inputSolverPathName, &
                                         inputModelFileName,  &
                                         inputPathFileName,   &
                                         inputBlockFileName,  &
                                         batchFileCommands,   &
                                         HPC_Arch)
    
        implicit none;
        
        !Declarative zone        
        character(len = 200)              :: inputSolverPathName;      !==>TO DELETE?
        character(len = 200)              :: inputModelFileName;       !==>TO DELETE?
        character(len = 200)              :: inputPathFileName;        !==>TO DELETE?
        character(len = 200)              :: inputBlockFileName;       !==>TO DELETE?
        character(len = 200)              :: currentLine;              !==>TO DELETE?
        character(len = 100)              :: currCommandArg;           !==>TO DELETE?
        character(len = 100), allocatable :: batchFileCommands(:);
        integer(INT32)                    :: HPC_Arch;
        integer(INT32)                    :: k;
        integer(INT32)                    :: j;
        integer(INT32)                    :: i;
        integer(INT64)                    :: file_size;
        integer(INT64)                    :: count;
        integer(INT32)                    :: status;
        integer(INT32)                    :: indexPath;
        integer(INT32)                    :: fid;
        logical                           :: file_exists;        
        
        ! Body of the program        
        count = command_argument_count()
        if( count > 1 )then
            
            do i = 1,count
                call get_command_argument( i, currCommandArg);
                if( len_trim(currCommandArg) == 0 ) exit;
            
                if( index( trim(currCommandArg), '-pf=') /= 0 ) then !check for file containing input parameters
                    
                    indexPath = scan( currCommandArg, "=" );
            
                    ! Check if the entered file name "currcommandarg" exists                    
                    inquire( FILE = currCommandArg( ( indexPath + 1): ) , EXIST = file_exists );  !check if the input file currcommandarg exists.
                    inquire( FILE = currCommandArg( ( indexPath + 1): ) , SIZE  = file_size );    !check if the input file currcommandarg is emptY.
            
                    if( file_exists .AND. (file_size /= 0) ) then
                        open(NEWUNIT = fid, ACTION = 'read', FILE = currCommandArg( ( indexPath + 1): ), STATUS = 'old');
                        read(fid, '(A)', IOSTAT = status) inputSolverPathName;
                        read(fid, '(A)', IOSTAT = status) inputModelFileName;
                        read(fid, '(A)', IOSTAT = status) inputPathFileName;
                        read(fid, '(A)', IOSTAT = status) inputBlockFileName;
                        read(fid, '(I2)', IOSTAT = status) HPC_Arch;
                        close(fid);       
                    endif
              
                elseif( index( trim(currCommandArg), '-bf=') /= 0 ) then !check for file containing batch commands
                
                    indexPath = scan( currCommandArg, "=" );
                
                    ! Check if the entered file name "currcommandarg" exists and counts how many lines there are inside it                
                    inquire( FILE = currCommandArg( ( indexPath + 1): ) , EXIST = file_exists );  !check if the input file currcommandarg exists.
                    inquire( FILE = currCommandArg( ( indexPath + 1): ) , SIZE  = file_size );    !check if the input file currcommandarg is empty.
            
                    if( file_exists .AND. (file_size /= 0) ) then
                        open( NEWUNIT = fid, ACTION = 'read', FILE = currCommandArg( ( indexPath + 1): ), STATUS = 'old');
            
                        count = 0;
                        do
                            read(fid, '(A)', IOSTAT = status) currentLine;
                            if (status == IOSTAT_END) then
                                exit; ! end of file reached.abm
                            end if                
                        
                            count = count + 1;
                        end do  !end of do statement to the end of the file self%inputFilePathName    
            
                        close(fid); 
                        
                        ! Allocate data for data structure batchfilecommands                        
                        allocate( batchFileCommands(count), STAT = status );
                        if (status /= 0) then
                            print*, "Failed allocation of batchFileCommands array!";
                            print*, "Errore code: ", status;
                            pause;
                            stop;
                        end if  
            
                        batchFileCommands(:) = '*';
                        
                        ! Fill batchfilecommands with commands from input batch file                        
                        open( NEWUNIT = fid, ACTION = 'read', FILE = currCommandArg( ( indexPath + 1): ), STATUS = 'old');        
                        count = 1;
                        do
                            read(fid, '(A)', IOSTAT = status) currentLine;
                            if (status == IOSTAT_END) then
                                exit; ! End of file reached.abm
                            end if                
                        
                            batchFileCommands(count) = currentLine;
                            count = count + 1;
                        end do  ! end of do statement to the end of the file self%inputFilePathName    
            
                        close(fid); 
                     
                    endif !end of if statement on checking if input batch file exists and is not empty
                    
                endif !end of if statement about checking for file containing input parameters
            
            end do !end of loop statement over number of command line arguments
        
        else

            ! Body of the main program inputSolverPathName            
            
            ! Enter pathname containing the beasy solver            
            write(OUTPUT_UNIT, '("Enter pathname containing the beasy solver...:")', ADVANCE = "NO");
            read(OUTPUT_UNIT, *)inputSolverPathName;
            write(OUTPUT_UNIT, *)"Beasy solver pathname entered is...:", inputSolverPathName; 
            
            
            ! Enter file pathname containing all input files to load            
            write(OUTPUT_UNIT, '("Enter input file pathname containing model informations...:")', ADVANCE = "NO");
            read(OUTPUT_UNIT, *)inputModelFileName;
            write(OUTPUT_UNIT, *)"input file pathname entered is...:", inputModelFileName; 
            
            
            ! Enter file pathname containing all input files to load            
            write(OUTPUT_UNIT, '("Enter file pathname containing all input files to load...:")', ADVANCE = "NO");
            read(OUTPUT_UNIT, *)inputPathFileName;
            write(OUTPUT_UNIT, *)"folder pathname entered is...:", inputPathFileName; 
            
            
            ! Enter file pathname containing dump system size used            
            write(OUTPUT_UNIT, '("Enter file pathname containing dump system size used...:")', ADVANCE = "NO");
            read(OUTPUT_UNIT, *)inputBlockFileName;
            write(OUTPUT_UNIT, *)"folder pathname entered is...:", inputBlockFileName; 
            
            
            ! Select hpc architecture where ddm will be run            
            write(OUTPUT_UNIT, '("Select HPC architecture...:")');

            write(OUTPUT_UNIT, '("Shared Memory arch...............==> 1 ")');
            write(OUTPUT_UNIT, '("Distributed Memory arch..........==> 2 ")');
            read(OUTPUT_UNIT, *)HPC_Arch;
            write(OUTPUT_UNIT, *)"HPC arch...:", HPC_Arch;   
            
        endif ! end of if statement about input parameters
        
    end subroutine Utility_CommandParameters

    
    ! Procedure to collect data from command to performe boundary solution    
    function Utility_ZoneBoundarySolution( inputStringLine, inputParamsData ) result(zonesBoundSolArray)
    
        implicit none;        
        
        ! Declarative zone        
        character(len = *), intent(in) :: inputStringLine;
        character(len = :), allocatable:: param;             
        class(DDDMS_InputParams)       :: inputParamsData;
        integer(INT32), allocatable    :: zonesBoundSolArray(:);
        integer(INT32)                 :: indexPath;
        integer(INT32)                 :: indexPath2;
        integer(INT32)                 :: status;
        integer(INT32)                 :: i;
        integer(INT32)                 :: numericParam1;
        integer(INT32)                 :: numericParam2;   
        
        ! Body of the program        
        indexPath = scan( trim(inputStringLine), "," );

        param = inputStringLine( (indexPath + 1):len_trim(inputStringLine) );
        read(param(1:), '(I2)', IOSTAT = status) numericParam1;

        ! Option all zones: CALC_BS,ALL        
        if( param == "ALL" ) then
            
            allocate( zonesBoundSolArray( inputParamsData%numberOfZones ), STAT = status );
            if (status /= 0) then
                print*, "Failed allocation of DDM zonesBoundSolArray array!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                zonesBoundSolArray = 0;
            end if
            
            zonesBoundSolArray = [ (i, i = 1, inputParamsData%numberOfZones ) ];

        ! Option single zone: CALC_BS,N        
        elseif( ( numericParam1 >= 1 )                             .AND. &
                ( numericParam1 <= inputParamsData%numberOfZones ) .AND. &
                ( scan( param, ",") == 0) ) then
            
            allocate( zonesBoundSolArray( 1 ), STAT = status );
            if (status /= 0) then
                print*, "Failed allocation of DDM zonesBoundSolArray array!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                zonesBoundSolArray = 0;
            end if
            
            zonesBoundSolArray = [ numericParam1 ];

        ! Option range: CALC_BS,N-M -- (WITH N <= M)        
        elseif( scan( trim(param), "-" ) /= 0 ) then
                    
            indexPath2 = scan( trim(param), "-" );
            read( param( 1:(indexPath2 - 1) ),               '(I2)', IOSTAT = status ) numericParam1;
            read( param( (indexPath2 + 1):len_trim(param) ), '(I2)', IOSTAT = status ) numericParam2;
            
            if( (numericParam1 <= numericParam2) .AND. &
                (numericParam1 >= 1)             .AND. &
                (inputParamsData%numberOfZones >= numericParam2) ) then
                
                allocate( zonesBoundSolArray( numericParam2 - numericParam1 + 1 ), STAT = status );
                if (status /= 0) then
                    print*, "Failed allocation of DDM zonesBoundSolArray array!";
                    print*, "Errore code: ", status;
                    pause;
                    stop;
                else
                    zonesBoundSolArray = 0;
                end if
            
                zonesBoundSolArray = [ (i, i = numericParam1, numericParam2 ) ];
                
            else
                pause;
                stop;
                
            end if

        ! Option list: CALC_BS,N,M -- (WITH N <= M)        
         elseif( scan( trim(param), "," ) /= 0 ) then
                    
            indexPath2 = scan( trim(param), "," );
            read( param( 1:(indexPath2 - 1) ),               '(I2)', IOSTAT = status ) numericParam1;
            read( param( (indexPath2 + 1):len_trim(param) ), '(I2)', IOSTAT = status ) numericParam2;
            
            if( (numericParam1 <= numericParam2) .AND. &
                (numericParam1 >= 1)             .AND. &
                (inputParamsData%numberOfZones >= numericParam2) ) then
                
                allocate( zonesBoundSolArray( 2 ), STAT = status );
                if (status /= 0) then
                    print*, "Failed allocation of DDM zonesBoundSolArray array!";
                    print*, "Errore code: ", status;
                    pause;
                    stop;
                else
                    zonesBoundSolArray = 0;
                end if
            
                zonesBoundSolArray(1) = numericParam1;
                zonesBoundSolArray(2) = numericParam2;
                
            else
                pause;
                stop;
                    
            end if
                       
        else
            write(OUTPUT_UNIT, '("Parameter was not recognised!")');
            write(OUTPUT_UNIT, '("Please, check it!")');
                
        end if      
        
    end function Utility_ZoneBoundarySolution
 
    ! Procedure to convert a generic strint to upper case string    
    subroutine To_upper( str )
        character(*) :: str;
        integer      :: i;
 
        do i = 1, len( trim(str) )
            select case( str(i:i) )
                case("a":"z")
                str(i:i) = achar( iachar( str( i:i ) ) - 32 )
            end select
        end do 
                 
    end subroutine To_upper
    
    ! Procedure to print out nodes deformed coordinates data to paraview    
    subroutine Fortran2VTK( fileName, dataToPrint )
    
        implicit none;    
        
        ! Declarative zone        
        character(len = *), intent(in) :: fileName;
        integer(INT32)                 :: fid;
        real(REAL64),       intent(in) :: dataToPrint(:, :);
        integer(INT64)                 :: i;
        
        
        ! Body of the program        
        open( NEWUNIT = fid, ACTION = 'write', FILE = fileName, STATUS = 'unknown' );       
            write(fid, FMT = '("# vtk DataFile Version 2.0")');                         
            write(fid, FMT = '("InputFile4 Output")');
            write(fid, FMT = '("ASCII")');
            write(fid, FMT = '("DATASET UNSTRUCTURED_GRID")');
            write(fid, FMT = '("POINTS",I5,X,"float")')size(dataToPrint, 1);
            
            do i = 1,size(dataToPrint, 1)
                write(fid, FMT = '(3(E12.6,X))')dataToPrint(i, 2), dataToPrint(i, 3), dataToPrint(i, 4);    
            end do

            write(fid, FMT = '("CELLS",2(I5,X))')size(dataToPrint, 1), (2*size(dataToPrint, 1));
            do i = 1,size(dataToPrint, 1)
                write(fid, FMT = '(2(I5,X))')1, (i - 1);    
            end do

            write(fid, FMT = '("CELL_TYPES",I5,X)')size(dataToPrint, 1);
            do i = 1,size(dataToPrint, 1)
                write(fid, FMT = '(I5)')1;    
            end do
            
        close( fid );        
                 
    end subroutine Fortran2VTK

end module Utility_Module
    