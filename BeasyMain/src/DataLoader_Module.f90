module DataLoader_Module
    use, intrinsic :: iso_fortran_env;
    use DataStructure_Module;
    use Utility_Module;                                      
    implicit none;
    private;
    public :: DDDMS_DataLoaderFromFile;

    ! Dataloader
    ! Dataloader class definition           
    type :: DDDMS_DataLoaderFromFile
        character(len = :), allocatable :: inputFilePathName;
        
    contains
        procedure :: CheckData;
        procedure :: LoadData;
        final     :: DelDDDMS_DataLoaderFromFile;
        
    end type DDDMS_DataLoaderFromFile

    
    ! DDDMS_DATALOADERFROMFILE user-defined data type constructor definition    
    interface DDDMS_DataLoaderFromFile
        module procedure NewDDDMS_DataLoaderFromFile;       ! add constructor to DDDMS_DataLoaderFromFile generic interface
    end interface DDDMS_DataLoaderFromFile

contains    
    ! DDDMS_DATALOADERFROMFILE user-defined data type constructor implementation    
    type(DDDMS_DataLoaderFromFile) function NewDDDMS_DataLoaderFromFile( inputFileName )
    
        implicit none;
        
        !********************************************
        !DECLARATIVE ZONE
        !********************************************
        character(len = *), intent(in) :: inputFileName;
        
        !********************************************
        !BODY OF THE PROGRAM
        !******************************************** 
        NewDDDMS_DataLoaderFromFile%inputFilePathName = inputFileName;
        
    end function NewDDDMS_DataLoaderFromFile ! End of function NewDDDMS_DataLoaderFromFile

    
    ! DDDMS_DataLoaderFromFile user-defined data type destructor implementation    
    elemental subroutine DelDDDMS_DataLoaderFromFile(self)

        !********************************************
        !DECLARATIVE ZONE
        !********************************************
        type(DDDMS_DataLoaderFromFile), intent(inout) :: self;
        integer                                       :: deallocStatusErr;
        
        !********************************************
        !BODY OF THE PROGRAM
        !******************************************** 
        deallocate(self%inputFilePathName, STAT = deallocStatusErr);
        
    end subroutine DelDDDMS_DataLoaderFromFile !END OF SUBROUTINE DELDDDMS_DATALOADERFROMFILE

    
    ! LoadData procedure implementation    
    subroutine LoadData(self, dataContainer, inputParamsData)
    
        implicit none;
        
        !********************************************
        !DECLARATIVE ZONE
        !********************************************
        class(DDDMS_DataLoaderFromFile) :: self;
        class(DDDMS_DataContainer)      :: dataContainer;
        class(DDDMS_InputParams)        :: inputParamsData;
        integer(INT64)                  :: countLinesFile;
        integer(INT32)                  :: status;
        integer(INT64),  allocatable    :: IndepInterNodesArrayTEMP(:, :);
        integer(INT64)                  :: temporaryIntVal(8); 
        real(REAL64)                    :: temporaryDblVal(4);
        real(REAL64),    allocatable    :: tempA_Matrix(:, :);
        real(REAL64),    allocatable    :: tempB_Matrix(:, :);
        real(REAL64),    allocatable    :: recordRow(:);
        integer(INT64)                  :: k;
        integer(INT64)                  :: j;
        integer(INT64)                  :: i;
        integer(INT64)                  :: countDuplicates;
        integer(INT32)                  :: IndepInterfNodesCount;
        integer(INT64)                  :: scaleIndexVal;
        integer(INT64)                  :: CurrentNodeId;
        integer(INT32)                  :: iszAOL(2);
        integer(INT64)                  :: fid;
        integer(INT32)                  :: fid_2;
        integer(INT64)                  :: numRowBlocksInZone;
        integer(INT64)                  :: recLen;
        integer(INT64)                  :: currRecord;
        integer(INT64)                  :: Res;              
        character(len = 200)            :: currentLine;
        character(len = 200)            :: FileNameAMatrix;
        character(len = 200)            :: FileNameBMatrix;
        character(len = 200)            :: FileNameCoordData;
        character(len = 200)            :: interfNodes;
        character(len = 4)              :: seqstring;                
! APB
        integer(INT64),  allocatable    :: bccode(:), nodeid(:);
        real(REAL64),    allocatable    :: bcval(:);
        integer(INT64)                  :: index1,index2;      
        integer(INT64)                  :: ibs,idim;
        integer(INT32)                  :: ZoneID1,ZoneID2,ZoneIDm,ifindex1,ifindex2
        real(REAL64)                    :: tmp
! APB        
!-------------------------------
! Initialisation
!-------------------------------        
        if (inputParamsData%analysType .eq. 1) then     ! Problem dimension
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif
        
        !********************************************
        !BODY OF THE PROGRAM
        !******************************************** 
        
        ! Count how many lines there are in the file about nodes on the interface in the model        
        countLinesFile = 0;
        open(NEWUNIT = fid, ACTION = 'read', FILE = inputParamsData%interfNodesFileName, STATUS = 'old');                  

        do
            read(fid, '(A)', IOSTAT = status) currentLine;
            if (status == IOSTAT_END) then
                exit;! END OF FILE REACHED
            endif  
            countLinesFile = countLinesFile + 1;
        end do
            
        close(fid);        
        
        ! Allocate temporary array to store data regarding node id on the interface before those duplicates are eliminated        
        allocate(IndepInterNodesArrayTEMP( ( countLinesFile - 1 ), 4 ), STAT = status);
        if (status /= 0) then
            print*, "Failed allocation of array IndepInterNodesArrayTEMP!";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            IndepInterNodesArrayTEMP(:, :) = -1; ! Initialise "IndepInterNodesArrayTEMP" array to value (-1).
        end if
        
        allocate(dataContainer%InterfaceNP(( countLinesFile - 1 ),2), STAT=status);
        if (status /= 0) then
            print*, "Failed allocation of array InterfaceNP!";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            dataContainer%InterfaceNP(:, :) = 0;
        end if        
        
        ! Load data regarding nodes index on interface in the model        
        write(*,*)'Reading interface definition file......'
        countLinesFile = 1;
        open(NEWUNIT = fid, ACTION = 'read', FILE = inputParamsData%interfNodesFileName, STATUS = 'old');
        read(fid, '(a)', IOSTAT = status) currentLine; ! Read the first line in the file containing an useless value
        do
            read(fid, *, IOSTAT = status) ( temporaryIntVal(j), j = 1, 4 );
            if (status == IOSTAT_END) then
                exit; ! end of file reached                
            end if              
                IndepInterNodesArrayTEMP(countLinesFile, 1) = temporaryIntVal(1);
                IndepInterNodesArrayTEMP(countLinesFile, 2) = temporaryIntVal(2);
                IndepInterNodesArrayTEMP(countLinesFile, 3) = temporaryIntVal(3);
                IndepInterNodesArrayTEMP(countLinesFile, 4) = temporaryIntVal(4);
                
            countLinesFile = countLinesFile + 1;
        end do
        close(fid); 
                
        ! Remove duplicates regarding nodes on interface in the model        
        countDuplicates = 0;
        do i = 1, size(IndepInterNodesArrayTEMP, 1)
            do j = (i + 1), size(IndepInterNodesArrayTEMP, 1)
                if( ( IndepInterNodesArrayTEMP(i, 1) == IndepInterNodesArrayTEMP(j, 1) ) .AND. &
                    ( IndepInterNodesArrayTEMP(i, 1) /= -1 ) ) then
                        
                    IndepInterNodesArrayTEMP(j, :) = -1;
                    countDuplicates = countDuplicates + 1;
                end if
            end do
        end do        
        
        ! Allocate sub data structure "dataContainer%IndependenInterfaceNodesArray"        
        allocate( dataContainer%IndependenInterfaceNodesArray( ( size( IndepInterNodesArrayTEMP, 1 ) - countDuplicates ), &
                                                                 inputParamsData%numberOfZones ), STAT = status );
        if (status /= 0) then
            print*, "Failed allocation of array dataContainer%IndependenInterfaceNodesArray!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            dataContainer%IndependenInterfaceNodesArray(:, :) = 0;    
        end if
        
        
        ! Copy independent interface nodes data from temporary data structure "IndepInterNodesArrayTEMP" TO
        ! Definitive data structure "dataContainer%IndependenInterfaceNodesArray"        
        countDuplicates = 1;
        do i = 1, size(IndepInterNodesArrayTEMP, 1)
            ifindex1=0
            ifindex2=0
            if( IndepInterNodesArrayTEMP(i, 1) /= -1 ) then                
                ZoneID1=IndepInterNodesArrayTEMP(i, 2)
                ZoneID2=IndepInterNodesArrayTEMP(i, 4)
                do j = 1,size(inputParamsData%zormap)                    
                    zoneIDm=inputParamsData%zormap(j)
                    if (ZoneID1.eq.zoneIDm)then
                        ifindex1=j
                    endif
                    if (ZoneID2.eq.zoneIDm)then
                        ifindex2=j
                    endif
                enddo
               
                dataContainer%IndependenInterfaceNodesArray( countDuplicates, ifindex1) = IndepInterNodesArrayTEMP(i, 1);
                dataContainer%IndependenInterfaceNodesArray( countDuplicates, ifindex2) = IndepInterNodesArrayTEMP(i, 3);
                
                dataContainer%InterfaceNP(i,1)=IndepInterNodesArrayTEMP(i, 1);
                dataContainer%InterfaceNP(i,2)=IndepInterNodesArrayTEMP(i, 3);
                countDuplicates = countDuplicates + 1;
            end if
        end do             

        ! Deallocate temporary data structure "IndepInterNodesArrayTEMP" after interface node data has been copied in
        ! Definitive data structure "dataContainer%IndependenInterfaceNodesArray"
        if( allocated( IndepInterNodesArrayTEMP ) ) then
            deallocate( IndepInterNodesArrayTEMP );
            if(status /= 0) then
                print*, "failed deallocation of array independeninterfacenodesarraytemp!";
                print*, "Error code: ", status;
                pause;
                stop;
            end if
        end if

        ! Allocate data structure "datacontainer%zonesdata" containing data for all zones        
        allocate( dataContainer%ZonesData( inputParamsData%numberOfZones ), STAT = status );                
        if(status /= 0) then
            print*, "Failed allocation of array dataContainer%ZonesData!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
        
        ! Loop over number of zones in the model to allocate data sub-structure for data container
        ! Next release===>the following loop could be replaced with a multithreaded block                 
        do i = 1, inputParamsData%numberOfZones
            IndepInterfNodesCount = count( dataContainer%IndependenInterfaceNodesArray(:, i) /= 0 ); ! number of independent nodes lying on the interface for the current zone
            write(*,*)'Zone #:',i            
            
            ! Retrive the name of the file containing a matric data            
            open(NEWUNIT = fid, ACTION = 'read', FILE = self%inputFilePathName, STATUS = 'old');       
            do j = 1, (6 + (i - 1)*4) 
                read(fid, '(A)', IOSTAT = status) FileNameAMatrix;      
            end do            
            
            ! Retrive the name of the file containing B matrix data            
            read(fid, '(A)', IOSTAT = status) FileNameBMatrix;     

            ! Retrive the name of the file containing nodes coordinates for the current zone            
            read(fid, '(A)', IOSTAT = status) FileNameCoordData;    

            close( fid );
            
            
            ! Open file containing nodes coordinates and boundary conditions for the current zone            
            write(*,*)'Reading nodal coordinates file....'
            
            open(NEWUNIT = fid, ACTION = 'read', FILE = FileNameCoordData, status = 'old'); 
            read(fid, '(I8)', IOSTAT = status) temporaryIntVal(1); ! Variable containing the number of nodes belonging to the current zone
            
            ! Allocate and initialise sub matrix "A00" for the current zone            
            allocate( dataContainer%ZonesData( i )%A00( idim*( temporaryIntVal(1) - IndepInterfNodesCount ), &
                                                        idim*( temporaryIntVal(1) - IndepInterfNodesCount ) ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix A00!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
               dataContainer%ZonesData( i )%A00 = 0.0; 
            end if           
            
            ! Allocate and initialise sub matrix "AOL" for the current zone
            allocate( dataContainer%ZonesData( i )%A0L( idim*( temporaryIntVal(1) - IndepInterfNodesCount ), &
                                                        idim*IndepInterfNodesCount ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix A0L!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%A0L = 0.0;    
            end if                        
            
            ! allocate and initialise sub matrix "AL0" for the current zone            
            allocate( dataContainer%ZonesData( i )%AL0( idim*IndepInterfNodesCount, &
                                                        idim*( temporaryIntVal(1) - IndepInterfNodesCount ) ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix AL0!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
               dataContainer%ZonesData( i )%AL0 = 0.0; 
            end if            

            ! Allocate and initialise sub matrix "ALL" for the current zone            
            allocate( dataContainer%ZonesData( i )%ALL( idim*IndepInterfNodesCount, &
                                                        idim*IndepInterfNodesCount ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix ALL!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
               dataContainer%ZonesData( i )%ALL = 0.0; 
            end if

            
            ! Allocate and initialise sub matrix "B00" for the current zone            
            allocate( dataContainer%ZonesData( i )%B00( idim*( temporaryIntVal(1) - IndepInterfNodesCount ), &
                                                        idim*( temporaryIntVal(1) - IndepInterfNodesCount ) ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix B00!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%B00 = 0.0;    
            end if

            ! Allocate and initialise sub matrix "B0L" for the current zone            
            allocate( dataContainer%ZonesData( i )%B0L( idim*( temporaryIntVal(1) - IndepInterfNodesCount ), &
                                                        idim*IndepInterfNodesCount ), STAT = status );
                                                        
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix B0L!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%B0L = 0.0;    
            end if
            
            ! Allocate and initialise sub matrix "BL0" for the current zone            
            allocate( dataContainer%ZonesData( i )%BL0( idim*IndepInterfNodesCount, &
                                                        idim*( temporaryIntVal(1) - IndepInterfNodesCount ) ), STAT = status );!--->REMEMBER TO DEALLOCATE DATA TYPE.
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix BL0!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%BL0 = 0.0;    
            end if
            
            ! Allocate and initialise sub matrix "BLL" for the current zone            
            allocate( dataContainer%ZonesData( i )%BLL( idim*IndepInterfNodesCount, &
                                                        idim*IndepInterfNodesCount ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix BLL!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%BLL = 0.0;    
            end if
            
            ! Allocate and initialise sub matrix "DL" for the current zone
            allocate( dataContainer%ZonesData( i )%DL( ( idim*IndepInterfNodesCount ), &
                                                       ( idim*IndepInterfNodesCount ) ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix DL!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%DL = 0.0;    
            end if
            
            ! Allocate and initialise sub vector "CB" for the current zone            
            allocate( dataContainer%ZonesData( i )%CB( idim*IndepInterfNodesCount ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-matrix DL!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%CB = 0.0;    
            end if
            
            ! Allocate and initialise sub vector "Y" for the current zone            
            allocate( dataContainer%ZonesData( i )%Y( idim*temporaryIntVal(1) ), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-vector Y!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%y = 0.0;    
            end if            
            
            ! Allocate and initialise sub vector "XB" for the current zone
            allocate( dataContainer%ZonesData( i )%XB( idim*( temporaryIntVal(1) - IndepInterfNodesCount ) ), STAT = status );!--->REMEMBER TO DEALLOCATE DATA TYPE.
            if(status /= 0) then
                print*, "Failed allocation of sub-vector XB!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%XB = 0.0;    
            end if
            
            ! Allocate and initialise sub vector "ULI" for the current zone               
            allocate( dataContainer%ZonesData( i )%ULI( idim*( IndepInterfNodesCount)), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-vector ULI!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%ULI = 0.0;    
            end if
                        
            ! Allocate and initialise sub vector "resvec" for the current zone            
            allocate( dataContainer%ZonesData( i )%resVec( idim*( temporaryIntVal(1))), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-vector resVec!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                dataContainer%ZonesData( i )%resVec = 0.0;    
            end if            
                       
            ! Allocate and initialise sub vector "coordsnodesinzone" to store nodes data for the current zone            
            allocate( dataContainer%ZonesData( i )%CoordsNodesInZone( temporaryIntVal(1), 4 ), STAT = status );!--->REMEMBER TO DEALLOCATE DATA TYPE.
            if(status /= 0) then
                print*, "Failed allocation of sub-vector CoordsNodesInZone!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if            
                       
            ! Load nodes coordinates data for the current zone in the sub-structure vector coordsnodesinzone
            countLinesFile = 1;
            do j = 1, temporaryIntVal(1) 
                read(fid, '(I8,3E20.8)', IOSTAT = status) temporaryIntVal(2), temporaryDblVal(1), &
                                                          temporaryDblVal(2), temporaryDblVal(3);  
                dataContainer%ZonesData( i )%CoordsNodesInZone( countLinesFile, 1) = temporaryIntVal(2);
                dataContainer%ZonesData( i )%CoordsNodesInZone( countLinesFile, 2) = temporaryDblVal(1);
                dataContainer%ZonesData( i )%CoordsNodesInZone( countLinesFile, 3) = temporaryDblVal(2);
                dataContainer%ZonesData( i )%CoordsNodesInZone( countLinesFile, 4) = temporaryDblVal(3);                                                                          
                countLinesFile = countLinesFile + 1;                
            end do ! end of loop over loading nodes coordinates data            
                        
            !  Allocate and initialise sub vector "bcstype" to store boundary conditions applied for the current zone            
            allocate( dataContainer%ZonesData( i )%BCsType( temporaryIntVal(1), idim), STAT = status );
            if(status /= 0) then
                print*, "Failed allocation of sub-vector BCsType!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if

            allocate(nodeid(1:idim), STAT = status )
            if(status /= 0) then
                print*, "Failed allocation of vector nodeid!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if
            
            allocate(bccode(1:idim), STAT = status )
            if(status /= 0) then
                print*, "Failed allocation of vector bccode!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if
            
            allocate(bcval(1:idim), STAT = status )
            if(status /= 0) then
                print*, "Failed allocation of vector bcval!";
                print*, "Errore code: ", status;
                pause;
                stop;
            end if                        
            
            ! Load boundary conditions data for the current zone            
            read(fid, '(I8)', IOSTAT = status) temporaryIntVal(2); !variable contains the number of boundary conditions applied
            scaleIndexVal = minval( dataContainer%ZonesData( i )%CoordsNodesInZone( :, 1) );
            
            ! Read the zonal node ID and boundary condition data from *.z00.txt file            
            do j = 1, temporaryIntVal(2), idim
                do k =1,idim
                    read(fid, '(2I8,E20.8)', IOSTAT = status) nodeid(k), bccode(k), bcval(k);
                enddo                
            
                if (status == IOSTAT_END) then
                    stop;
                endif  
            
                CurrentNodeId = int( nodeid(1) - scaleIndexVal + 1 );
                    
                do k =1,idim
                    if (bccode(k).eq.5 .or. bccode(k).eq.8 .or. bccode(k).eq.39 .or. bccode(k).eq.1 .or. &
                        bccode(k).eq.2 .or. bccode(k).eq.3 .or. bccode(k).eq.4) then                    
                        dataContainer%ZonesData(i)%BCsType( CurrentNodeId,1) = bccode(k);
                        if (inputParamsData%analysType .ne. 1) then
                            dataContainer%ZonesData(i)%y((idim*CurrentNodeId-2)) = bcval(k);
                        else
                            dataContainer%ZonesData(i)%y((CurrentNodeId)) = bcval(k);
                        endif
                    elseif (bccode(k).eq.6 .or. bccode(k).eq.9 .or. bccode(k).eq.40) then                    
                        dataContainer%ZonesData(i)%BCsType( CurrentNodeId,2) = bccode(k);
                        dataContainer%ZonesData(i)%y((idim*CurrentNodeId-1)) = bcval(k);
                    elseif (bccode(k).eq.7 .or. bccode(k).eq.10 .or. bccode(k).eq.41) then                    
                        dataContainer%ZonesData(i)%BCsType( CurrentNodeId,3) = bccode(k);
                        dataContainer%ZonesData(i)%y((idim*CurrentNodeId-0)) = bcval(k);
                    endif
                enddo                    
            enddo 
            
            
            ! Close file containing nodes coordinates and boundary conditions for the current zone
            close(fid);
            
            ! Calculate no. of blocks from which the matrices are to be read
            numRowBlocksInZone = ceiling( dble( idim*size( dataContainer%ZonesData( i )%CoordsNodesInZone( :, 1 ) , 1) )/dble( inputParamsData%blockSizeValue ) );
            
            
            ! Allocate temporary data structure tempa_matrix for the current zone            
            allocate( tempA_Matrix( ( inputParamsData%blockSizeValue*numRowBlocksInZone ), &
                                    ( inputParamsData%blockSizeValue*numRowBlocksInZone ) ), STAT = status );!--->REMEMBER TO DEALLOCATE TEMPA_MATRIX DATA TYPE.
            if(status /= 0) then
                print*, "Failed allocation of array tempA_Matrix!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                tempA_Matrix(:, :) = 0.0;
            end if

            ! Allocate temporary data structure tempb_matrix for the current zone
            allocate( tempB_Matrix( ( inputParamsData%blockSizeValue*numRowBlocksInZone ), &
                                    ( inputParamsData%blockSizeValue*numRowBlocksInZone ) ), STAT = status );!--->REMEMBER TO DEALLOCATE DATA TYPE TEMPB_MATRIX.
            if (status /= 0) then
                print*, "Failed allocation of array tempB_Matrix!";
                print*, "Errore code: ", status;
                pause;
                stop;
            else
                tempB_Matrix(:, :) = 0.0;
            end if

            ! Allocate data structure recordrow for the current zone
            allocate( recordRow( inputParamsData%blockSizeValue*inputParamsData%blockSizeValue + 2 ), STAT = status );
            if(status /= 0) then
                print*, "Allocation of array recordRow failed!.";
                print*, "Error code: ", status;
                pause;
                stop;
            else
                recordRow(:) = 0.0;
            end if
                        
            ! Load data concerning matrix a and b for the current zone from unformatted files            
            inquire( IOLENGTH = recLen ) recordRow;            

            open(NEWUNIT = fid,   ACTION = 'read', FILE = FileNameAMatrix, STATUS = 'old', ACCESS = 'direct', FORM = 'unformatted', RECL = recLen);
            open(NEWUNIT = fid_2, ACTION = 'read', FILE = FileNameBMatrix, STATUS = 'old', ACCESS = 'direct', FORM = 'unformatted', RECL = recLen);
            
            write(*,*)'Reading A and B matrices....'           
            currRecord = 0;
            ibs=inputParamsData%blockSizeValue;
            do k = 1,numRowBlocksInZone
                do j = 1,numRowBlocksInZone
                    currRecord = currRecord + 1;
                
                    read(fid, REC = currRecord, IOSTAT = status) recordRow;                    
                    if( status /= 0) then
                        pause;
                        exit;
                    endif                    
                    
!                    write(*,*)'Copying record for iteration: k=', k, 'and j=',j, 'Zone ID:',i
                                  
                   tempA_Matrix(((k-1)*ibs+1):(k*ibs),((j-1)*ibs+1):(j*ibs)) &
                                = reshape(recordRow(3:((ibs*ibs)+2)), &
                                  [ibs,ibs] );
                    
                     
                    read(fid_2, REC = currRecord, IOSTAT = status) recordRow;                    
                    if( status /= 0) then
                        exit;
                    endif
                    
                   tempB_Matrix(((k-1)*ibs+1):(k*ibs),((j-1)*ibs+1):(j*ibs)) &
                                = reshape(recordRow(3:((ibs*ibs)+2)), &
                                  [ibs,ibs] );                                  
                end do     
            end do

            close(fid);                 
            close(fid_2); 
            
            write(*,*)'Assigning A and B matrices data ....'
            write(*,*)

            ! Fill A00 sub-matrix for the current zone            
            write(*,*)'Creating A00 ....'
            index1=idim*(size(dataContainer%ZonesData( i )%CoordsNodesInZone( :, 1 ) , 1)- IndepInterfNodesCount);
            
            dataContainer%ZonesData( i )%A00(:, :) = tempA_Matrix( 1:index1, 1:index1);                                                                 
            
            ! Fill A0L sub-matrix for the current zone
            write(*,*)'Creating A0L ....'
			iszAOL=shape(dataContainer%ZonesData( i )%A0L)
	    	index2=index1+iszAOL(2)                                                             
            dataContainer%ZonesData( i )%A0L(:, :) = tempA_Matrix( 1:index1, ((index1)+1):index2);

            ! Fill AL0 sub-matrix for the current zone
            write(*,*)'Creating AL0 ....'
            dataContainer%ZonesData( i )%AL0(:, :) = tempA_Matrix( ((index1)+1):index2, 1:index1);
            
            ! Fill ALL sub-matrix for the current zone
            write(*,*)'Creating ALL ....'
            dataContainer%ZonesData( i )%ALL(:, :) = tempA_Matrix( ((index1)+1):index2, ((index1)+1):index2);

            ! Fill B00 sub-matrix for the current zone
            write(*,*)'Creating B00 ....'
            dataContainer%ZonesData( i )%B00(:, :) = tempB_Matrix( 1:index1, 1:index1);
            
            ! Fill B0L sub-matrix for the current zone
            write(*,*)'Creating B0L ....'
            dataContainer%ZonesData( i )%B0L(:, :) = tempB_Matrix( 1:index1, ((index1)+1):index2);
                        
            ! Fill BL0 sub-matrix for the current zone
            write(*,*)'Creating BL0 ....'
            dataContainer%ZonesData( i )%BL0(:, :) = tempB_Matrix(((index1)+1):index2, 1:index1);
            
            ! Fill BLL sub-matrix for the current zone
            write(*,*)'Creating BLL ....'
            dataContainer%ZonesData( i )%BLL(:, :) = tempB_Matrix( ((index1)+1):index2, ((index1)+1):index2);
            
            if( allocated(nodeid)) then
                deallocate(nodeid);
                if(status /= 0) then
                    print*, "Failed deallocation of array nodeid";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
            
            if( allocated(bccode)) then
                deallocate(bccode);
                if(status /= 0) then
                    print*, "Failed deallocation of array bccode";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
            
            if( allocated(bcval)) then
                deallocate(bcval);
                if(status /= 0) then
                    print*, "Failed deallocation of array bcval";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
            
            ! Deallocate temporary data structure recordrow for the current zone            
            if( allocated( recordRow ) ) then
                deallocate( recordRow );
                if(status /= 0) then
                    print*, "Failed deallocation of array recordRow";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if

            ! Deallocate temporary data structure tempa_matrix for the current zone
            if( allocated( tempA_Matrix ) ) then
                deallocate( tempA_Matrix );
                if(status /= 0) then
                    print*, "Failed deallocation of array tempA_Matrix";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if

            ! Deallocate temporary data structure tempb_matrix for the current zone            
            if( allocated( tempB_Matrix ) ) then
                deallocate( tempB_Matrix );
                if(status /= 0) then
                    print*, "Failed deallocation of array tempA_Matrix";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
            end if
        end do ! End of loop statement over number of zones
        
    end subroutine LoadData ! End of subroutine loaddata

    ! Checkdata procedure implementation    
    logical function CheckData(self)
    
        implicit none;        
        !********************************************
        ! Declarative zone
        !********************************************
        class(DDDMS_DataLoaderFromFile) :: self;
        character(len = 200)            :: currentLine;
        integer(INT32)                  :: indexPath;
        integer(INT32)                  :: countLinesFile;
        integer(INT32)                  :: status;
        integer(INT64)                  :: file_size;
        integer(INT64)                  :: NumOfZones;
        integer(INT32)                  :: i;
        integer(INT32)                  :: fid;
        logical                         :: file_exists; 
        
        !********************************************
        ! Body of the program
        !******************************************** 
        
        !********************************************************************************************************
        ! Check if the entered file name self%inputfilepathname exists and counts how many lines there are inside it
        !********************************************************************************************************
        inquire(FILE = self%inputFilePathName , EXIST = file_exists); !check if the input file "SELF%inputFilePathName" exists.
        inquire(FILE = self%inputFilePathName , SIZE  = file_size);   !check if the input file "SELF%inputFilePathName" is empty.
        
        indexPath = scan( self%inputFilePathName, "\", BACK = .TRUE.);
        
        if( file_exists .AND. (file_size /= 0) ) then
            CheckData = .TRUE.;
            open(NEWUNIT = fid, ACTION = 'read', FILE = self%inputFilePathName, STATUS = 'old');
        
            do
                read(fid, '(A)', IOSTAT = status) currentLine;
                if (status == IOSTAT_END) then
                    exit; ! end of file reached.abm
                end if
                
                if( ( index( trim(currentLine), '.txt') /= 0 ) .OR. &
                    ( index( trim(currentLine), '.c')   /= 0 ) .OR. & 
                    ( index( trim(currentLine), '.d')   /= 0 ) )then !the name of file was found
                    
                    inquire(FILE = self%inputFilePathName(1:indexPath)//currentLine , EXIST = file_exists); ! Check if the input file "SELF%inputFilePathName" exists.
                    inquire(FILE = self%inputFilePathName(1:indexPath)//currentLine , SIZE  = file_size);   ! Check if the input file "SELF%inputFilePathName" is empty.
                    
                    if( ( .not. file_exists) .AND. (file_size == 0) )then
                        print*, "File ", currentLine, " does not exist .";
                        print*, "Please, check input files path name folder ";
                        CheckData = .FALSE.;                        
                    end if
                end if
            end do  ! End of do statement to the end of the file self%inputFilePathName    
            
            close(fid);
        else
            write(OUTPUT_UNIT, *) "File name entered does not exist";
            write(OUTPUT_UNIT, *) "or its size is equal to zero!";
            pause;
            stop;
        end if
        
    end function CheckData ! End of function checkdata

end module DataLoader_Module