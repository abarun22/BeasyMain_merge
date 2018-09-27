module GlobalSchurComplementAssembler_Module
    use, intrinsic :: iso_fortran_env;
    use DataStructure_Module;
    use Utility_Module;                                        !==> TO DELETE?
    implicit none;
    private;
    public DDDMS_GlobalSchurComplementAssembler;

    ! DDDMS_GlobalSchurComplementAssembler definition        
    ! GlobalSchurComplementAssembler
    type :: DDDMS_GlobalSchurComplementAssembler
        
    contains
        procedure :: Assembly;
        
    end type DDDMS_GlobalSchurComplementAssembler    
    
    ! DDDMS_GlobalSchurComplementAssembler user-defined data type constructor definition    
    interface DDDMS_GlobalSchurComplementAssembler
        module procedure NewDDDMS_GlobalSchurComplementAssembler;       ! add constructor to DDDMS_DATALOADERFROMFILE generic interface
    end interface DDDMS_GlobalSchurComplementAssembler

contains
    
    ! DDDMS_GlobalSchurComplementAssembler user-defined data type constructor implementation    
    type(DDDMS_GlobalSchurComplementAssembler) function NewDDDMS_GlobalSchurComplementAssembler( )
    
        implicit none;
        
        !********************************************
        !DECLARATIVE ZONE
        !********************************************
        
        !********************************************
        !BODY OF THE PROGRAM
        !******************************************** 
        
    end function NewDDDMS_GlobalSchurComplementAssembler ! end of function NewDDDMS_GlobalSchurComplementAssembler    
    
    ! Assembly procedure implementation    
    subroutine Assembly(self, inputData, inputParamsData)

        
        ! Declarative zone        
        class(DDDMS_GlobalSchurComplementAssembler) :: self;
        class(DDDMS_DataContainer)                  :: inputData;
        class(DDDMS_InputParams)                    :: inputParamsData;
        integer(INT32)                              :: numOfZones;
        integer(INT64)                              :: i;
        integer(INT64)                              :: j;
        integer(INT64)                              :: k;
        integer(INT32)                              :: status;
        integer(INT64)                              :: RowCounter;
        integer(INT64)                              :: ColCounter;        
        integer(INT64)                              :: LocRowX;
        integer(INT64)                              :: LocRowY;
        integer(INT64)                              :: LocRowZ;
        integer(INT64)                              :: GloRowX;
        integer(INT64)                              :: GloRowY;
        integer(INT64)                              :: GloRowZ;
        integer(INT64)                              :: LocColX;
        integer(INT64)                              :: LocColY;
        integer(INT64)                              :: LocColZ;
        integer(INT64)                              :: GloColX;
        integer(INT64)                              :: GloColY;
        integer(INT64)                              :: GloColZ;
        character(len = 4)                          :: seqstring;        
        
! APB
        integer(INT64)                              :: ifnode;
        integer(INT64)                              :: idim;
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
        
        !Allocate and initialise matrix "gGlob" for the current zone        
        allocate( inputData%gGlob( idim*size( inputData%IndependenInterfaceNodesArray, 1 ) ) , STAT = status );
        if(status /= 0) then
            print*, "Failed allocation of gGlob!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            inputData%gGlob(:) = 0.0; 
        end if
                
        ! Allocate and initialise matrix "SGlob" for the current zone        
        allocate( inputData%SGlob( idim*size( inputData%IndependenInterfaceNodesArray, 1 ), &
                                   idim*size( inputData%IndependenInterfaceNodesArray, 1 ) ), STAT = status );
        if(status /= 0) then
            print*, "Failed allocation of SGlob!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            inputData%SGlob(:, :) = 0.0; 
        end if
        
        numOfZones = size( inputData%ZonesData, 1 );
        
        do k = 1,numOfZones
            RowCounter = 0;             
            
            do i = 1, size( inputData%IndependenInterfaceNodesArray(:, k) )
            
                if( inputData%IndependenInterfaceNodesArray(i, k) /= 0 ) then                                
                    
                    ! Index for local rows                    
                    LocRowX = int( idim*RowCounter + 1 );
                    if (idim .ne. 1) then
                        LocRowY = int( idim*RowCounter + 2 );
                        LocRowZ = int( idim*RowCounter + 3 );
                    endif

                    ! Index for global rows                                        
                    GloRowX = int( idim*( i - 1 ) + 1 );
                    if (idim .ne. 1) then
                        GloRowY = int( idim*( i - 1 ) + 2 );
                        GloRowZ = int( idim*( i - 1 ) + 3 );
                    endif
                    
                    inputData%gGlob( GloRowX ) = inputData%gGlob( GloRowX ) + inputData%ZonesData( k )%CB( LocRowX );  
                    if (idim .ne. 1) then
                        inputData%gGlob( GloRowY ) = inputData%gGlob( GloRowY ) + inputData%ZonesData( k )%CB( LocRowY );
                        inputData%gGlob( GloRowZ ) = inputData%gGlob( GloRowZ ) + inputData%ZonesData( k )%CB( LocRowZ );
                    endif                   
                    !=====================================================================        
                    ColCounter = 0;    
                    do j = 1, size( inputData%IndependenInterfaceNodesArray(:, k) )
                    
                        if( inputData%IndependenInterfaceNodesArray(j, k) /= 0 ) then                             
                            
                            ! Index for local columns                            
                            LocColX = int( idim*ColCounter + 1 );
                            if (idim .ne. 1) then
                                LocColY = int( idim*ColCounter + 2 );
                                LocColZ = int( idim*ColCounter + 3 );
                            endif                            
                            
                            ! Index for global columns                            
                            GloColX = int( idim*( j - 1 ) + 1 );
                            if (idim .ne. 1) then
                                GloColY = int( idim*( j - 1 ) + 2 );
                                GloColZ = int( idim*( j - 1 ) + 3 );                                    
                            endif                            
                            !=================================================================================================================================
                            inputData%SGlob( GloRowX, GloColX ) = inputData%SGlob( GloRowX, GloColX ) + inputData%ZonesData( k )%DL( LocRowX, LocColX );  
                            if (idim .ne. 1) then
                                inputData%SGlob( GloRowX, GloColY ) = inputData%SGlob( GloRowX, GloColY ) + inputData%ZonesData( k )%DL( LocRowX, LocColY );  
                                inputData%SGlob( GloRowX, GloColZ ) = inputData%SGlob( GloRowX, GloColZ ) + inputData%ZonesData( k )%DL( LocRowX, LocColZ );  
                            
                                inputData%SGlob( GloRowY, GloColX ) = inputData%SGlob( GloRowY, GloColX ) + inputData%ZonesData( k )%DL( LocRowY, LocColX );  
                                inputData%SGlob( GloRowY, GloColY ) = inputData%SGlob( GloRowY, GloColY ) + inputData%ZonesData( k )%DL( LocRowY, LocColY );                         
                                inputData%SGlob( GloRowY, GloColZ ) = inputData%SGlob( GloRowY, GloColZ ) + inputData%ZonesData( k )%DL( LocRowY, LocColZ ); 
                            
                                inputData%SGlob( GloRowZ, GloColX ) = inputData%SGlob( GloRowZ, GloColX ) + inputData%ZonesData( k )%DL( LocRowZ, LocColX );
                                inputData%SGlob( GloRowZ, GloColY ) = inputData%SGlob( GloRowZ, GloColY ) + inputData%ZonesData( k )%DL( LocRowZ, LocColY );  
                                inputData%SGlob( GloRowZ, GloColZ ) = inputData%SGlob( GloRowZ, GloColZ ) + inputData%ZonesData( k )%DL( LocRowZ, LocColZ );  
                            endif
                            ColCounter = ColCounter + 1;
                            !=================================================================================================================================
                        end if
                    end do
                
                    RowCounter = RowCounter + 1; 
                end if
            end do            
        end do ! End of loop statement on number of zones in the model        
        
        write(*,*)'Assembly of SGlob and gGlob done'
        
    end subroutine Assembly ! End of subroutine deldddms_dataloaderfromfile
    
end module GlobalSchurComplementAssembler_Module