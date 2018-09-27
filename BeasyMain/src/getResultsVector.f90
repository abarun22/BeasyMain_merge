module getResultsVector_Module
    use, intrinsic :: iso_fortran_env;                
    use DataStructure_Module
    implicit none;
    private;
    public :: DDDMS_getResultsVector;    
    
    type :: DDDMS_ResultsVar
        real(REAL64), allocatable :: ResultsVar(:)
    end type DDDMS_ResultsVar
    
    type :: DDDMS_getResultsVector
        type(DDDMS_ResultsVar), allocatable :: ZonalResults(:);
    contains 
        procedure :: getResultsVector;
    end type DDDMS_getResultsVector    
    
contains   
    
!
!*********************************************************************************************
    subroutine getResultsVector(self,resultsData, inputParamsData)
!*********************************************************************************************
!   Objective: Stack and output results vector in the following sequence
!              Zone 1: Boundary results
!                      Interface results
!              Zone 2: Boundary results
!                      Interface results
!                   .
!                   .
!                   .    
!                   .
!              Zone n: Boundary results
!                      Interface results    
!*********************************************************************************************
!   
#ifdef __WIN64
        use ifport
#endif
        use iso_fortran_env;
        implicit none;                
    
!   Declarations    
        class(DDDMS_getResultsVector)             :: self;
        class(DDDMS_DataContainer)                :: resultsData;
        class(DDDMS_InputParams)                  :: inputParamsData;
        character*100                             :: dirname,filename
        integer(INT32)                            :: i;
        integer(INT32)                            :: k;
        integer(INT32)                            :: iszXB;
        integer(INT32)                            :: iszUL;
        integer(INT32)                            :: iszresv;
        integer(INT32)                            :: j;        
        integer(INT32)                            :: nzones;
        integer(INT32)                            :: status,istat,nunit;
        integer(INT64)                            :: idim;
        
        if (inputParamsData%analysType .eq. 1) then
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif
    
! Output results vector in the file 'resultvector.txt' in the folder Output_file (in source directory)
        istat=getcwd(dirname)

#ifdef __WIN64
        dirname=trim(dirname)//'\Output_file'
        filename=trim(dirname)//'\resultvector.txt'
#elif __LINUX__
        dirname=trim(dirname)//'/Output_file'
        filename=trim(dirname)//'/resultvector.txt'
#endif

! open(NEWUNIT=nunit, FILE=filename, STATUS='UNKNOWN',ERR=10)
	OPEN(NEWUNIT=nunit, FILE=filename,  ACTION="WRITE", STATUS = "UNKNOWN",ERR=10)

        
! Write results for all zones in the output file
        nzones=size(resultsData%ZonesData);           
        write(*,*)'Stacking results in results vector.....'
        do i = 1,nzones            
            iszXB=size(resultsData%ZonesData(i)%XB(:))        
            iszUL=size(resultsData%ZonesData(i)%ULI(:))            
            iszresv=(iszXB+iszUL)
            allocate( self%ZonalResults(i)%ResultsVar(iszresv), STAT = status)            
            self%ZonalResults(i)%ResultsVar(1:iszXB)=resultsData%ZonesData(i)%XB(:)
            self%ZonalResults(i)%ResultsVar(iszXB+1:iszresv)=resultsData%ZonesData(i)%ULI(:)            
            if (idim.eq.1)then
                do k = 1,iszresv
                    write(nunit,*)self%ZonalResults(i)%ResultsVar(k)
                enddo
            else
#ifdef __WIN64	
                write(nunit,*)self%ZonalResults(i)%ResultsVar
#elif __LINUX__
				do k = 1,(iszresv/idim)
	           		write(nunit,*)self%ZonalResults(i)%ResultsVar(idim*(k-1)+1:k*idim)
		    	enddo
#endif
            endif            
        enddo   
        goto 20
        
 10     Continue
        write(OUTPUT_UNIT,*)'Error opening file......',filename

 20     continue
        write(OUTPUT_UNIT,*)'Results vector outputted in file:',filename
 
    end subroutine
!*********************************************************************************************
    
end module getResultsVector_Module