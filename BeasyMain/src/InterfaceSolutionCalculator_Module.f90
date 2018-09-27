module InterfaceSolutionCalculator_Module
    use, intrinsic :: iso_fortran_env;
    use DataStructure_Module;
    use Utility_Module;                                        
    implicit none;
    private;
    public :: DDDMS_InterfaceSolutionCalculator;
    public :: DDDMS_SMInterfaceSolutionCalculator;
    
    ! DDDMS_InterfaceSolutionCalculator definition    
    ! Abstract DDDMS_InterfaceSolutionCalculator user-defined data type definition    
    type, abstract :: DDDMS_InterfaceSolutionCalculator
        
    contains 
        procedure(CalculateABTRACT), deferred :: Dense_InterfaceSolution;
        procedure(CalculateABTRACT), deferred :: Sparse_InterfaceSolution        
    end type DDDMS_InterfaceSolutionCalculator
    
    ! DDDMS_InterfaceSolutionCalculator abstract user-defined data type procedure definition    
    abstract interface
        subroutine CalculateABTRACT(self, inputData, inputParamsData)
            import;
            class(DDDMS_InterfaceSolutionCalculator) :: self;
            class(DDDMS_DataContainer)               :: inputData;
            class(DDDMS_InputParams)                 :: inputParamsData;
        end subroutine CalculateABTRACT
    end interface
    
    ! Inherit  DDDMS_InterfaceSolutionCalculator user-defined data type
    type, extends(DDDMS_InterfaceSolutionCalculator) :: DDDMS_SMInterfaceSolutionCalculator        
    contains 
        procedure       :: Dense_InterfaceSolution => calculateDense;
        procedure       :: Sparse_InterfaceSolution => calculateSparse;
        
    end type DDDMS_SMInterfaceSolutionCalculator

    ! DDDMS_SMInterfaceSolutionCalculator user-defined data type constructor definition    
    interface DDDMS_SMInterfaceSolutionCalculator
        module procedure NewDDDMS_SMInterfaceSolutionCalculator; ! add constructor to DDDMS_SMInterfaceSolutionCalculator generic interface
    end interface DDDMS_SMInterfaceSolutionCalculator

contains

    ! DDDMS_SMInterfaceSolutionCalculator user-defined data type constructor implementation    
    type(DDDMS_SMInterfaceSolutionCalculator) function NewDDDMS_SMInterfaceSolutionCalculator(self)
    
        implicit none;
        
        !Declarative zone        
        class(DDDMS_SMInterfaceSolutionCalculator) :: self;
        
        ! Body of the program        
        
    end function NewDDDMS_SMInterfaceSolutionCalculator
    
    ! "calculateDense" procedure implementation    
    subroutine calculateDense(self, inputData, inputParamsData)
    
        implicit none;
        
        ! Declarative zone
        class(DDDMS_SMInterfaceSolutionCalculator) :: self;
        class(DDDMS_DataContainer)                 :: inputData;
        class(DDDMS_InputParams)                   :: inputParamsData;
        integer(INT64), allocatable                :: pivot(:);
        integer(INT32)                             :: info;
        integer(INT32)                             :: status;
! APB
        integer(INT64)                             :: idim;
        integer(INT64)                             :: ic1,ic2,elsc
! APB
         
!-------------------------------
! Initialisation
!-------------------------------        
        write(*,*)'Solution computed using MKL Dense solver (blas/lapack) routines....'
        
        if (inputParamsData%analysType .eq. 1) then     ! Problem dimension
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif

        ! Body of the program        
        ! Allocate the integer vector "pivot" containing the pivot indices
        allocate( pivot( size( inputData%SGlob, 1 ) ), STAT = status );
        if(status /= 0) then
            print*, "Allocation of array pivot failed!.";
            print*, "Error code: ", status;
            pause;
            stop;
        else
            pivot = 0.0;
        end if        
        
        ! Find the solution of the equations (sGlobal)A*x = (gGlobal)B using the lapack routine dgesv        
        write(*,*)'DGESV: Solving system of equations....'
       	call DGESV(size( inputData%SGlob, 1 ), & !is the order n of matrix A and the number of rows of matrix B.
                   1,                          & !is the number of right-hand sides; that is, the number of columns of matrix B.
                   inputData%SGlob,            & !S the general matrix A to be factored.
                   size( inputData%SGlob, 1 ), & !is the leading dimension of the array specified for A.
                   pivot,                      & 
                   inputData%gGlob,            & !Is the general matrix B, containing the nrhs right-hand sides of the system. 
                   size( inputData%gGlob, 1 ), & !Is the leading dimension of the array specified for B.
                   info);                           

        if( info > 0 ) then
            write(*,*)'The diagonal element of the triangular factor of inputData%SGlob,';
            write(*,*)'U(',info,',',info,') is zero, so that';
            write(*,*)'inputData%SGlob is singular; the solution could not be computed.';
            pause;
            stop;
        end if
        
        write(*,*)'DGESV: Solution phase completed'        
        
        ! Allocate and initialise vector INPUTDATA%U_INTERFACESOLUTION        
        allocate( inputData%U_InterfaceSolution( idim*size( inputData%IndependenInterfaceNodesArray, 1 ) ) , STAT = status );
        if(status /= 0) then
            print*, "Failed allocation of U_InterfaceSolution!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            inputData%U_InterfaceSolution = 0.0; 
        end if               
                
        ! Fill the vector gGlob using the solution of the equations A*X = B        
        inputData%U_InterfaceSolution = inputData%gGlob;               
        
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
    end subroutine calculateDense    
!
!****************************************************************************************
    subroutine calculateSparse(self, inputData, inputParamsData)
!****************************************************************************************
!   Objective: Usage of pardiso sparse routines to obtain the interface solution 
!              for system of equations  SGlob*UL=gGLOB
!****************************************************************************************
        implicit none;
        
!   Declarations
        class(DDDMS_SMInterfaceSolutionCalculator) :: self;
        class(DDDMS_DataContainer)                 :: inputData;
        class(DDDMS_InputParams)                   :: inputParamsData;
        integer(INT32)                             :: inz,i,j,k;
        integer(INT32)                             :: status;
        integer(INT64)                             :: idim,pt(64),szA(2),nrow,ncol,iparm(64),mtype;        
        integer(INT64),  allocatable               :: iia(:),ia(:),ja(:),perm(:);        
        integer(INT64)                             :: maxfct,mnum,phase,nrhs,msglvl,error,idum(1)
        real(REAL64),    allocatable               :: a(:);
        real(REAL64)                               :: ddum(1)
        integer(INT64)                             :: ic1,ic2,elsc
        
        write(*,*)'Solution computed using sparse MKL pardiso routines.....'
!
!-------------------------------
! Initialisation
!-------------------------------        
!
        inz=0
        maxfct=1
        mnum=1
        mtype=11
        nrhs=1
        msglvl=1
        error=0
        pt(1:64)=0
        
        if (inputParamsData%analysType .eq. 1) then     ! Problem dimension
            idim=inputParamsData%analysType
        elseif (inputParamsData%analysType .eq. 2) then
            idim=3
        endif
!
!---------------------------------------------------------------        
! Conversion of matrix SGlob in rectangular format to CSR format
!---------------------------------------------------------------        
!
        write(*,*)'Converting matrix SGlob from rectangular to CSR format...'
        szA=shape(inputData%SGlob)
        nrow=szA(1)
        ncol=szA(2)
        do i = 1,nrow
            do j = 1,ncol
                if (inputData%SGlob(i,j).ne.0) then
                    inz=inz+1
                endif
            enddo
        enddo
        
        allocate(a(1:inz), STAT = status )
        if(status /= 0) then
            print*, "Failed allocation of vector a!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
        
        allocate(ja(1:inz), STAT = status )
        if(status /= 0) then
            print*, "Failed allocation of vector ja!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
        
        allocate(iia(1:inz), STAT = status )
        if(status /= 0) then
            print*, "Failed allocation of vector iia!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
            
        allocate(ia(1:(nrow+1)), STAT = status )
        if(status /= 0) then
            print*, "Failed allocation of vector ia!";
            print*, "Errore code: ", status;
            pause;
            stop;
        end if
        
        allocate(perm(1:nrow), STAT = status )
        if(status /= 0) then
            print*, "Failed allocation of vector perm!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            perm(1:nrow)=0
        end if
        
        allocate( inputData%U_InterfaceSolution(1:nrow) , STAT = status );
        if(status /= 0) then
            print*, "Failed allocation of U_InterfaceSolution!";
            print*, "Errore code: ", status;
            pause;
            stop;
        else
            inputData%U_InterfaceSolution = 0.0; 
        end if               

! Form 'a' vector and the indices of SGlob for the non-zero terms
        inz=0
        do i = 1,nrow
            do j = 1,ncol
                if (inputData%SGlob(i,j).ne.0) then
                    inz=inz+1
                    a(inz)=inputData%SGlob(i,j)
                    iia(inz)=i
                    ja(inz)=j
                endif
            enddo
        enddo

! Form 'ia' vector with pointers to first non zero term in the column of each row
        do i = 1,nrow
            do j = 1,ncol
                if (inputData%SGlob(i,j).ne.0) then
                    exit
                endif                
            enddo
            do k = 1,inz
                if (i.eq.iia(k).and.j.eq.ja(k)) then
                    ia(i)=k
                    exit
                endif
            enddo
        enddo
        ia(nrow+1)=inz+1
        
        write(*,*)'Setting the working array iparm()'
! Initialize pardiso solver
        iparm(1) = 0 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(3) = 1 ! numbers of processors
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(7) = 0 ! not in use
        iparm(8) = 9 ! numbers of iterative refinement steps
        iparm(9) = 0 ! not in use
        iparm(10) = 13 ! handle non-symmetric matrices
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(12) = 0 ! not in use
        iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(15) = 0 ! not in use
        iparm(16) = 0 ! not in use
        iparm(17) = 0 ! not in use
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations

! Analysis, Numerical Factorization, Solution and Iterative refinement
        write(*,*)'Performing Analysis, Numerical Factorization, Solution and Iterative refinement.....'
        phase=13     
        CALL pardiso_64 (pt, maxfct, mnum, mtype, phase, nrow, a, ia, ja,  &
                         perm, nrhs, iparm, msglvl, inputData%gGlob, inputData%U_InterfaceSolution, error)
        
! Release internal memory        
        write(*,*)'Releasing internal memory acquired for computation.....'
        phase=-1
        CALL pardiso_64 (pt, maxfct, mnum, mtype, phase, nrow, ddum, idum, idum,  &
                         idum, nrhs, iparm, msglvl, ddum, ddum, error)
                         
        if( allocated(a)) then
                deallocate(a);
                if(status /= 0) then
                    print*, "Failed deallocation of array a";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
        end if
        if( allocated(iia)) then
                deallocate(iia);
                if(status /= 0) then
                    print*, "Failed deallocation of array iia";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
        end if
        if( allocated(ia)) then
                deallocate(ia);
                if(status /= 0) then
                    print*, "Failed deallocation of array ia";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
        end if
        if( allocated(ja)) then
                deallocate(ja);
                if(status /= 0) then
                    print*, "Failed deallocation of array ja";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
        end if
        if( allocated(perm)) then
                deallocate(perm);
                if(status /= 0) then
                    print*, "Failed deallocation of array perm";
                    print*, "Error code: ", status;
                    pause;
                    stop;
                end if
        end if        
        
    end subroutine calculateSparse
    
end module InterfaceSolutionCalculator_Module