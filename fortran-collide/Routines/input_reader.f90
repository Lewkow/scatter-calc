
SUBROUTINE collision_input_reader
  USE physics_constants, ONLY : TOAMU, TOEV, TODEG
  USE collision_inputs
	USE current_energy, ONLY : Channel
  USE mpi_info

	!-------------------------
	! Internal	
	!-------------------------

	INTEGER						:: i	! counter
	
	REAL(KIND=8)			:: DE	! Differential Energy

	WRITE(*,*)
	WRITE(*,*) "        $@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$"
	WRITE(*,*) "                    ASTROSCATT"
	WRITE(*,*) "        $@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$@$"

	OPEN(11, FILE="../Inputs/Keys.in", STATUS="old", ACTION="read")

	READ(11,*)								!1 
	READ(11,*)								!2
	READ(11,*)								!3
	READ(11,*)				  			!4	
	READ(11,*)								!5
	READ(11,*)								!6
	READ(11,*)								!7
	READ(11,*)								!8
	READ(11,*)								!9
	READ(11,*) Proj						!10
	READ(11,*) Targ						!11
	READ(11,*) IDNT						!12
	READ(11,*) FRAME					!13
	READ(11,*) CALC_TYPE			!14
	READ(11,*)								!15
	READ(11,*)								!16
	READ(11,*)								!17
	READ(11,*)								!18
	READ(11,*) TCS						!19
	READ(11,*) NUM_TCS			
	READ(11,*) DFCS						
	READ(11,*) DCS		  			!20
	READ(11,*) DC3D						!21
	READ(11,*) HS_TEST				!22
	READ(11,*) PROB_DEN				!23
	READ(11,*) PROB_DEN_3D		!24
	READ(11,*) PROB_DEN_FIT
	READ(11,*) AVE_EN_LOSS		
	READ(11,*) PHASE_WRITE		!25
	READ(11,*) PHASE_TIMING		!26
	READ(11,*) POT_TEST				!27
	READ(11,*)              	!28
	READ(11,*)								!29
	READ(11,*)								!30
	READ(11,*) POT_I					!31
	READ(11,*) POT_F					!32
	READ(11,*) POT_N					!33
	READ(11,*)								!34
	READ(11,*)								!35
	READ(11,*)								!36
	READ(11,*) ThI		 				!37
	READ(11,*) ThF						!38
	READ(11,*) NTh						!39
	READ(11,*)								!40
	READ(11,*) 								!41
	READ(11,*) 								!42
	READ(11,*) NumEngy				!43
	READ(11,*)								!44
	READ(11,*)								!45
	READ(11,*)								!46
	READ(11,*) EngyI					!47
	READ(11,*) EngyF					!48
	READ(11,*) EngyN					!49
	READ(11,*)								!50
	READ(11,*)								!51
	READ(11,*)								!52

	IF (NumEngy /= 0) THEN
		ALLOCATE(Engy(NumEngy))
		DO i=1,NumEngy
			READ(11,*) Engy(i)
		END DO
	ELSE
		NumEngy = EngyN
		ALLOCATE(Engy(NumEngy))
		DE = (EngyF - EngyI)/DBLE(NumEngy-1)
		DO i=1,NumEngy
			Engy(i) = EngyI + (i-1)*DE		
		END DO
	END IF

	CLOSE(11)

	!! Set Channel parameter if the collision being 
	!! calculated has multiple potential channels

	IF (Proj .EQ. 'H' .AND. Targ .EQ. 'H') THEN
		Channel = 'S'
	END IF

	IF ( (Proj .EQ. 'He4' .AND. Targ .EQ. 'O') .OR. (Proj .EQ. 'O' .AND. Targ .EQ. 'He4') ) THEN
		Channel = 'S'
	END IF

	CALL depth_finder
	CALL mass_finder

  M1    	= M1*TOAMU
  M2    	= M2*TOAMU
  MU    	= M1*M2/(M1+M2)
	ThI 		= ThI/TODEG
	ThF 		= ThF/TODEG
	Engy(:) = Engy(:)/TOEV

	66 FORMAT(A,I5)
	67 FORMAT(A,ES8.2)
 
  WRITE(*,*) 
	WRITE(*,*) "--------------------------------------------------"
	WRITE(*,66) "MPI tasks being used   = ", numrank	
	WRITE(*,66) "Number of Energies     = ", NumEngy	
	WRITE(*,*) "--------------------------------------------------"
  WRITE(*,67) "Projectile   (amu) M1  = ", M1/TOAMU
  WRITE(*,67) "Target       (amu) M2  = ", M2/TOAMU
  WRITE(*,67) "Reduced Mass (amu) MU  = ", MU/TOAMU
	WRITE(*,67) "Well Depth   (eV) Depth= ", Depth*TOEV	
	WRITE(*,*) "--------------------------------------------------"
  WRITE(*,*) 

	IF (TCS == 1) THEN
		WRITE(*,*) "--------------------------------------------------"
	  WRITE(*,*) "%      Total Cross Section Being Computed        %"
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,*)
	END IF
	
	IF (DCS == 1) THEN
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,*) "%   Differential Cross Section Being Computed    %"
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,'(A,F6.2)') "Initial Angle CM(deg) ThI= ", ThI*TODEG
		WRITE(*,'(A,F6.2)') "Final Angle   CM(deg) ThF= ", ThF*TODEG 
		WRITE(*,'(A,I6)') "Number of Angles      NTh= ", NTh
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,*)
	END IF 
	
	IF (HS_TEST == 1) THEN
		WRITE(*,*) "--------------------------------------------------"
	  WRITE(*,*) "%   Hard Sphere Accuracy Test Being Performed    %"
		WRITE(*,*) "--------------------------------------------------"
		WRITE(*,*)
	END IF

	WRITE(*,*) "--------------------------------------------------"
	IF ( (CALC_TYPE == 0) .OR. (CALC_TYPE == 99) ) THEN
		WRITE(*,*) "Partial Wave Calculation Being Used"
	END IF
	IF ( (CALC_TYPE == 1) .OR. (CALC_TYPE == 99) ) THEN
		WRITE(*,*) "Eikonal Approximation Calculation Being Used"
	END IF

	IF (FRAME == 0) THEN
		WRITE(*,*) "FRAME    -> Center of Mass"
	ELSE
		WRITE(*,*) "FRAME    -> Laberatory"
	END IF

	IF (IDNT == 0) THEN
		WRITE(*,*) "Symmetry -> Non-identical Particles"
	ELSE IF (IDNT == 1) THEN
		WRITE(*,*) "Symmetry -> Fermions (odd spin)"
	ELSE
		WRITE(*,*) "Symmetry -> Bosons (even spin)"
	END IF
	WRITE(*,*) "--------------------------------------------------"
  WRITE(*,*) 
  
END SUBROUTINE collision_input_reader

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE collision_input_bcast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Root process calls collision input reader to 
! parse the input keys file. The inputs are then
! packaged into a buffer, dependent on data type, 
! and is then broadcast to the rest of the MPI
! processes being used. 
!
!	The buffers must be changed as the input file and
! its arguments are modified for future tests or
! modules. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE mpi_info  
  USE collision_inputs

	IMPLICIT NONE

  INCLUDE 'mpif.h'

	INTEGER, PARAMETER									:: real_size = 10
	INTEGER, PARAMETER									:: int_size  = 19
 
  REAL(KIND=8), DIMENSION(real_size)	:: buff_real
  INTEGER, DIMENSION(int_size)				:: buff_int

  IF (rank == 0) THEN
    CALL collision_input_reader

    buff_real(1)  = M1
    buff_real(2)  = M2
    buff_real(3)  = MU
    buff_real(4)  = ThI
    buff_real(5)  = ThF
		buff_real(6)  = POT_I 
		buff_real(7)  = POT_F 
		buff_real(8)  = EngyI
		buff_real(9)  = EngyF
		buff_real(10) = Depth
 
    buff_int(1)  = TCS
		buff_int(2)	 = DFCS
    buff_int(3)  = DCS
    buff_int(4)  = FRAME
    buff_int(5)  = IDNT
    buff_int(6)  = NTh
		buff_int(7)  = DC3D
		buff_int(8)  = HS_TEST
		buff_int(9)  = PROB_DEN
		buff_int(10) = PHASE_WRITE
		buff_int(11) = PHASE_TIMING
		buff_int(12) = POT_TEST
		buff_int(13) = POT_N
		buff_int(14) = EngyN
 		buff_int(15) = PROB_DEN_3D
		buff_int(16) = CALC_TYPE
		buff_int(17) = AVE_EN_LOSS 
		buff_int(18) = NUM_TCS
		buff_int(19) = PROB_DEN_FIT
	END IF !rank=0
  
	!---------------------------
	! Broadcast buffer arrays	
	!---------------------------

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( buff_real, real_size, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( buff_int,  int_size,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Proj, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Targ, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
	!-----------------------------------------
	! Broadcast size of Engy array to allocate
	!-----------------------------------------

	CALL MPI_BCAST( NumEngy ,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

	!----------------------------------	
	! Allocate Engy array if not root
	! Assign module inputs from buffer
	!----------------------------------	

	IF (rank .NE. 0) THEN
		ALLOCATE(Engy(NumEngy))
	  M1    				= buff_real(1)
    M2    				= buff_real(2)
    MU    				= buff_real(3)
    ThI   				= buff_real(4)
    ThF   				= buff_real(5)
		POT_I 				= buff_real(6)
		POT_F 				= buff_real(7)
		EngyI 				= buff_real(8)
		EngyF 				= buff_real(9)
		Depth					= buff_real(10)

    TCS   				= buff_int(1)
		DFCS					= buff_int(2)
    DCS   				= buff_int(3)
    FRAME 				= buff_int(4)
    IDNT  				= buff_int(5)
    NTh   				= buff_int(6)
		DC3D					= buff_int(7)
		HS_TEST 			= buff_int(8)
		PROB_DEN 			= buff_int(9)
		PHASE_WRITE 	= buff_int(10)
		PHASE_TIMING 	= buff_int(11)
		POT_TEST 		 	= buff_int(12)
		POT_N					= buff_int(13)
		EngyN					= buff_int(14)	
		PROB_DEN_3D		= buff_int(15)
		CALC_TYPE			= buff_int(16)
		AVE_EN_LOSS		= buff_int(17)
		NUM_TCS				= buff_int(18)
		PROB_DEN_FIT	= buff_int(19)
	END IF 

	!--------------------------------
	! Broadcast Engy array from root
	!--------------------------------

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Engy, NumEngy, MPI_REAL8,  0, MPI_COMM_WORLD, ierr ) 
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE collision_input_bcast

