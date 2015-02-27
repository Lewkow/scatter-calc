
SUBROUTINE ave_energy_loss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate the average percent energy loss
! per collision. 
!
! FORM: <dE>/E = (4*pi*m1*m2)/((m1+m2)^2*TC) *
!	int_{0}^{pi} sin(theta)*DCS*(1-cos(theta)) d theta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : M1, M2
	USE physics_constants, ONLY : PI, TOEV
	USE dcs, ONLY : Amp
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'	

	!--------------------------
	! Internal	
	!--------------------------

	REAL(KIND=8)										:: dE   	! Percent Energy Loss
	REAL(KIND=8)										:: TC 		! Total Cross Section
	REAL(KIND=8)										:: C1			! Const non-integrand term in return
	REAL(KIND=8)										:: C2			! Const integrand term in return
	REAL(KIND=8)										:: T			! Theta
	REAL(KIND=8)										:: DT			! dTheta
	REAL(KIND=8)										:: tstart	! dTheta
	REAL(KIND=8)										:: tend		! dTheta

	REAL(KIND=8),DIMENSION(numrank)	:: buff 	! Buffer to communicate

	INTEGER													:: i			! Counter
	INTEGER													:: N			! Number of point for integral
	INTEGER													:: chunk 	! Chunk size for each rank
	INTEGER													:: my_st 	! Starting counter for each rank
	INTEGER													:: my_fn 	! Ending counter for each rank

	N = 100000

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL total_cross_section_MPI(TC)

	C1 = 4.0*PI*M1*M2/(((M1+M2)**2.0)*TC)
	C2 = 0.0
	DT = PI/(N-1)	

	chunk = N/numrank
	my_st = rank*chunk
	my_fn = (rank+1)*chunk - 1

	IF (rank == (numrank-1)) THEN
		my_fn = N
	END IF

	IF (rank == 0) THEN
		tstart = MPI_WTIME()
	END IF	
	
	DO i=my_st,my_fn
		T = 0.0 + (i-1)*DT
		CALL pw_amp(T)
		C2 = C2 + SIN(T)*Amp*(1.0-COS(T))*DT
	END DO
	
	CALL MPI_GATHER(C2,1,MPI_REAL8,buff,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (rank == 0) THEN
		tend = MPI_WTIME()
		dE   = 0.0

		DO i=1,numrank
			dE = dE + buff(i)
		END DO
	
		dE = dE*C1	

		OPEN(UNIT=56, FILE="../Data/AveEngyRelax.dat", ACCESS="APPEND")
		WRITE(56,*) E*TOEV, dE
		CLOSE(56)

		WRITE(*,"(A,F5.2,A)") "Average Energy Loss Calculation Complete in ",tend-tstart, " sec"

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE ave_energy_loss

!################################################################
!################################################################

SUBROUTINE energy_loss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculate the energy loss per unit time when
! multiplied by density and velocity 
!
! FORM: <dE>/E = (4*pi*m1*m2)/(m1+m2)^2 *
!	int_{0}^{pi} sin(theta)*DCS*(1-cos(theta)) d theta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : M1, M2
	USE physics_constants, ONLY : PI, TOEV
	USE dcs, ONLY : Amp
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'	

	!--------------------------
	! Internal	
	!--------------------------

	REAL(KIND=8)										:: dE   	! Percent Energy Loss
	REAL(KIND=8)										:: C1			! Const non-integrand term in return
	REAL(KIND=8)										:: C2			! Const integrand term in return
	REAL(KIND=8)										:: T			! Theta
	REAL(KIND=8)										:: DT			! dTheta
	REAL(KIND=8)										:: tstart	! dTheta
	REAL(KIND=8)										:: tend		! dTheta

	REAL(KIND=8),DIMENSION(numrank)	:: buff 	! Buffer to communicate

	INTEGER													:: i			! Counter
	INTEGER													:: N			! Number of point for integral
	INTEGER													:: chunk 	! Chunk size for each rank
	INTEGER													:: my_st 	! Starting counter for each rank
	INTEGER													:: my_fn 	! Ending counter for each rank

	N = 100000

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	C1 = 4.0D0*PI*M1*M2/((M1+M2)**2.0)
	C2 = 0.0D0
	DT = PI/REAL(N-1)	

	chunk = N/numrank
	my_st = rank*chunk
	my_fn = (rank+1)*chunk - 1

	IF (rank == (numrank-1)) THEN
		my_fn = N
	END IF

	IF (rank == 0) THEN
		tstart = MPI_WTIME()
	END IF	
	
	DO i=my_st,my_fn
		T = REAL(i-1)*DT
		CALL pw_amp(T)
		C2 = C2 + SIN(T)*Amp*(1.0D0-COS(T))*DT
	END DO
	
	CALL MPI_GATHER(C2,1,MPI_REAL8,buff,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (rank == 0) THEN
		tend = MPI_WTIME()
		dE   = 0.0

		DO i=1,numrank
			dE = dE + buff(i)
		END DO
	
		dE = dE*C1	

		OPEN(UNIT=56, FILE="../Data/EngyRelax.dat", ACCESS="APPEND")
		WRITE(56,*) E*TOEV, dE
		CLOSE(56)

		WRITE(*,"(A,F5.2,A)") "Energy Loss Calculation Complete in ",tend-tstart, " sec"

	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE energy_loss

!################################################################
!################################################################

SUBROUTINE ave_percent_energy_loss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calculates the average percent energy loss for several angles btwn
! 0 and PI using: 
!
! dE(t) = int_{0}^{t} (1-gam)*rho*dt
! gam   = ( 1 + 2*p*cos(t) + p^2 )/( 1 + p )^2
! rho   = 2*pi*DCS/TCS
! p     = mp/mt
!
! t     = scattering angle in CM frame
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE current_energy, ONLY : E
	USE collision_inputs, ONLY : M1, M2
	USE physics_constants, ONLY : PI, TOEV
	USE dcs, ONLY : Amp

	IMPLICIT NONE

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: x,y,grand,gg,pd
	REAL(KIND=8)				:: full, t_i, t_f, dt, t, p, t_now, tot, gam, rho
	REAL(KIND=8)				:: E_cm, E_lab, v
	INTEGER							:: N_t, i, j	

	OPEN(UNIT=77,FILE="../Data/AvePercentEngyLoss.dat",ACCESS="APPEND")
	
	p   = M1/M2											! projectile-target mass ratio
	N_t = 10000											! number of angles to calculate for
	t_i = 0.0D0											! initial scattering angle CM frame [rad]
	t_f = PI												! final scattering angle CM frame [rad]
	dt  = (t_f - t_i)/REAL(N_t-1)		! angle differential [rad]

	ALLOCATE(x(N_t),y(N_t),grand(N_t),gg(N_t),pd(N_t))

	DO i=1,N_t
		t        = t_i + (i-1)*dt
		gg(i)    = 1.0D0-(1.0D0+2.0D0*p*COS(t)+p*p)/((1.0D0+p)*(1.0D0+p)) 
		CALL prob_density(t,rho) 
		pd(i)    = rho	
		grand(i) = gg(i)*pd(i)*dt
	END DO	

	DO i=1,N_t
		t 	= t_i + (i-1)*dt
		tot = 0.0D0
		DO j=1,i
			tot = tot + grand(j)
		END DO !j
		x(i) = t
		y(i) = tot
	END DO !i	

	full = y(N_t)
	DO i=1,N_t
		y(i) = y(i)/full
		WRITE(77,*) x(i)*180.0D0/PI, gg(i), pd(i), y(i)
	END DO
	
	DEALLOCATE(x,y,grand,gg,pd)
	CLOSE(77)

END SUBROUTINE ave_percent_energy_loss

