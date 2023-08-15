!! version 3 is a copy of 2
!! version 3 has no tolerance criteria
!! program ends after all newton iterations are over

PROGRAM MAIN
!find the equilibrium SRO parameter based on flux calculated
! the goal is to find where the flux is zero (approximately)

!the user can specify nSRO parameters, the flux is evaluated at these 
!SRO parameters

!composition has to be provided in RMC
!! number of SRO parameters set as parameter apriori !!
!! user defined parameters to be taken from newton_input.txt (argument input)
!! input structure
!! <z initial>
!! tol (tolerance)
!! maxiter

    IMPLICIT NONE
    INTEGER :: nSRO,rmdir_flag
    DOUBLE PRECISION :: x, tol, dp
    double precision, dimension(:),allocatable :: z,f,p,fl,fr,f_all,dz
    double precision, dimension(:,:),allocatable :: g
    INTEGER :: iter !number of iterations required
    INTEGER :: maxiter !max number of iterations allowed
    INTEGER :: i, j, frustration_level, threads, sleep_time
    integer :: t1,t2,t3,k,con_flag
    integer,dimension(1:11) :: SRO_flag	!! determines which clusters are constrained
    LOGICAL :: IsPrint=.FALSE.
    character(len = 50) :: get_arg,filename,filename_J
    character(len = 50) :: check_file,out_file, newt_inp
    character(len = 50) :: check_file_J,out_file_J
    character(len = 100) :: command
    call get_command_argument(1,get_arg)
    read(get_arg,*) threads
    call get_command_argument(2,get_arg)
    newt_inp = trim(get_arg)
    open(unit = 40,file = newt_inp,status = "old")
    read(40,*) nSRO
    allocate( z(nSRO), dz(nSRO), f(nSRO), fl(nSRO), fr(nSRO), g(nSRO,nSRO), p(nSRO), f_all(1:11) )
    read(40,*) x
    read(40,*) z(:)
    read(40,*) tol
    read(40,*) maxiter
    read(40,*) rmdir_flag
    read(40,*) SRO_flag(:)
    !! 0 to not remove folders (only for 1 newton iteration!!)
    !! 1 to keep final files
    !! otherwise always delete 
    close(40)
    print *, "Inputs"
    print *, "Number of SRO parameters: ",nSRO
    print *, "Fraction of A atoms: ", x
    print *, "Initial SRO values: ", z(:)
    print *, "Tolerance of flux: ", tol
    print *, "Maximum iterations: ", maxiter
    print *, "Number of threads per RMC simulations:", threads
    print *, "Output flag:", rmdir_flag
    print *, "SRO identifiers:", SRO_flag
!~    tol=1.e2 !flux
    write(filename,"(A)") "RMC_flux"	!! filename starts with R - necessary condition!!
    write(check_file,"(A)") trim(filename)//"/job_is_complete"
    write(out_file,"(A)") trim(filename)//"/current_flux.txt"
    !! Begin calculations
    print *, "Begin calculations"
!~     open(unit = 20,file = "cal.txt",status="replace")
!~ 	write(20,*) trim(filename), " 0.5 0.5", z(:)
!~ 	close(20)
!~     f = flux()
!~     write(command,"(A,A)") "./cleanup.sh ",filename 
!~     call system(command)
!~     g = flux_Jacobian()
!~     print *, g(1,:)
!~     print *, g(2,:)
!~     print *,"Done"
!~    stop
!~    CALL Newton(z,tol,maxiter,iter,errorstatus)
	con_flag = 0	!! Convergence flag
	do i = 1,maxiter
		call system("echo `date`")
		call newton_iteration
		if(i > 1 .AND. ALL(ABS(f)<tol) .AND. con_flag == 0) then
			con_flag = 1
			print *, "Converged"
		end if
		dz = z
		z = z - 0.5*matmul(inv(g),f)
		print *, "New z:", z(:)
		do j = 1,nSRO
			if(z(j) > 1) then
				z(j) = 1
			end if
			if(z(j) < 0) then
				z(j) = 0
			end if
		end do
		dz = dz - z
		if(i > 1 .AND. ALL(ABS(dz)<0.001)) then
			print *, "No change"
!~ 			call system("echo `date`")
!~ 			deallocate( dz,z,f,p,g,fl,fr,f_all )
!~ 			stop
		end if
		if(rmdir_flag /= 0) then
			call file_removal
		end if
	end do
	print *, "End of calculations"
	print *, "Final z:", z(:)
	call system("echo `date`")
	deallocate( dz,z,f,p,g,fl,fr,f_all )
    contains
    subroutine newton_iteration
		implicit none
		integer :: ii,jj
		logical :: rmc_status
		!! Begin by setting up the cal.txt file for a single newton iteration
		open(unit = 20,file = "cal.txt",status="replace")
		write(20,*) trim(filename),x,(1-x), z(:)
		close(20)
		open(unit = 20,file = "cal.txt",status = "old",position = "append")
		do ii = 1,nSRO
			p = z
			p(ii) = z(ii) + min(0.01,0.1*z(ii))
			write(filename_J,"(A,I0)") trim(filename)//"_",(2*ii-1)
			write(20,*) trim(filename_J),x,(1-x),p(:)
			p(ii) = z(ii) - min(0.01,0.1*z(ii))
			write(filename_J,"(A,I0)") trim(filename)//"_",(2*ii)
			write(20,*) trim(filename_J),x,(1-x),p(:)
		end do
		close(20)
		!! cal.txt is properly set up
		!! Begin RMC
		write(command,"(A,I0)") "./flux_master.sh cal.txt ",threads
		call system(command)
		print *, "Running all RMC simulations in parallel"
		!! RMC_flux is queued first - technically it should finish 1st
		!! Others should finish within a lag of about 5-10 seconds
		sleep_time = 60
		frustration_level = 0
		rmc_status = .FALSE.
		do while(.NOT. rmc_status)
			call sleep(sleep_time)
			frustration_level = frustration_level + 1
!~ 			print "(A,I0,A)", "Wait time: ", frustration_level, " minutes"
!~ 			if(frustration_level == 2*nSRO) then
!~ 				print *, "Patiently waiting"
!~ 			elseif(frustration_level == 5*nSRO) then
!~ 				print *, "Losing patience"
!~ 			elseif(frustration_level == 11*nSRO) then
!~ 				print *, "Please do better than this!!"
!~ 			end if
			inquire(file = check_file, exist = rmc_status)
		end do
		print "(A,I0,A)", "Total Wait time: ", frustration_level, " minutes"
		open(unit = 30,file = out_file,status = "old")
		read(30,*) f_all(:)
		close(30)
		k = 1
		do jj = 1,11
			if(SRO_flag(jj) == 1) then
				f(k) = f_all(jj)
				k = k + 1
			end if
		end do
		print *, "Flux :", f_all(:)
		print *, "Flux in analysis:", f(:)
		print *, "FLux computation complete"
		!! setup for the RMC simulations for the Jacobian matrix
		!! this is to be done in series
		do ii = 1,nSRO
			!! setup check_file an out_file TWICE for each SRO
			rmc_status = .FALSE.
			write(check_file_J,"(A,I0,A)") trim(filename)//"_",(2*ii-1),"/job_is_complete"
			do while(.NOT. rmc_status)
				call sleep(5)
				print *, "Waiting for fr"
				inquire(file = check_file_J, exist = rmc_status)
			end do
			write(out_file_J,"(A,I0,A)") trim(filename)//"_",(2*ii-1),"/current_flux.txt"
			open(unit = 30,file = out_file_J,status = "old")
			read(30,*) f_all(:)
			k = 1
			do jj = 1,11
				if(SRO_flag(jj) == 1) then
					fr(k) = f_all(jj)
					k = k + 1
				end if
			end do
			print *, "Fluxr :", f_all(:)
			print *, "Fluxr in analysis:", fr(:)
			close(30)
			rmc_status = .FALSE.
			write(check_file_J,"(A,I0,A)") trim(filename)//"_",(2*ii),"/job_is_complete"
			do while(.NOT. rmc_status)
				call sleep(5)
				print *, "Waiting for fl"
				inquire(file = check_file_J, exist = rmc_status)
			end do
			write(out_file_J,"(A,I0,A)") trim(filename)//"_",(2*ii),"/current_flux.txt"
			open(unit = 30,file = out_file_J,status = "old")
			read(30,*) f_all(:)
			k = 1
			do jj = 1,11
				if(SRO_flag(jj) == 1) then
					fl(k) = f_all(jj)
					k = k + 1
				end if
			end do
			close(30)
			print *, "Fluxl :", f_all(:)
			print *, "Fluxl in analysis:", fl(:)
			g(:,ii) = (fr-fl)/min(0.02,0.2*z(ii))
		end do
		print *, "Jacobian computation complete"
		
	end subroutine
	
	subroutine file_removal
		implicit none
		integer :: ii
		!! remove all the RMC simulation folders (keep the parent directory clean)
		write(command,"(A,A)") "./cleanup.sh ",trim(filename) 
		call system(command)
		do ii = 1,nSRO
			write(filename_J,"(A,I0)") trim(filename)//"_",(2*ii-1)
			write(command,"(A,A)") "./cleanup.sh ", trim(filename_J)
			call system(command)
			write(filename_J,"(A,I0)") trim(filename)//"_",(2*ii)
			write(command,"(A,A)") "./cleanup.sh ", trim(filename_J)
			call system(command)
		end do
		
	end subroutine
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
		DOUBLE PRECISION, dimension(:,:), intent(in) :: A
		DOUBLE PRECISION, dimension(size(A,1),size(A,2)) :: Ainv

		DOUBLE PRECISION, dimension(size(A,1)) :: work  ! work array for LAPACK
		integer, dimension(size(A,1)) :: ipiv   ! pivot indices
		integer :: n, info

		!External procedures defined in LAPACK
		external DGETRF
		external DGETRI

		! Store A in Ainv to prevent it from being overwritten by LAPACK
		Ainv = A
		n = size(A,1)

		! DGETRF computes an LU factorization of a general M-by-N matrix A
		! using partial pivoting with row interchanges.
		call DGETRF(n, n, Ainv, n, ipiv, info)

		if (info /= 0) then
			stop 'Matrix is numerically singular!'
		end if

		! DGETRI computes the inverse of a matrix using the LU factorization
		! computed by DGETRF.
		call DGETRI(n, Ainv, n, ipiv, work, n, info)

		if (info /= 0) then
			stop 'Matrix inversion failed!'
		end if
   end function
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    SUBROUTINE PrintMatrix(A)
       IMPLICIT NONE
       DOUBLE PRECISION, DIMENSION(:,:) :: A
       INTEGER :: j
      
       WRITE(6,*) '-----'
       DO j=1,SIZE(A,1) !print row
          WRITE(6,*) A(j,:)
       END DO
     END SUBROUTINE
     !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END PROGRAM
