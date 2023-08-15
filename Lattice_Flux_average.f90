!! utility code to average the flux present in current_flux_<i>.txt

!! use as ./LF_averager <nSRO> <# of threads>

program Lattice_Flux_average

	implicit none
	character(len = 50) :: get_arg,filename
	integer :: i,nSRO,threads
	real,dimension(:),allocatable :: f,f_avg
	call get_command_argument(1, get_arg)
	open(unit = 50,file = get_arg,status = "old")
	read(50,*) nSRO
	close(50)
	nSRO=11 !! force all flux to be computed
	call get_command_argument(2, get_arg)
	read(get_arg,*) threads
	
	allocate( f(nSRO), f_avg(nSRO) )
	f_avg = 0
	do i = 0,threads-1
		write(filename,"(A,I0,A)") "current_flux_",i,".txt"
		open(unit = 50,file = filename,status = "old")
		read(50,*) f(:)
		f_avg = f_avg + f
		close(50)
	end do
	f_avg = f_avg/threads
	write(6,*) f_avg
	deallocate( f, f_avg )
end program
