!! Fortran code to calculate the flux
!! code will take in the occupation xyz from RMC
!! code will calculate the clusters in the lattice
!! code will perform cluster perturbations in the lattice and find flux
!! NEVER accept any change to the lattice configuration
!! final: code will provide a flux vector (each cluster will give a flux count)
!! run code as ./LF <RMC_xyz_file> <lattice_flux_input_file>

program Lattice_Flux
	!! declaration of variables
	character(len = 50) :: get_arg, lf_inp
	integer :: row,col,i,j,nA,c,n_seed
	integer, dimension(0:11) :: clusters
	integer, dimension(1:8) :: dt_seed
	real, dimension(0:11) :: wA
	real, dimension(1:11) :: flux
	! clusters array =  [nA,1nn,2nn,3nn,4nn,Ltrip,Itrip,L2,L3,quad,Q2,Q3]
	integer, dimension(:), allocatable :: seeds
	integer, dimension(:,:), allocatable :: lattice, lattice_cl
	real :: rn,Kb,T,beta,u_sys,beta_delmu, avg_cov
	character(len = 50) :: wdata,rmc_xyz
	character(len = 100) :: full_line
	real :: u_new,beta_del_u,bias,beta_u0,beta_u1
	integer,dimension(0:11) :: clusters0,clusters1
	integer :: r0,r1,c0,c1,m,nSRO
	! m is the number of mirrored rows/columns
	
	!! Random seeding implementation
	call random_seed(size = n_seed)
	allocate( seeds(1:n_seed) )
	do i = 1,n_seed
		seeds(i) = 1396578+i
	end do
	call random_seed(get=seeds)
	call date_and_time(values = dt_seed)
	seeds(n_seed)=dt_seed(8); seeds(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
	call random_seed(put=seeds)
	deallocate(seeds)
	
	!!input filename
	call get_command_argument(1,get_arg)
	rmc_xyz = trim(get_arg)
	call get_command_argument(2,get_arg)
	lf_inp = trim(get_arg)
	open(unit = 30,file = lf_inp,status = "old")
	read(30,"(A)") wdata
	read(30,*) row
	read(30,*) col
	close(30)
	open(unit = 30,file = "../newton_input.txt",status = "old")
	read(30,*) nSRO
	close(30)
	
	!! Total number of lattice sites = row*col
	c = 4							! Coordination number		
	Kb = 8.617333262e-5
	T = 1/Kb
	open(unit = 3,file = wdata,status = "old")
	read(3,*) 						!! Used to jump comment line in data file
	read(3,*) bias					!! Bias term of CEM
	read(3,*) wA(0)					!! Singlet - binding energy
	read(3,*) wA(1)					!! 1-NN interaction
	read(3,*) wA(2)					!! 2-NN interaction
	read(3,*) wA(3)					!! 3-NN interaction
	read(3,*) wA(4)					!! 4-NN interaction
	read(3,*) wA(5)					!! L-triplet interaction
	read(3,*) wA(6)					!! I-triplet interaction
	read(3,*) wA(7)					!! L2 interaction
	read(3,*) wA(8)					!! L3 interaction
	read(3,*) wA(9)					!! Quadruplet interaction
	read(3,*) wA(10)				!! Q2 interaction
	read(3,*) wA(11)				!! Q3 interaction
	close(3)
	beta = 1/Kb/T
	wA = beta*wA
	m = 7
	allocate( lattice(0:row-1,0:col-1) )
	allocate( lattice_cl(0:row+2*m-1,0:col+2*m-1) )
	
	flux = 0
	call read_xyz
!~ 	call show_lattice
!~ 	r0 = m
!~ 	r1 = row+m-1
!~ 	c0 = m
!~ 	c1 = col+m-1
	
	call find_flux
!~ 	open(unit = 10,file = "fast_flux_results.txt",status="old",position="append")
!~ 	write(10,*) clusters(1)/105625.0
	write(6,*) flux(:)
!~ 	close(10)
	
	contains
	subroutine system_energy() 						! calculates the energy of the system (non-dimensionalized)
		implicit none
		integer :: ii
		u_new = bias
		do ii = 0,11
			u_new = u_new + wA(ii)*clusters(ii)
		end do
	end subroutine
	
	subroutine find_flux()							! find the flux by attempting different hops
		implicit none
		integer :: ii,jj, atom_count
		flux = 0
		atom_count = 0
		do ii = 0,row-1
			do jj = 0,col-1
				if(lattice(ii,jj) == 1) then
!~ 					print *, "Flux for atom count:", atom_count
!~ 					print *, "Position:",ii,jj
					atom_count = atom_count + 1
					r0 = ii+m-4
					r1 = ii+m+4
					c0 = jj+m-4
					c1 = jj+m+4
					call find_clusters
					clusters0 = clusters
					call system_energy
					beta_u0 = u_new
					lattice(ii,jj) = 0
					! 1-NN hops
					if(ii < row-1 .AND. lattice(ii+1,jj) == 0) then
						lattice(ii+1,jj) = 1
						call hop
						lattice(ii+1,jj) = 0
					end if
					if(ii > 0 .AND. lattice(ii-1,jj) == 0) then
						lattice(ii-1,jj) = 1
						call hop
						lattice(ii-1,jj) = 0
					end if
					if(jj < col-1 .AND. lattice(ii,jj+1) == 0) then
						lattice(ii,jj+1) = 1
						call hop
						lattice(ii,jj+1) = 0
					end if
					if(jj > 0 .AND. lattice(ii,jj-1) == 0) then
						lattice(ii,jj-1) = 1
						call hop
						lattice(ii,jj-1) = 0
					end if
					lattice(ii,jj) = 1
				end if
			end do
		end do
						
	end subroutine
	
	subroutine hop()								! add the flux of a specific hop
		implicit none
		integer :: f
		real :: factor
		call find_clusters
		call system_energy
		clusters1 = clusters
!~ 		print *, clusters0(:)
!~ 		print *, clusters1(:)
		beta_u1 = u_new
		beta_del_u = beta_u1 - beta_u0
		factor = (min(1.0,exp(-beta_del_u)))
		do f = 1,11
			flux(f) = flux(f) + (clusters1(f)-clusters0(f))*factor
		end do
		
	end subroutine
	
	subroutine lattice_cl_init()					! initialize lattice_Cl
		implicit none
		integer :: ii
		lattice_cl(m:row+m-1,m:col+m-1) = lattice
		do ii = 1,m
			lattice_cl(m-ii,m-ii+1:) = lattice_cl(m+row-ii,m-ii+1:)
			lattice_cl(m+row+ii-1,m-ii+1:) = lattice_cl(m+ii-1,m-ii+1:)
			lattice_cl(m-ii+1:,m+col+ii-1) = lattice_cl(m-ii+1:,m+ii-1)
			lattice_cl(m-ii+1:,m-ii) = lattice_cl(m-ii+1:,m+col-ii)
			
			lattice_cl(m-ii,m-ii) = lattice_cl(m-ii,m+col-ii)
			lattice_cl(m-ii,m+col+ii-1) = lattice_cl(m-ii,m+ii-1)
			lattice_cl(m+row+ii-1,m-ii) = lattice_cl(m+ii-1,m-ii)
			lattice_cl(m+row+ii-1,m+col+ii-1) = lattice_cl(m+row+ii-1,m+ii-1)
		end do
	end subroutine
	
	subroutine clusters_info()						! print cluster information for current lattice
		implicit none
		print *
		print *, "Clusters information:-"
		print "(A20,A2,I20)", "Singlets",": ",		 clusters(0)
		print "(A20,A2,I20)", "1-NN clusters", ": ", clusters(1)
		print "(A20,A2,I20)", "2-NN clusters", ": ", clusters(2)
		print "(A20,A2,I20)", "3-NN clusters", ": ", clusters(3)
		print "(A20,A2,I20)", "4-NN clusters", ": ", clusters(4)
		print "(A20,A2,I20)", "L1-triplet clusters", ": ", clusters(5)
		print "(A20,A2,I20)", "I-triplet clusters", ": ", clusters(6)
		print "(A20,A2,I20)", "L2-triplet clusters", ": ", clusters(7)
		print "(A20,A2,I20)", "L3-triplet clusters", ": ", clusters(8)
		print "(A20,A2,I20)", "Quadruplet-1 clusters", ": ", clusters(9)
		print "(A20,A2,I20)", "Quadruplet-2 clusters", ": ", clusters(10)
		print "(A20,A2,I20)", "Quadruplet-3 clusters", ": ", clusters(11)
	end subroutine
	
	subroutine read_xyz()							! read an xyz file (RMC output) and create lattice
		implicit none
		real :: ii,jj,kk
		integer :: l
		character :: atom_type
		
		open(unit = 4,file = rmc_xyz,status = "old")
		read (4,"(A)") full_line
		read (4,"(A)") full_line
		nA = 0
		lattice = 0
		do l = 1,row*col
			read(4,"(A)") full_line
			read(full_line,*) atom_type,ii,jj,kk
			if(atom_type == 'A') then
				lattice(nint(ii),nint(jj)) = 1
				nA = nA + 1
			else
				lattice(nint(ii),nint(jj)) = 0
			end if
		end do
		close(4)
!~ 		print *, "Lattice created from xyz file"
	end subroutine
	
	subroutine show_lattice()						! print thr current square lattice with atom ids
		implicit none
		integer :: ii
		do ii = 0,row-1
			print *, lattice(ii,:)
			print *
			print *
		end do
	end subroutine
	
	subroutine find_clusters()						! all cluster subroutines in one subroutine
		implicit none
		clusters = 0
		call lattice_cl_init
		clusters(0) = nA
		call cluster_1nn
		call cluster_2nn
		call cluster_3nn
		call cluster_4nn
		call cluster_Ltrip
		call cluster_Itrip
		call cluster_L2
		call cluster_L3
		call cluster_quad
		call cluster_Q2
		call cluster_Q3
	end subroutine
	
	subroutine cluster_1nn()						! count 1-NN clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii,jj+1) == 1) then
						clusters(1) = clusters(1) + 1
					end if
					if(lattice_cl(ii+1,jj) == 1) then
						clusters(1) = clusters(1) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_2nn()						! count 2-NN clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+1,jj+1) == 1) then
						clusters(2) = clusters(2) + 1
					end if
					if(lattice_cl(ii+1,jj-1) == 1) then
						clusters(2) = clusters(2) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_3nn()						! count 3-NN clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj) == 1) then
						clusters(3) = clusters(3) + 1
					end if
					if(lattice_cl(ii,jj+2) == 1) then
						clusters(3) = clusters(3) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_4nn()						! count 4-NN clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj+1) == 1) then
						clusters(4) = clusters(4) + 1
					end if
					if(lattice_cl(ii+1,jj+2) == 1) then
						clusters(4) = clusters(4) + 1
					end if
					if(lattice_cl(ii+2,jj-1) == 1) then
						clusters(4) = clusters(4) + 1
					end if
					if(lattice_cl(ii+1,jj-2) == 1) then
						clusters(4) = clusters(4) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_Ltrip()						! count L-triplet clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+1,jj+1) == 1 .AND. lattice_cl(ii-1,jj+1) == 1) then
						clusters(5) = clusters(5) + 1
					end if
					if(lattice_cl(ii+1,jj+1) == 1 .AND. lattice_cl(ii+1,jj-1) == 1) then
						clusters(5) = clusters(5) + 1
					end if
					if(lattice_cl(ii-1,jj+1) == 1 .AND. lattice_cl(ii-1,jj-1) == 1) then
						clusters(5) = clusters(5) + 1
					end if
					if(lattice_cl(ii+1,jj-1) == 1 .AND. lattice_cl(ii-1,jj-1) == 1) then
						clusters(5) = clusters(5) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_Itrip()						! count I-triplet clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+1,jj+1) == 1 .AND. lattice_cl(ii-1,jj-1) == 1) then
						clusters(6) = clusters(6) + 1
					end if
					if(lattice_cl(ii+1,jj-1) == 1 .AND. lattice_cl(ii-1,jj+1) == 1) then
						clusters(6) = clusters(6) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_L2()							! count L2 clusters
	implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj+1) == 1 .AND. lattice_cl(ii+1,jj+2) == 1) then
						clusters(7) = clusters(7) + 1
					end if
					if(lattice_cl(ii-2,jj+1) == 1 .AND. lattice_cl(ii-1,jj+2) == 1) then
						clusters(7) = clusters(7) + 1
					end if
					if(lattice_cl(ii+2,jj-1) == 1 .AND. lattice_cl(ii+1,jj-2) == 1) then
						clusters(7) = clusters(7) + 1
					end if
					if(lattice_cl(ii-2,jj-1) == 1 .AND. lattice_cl(ii-1,jj-2) == 1) then
						clusters(7) = clusters(7) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_L3()							! count L3 clusters
	implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj+1) == 1 .AND. lattice_cl(ii+2,jj-1) == 1) then
						clusters(8) = clusters(8) + 1
					end if
					if(lattice_cl(ii-2,jj+1) == 1 .AND. lattice_cl(ii-2,jj-1) == 1) then
						clusters(8) = clusters(8) + 1
					end if
					if(lattice_cl(ii+1,jj+2) == 1 .AND. lattice_cl(ii-1,jj+2) == 1) then
						clusters(8) = clusters(8) + 1
					end if
					if(lattice_cl(ii+1,jj-2) == 1 .AND. lattice_cl(ii-1,jj-2) == 1) then
						clusters(8) = clusters(8) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_quad()						! count quadruplet clusters
		implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+1,jj+1) == 1 .AND. lattice_cl(ii+1,jj-1) == 1 .AND. lattice_cl(ii+2,jj) == 1) then
						clusters(9) = clusters(9) + 1
					end if
				end if
			end do
		end do
	end subroutine
	
	subroutine cluster_Q2()							! count q2 clusters
	implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj+1) == 1 .AND. lattice_cl(ii+2,jj-1) == 1 .AND. lattice_cl(ii+3,jj) == 1) then
						clusters(11) = clusters(11) + 1
					end if
					if(lattice_cl(ii-2,jj+1) == 1 .AND. lattice_cl(ii-2,jj-1) == 1 .AND. lattice_cl(ii-3,jj) == 1) then
						clusters(11) = clusters(11) + 1
					end if
					if(lattice_cl(ii+1,jj+2) == 1 .AND. lattice_cl(ii-1,jj+2) == 1 .AND. lattice_cl(ii,jj+3) == 1) then
						clusters(11) = clusters(11) + 1
					end if
					if(lattice_cl(ii+1,jj-2) == 1 .AND. lattice_cl(ii-1,jj-2) == 1 .AND. lattice_cl(ii,jj-3) == 1) then
						clusters(11) = clusters(11) + 1
					end if
				end if
			end do
		end do
		
	end subroutine
	
	subroutine cluster_Q3()							! count q3 clusters
	implicit none
		integer :: ii,jj
		
		do ii = r0,r1
			do jj = c0,c1
				if(lattice_cl(ii,jj) == 1) then
					if(lattice_cl(ii+2,jj) == 1 .AND. lattice_cl(ii+2,jj+2) == 1 .AND. lattice_cl(ii,jj+2) == 1) then
						clusters(10) = clusters(10) + 1
					end if
				end if
			end do
		end do
		
	end subroutine
end program
