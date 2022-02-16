!  code fragment to read parameters from a file. nah 28.8.08

!   included in alphabranch3

!  variables declared in module tree_pres3

! !  integer, parameter :: lng = selected_real_kind(12,99) in f90_tools.f90
!   real(lng), save :: Period, trm_rst, ff1, ff2, ff3, r_root, r_min, expo, asym, lrr, deltaR
!   integer, save   :: Maxgen, tmstps, maxStat, Ctmpstps
!
  character :: parafile*100 = 'parameters.dat'

  write (*,*) parafile
  open (1, file = parafile) ! Deals with empty spaces.

	read (1, '(///)') !  3 slashes skips 4 lines.
	read (1,*) Period ! List directed input.
	write (*,*) 'Period =', Period
	
	read (1,*) ; read (1,*) ! Skips 2 lines.
	read (1,*) trm_rst
	write (*,*) 'trm_rst =', trm_rst
	
	read (1,'(////)') ! Skips 5 lines & f1, f2, f3.
	read (1,*) ff1, ff2, ff3  ! List directed input.
	write (*,*) 'ff1 = ', ff1, 'ff2 = ', ff2, 'ff3 = ', ff3
	
	read (1,'(1x)') ; read (1,'(1x)') ! Skips 2 lines - another alternative
	read (1,*) r_root
	write (*,*) 'r_root =', r_root
	
	read (1,*) ; read (1,*)
	read (1,*) r_min
	write (*,*) 'r_min =', r_min
	
	read (1,*) ; read (1,*)
	read (1,*) alpha
	write (*,*) 'alpha =', alpha
	
	read (1,*) ; read (1,*)
	read (1,*) beta
	write (*,*) 'beta =', beta
	
	read (1,*) ; read (1,*)
	read (1,*) lrr
	write (*,*) 'lrr =', lrr
	
	read (1,*) ; read (1,*)
	read (1,*) Maxgen
	write (*,*) 'Maxgen =', Maxgen
	
	read (1,*) ; read (1,*)
	read (1,*) tmstps
	write (*,*) 'tmstps =', tmstps

	read (1,*) ; read (1,*)
	read (1,*) Ctmstps             !  used by Creaddata.C only
	write (*,*) 'Ctmstps =', Ctmstps

	read (1,*) ; read (1,*)
	read (1,*) maxStat
	write (*,*) 'maxStat =', maxStat
	
	read (1,*) ; read (1,*)
	read (1,*) deltaR
	write (*,*) 'deltaR =', deltaR
	
close (1)

