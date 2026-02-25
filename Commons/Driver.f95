program Driver
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File Driver.f95
!
! J.-M. Petit  Observatoire de Besac
! Version 1:  May 2013
!
! The purpose of this program is to compare models of the Kuiper Belt to
! reality by actually comparing what the surveys (say CFEPS or OSSOS)
! would have found if the model was a good representation of the real
! world, to the objects that were actually found.
!
! The workflow of the survey simulator driver is as follows:
!
!     Loop (some condition):
!         call GiMeObj(arg_list_1)
!         Check model ended:
!             set exit condition
!         call Detos1(arg_list_2)
!         Check detection and tracking:
!             store results
!
! The GiMeObj routine is in charge of providing a single new object at
! each call. The Detos1 routine determines is the proposed object would
! have been detected by the survey.
!
! The survey simulator expects orbital elements with respect to ecliptic
! reference frame.
!
! Logical unit numbers 7 to 19 are reserved to use by Driver.f and
! SurveySubs.f and should not be used by GiMeObj or any other routine
! that you may add to the driver. Please use logical unit numbers
! starting from 20.
!
! As currently written, the driver includes a file 'GiMeObj.f'
! containing the definition of the model. The correct way to use this
! feature is to have one's GiMeObj routine in a file <whatever.f> and
! create a symobolic link: 
!
!     ln -s <whatever.f> GiMeObj.f
!
! For more information, please refer to README files in the source directory.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  use gimeobjut

  implicit none

  integer, parameter :: screen = 6, keybd = 5, verbose = 9
  integer :: lun_h
  type(t_orb_m) :: o_m
! color array NEEDS to be length 10 or more!
  real (kind=8) :: h, epoch, hmax, rn_iter
  integer :: ierr, seed, n_iter, nmax, nchar, values(8), seed_mod
  character(80) :: det_outfile, line
  character(100) :: comments
  character(10) :: time
  character(8) :: date
  character(5) :: zone
  logical :: keep_going, rec

  lun_h = 10
  keep_going = .true.

! Get arguments
! Seed for random number generator
  read (keybd, *, err=9999) seed
  seed_mod = seed
! -maximum number of trials (<0) or calibrated number of objects
  read (keybd, *, err=9999) nmax
! Maximum value of H to draw
  read (keybd, *, err=9999) hmax
! Do we want to record the objects inside GiMeObj ?
  read (keybd, *, err=9999) rec
! Name for the model objects outfile
  read (keybd, '(a)', err=9999) det_outfile
  det_outfile = strip_comment(det_outfile)

!  print *, n_track_max

! Open output files and write header
  open (unit=lun_h, file=det_outfile, status='new', err=9500)

  call date_and_time(date, time, zone, values)
  write (lun_h, '(a17,a23,2x,a5)') '# Creation time: ', &
       date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'  &
       //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
  write (lun_h, '(''#'')')
  write (lun_h, '(''# Seed: '', i10)') seed
  write (lun_h, '(''#'')')
  write (lun_h, '(a,a)') &
       '#   a        e        i      Omega    omega      M', &
       '        H       epoch      comment '
  
! Initialize counters
  n_iter = 0
  rn_iter = 0.d0

! Open and read in object distribution
100 continue
      
!     Main loop: loop on objects
  if (keep_going) then

!        Select object.
     nchar = 0
     call GiMeObj (seed_mod, nmax, hmax, rec, &
          o_m, epoch, h, comments, nchar, ierr)

     if (ierr .eq. -10) then       ! Something wrong with this object, go to next one
        goto 100
     else if (ierr .eq. -20) then     ! Something very wrong happend, stop
        write (screen, '(a)') 'GiMeObj returned -20, stopping.'
        goto 2010
     else if (ierr .eq. 100) then     ! Reached end of model, prepare to stop
        keep_going = .false.
     end if

!        Count number of iterations; may exceed 2**31, so use a real*8 helper
     n_iter = n_iter + 1
     if (n_iter .gt. 2000000000) then
        rn_iter = rn_iter + dble(n_iter)
        n_iter = 0
        if (rn_iter .ge. 9.d9) goto 2000
     end if

     write (lun_h, 9000) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
          o_m%peri/drad, o_m%m/drad, h, epoch, comments(1:nchar)

     goto 100
!     end of the if( keep_going ) loop    
  end if

2000 continue
  write (lun_h, '(''#'')')
  write (lun_h, '(''# Total number of objects:   '', f11.0)') &
       rn_iter + dble(n_iter)
  close (lun_h)

  call exit (0)

9000 format (f9.4,1x,5(f8.4,1x),f6.2,1x,f13.5,1x,a)

9500 continue
  write (screen, *) 'File --- "', det_outfile, '" already exists. '
  goto 9502

9502 continue
  write (screen, *) 'Make sure "', det_outfile, &
       '" do not exist and restart ModelDriver.'
  call exit(-1)

9999 continue
  write (screen, *) 'Usage: ModelDriver < input'
  write (screen, *)
  write (screen, *) 'This will read the following from the keyboard:'
  write (screen, *) '<seed>'
  write (screen, *) '<nmax>'
  write (screen, *) '<Hmax>'
  write (screen, *) 'rec'
  write (screen, *) '<det_file>'
  write (screen, *) 'where:'
  write (screen, *) &
       '<seed>: integer used as seed for the random number generator'
  write (screen, *) '<nmax>: -maximum number of trials if < 0'
  write (screen, *) '        Calibrated number of object if >= 0'
  write (screen, *) '<Hmax>: maximum value of H up to which to draw'
  write (screen, *) &
       '<rec>: do we want to record the model objects inside GiMeObj'
  write (screen, *) '<det_file>: output file for model objects'

  call exit(0)

  2010 continue
  write (lun_h,*) "# Failed while generating model object."
  call exit(-20)


end program Driver
