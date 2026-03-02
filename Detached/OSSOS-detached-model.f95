module gimeobjut

  use datadec
  use elemutils
  use rot
  use modelutils
  use ioutils

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! First define variables so they are accessible from a Python wrapper
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! Values of time and planet positions for all models
  real (kind=8) :: lambdaN, epoch_m
  common /com_time/ epoch_m, lambdaN
  data lambdaN /5.489d0/, epoch_m /2453157.5d0/

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Generic model routines for Survey Simulator, version 2.0 for OSSOS
!
! Calling sequence in Driver.f95 (survey simulator driver) is:
!
! Loop (some condition):
!     call GiMeObj(arg_list_1)
!     Check model ended:
!         set exit condition
!     call Detos1(arg_list_2)
!     Check detection and tracking:
!         store results
!
! where arg_list_1 is
! (seed, nmax, hmax, rec, o_m, epoch, h, commen, nchar, ierr)
! with:
!
! INPUT
!     seed  : Random number generator seed (I4)
!     nmax  : -Maximum number of model objects, if < 0 (I4)
!     hmax  : Maximum value of H for model objects [mag] (R8)
!     rec   : Do we want to record the objects ? (logical)
!
! OUTPUT
!     o_m   : orbital elements of object (t_orb_m)
!     epoch : Time of elements [JD] (R8)
!     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code (I4)
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The model subroutines can access files using logical unit numbers from
! 10 to 15. This range in reseved for them and won't be used by the
! drivers nor SurveySubs routines.
!
! It is good practice that when first started, the GiMeObj routine
! writes a file describing the model used, the version and the date of
! the routine.
!
! Since this routine is called once for every object created, it needs
! to get all the required parameters once when it is called the first
! time, then save these values for future use.
!
! The following routine gives a working example of a model routine. It
! is probably worth reading it through.
!
! The survey simulator expects orbital elements with respect to ecliptic
! reference frame.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File generated on 2026-02-27
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (seed, nmax, hmax, rec, o_m, epoch, h, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model parametric model of the
! outer/detached population.
!
! Version 1.0
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1.0 : February 2026
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     seed  : Random number generator seed (I4)
!     nmax  : -Maximum number of model objects, if < 0 (I4)
!     hmax  : Maximum value of H for model objects [mag] (R8)
!     rec   : Do we want to record the objects ? (logical)
!
! OUTPUT
!     o_m   : orbital elements of object (t_orb_m)
!     epoch : Time of elements [JD] (R8)
!     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code (I4)
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The user can fill the 100-character 'commen' string any way they
! wish; this comment string will be printed by the driver on the output
! line of each detection.  Examples of the comment might be resonance name
! and libration amplitude, or the name of a component in the GiMeObj model
! that the object responds to. The nchar variable (passed back to Driver)
! allows the user to print only the first nchar characters of this string.
!
! This routine uses logical unit 10 to access the file containing the model.
!
! The model uses the following indices in incdism to define the various
! distributions. Remember that only indices from 1 to 10 are allowed.
!    1:
!    2: qhot
!    3: offgau
!    4:
!    8:
!    9:
!   10:
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!f2py intent(in) nmax
!f2py intent(in) hmax
!f2py intent(in) rec
!f2py intent(in,out) seed
!f2py intent(out) o_m
!f2py intent(out) epoch
!f2py intent(out) h
!f2py intent(out) commen
!f2py intent(out) nchar
!f2py intent(out) ierr
    implicit none

! Calling arguments
    real (kind=8), intent(in) :: hmax
    integer (kind=4), intent(in) :: nmax
    logical, intent(in) :: rec
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(out) :: ierr, nchar
    type(t_orb_m), intent(out) :: o_m
    real (kind=8), intent(out) :: epoch, h
    character(100), intent(out) :: commen

! Some values better set up as parameters
    integer, parameter :: &
         lun_m = 20,            &! Logical unit number for data file reading
         lun_ll = 21             ! Logical unit number for logging
    real (kind=8), parameter :: &
         Pi = 3.141592653589793238d0, &! Pi
         TwoPi = 2.0d0*Pi,      &! 2*Pi
         drad = Pi/180.0d0       ! Degree to radian convertion: Pi/180

! Internal storage
    type(t_v3d) :: p, opos, ovel
    type(func_holder) :: fp
    character(5) :: zone         ! Time zone
    character(8) :: date         ! Date of execution
    character(10) :: time        ! Time of execution
    integer (kind=4), save :: &
         values(8),             &! Date and time of execution
         flag,                  &! Tell invar_ecl_osc which direction to go
         i,                     &! Dummy index
         comp,                  &! Index of the component
         n_h,                   &! Number of parameters for H distributions
         n_q,                   &! Number of parameters for q distribution
         nparam,                &! Number of parameters for the function called
                                 ! by routine incdism
         n_calib,               &! Calibrated number of objects at H_calib
         n_iter,                &! Number iterations
         n_hits                  ! Number of draws wit H <= H_calib
    real (kind=8), save :: &
         h_calib,               &! H value for calibration
         h_max,                 &! max(h_calib, hmax)
         h_params(60),          &! Parameters for the H distribution
         q_params(60),          &! Parameters for the q distribution
         inc_f, node_f, peri_f, &! Free inclination and node and arg of peri
         i_ref, om_ref,         &! Coordinates of forced plane
         param(50),             &! Temporary storage for distribution parameters
         q,                     &! Perihelion distance
         random,                &! Random number
         r,                     &! Distance of object to Sun
         colorH(10),            &! Color parameters of model for hot component
         ra, dec,               &!
         delta, or, mag, alpha, &!
         rn_iter,               &! Number iterations
         beta_ah,               &! Index of the a distribution
         ah_min, ah_max,        &! Lower and upper limit of a distribution
         a0sl, a1sl              ! Intermediate values for a distribution
    logical, save :: &
         bool,                  &! Dummy logical value
         first                   ! Tells if first call to routine
    real (kind=8) :: &
         epoch_m,               &! Epoch of elements [JD]
         lambdaN                 ! Longitude of Neptune at epoch
! Lightcurve and opposition surge effect parameters
    real (kind=8), save :: gb0, ph0, period0, amp0

! Place some variables in common block so they can be accessed directly
! by a Python program.
    common /com_time/ epoch_m, lambdaN

! Sets initial values
    data &
         first /.true./,        &! First call
         gb0     / 0.15d0/,     &! Opposition surge effect
         ph0     / 0.00d0/,     &! Initial phase of lightcurve
         period0 / 0.60d0/,     &! Period of lightcurve
         amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

! Calibrated number of objects
    data &
         n_calib /34500/,          &! Number of object at H_calib
         h_calib /8.66d0/

! This is the first call
    if (first) then
! Reads in other parameters describing the model.
!       open (unit=lun_m, file=filena, status='old', err=1000)
!       read (lun_m, *) n_h
       n_h = 4
! hot
!       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
       h_params(1:4) = [8.3d0, 0.316d0, 0.4d0, 10.5d0]
       h_max = max(h_calib, hmax)
       h_params(4:4:4) = h_max
       comp = 1
!       read (lun_m, *) (param(comp*10+i),i=1,2)
       param(comp*10+1:comp*10+2) = [0.0d0, 18.0d0]
       do i = 1, 2
          param(comp*10+i) = param(comp*10+i)*drad
       end do
!       read (lun_m, *) n_q, q_params(1)
       n_q = 7
!       do i = 2, n_q, 2
!          read (lun_m, *) q_params(i), q_params(i+1)
!       end do
       q_params(1:n_q) = [34.0d0, 0.080d0, 35.5d0, 0.580d0, 39.0d0, &
            0.340d0, 49.0d0]
       random = sum(q_params(2:n_q:2))
       q_params(2:n_q:2) = q_params(2:n_q:2)/random
!       read (lun_m, *) beta_ah, ah_min, ah_max
       beta_ah = -3.0d0
       ah_min = 48.0d0
       ah_max = 250.0d0
       a0sl = ah_min**(beta_ah + 1.d0)
       a1sl = ah_max**(beta_ah + 1.d0)
!       read (lun_m, *) (colorH(i),i=1,10)      ! read color array
       colorH(1:10) = [0.60d0, 0.0d0, -0.5d0, -1.0d0, 1.5d0, &
            1.2d0,  0.8d0, -0.1d0, -0.5d0, 0.0d0]
!       read (lun_m, *) log          ! logical variable, turn on drawing log
!       close (lun_m)
! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a)') '# File: ModelUsed.dat'
       call date_and_time(date, time, zone, values)
       write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ', &
            date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' &
            //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a)') '# Detached population model.'
       write (lun_ll, '(a)') '# Version OSSOS 1.0, 2026-02-27'
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,i10)') '# Seed:', seed
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,f13.5)') '# Epoch:', epoch_m
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for hot: ', (colorH(i),i=1,10)
       write (lun_ll, '(''#'')')
       comp = 1
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Parameters for inclination distribution: ', &
            (param(comp*10+i)/drad, i=1,2)
       comp = 2
       write (lun_ll, '(a,1x,i3,1x,f6.3)') &
            '# Parameters for q distribution: ', &
            n_q, q_params(1)
       do i = 2, n_q, 2
          write (lun_ll, '(a,2(1x,f6.3))') &
               '#                                ', &
               q_params(i), q_params(i+1)
       end do
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# Parameters for a distribution: ', &
            beta_ah, ah_min, ah_max
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Hot H-dist. parameter:          ', (h_params(i+0*n_h), i=1,n_h)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(2(a,f6.3,/),a,f6.3)') &
            '# lightcurve initial phase: ', ph0, &
            '# lightcurve period [day]: ', period0, &
            '# lightcurve amplitude: ', amp0
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,a)') &
            '#   a        e        i      Omega    omega      M', &
            '        H       epoch        dist    comment '
       close (lun_ll)
! Initialize counters
       n_hits = 0
       n_iter = 0
       rn_iter = 0.0d0
! Change "first" so this is not called anymore
       first = .false.
    end if
!
! Only 1 component, detached
    commen = 'detached_'
    nchar = 9

1100 continue
!
! For this component, draw 'a'
    random=ran_3(seed)
    o_m%a = (a0sl + (a1sl-a0sl)*random)**(1.0d0/(beta_ah+1.0d0))
!
! Index of distributions fors incdism
!   3: offgau
!
! Draw 'q'
    q = det_q_n(seed, n_q, q_params)
    if (q .ge. o_m%a) goto 1100
!
! Determination of "e"
    o_m%e = 1.d0 - q/o_m%a
!
! Now select inclination cell
1150 continue
    nparam = 2
    fp%f_ptr => offgau
    call incdism (seed, nparam, param(10+1), 0.0d0*drad, 70.d0*drad, inc_f, &
         3, ierr, fp)
! Put in a cut for instability at low q and low i
    if (q .lt. 37.d0-inc_f/drad*0.2d0) goto 1150
!
! H-mag distribution: Exponential cutoff and broken exponential law
    h = H_draw_hot_6(seed, n_h, h_params(0*n_h+1))
!
! Angles: uniform distribution on allowable values
    random=ran_3(seed)
    o_m%node = random*TwoPi
    random=ran_3(seed)
    o_m%peri = random*TwoPi
    random=ran_3(seed)
    o_m%m = random*TwoPi
!
! Set up epoch for orbial elements
    epoch = epoch_m

!
! The model above gives orbital elements with respect to the forced
! plane reference frame (orientation depending on 'a')
! The survey simulator expects the orbital elements with respect to the
! ecliptic, so convert them.
    o_m%inc = inc_f
    node_f = o_m%node
    peri_f = o_m%peri
    call forced_plane_damp(o_m%a, inc_f/drad, i_ref, om_ref)
    flag = 1
    call ref_ecl_osc (flag, o_m, o_m, i_ref*drad, om_ref*drad, ierr)
!
! Store object if user requested
! Normally, we should explicitly open a file and write to its end, it
! seems like the pointer to the file is not retained from one call to
! the other, so simply use the default file assigned to the logical
! unit. In this case, the output file will be something like "fort.11"
    if (rec) then
       call pos_cart(o_m, p)
       open (unit=lun_ll, file='ModelUsed.dat', access='append', &
            status='old')
       write(lun_ll,101) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
            o_m%peri/drad, o_m%m/drad, h, epoch, sqrt(p%x**2+p%y**2+p%z**2), &
            commen(1:nchar)
       close (lun_ll)
101    format(f9.4,1x,5(f8.4,1x),f6.2,1x,f13.5,1x,f9.4,1x,a9)
    end if

    n_iter = n_iter + 1
    if (n_iter .gt. 2000000000) then
       rn_iter = rn_iter + dble(n_iter)
       n_iter = 0
    end if

    if (h .le. h_calib) then
       n_hits = n_hits + 1
    end if

    if (((nmax .lt. 0) .and. (rn_iter+dble(n_iter) .ge. dble(-nmax))) .or. &
         ((nmax .ge. 0) .and. ((rn_iter .ge. 1.d12) .or. (n_hits .ge. n_calib)))) &
         then
       ierr = 100
       open (unit=lun_ll, file='ModelUsed.dat', access='append', &
            status='old')
       write(lun_ll, '(''#'')')
       write(lun_ll, '(a,f5.2,a,f13.0)') &
            '# Total number of objects up to H =         ', h_max, ': ', &
            rn_iter + dble(n_iter)
       write(lun_ll, '(a,f5.2,a,i10)') &
            '# Number of objects brighter than H_calib = ', h_calib, &
            ': ', n_hits
       close (lun_ll)
       return
    end if

! Prepare return code
    ierr = 0

    return

1000 continue
! If we get here, there is something really wrong, better return with
! panic code.
    ierr = -20
    return

  end subroutine GiMeObj

  real (kind=8) function det_q_n(seed, np, param)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws "q" according to np-1 uniform distributions over np-1
! consecutive ranges, with given fraction in each range.
!
! BEWARE: sum(param(2:np:2) MUST be 1.0d0
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     np    : number of parameters    
!     param : parameters for the distribution (4*R8)
!             param(1) = q_min
!             param(2) = fraction of range 1
!             param(3) = end of range 1
!             param(4) = fraction of range 2
!             param(5) = end of range 2
!             ...
!
! OUTPUT
!     det_q_2: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Set of F2PY directives to create a Python module
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) param
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: param(np)
    real (kind=8) :: random, s1
    integer (kind=4) :: i

    random = ran_3(seed)
    s1 = 0.0d0
    do i = 2, np, 2
       if ((random .ge. s1) .and. (random .le. s1+param(i))) then
          det_q_n = param(i-1) + (param(i+1)-param(i-1))*(random-s1)/param(i)
          exit
       end if
       s1 = s1 + param(i)
    end do

    return
  end function det_q_n

end module gimeobjut
