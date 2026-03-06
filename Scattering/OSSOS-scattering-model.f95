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
! File generated on 2026-03-02
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (seed, nmax, hmax, rec, o_m, epoch, h, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model parametric model of the
! outer/detached population.
!
! Version 1.0 draws according to a distributino of the form
! P(a) x P(q) x P(i_free) x P(H_r).
!
! P(a) = \left(\frac{a - a_{min}}{a_{max} - a_{min}}\right)^\alpha
! with a_{min} = 30.0d0, a_{max} = 700.0d0 and \alpha = 0.9d0
!
! P(q) is defined as follows.
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is
! $q_{min} = max(10, q_3 (a/a_3)^\beta)$,
! with $\beta = log10(q_3/q_4)/log10(a_3/a_4)$,
! $a_3 = 100$, $q_3 = 10$, $a_4 = 700$ and $q_4 = 40$
! Given this,
! - $P(q)$ is uniform between 10 and 20 with total fraction fr_1 = 0.040d0;
!   if q_{min} > 10, then fr_1 is reduced according to reduced range
! - $P(q)$ is uniform between 20 and 32 with total fraction fr_2 = 0.20d0;
!   if q_{min} > 20, then fr_2 is reduced according to reduced range
! - for the rest, $P(xq)$ is increasing as a power of xq with remaining fraction;
!   xq in range [fr_1+fr_2; 1] corresponds linearly to q in range
!   [max(32, q_{min}); q_{max}]
! with:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = 0$, $xq_1 = 1$, $q_0 = q_{min}$,
! q_{max} = 33.0d0*(a/43.0d0)^\alpha,
! \alpha = Log(33.0d0/46.0d0)/Log(43.0d0/425.0d0)
! All this has been determined empirically based on stability diagrams.
!
! P(i_free) is the usual Brown function of width 19°
!
! P(H_r) is the analytical size distribution for
! hot from Petit et al. (2023), ApJL, 947:L4. Implementation: =H_draw_hot_6=.
!
! The other angles follow a factorized uniform probability.
!
! Input files and paramters for the model are hardcoded to avoid misuse.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Institut UTINAM, UMR 6213 CNRS-UMLP, OSU THETA, Besançon, France
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
         nparam,                &! Number of parameters for the function called
                                 ! by routine incdism
         n_calib,               &! Calibrated number of objects at H_calib
         n_iter,                &! Number iterations
         n_hits                  ! Number of draws wit H <= H_calib
    real (kind=8), save :: &
         h_calib,               &! H value for calibration
         h_max,                 &! max(h_calib, hmax)
         h_params(60),          &! Parameters for the H distribution
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
         n_calib /66000/,          &! Number of object at H_calib
         h_calib /8.66d0/

! This is the first call
    if (first) then
! Reads in other parameters describing the model.
!       open (unit=lun_m, file=filena, status='old', err=1000)
!       read (lun_m, *) n_h
       n_h = 3
! hot
!       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
       h_params(1:3) = [8.3d0, 0.316d0, 13.5d0]
       h_max = max(h_calib, hmax)
       h_params(3:3:3) = h_max
       comp = 1
!       read (lun_m, *) (param(comp*10+i),i=1,2)
       param(comp*10+1:comp*10+2) = [0.0d0, 19.0d0]
       do i = 1, 2
          param(comp*10+i) = param(comp*10+i)*drad
       end do
       comp = 2
!       read (lun_m, *) (param(comp*10+i),i=1,3)
       param(comp*10+1:comp*10+3) = [30.0d0, 700.0d0, 0.90d0]
       comp = 3
!       read (lun_m, *) (param(comp*10+i),i=1,3)
       param(comp*10+1:comp*10+3) = [0.040d0, 0.200d0, 2.50d0]
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
       write (lun_ll, '(a)') '# Scatering population model.'
       write (lun_ll, '(a)') '# Version OSSOS 1.0, 2026-02-02'
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
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# Parameters for a distribution: ', &
            (param(comp*10+i), i=1,3)
       comp = 3
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# Parameters for q distribution: ', &
            (param(comp*10+i), i=1,3)
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
    commen = 'scattering_'
    nchar = 11

1100 continue
!
    ! For this component, draw 'a'
    o_m%a = scat_a_2(seed, 3, param(20+1))
!
! Index of distributions for incdism
!   3: offgau
!
! Draw 'q'
    q = scat_q_4(seed, 3, param(30+1), o_m%a)
    if (q .ge. 0.82d0*o_m%a) goto 1100
!
! Determination of "e"
    o_m%e = 1.d0 - q/o_m%a
!
! Now select inclination cell
1150 continue
    nparam = 2
    fp%f_ptr => offgau
    call incdism (seed, nparam, param(10+1), 0.0d0*drad, 110.d0*drad, inc_f, &
         3, ierr, fp)
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

  real (kind=8) function scat_q_4 (seed, np, p, a)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering q
! distribution, described below.
! The first values have been eye-balled on the scattering q distribution.
!
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is
! $q_{min} = max(10, q_3 (a/a_3)^\beta)$, with $\beta = log10(q_3/q_4)/log10(a_3/a_4)$,
! $a_3 = 100$, $q_3 = 10$, $a_4 = 700$ and $q_4 = 40$
! Given this,
! - $P(q)$ is uniform between 10 and 20 with total fraction fr_1; if q_{min} > 10,
!   then fr_1 is reduced according to reduced range
! - $P(q)$ is uniform between 20 and 32 with total fraction fr_2; if q_{min} > 20,
!   then fr_2 is reduced according to reduced range
! - for the rest, $P(xq)$ is increasing as a power of xq with remaining fraction;
!   xq in range [fr_1+fr_2; 1] corresponds linearly to q in range
!   [max(32, q_{min}); q_{max}]
! with:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = 0$, $xq_1 = 1$, $q_0 = q_{min}$.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2023
!             In this version I set the various points that define $\alpha$
!             and $\beta$. I can vary fr_1 and fr_2 only
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): fr_1
!             p(2): fr_2
!             p(2): al, the exponent of the increasing distribution of xq
!     a     : Value of semimajor axis
!
! OUTPUT
!     scat_q_4: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) a
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), a
    integer (kind=4) :: i
    real (kind=8), save :: random, a1, qa1, a2, qa2, a3, qa3, a4, qa4
    real (kind=8), save :: q0, q1, q2, qm, al, al1
    real (kind=8), save :: alpha, beta, fr1, fr2, xq, qmax, qmin, f1, f2
    logical, save :: first

    data a1 /43.0d0/, qa1 /33.0d0/, a2 /425.0d0/, qa2 /46.0d0/, &
         a3 /100.0d0/, qa3 /10.0d0/, a4 /700.0d0/, qa4 /40.0d0/
    data q0 /10.0d0/, q1 /20.0d0/, q2 /32.0d0/
    data fr1 /0.1d0/, fr2 /0.2d0/, al /2.5d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          fr1 = p(1)
       end if
       if (np .ge. 2) then
          fr2 = p(2)
       end if
       if (np .ge. 3) then
          al = p(3)
       end if
       alpha = log10(qa1/qa2)/log10(a1/a2)
       beta = log10(qa3/qa4)/log10(a3/a4)
       al1 = 1.0d0/al
       first = .false.
    end if

    qmax = qa1*(a/a1)**alpha
    qmin = max(q0, qa3*(a/a3)**beta)
    f1 = min(fr1, max(0.0d0, fr1*(q1-qmin)/(q1-q0)))
    f2 = min(fr2, max(0.0d0, fr2*(q2-qmin)/(q2-q1)))
    random=ran_3(seed)
    if (random .le. f1) then
       scat_q_4 = q0 + random*(q1-q0)/f1
    else if (random .le. f1+f2) then
       scat_q_4 = q1 + (random-f1)*(q2-q1)/f2
    else
       xq = ((random-f1-f2)/(1.0d0-f1-f2))**al1
       qm = max(q2, qmin)
       scat_q_4 = qm + xq*(qmax-qm)
    end if

    return
  end function scat_q_4

  real (kind=8) function scat_a_2 (seed, np, p)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering a
! distribution, represented by a power-law of the shifted semimajor axis.
! The first parameter values have been eye-balled on the scattering a distribution.
!
! The cumulative probability of $a$ is
! \begin{displaymath}
! P(a) = \left(\frac{a - a_{min}}{a_{max} - a_{min}}\right)^\alpha
! \end{displaymath}
! where $a_{min}$ and $a_{max}$ are the smallest and largest semimajor axis
! and $\alpha$ is the index that will set the global shape.
! 
! Once a probablity $P(a)$ is drawn in [0; 1], one can easily obtain $a$ with:
! \begin{displaymath}
! a = a_{min}} + (a_{max} - a_{min})P(a)^{1/\alpha}
! \end{displaymath}
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): minimum value of a
!             p(2): maximum value of a
!             p(3): exponential slope of distribution
!
! OUTPUT
!     scat_a_2: Random value of a (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*)
    integer (kind=4) :: i
    real (kind=8), save :: amini, amaxi, alpha, al1
    real (kind=8) :: random
    logical, save :: first

    data amini /30.0d0/, amaxi /700.0d0/, alpha /0.5d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          amini = p(1)
       end if
       if (np .ge. 2) then
          amaxi = p(2)
       end if
       if (np .ge. 3) then
          alpha = p(3)
       end if
       al1 = 1.0d0/alpha
!       print *, 'Scat_a_2'
!       print *, amini, amaxi, alpha
!       print *, al1
       first = .false.
    end if

    random=ran_3(seed)
    scat_a_2 = amini + (amaxi - amini)*random**al1

    return
  end function scat_a_2

end module gimeobjut
