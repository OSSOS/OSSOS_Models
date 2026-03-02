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
! File generated on 2026-02-24
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (seed, nmax, hmax, rec, o_m, epoch, h, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model parametric model of the main
! classical belt.
!
! Version 1.0 has 4 components (cold, low-a; cold, high-a; hot; kernel embeded
! inside cold, low-a). For each component, it assumes a distribution of the
! form P(a, q) x P(i_free) x P(H_r) where P(H_r) and P(i_free) are analytical
! functions and P(a, q) is the marginalized debiased distribution from
! Bias-aq-i-H-v6.0-All_Surveys.
! For the hot population, removes the low inclination objects in the secular
! resonance: - a < 42.5 and i_free < 5.5 + (a - 42.5)/(40 - 42.5)*(14. - 5.)
! Debiasing was done using the free inclination and inclination and a-dependent
! forced plane derived from double averaging (Huang et al. (2022), ApJS, 259:54).
!
! The input model for (a, q) distribution is a smooothed version of the
! actual debiasing of real detections using KDE, with one set for cold component
! (low a, hight a, kernel) and another set for hot component.
!
! P(i_free) is an analytical function, one for each of the (cold, low-a + kernel),
! (cold, hight-a) and (hot). P(i_free) for hot is the usual Brown function, with
! width 16.5°
!
! P(H_r) is the analytical size distribution for cold from Kavelaars et al.
! (2021), ApJL, 920:28 for cold, low-a; kernel and cold, high-a, and the one for
! hot from Petit et al. (2023), ApJL, 947:L4. Implementation: =H_draw_hot_5=.
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
!     hmax  : Maximum value of H for model objects - cannot be less than
!             H_calib [mag] (R8)
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
!    1: offgau
!    2: cold_inc_sq
!    3: warm_inc_sq
!    4: cold_inc_sq
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
         lun_ll = 21,           &! Logical unit number for logging
         n_tab_max = 20000,     &! max number of cells in model
         n_comp = 4,            &! number of components
         nf_max = 4              ! max number of files for model
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
         i, i1, i2,             &! Dummy indices
         naq_tab(n_comp),       &! Actual number of object read (a,q)
         comp,                  &! Index of the component
         nparam,                &! Number of parameters for the function called
                                 ! by routine incdism
         n_h,                   &! Number of parameters for H distributions
         nfiles,                &! Number of input files for model
         n_calib,               &! Calibrated number of objects at H_calib
         n_iter,                &! Number iterations
         n_hits                  ! Number of draws wit H <= H_calib
    real (kind=8), save :: &
         fr_cl, fr_ch,          &! Fractions of cold, low-a, cold, high-a
         fr_h, fr_k,            &! hot and warm
         hcut,                  &! Largest Hx value in debiased model
         h_calib,               &! H value for calibration
         h_max,                 &! max(hmax, h_calib)
         h_params(60),          &! Parameters for the H distribution
         inc_f, node_f, peri_f, &! Free inclination and node and arg of peri
         i_ref, om_ref,         &! Coordinates of forced plane
         param(50),             &! Temporary storage for distribution parameters
         q,                     &! Perihelion distance
         random,                &! Random number
         r,                     &! Distance of object to Sun
         colorC(10),            &! Color parameters of model for cold component
         colorH(10),            &! Color parameters of model for hot component
         ra, dec,               &!
         delta, or, mag, alpha, &!
         rn_iter,               &! Number iterations
         a_mins(n_comp),        &! a lower limits of cells
         q_mins(n_comp),        &! q lower limits of cells
         a_min, q_min,          &! (a, q) lower limits of current cell
         a_step, q_step,        &! (a, q) sizes of cell
         a_tab(n_tab_max,n_comp),&! Semi-major axis array
         q_tab(n_tab_max,n_comp),&! Pericenter distance array
         waq_tab(0:n_tab_max,n_comp)! Bias array (a,q distribs)
    logical, save :: &
         bool,                  &! Dummy logical value
         hot,                   &! Is this the hot component ?
         first                   ! Tells if first call to routine
    character(100) :: debiasfile(nf_max)
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
         n_calib /53000/,          &! Number of object at H_calib
         h_calib /8.66d0/

! This is the first call
    if (first) then
! Reads in other parameters describing the model.
!       open (unit=lun_m, file=filena, status='old', err=1000)
!       read (lun_m, *) nfiles
       nfiles = 2
! Read in model file names
! Hot first
! Cold second
!       do i = 1, nfiles
!          read (lun_m, '(a)') debiasfile(i)    ! filename of models
!       end do
       debiasfile(1) = 'Hot_aq_density_Bias-aq-i-H-v6.3-All_Surveys_0.20_cut.dat'
       debiasfile(2) = 'Cold_aq_density_Bias-aq-i-H-v6.3-All_Surveys_0.20_cut.dat'
!       read (lun_m, *) n_h
       n_h = 4
! cold, low-a
!       read (lun_m, *) (h_params(i+1*n_h), i=1,n_h)
! cold, high a
!       read (lun_m, *) (h_params(i+2*n_h), i=1,n_h)
! hot
!       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
! kernel
!       read (lun_m, *) (h_params(i+3*n_h), i=1,n_h)
       h_params(1:16) = [10.0d0, 1.0d0, 0.4d0, 10.0d0, &
            10.0d0, 1.0d0, 0.4d0, 10.0d0, &
            10.0d0, 1.0d0, 0.4d0, 10.0d0, &
            10.0d0, 1.0d0, 0.4d0, 10.0d0]
       h_max = max(hmax, h_calib)
       h_params(1:16:4) = h_max
       h_params(4:16:4) = h_max
!       read (lun_m, *) fr_cl, fr_ch, fr_h, fr_k
       fr_cl = 0.280d0
       fr_ch = 0.290d0
       fr_h = 1.000d0
       fr_k = 0.090d0
       random = fr_cl + fr_ch + fr_h + fr_k
       fr_cl = fr_cl/random
       fr_ch = fr_ch/random
       fr_h = fr_h/random
       fr_k = fr_k/random
       comp = 1
!       read (lun_m, *) (param(comp*10+i),i=1,2)
       param(comp*10+1:comp*10+2) = [0.0d0, 16.50d0]
       do i = 1, 2
          param(comp*10+i) = param(comp*10+i)*drad
       end do
       comp = 2
!       read (lun_m, *) (param(comp*10+i),i=1,2)
       param(comp*10+1:comp*10+2) = [4.90d0, 0.70d0]
       param(comp*10+1) = param(comp*10+1)*drad
       param(4*10+1) = param(comp*10+1)
       param(4*10+2) = param(comp*10+2)
       comp = 3
!       read (lun_m, *) (param(comp*10+i),i=1,2)
       param(comp*10+1:comp*10+2) = [12.0d0, 0.50d0]
       param(comp*10+1) = param(comp*10+1)*drad
!       read (lun_m, *) (colorC(i),i=1,10)      ! read color array
!       read (lun_m, *) (colorH(i),i=1,10)      ! read color array
       colorC(1:10) = [0.90d0, 0.0d0, -0.5d0, -1.0d0, 1.5d0, &
            1.2d0,  0.8d0, -0.1d0, -0.5d0, 0.0d0]
       colorH(1:10) = [0.60d0, 0.0d0, -0.5d0, -1.0d0, 1.5d0, &
            1.2d0,  0.8d0, -0.1d0, -0.5d0, 0.0d0]
!       read (lun_m, *) log          ! logical variable, turn on drawing log
!       close (lun_m)
! Reads in the '(a,q,w)' arrays for the model.
       call set_proba_tabs(debiasfile, nfiles, hcut, a_tab, q_tab, &
            waq_tab, a_mins, q_mins, naq_tab, a_step, q_step, ierr)
       if (ierr .lt. 0) then
          return
       end if
! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a)') '# File: ModelUsed.dat'
       call date_and_time(date, time, zone, values)
       write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ', &
            date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' &
            //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a)') '# Classical population model.'
       write (lun_ll, '(a)') '# Version OSSOS 1.0, 2026-02-24'
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,i10)') '# Seed:', seed
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,f13.5)') '# Epoch:', epoch_m
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a)') '# Debiased main belt model files:'
       do i = 1, nfiles
          call read_file_name(debiasfile(i), i1, i2, bool, len(debiasfile(i)))
          write (lun_ll, '(a,a)') '#     ', debiasfile(i)(i1:i2)
       end do
       write (lun_ll, '(''#'')')
       write (lun_ll, '(''# a_step:  '',f4.2)') a_step
       write (lun_ll, '(''# q_step:  '',f4.2)') q_step
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,f5.2)') '# Hcut:  ', hcut
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for cold:', (colorC(i),i=1,10)
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for hot: ', (colorH(i),i=1,10)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(4(a,f5.3))') &
            '# Fractions of (cold, low-a): ', fr_cl, &
            '; (cold, high-a): ', fr_ch, &
            '; hot: ', fr_h, &
            '; kernel: ', fr_k
       write (lun_ll, '(''#'')')
       comp = 1
       write (lun_ll, '(a,2(1x,f5.2))') &
            '# Parameters for hot inclination distribution: ', &
            param(comp*10+1)/drad, param(comp*10+2)/drad
       comp = 2
       write (lun_ll, '(a,2(1x,f5.2))') &
            '# Parameters for cold inclination distribution: ', &
            param(comp*10+1)/drad, param(comp*10+2)
       comp = 3
       write (lun_ll, '(a,2(1x,f5.2))') &
            '# Parameters for warm inclination distribution: ', &
            param(comp*10+1)/drad, param(comp*10+2)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Cold, low-a H-dist. parameter:  ', (h_params(i+1*n_h), i=1,n_h)
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Cold, high-a H-dist. parameter: ', (h_params(i+2*n_h), i=1,n_h)
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Hot H-dist. parameter:          ', (h_params(i+0*n_h), i=1,n_h)
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Kernel H-dist. parameter:       ', (h_params(i+3*n_h), i=1,n_h)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,f5.2)') '# Phase effect parameter G:', gb0
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
! Select component the object belongs to.
    random=ran_3(seed)
    if (random .lt. fr_cl) then
! Cold low-a component
       commen = 'coldl_'
       nchar = 6
       comp = 2
    else if (random .lt. fr_cl+fr_ch) then
! Cold high-a component
       commen = 'coldh_'
       nchar = 6
       comp = 3
    else if (random .lt. fr_cl+fr_ch+fr_h) then
! Hot component
       commen = 'hot_'
       nchar = 4
       comp = 1
    else
! Kernel component
       commen = 'kernel_'
       nchar = 7
       comp = 4
    end if
!
! For this component, select the (a, q) cell
    random=ran_3(seed)
    i = find_index(waq_tab(0,comp), naq_tab(comp), random)
    o_m%a = a_tab(i,comp)
    q = q_tab(i,comp)
    nchar = nchar+3
!
! Determine boundaries of cell
    a_min = int(o_m%a/a_step)*a_step
    q_min = int(q/q_step)*q_step
1100 continue
!
! Randomize 'a'
    random=ran_3(seed)
    o_m%a = a_min + a_step*random
!
! Randomize 'q'
    random=ran_3(seed)
    q = q_min + q_step*random
    if (q .ge. o_m%a) goto 1100
!
! Determination of "e"
    o_m%e = 1.d0 - q/o_m%a
!
! Now select inclination cell
1150 continue
    param(1) = 0.0d0
! comp is used as index of the distribution; 1: cold, low a; 2: cold,
! high a; 3: hot
    if (comp .eq. 2) then
       nparam = 2
       fp%f_ptr => cold_inc_sq
       call incdism (seed, nparam, param(comp*10+1), 0.0d0, 5.d0*drad, inc_f, &
            comp, ierr, fp)
    else if (comp .eq. 3) then
       nparam = 2
       fp%f_ptr => warm_inc_sq
       call incdism (seed, nparam, param(comp*10+1), 0.0d0, 12.d0*drad, inc_f, &
            comp, ierr, fp)
    else if (comp .eq. 1) then
       nparam = 2
       fp%f_ptr => offgau
       call incdism (seed, nparam, param(comp*10+1), 0.0d0*drad, 50.d0*drad, inc_f, &
            comp, ierr, fp)
    else if (comp .eq. 4) then
       nparam = 2
       fp%f_ptr => cold_inc_sq
       call incdism (seed, nparam, param(comp*10+1), 0.0d0, 5.d0*drad, inc_f, &
            comp, ierr, fp)
    else
       print *, 'ERROR: unknown component '//commen
       stop
    end if
! Put in a cut for instability at low q and low i
    if (q .lt. 37.d0-inc_f/drad*0.2d0) goto 1150
! Checking for secular resonant region for hot component
    if (comp .eq. 1) then
       if (o_m%a .lt. 42.5d0) then
          if (inc_f .lt. &
               5.5d0*drad+(o_m%a-42.5d0)/(40.d0-42.5d0)*(14.d0-5.d0)*drad) &
               goto 1150
       end if
    end if
!
! H-mag distribution: Exponential cutoff and broken exponential law
    if ((comp .eq. 2) .or. (comp .eq. 3) .or. (comp .eq. 4)) then
       h = H_dist_cold_2(seed, n_h, h_params((comp-1)*n_h+1))
    else if (comp .eq. 1) then
       h = H_draw_hot_5(seed, n_h, h_params((comp-1)*n_h+1))
    end if
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

  subroutine set_proba_tabs(debfile, nfiles, hcut, a_tab, q_tab, &
       waq_tab, a_mins, q_mins, naq_tab, a_step, q_step, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This subroutine reads in a detection file with biases and returns the
! distribution function arrays (a, q, 1/bias) for the 4 components.
!
! The first file corresponds to the hot component, comp <- 1.
!
! The second file corresponds to the cold component. If a >= 44.5, comp <- 3;
! if a < 44.5 and not in kernel, comp <- 2;
! if (43.8 <= a <= 44.4) and (40.1 <= q <= 41.5), comp <- 4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1.0 : May 2018
! Version 1.1 : May 2018
! Version 1.2 : Decembre 2019 - reading only (a, q, bias)
! Version 2.0 : December 2021 - 4 components
! Version 2.1 : April 2022 - 4 components, 2 model files    
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     debfile: Object detection (a, q) with biases file name (CH)
!     nfiles: Number of files to read
!
! OUTPUT
!     hcut  : Largest Hx magnitude used in model
!     a_tab : Semi-major axis array [au] ((n,2)*R8)
!     q_tab : Pericenter distance array [au] ((n,2)*R8)
!     waq_tab: 1./bias array  for (a, q) ((n,2)*R8)
!     a_mins: a lower limits of cells [au] (n_comp*R8)
!     q_mins: q lower limits of cells [au] (n_comp*R8)
!     naq_tab: Number of entries in (a,q) arrays (2*I4)
!     a_step: Cell size in 'a' [au] (R8)
!     q_step: Cell size in 'q' [au] (R8)
!     ierr  : return code (I4)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) debfile
!f2py intent(in) nfiles
!f2py intent(out) hcut
!f2py intent(out) a_tab
!f2py intent(out) q_tab
!f2py intent(out) waq_tab
!f2py intent(out) a_mins
!f2py intent(out) q_mins
!f2py intent(out) naq_tab
!f2py intent(out) a_step
!f2py intent(out) q_step
!Cf2py intent(out) ierr

    implicit none
! Some values better set up as parameters
    integer (kind=4), parameter :: &
         lun_i = 22,              &! Logical unit number for model input
         n_tab_max = 20000,       &! Maximum size of input file
         n_comp = 4,              &! Number of components
         nw_max = 50               ! Maximum number of words
    real (kind=8), parameter :: &
         Pi = 3.141592653589793238d0, &! Pi
         TwoPi = 2.0d0*Pi,            &! 2*Pi
         drad = Pi/180.0d0,           &! Degree to radian convertion: Pi/180
         akl = 43.8d0, akh = 44.4d0,  &! a limits of kernel
         qkl = 40.1d0, qkh = 41.5d0,  &! q limits of kernel
         a_gap = 44.5d0                ! gap
! Calling arguments
    character(100), intent(in) :: debfile(n_comp)
    integer (kind=4), intent(in) :: &
         nfiles                    ! Number of input files
    real (kind=8), intent(out) :: hcut ! Largest Hx magnitude used in model
    real (kind=8), intent(out) :: &
         a_tab(n_tab_max,n_comp), &! Semi-major axis array
         q_tab(n_tab_max,n_comp), &! Pericenter distance array
         waq_tab(0:n_tab_max,n_comp),&! Bias array (a, q)
         a_mins(n_comp),          &! a lower limits of cells
         q_mins(n_comp),          &! q lower limits of cells
         a_step, q_step            ! (a, q) sizes of cell
    integer (kind=4), intent(out) :: &
         naq_tab(n_comp),         &! Number of entries in arrays (a,q distribs)
         ierr                      ! Error code
! Internal storage
    integer (kind=4) :: &
         i1, i2, flag, comp, nw, lw(nw_max)
    real (kind=8) :: &
         a, q, w, mean, std, tmp
    logical :: finished
    character(400) :: line
    character(20) :: word(nw_max)

! Initializes some values
    a_step = 0.0d0
    q_step = 0.0d0
    do comp = 1, n_comp
       a_mins(comp) = 0.0d0
       q_mins(comp) = 0.0d0
       naq_tab(comp) = 0
       waq_tab(0,comp) = 0.0d0
    end do

    ierr = 0
! Reads in the '(a,q,w)' arrays for the hot model.
    call read_file_name(debfile(1), i1, i2, finished, len(debfile(1)))
    open (unit=lun_i, file=debfile(1)(i1:i2), status='old', err=1000)
1150 continue
    do i1 = 1, len(line)
       line(i1:i1) = ' '
    end do
    read (lun_i, '(a)', err=1151, end=1151) line
    if (line(1:1) .eq. '#') then
       if (line(3:9) .eq. 'a_min =') read(line(10:), *) a_mins(1)
       if (line(3:9) .eq. 'q_min =') read(line(10:), *) q_mins(1)
       if (line(3:10) .eq. 'a_step =') read(line(11:), *) a_step
       if (line(3:10) .eq. 'q_step =') read(line(11:), *) q_step
       if (line(3:8) .eq. 'Hcut =') read(line(9:), *) hcut
       goto 1150
    end if
    read(line, *) a, q, w
!
! Determine component. Here hot = 1
    comp = 1
    naq_tab(comp) = naq_tab(comp) + 1
    a_tab(naq_tab(comp), comp) = a
    q_tab(naq_tab(comp), comp) = q
    waq_tab(naq_tab(comp), comp) = w
    if (naq_tab(comp) .ge. n_tab_max) then
       print *, &
            'WARNING: input file size exceeds maximum allowed ' &
            //'size', n_tab_max, &
            '. Restricting input to that size.'
       goto 1151
    end if
    goto 1150
1151 continue
    close (lun_i)

! Reads in the '(a,q,w)' arrays for the cold model.
    call read_file_name(debfile(2), i1, i2, finished, len(debfile(2)))
    open (unit=lun_i, file=debfile(2)(i1:i2), status='old', err=1000)
2150 continue
    do i1 = 1, len(line)
       line(i1:i1) = ' '
    end do
    read (lun_i, '(a)', err=2151, end=2151) line
    if (line(1:1) .eq. '#') then
       if (line(3:9) .eq. 'a_min =') read(line(10:), *) a_mins(2)
       if (line(3:9) .eq. 'q_min =') read(line(10:), *) q_mins(2)
       goto 2150
    end if
    read(line, *) a, q, w
!
! Determine component. Here hot = 1
    if (a .ge. a_gap) then
       comp = 3
    else if ((akl .le. a) .and. (a .le. akh) .and. (qkl .le. q) .and. (q .le. qkh)) then
       comp = 4
    else
       comp = 2
    end if
    naq_tab(comp) = naq_tab(comp) + 1
    a_tab(naq_tab(comp), comp) = a
    q_tab(naq_tab(comp), comp) = q
    waq_tab(naq_tab(comp), comp) = w
    if (naq_tab(comp) .ge. n_tab_max) then
       print *, &
            'WARNING: input file size exceeds maximum allowed ' &
            //'size', n_tab_max, &
            '. Restricting input to that size.'
       goto 1151
    end if
    goto 2150
2151 continue
    close (lun_i)
    a_mins(3) = a_gap
    a_mins(4) = akl
    q_mins(3) = q_mins(2)
    q_mins(4) = qkl

! Generates the probability array
    do i2 = 1, n_comp
       if (naq_tab(i2) .gt. 0) then
          mean = 0.d0
          std = 0.d0
          do i1 = 1, naq_tab(i2)
             mean = mean + waq_tab(i1,i2)
             std = std + waq_tab(i1,i2)**2
          end do
          mean = mean/dble(naq_tab(i2))
          std = sqrt(std/dble(naq_tab(i2)) - mean**2)
          do i1 = 1, naq_tab(i2)
             if (waq_tab(i1,i2) .gt. mean+3.d0*std) then
                waq_tab(i1,i2) = mean+3.d0*std
             end if
          end do
       end if
    end do
    do i2 = 1, n_comp
       if (naq_tab(i2) .gt. 0) then
          do i1 = 2, naq_tab(i2)
             waq_tab(i1,i2) = waq_tab(i1-1,i2) + waq_tab(i1,i2)
          end do
          tmp = waq_tab(naq_tab(i2),i2)
          do i1 = 0, naq_tab(i2)
             waq_tab(i1,i2) = waq_tab(i1,i2)/tmp
          end do
       end if
    end do

    return

1000 continue
! If we get here, there is something really wrong, better return with
! panic code.
    ierr = -20
    return
  end subroutine set_proba_tabs

  integer (kind=4) function find_index(tabl, n, val)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This function returns the index in 'tabl' corresponding to 'val'.
! 'tabl' contains 'n' elements. 'tabl' is assumed to be sorted in increasing
! order. Index 'i' is such that 'tabl(i-1) < val <= tabl(i)'
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1.0 : April 2018
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     tabl   : Array of sorted values (n*R8)
!     n     : Number of elements in 'tabl' (I4)
!     val   : Value we are looking for in 'tabl' (R8)
!
! OUTPUT
!     find_index: index in 'tabl' array (I4)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) tabl
!f2py intent(in) n
!f2py intent(in) val
    implicit none
    integer (kind=4), intent(in) :: n
    real (kind=8), intent(in) :: tabl(0:n), val
    integer (kind=4) :: i, ilo, ihi

    ilo = 0
    ihi = n
    if (val .lt. tabl(0)) then
       find_index = -1
       return
    end if
    if (val .ge. tabl(n)) then
       find_index = n
       return
    end if
1000 continue
    if (ihi - ilo .gt. 1) then
       i = (ihi + ilo)/2
       if (tabl(i) .lt. val) then
          ilo = i
       else if (tabl(i) .gt. val) then
          ihi = i
       else
          find_index = i
          return
       end if
       goto 1000
    end if
    if (tabl(i) .lt. val) i = i + 1
    find_index = i

    return
  end function find_index

  real (kind=8) function cold_inc_sq (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     cold_inc_sq: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, a0, p
    logical, save :: first

    data a0 /5.0d0/, p /1.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          p = param(2)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if ((angle .gt. 0.0d0) .and. (angle .lt. a0)) then
       cold_inc_sq = angle**p*(a0**p - angle**p)
    else
       cold_inc_sq = 0.0d0
    end if

    return
  end function cold_inc_sq

  real (kind=8) function warm_inc_sq (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for warm objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     warm_inc_sq: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, a0, p
    logical, save :: first

    data a0 /12.0d0/, p /1.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          p = param(2)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if ((angle .gt. 0.0d0) .and. (angle .lt. a0)) then
       warm_inc_sq = angle**p*(a0**p - angle**p)
    else
       warm_inc_sq = 0.0d0
    end if

    return
  end function warm_inc_sq

end module gimeobjut
