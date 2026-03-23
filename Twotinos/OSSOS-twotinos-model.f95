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
! File generated on 2026-03-19
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (seed, nmax, hmax, rec, o_m, epoch, h, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model parametric model of the
! outer/detached population.
!
! Version 1.0 draws according to a distributino of the form
! P(a) x P(e) x P(i_free) x P(H_r) x P(Res_amp) x P(\Phi_21)
!
! P(a) is linearly increasing from 0 at 47.45 to 1 at 47.8 with 80% of the
! objects and uniform from 47.8 to 48.1 with 20% of the population.
! 
! P(i_free) is a combination of 20% of the usual Brown function of width 1.7°
! (cold component) and 80% with width 10°, a somewhat not-so excited version of
! the hot component. Refered to the invariable plane.
!
! P(e) depends on whether we consider a symmetric or asymmetric island.
! - For asymmetric island Gaussian centered on 0.275, width 0.045, limited
!   to range [0.1; 0.45]; trimmed if q = a(1-e) < 22 au.
! - For symmetric island Gaussian centered on 0.21, width 0.14, limited to range
!   [0.07; 0.35]; trimmed if q = a(1-e) < 22 au.
!
! P(Res_amp) depends on the libration island, and on e. Have a look at the code
! for the precise algorithm.
!
! P(\Phi_21) is uniform on range of width Res_amp, centered on value depending
! on the island, and possibly the eccentricity.
!
! P(H_r) is the analytical size distribution for
! hot from Petit et al. (2023), ApJL, 947:L4. Implementation: =H_draw_hot_6=.
!
! The other angles (mean anomaly and logitude of node) follow a factorized
! uniform probability. The argument of perihelion is computed from phi21,
! M, node and the mean longitude of Neptune.
!
! Paramters for the model are hardcoded to avoid misuse.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Institut UTINAM, UMR 6213 CNRS-UMLP, OSU THETA, Besançon, France
! Version 1.0 : March 2026
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
         i, j,                  &! Dummy indices
         comp,                  &! Index of the component
         kozflag,               &! Flag for kozai resonance
         n_h,                   &! Number of parameters for H distributions
         nparam,                &! Number of parameters for the function called
                                 ! by routine incdism
         n_calib,               &! Calibrated number of objects at H_calib
         n_iter,                &! Number iterations
         n_hits,                &! Number of draws wit H <= H_calib
         lts                     !sym = 0, leading = 1, trailing = 2
    real (kind=8), save :: &
         resamp,                &! Resonance amplitude
         fsym,                  &! symmetric fraction
         phi21,                 &! Resonant angle
         h_calib,               &! H value for calibration
         h_max,                 &! max(h_calib, hmax)
         h_params(60),          &! Parameters for the H distribution
         hmin,                  &! Minimum value of the Hr distribution
         inc_f, node_f, peri_f, &! Free inclination and node and arg of peri
         i_ref, om_ref,         &! Coordinates of forced plane
         params(50),            &! Temporary storage for distribution parameters
         incmin,                &! Lower limit of inclination distribution
         incmax,                &! Upper limit of inclination distribution
         q,                     &! Perihelion distance
         ecas,                  &! mean eccentricity of asymmetric
         esas,                  &! sigma eccentricity of asymmetric
         emin,                  &! minimum eccentricity of asymmetric
         emax,                  &! maximum eccentricity of asymmetric
         emins,                 &! mean eccentricity
         emaxs,                 &! standard deviation of eccentricity
         emax52,                &! Maximum eccentricity of plutinos
         fleading,              &! fraction of asymm lib that are in leading island
         f1,                    &! fraction of low-inclination component
         libampc,               &! center of the libration amp dist
         fkoz,                  &! Fraction of kozai resonators
         random,                &! Random number
         r,                     &! Distance of object to Sun
         colorH(10),            &! Color parameters of model for hot component
         ra, dec,               &!
         delta, or, mag, alpha, &!
         rn_iter,               &! Number iterations
         libc, ymax, ampmax,    &!
         mslope, binter, temp,  &!
         temp2, libcl, libcf,   &!
         amplowest, ldl
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
         amp0    / 0.00d0/,     &! Amplitude of lightcurve (peak-to-peak)
         emax52  /0.32d0/        ! Maximum eccentricity of twotinos

! Calibrated number of objects
    data &
         n_calib /3200/,          &! Number of object at H_calib
         h_calib /8.66d0/

! This is the first call
    if (first) then
! Reads in other parameters describing the model.
!       open (unit=lun_m, file=filena, status='old', err=1000)
!       read (lun_m, *) n_h
       n_h = 4
!       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
       h_params(1:4) = [10.7d0, 1.0d0, 0.4d0, 10.7d0]
       h_max = max(h_calib, hmax)
       h_params(1:4:4) = h_max
       h_params(4:4:4) = h_max
!       read (lun_m, *) hmin
       hmin = 6.0d0
!       read(lun_m,*) sg1d, sg2d, f1
       params(1*10+1) = 1.7d0
       params(2*10+1) = 10.0d0
       f1 = 0.2d0
       params(1*10+1) = params(1*10+1)*drad
       params(2*10+1) = params(2*10+1)*drad
!       read(lun_m,*) emin
       emin = 0.10d0
!       read(lun_m,*) emax
       emax = 0.45d0
!       read(lun_m,*) ecas
       ecas = 0.275d0
!       read(lun_m,*) esas
       esas = 0.045d0
!       read(lun_m,*) emins
       emins = 0.070d0
!       read(lun_m,*) emaxs
       emaxs = 0.350d0
!       read(lun_m,*) fleading
       fleading = 0.5d0
!       read(lun_m,*) fsym
       fsym = 0.5d0
!       read (lun_m, *) (colorH(i),i=1,10)
       colorH(1:10) = [0.60d0, 0.0d0, -0.5d0, -1.0d0, 1.5d0, &
            1.2d0,  0.8d0, -0.1d0, -0.5d0, 0.0d0]
!       read (lun_m, *) nlog
!       close(lun_m)

! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a)') '# File: ModelUsed.dat'
       call date_and_time(date, time, zone, values)
       write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ', &
            date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' &
            //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a)') '# Twotino population model.'
       write (lun_ll, '(a)') '# Version OSSOS 1.0, 2026-03-20'
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,i10)') '# Seed:', seed
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,f13.5)') '# Epoch:', epoch_m
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for hot: ', (colorH(i),i=1,10)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(''# Cold component: '',f9.3)') &
            params(1*10+1)/drad
       write (lun_ll, '(''# Hot component: '',f9.3)') &
            params(2*10+1)/drad
       write (lun_ll, '(''# f1:        '',f9.3)') f1
       write (lun_ll, '(''# emin:      '',f9.3)') emin
       write (lun_ll, '(''# emax:      '',f9.3)') emax
       write (lun_ll, '(''# ecas:      '',f9.3)') ecas
       write (lun_ll, '(''# esas:      '',f9.3)') esas
       write (lun_ll, '(''# emins:     '',f9.3)') emins
       write (lun_ll, '(''# emaxs:     '',f9.3)') emaxs
       write (lun_ll, '(''# fleading:  '',f9.3)') fleading
       write (lun_ll, '(''# fsym:      '',f9.3)') fsym
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Hot H-dist. parameter:          ', (h_params(i+0*n_h), i=1,n_h)
       write (lun_ll, '(''# hmin:      '',f9.3)') hmin
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

! Use invariable plane as reference plane
       i_ref = 5713.86d0/3600.d0
       om_ref = 387390.8d0/3600.d0

       incmin = (0.0d0*drad)
       incmax = (90.0d0*drad)

! Initialize counters
       n_hits = 0
       n_iter = 0
       rn_iter = 0.0d0
! Change "first" so this is not called anymore
       first = .false.
    end if

!  Algorithm from Charles (Ying-Tung Chen)
!
! pick 'a' and re-draw if a/e outside of bounds
    random=ran_3(seed)
    if (random .lt. 0.8d0) then
       o_m%a = 47.45d0 + 0.35d0*(random/0.8d0)**(1.0d0/2.0d0)
    else
       o_m%a = 47.8d0 + (1.d0-random)/0.2d0*0.3d0
    end if
!
! Pull an inclination
! J-M's code:
! incdism indices:
! 1: onecomp, sg1 (same as offgau, 0, sg1)
! 2: onecomp, sg2 (same as offgau, 0, sg2)
    nparam = 1
    incmin = (0.0d0*drad)
    incmax = (90.0d0*drad)
    fp%f_ptr => onecomp
    random = ran_3(seed)
    if (random .le. f1) then
       call incdism (seed, nparam, params(1*10+1), incmin, incmax, o_m%inc, &
            1, ierr, fp)
    else
       call incdism (seed, nparam, params(2*10+1), incmin, incmax, o_m%inc, &
            2, ierr, fp)
    end if

! node and M picked randomly.  
    random=ran_3(seed)
    o_m%node = random*TwoPi
    random=ran_3(seed)
    o_m%m = random*TwoPi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! have to pick whether it's symmetric or asymmetric, and which island
! it orbits
    random=ran_3(seed)
    if (random .gt. fsym) then
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1051   continue

! new asymmetric rules:
! pick an e.
! maximum libration amplitude for that given e is roughly:
! 81.6 - 1.5631/sinh(e/1.636) based solely off a rough fit to the figure
! in Nesvorny & Roig 2001
! amplitude is then uniformly selected between 0 and that value
! based off the same figure, the libration center is given by:
! 72.75 - 33.1*x - 15.32*x**2 + 3.1616/x
!
! pick e
!          o_m%e = ran_3(seed)*(emax-emin)+emin
! give a slope distribution
133    continue
         o_m%e =  gasdev(ecas, esas, seed)
         if ((o_m%e .lt. emin) .or. (o_m%e .gt. emax)) goto 133
!        cutoff if Uranus-approaching. , reselect e 
          q = o_m%a*(1.0d0-o_m%e)
       if (q .lt. 22.000) goto 1051         
!
! pick a libration center
       if (o_m%e .lt. 0.05) then
          libc = 140d0 + (1d0-ran_3(seed))*10d0
          ampmax = 15d0
       else
          libcl = 133.157d0 -315.272d0*o_m%e+377.082d0*o_m%e**2 -0.491667/o_m%e
!        get the libc
          random=ran_3(seed)
          libc = libcl+ (1-sin(0.5d0*Pi*random))*20
          ampmax = -403.632d0 +9.09917d0*libc-0.0442498d0*libc*libc &
               -0.0883975d0/libc
          amplowest = 79.031d0*exp(-(libc-121.3435d0)**2/(2*15.51349d0**2))
          lts = 1
          random = ran_3(seed)
          if (random .gt. fleading) then
             libc = 360.0d0-libc
             lts = 2 !trailling
          endif

       end if

! Draws a number with probability linearly decreasing from 1 at 0 to 1/4 at 1.
! Can be done as follow by trial and error
!222    random =  ran_3(seed)
!       ldl = 1.0d0 - 0.75d0*random
!       if (ran_3(seed) .gt. ldl) goto 222
! or more efficiently as
       ldl = ran_3(seed)
       random = (8.0d0 - sqrt(64.0d0-60.0d0*ldl))/6.0d0
       resamp = ampmax - (random)*(ampmax - amplowest)
!
!cccccccccccccccccccccc Now choose phi
       phi21 = (libc + 2*resamp*(ran_3(seed) - 0.5))*drad
!
!             Set argument of pericenter based on resonant angle
       o_m%peri = phi21 - 2.0d0*o_m%m + lambdaN  - o_m%node
       libcf = libc
       call zero2pi(o_m%peri)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! symmetric
    else
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! symmetric, libration centers are all 180
!            pick libration amplitudes uniformly from 125-165
!            pick e uniformly from 0.05 -> 0.3 (Chiang&Jordan02)
       lts = 0
!
! pick a libration center
       libampc=180.0d0

1052   continue
! pick e
          o_m%e = gasdev((emaxs+emins)/2.0d0, (emaxs-emins)/2.0d0, seed)
          if ((o_m%e .lt. emins) .or. (o_m%e .gt. emaxs)) goto 1052
!        cutoff if Uranus-approaching.  
          q = o_m%a*(1.0d0-o_m%e)
       if (q .lt. 22.000) goto 1052       
!
! pick a libration amplitude
       resamp = ran_3(seed)*40d0 + 125d0
!
!cccccccccccccccccccccc Now choose phi
       random=ran_3(seed)
       phi21 = libampc*drad + (2*random-1)*resamp*drad
!
!             Set argument of pericenter based on resonant angle
       o_m%peri = phi21 - 2.0d0*o_m%m + lambdaN  - o_m%node
       libcf = libampc
       call zero2pi(o_m%peri)
    end if
!
! H-mag distribution
1888 continue
    h = H_draw_hot_6(seed, n_h, h_params(0*n_h+1))
    if (h .lt. hmin) goto 1888
!
! Set up epoch for orbial elements
    epoch = epoch_m
!
! Stores informations on resonance
! This will get written out to the detections/tracked files
    commen = '02:01 '
    kozflag = 0
    write (commen(7:17), '(f10.6,1x)') phi21
    write (commen(18:28), '(f10.5,1x)') resamp
    write (commen(29:29), '(i1)') lts
    write (commen(31:41), '(f10.6,1x)') libcf
    nchar = 41
!
! The model above gives orbital elements with respect to invariable
! plane reference frame.
! The survey simulator expects the orbital elements with respect to the
! ecliptic, so convert them.
    inc_f = o_m%inc
    node_f = o_m%node
    peri_f = o_m%peri
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

end module gimeobjut
