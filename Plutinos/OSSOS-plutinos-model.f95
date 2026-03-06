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
  data lambdaN /5.83917354547d0/, epoch_m /2456505.5d0/

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Generic model routines for Survey Simulator, version 2.0 for OSSOS
!
! Calling sequence in SurveySimulator.f (survey simulator driver) is:
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
! File generated on 2026-03-04
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
         n_hits                  ! Number of draws wit H <= H_calib
    real (kind=8), save :: &
         resamp,                &! Resonance amplitude
         phi32,                 &! Resonant angle
         h_calib,               &! H value for calibration
         h_max,                 &! max(h_calib, hmax)
         h_params(60),          &! Parameters for the H distribution
         hmin,                  &! Minimum value of the Hr distribution
         inc_f, node_f, peri_f, &! Free inclination and node and arg of peri
         i_ref, om_ref,         &! Coordinates of forced plane
         param(50),             &! Temporary storage for distribution parameters
         incmin,                &! Lower limit of inclination distribution
         incmax,                &! Upper limit of inclination distribution
         q,                     &! Perihelion distance
         ecent,                 &! mean eccentricity
         ew,                    &! standard deviation of eccentricity
         sg2deg,                &! standard deviation in inclination (degree)
         sg2,                   &! standard deviation in inclination (radian)
         libampmin,             &! min of the libration amp dist
         libampc,               &! center of the libration amp dist
         libampmax,             &! max of the libration amp dist
         libamptail,            &! tail of the libration amp dist
         libamplow,             &! fraction of population in low libamp part
         libampmid,             &! fraction of population in medium libamp part
         capomc,                &! middle point of Omega distribution
         capomfrac,             &! Fraction of population with Omega < capomc
         fkoz,                  &! Fraction of kozai resonators
         h0s10,                 &! normalizing factor at low end of H-distrib
         h1s10,                 &! normalizing factor at high end of H-distrib
         p1deg,                 &! 1st gamma inc dist param
         p2deg,                 &! 2nd gamma inc dist param
         random,                &! Random number
         r,                     &! Distance of object to Sun
         colorH(10),            &! Color parameters of model for hot component
         ra, dec,               &!
         delta, or, mag, alpha, &!
         rn_iter,               &! Number iterations
         a_cent,                &! center of a distribution
         a_min, a_max,          &! Lower and upper limit of a distribution
         alp                     ! 
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
         a_cent /39.42d0/,      &! Central a of resonance
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
!       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
       h_params(1:4) = [12.0d0, 1.0d0, 0.4d0, 12.0d0]
       h_max = max(h_calib, hmax)
       h_params(4:4:4) = h_max
!       read (lun_m, *) hmin
       hmin = 6.0d0
!       read(lun_m,*) p1deg
       p1deg = 2.0d0
!       read(lun_m,*) p2deg
       p2deg = 7.5d0
!       read(lun_m,*) ew
       ew = 0.06d0
!       read(lun_m,*) ecent
       ecent = 0.19d0
!       read(lun_m,*) libampmin
       libampmin = 10.0d0
!       read(lun_m,*) libampc
       libampc = 50.0d0
!       read(lun_m,*) libampmax
       libampmax = 120.0d0
!       read(lun_m,*) libamptail
       libamptail = 160.0d0
!       read(lun_m,*) libamplow
       libamplow = 0.14d0
!       read(lun_m,*) libampmid
       libampmid = 0.76d0
!       read(lun_m,*) capomfrac
       capomfrac = 0.48d0
       capomc = 120.0d0
!       read (lun_m, *) (color(i),i=1,10)
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
       write (lun_ll, '(a)') '# Plutino population model.'
       write (lun_ll, '(a)') '# Version OSSOS 1.0, 2026-03-04'
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,i10)') '# Seed:', seed
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,1x,f13.5)') '# Epoch:', epoch_m
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for hot: ', (colorH(i),i=1,10)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(''# p1deg:     '',f9.3)') p1deg
       write (lun_ll, '(''# p2deg:     '',f9.3)') p2deg
       write (lun_ll, '(''# ew:        '',f9.3)') ew
       write (lun_ll, '(''# ecent:     '',f9.3)') ecent
       write (lun_ll, '(''# libampmin: '',f9.3)') libampmin
       write (lun_ll, '(''# libampc:   '',f9.3)') libampc
       write (lun_ll, '(''# libampmax: '',f9.3)') libampmax
       write (lun_ll, '(''# libamptail:'',f9.3)') libamptail
       write (lun_ll, '(''# libamplow: '',f9.3)') libamplow
       write (lun_ll, '(''# libampmid: '',f9.3)') libampmid
       write (lun_ll, '(''# capomfrac: '',f9.3)') capomfrac
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Hot H-dist. parameter:          ', (h_params(i+0*n_h), i=1,n_h)
       write (lun_ll, '(''# hmin:      '',f9.3)') hmin
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

! Neptune plane : i = 1.770481500357105E+00, \Omega = 1.317941625638794E+02
       i_ref = 1.770481500357105d+00
       om_ref = 1.317941625638794d+02

       incmin = (0.0d0*drad)
       incmax = (90.0d0*drad)

! Initialize counters
       n_hits = 0
       n_iter = 0
       rn_iter = 0.0d0
! Change "first" so this is not called anymore
       first = .false.
    end if

! Select the resonant argument
! (done first as this is the same for Kozai and non-Kozai)
!
!             phi32 = 3*lambda - 2*lambdaN - longperi
!
! Assume a triangular lib amp distribution
! starts at 0, ramps up linearly to libampc
! then decreases linearly to 150 degrees
!
! KV: modified from Sam's model by eliminating
! the lower limit (was 10 degrees) and increasing
! the upper limit (was 140 degrees)
! the OSSOS data requires these changes
!
! choose libration amplitude
    random = ran_3(seed)
    if (random .le. libamplow) then
       resamp = libampmin + random*(libampc-libampmin)/libamplow
    elseif (random .le. libamplow+libampmid) then
       resamp = libampc + (random-libamplow)*(libampmax-libampc)/libampmid
    else
       resamp = libamptail - (libamptail-libampmax)*sqrt(1.d0-(random-libamplow-libampmid)/(1.0d0-libamplow-libampmid))
    end if

! choose value of phi32

    random=ran_3(seed)
    phi32 = Pi + sin(2.0d0*Pi*random)*resamp*drad

! pick mean anomaly randomly over two full orbits to cover
! the two orbits per resonant cycle
    random=ran_3(seed)
    o_m%m = random*(2d0*TwoPi)

    kozflag=0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! NON-KOZAI PLUTINO SECTION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1051 continue
    random=ran_3(seed)
    o_m%e = gasdev(ecent,ew,seed)
    if ((o_m%e .lt. 0.0d0) .or. (o_m%e .gt. 0.4d0)) goto 1051

! We want a on a triangular distribution whose limits depend on e to reproduce
! the stability of non-kozai plutinos found in Sam's paper on Kozai for FiDO.
    if (o_m%e .le. 0.2d0) then
       a_min = 39.3d0 + (39.15d0-39.3d0)*o_m%e/0.2d0
       a_max = 39.5d0 + (39.68d0-39.5d0)*o_m%e/0.2d0
    else
       a_min = 39.15d0 + (a_cent-39.15d0)*(o_m%e-0.2d0)/0.2d0
       a_max = 39.68d0 + (a_cent-39.68d0)*(o_m%e-0.2d0)/0.2d0
    end if
    random=ran_3(seed)
    alp = 0.4d0
    if (random .lt. 0.5d0) then
       o_m%a= a_min + (a_cent-a_min)*(random/0.5d0)**(1.0d0/(alp+1.0d0))
    else
       o_m%a= a_max + (a_cent-a_max)*((1.d0-random)/0.5d0)**(1.0d0/(alp+1.0d0))
    end if
! KV: commented out the resonance phase space-based a for two reasons:
! 1)  semimajor axis really doesn't matter at all for detectability
!     so there's not much point. Even the above isn't strictly 
!     necessary
! 2)  using the below procedure means you aren't *really* setting
!     the eccentricity distribution to the specified gaussian 
!     because it is going to force a re-draw much more often for
!     low eccentricity draws (there *is* less phase space at low-e
!     in the resonance, so it's not unreasonable, but I don't like
!     that it means you aren't using the eccentricity distribution
!     you say you're using. I would prefer to just come up with a
!     better eccentricity distribution model.)
! center and shape from stability plots in Tiscareno paper
!      if (a .gt. (39.45d0 + 4.0d0/3.0d0*(e - 0.01d0)) ) goto 1051
!      if (a .lt. (39.45d0 - 4.0d0/3.0d0*(e - 0.01d0)) ) goto 1051
!      a = (39.45d0 +(random-0.5d0)*4.0d0/3.0d0*(e - 0.01d0))
!      if (a .lt. (39.45d0 - 4.0d0/3.0d0*(e - 0.01d0)) ) goto 1051


!        cutoff if Uranus-approaching.  
    q = o_m%a*(1.0d0-o_m%e)
    if (q .lt. 22.d0) goto 1051

! Pull an inclination
! Modified to a gamma distribution because a gaussian wasn't a good
! fit to the data and this is a better fit that still only has two 
! parameters
    random = ran_3(seed)
9595 continue
    o_m%inc = rand_gamma(p1deg, p2deg, seed)
    o_m%inc = o_m%inc*drad
    if(o_m%inc .lt. incmin .or. o_m%inc .gt. incmax) goto 9595

! node picked randomly.  
    random=ran_3(seed)
    if (random .le. capomfrac) then
       o_m%node = 0.0d0 + random*(capomc-0.0d0)/capomfrac
    else
       o_m%node = capomc + (random-capomfrac)*(360.0d0-capomc)/(1.0d0-capomfrac)
    end if
    o_m%node = (o_m%node-om_ref)*drad

! Set argument of pericenter based on resonant angle
    o_m%peri = 0.5d0*(phi32 - 3.d0*o_m%m + 2.d0*lambdaN)  - o_m%node
    call zero2pi(o_m%peri)

    goto 1999

1999 continue

!
! H-mag distribution
1888 continue    
    h = H_draw_hot_6(seed, n_h, h_params(0*n_h+1))
    if (h .lt. hmin) goto 1888
!
! Set up epoch for orbital elements
    epoch = epoch_m
!
! Stores informations on resonance
! This will get written out to the detections/tracked files
    commen = '03:02_'
    write (commen(7:17), '(f10.6,1x)') phi32/drad
    write (commen(18:28), '(f10.5,1x)') resamp
    write (commen(29:29), '(i1)') kozflag
    nchar = 29
    do i = 7, nchar
       if (commen(i:i) .eq. ' ') commen(i:i) = '_'
    end do

! Neptune plane : i = 1.770481500357105E+00, \Omega = 1.317941625638794E+02
!
! The model above gives orbital elements with respect to Neptune
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

  real (kind=8) function gasdev(x0,sigma,rs)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!  Given a center, width and seed return a value drawn from a gaussian
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    implicit none

    real (kind=8), intent(in) :: x0, sigma
    integer (kind=4), intent(inout) :: rs
    integer (kind=4) :: iset
    real (kind=8) :: v1, v2, rsq, gset
    real (kind=8) :: fac, random
    SAVE :: iset,gset
      
    if  (iset.eq.0) then
12     continue
          random = ran_3(rs)
          v1 = 2d0*random - 1d0
          random = ran_3(rs)
          v2 = 2d0*random - 1d0
          rsq = v1*v1+v2*v2
       if (rsq.ge.1.0) goto 12 

       fac=sqrt(-2d0*log(rsq)/rsq)
       gset=v1*fac
       iset=1
       gasdev=v2*fac*sigma+x0
       return 
    else
       iset=0
       gasdev=gset*sigma+x0
       return
    endif
  end function gasdev

  real (kind=8) function rand_gamma(shape, scale, rs)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ## Implementation based on "A Simple Method for Generating Gamma Variables"
! ## by George Marsaglia and Wai Wan Tsang.
! ## ACM Transactions on Mathematical Software
! ## Vol 26, No 3, September 2000, pages 363-372.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    implicit none

    real (kind=8), intent(in) :: shape,scale
    integer (kind=4), intent(inout) :: rs
    real (kind=8) :: u, v, d, c, x, xsq, g

    IF (shape <= 0.0d0) THEN
       WRITE(*,*) "Shape PARAMETER must be positive"
    END IF
    IF (scale <= 0.0d0) THEN
       WRITE(*,*) "Scale PARAMETER must be positive"
    END IF
!
    IF (shape >= 1.0d0) THEN
       d = SHAPE - 1.0d0/3.0d0
       c = 1.0d0/(9.0d0*d)**0.5
       DO while (.true.)
          x = gasdev(0d0,1d0,rs)
          v = 1.0 + c*x
          DO while (v <= 0.0d0)
             x = gasdev(0d0,1d0,rs)
             v = 1.0d0 + c*x
          END DO

          v = v*v*v
          u =ran_3(rs)
          xsq = x*x
          IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR. &
               (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) )then
             rand_gamma=scale*d*v
             RETURN
          END IF

       END DO
    ELSE
       rand_gamma = 0d0
       RETURN
    END IF

  end function rand_gamma

end module gimeobjut
