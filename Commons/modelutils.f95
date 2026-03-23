module modelutils

  use datadec

  interface
     function func(nparam, param, inc)
       real (kind=8) :: func
       integer (kind=4), intent(in) :: nparam
       real (kind=8), intent(in) :: param(*), inc
     end function func
  end interface

  type func_holder
     procedure(func), pointer, nopass :: f_ptr => null()
  end type func_holder

contains
  subroutine incdism (seed, nparam, param, incmin, incmax, inc, &
       dist, ierr, func)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a variable according to probability
! density \verb|func| with parameters \verb|param|. Same as previous
! routine, but can remember up to 10 different distributions at a time.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     incmin: Minimum inclination (R8)
!     incmax: Maximum inclination (R8)
!     dist  : index of the selected distribution (I4)
!     func  : probability density function
!
! OUTPUT
!     inc   : Inclination (R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : wrong input data
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) incmin
!f2py intent(in) incmax
!f2py intent(in) dist
!f2py intent(out) inc
!f2py intent(out) ierr
    implicit none

    integer (kind=4), intent(in) :: nparam, dist
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(out) :: ierr
    real (kind=8), intent(in) :: param(nparmax), incmin, incmax
    real (kind=8), intent(out) :: inc
    type(func_holder), intent(in) :: func
    integer :: i, ilo, ihi, di
    real (kind=8) :: random
    integer (kind=4), parameter :: np = 16384, nd = 10
    real (kind=8), save :: proba(0:np,nd), inctab(0:np,nd)
    logical, save :: first(nd)

    data first /.true.,.true.,.true.,.true.,.true., &
         .true.,.true.,.true.,.true.,.true./

    ierr = 0
    di = min(nd, dist)
    if (first(di)) then
       inctab(0,di) = incmin
       proba(0,di) = 0.d0
       do i = 1, np
          inctab(i,di) = incmin + dfloat(i)*(incmax-incmin)/dfloat(np)
          proba(i,di) = func%f_ptr(nparam, param(1:nparam), inctab(i,di)) + proba(i-1,di)
       end do
       do i = 1, np
          proba(i,di) = proba(i,di)/proba(np,di)
       end do
       first(di) = .false.
    end if

    random = ran_3(seed)
    inc = interp(proba(0,di), inctab(0,di), random, np+1)

    return
  end subroutine incdism

  real (kind=8) function Variably_tapered(h, params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the number of objects brighter or equal to H
! following an exponentially tapered exponential with parameters in
! params.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     params: parameters for the distribution (4*R8)
!
! OUTPUT
!     Variably_tapered : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) h
!f2py intent(in) params
    implicit none

    real (kind=8), intent(in) :: params(4), h

    Variably_tapered = 10.d0**(params(3)*3.d0*(h-params(1))/5.d0) &
         *exp(-10.d0**(-params(4)*3.d0*(h-params(2))/5.d0))

    return
  end function Variably_tapered

  real (kind=8) function Variably_tapered_diff(h, params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the differential number of objects at mag H per unit
! mag following an exponentially tapered exponential with parameters in
! params.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     params: parameters for the distribution (4*R8)
!
! OUTPUT
!     Variably_tapered_diff : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) h
!f2py intent(in) params
    implicit none

    real (kind=8), intent(in) :: params(4), h

    Variably_tapered_diff = (log(10.0d0)*3.0d0/5.0d0) &
         *10.0d0**(params(3)*3.0d0*(h-params(1))/5.0d0) &
         *(params(3) + params(4)*10.0d0**(-params(4)*3.0d0*(h-params(2))/5.0d0)) &
         *exp(-10.0d0**(-params(4)*3.0d0*(h-params(2))/5.0d0))

    return
  end function Variably_tapered_diff

  real (kind=8) function onecomp (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density of Brown.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2007
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     onecomp: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t1, t3

    if (nparam .ne. 1) stop
    s1 = param(1)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -angle**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    onecomp = dsin(angle)*fe

    return
  end function onecomp

  real (kind=8) function offgau (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density as a non-zero centered gaussian.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     offgau: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t0, t1, t3

    if (nparam .ne. 2) stop
    t0 = param(1)
    s1 = param(2)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -(t0-angle)**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    offgau = dsin(angle)*fe

    return
  end function offgau

  real (kind=8) function H_dist_cold_2(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the cold belt H_r
! distribution, represented by an exponentially tapered exponential,
! with parameters I've fitted on the OSSOS cold belt data.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
! Version 2 : December 2021- updated parameter values.
! Version 3 : January 2022 - forcing slope at small sizes.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_dist_cold : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(4)
    integer (kind=4) :: i
    real (kind=8) :: params(4), random, h_min, h_max, n, c
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /4.6d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       proba(0) = 1.d-10
       n = Variably_tapered(hparam(1), params)
       c = hparam(2)
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
          if (htab(i) .lt. hparam(1)) then
             proba(i) = Variably_tapered(htab(i), params)
          else
             proba(i) = n &
                  + n*c*(10.0d0**(hparam(3)*(htab(i)-hparam(1))) - 1.0d0)
            end if
       end do
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_dist_cold_2 = interp(proba, htab, random, np+1)

    return
  end function H_dist_cold_2

  subroutine H_diff_hot_5(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the differential distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
! This version provides a continuous differential function, except for the divot.
!
! This version has modified parameters to obtain a cumulative distribution
! that looks like, and has the same normalisation as H_dist_hot_3.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(*)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: a1, a2, a3, sl1, sl2, sl3, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.7d0/
    data h1 /2.5d0/, h2 /5.8d0/
    data scale /2.2d0/

    h_max = hparam(nparam)
    h3 = hparam(1)
    a2 = scale*Variably_tapered_diff(h2, params)
    a3 = hparam(2)*scale*Variably_tapered_diff(h3, params)
    sl3 = hparam(3)
    a1 = a2*10.0d0**(sl2*(h1-h2))
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
    ! for the appropriate formula.
    do i = 0, np
       if (hs(i) .le. h_min) then
          dist(i) = 0.0d0
       else if (hs(i) .le. h1) then
          dist(i) = a1*10.0d0**(sl1*(hs(i)-h1))
       else if (hs(i) .le. h2) then
          dist(i) = a2*10.0d0**(sl2*(hs(i)-h2))
       else if (hs(i) .le. h3) then
          dist(i) = scale*Variably_tapered_diff(hs(i), params)
       else if (hs(i) .le. h_max) then
          dist(i) = a3*10.0d0**(sl3*(hs(i)-h3))
       else
          dist(i) = 0.0d0
       end if
    end do
    
    return
  end subroutine H_diff_hot_5

  real (kind=8) function H_draw_hot_5(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a differntial exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_draw_hot_5 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(*)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_diff_hot_5(nparam, hparam, np, htab, proba)
       proba(0) = 0.0d0
       do i = 1, np
          proba(i) = (proba(i)+proba(i-1))*(htab(i)-htab(i-1))/2.0d0 + proba(i-1)
       end do
       proba(0) = 1.0d-10
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_5 = interp(proba, htab, random, np+1)

  end function H_draw_hot_5

  subroutine H_diff_hot_6(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the differential distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
! This version provides a continuous differential function, except for the divot.
!
! This version has modified parameters to obtain a cumulative distribution
! that looks like, and has the same normalisation as H_dist_hot_3.
!
! Here, there is a simple divot at given mag, and then the distribution
! continues according to the expoential taper, simply rescaled by the contrast.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_5
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at divot
!             hparam(3): end of distribution
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(*)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: a1, a2, a3, sl1, sl2, sl3, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.7d0/
    data h1 /2.5d0/, h2 /5.8d0/
    data scale /2.2d0/

!    print *, nparam
!    print *, hparam(1:nparam)
    h_max = hparam(nparam)
    h3 = hparam(1)
    a2 = scale*Variably_tapered_diff(h2, params)
    a1 = a2*10.0d0**(sl2*(h1-h2))
!    print *, a1, a2, a3, h3, h_max
!    print *, scale*Variably_tapered_diff(h3, params), hparam(2)*scale*Variably_tapered_diff(h3, params)
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
    ! for the appropriate formula.
    do i = 0, np
       if (hs(i) .le. h_min) then
          dist(i) = 0.0d0
       else if (hs(i) .le. h1) then
          dist(i) = a1*10.0d0**(sl1*(hs(i)-h1))
       else if (hs(i) .le. h2) then
          dist(i) = a2*10.0d0**(sl2*(hs(i)-h2))
       else if (hs(i) .le. h3) then
          dist(i) = scale*Variably_tapered_diff(hs(i), params)
       else if (hs(i) .le. h_max) then
          dist(i) = hparam(2)*scale*Variably_tapered_diff(hs(i), params)
       else
          dist(i) = 0.0d0
       end if
    end do
    
    return
  end subroutine H_diff_hot_6

  real (kind=8) function H_draw_hot_6(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a differntial exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_5
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at divot
!             hparam(3): end of distribution
!
! OUTPUT
!     H_draw_hot_6 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(*)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_diff_hot_6(nparam, hparam, np, htab, proba)
       proba(0) = 0.0d0
       do i = 1, np
          proba(i) = (proba(i)+proba(i-1))*(htab(i)-htab(i-1))/2.0d0 + proba(i-1)
       end do
       proba(0) = 1.0d-10
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_6 = interp(proba, htab, random, np+1)

  end function H_draw_hot_6

  real (kind=8) function interp (x, y, val, n)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function linearly interpolates the function y(x) at value x=val.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     x     : Abscissa of the function, sorted in ascending order (n*R8)
!     y     : Values of the function (n*R8)
!     val   : Value of x at which to interpolate (R8)
!     n     : Size of x and y arrays (I4)
!
! OUTPUT
!     interp: Interpolated value (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) n
!f2py intent(in), depend(n) :: x
!f2py intent(in), depend(n) :: y
!f2py intent(in) val
    implicit none

    integer (kind=4), intent(in) ::  n
    real (kind=8), intent(in) :: x(*), y(*), val
    integer ::  ilo, ihi, i

    if (val .le. x(1)) then
       interp = y(1)
    else if (val .ge. x(n)) then
       interp = y(n)
    else
       ilo = 1
       ihi = n
1000   continue
       if (ihi - ilo .gt. 1) then
          i = (ihi + ilo)/2
          if (x(i) .lt. val) then
             ilo = i
          else if (x(i) .gt. val) then
             ihi = i
          else
             interp = y(i)
             return
          end if
          goto 1000
       end if
       interp = y(ilo) + (y(ihi) - y(ilo))*(val - x(ilo))/(x(ihi) - x(ilo))
    end if

    return
  end function interp

  real (kind=8) function ran_3(idum)
!f2py intent(in,out) idum
    INTEGER (KIND=4), intent(inout) :: idum
    INTEGER (KIND=4), parameter :: MBIG=1000000000, MSEED=161803398, MZ=0
    REAL (KIND=8), parameter :: FAC=1.d0/MBIG
    INTEGER :: i,ii,k,mj,mk
    INTEGER (KIND=4), save :: iff,inext,inextp,ma(55)
    data iff /0/
    if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=abs(MSEED-abs(idum))
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
       end do
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
       end do
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran_3=mj*FAC
    return
  end function ran_3

  subroutine zero2pi (var)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This function resets variable 'var' to be between 0 and 2pi
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! B. Gladman  UBC
! Version 1 : January 2007
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT/OUPUT
!     var   : Variable to reset to be between 0 and 2*Pi (R8)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) var
    implicit none

! Calling arguments
    real (kind=8), intent(inout) :: var

771 if (var .gt. TwoPi) then
       var = var - TwoPi
       goto 771
    endif
772 if (var .lt. 0.0d0) then
       var = var + TwoPi
       goto 772
    endif
    return
  end subroutine zero2pi

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

    data iset /0/
      
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

end module modelutils
