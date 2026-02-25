module datadec

  ! define length of array parameters
  integer, parameter :: nw_max = 20, nparmax = 20

  ! define some useful constants
  real (kind=8), parameter :: Pi = 3.141592653589793238d0, drad = Pi/180.0D0, &
       TwoHours = 2.d0/24.d0, TwoPi = 2.0d0*Pi, eps = 1.d-14, &
       Flatten = 0.003352813d0, & ! flattening of earth, 1/298.257
       Equat_Rad = 6378137.0d0 ! equatorial radius of earth, meters

  real (kind=8), parameter :: km2AU = 149597870.700d0


  real (kind=8), parameter :: gmb = 1.d0+1.d0/6023600.0d0+1.d0/408523.71d0 &
       +1.d0/328900.56d0+1.d0/3098708.0d0+1.d0/1047.3486d0+1.d0/3497.898d0 &
       +1.d0/22902.98d0+1.d0/19412.24d0+1.d0/1.35d8

  ! Internal variables
  real (kind=8) :: om_lim_low, om_lim_high
  common /om_lim_com/ om_lim_low, om_lim_high

  ! define data type to represent survey efficiency and pointings, and objects
  type t_orb_m
     real (kind=8) :: a, e, inc, node, peri, m
  end type t_orb_m

  type t_orb_p
     real (kind=8) :: a, e, inc, node, peri, tperi
  end type t_orb_p

  type t_v3d
     real (kind=8) :: x, y, z
  end type t_v3d

end module datadec
