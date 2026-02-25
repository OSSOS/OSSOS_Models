module ioutils

  use datadec

contains

  subroutine trim (base_name, i1, i2, finished, len)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine trims the input character string 'base_name' and returns
! the 
!
! Author: Jean-Marc Petit
! version 1: December 2019
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     base_name: Input string with file name in it (CH)
!     len   : Length of input string (I4)
!
! OUPUT
!     i1    : Index of first character of file name (I4)
!     i2    : Index of last character of file name (I4)
!     finished: Whether it reached the end of the string when searching
!               for first character
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) base_name
!f2py intent(out) i1
!f2py intent(out) i2
!f2py intent(out) finished
!f2py intent(in) len
    implicit none

    integer, intent(in) :: len
    integer, intent(out) :: i1, i2
    character(*), intent(in) :: base_name
    logical, intent(out) :: finished

    finished = .false.
    i1 = 1
100 continue
    if ((base_name(i1:i1) .eq. char(0)) &
         .or. (base_name(i1:i1) .eq. char(9)) &
         .or. (base_name(i1:i1) .eq. ' ')) then
       i1 = i1 + 1
       if (i1 .eq. len) then
          finished = .true.
          return
       end if
       goto 100
    end if
101 continue

    i2 = len
110 continue
    if ((base_name(i2:i2) .eq. char(0)) &
         .or. (base_name(i2:i2) .eq. char(9)) &
         .or. (base_name(i2:i2) .eq. ' ')) then
       if (i2 .eq. i1) goto 111
       i2 = i2 - 1
       goto 110
    end if
111 continue

    return
  end subroutine trim

  subroutine read_file_name (base_name, i1, i2, finished, len)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine trims the input character string 'base_name' and returns
! the 
!
! Author: Jean-Marc Petit
! version 1: September 2017
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     base_name: Input string with file name in it (CH)
!     len   : Length of input string (I4)
!
! OUPUT
!     i1    : Index of first character of file name (I4)
!     i2    : Index of last character of file name (I4)
!     finished: Whether it reached the end of the string when searching
!               for first character
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) base_name
!f2py intent(out) i1
!f2py intent(out) i2
!f2py intent(out) finished
!f2py intent(in) len
    implicit none

    integer, intent(in) :: len
    integer, intent(out) :: i1, i2
    character(*), intent(in) :: base_name
    logical,intent(out) :: finished

    finished = .false.
    i1 = 1
100 continue
    if ((base_name(i1:i1) .eq. char(0)) &
         .or. (base_name(i1:i1) .eq. char(9)) &
         .or. (base_name(i1:i1) .eq. ' ')) then
       i1 = i1 + 1
       if (i1 .eq. len) then
          finished = .true.
          return
       end if
       goto 100
    end if
101 continue

    i2 = i1 + 1
110 continue
    if ((base_name(i2:i2) .ne. char(0)) &
         .and. (base_name(i2:i2) .ne. char(9)) &
         .and. (base_name(i2:i2) .ne. ' ')) then
       if (i2 .eq. len) goto 111
       i2 = i2 + 1
       goto 110
    end if
    i2 = i2 - 1
111 continue

    return
  end subroutine read_file_name

  character(100) function strip_comment(str)
    character (*), intent(in) :: str
    integer c_idx, i1, i2
    logical finished
    strip_comment = str
    c_idx = index(strip_comment, '!')
    if (c_idx /= 0) then
       strip_comment = strip_comment(1:c_idx-1)
    end if
    call trim(strip_comment, i1, i2, finished, len(strip_comment))
    strip_comment = strip_comment(i1:i2)

  end function strip_comment

end module ioutils

