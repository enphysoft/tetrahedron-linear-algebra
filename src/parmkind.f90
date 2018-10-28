module parmkind ! Type kind paramters 
  implicit none
  private
  public :: rkind, ikind, ckind, success, i,j,k, ip, jp, kp,ihalf, xrand 
  integer   ,parameter :: rdigits=8
  integer   ,parameter :: idigits=8
  integer   ,parameter :: rkind = 8 ! selected_real_kind(rdigits)
  integer   ,parameter :: ikind = 4 ! selected_int_kind (idigits)
  integer   ,parameter :: ckind = selected_char_kind("ascii")
  logical              :: success, message 
  integer(ikind)       :: i,j,k, ip, jp, kp 
  integer(ikind)       :: ihalf
  real   (rkind)       :: xrand
end module
