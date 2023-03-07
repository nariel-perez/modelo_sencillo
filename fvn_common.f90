module fvn_common
! This module contains routines that are used by more than one fvn submodule

implicit none
! Kind Definition Module integrated into fvn_common
integer, parameter :: ip_kind = kind(1)
integer, parameter :: sp_kind = kind(1.0E0)
integer, parameter :: dp_kind = kind(1.0D0)

! We define pi and i for the module
real(kind=dp_kind),parameter :: fvn_pi = 3.141592653589793_dp_kind
real(kind=dp_kind),parameter :: fvn_el = 0.5772156649015328_dp_kind
complex(kind=dp_kind),parameter :: fvn_i = (0._dp_kind,1._dp_kind)

! an integer variable that can be used to store the return status of different fvn subroutines
integer :: fvn_status

interface
    function d1mach(i)
        integer :: i
        double precision :: d1mach
    end function
    function r1mach(i)
        integer :: i
        real :: r1mach
    end function
end interface

end module fvn_common