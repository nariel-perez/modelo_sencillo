module fvn_interpol
use fvn_common
implicit none


! Utility procedure find interval
interface fvn_find_interval
    module procedure fvn_s_find_interval,fvn_d_find_interval
end interface fvn_find_interval

! Quadratic 1D interpolation
interface fvn_quad_interpol
    module procedure fvn_s_quad_interpol,fvn_d_quad_interpol
end interface fvn_quad_interpol

! Quadratic 2D interpolation
interface fvn_quad_2d_interpol
    module procedure fvn_s_quad_2d_interpol,fvn_d_quad_2d_interpol
end interface fvn_quad_2d_interpol

! Quadratic 3D interpolation
interface fvn_quad_3d_interpol
    module procedure fvn_s_quad_3d_interpol,fvn_d_quad_3d_interpol
end interface fvn_quad_3d_interpol

! Akima interpolation
interface fvn_akima
    module procedure fvn_s_akima,fvn_d_akima
end interface fvn_akima

! Akima evaluation
interface fvn_spline_eval
    module procedure fvn_s_spline_eval,fvn_d_spline_eval
end interface fvn_spline_eval

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Quadratic interpolation of tabulated function of 1,2 or 3 variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fvn_s_find_interval(x,i,xdata,n)
      implicit none
      ! This routine find the indice i where xdata(i) <= x < xdata(i+1)
      ! xdata(n) must contains a set of increasingly ordered values 
      ! if x < xdata(1) i=0 is returned
      ! if x > xdata(n) i=n is returned
      ! special case is where x=xdata(n) then n-1 is returned so 
      ! we will not exclude the upper limit
      ! a simple dichotomy method is used

      real(kind=sp_kind), intent(in) :: x
      integer(kind=sp_kind), intent(in) :: n
      real(kind=sp_kind), intent(in), dimension(n) :: xdata
      integer(kind=sp_kind), intent(out) :: i

      integer(kind=sp_kind) :: imin,imax,imoyen

      ! special case is where x=xdata(n) then n-1 is returned so 
      ! we will not exclude the upper limit
      if (x == xdata(n)) then
            i=n-1
            return
      end if

      ! if x < xdata(1) i=0 is returned
      if (x < xdata(1)) then
            i=0
            return
      end if

      ! if x > xdata(n) i=n is returned
      if (x > xdata(n)) then
            i=n
            return
      end if

      ! here xdata(1) <= x <= xdata(n)
      imin=0
      imax=n+1

      do while((imax-imin) > 1)
            imoyen=(imax+imin)/2
            if (x >= xdata(imoyen)) then
                  imin=imoyen
            else
                  imax=imoyen
            end if
      end do

      i=imin

end subroutine


subroutine fvn_d_find_interval(x,i,xdata,n)
      implicit none
      ! This routine find the indice i where xdata(i) <= x < xdata(i+1)
      ! xdata(n) must contains a set of increasingly ordered values 
      ! if x < xdata(1) i=0 is returned
      ! if x > xdata(n) i=n is returned
      ! special case is where x=xdata(n) then n-1 is returned so 
      ! we will not exclude the upper limit
      ! a simple dichotomy method is used

      real(kind=dp_kind), intent(in) :: x
      integer(kind=sp_kind), intent(in) :: n
      real(kind=dp_kind), intent(in), dimension(n) :: xdata
      integer(kind=sp_kind), intent(out) :: i

      integer(kind=sp_kind) :: imin,imax,imoyen

      ! special case is where x=xdata(n) then n-1 is returned so 
      ! we will not exclude the upper limit
      if (x == xdata(n)) then
            i=n-1
            return
      end if

      ! if x < xdata(1) i=0 is returned
      if (x < xdata(1)) then
            i=0
            return
      end if

      ! if x > xdata(n) i=n is returned
      if (x > xdata(n)) then
            i=n
            return
      end if

      ! here xdata(1) <= x <= xdata(n)
      imin=0
      imax=n+1

      do while((imax-imin) > 1)
            imoyen=(imax+imin)/2
            if (x >= xdata(imoyen)) then
                  imin=imoyen
            else
                  imax=imoyen
            end if
      end do

      i=imin

end subroutine


function fvn_s_quad_interpol(x,n,xdata,ydata)
      implicit none
      ! This function evaluate the value of a function defined by a set of points
      ! and values, using a quadratic interpolation
      ! xdata must be increasingly ordered
      ! x must be within xdata(1) and xdata(n) to actually do interpolation
      ! otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: n
      real(kind=sp_kind), intent(in), dimension(n) :: xdata,ydata
      real(kind=sp_kind), intent(in) :: x
      real(kind=sp_kind) ::  fvn_s_quad_interpol

      integer(kind=sp_kind) :: iinf,base,i,j
      real(kind=sp_kind) :: p

      call fvn_s_find_interval(x,iinf,xdata,n)

      ! Settings for extrapolation
      if (iinf==0) then
            ! TODO -> Lower bound extrapolation warning
            iinf=1
      end if

      if (iinf==n) then
            ! TODO -> Higher bound extrapolation warning
            iinf=n-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (iinf==1) then
            base=0
      else
            base=iinf-2
      end if

      ! The three points we will use are :
      ! xdata/ydata(base+1),xdata/ydata(base+2),xdata/ydata(base+3)

      ! Straight forward Lagrange polynomial
      fvn_s_quad_interpol=0.
      do i=1,3
            !  polynome i
            p=ydata(base+i)
            do j=1,3
                  if (j /= i) then
                        p=p*(x-xdata(base+j))/(xdata(base+i)-xdata(base+j))
                  end if
            end do
            fvn_s_quad_interpol=fvn_s_quad_interpol+p
      end do

end function


function fvn_d_quad_interpol(x,n,xdata,ydata)
      implicit none
      ! This function evaluate the value of a function defined by a set of points
      ! and values, using a quadratic interpolation
      ! xdata must be increasingly ordered
      ! x must be within xdata(1) and xdata(n) to actually do interpolation
      ! otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: n
      real(kind=dp_kind), intent(in), dimension(n) :: xdata,ydata
      real(kind=dp_kind), intent(in) :: x
      real(kind=dp_kind) ::  fvn_d_quad_interpol

      integer(kind=sp_kind) :: iinf,base,i,j
      real(kind=dp_kind) :: p

      call fvn_d_find_interval(x,iinf,xdata,n)

      ! Settings for extrapolation
      if (iinf==0) then
            ! TODO -> Lower bound extrapolation warning
            iinf=1
      end if

      if (iinf==n) then
            ! TODO Higher bound extrapolation warning
            iinf=n-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (iinf==1) then
            base=0
      else
            base=iinf-2
      end if

      ! The three points we will use are :
      ! xdata/ydata(base+1),xdata/ydata(base+2),xdata/ydata(base+3)

      ! Straight forward Lagrange polynomial
      fvn_d_quad_interpol=0.
      do i=1,3
            !  polynome i
            p=ydata(base+i)
            do j=1,3
                  if (j /= i) then
                        p=p*(x-xdata(base+j))/(xdata(base+i)-xdata(base+j))
                  end if
            end do
            fvn_d_quad_interpol=fvn_d_quad_interpol+p
      end do

end function


function fvn_s_quad_2d_interpol(x,y,nx,xdata,ny,ydata,zdata)
      implicit none
      ! This function evaluate the value of a two variable function defined by a
      ! set of points and values, using a quadratic interpolation
      ! xdata and ydata must be increasingly ordered
      ! the couple (x,y) must be as x within xdata(1) and xdata(nx) and
      ! y within ydata(1) and ydata(ny) to actually do interpolation
      ! otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: nx,ny
      real(kind=sp_kind), intent(in) :: x,y
      real(kind=sp_kind), intent(in), dimension(nx) :: xdata
      real(kind=sp_kind), intent(in), dimension(ny) :: ydata
      real(kind=sp_kind), intent(in), dimension(nx,ny) :: zdata
      real(kind=sp_kind) :: fvn_s_quad_2d_interpol

      integer(kind=sp_kind) :: ixinf,iyinf,basex,basey,i
      real(kind=sp_kind),dimension(3) :: ztmp
      !real(kind=4), external :: fvn_s_quad_interpol

      call fvn_s_find_interval(x,ixinf,xdata,nx)
      call fvn_s_find_interval(y,iyinf,ydata,ny)

      ! Settings for extrapolation
      if (ixinf==0) then
            ! TODO -> Lower x  bound extrapolation warning
            ixinf=1
      end if

      if (ixinf==nx) then
            ! TODO -> Higher x bound extrapolation warning
            ixinf=nx-1
      end if

      if (iyinf==0) then
            ! TODO -> Lower y  bound extrapolation warning
            iyinf=1
      end if

      if (iyinf==ny) then
            ! TODO -> Higher y bound extrapolation warning
            iyinf=ny-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (ixinf==1) then
            basex=0
      else
            basex=ixinf-2
      end if

      if (iyinf==1) then
            basey=0
      else
            basey=iyinf-2
      end if

      ! First we make 3 interpolations for x at y(base+1),y(base+2),y(base+3)
      ! stored in ztmp(1:3)
      do i=1,3
            ztmp(i)=fvn_s_quad_interpol(x,nx,xdata,zdata(:,basey+i))
      end do

      ! Then we make an interpolation for y using previous interpolations
      fvn_s_quad_2d_interpol=fvn_s_quad_interpol(y,3,ydata(basey+1:basey+3),ztmp)
end function


function fvn_d_quad_2d_interpol(x,y,nx,xdata,ny,ydata,zdata)
      implicit none
      ! This function evaluate the value of a two variable function defined by a
      ! set of points and values, using a quadratic interpolation
      ! xdata and ydata must be increasingly ordered
      ! the couple (x,y) must be as x within xdata(1) and xdata(nx) and
      ! y within ydata(1) and ydata(ny) to actually do interpolation
      ! otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: nx,ny
      real(kind=dp_kind), intent(in) :: x,y
      real(kind=dp_kind), intent(in), dimension(nx) :: xdata
      real(kind=dp_kind), intent(in), dimension(ny) :: ydata
      real(kind=dp_kind), intent(in), dimension(nx,ny) :: zdata
      real(kind=dp_kind) :: fvn_d_quad_2d_interpol

      integer(kind=sp_kind) :: ixinf,iyinf,basex,basey,i
      real(kind=dp_kind),dimension(3) :: ztmp
      !real(kind=8), external :: fvn_d_quad_interpol

      call fvn_d_find_interval(x,ixinf,xdata,nx)
      call fvn_d_find_interval(y,iyinf,ydata,ny)

      ! Settings for extrapolation
      if (ixinf==0) then
            ! TODO -> Lower x  bound extrapolation warning
            ixinf=1
      end if

      if (ixinf==nx) then
            ! TODO -> Higher x bound extrapolation warning
            ixinf=nx-1
      end if

      if (iyinf==0) then
            ! TODO -> Lower y  bound extrapolation warning
            iyinf=1
      end if

      if (iyinf==ny) then
            ! TODO -> Higher y bound extrapolation warning
            iyinf=ny-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (ixinf==1) then
            basex=0
      else
            basex=ixinf-2
      end if

      if (iyinf==1) then
            basey=0
      else
            basey=iyinf-2
      end if

      ! First we make 3 interpolations for x at y(base+1),y(base+2),y(base+3)
      ! stored in ztmp(1:3)
      do i=1,3
            ztmp(i)=fvn_d_quad_interpol(x,nx,xdata,zdata(:,basey+i))
      end do

      ! Then we make an interpolation for y using previous interpolations
      fvn_d_quad_2d_interpol=fvn_d_quad_interpol(y,3,ydata(basey+1:basey+3),ztmp)
end function


function fvn_s_quad_3d_interpol(x,y,z,nx,xdata,ny,ydata,nz,zdata,tdata)
      implicit none
      ! This function evaluate the value of a 3 variables function defined by a
      ! set of points and values, using a quadratic interpolation
      ! xdata, ydata and zdata must be increasingly ordered
      ! The triplet (x,y,z) must be within xdata,ydata and zdata to actually 
      ! perform an interpolation, otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: nx,ny,nz
      real(kind=sp_kind), intent(in) :: x,y,z
      real(kind=sp_kind), intent(in), dimension(nx) :: xdata
      real(kind=sp_kind), intent(in), dimension(ny) :: ydata
      real(kind=sp_kind), intent(in), dimension(nz) :: zdata
      real(kind=sp_kind), intent(in), dimension(nx,ny,nz) :: tdata
      real(kind=sp_kind) :: fvn_s_quad_3d_interpol

      integer(kind=sp_kind) :: ixinf,iyinf,izinf,basex,basey,basez,i,j
      !real(kind=4), external :: fvn_s_quad_interpol,fvn_s_quad_2d_interpol
      real(kind=sp_kind),dimension(3,3) :: ttmp

      call fvn_s_find_interval(x,ixinf,xdata,nx)
      call fvn_s_find_interval(y,iyinf,ydata,ny)
      call fvn_s_find_interval(z,izinf,zdata,nz)

      ! Settings for extrapolation
      if (ixinf==0) then
            ! TODO -> Lower x  bound extrapolation warning
            ixinf=1
      end if

      if (ixinf==nx) then
            ! TODO -> Higher x bound extrapolation warning
            ixinf=nx-1
      end if

      if (iyinf==0) then
            ! TODO -> Lower y  bound extrapolation warning
            iyinf=1
      end if

      if (iyinf==ny) then
            ! TODO -> Higher y bound extrapolation warning
            iyinf=ny-1
      end if

      if (izinf==0) then
            ! TODO -> Lower z bound extrapolation warning
            izinf=1
      end if

      if (izinf==nz) then
            ! TODO -> Higher z bound extrapolation warning
            izinf=nz-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (ixinf==1) then
            basex=0
      else
            basex=ixinf-2
      end if

      if (iyinf==1) then
            basey=0
      else
            basey=iyinf-2
      end if

      if (izinf==1) then
            basez=0
      else
            basez=izinf-2
      end if

      ! We first make 9 one dimensional interpolation on variable x.
      ! results are stored in ttmp
      do i=1,3
            do j=1,3
                  ttmp(i,j)=fvn_s_quad_interpol(x,nx,xdata,tdata(:,basey+i,basez+j))
            end do
      end do

      ! We then make a 2 dimensionnal interpolation on variables y and z
      fvn_s_quad_3d_interpol=fvn_s_quad_2d_interpol(y,z, &
            3,ydata(basey+1:basey+3),3,zdata(basez+1:basez+3),ttmp)
end function


function fvn_d_quad_3d_interpol(x,y,z,nx,xdata,ny,ydata,nz,zdata,tdata)
      implicit none
      ! This function evaluate the value of a 3 variables function defined by a
      ! set of points and values, using a quadratic interpolation
      ! xdata, ydata and zdata must be increasingly ordered
      ! The triplet (x,y,z) must be within xdata,ydata and zdata to actually 
      ! perform an interpolation, otherwise extrapolation is done
      integer(kind=sp_kind), intent(in) :: nx,ny,nz
      real(kind=dp_kind), intent(in) :: x,y,z
      real(kind=dp_kind), intent(in), dimension(nx) :: xdata
      real(kind=dp_kind), intent(in), dimension(ny) :: ydata
      real(kind=dp_kind), intent(in), dimension(nz) :: zdata
      real(kind=dp_kind), intent(in), dimension(nx,ny,nz) :: tdata
      real(kind=dp_kind) :: fvn_d_quad_3d_interpol

      integer(kind=sp_kind) :: ixinf,iyinf,izinf,basex,basey,basez,i,j
      !real(kind=8), external :: fvn_d_quad_interpol,fvn_d_quad_2d_interpol
      real(kind=dp_kind),dimension(3,3) :: ttmp

      call fvn_d_find_interval(x,ixinf,xdata,nx)
      call fvn_d_find_interval(y,iyinf,ydata,ny)
      call fvn_d_find_interval(z,izinf,zdata,nz)

      ! Settings for extrapolation
      if (ixinf==0) then
            ! TODO -> Lower x  bound extrapolation warning
            ixinf=1
      end if

      if (ixinf==nx) then
            ! TODO -> Higher x bound extrapolation warning
            ixinf=nx-1
      end if

      if (iyinf==0) then
            ! TODO -> Lower y  bound extrapolation warning
            iyinf=1
      end if

      if (iyinf==ny) then
            ! TODO -> Higher y bound extrapolation warning
            iyinf=ny-1
      end if

      if (izinf==0) then
            ! TODO -> Lower z bound extrapolation warning
            izinf=1
      end if

      if (izinf==nz) then
            ! TODO -> Higher z bound extrapolation warning
            izinf=nz-1
      end if

      ! The three points we will use are iinf-1,iinf and iinf+1 with the
      ! exception of the first interval, where iinf=1 we will use 1,2 and 3
      if (ixinf==1) then
            basex=0
      else
            basex=ixinf-2
      end if

      if (iyinf==1) then
            basey=0
      else
            basey=iyinf-2
      end if

      if (izinf==1) then
            basez=0
      else
            basez=izinf-2
      end if

      ! We first make 9 one dimensional interpolation on variable x.
      ! results are stored in ttmp
      do i=1,3
            do j=1,3
                  ttmp(i,j)=fvn_d_quad_interpol(x,nx,xdata,tdata(:,basey+i,basez+j))
            end do
      end do

      ! We then make a 2 dimensionnal interpolation on variables y and z
      fvn_d_quad_3d_interpol=fvn_d_quad_2d_interpol(y,z, &
            3,ydata(basey+1:basey+3),3,zdata(basez+1:basez+3),ttmp)
end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Akima spline interpolation and spline evaluation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Single precision
subroutine fvn_s_akima(n,x,y,br,co)
    implicit none
    integer, intent(in)  :: n
    real, intent(in) :: x(n)
    real, intent(in) :: y(n)
    real, intent(out) :: br(n)
    real, intent(out) :: co(4,n)
    
    real, allocatable :: var(:),z(:)
    real :: wi_1,wi
    integer :: i
    real :: dx,a,b

    ! br is just a copy of x
    br(:)=x(:)
    
    allocate(var(n+3))
    allocate(z(n))
    ! evaluate the variations
    do i=1, n-1
        var(i+2)=(y(i+1)-y(i))/(x(i+1)-x(i))
    end do
    var(n+2)=2.e0*var(n+1)-var(n)
    var(n+3)=2.e0*var(n+2)-var(n+1)
    var(2)=2.e0*var(3)-var(4)
    var(1)=2.e0*var(2)-var(3)
  
    do i = 1, n
    wi_1=abs(var(i+3)-var(i+2))
    wi=abs(var(i+1)-var(i))
    if ((wi_1+wi).eq.0.e0) then
        z(i)=(var(i+2)+var(i+1))/2.e0
    else
        z(i)=(wi_1*var(i+1)+wi*var(i+2))/(wi_1+wi)
    end if
    end do
    
    do i=1, n-1
        dx=x(i+1)-x(i)
        a=(z(i+1)-z(i))*dx                      ! coeff intermediaires pour calcul wd
        b=y(i+1)-y(i)-z(i)*dx                   ! coeff intermediaires pour calcul wd
        co(1,i)=y(i)
        co(2,i)=z(i)
        !co(3,i)=-(a-3.*b)/dx**2                ! méthode wd
        !co(4,i)=(a-2.*b)/dx**3                 ! méthode wd
        co(3,i)=(3.e0*var(i+2)-2.e0*z(i)-z(i+1))/dx   ! méthode JP Moreau
        co(4,i)=(z(i)+z(i+1)-2.e0*var(i+2))/dx**2  !
        ! les coefficients donnés par imsl sont co(3,i)*2 et co(4,i)*6
        ! etrangement la fonction csval corrige et donne la bonne valeur ...
    end do
    co(1,n)=y(n)
    co(2,n)=z(n)
    co(3,n)=0.e0
    co(4,n)=0.e0

    deallocate(z)
    deallocate(var)

end subroutine

! Double precision
subroutine fvn_d_akima(n,x,y,br,co)

    implicit none
    integer, intent(in)  :: n
    double precision, intent(in) :: x(n)
    double precision, intent(in) :: y(n)
    double precision, intent(out) :: br(n)
    double precision, intent(out) :: co(4,n)
    
    double precision, allocatable :: var(:),z(:)
    double precision :: wi_1,wi
    integer :: i
    double precision :: dx,a,b
    
    ! br is just a copy of x
    br(:)=x(:)

    allocate(var(n+3))
    allocate(z(n))
    ! evaluate the variations
    do i=1, n-1
        var(i+2)=(y(i+1)-y(i))/(x(i+1)-x(i))
    end do
    var(n+2)=2.d0*var(n+1)-var(n)
    var(n+3)=2.d0*var(n+2)-var(n+1)
    var(2)=2.d0*var(3)-var(4)
    var(1)=2.d0*var(2)-var(3)
  
    do i = 1, n
    wi_1=dabs(var(i+3)-var(i+2))
    wi=dabs(var(i+1)-var(i))
    if ((wi_1+wi).eq.0.d0) then
        z(i)=(var(i+2)+var(i+1))/2.d0
    else
        z(i)=(wi_1*var(i+1)+wi*var(i+2))/(wi_1+wi)
    end if
    end do
    
    do i=1, n-1
        dx=x(i+1)-x(i)
        a=(z(i+1)-z(i))*dx                      ! coeff intermediaires pour calcul wd
        b=y(i+1)-y(i)-z(i)*dx                   ! coeff intermediaires pour calcul wd
        co(1,i)=y(i)
        co(2,i)=z(i)
        !co(3,i)=-(a-3.*b)/dx**2                ! méthode wd
        !co(4,i)=(a-2.*b)/dx**3                 ! méthode wd
        co(3,i)=(3.d0*var(i+2)-2.d0*z(i)-z(i+1))/dx   ! méthode JP Moreau
        co(4,i)=(z(i)+z(i+1)-2.d0*var(i+2))/dx**2  !
        ! les coefficients donnés par imsl sont co(3,i)*2 et co(4,i)*6
        ! etrangement la fonction csval corrige et donne la bonne valeur ...
    end do
    co(1,n)=y(n)
    co(2,n)=z(n)
    co(3,n)=0.d0
    co(4,n)=0.d0

    deallocate(z)
    deallocate(var)

end subroutine

!
! Single precision spline evaluation
!
function fvn_s_spline_eval(x,n,br,co)
    implicit none
    real, intent(in) :: x           ! x must be br(1)<= x <= br(n+1) otherwise value is extrapolated
    integer, intent(in) :: n        ! number of intervals
    real, intent(in) :: br(n+1)     ! breakpoints
    real, intent(in) :: co(4,n+1)   ! spline coeeficients
    real :: fvn_s_spline_eval
    
    integer :: i
    real :: dx
    
    if (x<=br(1)) then
        i=1
    else if (x>=br(n+1)) then
        i=n
    else
    i=1
    do while(x>=br(i))
        i=i+1
    end do
    i=i-1
    end if
    dx=x-br(i)
    fvn_s_spline_eval=co(1,i)+co(2,i)*dx+co(3,i)*dx**2+co(4,i)*dx**3

end function

! Double precision spline evaluation
function fvn_d_spline_eval(x,n,br,co)
    implicit none
    double precision, intent(in) :: x           ! x must be br(1)<= x <= br(n+1) otherwise value is extrapolated
    integer, intent(in) :: n        ! number of intervals
    double precision, intent(in) :: br(n+1)     ! breakpoints
    double precision, intent(in) :: co(4,n+1)   ! spline coeeficients
    double precision :: fvn_d_spline_eval
    
    integer :: i
    double precision :: dx
    
    
    if (x<=br(1)) then
        i=1
    else if (x>=br(n+1)) then
        i=n
    else
    i=1
    do while(x>=br(i))
        i=i+1
    end do
    i=i-1
    end if
    
    dx=x-br(i)
    fvn_d_spline_eval=co(1,i)+co(2,i)*dx+co(3,i)*dx**2+co(4,i)*dx**3

end function


end module fvn_interpol