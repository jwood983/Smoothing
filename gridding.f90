! This program will smooth a set of data (read in from a file) onto a grid
!  with a size that depends on the desired resolution. The file is expected
!  to be a 2-column data (x,y) with values in the range [0,1]. Values out
!  of that range will not be included in the smoothing (and may cause a
!  segmentation fault near the boundaries).

! preprocessor directives (quartic = 1, epan = 2, gaussian = 3, tricubic = 4)
#define WEIGHT 1

program gridding
    implicit none
    integer, parameter :: wp = selected_real_kind(15,307)
    real(wp) :: E
    real(wp), allocatable, dimension(:,:) :: Map
    real(wp), allocatable, dimension(:) :: x, y, xgrid, ygrid
    real(wp), dimension(10) :: resolution = -1.0
    real(wp) :: h, dij, u, v, xp, yp, thick, radius
    integer :: i, j, k, n, ii, jj, numpart, numgrid, ierr, lun, fileread, partStep, numres
    character(len=100) :: filename

    lun = 10

! open input file & error check
    namelist / particles / filename, resolution
    open(unit=lun,file="init.nml",iostat=ierr)
    if(ierr/=0) then
        print*,"Unable to open 'init.nml', does it exist?"
        stop
    endif
    read(lun,nml=particles,iostat=ierr)
    if(ierr/=0) then
       print *,"Unable to read 'particles' namelist, does it exist in file?"
    endif
! error check resolution
    if(any(resolution > 0.5)) then
       print *,"Error: resolution above 0.5 causes segmentation fault error"
       stop
    endif
    close(lun)

! open file with particle position data
    open(unit=lun,file=filename,iostat=ierr)
    if(ierr/=0) then
       print *,"Error, unable to read file '",filename,"'"
       stop
    endif
    
! count number of particles
    i=1
    do
       read(lun,*,iostat=ierr) xp,yp,E
       if(ierr/=0) exit
       i = i +1
    enddo
    rewind(lun)
    numpart = i-1
    
! allocate particle positions and then read them in
    allocate(x(numpart),y(numpart))
    do i=1,numpart
       read(lun,*) x(i),y(i)
    enddo

    numres = count(resolution > 0.0)
! run through resolutions
    do k=1,numres
        h = resolution(k)
        numgrid = int(2d0/h)
! allocate grid space
        print '(a,f8.5,a,i0,a)',"Running h=",h," --> grid size of ",numgrid,"^2"
        call AllocateGrid(numgrid,Map,xgrid,ygrid)

! step size (set > 1  for more diffuse display)
        partStep = 1

! loop through the particles from the file
        do n=1,numpart,partStep
! find the closest cells to the particle position
            xp = x(n); yp = y(n)
            ii = FindBox(numgrid,xp,xgrid)
            jj = FindBox(numgrid,yp,ygrid)

! loop through the adjacent cells of centroid (ii,jj) and measure "leakage"
            do j=jj-1,jj+1
                do i=ii-1,ii+1
                    u = xp - xgrid(i)
                    v = yp - ygrid(j)
                    dij = sqrt(u*u + v*v)
                    if(dij <= h) then
                        map(i,j) = map(i,j) + kernel(dij,h)
                    endif
                enddo !- i
            enddo !- j
        enddo !- n

! write the map to a file
        call WriteOut(numgrid,xgrid,ygrid,map,h)

    enddo !- k

 contains

    !> @brief allocate grid positions and map
    subroutine AllocateGrid(n,M,x,y)
        integer, intent(in) :: n
        real(wp), allocatable,intent(out) :: M(:,:), x(:), y(:)
        real(wp) :: xi, yi, dx
        integer :: ierr, i

! deallocate if necessary
        if(allocated(M)) deallocate(M)
        if(allocated(x)) deallocate(x)
        if(allocated(y)) deallocate(y)
        allocate(M(n,n),x(n),y(n),stat=ierr)
        if(ierr/=0) then
            print *,"unable to allocate M,x,or y"
            stop
        endif

        xi = -0.5
        yi = -0.5
        dx = 2.0/real(n+1)
        do i=1,n
            x(i) = xi + (real(i) - 0.5)*dx
            y(i) = yi + (real(i) - 0.5)*dx
        enddo
        M(:,:) = 0.0
    end subroutine AllocateGrid

    !> @brief determine closest cell to point
    pure function FindBox(n,rp,rg) result(i)
        integer, intent(in) :: n
        real(wp), intent(in) :: rp
        real(wp), dimension(n), intent(in) :: rg
        integer:: i
        real(wp) :: r, rmin
        integer :: j,jcell
        rmin = 100d0
        do j=1,n
            r = abs(rp - rg(j))
            if(r < rmin) then
                rmin = r
                jcell = j
            endif
        enddo
        i = jcell
    end function FindBox

    !> @brief tricubic interpolation
    !! @info 0.864... is 70/81
    real(wp) function kernel(d,h)
        real(wp), intent(in) :: d, h
        real(wp) :: u, coef, a, b
        u = d/h
#if WEIGHT == 1
! Quartic biweight
        a = 0.9375
        coef = a/h
        b = 1d0 - u*u
        kernel = coef*b*b
#elif WEIGHT == 2
! Epanechnikov
        a = 0.75
        coef = a/h
        kernel = coef*(1d0 - t*t/5d0)
#elif WEIGHT == 3
! Gaussian
        a = 0.39894228 ! 1/sqrt(2*pi)
        coef = a/h
        kernel = coef*exp(-0.5*t*t)
#elif WEIGHT == 4
! Tricubic
        a = 0.864197531
        coef = a/h
        b = 1d0 - u*u*u
        kernel = coef*b*b*b
#else
#error "Unknown kernel"
#endif
    end function kernel

    !> @brief write data to file
    subroutine WriteOut(ng,x,y,M,h)
        integer, intent(in) :: ng
        real(wp), dimension(ng),intent(in) :: x,y
        real(wp), dimension(ng,ng), intent(in) :: M
        real(wp), intent(in) :: h
        integer :: i,j, ierr, lun
        character(len=80) :: filename
        character(len=8) :: nchar
        
        write(nchar,'(f8.5)') h
        filename = 'grid_'//trim(adjustl(nchar))//'.dat'

        lun = 10

! open output file & error check
        open(unit=lun,file=trim(filename),iostat=ierr)
        if(ierr/=0) then
            print *,"Error:unable to open file",filename
            stop
        endif

! this writes to a gnuplot-readable format
        do j=1,ng
            do i=1,ng
                write(lun,*) x(i),y(j),M(i,j)
            enddo
            write(lun,*) " "
        enddo
        close(lun)
    end subroutine WriteOut


end program gridding
