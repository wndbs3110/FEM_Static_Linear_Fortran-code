
module Beam_Input_3D    
    implicit none
    
    contains
    
    
    subroutine Make_Input
    
    
    
   

    !-----------------------------------------------------------------------------    

    
        
    double PRECISION, parameter :: length=20.0d0
    double PRECISION, parameter :: width=4.0d0, thick=1.0d0
    double PRECISION, parameter :: el_size=0.50d0
    integer, parameter :: n_x=int(length/el_size)+1, n_y=int(width/el_size)+1, n_z=int(thick/el_size)+1
    double PRECISION, parameter :: force=-3.0d4
    double PRECISION, parameter :: Young=2.0d11, Poisson=0.3d0
    
    
    
    !-----------------------------------------------------------------------------
    
    integer,parameter :: n_xnod=n_x
    integer,parameter :: n_ynod=n_y
    integer,parameter :: n_znod=n_z
    
    double precision :: x_coord(n_xnod)
    double precision :: y_coord(n_ynod)
    double precision :: z_coord(n_znod)

    
    integer,parameter :: n_xnel=n_x-1
    integer,parameter :: n_ynel=n_y-1
    integer,parameter :: n_znel=n_z-1
    

    
    integer,parameter :: n_nod=n_xnod*n_ynod*n_znod
    integer,parameter :: n_el= n_xnel*n_ynel*n_znel
    integer :: u(n_nod), v(n_nod), w(n_nod)
    
    double PRECISION :: coord(3*n_nod)
    double PRECISION :: fx(n_el), fy(n_el), fz(n_el), px(n_nod), py(n_nod), pz(n_nod)   
    
    integer :: k
    integer :: i, j, m
    
    fx(:)=0; fy(:)=0; fz(:)=0
    px(:)=0; py(:)=0; pz(:)=0
    py(n_znod*n_ynod*(n_xnod-1)+1:n_nod)=force/DBLE(n_ynod*n_znod)
    
    
    !-----------------------------------------------------------------------------------
    
    open(unit=1, file="input.txt", form="formatted")
    
    
    k=1
    do i=1,n_xnod
        x_coord(i)=(i-1)*dble(length)/dble(n_xnel)
        do j=1,n_ynod
            y_coord(j)=(j-1)*dble(width)/dble(n_ynel)
            do m=1,n_znod
                z_coord(m)=(m-1)*dble(thick)/dble(n_znel)
                coord(3*k)= z_coord(m)
                coord(3*k-2)= x_coord(i) 
                coord(3*k-1)  = y_coord(j)
                k=k+1
            end do
            
        end do
    end do
    
    
    !-----------------------------------------------------------------------------------
    
    u(1:n_ynod*n_znod)=1; u(n_ynod*n_znod+1:n_nod)=0
    v(1:n_ynod*n_znod)=1; v(n_ynod*n_znod+1:n_nod)=0
    w(1:n_ynod*n_znod)=1; w(n_ynod*n_znod+1:n_nod)=0
    !-----------------------------------------------------------------------------------
    
    
    
    write(1,"(a)") "Nodes"
    write(1,"(i8)") n_nod
    write(1,"(a)") "#     x   y   z    u v w   px  py pz"
    do i=1,n_nod
        write(1,"(i5, 3f10.2, 3i5, 3f10.1)") i, coord(3*i-2), coord(3*i-1), coord(3*i), u(i), v(i), w(i), px(i), py(i), pz(i)
    end do
    write(1,*)
    write(1,"(a)") "Elements"
    write(1,"(i5)") n_el
    write(1,"(a)") "#     connectivity   fx  fy  fz"
    
    
    k=1
    do i=1,n_xnod-1
        
        do j=1,n_ynod-1
            do m=1,n_znod-1
        write(1, "(i5, 8i5, 3f7.1)") k,n_ynod*(i-1)*n_znod+n_znod*(j)+m+1, n_ynod*(i-1)*n_znod+n_znod*(j)+m,&     ! 1, 2 node
									n_ynod*(i-1)*n_znod+n_znod*(j-1)+m, n_ynod*(i-1)*n_znod+(j-1)*n_znod+m+1,&		    ! 3, 4 node
                                    n_ynod*i*n_znod+n_znod*j+m+1, n_ynod*(i)*n_znod+n_znod*j+m,&                ! 5, 6 node
                                    n_ynod*(i)*n_znod+n_znod*(j-1)+m, n_ynod*i*n_znod+n_znod*(j-1)+m+1,fx(k), fy(k), fz(k)  ! 7, 8 node
            k=k+1
            end do
        end do
        
    end do
    
    write(1,*)
    write(1,"(a)") "Young"
    write(1,"(es12.2)") Young
    write(1,"(a)") "Poisson"
    write(1,"(es10.4)")Poisson
    
    close(1)
    
    end subroutine Make_Input
    
    end module Beam_Input_3D
    
    
    
    
    
    
    
    