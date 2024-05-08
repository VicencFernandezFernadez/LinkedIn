! Aquest programa genera la Resolució de EDP, equacions el·líptiques, equació de Poisson, equació de la calor
! Vicenç Fernàndez Fernàndez
! 19/12/2023

program practica9
    implicit none
    ! Constants del programa
    double precision :: Lx, Ly, h, w, preci
    integer :: Ny, Nx
    ! Coordendes
    double precision :: x, y
    ! Condicions de contorn
    double precision :: Tx0, T0y, TxLy, TLxy
    ! Matriu on emmagatzamem els canvis de temperatura
    double precision, dimension(:,:), allocatable :: Temp
    ! Variables per l'apartat 3
    double precision :: xi, yi, Ti(3)
    ! Integer per crear els bucles do
    integer :: i, j
    ! Creem els comuns per introduirlos dins de la subrutina
    common/longituds/ Lx, Ly
    common/temperatura/ Tx0, T0y, TxLy, TLxy, w

    ! CONSTANTS DEL PROGRAMA
    Lx = 33.5d0
    Ly = 45.5d0
    h = 0.25d0   
    w = 1.35d0
    Nx = int(Lx/h) + 1
    Ny = int(Ly/h) + 1
    preci = 1.d-4
    ! -----------------------

    ! CONDICIONS DE CONTORN
    Tx0 = 0.5d0 ! T en (x,0)
    T0y = 17.d0 ! T en (0,y)
    TxLy = 11.2d0 ! T en (x,Ly)
    TLxy = 25.3d0 ! T en (Lx,y)
    ! -----------------------


    ! ________ APARTAT 3_______
    ! Estudia la convergencia en el punt (7.5,23.5) amb els dos metodes
    ! Amb tres matrius amb Tinterior=10,120,1040ºC

    ! Punt d'estudi
    xi = 7.5d0
    yi = 23.5d0

    Ti = [10.d0, 120.d0, 1040.d0] ! Creem un vector que conté les 3 temperatures a estudiar

    allocate(Temp(Nx,Ny))

    open(1,file='pre9-23-24-res1.dat')

    ! Metode de Gauss-Seidel per les tres temperatures incials
    do i=1,3
        write(1,*) 'Metode Gauss-Seidel'
        call Poisson(2, Nx, Ny, h, Ti(i), preci, xi, yi, 1, Temp)
        write(1,*)
        write(1,*)
    end do

    ! Metode de Sobrerelaxació per les tres temperatures incials
    do i=1,3
        write(1,*) 'Metode Sobrerelaxació'
        call Poisson(3, Nx, Ny, h, Ti(i), preci, xi, yi, 1, Temp)
        write(1,*)
        write(1,*)
    end do
    
    close(1)

    ! ________ APARTAT 4_______
    ! Genera un dat que emagatzema la temperatura a cada punt de la malla

    open(2,file='mapa_temperatures.dat')

    call Poisson(3, Nx, Ny, h, Ti(2), preci, xi, yi, 0, Temp)

    do i = 1, Nx
        x = (i-1)*h
        do j = 1, Ny
            y = (j-1)*h
            write(2,*) x, y, Temp(i,j)
        end do
        write(2,*)
    end do


    close(2)

    !________ APARTAT 5_______
    ! Genera un dat que emagatzema la temperatura a cada punt de la malla sense fonts de calor

    open(3,file='mapa_temperatures_no_fonts.dat')

    ! Cridem a la funció sense fonts de calor
    call poisson_no_fonts(3, Nx, Ny, h, Ti(2), preci, xi, yi, 0, Temp)

    do i = 1, Nx
        x = (i-1)*h
        do j = 1, Ny
            y = (j-1)*h
            write(3,*) x, y, Temp(i,j)
        end do
        write(3,*)
    end do


    close(3)


end program practica9

! ---------------------------------------------------------------------------------------------------
! SUBRUTINES

! FONTS DE CALOR
subroutine rho(x, y, rho_k)
    implicit 
    ! Coordenades de la malla
    double precision :: x, y, r 
    ! Valor de cada fogó i la suma total
    double precision :: rho02, rho03, rho1, rho_k
    ! Constants per la rho2 i rho3
    double precision :: rho2, rho3

    !Definicio de rho1
    if ( x.ge.18.and.x.le.22.and.y.ge.29.and.y.le.35 ) then
        rho1 = 3.d0
    else
        rho1  = 0d0
    end if

    !Definicio de rho2
    rho02 = 10.d0
    r = dsqrt((x-8.d0)**2 + (y-22.5d0)**2)
    rho2 = rho02*dexp(-((r-5d0) / 0.3d0)**2)
    
    ! Definicio de rho3
    rho03 = 6.d0
    r = dsqrt((x-22.d0)**2 + (y-10.5d0)**2)
    rho3 = rho03*dexp(-((r-4d0) / 0.8d0)**2)

    ! Rho total
    rho_k = rho1 + rho2 + rho3

    return
end subroutine rho


! SUBRUTINA POISSON
! Subrutina que resol la equacició amb un dels tres metodes depenent de icontrol. Si:
! icontrol = 1 --> Metode Jacobi
! icontrol = 2 --> Metode Gauss-Seidel
! icontrol = 2 --> Metode sobrerelaxació
subroutine poisson(icontrol, Nx, Ny, h, Ti, epsilon, xi, yi, escriu, T)
    implicit none
    double precision:: T0y, TxLy, TLxy, Tx0, h, Ti, error, epsilon, x, y, rho_k, w, xi, yi
    double precision:: T(Nx,Ny), Tk(Nx,Ny)
    integer:: icontrol, escriu, Nx, Ny, i, j, k
    common/temperatura/Tx0, T0y, TxLy, TLxy, w

    ! CONDICIONS DE CONTORN
    do i = 1, Nx
        T(i, 1) = Tx0
        T(i, Ny) = TxLy
    end do
 
    do i = 1, Ny
        T(1, i) = T0y
        T(Nx, i) = TLxy
    end do

    ! AFEGIM LES T inicials
    do i = 2, Nx-1
        do j = 2, Ny-1
            T(i, j) = Ti
        end do
    end do

 
    error = 1.d0

    ! INICIEM EL NOMBRE DE ITERACIONS
    k = 0   

    ! INICIEM UN BUCLE
    ! Fins que l'error sigui menor a la precisió desitjada
    do while (error>epsilon)
        Tk = T      !La matriu Tk sera T  anterior 

        ! Recorrem tota la matriu excepte on hem definit les CC
        do i = 2, Nx-1
            x = h*(i-1)    
            do j = 2, Ny-1
                y = h*(j-1)

                !Per cada punt trobem el valor de rho
                call rho(x, y, rho_k) 

                ! ESCOLLIM EL METODE A UTILITZAR
                if (icontrol.eq.1) then     !Jacobi

                    T(i,j) = (Tk(i+1,j) + Tk(i-1,j) + Tk(i,j+1) + Tk(i, j-1) + rho_k*h**2)/4.d0

                else if (icontrol.eq.2) then !Gauss-Seidel

                    T(i,j) = (Tk(i+1,j) + T(i-1,j) + Tk(i,j+1) + T(i, j-1) + rho_k*h**2)/4.d0

                else if (icontrol.eq.3) then    !Sobrerelaxació

                    T(i,j) = Tk(i,j) + w*(Tk(i+1,j) + T(i-1,j) + Tk(i,j+1) + T(i, j-1) + rho_k*h**2 - 4*Tk(i,j))/4.d0
                
                endif

                !Si escriu=1 guardarem els valors de T en el punt (xi,yi)
                if ( escriu.eq.1 ) then
                    if ((x.eq.xi).and.(y.eq.yi)) then
                        write(1,*) k, T(i,j)
                    end if
                end if
                
            end do
        end do

        
        error = maxval(abs(T - Tk))

        ! Sumem nova iteració
        k = k+1
    end do

    return

end subroutine

! SUBRUTINA POISSON_NO_FONTS
! Subrutina que resol la equacició amb un dels tres metodes depenent de icontrol. Si:
! icontrol = 1 --> Metode Jacobi
! icontrol = 2 --> Metode Gauss-Seidel
! icontrol = 2 --> Metode sobrerelaxació
! A més no te en compte les fonts de calor
subroutine poisson_no_fonts(icontrol, Nx, Ny, h, Ti, epsilon, xi, yi, escriu, T)
    implicit none
    double precision:: T0y, TxLy, TLxy, Tx0, h, Ti, error, epsilon, x, y, rho_k, w, xi, yi
    double precision:: T(Nx,Ny), Tk(Nx,Ny)
    integer:: icontrol, escriu, Nx, Ny, i, j, k
    common/temperatura/Tx0, T0y, TxLy, TLxy, w

    ! CONDICIONS DE CONTORN
    do i = 1, Nx
        T(i, 1) = Tx0
        T(i, Ny) = TxLy
    end do
   
    do i = 1, Ny
        T(1, i) = T0y
        T(Nx, i) = TLxy
    end do

    ! AFEGIM LES T inicials
    do i = 2, Nx-1
        do j = 2, Ny-1
            T(i, j) = Ti
        end do
    end do

    
    error = 1.d0

    ! INICIEM EL NOMBRE DE ITERACIONS
    k = 0   

    ! INICIEM UN BUCLE
    ! Fins que l'error sigui menor a la precisió desitjada
    do while (error>epsilon)
        Tk = T      !La matriu Tk sera T de la iter anterior (i-1)

       ! Recorrem tota la matriu excepte on hem definit les CC
        do i = 2, Nx-1
            x = h*(i-1)    
            do j = 2, Ny-1
                y = h*(j-1)

                
                ! ESCOLLIM EL METODE A UTILITZAR
                if (icontrol.eq.1) then     !Jacobi

                    T(i,j) = (Tk(i+1,j) + Tk(i-1,j) + Tk(i,j+1) + Tk(i, j-1) )/4.d0

                else if (icontrol.eq.2) then !Gauss-Seidel

                    T(i,j) = (Tk(i+1,j) + T(i-1,j) + Tk(i,j+1) + T(i, j-1))/4.d0

                else if (icontrol.eq.3) then    !Sobrerelaxació

                    T(i,j) = Tk(i,j) + w*(Tk(i+1,j) + T(i-1,j) + Tk(i,j+1) + T(i, j-1) - 4*Tk(i,j))/4.d0
                
                endif

                !Si escriu=1 guardarem els valors de T en el punt (xi,yi)
                if ( escriu.eq.1 ) then
                    if ((x.eq.xi).and.(y.eq.yi)) then
                        write(1,*) k, T(i,j)
                    end if
                end if
                
            end do
        end do

       
        error = maxval(abs(T - Tk))

        ! Sumem nova iteració
        k = k+1
    end do

    return

end subroutine
