! Aquest programa Métodes de Montecarlo (cru, sampleig d’importància), nombres aleatoris
! Vicenç Fernàndez Fernàndez
! 21/11/2023

program practica6
    implicit none
    double precision :: pi, a, b, integral, error, error_real
    double precision :: resultat, L
    double precision :: xhigh, xlow, cotasup, M
    double precision, dimension(:), allocatable:: nums
    integer :: iseed, i, ndades
    external :: I1, p, I2, gxpx

    ! Generem constants que estaran presents a tota la pràctica
    common/dades/pi,L
    pi= dacos(-1.d0)
    L= 50.d0

    ! Generem la seed de nombres aleatoris
    ISEED=20584933
    CALL SRAND(ISEED)

    ! Límit d'integració
    a = 0.d0
    b = 2.d0*pi
    resultat= pi**2*(2.d0*pi**2+(3.d0/2.d0)) ! Resultat real

    ! APARTAT 1
    !  Fent servir el mètode de Montecarlo cru per a calcular la següent integral definida I1
    open(1,file='P6-23-24-res.dat')
    write(1,* ) "#INDEX 0: INEGRAL i ERRORS"
    do i= 10000, 1000000, 10000
        call montecarlocru(i, a, b, integral, error, I1)
        error_real = dabs(resultat-integral) ! Calculem l'error a partir del resultat
        write(1,*) i, integral, error, error_real
    end do


    ! APARTAT B
    !  Genera 1000000 valors aleatoris amb densitat de probabilitat p(x)
    cotasup = 1/L
    xlow = 0.d0
    xhigh = 2*L
    ndades = 1000000
    allocate (nums(ndades))
    call accepta(ndades, nums, xlow, xhigh, cotasup, p)

    write(1,*)
    write(1,*)

    ! APARTAT C
    ! Amb els nombres generats calcula l'integral I2
    write(1,*)"#INDEX 1: I2"
    do i=10000,1000000,10000
        call montecarloconnums(nums, i, integral, error, I2)
        write(1,*) i, integral, error
    end do


    write(1,*)
    write(1,*)

    ! APARTAT 2
    ! Programa un algorisme d’encert/errada per tal de trovar el valor d’I2 i una estimació d'error
    write(1,*)" # INDEX 2: ENCERT/ERRADA"
    a= 0.d0
    b= 2.d0*L
    M= 1.d0/L
    do i = 10000, 300000, 10000
        call encerterrada(i, a, b, M, gxpx, integral, error)
        write(1,*) i, integral, error
    end do


end program practica6

! ---------------------------------------------------------------------------------------------------
! FUNCIONS I SUBRUTINTES !
! ---------------------------------------------------------------------------------------------------

double precision function I1(x)
    implicit none
    double precision :: fx, x
    fx=x**3*dcos(x)**2
    return
end

double precision function p(x)
    implicit none
    double precision:: x, pi, L
    common/dades/pi, L

    p = (1/L)*sin(pi*(x-2*L)/(2*L))**2
    return
end function

double precision function I2(x)
    implicit none
    double precision:: x, pi, L, fx
    common/dades/pi, L

    fx= ((1/L)*sin(pi*(x-2*L)/L)**2)
    return
end


double precision function gxpx(x)
    implicit none 
    double precision :: x, pi, L, fx
    common/dades/pi,L
    fx=  ((1/L)*sin(pi*(x-2*L)/L)**2) * (1/L)*sin(pi*(x-2*L)/(2*L))**2 
    return
end


! SUBRUTINES UTILITZADES

! SUBRUTINA MONTECARLO
! Aquesta subrutina es basa en fer el canvi de l'integral de f(x) de a a b 
! a una altre funcio h(t)=(b − a)f((b − a)t + a) segons una distrubucio uniforme
! El resultat be donat per el promig de h(t).
! HN = 1/N ∑ h(t)

subroutine montecarlocru(N, a, b, integral, error, f)
    implicit none
    integer :: N, i
    double precision :: integral, error, x, f
    double precision :: a, b 
    external :: f

    ! Iniciem els dos sumatoris en 0
    integral = 0.d0
    error = 0.d0

    do i=0, N
        x= (b-a)*rand() + a ! Fem el canvi de variable en les x
        integral= integral + f(x) ! Sumem el valor de cada x
        error = error + f(x)**2
    end do

    ! Fem la mitjana per trobar el valor de l'integral i l'error
    integral = ((b-a)/dble(N)) * integral
    error = dsqrt((error*(b-a)**2/dble(N) - integral**2)/dble(N))

    return
end

! SUBRUITINA ACCEPTAREBUIG
subroutine accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)
    implicit none
    double precision :: numeros(ndades)
    integer :: i, ndades, iseed
    double precision :: xlow, xhigh
    double precision :: x1, x2, x, p, fx, cotasup, funcio
    double precision :: valormitja, variancia, desviacio
    external :: funcio

    do i=1, ndades
        fx = -1
        p = 0
        do while (fx.lt.p)
            x1= rand()
            x2= rand()
            x= x1*(xhigh-xlow)+xlow
            p= cotasup*x2
            fx=funcio(x)
        enddo
        numeros(i)=x
    enddo
end

!! MONTECARLO AMB NUMS ALEATORIS
! Donda una funió f(x) i una densitat de probabilitat g(x)
! Per calcular la integral de la funció a partir de nombres distrubits per g(x)
! ∫f(x)= 1/N ∑ f(x)/g(x)
! Tot i que a la nostre subrutina demanem directament  f(x)/g(x) com a funcio
subroutine montecarloconnums(numeros, N, integral, error, f)
    implicit none
    integer :: N, j
    double precision :: integral, error, x, f
    double precision :: a, b 
    double precision :: numeros(N)
    external :: f

    ! Iniciem els sumatoris a 0
    integral = 0.d0
    error = 0.d0
    do j=1, N
        x= numeros(j) ! Cridem al vector amb els numeros distribuits segons g(x)
        integral= integral + f(x)
        error = error + f(x)**2
    end do

    ! Fem la mktja per donar el valor de l'integral i l'error
    integral = (integral/dble(N))
    error = dsqrt((error/dble(N) - integral**2)/dble(N))

    return
end


subroutine encerterrada(N, a, b, M, fun, integral, sigma)
    implicit none
    integer N, i, ndins
    double precision  a, b, M, p, fun, x, integral, sigma

    ndins = 0

    do i = 1, N

        x = (b-a)*rand(0) + a

        p = M*rand(0)

        if (fun(x).ge.p) then
            ndins = ndins + 1
        endif
    enddo

    integral = M*(b-a)*dble(ndins)/dble(N)
    sigma = (M*(b-a)/DSQRT(dble(N))*DSQRT(((dble(ndins)/dble(N))*(1.d0-dble(ndins)/dble(N)))))

end 
