! Aquest programa soluciona dos problemes mitjançant la aplicació de algoritmes computacionals
! El primer és tracta d'una equació diferencial
! El segon la resolució d'una integral
! Vicenç Fernàndez Fernàndez
! 26/01/2024

program examen2024
    implicit none

    ! VARIABLES PROBLEMA 1

    ! constants del problema 1
    double precision :: omega, alpha
    ! Funcin vector que conte (z(i),phi(i)), funcout retorna(z(i+1), phi(i+1))
    double precision :: funcin(2), funcout(2)
    ! Variable que marca el temps i el pas utilitzat
    double precision :: t, dt
    ! Energia
    double precision :: Energia, E
    ! utilitzat per els do's
    integer :: i
    ! Funció energia
    external :: E


    ! VARIABLES PROBLEMA 2

    ! Constants del problema 2
    double precision :: x0, gamma
    ! Límits integració I1
    double precision :: x1, x2
    ! Resultat de la integral
    double precision :: resultat
    ! Funcions a integrar
    external :: f, h
    ! vector que emagatzema nombres aleatoris
    double precision, dimension(:), allocatable :: gauss_nums
    ! Resultat de la integral de montecarlo amb el seu error
    double precision :: integral, error
    ! Valor per fer els nombres aleatoris
    integer :: ISEED

    ! Definim els commons utilitzats per els dos problemes
    common/constants/ omega, alpha
    common/constants2/ x0, gamma

    ! dat on ecribim en tot el programa
    open(1,file="Exa-jan-24-res1.dat")

    ! ***** DEFINIM CONSTANTS *****
    omega = 1.d0
    alpha = 2.5d0
    ! *****************************

    !! --------- PROBLEMA 1  ---------
     
    ! Definim les condicions inicials
    funcin(1)=0.3d0 ! z(0)
    funcin(2)=0.d0 ! phi(0)

    dt=(10.d0/1000.d0) ! Definim el pas
    t=0.d0 ! Iniciem el comptador del temps 0

    ! Busquem l'espai fasic per z(0)=0.3
    write(1,*)"#INDEX 0: t, z(t), phi(t). Per z(0)=0.3"
    do i = 0, 1000
        t=dt*dble(i)
        call RungeKutta4order(t,dt,funcin,funcout,3)
        write(1,*) t, funcout(1), funcout(2)
        funcin=funcout
    end do

    write(1,*)
    write(1,*)

    ! Busquem l'espai fasic per z(0)=0.9
    funcin(1)=0.9d0
    funcin(2)=0.d0
    t=0.d0
    write(1,*)"#INDEX 1: t, z(t), phi(t). Per z(0)=0.9"
    do i = 0, 1000
        t=dt*dble(i)
        call RungeKutta4order(t,dt,funcin,funcout,3)
        write(1,*) t, funcout(1), funcout(2)
        funcin=funcout
    end do

    ! ------------------------------------------
    write(1,*)
    write(1,*)
    ! ------------------------------------------

    ! Apartat b
    ! Escrbim l'energia a un fitxer mitançant una funció externa que necesita z i phi en
    ! un temps donat. En el nostre cas utilitzem t=0s i t=10s.

    write(1,*)"#INDEX 2: t, E(t). Per z(0)=0.3"
    funcin(1)=0.3d0 ! z(0)
    funcin(2)=0.d0 ! phi(0)
    ! Calculem la energia per t=0s
    Energia=E(funcin(1),funcin(2))
    write(1,*) 0.d0, Energia

    do i = 0, 1000
        call RungeKutta4order(t,dt,funcin,funcout,3)
        t=dt*dble(i)
        if (dabs(t-10).lt.1.d-8) then ! calculem l'energia per t=10s
            Energia=E(funcout(1),funcout(2))
            write(1,*) t, Energia
        end if
        funcin=funcout
    end do

    write(1,*)
    write(1,*)

    ! Repetim el proces per z(0)=0.9
    write(1,*)"#INDEX 3: t, E(t). Per z(0)=0.9"
    funcin(1)=0.9d0 
    funcin(2)=0.d0
    Energia=E(funcin(1),funcin(2))
    write(1,*) 0.d0, Energia

    do i = 0, 1000
        call RungeKutta4order(t,dt,funcin,funcout,3)
        t=dt*dble(i)
        if (dabs(t-10).lt.1.d-8) then
            Energia=E(funcout(1),funcout(2))
            write(1,*) t, Energia
        end if
        funcin=funcout
    end do
    
    ! ------------------------------------------
    write(1,*)
    write(1,*)
    ! ------------------------------------------


    !! --------- PROBLEMA 2  ---------

    ! ***** DEFINIM CONSTANTS *****
    x0 = 1.d0
    gamma = 0.5d0
    ! *****************************

    ! Apartat a: calcula la integral I1, utilitzant un mètode d'integració tencat. 
    ! En el meu cas he decidit calcular-la amb el mètode simpson 3/8

    ! Límits integració
    x1= 0.d0
    x2 =dacos(-1.d0)

    ! Cridem a la funció
    call simpsontresvuit(x1, x2, 10, f, resultat)
    
    ! Escribim resultat al fitxer
    write(1,*)"#INDEX 4: resultat I1"
    write(1,'(F20.10)') resultat

    ! ------------------------------------------
    write(1,*)
    write(1,*)
    ! ------------------------------------------

    ! Apartat b: calcula I2 mitjançant Montecralo amb sampleig

    ! Generem la base dels nombres aleatoris
    ISEED = (20584933)
    call SRAND(ISEED)

    ! Per crear la distrubció de numeros he escollit una gausiana ja que el metode accepta-rebuig rquereix 
    ! extrems d'integració definits. Aixi doncs fem una gausiana amb valor mitjà 0 i variancia 2.
    
    ! Cridem a la funció box-muller que ens crea un llistat de numeros amb distribució gaussina
    ! I calulem per a cada N el valor de l'integral mitjaçant montecarlo amb smapleig.

    write(1,*)"#INDEX 5: N, integral, error"
    do i=10000,100000,10000
        allocate(gauss_nums(i))
        call boxmuller(i,0.d0,dsqrt(2.d0),gauss_nums)
        call montecarloconnums(gauss_nums, i, integral, error, h)
        write(1,*) i, integral, error
        deallocate(gauss_nums)
    end do


end program examen2024

! --------------------------------------------------------------------------------------

!! --------- FUNCIONS I SUBRUTINES  ---------

! FUNCIONS I SUBRUTINES APARTAT 1

! SUBRUTINA EDO
! Rep com a input (z(0),phi(0)) i retorn (z',phi')
! Es util per poder aplicar RK4
subroutine edo(ndim,t, yin, dyout)
    implicit none
    ! dimensions dels vectors d'entrada i sortida
	integer ndim
    ! 
	double precision :: z, phi, t
    ! yin: conté (z(0),phi(0)) dyout: retorn (z',phi')
	double precision, dimension(ndim) :: yin, dyout
    ! constants del problema
    double precision :: omega, alpha

    common/constants/ omega, alpha
    z = yin(1)
    phi = yin(2)
   
    dyout(1) = - omega * dsqrt(1-(z**2)) * dsin(phi)
    dyout(2) = alpha*z + omega*(z/dsqrt(1-(z**2)))* dcos(phi)
      
end subroutine edo

! RK4
subroutine RungeKutta4order(x,dx,funcin,funcout,nequs)
	implicit none
	integer i, nequs
	! Funcin son dadas inicials i funcout: dades finals
	double precision funcin(nequs), funcout(nequs)
	double precision dx, x, ytemp(nequs), k1(nequs), k2(nequs), k3(nequs), k4(nequs)
   
	! Pas k1
	call EDO(nequs,x,funcin,k1)

	! Pas k2
	do i=1,nequs
	  ytemp(i)= funcin(i) + 0.5d0*dx*k1(i)
	end do
	call EDO(nequs,x+0.5d0*dx, ytemp, k2)
  
	! Pas k3
	do i=1,nequs
	  ytemp(i)= funcin(i) + 0.5d0*dx*k2(i)
	end do
	call EDO(nequs,x+0.5d0*dx,ytemp,k3)
  
	! Pas k4
	do i=1,nequs
	  ytemp(i)= funcin(i) + dx*k3(i)
	end do
	call EDO(nequs,x+dx,ytemp,k4)
	

	do i=1,nequs
	  funcout(i)= funcin(i) + dx/6d0*(k1(i) + 2d0*k2(i) + 2d0*k3(i) + k4(i))
	end do

	return
end 

! FUNCIÓ ENERGIA
! Ens calcula l'energia en un punt en el temps
! Necesita com a input z(t) i phi(t)
double precision function E(z,phi)
    implicit none
    double precision z, phi 
    double precision :: omega, alpha
    common/constants/ omega, alpha
    E = (alpha*(z**2))/2.d0 - dsqrt(1-(z**2))*cos(phi)
    return

end


! FUNCIONS I SUBRUTINES APARTAT 2

! Funció que conte el que s'integra a I1
double precision function f(x)
    implicit none
    ! constants del problema
    double precision :: x0, gamma
    ! punts on avaluem la funció
    double precision :: x

    common/constants2/ x0, gamma

    f = (dsin(x)**2)/((x-x0)**2+gamma**2) 
    return 
end

! Funció que conté h(t)= f(x)/p(x)
! f(x) --> Funció a integrar I2
! p(x) --> Densitat probabilitat gaussiana on extreiem nombres aleatoris
double precision function h(t)
    implicit none
    double precision :: f, p, x0, gamma, pi, t
    common/constants2/x0, gamma
    pi=dacos(-1.d0)

    f = (dsin(t)**2)/((t-x0)**2+gamma**2) 
    p= dexp(-t**2/4.d0)/(2*dsqrt(pi))
    h = f/p
end 

! REGLA DE SIMPSON 3/8. Calcula l'integral de f(x) en a,b tal que:
! ∫ f(x)dx ≈ h/3 [f(a)+2∑f(x_(multiplede3))+3∑f(x_restadenums)+f(b)]

subroutine simpsontresvuit(x1, x2, k, f, resultat)
	implicit none
	! Límts d'integració
	double precision :: x1, x2  
	! Per calcular el numero de subintervals
	integer :: k        
	! Resultat de l'integral
	double precision :: resultat
	! Longitud del subinterval
	double precision :: h
	! I: integer per fer dos, n nombre subintervals
	integer :: i, n
	! xk= a+ n*h
	double precision :: x 
	! Funció a integrar
	double precision :: f 
  
	n = 3**k  ! Calcula el número de subintervals 
  
	h = (x2 - x1) / dble(n)  ! Calcula la longitud de cada subinterval
  
	resultat = f(x1) + f(x2)   ! Suma els valors en los extrems
  
	! Calcul dels punts intermitjos
	do i = 1, n - 1
	  if (mod(i, 3) == 0) then ! Per multiples de 3
		x= x1 + dble(i) * h
		resultat = resultat + 2.0d0 * f(x)
	  else 	! Per la resta de nombres
		x= x1 + dble(i) * h
		resultat = resultat + 3.0d0 * f(x)
	  end if
	end do
  
    resultat = 3.0d0 * h * resultat / 8.0d0 ! Resultat final
	return
end 


!! MONTECARLO AMB SAMPLEIG D'IMPORTÀNCIA
! Donda una funió f(x) i una densitat de probabilitat g(x)
! Per calcular la integral de la funció a partir de nombres distrubits per g(x)
! ∫f(x)= 1/N ∑ f(x)/g(x)
! Tot i que a la nostre subrutina demanem directament h(t)=f(x)/g(x) com a funcio
subroutine montecarloconnums(numeros, N, integral, error, h)
    implicit none
    ! N: nombre de pasos a realitzar i j:utilitzat per fer do's
    integer :: N, j
    ! Resultat de la integral i el seu error
    double precision :: integral, error
    ! h=f(x)/g(x) i t:variable utilitzada per calcular h(t)
    double precision :: t, h
    !Límits d'integració
    double precision :: a, b 
    ! Vector que  conté els nombres aleatoris distribuits segons g(x)
    double precision :: numeros(N)
    ! Funció a integrar
    external :: h

    ! Iniciem els sumatoris a 0
    integral = 0.d0
    error = 0.d0
    do j=1, N
        t= numeros(j) ! Cridem al vector amb els numeros distribuits segons g(x)
        integral= integral + h(t)
        error = error + h(t)**2
    end do

    ! Fem la mitja per donar el valor de l'integral i l'error
    integral = (integral/dble(N))
    error = dsqrt((error/dble(N) - integral**2)/dble(N))

    return
end


! GAUSSIANA NUMS ALEATORIS
! Coneixent el valor mitjà i la desviació podem crear una distribució gaussiana seguint:
! x_n = p*σ + μ
! Amb p= r*Φ, r=sqrt(-2log*rand(0,1)), Φ=2*pi*rand(0,1)
subroutine boxmuller(N,valormitja,variancia,gauss_nums)
	implicit none
	integer N, i
	! Inputs necesaris per saber el tipus de gaussina
	double precision valormitja, variancia
	! Vector on emmagatzmem els resultats
	double precision :: gauss_nums(N)
	! Nombre aleatori entre 0 i 1
	double precision :: x1, x2
	! Altres constants del programa
	double precision:: z1, z2, pi

	pi=dacos(-1.d0)
	do i=1,N-1,2
		! Creem dos numeros aleatoris
		x1=rand() 
		x2=rand()
		! p= r*Φ, r=sqrt(-2log*rand(0,1)), Φ=2*pi*rand(0,1)
		z1=dsqrt(-2.d0*dlog(x1))*dcos(2.d0*pi*x2)
		z2=dsqrt(-2.d0*dlog(x1))*dsin(2.d0*pi*x2)
		! Fem x_n = p*σ + μ i introduïm el resultat a un vector
		gauss_nums(i)=z1*variancia+valormitja
		gauss_nums(i+1)=z1*variancia+valormitja
	enddo

	return
end