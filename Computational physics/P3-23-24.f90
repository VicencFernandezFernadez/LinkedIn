! Aquest programa genera les dades per la pràctica 3
! Vicenç Fernàndez Fernàndez
! 24/10/2023
program practica3
    implicit none
    double precision :: valorf, valordf, E, pi 
    integer :: i, ndates, nitera
    double precision :: valorsE120(120), funci120(120), dfunci120(120), t(100)
    double precision :: C, B, valorarrel, preci, a = 189.857d0, exc = 0.995086d0 
    external :: fun, fun1
    double precision :: xini, x1, xe, ye, fu, dfu
    

    open(1,file='P3-23-24-res.dat')
    write(1,*)"#INDEX 0: E, D(E) i D'(E)"
    
    
    pi= 4.d0*atan(1.d0)
    ! Ndates= 120
    do i = 0, 119
        valorsE120(i+1) = i * (2d0*pi/120d0) ! vector de 120 punts equiespaiats entre [0, 2π]
        CALL funDE(valorsE120(i+1), valorf)
        funci120(i+1) = valorf ! vector de 120 punts del resultat de D(x)
    end do
    
    
    ! Escribim el resultat en dos arxius E, D(E), D'aprox(E) 
    do ndates = 1,120 
        CALL derivataula(ndates,valorsE120,funci120,dfunci120) 
        write(1,'(F20.12,F20.12,F20.12)') valorsE120(ndates), funci120(ndates), dfunci120(ndates)
    end do

    ! Escrbim dos espais per separar el següent bloc
    write(1,*)
    write(1,*)

    ! Agafem un interval 
    preci = 1.d-9
    C=0.2d0
	B=6.1d0
    call biseccio(fun,C,B,preci,nitera,valorarrel) 
    write(1,*) "#INDEX 1: Emax I D(Emax) "
    call funDE(valorarrel, valorf)
    write(1,'(F20.12,F20.12)') valorarrel, valorf

    write(1,*)
    write(1,*)
    

    write(1,*)'#INDEX 2: T, E, X(E) I Y(E)'
    do i=1,100
		t=i*((2526.5d0)/100.d0)
		xini=(pi/6.d0)
		x1=xini
		call NewtonRaphson(xini,1.d-11,fun1,t,nitera,valorarrel)
		xe=a*(dcos(valorarrel)-exc)
		ye=a*dsqrt(1.d0-(exc**2))*dsin(valorarrel)
		write(1,*) t, valorarrel, xe, ye
	enddo
	
end program practica3



! Subrutina que en es retorna el valor de D(E)= sqrt(x^2+y^2)
! x(E) = a(cos(E) - e) i y(E)=a sqrt(1-e^2)sin(E)
subroutine funDE(E,valorf)
    implicit none
    double precision ::  valorx, valory ! resultat de X(E) i Y(E)
    double precision :: a = 189.857d0, exc = 0.995086d0 
    double precision :: valorf, valordf
    double precision :: E
    
    valorx= a*(dcos(E)-exc)
    valory= a*sqrt(1-exc**2)*dsin(E)
    valorf=sqrt(valorx**2+valory**2)

    return 
end

! Subrutina que ens retorna el valor de F(E) i F'(E)
! F(E)=sin(2E)(1-e^2) - (cos(E)(2-))
subroutine fun(E,F,df)
    implicit none
    double precision :: F, E, exc= 0.995086d0, df
    F=dsin(2.d0*E)*(1.d0-exc**2)-(dcos(E)*(2.d0-exc**2)-exc)*dsin(E)
    df=1
    return 
end

! Subrutina que busca els zeros d'una funció externa
! Utilitzem el metode de bisecció que es basa en que si una funció continua cambia de signe 
! entre dos punts hi ha almenys un zero.
! Necesitem dos punts inicials A,B amb F(a)F(b)<0
! La subrutina s'atura al trobar el zero de manera exacta o amb la precisió desitjada
subroutine biseccio(fun,A,B,preci,nitera,valorarrel)
    implicit none
    double precision :: A,B ! punts inicials de bisecció
    double precision :: C ! punt intermig del interval
    external :: fun ! subroutina que retorna el valor de f(x) i f′(x)
    double precision ::  preci,valorarrel ! precisió desitjada i valor de l'arrel
    integer :: nitera ! nombre d’iteracions per aconseguir la precisió
    double precision :: FA,FB,FC,dif,dfA,dfB,dfC ! variables que ens serveixen per guardar F(x) i df(x)
    integer :: i,maxiter ! maxiter sera el maxim d'itereccions que podrem arribar a fer
    
    ! cridem la subrutina ja que necesitem f(a) i f(b)
    call fun(A,FA,dfA)
    call fun(B,FB,dfB)

    ! comprova que compleixi la condició F(a)F(b)<0
    if (FA*FB.ge.0) then 
        write(*,'(A40,F20.12)') 'La funcio no cambia de signe en aquest interval'
        return
    ! tambe comprovem que els extrems del interval no siguin zeros de la funció
    else if (FA.eq.0) then
        valorarrel=A
        write(*,'(A40,F20.12)') 'La funcio conte una arrel exacta per x=', A
        return
    else if (FB.eq.0) then
        write(*,'(A40,F20.12)') 'La funcio conte una arrel exacta per x=', B
        valorarrel=B
        return
    endif

    ! com l'error disminueix un factor 2 cada itereccio utilitzem aquesta formula per calcular el maxim d'iteracions
    maxiter=nint(log((B-A)/preci)/log(2.d0))+1
    do nitera= 1,maxiter
        C=(A+B)/2.0d0 ! Punt intermig del interval
        ! cridem la subrutina ja que necesitem f(a),f(b) i f(c)
        call fun(A,FA,dfA)
        call fun(B,FB,dfB)
        call fun(C,FC,dfC)
    
       
        ! comprovem que f(c) no sigui arrel
        if (FC.eq.0) then 
            valorarrel=C
            write(*,'(A40,F20.12)') 'La funcio conte una arrel exacta per x=', C
            return
            stop
        ! redefineix el subinterval
        else if (FC*FA.lt.0) then
            B=C 
        else 
            A=C
        endif

        ! si el interval és més petit a la precisio hem arribat a una arrel
        dif=B-A
        if (dif.lt.preci) then
            valorarrel= C
            write(*,"(A40, I10)") 'Nombre iteracions =',nitera
            write(*,'(A40,F20.12)') 'Solucio aproximada x=', C 
            write(*,'(A40,F20.12)') 'Error <',dif
            return
            stop
        end if
    
    end do

    write(*,*)'Hem arribat al màxim de iteracions'
    return
end


! Subrutina que calcula la primera derivada d'una funció dins de l'interval [a,b]
! Per calcular la derivada la subrutina fa:
! si x=a --> f'(a)= (f(a+h)-f(a))/h
! si x ∈ (a,b) --> f'(x)= (f(x + h) − f(x − h))/2h
! si x=b --> f'(b)=(f(b) − f(b − h))/h
! El que fa la subrutina es que rebre dos vectors:
! - Un amb els valors de la variable equiespaiats,valorsx(ndates)
! - L’altre amb els valors corresponents de la funció f(x), funci(ndates),
! retorna un vector amb la derivada calculada numèricament f'(x), dfunci(ndates)
subroutine derivataula(ndates,valorsx,funci,dfunci)
    implicit none
    integer :: ndates ! nombre de arguments que contenen els vectors
    ! VARIABLES D'ENTRADA
    double precision :: valorsx(ndates) ! vector que rep valors de x equiespaiats
    double precision :: funci(ndates) ! vector corresponent a la solució de f(x)
    ! OUTPUT
    double precision :: dfunci(ndates) ! vector amb la derivada calculada numèricament f'(x)
    double precision :: h ! correspon a la distancia entre dues x equiespaiades h= xk+1 - xk
    integer :: i

    h=valorsx(2)-valorsx(1)
    ! Calculem el cas si x=a
    dfunci(1)=(funci(2)-funci(1))/h
    ! Calculem el cas si x ∈ (a,b)
    do i=2,ndates-1
        dfunci(i)= (funci(i+1)-funci(i-1))/2d0*h
    end do
    ! Calculem el cas si x=b
    dfunci(ndates)=(funci(ndates)-funci(ndates-1))/h

    return
end


! Subrutina que busca els zeros d'una funció externa
! Comencem amb un valor x0 a prop de la funció
! Construim un polinomi de Taylor de grau 1 --> T1(x; x0) = f(x0) + (x − x0)f'(x0)
! T1(x1, x0)=0 <-> x1=x0- f(xo)/f'(x0)
! Si  |∆| = |x1 − x0| < precisió hem trobat l'arrel sino redefinim x0=x1
subroutine NewtonRaphson(xini,eps,fun1,t,niter,valorarrel)
	implicit none
	double precision xini, x1, eps, valorarrel, fx, dfx, pmod, t
	external fun1
	integer niter
	character*2 sn
	niter=0
	valorarrel=99999999
	x1=xini
	do while (abs(xini-valorarrel).ge.eps)
		niter=niter+1
!pmod i el if verifiquen quantens iteracions hem fet i ens demanen si volem continuar cada 
! 2000 iteracións per pantalla (escribim SI o NO)
		pmod=mod(niter,2000)
		if (pmod.eq.0) then
			write(*,*)"Aquesta iteració", valorarrel
			write(*,*)"La iteració anterior", xini
			write(*,*)"2000 iteracións. Vols continuar? SI/NO"
1			read*, sn
			if (sn.eq."NO".or.sn.eq."no") then
				return
			else if (sn.ne."SI".or.sn.eq."si") then
				write(*,*)"No es un valor valid, entra SI o NO"
				goto 1
			endif
		endif
		xini=x1
		call fun1(xini,fx,dfx,t)
		valorarrel=xini-fx/dfx
		x1=valorarrel
	enddo
	end



subroutine fun1(x,fu,dfu,t)
	implicit none 
	double precision x, fu, dfu, t, e, a
    integer :: i
	e=0.995086
    
    do i= 1,100
        fu=x+e*dsin(x)-(2*(dacos(-1.d0))/2526.5)*t
        dfu=e*dcos(x)+1
    end do
	return
end


