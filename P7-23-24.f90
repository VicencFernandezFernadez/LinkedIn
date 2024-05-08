! Aquest programa genera totes les dades per a la pràctica 7
! Vicenç Fernàndez Fernàndez
! 28 / 11 / 2023

program practica7
    implicit none
    ! Dades que utilitzarem com a common més wn i Tn que no canvien al llarg de la practica
    double precision :: g, l, pi, m, wn, Tn
    ! Vectors que tindran els resultats de euler i l'algoritme
    double precision, dimension(:), allocatable:: phi, dphi, phi_a, dphi_a
    ! Intervals de temps on avaluarem les subrutines
    double precision :: tmax, tmin, dt
    ! Nombre de iteracions
    integer :: npas, i, j, pasos(4)
    ! Condicions inicials
    double precision  :: phi_0, dphi_0
    ! Variables de l'energia aprtat c
    double precision :: Etotal, Etotal2, Ekin, Epot, cin, pot
    ! Funcio del pendol simple
    external :: fun

    common/dades/g,l,pi,m

    ! Dades utilitzades per tot el programa
    ! _____________________
    g=10.44d0
    l= 1.07d0
    m= 0.98d0
    pi= dacos(-1.d0)
    wn = dsqrt(g/l)
    Tn= 2.d0 * pi/ wn
    ! _________________________

    ! --------------- APARTAT A -----------------
    tmin=0.d0
    tmax= 7.d0* Tn
    npas=1500

    allocate(phi(npas), dphi(npas), phi_a(npas),dphi_a(npas))

    open(1,file='P7-23-24-resf.dat')
   
    call euler(npas,tmin,tmax,phi,dphi,0.025d0,0.d0,dt,fun)

    write(1,*)"#INDEX 0: EULER"
 
    do i=1,npas
        write(1,*) dt*i, phi(i)
    end do

    write(1,*)
    write(1,*)

    call metodesegon(npas, tmin, tmax, phi_a, dphi_a, 0.025d0, 0.d0, dt, fun)

    write(1,*)"#INDEX 1: RK2"
    do i=1,npas
        write(1,*) dt*i, phi_a(i)
    end do

    ! --------------- APARTAT B -----------------
    write(1,*)
    write(1,*)

    phi_0= pi-0.15d0
    call euler(npas,tmin,tmax,phi,dphi,phi_0,0.d0,dt,fun)

    write(1,*)"#INDEX 2: EULER"
 
    do i=1,npas
        write(1,*) dt*i, phi(i), dphi(i)
    end do

    write(1,*)
    write(1,*)

    call metodesegon(npas, tmin, tmax, phi_a, dphi_a, phi_0, 0.d0, dt, fun)
    write(1,*)"#INDEX 3: RK2"
    do i=1,npas
        write(1,*) dt*i, phi_a(i), dphi(i)
    end do

    ! --------------- APARTAT C -----------------
    write(1,*)
    write(1,*)

    phi_0= pi-0.025d0

    call euler(npas,tmin,tmax,phi,dphi,phi_0,0.12d0,dt,fun)
    call metodesegon(npas,tmin,tmax,phi_a,dphi_a,phi_0,0.12d0,dt,fun)

    write(1,*)"#INDEX 4: EPOTEN/ ETOTAL"
   
    do i=1,npas
        Etotal= Epot(phi(i)) + Ekin (dphi(i))
        Etotal2= Epot(phi_a(i)) + Ekin (dphi_a(i))
        write(1,*) dt*i, Ekin (dphi(i)), Etotal,Ekin (dphi_a(i)), Etotal2
    end do
    
    ! --------------- APARTAT D -----------------
    tmin= 0.d0
    tmax= 15.d0*Tn
    npas=6000

    deallocate(phi_a, dphi_a)
    allocate(phi_a(npas),dphi_a(npas))

    dphi_0= 2.d0*dsqrt(g/l)+ 0.05d0
    call metodesegon(npas,tmin,tmax,phi_a,dphi_a,0.d0,dphi_0,dt,fun)

    write(1,*)
    write(1,*)
    write(1,*)"#INDEX 5: Espai fàsic +0.05"
  
    do i=1,npas
        write(1,*) phi_a(i), dphi_a(i)
    end do


    dphi_0= 2.d0*dsqrt(g/l) - 0.05d0
    call metodesegon(npas,tmin,tmax,phi_a,dphi_a,0.d0,dphi_0,dt,fun)
    write(1,*)
    write(1,*)
    write(1,*)"#INDEX 6: Espai fàsic  -0.05"

    do i=1,npas
        write(1,*) phi_a(i), dphi_a(i)
    end do

    ! ------- APARTAT E ----------
    pasos=[300,550,1000,20000]
    tmin= 0.d0
    tmax= 11.d0*Tn

    do i=1,4
        do i=1,4
        write(1,*)
        write(1,*)
        write(1,*)"#INDEX", 6+i, pasos(i)
        deallocate(phi_a, dphi_a)
        allocate(phi_a(pasos(i)),dphi_a(pasos(i)))
        do j=1, pasos(i)
            call metodesegon(pasos(i),tmin,tmax,phi_a,dphi_a,2.87d0,0.0d0,dt,fun)
            cin=Ekin(dphi_a(j))
            pot=Epot(phi_a(j))
            Etotal=cin+pot
            write(1,*) dt*j,Etotal
        end do
    end do
end program practica7

! ---------------------------------------------------------------------------------------------------
! FUNCIONS I SUBRUTINTES !
! ---------------------------------------------------------------------------------------------------


!! Método de Euler
! El metode de Euler per obtenir el resultat d'una EDO es basa en:
! y_(n+1)=y_n+ h·f_n
! Si apliquem per el pendol el metodo es descriu:
! phi_n = phi_(n-1)+ ∆t·dphi_(n-1)
! dphi_n = dphi_(n-1)+ ∆t·(−g/l)·sin(θ)
subroutine euler(npas, tmin, tmax, y, dy, y0, dy0, h, fun)
	implicit none
	double precision :: tmin, tmax, h  ! Temps inicial i final on determinem intervals h en funcio de nombre iteracions 
    double precision :: y0, dy0 ! y_0 i y'_0
    double precision :: y(npas), dy(npas) ! Vector on guardem tots els resultats de y_n i y'_n
	double precision :: pi, m, l, g, fun ! Common/dades
	integer :: npas, i ! Nombre de iteracions del metode
	external :: fun ! Funcio externa a calcular la EDO

	common/dades/g,l,pi,m

	h=(tmax-tmin)/dble(npas) ! Calcul del interval h

    ! Introduim el primer valor y_0 i y'_0 al vector dels resultats
	y(1)= y0 
	dy(1)= dy0

    ! Anem introduint dins dels vectors els resultats obtinguts
	do i=2,npas
		y(i)=y(i-1)+h*dy(i-1)
		dy(i)=dy(i-1)+h*fun(y(i-1))
	enddo
	return 
end


subroutine metodesegon(npas, tmin, tmax, y, dy, y0, dy0, h, fun)
    implicit none 
    double precision :: tmin, tmax, h  ! Temps inicial i final on determinem intervals h en funcio de nombre iteracions 
    double precision :: y0, dy0 ! y_0 i y'_0
    double precision :: y(npas), dy(npas) ! Vector on guardem tots els resultats de y_n i y'_n
    double precision :: ky1, ky2, kdy1, kdy2
	double precision :: pi, m, l, g, fun ! Common/dades
	integer :: npas, i ! Nombre de iteracions del metode
	external :: fun ! Funcio externa a calcular la EDO

	common/dades/g,l,pi,m

	h=(tmax-tmin)/dble(npas) ! Calcul del interval h

    ! Introduim el primer valor y_0 i y'_0 al vector dels resultats
	y(1)= y0 
	dy(1)= dy0

    ! Anem introduint dins dels vectors els resultats obtinguts
	do i=2,npas
        ky1= fun(y(i-1))
        kdy1 = dy(i-1)

        ky2 = dy(i-1) + (h/2.d0)*ky1
        kdy2 = fun(y(i-1) + (h/2.d0)*kdy1)

		y(i)=y(i-1)+h*ky2
		dy(i)=dy(i-1)+h*kdy2
	enddo
	return 
end

double precision function fun(x)
	implicit none
	double precision x, fu
	double precision pi, m, l, g
	common/dades/g,l,pi,m
	fu=-(g/l)*dsin(x)
	return
end


double precision function Ekin(dphi)
    implicit none 
    double precision k, dphi
    double precision pi, m, l, g
    common/dades/g,l,pi,m
    K= (1.d0/2.d0)*m*(dphi**2)*l**2
    return
end

double precision function Epot(phi)
    implicit none
    double precision V, phi
    double precision pi, m, l, g
    common/dades/g,l,pi,m
    V=-m*g*l*dcos(phi)
    return
end
