! Aquest programa genera les dades per a la pràctica 8
! 12/12/2023
! Vicenç Fernàndez Fernàndez

program practica8
    implicit none
    ! constants del programa
    double precision :: dx, xf, x0, L, y0, dy0, E, x,  beta, alfa, delta, V0
    ! variable apartat 1
    double precision ::  Ei(4)
    double precision, dimension(:), allocatable :: yin, yout, phi1
    double precision, dimension(:), allocatable :: phi, valors_dy, yout_norm, phiout
    ! variables per la sub de tir
    double precision :: E3, phi3, error, E1(3), E2(3)
    ! variables per l'apartat 3
    double precision :: prob1, prob2, prob3 !probabilitatz
    double precision, dimension(:), allocatable :: yout_norm1, yout_norm2, yout_norm3 !vectors normalitzts
    ! integrs del programa
    integer N, i, j, nequ
    !Subrutines utilitzades
    external EDO, RungeKutta4order, inte 

    common/dades/ E, alfa, delta, beta, V0

    ! constants que son utils per tot el programa

    ! -------------------------------------
    alfa= 2.d0
    delta= 0.05d0
    beta=0 ! canvia en el aprtat 3
    V0=-20.d0

    L= 8.d0 !longitud de la caixa
    x0=-L/2 ! intervals 
    xf=L/2
    !Condicions inicials
    y0=0.d0 
    dy0=2*10.d0**(-6.d0)

    N=400 ! nombre de pasos
    nequ=2 ! nombre d'equacions
    dx = (xf-x0)/dble(N) ! interval
    ! -------------------------------------


    ! ____________ APARTAT 1 ____________
    ! Aquest apartat integra l'equació i obte les dades per fer una figura de les solucions
    ! en el interval [-L/2 : L/14]

    open(1,file="dades.dat") ! obrim el dat on posarem els resultats
    
    allocate(yin(nequ), yout(nequ), phi(N+1), valors_dy(N+1))
   
    Ei=[-21.d0,-20.5d0,-14.d0, -13.d0] 

    do i=1,4
        E=Ei(i)
        call inte(N,x0,xf,y0,dy0,nequ,phi,valors_dy) ! integrem l'equció
        do j=1, N+1
            x=(j-1)*dx+ x0
            if (x.le.L/14.d0) then ! si la x esta l'interval afegim la seva solució
                write(1,*) x, phi(j)
            endif
        end do
        write(1,*)
        write(1,*)
    end do


! __________  APARTAT 2 ________________
! Apartat a: mitjançant el metòde de tir calcula els primers tres autovalors
! I amagatzema les dades per fer una figura  mostrant la convergencia

! Apartat b: calcula els autovectors normalitzats a partir dels aurovalors

    error= 1.d-6 ! error del metode
    E1=[-21.d0,-14.d0,-8.d0]
    E2=[-20.5d0,-13.d0, -7.5d0]

    allocate(yout_norm(N), phiout(N+1))

    open(3,file='normal.dat') ! dat on escribim els autovectors normalitzats

    do i=1,3
        ! apartat a
        call tir(N, error, E1(i), E2(i), E3, nequ, phiout, y0, dy0,x0, xf)
        write(2,*)
        write(2,*)
        ! aprtat b
        call normalitzacio(N,dx,phiout,yout_norm)
        do j=1,n+1
            x=(j-1)*dx+ x0
            write(3,*) x, yout_norm(j)
        end do
        write(3,*)
        write(3,*)
    end do


! __________ APARTAT 3 ___________
! Compara els autovectors normalitzats per diferents valors de beta en l'estat fonamental

    open(4,file='normabeta.dat')

    E1=-21.d0
    E2=-20.d0

    allocate(yout_norm1(N), yout_norm2(N), yout_norm3(N))
    
    beta=0.d0
    call tir(N, error, E1, E2, E3, nequ, phiout, y0, dy0,x0, xf)
    call normalitzacio(N,dx,phiout,yout_norm1)
    do j=1,n+1
        x=(j-1)*dx+ x0
        write(4,*) x, yout_norm1(j)       
    end do
        
    write(4,*)
    write(4,*)

    beta=1.d0
    call tir(N, error, E1, E2, E3, nequ, phiout, y0, dy0,x0, xf)
    call normalitzacio(N,dx,phiout,yout_norm2)
    do j=1,n+1
        x=(j-1)*dx+ x0
        write(4,*) x, yout_norm2(j)       
    end do
        
    write(4,*)
    write(4,*)

    beta=5.d0
    call tir(N, error, E1, E2, E3, nequ, phiout, y0, dy0,x0, xf)
    call normalitzacio(N,dx,phiout,yout_norm3)
    do j=1,n+1
        x=(j-1)*dx+ x0
        write(4,*) x, yout_norm3(j)     
    end do
        
    write(4,*)
    write(4,*)
    

    ! apartat b
    ! calcula la probabilitat de que un electro es trobi en x[-1.3:1.3]
    open(5,file="P8-23-24-res1.dat")
    prob1 = 0.d0
    prob2 = 0.d0
    prob3 = 0.d0
    do i=1, N+1
        x=(i-1)*dx+x0
        if (x.le.1.3d0.and.x.ge.-1.3d0) then
            prob1=prob1+dx*yout_norm1(i)**2.d0
            prob2=prob2+dx*yout_norm2(i)**2.d0
            prob3=prob3+dx*yout_norm3(i)**2.d0
        end if
    end do
    write(5,*)'Probabilitat per beta = 0 eVA^{-2} ', prob1*100, '%'
    write(5,*)'Probabilitat per beta = 1 eVA^{-2} ', prob2*100, '%'
    write(5,*)'Probabilitat per beta = 5 eVA^{-2} ', prob3*100, '%'
   

end program practica8

! ---------------------------------------------------------------------------------------------------
! SUBRUTINES

! EDO
subroutine EDO(nequ,x,yin,dyout)
    implicit none
    integer nequ
    double precision x, yin(nequ), dyout(nequ), V, E, alfa, delta, beta, V0
    common/dades/ E, alfa, delta, beta, V0


    dyout(1)=yin(2)
    V=(V0* dsinh(alfa/delta)/(dcosh(alfa/delta)+dcosh(x/delta))) + beta*x**2
    

    dyout(2) = (2.d0/7.6199d0)*(V-E)*yin(1)

end 

! RK4
subroutine RungeKutta4order(x,dx,funcin,funcout,nequs)
	implicit none
	integer i, nequs
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

! INTEGRACIÓ
subroutine inte(N,x0, xf, y0, dy0, nequ, valors_y, valors_dy)
    implicit none
    integer i, N, nequ
    double precision E, x0, xf, dx, y0, dy0, x
    double precision, dimension(nequ) :: yin, yout
    double precision, dimension(N+1) :: valors_y, valors_dy
    
  

    dx=(xf-x0)/dble(N)
    valors_y(1)=y0
    valors_dy(1)=dy0
    yin(1)=y0
    yin(2)=dy0
   

    do i=1, N
        x=(i-1)*dx+x0
        call RungeKutta4order(x,dx,yin,yout,nequ)
        yin=yout
        valors_y(i+1)=yin(1)
        valors_dy(i+1)=yin(2)
    end do

end subroutine inte

! METODE TIR
subroutine tir(N, error,E1, E2, E3, nequ, phiout, y0, dy0,x0, xf)
    implicit none
    double precision error, E1, E2, E3, phi3, x, y0, dy0, dx, E, x0, xf, alfa, delta, beta, V0
    integer i, N, nequ, count
    external RungeKutta4order, EDO
    double precision, dimension(nequ) :: yin, yout
    double precision, dimension(N+1) ::valors_phi1, valors_phi2, phiout
    common/dades/ E, alfa, delta, beta, V0

    phi3=100.d0*error
    valors_phi1(1)=y0
    valors_phi2(1)=y0
    phiout(1)=y0
    dx=(xf-x0)/dble(N)

    count=1
    do while (abs(phi3).gt.(error)) 
        yin(1)=y0
        yin(2)=dy0
        
        E=E1
        do i=1, N
            x=(i-1)*dx+x0
            
            call RungeKutta4order(x,dx,yin,yout,nequ)
            yin=yout
            valors_phi1(i+1)=yin(1)
        
        enddo

        yin(1)=y0
        yin(2)=dy0
    

        E=E2
        do i=1, N
            x=(i-1)*dx+x0
            call RungeKutta4order(x,dx,yin,yout,nequ)
            yin=yout
            valors_phi2(i+1)=yin(1)
        
        enddo
      
        E3=(E1*valors_phi2(N+1)-E2*valors_phi1(N+1))/(valors_phi2(N+1)-valors_phi1(N+1))

        open(2,file='tir.dat')
        

        yin(1)=y0
        yin(2)=dy0
    
        E=E3
        do i=1, N
            x=(i-1)*dx+x0
            call RungeKutta4order(x,dx,yin,yout,nequ)
            yin=yout
            phiout(i+1)=yin(1)
        enddo
        
        
        phi3=phiout(N+1)
        write(2,*) count, E3, phi3

        E1=E2
        E2=E3
        count=count+1
     

    enddo

end subroutine tir

! NORMALITZACIÓ
subroutine normalitzacio(npasos,dx,yin,yout_norm)
	implicit none
	integer npasos, i
	double precision dx, yin(npasos), yout_norm(npasos), area

	area=0.d0
		do i=1,npasos
			area= area + yin(i)**2*dx
		end do

		do i=1,npasos
			yout_norm(i) = yin(i)/DSQRT(area)
		end do
	return
end 