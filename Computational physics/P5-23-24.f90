! Aquesta pràctica genera totes les dades per la P5
! Vicenç Fernàndez Fernàndez
! 14/11/2023

program practica4
    implicit none
    double precision :: xlow, xhigh, cotasup ! intervals de la funció i cotasuperior
    double precision :: l, pi, boxsize, sig
    integer :: ndades, nbox, i,  iseed
    double precision, dimension(:), allocatable::numeros, xhis, vhis, errhis
    external :: funcio, funcio2
    double precision :: x1, x2, resultat
    double precision :: variancia, valormitja, desviacio

    ISEED=20584933
    CALL SRAND(ISEED)

    ! EXERCICI 1
    ! PART A
    ! Genera un histograma normalitzat amb 80 caixes entre [-3L,3L] a partir de 40000 valors de x

    l=3.d0
	pi=dacos(-1d0)
	cotasup=2.d0

	ndades=40000 !nombre de dades
	nbox=80 !nombre caixes
	xlow=-3*l 
	xhigh=3*l

    ! introduim els valors dintre dels vectors amb alocate
	allocate(numeros(ndades),xhis(nbox),vhis(nbox),errhis(nbox))

    ! cridem a les subrutines
	call accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)
	call histograma(ndades,numeros,xlow,xhigh,nbox,xhis,vhis,errhis,boxsize)

    ! obrim i escribim al dat les dades per fer l'histograma
    open(1,file="P5-23-24-res.dat")
    write(1,*)"#INDEX 0: PUNTS, VALORS I ERRORS. FIG 1"
    do i=1,nbox
		write(1,*) xhis(i), vhis(i), errhis(i)
	end do


    write(1,*)
    write(1,*)

    deallocate(numeros,errhis,vhis,xhis)

    ! PART B
    ! Aquest apartat calcula la probabilitat de trobar un electro entre -3L i L/2
    write(1,*) "INDEX 1: INTEGRALS"
    x1= -3.d0*l
    x2= L/2.d0
    call simpsontresvuit(x1, x2, 10, funcio, resultat)
    write(1,*) 'P1=', resultat

    ! Comprovem si la funció esta normalizada
    x1= -3.d0*l
    x2= 3.d0*l
    call simpsontresvuit(x1, x2, 10, funcio, resultat)
    write(1,*) 'La funció esta normalitzada ja que el resultat és', resultat


    write(1,*)
    write(1,*)

    ! EXERCICI 2
    ! Genera un histograma normalitzat amb 70 caixes entre [-3SIG,3SIG] a partir de 10000 valors de x
    sig=4.d0
	pi=dacos(-1d0)
	cotasup=2.d0

	ndades=10000
	nbox=70
	xlow=-3*sig
	xhigh=3*sig

    allocate(numeros(ndades),xhis(nbox),vhis(nbox),errhis(nbox))
	call accepta(ndades,numeros,xlow,xhigh,cotasup,funcio2)
	call histograma(ndades,numeros,xlow,xhigh,nbox,xhis,vhis,errhis,boxsize)

    write(1,*)"#INDEX 2: PUNTS, VALORS I ERRORS. FIG 2"
    do i=1,nbox
		write(1,*) xhis(i), vhis(i), errhis(i)
	end do


    write(1,*)
    write(1,*)

    ! PART B
    ! Calcula el valor mitja, variancia i desviació estandard
    write(1,*)"INDEX 3: VALOR  MITJÀ, VARIÀNCIA I DESVIACIÓ ESTÀNDARD"
    valormitja = 0
    do i= 1, ndades
        valormitja = valormitja + numeros(i)
    end do
    valormitja = valormitja/dble(ndades)
    write(1,*)'Valor mitjà =', valormitja

    variancia = 0
    do i= 1, ndades
        variancia = variancia + (numeros(i)-valormitja)**2
    end do
    variancia=variancia/dble(ndades)
    write(1,*)'variancia =', variancia

    desviacio = dsqrt(variancia)
    write(1,*)'desviacio =', desviacio

    close(1)




end program practica4

!---------------------------------------------------------------------------------
! SUBRUTINES 
! -------------------------------------------------------------------------------

! p(X) PER EXCERCI 1
subroutine funcio(x,px)
    implicit none
    double precision :: px, l=3.d0, x
    px =(5.d0/(324.d0*l))*((x/l)**2)*(9.d0-(x/l)**2)
	return
end 

! G(X) PER EXRCICI 2
subroutine funcio2(x,gx)
    implicit none 
    double precision :: x, gx, sig=4, pi=dacos(-1.d0)
    gx= (dexp((-x**2)/(2.d0*sig**2)))/ (dsqrt(2.d0*pi*(sig**2)))
    return
end

! SUBRUITINA ACCEPTAREBUIG
subroutine accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)
    implicit none
    double precision :: numeros(ndades)
    integer :: i, ndades, iseed
    double precision :: xlow, xhigh
    double precision :: x1, x2, x, p, fx, cotasup
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
			call funcio(x,fx)
		enddo
		numeros(i)=x
	enddo

    valormitja = 0
    do i= 1, ndades
        valormitja = valormitja + numeros(i)
    end do
    valormitja = valormitja/dble(ndades)

    variancia = 0
    do i= 1, ndades
        variancia = variancia + (numeros(i)-valormitja)**2
    end do
    variancia=variancia/dble(ndades)

    desviacio = dsqrt(variancia)



end 

! HISTOGRAMMA EXTRET DEL CAMPUS
SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS,&
    BOXSIZE)
    implicit none

! INPUT/OUTPUT VARIABLES

      INTEGER NDAT,NBOX
      DOUBLE PRECISION XDATA(NDAT),XA,XB
      DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
      INTEGER IERR
      INTEGER I,IBOX,ICOUNT
      DOUBLE PRECISION BOXSIZE

      IF (XA.GE.XB) THEN 
         IERR=1
         RETURN
      ENDIF
! BOX SIZE
        BOXSIZE=(XB-XA)/NBOX

! COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
      ICOUNT=0

! SETS ALL TO ZERO
      DO I=1,NBOX
         VHIS(I)=0
         ERRHIS(I)=0
      ENDDO

! WE RUN THROUGH THE DATASET
      DO I=1,NDAT
! CHECKS IF DATA LIES WITHIN XA,XB
        IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
           IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
! PUTS XB INTO THE LAST BOX, IF NEEDED
           IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

           VHIS(IBOX)=VHIS(IBOX)+1
           ICOUNT=ICOUNT+1
        ENDIF
       ENDDO

       IF (ICOUNT.EQ.0) THEN 
          IERR=2
          RETURN
       ENDIF

       IERR=0
       PRINT*,"ACCEPTED:",ICOUNT," OUT OF:",NDAT

       DO I=1,NBOX
! CENTRAL VALUE OF THE BAR
          XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
!  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
          ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT)) /BOXSIZE / SQRT(DBLE(ICOUNT))
! NORMALIZED VALUE OF THE BAR
          VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
       ENDDO
END


! REGLA DE SIMPSON 3/8 DE LA PRACTICA 4
subroutine simpsontresvuit(x1, x2, k, f, resultat)
	implicit none
	double precision :: x1, x2  ! Límites de integración
	integer :: k        ! Potencia de 3 para el número de subintervalos
	double precision :: resultat  ! Resultado de la integral
	double precision :: h, sum
	integer :: i, n
	double precision :: fx,  fx2, fx1,x ! f(x) utilizados para la subrutina
	external :: f ! Función externa a integrar
  
	n = 3**k  ! Calcula el número de subintervalos (3 elevado a la k)
  
	h = (x2 - x1) / dble(n)  ! Calcula la longitud de cada subintervalo
  
    call f(x1,fx1)
    call f(x2,fx2)
	sum = fx1 + fx2   ! Suma los valores en los extremos del intervalo
  
	! Calculo de los puntos intermedios del intervalo
	do i = 1, n - 1
	  if (mod(i, 3) == 0) then ! Para multiplos de 3
		x= x1 + dble(i) * h
        call f(x,fx)
		sum = sum + 2.0d0 * fx
	  else
		x= x1 + dble(i) * h
        call f(x,fx)
		sum = sum + 3.0d0 * fx
	  end if
	end do
  
    resultat = 3.0d0 * h * sum / 8.0d0 ! Calcula el resultado final
	return
end 