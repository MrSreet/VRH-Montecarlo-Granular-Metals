      module prceldas
      integer(kind=4) npar,npar1,npar2,npar3,npar4,npar5,npar6
      integer(kind=4) npold1,npold2,npold3,npold4,npold5,npold6
      integer(kind=4), allocatable, dimension(:) :: isit1,isit2
      real(kind=8), allocatable, dimension(:):: pr1,pr2,pr3,pr4,pr5,pr6
      real(kind=8) :: pr1ct,pr2ct,pr3ct,pr4ct,pr5ct,deltat
      end module
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!           Relajacion por el metodo de MonteCarlo                     C
!                                                                      C
!   Authors:   Cayetano Hernández Sánchez                              C
!                                                                      C
!   Version:   08.03.13                                                C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      PROGRAM      MC_EG   ! Monte Carlo Electron glass
      use prceldas
!$    use omp_lib
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER      NPMAX,ISEED1,ISEED2,nprom
      PARAMETER   (NPMAX=100,nprom=100000)
!
      real*4       Ts,TT(2),T1,T2,timefin
      real*8       T,xmonte, CONDUC, NTERM, CON, CON2
      REAL*8       RMU,F,PHI,FK,FMIN,B,X,Y,Z,AMIN,rmax, inten
      real*8       xmu,work,rloc,ttime,e
      INTEGER*4    NP,ID,NREALIZ,I,J,IREALIZ,NPC,ireflet
      integer*4    NC1,iu,ind
      integer(kind=8) :: nexit
      LOGICAL      SWCOU
      CHARACTER*3  SW
      character*1  letra
      character*20 filestore,fileinput
      DIMENSION    F(NPMAX,NPMAX),et(nprom),etot(nprom),atime(nprom)
      DIMENSION    PHI(NPMAX),X(NPMAX),Y(NPMAX),Z(NPMAX),AA(NPMAX)
      DIMENSION    NC1(NPMAX),E(NPMAX),ind(nprom)
!$    integer*4 proc_num, thread_num
!
      NP         = npmax          ! number of points ejmplo al inicio 100/ luego 400
      ID        = 2             ! dimension
      B         = 0.D0          ! phi=[-B/2,B/2] rango del desorden B = 0 luego 1
      SW        = 'rnd'         ! rnd or ltc lattice
      FK        = 0.D0         ! ocupation inicial
      RMU       = 0.0D0         ! chemical potential
      SWCOU     = .TRUE.
      rloc      = 1.d0          ! localization length
      prcorte   = 1.d-12        ! corte en prob para las rf
      rmax      = 25.d0         ! largest pair distance considered
      rmax=dmin1(rmax,-0.5d0*rloc*dlog(prcorte))
      nstep     = 1200 ! numero de pasos montecarlo entre promedios es variable y empieza en 10, 
      NTERM     = (NSTEP / 5.0D0) * NP
      npas      = 100  ! numero de promedios (pasos totales=nstep*npas)
      Nocup     = 1             ! num de ocup iniciales por muestra
      NREALIZ   = 100            ! number of realizations 8ANTES A 1)
      TEMPIN      = 8.D0     !, 2.D0, 2.5D0, 3.D0, 3.5D0       ! TEMPERATURA INERSA
      TEMP      = 1.0 / TEMPIN ! TEMPERATURA I.E INVERSO DEL INVERSO
!      TEMP      = 1.D0
      Efield    = TEMP/7.d0    ! Campo electrico
      prcorte   = 0.99d0
      inputfile = 0         ! 0 no input, 1 fichero de imput
      iu        = 80      ! unidad donde se guardan las configuraciones
      iinp      = 81      ! unidad donde se leen las configuraciones
      filestore = 'conf.out'
      fileinput = 'conf.inp'
      ISEED1    = 1210101
      ISEED2    = 1203333
      nblock    = 1
      CON = 0.0
      CON2 = 0.0
!
      if(npas.gt.nprom) then
            write(6,*) 'aumenta nprom',npas,nprom
            stop
      endif
      RL=DBLE(NP)**(1.D0/DBLE(ID))
      RL2=RL/2.d0
!      write(*,*) rl2
!
      IF(NP.NE.(NP/2)*2) WRITE(6,*) 'Cuidado NP no es par'
!$    proc_num = omp_get_num_procs ( )
!$    thread_num = omp_get_max_threads ( )
!$    write ( *, '(a)' ) ' '
!$    write ( *,'(a,i8)') ' Number of processor available=', proc_num
!$    write ( *,'(a,i8)') ' Number of threads available=', thread_num
!------------------------------------------------------------
! inicializa promedios sobre realizaciones
!
        do i=1,10000
        etot(i)=0.d0
        et(i)=0.d0
        atime(I)=0.d0
        enddo
!
      open(iu,file=filestore,status='unknown')
      write(iu,*) nrealiz,nocup
      
      
      
! comienza bucle sobre realizaciones
      xsalfin0 = 0
      xsalfin = 0 ! se pone a cero xsalfin aqui
      DO IREALIZ=1,NREALIZ !Como NREALIZ = 1, ESTE BUCLE NO HACE NADA 
       if(int(dble(irealiz)/10.d0)*10.eq.irealiz)write(*,*)' tirada= ',
     &        irealiz
!Importa muestra si inputfile=1
!  INICIALIZA RANDOM
       CALL srandom(ISEED1,ISEED2)
!
!------------------------------------------
! CONSTRUYE LA MUESTRA E INICIALIZA  TABLAS
       CALL NEWBUILD(NPMAX,NP,PHI,ID,SW,AMIN,X,Y,Z,B,F) ! no toques nada 
       call celdas(np,npmax,F,rmax,rloc,prcorte) !nada
       write(6,*) 'deltat=',deltat
       cintens=TEMP/(Efield*(RL**id)*deltat*dble(nstep*np))
!       write(6,*) temp/efield,deltat*dble(np)
!       write(6,*) cintens
!--------------------------
! Ocupaciones iniciales ! de momento mitad uno mitad cero, poner tods a 0. Hacer bucle. Quitar este bucle, nc1 = 0 \\\ definir nuava variable ecarga que sea igual 1 + phi. Subrutina Sum(E.Q^2_i, q es nc, poner aqui no en energy
        do iocup=1,nocup !El bucle se alarga mucho hacia abajo, pero no pasa nada, podriamos poner endo aqui abajo y no cambiaria nada 
         DO I=1,NP
           E(I) = PHI(i) + 1.d0
           NC1=0
         ENDDO

         
!NP es array de los factores de llenado, pero no entiendo
!   Empieza simulacion montecarlo.
!
         timefin=0.d0
         atime=0.d0
         nexit=0
        do ipas=1,npas
       if((ipas/1000)*1000.eq.ipas)write(*,*)' pasos= ', ipas
! inicializacion de promedios
        enermed=0.d0
        enermed2=0.d0
        xsaltot=0.d0
        ener = 0.d0
        XSAL0 = 0.0D0
!------------------
! paso montecarlo
        do imonte=1,nstep*np,nblock ! DE 1 A NSTEP*NP EN PASOS DE NBLOCK
        
         
         
         ! empieza el factor espacial
          do j=1,nblock
            call spacial(ipar)
            ind(j)=ipar
          enddo ! j
! acaba el factor espacial
!  --------------
! empieza el factor de energias
        do j=1,nblock
        ipar=ind(j)
!        i1=isit1(ipar) !CAMBIAR POR LO ESCRITO A CONTINUACION
!        i2=isit2(ipar)
         j1 = isit1(ipar)
         j2 = isit2(ipar)
         i1 = j1
         i2 = j2
         call random_number(xr)
         if (xr .gt. 0.5) then !que tengo que poner en lugar de xr? por qué invertir?         
            i1 = j2
            i2 = j1            
         end if
         
!        if (abs(nc1(i1) - nc1(i2)) .le. 1) then      !  Esto se quita del todo
           xn1=dble(1-2*nc1(i1))
           xn2=dble(1-2*nc1(i2))
!           enerdif=xn1*e(i1)+xn2*e(i2)  !Cambiar por energia de carga. (Abajo) !! formula de la pizarra
           dif2 = (2*nc1(i2) + 1)*e(i2)
           dif1 = (-2*nc1(i1) + 1)*e(i1)

           enerdif = dif1 + dif2
*           write(*,*) enerdif
           xsal=x(i2)-x(i1)
!           if(dabs(xsal) > RL2) xsal=RL-xsal incorrecto
           if(xsal >  RL2) xsal=xsal-RL
           if(xsal < -RL2) xsal=xsal+RL
!           xsal=xn1*xsal
           enerdif=enerdif-xsal*Efield
*           write(*,*) enerdif, xsal*Efield
           jump=0
           if(enerdif.lt.0.d0) then
             jump=1
           else
              call random_number(xmonte)
              xener=dexp(-enerdif/TEMP)
              if(xmonte.lt.1.d-4) then
                xener=xener*1.d4
              call random_number(xmonte)
              endif
              if(xmonte.lt.1.d-4) then
                xener=xener*1.d4
                call random_number(xmonte)
              endif
               if(xmonte.lt.1.d-4) then
                xener=xener*1.d4
                call random_number(xmonte)
              endif
              if(xmonte.lt.xener) then
               jump=1
              endif
           endif
              
           
              if(jump.eq.1)then
                 nc1(i1)=nc1(i1)-1 !ni esto?
                 nc1(i2)=1+nc1(i2)
                 ener=ener+enerdif
                 nexit=nexit+1
                 xsaltot=xsaltot+xsal
*                 write(*,*) nc1(i1), nc1(i2), ener, xsaltot
               endif
! acaba  el factor de energias
            END DO !end imonte
         
          
         if(ipas.eq.npas*0.2) then
         
            xsaltot = XSAL0
            
         end if
        
        enddo !end ipas

************************************************************************
        DO I=1,NP
           E(I) = PHI(i) + 1.d0
        ENDDO
! promedios del montecarlo
       enermed=enermed+ener
       enermed2=enermed2+ener**2
       ratio=dble(nexit)/dble(nstep*np)
        etot(ipas)=etot(ipas)+ener
        atime(ipas)=atime(ipas)+xsaltot*cintens
        et(ipas)=et(ipas)+ener**2
c       et(ipas)=et(ipas)+ratio
!        write(6,*) xsaltot
        xsalfin=xsalfin+xsaltot
        
        xsalfin0 = xsalfin0 + XSAL0
        
*        write(238,*) deltat*dble(nstep*np)*dble(ipas),xsalfin
*          write(*,*) deltat*dble(nstep*np)*dble(ipas),xsalfin, xsaltot
        enddo ! end ipas
!TERMINA MONTECARLO
! PROMEDIO SOBRE REALIZACIONES Y OCUPACIONES
!      guarda ultima configuracion
       call store(NP,id,b,amin,fk,rmu,rmax,rloc,temp,iseed1,iseed2,
     &       ttime,nc1,iu)
!
!---------------------------------------
! CALCULAMOS LOGARITMO DE CONDUCTIVUDAD PARA REALIZACIÓN PARA LUEGO CALCULAR EL ERROR


        VOLT = Efield * RL !POR QUE SOLO NOS INTERESA LA DIFERENCIA NO?
        INTEN = ((xsalfin - xsalfin0) / RL) /( NPAS* NSTEP* 0.8)

        CONDUC = INTEN / VOLT !fundamental, luego es la inversa de ohm
        
        CON = CON + LOG(CONDUC)
        CON2 = CON2 + LOG(CONDUC)**2







!---------------------------------------
!  ESCRIBE PROMEDIOS

        nsample=(irealiz-1)*nocup+iocup
        xsample=dble(nsample)
        open(16,file='prom.dat',status='unknown')
       atim=0.d0
       do i=1,npas
       ee=etot(i)/xsample
       ee2=et(i)/xsample
       de=dsqrt(ee2-ee*ee)/dble(np)
       atim=atim+atime(i)/xsample
       write(16,401) i*nstep,ee/dble(np),de,atime(I),nsample
 401    format(1x,i9,1x,e15.8,2(1x,e12.6),I6)
       enddo
       close(16)
       write(6,*) 'nexito=',dble(nexit)/(dble(nstep*np)*dble(npas))
       write(6,*) 'G*T=',atim/dble(npas)
       enddo  ! end iocup
c------------------------------------------------------------!
! PREPARA SEED PARA SIGUIENTE REALIZACION
       ISEED1=ISEED1+5
       IF(ISEED1.GT.31328) ISEED1=ISEED1-31328
       ISEED2=ISEED2+500
       IF(ISEED2.GT.30081) ISEED2=ISEED2-30081
c------------------------------------------------------------!
       ENDDO ! end irealiz
       

       close(iu)
       close(iinp)
       
       
* CALCULAMOS VARIANZA 

        VAR = ABS(CON2 / NREALIZ - (CON**2)/ NREALIZ**2)
        ER = SQRT(VAR/NREALIZ)
* TERMINADO BUCLE DE REALIZACIONES TOMAMOS PROMEDIOS DE ENERGÍA PARA CALCULAR CONDUCTIVIDAD EN FUNCION DE 1/T

       VOLT = Efield * RL !POR QUE SOLO NOS INTERESA LA DIFERENCIA NO?
       INTEN = ((xsalfin - xsalfin0) / RL) / (NREALIZ* NPAS* NSTEP* 0.8)

       
       CONDUC = INTEN / VOLT !voltaje /intensidad?
       
       
       
!       OPEN(10,file='MCCARGA-B1-10.DAT',access='append')



         OPEN(10,file='MCCARGAB-B0-10.DAT',access='append')
!           OPEN(10,file='MCCARGA-B0-40.DAT',access='append')
!          OPEN(10,file='MCCARGA-B0-30.DAT',access='append')
!         OPEN(10,file='MCCARGA-B0-20.DAT',access='append')
        WRITE(10,*) TEMPIN, LOG(CONDUC), xsalfin / RL, ER, NSTEP
        CLOSE (10)
      END

      include 'mclib.f'

