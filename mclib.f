!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!  Task:    Initialisation of the occupation numbers                   C
!                                                                      C
!  Input:   NP     - actual number of sites                            C
!           FK     - initial mean occupation number                    C
!                                                                      C
!  Output:  N - array of occupation numbers                            C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE NINIT(NP,FK,N)
      IMPLICIT  NONE
      INTEGER   N,NP,I,IC,NEL,NUM,NES
      REAL*8    FK,Rvec
      DIMENSION N(NP)
!
      NEL = FK*NP
      IF ((FK.LE.0.5D0).AND.(FK.GE.0.0D0)) THEN !ERA GT EN LUGAR DE GE, PERO HE PUESTO GE PARA PODER PONER FK = 0 COMO QUEREMOS
        DO I=1,NP
          N(I) = 0
        ENDDO
        IC  = 0
        DO WHILE (IC.NE.NEL)
          call random_number(rvec)
          NUM = NP*Rvec
          NUM = MOD(NUM,NP) + 1
          IF (N(NUM).EQ.0) THEN
            N(NUM) = 1
            IC     = IC + 1
          ENDIF
        ENDDO
      ELSE IF ((FK.GT.0.5D0).AND.(FK.LT.1.0D0)) THEN
        NES = (1 - FK)*NP
        DO I=1,NP
          N(I) = 1
        ENDDO
        IC  = 0
        DO WHILE (IC.NE.NES)
          call random_number(rvec)
          NUM = NP*Rvec
          NUM = MOD(NUM,NP) + 1
          IF (N(NUM).EQ.1) THEN
            N(NUM) = 0
            IC     = IC + 1
          ENDIF
        ENDDO
      ELSE
        PRINT *,'*** ERROR: UNREASONABLE FK = ', FK
        STOP
      ENDIF
!
C Dummy commands to neutralize an HP compiler bug:
      if(ic.lt.0) then
        print *,'*** error in ninit: ic =',ic
        stop
      endif
!
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Task:    Generation of  sites                                       C
C                                                                      C
C  Input:   NP    - actual number of sites                             C
C           ID    - dimension of the sample                            C
C           SW    -    'rnd'/'ltc'                                     C
C           AMIN  -  minimun value of distance to 'rnd' positions      C
C                                                                      C
C  Output:  X,Y,Z - arrays of randomly distributed in range            C
C                  (0..L) coordinates L=NP**(1/D) or a lattice         C
C                                                                      C
C                                                                      C
C  Authors: A. Neklioudov, A. Mobius, A. Diaz-Sanchez                  C
C                                                                      C
C  Version: 16.07.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE SITEGEN(NP,ID,SW,AMIN,X,Y,Z)
!
      IMPLICIT    NONE
      INTEGER     ID,NP,K,I,NPS,J,LI
      REAL*8      X,Y,Z,RL,RL2,RIJ,DIST,AMIN,AMINRU,Rvec
      CHARACTER*3 SW
      DIMENSION   X(NP),Y(NP),Z(NP)
!
      K = 0
      IF(SW .EQ. 'rnd') THEN     !---------------------
        RL=DBLE(NP)**(1.D0/DBLE(ID))
        RL2=RL/2.D0
        AMINRU=AMIN**2.D0
        DO I=1,NP
550       CONTINUE
          call random_number(rvec)
          X(I)=RL*Rvec
          IF(ID.GE.2) then
             call random_number(rvec)
             Y(I)=RL*Rvec
          ENDIF
          IF(ID.GE.3) then
             call random_number(rvec)
             Z(I)=RL*Rvec
          ENDIF
          DO J=1,I-1
            DIST=0.D0
            RIJ=DABS(X(J)-X(I))
            IF(RIJ.GT.RL2) RIJ=RL-RIJ
            DIST=DIST+RIJ**2.D0
            IF(ID.GE.2) THEN
              RIJ=DABS(Y(J)-Y(I))
              IF(RIJ.GT.RL2) RIJ=RL-RIJ
              DIST=DIST+RIJ**2.D0
              IF(ID.GE.3) THEN
                RIJ=DABS(Z(J)-Z(I))
                IF(RIJ.GT.RL2) RIJ=RL-RIJ
                DIST=DIST+RIJ**2.D0
              ENDIF
            ENDIF
            IF(DIST.LT.AMINRU) GOTO 550
          ENDDO
        ENDDO
      ELSE  !-------------------------------
        AMIN=1.D0
        NPS = NP**(1.0D0/DBLE(ID)) + 0.1D0
        IF (IABS(NPS**ID - NP).GT.0) THEN
          PRINT *,'*** ERROR: UNREASONABLE NP = ', NP
          STOP
        ENDIF
        IF(ID.EQ.1) THEN
          DO I=1,NP
            K    = K + 1
            X(K) = I
          ENDDO
        ELSE IF(ID.EQ.2) THEN
          DO I=1,NPS
            DO J=1,NPS
              K    = K + 1
              X(K) = I
              Y(K) = J
            ENDDO
          ENDDO
        ELSE IF(ID.EQ.3) THEN
          DO I=1,NPS
            DO J=1,NPS
               DO LI=1,NPS
                 K    = K + 1
                 X(K) = I
                 Y(K) = J
                 Z(K) = LI
               ENDDO
            ENDDO
          ENDDO
        ELSE
          PRINT *,'*** ERROR: UNREASONABLE ID = ',ID
          STOP
        ENDIF
      ENDIF
!
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Task:    Calculation of the interaction matrix according to the     C
C           periodic boundary conditions (see e.g. formula (10) MRD)   C
C                                                                      C
C  Input:   NPMAX  - maximum number of sites for array's initalisation C
C           NP     - actual number of sites                            C
C           ID     - dimension we're working with                      C
C           X,Y,Z  - coordinate arrays                                 C
C                                                                      C
C  Output:  F - array = 1/Rij where F(i=j)=0                           C
C                                                                      C
C                                                                      C
C  Authors: A. Neklioudov, A. Mobius                                   C
C                                                                      C
C  Version: 26.05.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DIST(NPMAX,NP,ID,X,Y,Z,F)
      IMPLICIT  NONE
      INTEGER(kind=4), intent(in) :: ID,NPMAX,NP
      REAL(kind=8),intent(in),dimension(npmax) :: X,Y,Z
      REAL(kind=8),intent(out),dimension(npmax,npmax) ::  F
      REAL(kind=8) :: RLX,RLY,RLZ,RL
      INTEGER(kind=4) :: i,j
!
      RL=dble(NP)**(1.D0/DBLE(ID))
      IF (ID.EQ.1) THEN
        DO I=1,NP
          F(I,I)=0.0D0
          DO J=I+1,NP
            RLX = DMIN1(DABS(X(I)-X(J)),RL-DABS(X(I)-X(J)))
            F(I,J)=1.D0/RLX
            F(J,I)=F(I,J)
          ENDDO
        ENDDO
      ELSE IF (ID.EQ.2) THEN
        DO I=1,NP
          F(I,I)=0.0D0
          DO J=I+1,NP
            RLX= DMIN1(DABS(X(I)-X(J)),RL-DABS(X(I)-X(J)))
            RLY= DMIN1(DABS(Y(I)-Y(J)),RL-DABS(Y(I)-Y(J)))
            F(I,J)=1.D0/DSQRT(RLX*RLX+RLY*RLY)
            F(J,I)=F(I,J)
          ENDDO
        ENDDO
      ELSE IF (ID.EQ.3) THEN
        DO I=1,NP
          F(I,I)=0.0D0
          DO J=I+1,NP
            RLX= DMIN1(DABS(X(I)-X(J)),RL-DABS(X(I)-X(J)))
            RLY= DMIN1(DABS(Y(I)-Y(J)),RL-DABS(Y(I)-Y(J)))
            RLZ= DMIN1(DABS(Z(I)-Z(J)),RL-DABS(Z(I)-Z(J)))
            F(I,J)=1.D0/DSQRT(RLX*RLX+RLY*RLY+RLZ*RLZ)
            F(J,I)=F(I,J)
          ENDDO
        ENDDO
      ELSE
        PRINT *,'*** ERROR: UNREASONABLE ID = ',ID
        STOP
      ENDIF
!
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Task:    Generation of random potentials on each site               C
C                                                                      C
C  Input:   NP - actual number of sites                                C
C           B  - phi = [-B/2;+B/2]                                     C
C                                                                      C
C  Output:  PHI - array of random potential on each site               C
C                 (rectangular distribution)                           C
C                                                                      C
C                                                                      C
C  Authors: A. Neklioudov, A. Mobius                                   C
C                                                                      C
C  Version: 26.05.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE POT(NP,B,PHI)
!
      IMPLICIT  NONE
      INTEGER   NP,I
      REAL*8    PHI,B,Rvec
      DIMENSION PHI(NP)
!
      DO I=1,NP
        call random_number(rvec)
        PHI(I)=B*(Rvec - 0.5D0)
      ENDDO
!
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Task:    Determination of the minimum value of interaction          C
C                                                                      C
C  Input:   NPMAX  - maximum number of sites                           C
C           NP     - actual number of sites                            C
C           F      - array of inverse distances                        C
C           ISRT   - matrix ISRT(i,j)                                  C
C                                | |___ number of site                 C
C                                |                                     C
C                                |___ number of site acc. to descent   C
C                                     inverse distance                 C
C                                                                      C
C  Output:  FMIN - minimum value of interaction                        C
C                                                                      C
C                                                                      C
C  Authors: A. Neklioudov, A. Mobius                                   C
C                                                                      C
C  Version: 31.05.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE MININT(NPMAX,NP,F,ISRT,FMIN)
!
      IMPLICIT  NONE
      INTEGER   NPMAX,NP,ISRT,I,J
      REAL*8    F,FMIN
      DIMENSION ISRT(NPMAX,NP),F(NPMAX,NP)
!
      J    = NP - 1
      FMIN = F(1,ISRT(J,1))
      DO I=1,NP
        FMIN = DMIN1(FMIN,F(I,ISRT(J,I)))
      ENDDO
!
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Task:    Calculation of total energy of the system formula (8)&(11) C
C           MRD and single-particle energies formula (13) MRD          C
C                                                                      C
C                                                                      C
C  Input:   NPMAX  - maximum number of sites                           C
C           NP     - actual number of sites                            C
C           F      - array of inverse distances between sites          C
C           N      - array of the sample filling factors               C
C           PHI    - array of potentials on each site                  C
C           FK     - mean filling factor of the sites                  C
C           RMU    - chemical potential of the system                  C
C           AA     - auxiliary array                                   C
C                                                                      C
C  Output:  HCA    - total energy in canonical ansemble                C
C           HGCA   - total energy in grand canonical ansemble          C
C           E      - array of single particle energies                 C
C                                                                      C
C                                                                      C
C  Authors: A. Neklioudov, A. Mobius                                   C
C                                                                      C
C  Version: 26.05.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE ENERGY(NPMAX,NP,F,N,PHI,FK,RMU,AA,HCA,HGCA,E)
      IMPLICIT  NONE
      INTEGER   N,NPMAX,NP,I,J
      REAL*8    F,PHI,AA,E,FK,HCA,HGCA,H2,PHII,RMU,RNI,H1
      DIMENSION N(NP),F(NPMAX,NP),PHI(NP),AA(NP),E(NP)
!
      HCA  = 0.0D0
      HGCA = 0.0D0
      DO I=1,NP
        AA(I) = N(I) - FK
      ENDDO
!
      DO I=1,NP
        PHII = PHI(I)
        RNI  = DBLE(N(I))
        H2   = 0.0D0
        DO J=1,NP
          H2 = H2 + F(J,I)*AA(J)
        ENDDO
        H1   = PHII*RNI + AA(I)*H2*0.5D0
        HCA  = HCA + H1
        HGCA = HGCA + H1 - RMU*RNI
        E(I) = PHII + H2
      ENDDO
!
      RETURN
      END

!=======================================================================
!   SUBROUTINE SRANDOM(is1,is2)
!------------------------------------------------------------------
!  subroutina para inicializar del generador de numeros aleatorios
!  intrinseco de fortran 90 utilizando solo 2 semillas, is1, is2
!-----------------------------------------------------------------------
      subroutine srandom(is1,is2)
      implicit none
      integer(kind=4), intent(inout) :: is1,is2
      integer(kind=4), allocatable :: iseed(:)
      integer(kind=4) :: n,i,itop
!
      itop=2**30
      if(is1 > itop) is1=is1-itop
      if(is2 > itop) is2=is2-itop
      call random_seed(size=n)
!      write(6,*) 'numero de seeds=',n
      allocate(iseed(n))
      do i=1,n,2
         iseed(i)=is1
      enddo
      if(n > 1) then
         do i=2,n,2
            iseed(i)=is2
         enddo
      endif
      call random_seed(put=iseed(1:n))
      deallocate(iseed)
!
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   Translation between extended and compressed representations of     C
C   the occupation vector                                              C
C                                                                      C
C                                                                      C
C   N          -  state (array of occupation numbers)   (input/output) C
C   NP         -  number of sites                       (input)        C
C   NCOMPR     -  occupation vector in compressed       (input/output) C
C                 representation                                       C
C   NPC        -  size of occupation vector in          (input)        C
C                 compressed representation                            C
C                 (NP/31 if MOD(NP,31)=0,                              C
C                 otherwise NP/31+1)                                   C
C   ISW        -  switch:   -1 = compression,           (input)        C
C                            1 = explosion                             C
C                                                                      C
C                                                                      C
C   Authors:   A. Mobius, A. Diaz-Sanchez                              C
C                                                                      C
C   Version:   13.05.97                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TRANSL(N,NP,NCOMPR,NPC,ISW)
      IMPLICIT  NONE
      INTEGER   NP,NPC,N,NCOMPR,ISW,J,I0,IM,NN,I
      DIMENSION N(NP),NCOMPR(NPC)
!
      IF(ISW.LT.0) THEN
         DO J=1,NPC
            I0=(J-1)*31
            IM=31
            IF(J.EQ.NPC) IM=NP-I0
            NN=0
            DO I=1,IM
               NN=2*NN+N(I0+I)
            ENDDO
            NCOMPR(J)=NN
         ENDDO
      ENDIF
!
      IF(ISW.GT.0) THEN
         DO J=1,NPC
            I0=(J-1)*31
            IM=31
            IF(J.EQ.NPC) IM=NP-I0
            NN=NCOMPR(J)
            DO I=IM,1,-1
               N(I0+I)=MOD(NN,2)
               NN=(NN-N(I0+I))/2
            ENDDO
         ENDDO
      ENDIF
!
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              Quick sorting procedure                                 C
C  from "Numerical recipes in FORTRAN", second ed., p. 330-3,          C
C  slightly modified.                                                  C
C                                                                      C
C  Input:  N    - number of elements in the array                      C
C          ARR  - array to be sorted                                   C
C          Both these quantities remain unchanged.                     C
C                                                                      C
C  Output: INDX - array with numbers of elements according to          C
C                 ascending order in array ARR                         C
C                                                                      C
C  Version: 02.06.97                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      SUBROUTINE INDEXX(N,ARR,INDX)
!
      IMPLICIT NONE
      INTEGER N,INDX(N),M,NSTACK
      REAL*8  ARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A
!
      DO J=1,NSTACK
        ISTACK(J)=0
      ENDDO
!
      DO J=1,N
        INDX(J)=J
      ENDDO
      JSTACK=0
      L=1
      IR=N
    1 IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          INDXT=INDX(J)
          A=ARR(INDXT)
          DO I=J-1,1,-1
            IF(ARR(INDX(I)).LE.A) GOTO 2
            INDX(I+1)=INDX(I)
          ENDDO
          I=0
    2     INDX(I+1)=INDXT
        ENDDO
        IF (JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        ITEMP=INDX(K)
        INDX(K)=INDX(L+1)
        INDX(L+1)=ITEMP
        IF(ARR(INDX(L+1)).GT.ARR(INDX(IR))) THEN
          ITEMP=INDX(L+1)
          INDX(L+1)=INDX(IR)
          INDX(IR)=ITEMP
        ENDIF
        IF(ARR(INDX(L)).GT.ARR(INDX(IR))) THEN
          ITEMP=INDX(L)
          INDX(L)=INDX(IR)
          INDX(IR)=ITEMP
        ENDIF
        IF(ARR(INDX(L+1)).GT.ARR(INDX(L))) THEN
          ITEMP=INDX(L+1)
          INDX(L+1)=INDX(L)
          INDX(L)=ITEMP
        ENDIF
        I=L+1
        J=IR
        INDXT=INDX(L)
        A=ARR(INDXT)
    3   CONTINUE
        I=I+1
        IF(ARR(INDX(I)).LT.A) GOTO 3
    4   CONTINUE
        J=J-1
        IF(ARR(INDX(J)).GT.A) GOTO 4
        IF(J.LT.I) GOTO 5
        ITEMP=INDX(I)
        INDX(I)=INDX(J)
        INDX(J)=ITEMP
        GOTO 3
    5   INDX(L)=INDX(J)
        INDX(J)=INDXT
        JSTACK=JSTACK+2
        IF(JSTACK.GT.NSTACK) THEN
          PRINT *,'*** NSTACK TOO SMALL IN INDEXX ***'
          STOP
        ENDIF
        IF(IR-I+1.GE.J-L) THEN
          ISTACK(JSTACK)=IR
          ISTACK(JSTACK-1)=I
          IR=J-1
        ELSE
          ISTACK(JSTACK)=J-1
          ISTACK(JSTACK-1)=L
          L=I
        ENDIF
      ENDIF
      GOTO 1
!
      END

!=======================================================================
!     subroutine celdas
!---------------------------------------------------------------------
!  subrutina que construye las celdas del factor para el MonteCarlo
!  tipo Efros. La informacion final se guarga en un modulo.
!  deltat(un sorteo)=1/zz
      subroutine celdas(np,npmax,F,rmax,rloc,prct,zz)
      use prceldas
      implicit real*8 (a-h,o-z)
      integer(kind=4) :: np,npmax
      real(kind=8), intent(in) :: rloc,prct
      real(kind=8), intent(in), dimension(npmax,npmax) :: F
      real(kind=8), dimension(:), Allocatable :: rf,work
      real(kind=8), intent(out) :: zz
      integer(kind=4),dimension(:),Allocatable:: iwork,ind
      real(kind=8) :: rmax,RLX
      integer(kind=4) :: i,j,ipar,Allocatestatus
!
! calculamos npar
      npar=0
      zz = 0.
      DO i=1,NP
         do j=i+1,np
            RLX = 1.d0/f(i,j)
            if(RLX.lt.rmax) then
               npar=npar+1
               zz = zz + dexp(-2.* RLX/rloc)
            endif
         enddo
      ENDDO
! dimensionamos vectores
      allocate(rf(npar),work(npar),stat=Allocatestatus)
      IF (AllocateStatus /= 0) STOP " Not enough memory *1*"
      if(allocated(isit1))then
         deallocate(isit1,isit2, stat = AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem *2*"
      endif
      allocate(isit1(npar),isit2(npar),stat=Allocatestatus)
      IF (AllocateStatus /= 0) STOP " Not enough memory *2*"
      allocate(iwork(npar),ind(npar),stat=Allocatestatus)
      IF (AllocateStatus /= 0) STOP " Not enough memory *3*"
! repetimos, guardamos y ordenamos
       ipar=0
         DO i=1,NP
           do j=i+1,np
             RLX = 1.d0/F(i,j)
             if(RLX.lt.rmax) then
               ipar=ipar+1
               if(ipar.gt.npar) then
                   write(6,*) 'increase nparmax'
                   stop
               endif
               rf(ipar)=RLX
               isit1(ipar)=i
               isit2(ipar)=j
             endif
           enddo
         ENDDO
!---------------------
! ordena
      call INDEXX(ipar,rf,ind)
      do i=1,ipar
        work(I)=rf(ind(i))
        iwork(i)=isit1(ind(i))
      enddo
      do i=1,ipar
        rf(i)=work(i)
!        write(66,*) i,rf(i)
        isit1(i)=iwork(i)
        iwork(i)=isit2(ind(i))
      enddo
      do i=1,ipar
        isit2(i)=iwork(i)
      enddo
!
      deallocate(iwork,work,ind)
!---------------------------------------------------
! ahora construimos las caja.
! calculamos npar1. Primero calculamos npar1 y despues repetimos.
      if(allocated(pr1))then
         deallocate(pr1,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr1"
      endif
      if(allocated(pr2))then
         deallocate(pr2,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr2"
      endif
      if(allocated(pr3))then
         deallocate(pr3,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr3"
      endif
      if(allocated(pr4))then
         deallocate(pr4,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr4"
      endif
      if(allocated(pr5))then
         deallocate(pr5,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr5"
      endif
      if(allocated(pr6))then
         deallocate(pr6,stat=AllocateStatus)
         IF (AllocateStatus /= 0) STOP " Deallocation problem pr6"
      endif
        ZZ=0.d0
        do i=1,npar
           ZZ=ZZ+Dexp(-2.d0*rf(I)/rloc)
        enddo
        deltat=1.d0/zz
        write(6,*) 'zz,npar=',zz,npar
        npar1=1
        tot=Dexp(-2.d0*rf(npar1)/rloc)/zz
        do while(tot.lt.prct)
        npar1=npar1+1
        fc=Dexp(-2.d0*rf(npar1)/rloc)/zz
        tot=tot+fc
        enddo
        npar1=npar1-1
        if(npar1.ne.0) then
          allocate(pr1(npar1),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *4*"
          npar1=1
          tot=Dexp(-2.d0*rf(npar1)/rloc)/zz
          do while(tot.lt.prct)
            pr1(npar1)=tot
            npar1=npar1+1
            fc=Dexp(-2.d0*rf(npar1)/rloc)/zz
            tot=tot+fc
          enddo
          npar1=npar1-1
          pr1ct=pr1(npar1)
        else
           pr1ct=0.d0
        endif
! calculamos npar2. Primero calculamos npar2 y despues repetimos.
        npold1=npar1
        zz=0.d0
        do i=1,npar-npold1
           ZZ=ZZ+Dexp(-2.d0*rf(i+npold1)/rloc)
        enddo
        npar2=1
        tot=Dexp(-2.d0*rf(npar2+npold1)/rloc)/zz
        do while(tot.lt.prct)
           npar2=npar2+1
           fc=Dexp(-2.d0*rf(npar2+npold1)/rloc)/zz
           tot=tot+fc
        enddo
        npar2=npar2-1
        if(npar2.ne.0) then
          allocate(pr2(npar2),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *5*"
          npar2=1
          tot=Dexp(-2.d0*rf(npar2+npold1)/rloc)/zz
          do while(tot.lt.prct)
             pr2(npar2)=tot
             npar2=npar2+1
             fc=Dexp(-2.d0*rf(npar2+npold1)/rloc)/zz
             tot=tot+fc
          enddo
          npar2=npar2-1
          pr2ct=pr2(npar2)
        else
           pr2ct=0.d0
        endif
! calculamos npar3. Primero calculamos npar3 y despues repetimos.
        npold2=npar1+npar2
        zz=0.d0
        do i=1,npar-npold2
           ZZ=ZZ+Dexp(-2.d0*rf(i+npold2)/rloc)
        enddo
        npar3=1
        tot=Dexp(-2.d0*rf(npar3+npold2)/rloc)/zz
        do while(tot.lt.prct)
           npar3=npar3+1
           fc=Dexp(-2.d0*rf(npar3+npold2)/rloc)/zz
           tot=tot+fc
        enddo
        npar3=npar3-1
        if(npar3.ne.0) then
          allocate(pr3(npar3),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *6*"
          npar3=1
          tot=Dexp(-2.d0*rf(npar3+npold2)/rloc)/zz
          do while(tot.lt.prct)
             pr3(npar3)=tot
             npar3=npar3+1
             fc=Dexp(-2.d0*rf(npar3+npold2)/rloc)/zz
             tot=tot+fc
          enddo
          npar3=npar3-1
          pr3ct=pr3(npar3)
        else
           pr3ct=0.d0
        endif
! calculamos npar4. Primero calculamos npar4 y despues repetimos.
        npold3=npar1+npar2+npar3
        zz=0.d0
        do i=1,npar-npold3
           ZZ=ZZ+Dexp(-2.d0*rf(i+npold3)/rloc)
        enddo
        npar4=1
        tot=Dexp(-2.d0*rf(npar4+npold3)/rloc)/zz
        do while(tot.lt.prct)
           npar4=npar4+1
           fc=Dexp(-2.d0*rf(npar4+npold3)/rloc)/zz
           tot=tot+fc
        enddo
        npar4=npar4-1
!       write(6,*) zz,npar4,1.d0-Dexp(-2.d0*rf(npar))/zz,prct
        if(npar4.ne.0) then
          allocate(pr4(npar4),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *7*"
          npar4=1
          tot=Dexp(-2.d0*rf(npar4+npold3)/rloc)/zz
          do while(tot.lt.prct)
             pr4(npar4)=tot
             npar4=npar4+1
             fc=Dexp(-2.d0*rf(npar4+npold3)/rloc)/zz
             tot=tot+fc
          enddo
          npar4=npar4-1
          pr4ct=pr4(npar4)
        else
           pr4ct=0.d0
        endif
! calculamos npar5. Primero calculamos npar5 y despues repetimos.
        npold4=npar1+npar2+npar3+npar4
        zz=0.d0
        do i=1,npar-npold4
           ZZ=ZZ+Dexp(-2.d0*rf(i+npold4)/rloc)
        enddo
        npar5=1
        tot=Dexp(-2.d0*rf(npar5+npold4)/rloc)/zz
        do while(tot.lt.prct)
           npar5=npar5+1
           fc=Dexp(-2.d0*rf(npar5+npold4)/rloc)/zz
           tot=tot+fc
        enddo
        npar5=npar5-1
!       write(6,*) zz,npar5,1.d0-Dexp(-2.d0*rf(npar))/zz,prct
        if(npar5.ne.0) then
          allocate(pr5(npar5),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *8*"
          npar5=1
          tot=Dexp(-2.d0*rf(npar5+npold4)/rloc)/zz
          do while(tot.lt.prct)
             pr5(npar5)=tot
             npar5=npar5+1
             fc=Dexp(-2.d0*rf(npar5+npold4)/rloc)/zz
             tot=tot+fc
          enddo
          npar5=npar5-1
          pr5ct=pr5(npar5)
        else
           pr5ct=0.d0
        endif
! calculamos npar6. Primero calculamos npar6 y despues repetimos.
        npold5=npar1+npar2+npar3+npar4+npar5
        zz=0.d0
        do i=1,npar-npold5
           ZZ=ZZ+Dexp(-2.d0*rf(i+npold5)/rloc)
        enddo
!       npar6=1
!       tot=Dexp(-2.d0*rf(npar6+npold5)/rloc)/zz
!       do while(tot.lt.prct)
!          npar6=npar6+1
!          fc=Dexp(-2.d0*rf(npar6+npold5)/rloc)/zz
!          tot=tot+fc
!       enddo
!       npar6=npar6-1
        npar6=npar-npold5
        if(npar6.ne.0) then
          allocate(pr6(npar6),stat=Allocatestatus)
          IF (AllocateStatus /= 0) STOP " Not enough memory *9*"
          tot=0.d0
          do npar6=1,npar-npold5
             fc=Dexp(-2.d0*rf(npar6+npold5)/rloc)/zz
             tot=tot+fc
             pr6(npar6)=tot
          enddo
          npar6=npar6-1
        endif
!        write(6,*) npar1,npar2,npar3,npar4,npar5,npar6
      return
      end

!----------------------------------------------------------------------
!     subroutine spacial
!----------------------------------------------------------------------
      subroutine spacial(ipar)
      use prceldas
      implicit none
      integer(kind=4), intent(out)::ipar
      real(kind=8) :: xmonte
!
      call random_number(xmonte)
      if(xmonte.lt.pr1ct) then          !primera caja de probabilidades
        call search(npar1,xmonte,pr1,ipar)
      else
        call random_number(xmonte)
        if(xmonte.lt.pr2ct) then        !segunda caja
          call search(npar2,xmonte,pr2,ipar)
          ipar=ipar+npold1
        else
          call random_number(xmonte)
          if(xmonte.lt.pr3ct) then      !tercera caja
            call search(npar3,xmonte,pr3,ipar)
            ipar=ipar+npold2
          else
            call random_number(xmonte)
            if(xmonte.lt.pr4ct) then    !cuarta caja
              call search(npar4,xmonte,pr4,ipar)
              ipar=ipar+npold3
            else
              call random_number(xmonte)
              if(xmonte.lt.pr5ct) then  !quinta caja
                call search(npar5,xmonte,pr5,ipar)
                ipar=ipar+npold4
              else ! sexta caja
                call random_number(xmonte)
                call search(npar6,xmonte,pr6,ipar)
                ipar=ipar+npold5
              endif ! quinta sexta caja
            endif ! cuarta caja
          endif  ! tercera caja
        endif   ! segunda caja
      endif    ! primera caja
!
      end


!-----------------------------------------------------------------------
!                                                                      C
!                     Build the sample                                 C
!                                                                      C
!                                                                      C
!   NPMAX      -  maximum number of sites               (input)        C
!   NP         -  actual number of sites                (input)        C
!   F          -  array of interactions                 (output)       C
!   PHI        -  array of potentials on each site      (output)       C
!   ID         -  dimension of the sample               (input)        C
!   SW         -    'rnd'/'ltc'                         (input)        C
!                                                                      C
!   X,Y,Z      - arrays of randomly distributed in      (output)       C
!                range(0..L) coordinates L=NP**(1/D)                   C
!                 or a lattice                                         C
!   B          - phi = [-B/2;+B/2]                      (input)        C
!                                                                      C
!                                                                      C
!     Authors:   A.M. Somoza, M. Ortuño                                C
!                                                                      C
!     Version:   20.07.2004                                            C
!                                                                      C
!-----------------------------------------------------------------------
      SUBROUTINE NEWBUILD(NPMAX,NP,PHI,ID,SW,AMIN,X,Y,Z,B,F)
      IMPLICIT NONE
      INTEGER(kind=4), intent(in) :: NPMAX,NP,ID
      REAL(kind=8), intent(in) ::  B,AMIN
      REAL(kind=8), intent(out), dimension(npmax) :: PHI,X,Y,Z
      REAL(kind=8), intent(out), dimension(npmax,npmax) :: f
      CHARACTER*3, intent(in) :: SW
      CALL SITEGEN(NP,ID,SW,AMIN,X,Y,Z)
      CALL DIST(NPMAX,NP,ID,X,Y,Z,F)
      CALL POT(NP,B,PHI)
      RETURN
      END


!-----------------------------------------------------------------------
!                                                                      C
!           Search for the first element of a list (x) greater         C
!             or equal to  value                                       C
!                                                                      C
!                                                                      C
!   N      - Dimension of X                                            C
!   value   -  actual value to search                (input)           C
!   x       -  ordered array with the list           (input)           C
!   ii      -  number of the first element >= to value   (output)      C
!                                                                      C
!     Authors:   A.M. Somoza, M. Ortuño                                C
!                                                                      C
!     Version:   20.07.2004                                            C
!                                                                      C
!-----------------------------------------------------------------------
      subroutine search(n,value,x,ii)
      implicit none
      integer(kind=4), intent(in):: n
      real(kind=8), intent(in),dimension(n):: x
      real(kind=8), intent(in):: value
      integer(kind=4), intent(out):: ii
      integer(kind=4) :: imin,imax,idif
!
      imin=0
      imax=n
      idif=imax-imin
      do while(idif.gt.1)
         if(value.gt.x(imin+idif/2)) then
            imin=imin+idif/2
         else
            imax=imin+idif/2
         endif
         idif=imax-imin
      enddo
      ii=imax
      end


!----------------------------------------------------------------
! escribe una configuracion.
      subroutine store(NP,id,b,amin,fk,rmu,rmax,rloc,temp,
     &            iseed1,iseed2,time,nc1,iu)
      integer  Np,id,iseed1,iseed2,nc1,iu,nfil,i
      real*8   b,amin,fk,rmu,rmax,rloc,temp,time
      dimension nc1(np)
!
      write(iu,10) np,id,b
      write(iu,20) amin,fk,rmu
      write(iu,30) rmax,rloc,temp
      write(iu,40) iseed1,iseed2,time
      write(iu,*)
      nfil=np/16
      if(np.ne.(np/16)*16) nfil=nfil+1
      do i=1,nfil-1
        write(iu,50) nc1((i-1)*16+1:i*16)
      enddo
        write(iu,50) nc1((nfil-1)*16+1:np)
 10   format(1x,I6,1x,I3,1x,f8.4)
 20   format(1x,f8.4,1x,f8.4,1x,f8.4)
 30   format(1x,f8.4,1x,f8.4,1x,f10.6)
 40   format(1x,i6,1x,i6,1x,e15.8)
 50   format(1x,16(i2,1x))
      end


!----------------------------------------------------------------
! lee una configuracion, para continuar el montecarlo a partir de ahi.
      subroutine import(NP,id,b,amin,fk,rmu,rmax,rloc,temp,
     &            iseed1,iseed2,time,nc1,iu)
      integer  Np,id,iseed1,iseed2,nc1,iu,nfil,i,Np1,id1
      real*8   b,amin,fk,rmu,rmax,rloc,temp,time
      real*8   b1,am1,fk1,rmu1,rmax1,rloc1,temp1
      dimension nc1(np)
!
      read(iu,10) np1,id1,b1
      read(iu,20) am1,fk1,rmu1
      read(iu,30) rmax1,rloc1,temp1
      if(np1.ne.np.or.id1.ne.id.or.dabs(b1-b).gt.1.d-6.or.
     &   dabs(amin-am1).gt.1.d-6.or.dabs(fk1-fk).gt.1.d-6.or.
     &   dabs(rmu1-rmu).gt.1.d-6.or.dabs(rmax1-rmax).gt.1.d-6.or.
     &   dabs(rloc1-rloc).gt.1.d-6.or.dabs(temp1-temp).gt.1.d-6) then
         write(6,*) 'la muestra no corresponde a los parametros'
         write(6,*) 'CHECKEA parametros o muestra'
         stop
      endif
      read(iu,40) iseed1,iseed2,time
      read(iu,*)
      nfil=np/16
      if(np.ne.(np/16)*16) nfil=nfil+1
      do i=1,nfil-1
        read(iu,50) nc1((i-1)*16+1:i*16)
      enddo
        read(iu,50) nc1((nfil-1)*16+1:np)
 10   format(1x,I6,1x,I3,1x,f8.4)
 20   format(1x,f8.4,1x,f8.4,1x,f8.4)
 30   format(1x,f8.4,1x,f8.4,1x,f10.6)
 40   format(1x,i6,1x,i6,1x,e15.8)
 50   format(1x,16(i2,1x))
      return
      end

! Guarda una configuracion en formato comprimido en la unidad iu
      subroutine storeconf(ttime,nc1,np,nsal,npc,niter,dt,ed,iu)
      implicit real*8 (a-h,o-z)
      dimension nc1(np),nsal(npc)
      call TRANSL(Nc1,NP,Nsal,NPC,-1)
      write(iu,21) ttime,niter
      write(iu,31) dt,ed
 21   format(1x,'time=',f30.16,' niter=',i8)
 31   format(' deltat=',e21.15,' enerdif=',e21.15)
      do i=1,npc
      write(iu,*) nsal(i)
      enddo
      return
      end

! Lee una configuracion en formato comprimido en la unidad iu
      subroutine readconf(ttime,nc1,np,nsal,npc,niter,dt,ed,iu)
      implicit real*8 (a-h,o-z)
      dimension nc1(np),nsal(npc)
      read(iu,20) ttime,niter
      read(iu,30) dt,ed
 20   format(6x,f30.16,7x,i8)
 30   format(8x,e21.15,9x,e21.5)
      do i=1,npc
      read(iu,*) nsal(i)
      enddo
      call TRANSL(Nc1,NP,Nsal,NPC,1)
      end

! Copia una configuracion comprimida de la unidad iu a la unidad iinp
      subroutine copyconf(ttime,nc1,np,nsal,npc,niter,dt,ed,iu,iinp)
      implicit real*8 (a-h,o-z)
      dimension nc1(np),nsal(npc)
      read(iinp,20) ttime,niter
      read(iinp,30) dt,ed
      write(iu,21) ttime,niter
      write(iu,31) dt,ed
 20   format(6x,f30.16,7x,i8)
 30   format(8x,e21.15,9x,e21.5)
 21   format(1x,'time=',f30.16,' niter=',i8)
 31   format(' deltat=',e21.15,' enerdif=',e21.15)
      do i=1,npc
      read(iinp,*) nsal(i)
      enddo
      do i=1,npc
      write(iu,*) nsal(i)
      enddo
      end


      subroutine setneigh(np,npmax,rmax,f,nveci,neigh,
     &              work,iwork,ind)
      implicit real*8 (a-h,o-z)
      dimension nveci(np),neigh(npmax,np),f(npmax,np),work(np)
      dimension iwork(np),ind(np)
!
      rinv=1.d0/rmax
      do i=1,np
!  recopila todos los vecinos
         ii=0
         do j=1,np
            if(j.ne.i.and.f(i,j).gt.rinv) then
            ii=ii+1
            work(ii)=1.d0/f(i,j)
            iwork(ii)=j
            endif
            nveci(i)=ii
         enddo
! ordena
        call INDEXX(ii,work,ind)
        do j=1,ii
        neigh(i,j)=iwork(ind(j))
        enddo
      enddo
!
      end
