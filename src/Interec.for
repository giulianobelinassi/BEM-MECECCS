************************************************************************
*                                                                      *
*       ESTA SUBROTINA CALCULA OS VALORES DOS DESLOCAMENTOS            *
*                                                                      *
*                        NOS PONTOS INTERNOS                           *
*                                                                      *
************************************************************************
*
      SUBROUTINE INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,ZSSOL,
     $  NE,NX,NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,NBE,RHO,
     $  ETAS,GI,OME,NP)
*
        USE omp_lib

        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
#ifdef USE_GPU
        INCLUDE 'kernels/Interec1_cu.fd'
#endif
        COMMON INP,INQ,IPR,IPS,IPT
        
        REAL, DIMENSION(:), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(L),  INTENT(IN) :: CXI, CYI, CZI
        REAL CO(4,3), DELTA(3,3)
        COMPLEX, DIMENSION(3*NBE), INTENT(INOUT) :: ZDFI, ZFI
        COMPLEX, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: ZDSOL, ZSSOL 
        INTEGER, INTENT(IN) ::  CONE(N,4),KODE(3*NBE)
        DIMENSION ZD(3,3,3),ZS(3,3,3)
        REAL, INTENT(IN) :: ETAS(3,NX)
        
        REAL, INTENT(IN) :: GI(NPG), OME(NPG)
        INTEGER stats1, stats2

        DOUBLE PRECISION :: t1, t2

#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif

#ifdef TEST_CUDA
#undef USE_CPU
#undef USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX, DIMENSION(3*L) :: ZDSOLP
#endif

#ifdef USE_CPU
        COMPLEX, DIMENSION(3,3) :: ZHELEM, ZGELEM 
#endif

*
* REARRANJA OS VETORES ZFI AND ZDFI PARA ARMAZENAR TODOS OS VALORES DOS
*
* DESLOCAMENTOS EM ZFI E TODOS OS VALORES DE FORÇAS DE SUPERFÍCIE EM ZDFI
*
        NN=3*NBE
        DO 20 I=1,NN

            IF (KODE(I) == 0) THEN
                ZCH = ZFI(I)*ZGE
                ZFI(I)=ZDFI(I)
                ZDFI(I)=ZCH
            ENDIF
  20    CONTINUE
*
* CALCULA OS VALORES DOS DESLOCAMENTOS EM PONTOS INTERNOS
*
      IF(L <= 0) RETURN
*
        ALLOCATE(ZDSOL(3*L), STAT = stats1)
        ALLOCATE(ZSSOL(9*L), STAT = stats2)
        IF (stats1 /= 0 .OR. stats2 /= 0) THEN
            PRINT*, "MEMORIA INSUFCIENTE"
            STOP
        ENDIF


#ifdef USE_CPU

        ZDSOL = 0
        t1 = OMP_GET_WTIME()

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,JJ,K,KK,CO,ZHELEM,ZGELEM)
!$OMP& REDUCTION(+:ZDSOL)
        DO J=1,NBE
*
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
! Código abaixo é desnecessário, pois a função Normvec
! realiza estes calculos apenas uma vez.
C            A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
C     $          (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
C            B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
C     $          (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
C            C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
C     $          (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
C            R=DSQRT(A*A+B*B+C*C)
C      
C            ETA(1)=A/R
C            ETA(2)=B/R
C            ETA(3)=C/R
            CO(1,1)=CX(N1)
            CO(1,2)=CY(N1)
            CO(1,3)=CZ(N1)
            CO(2,1)=CX(N2)
            CO(2,2)=CY(N2)
            CO(2,3)=CZ(N2)
            CO(3,1)=CX(N3)
            CO(3,2)=CY(N3)
            CO(3,3)=CZ(N3)
            CO(4,1)=CX(N4)
            CO(4,2)=CY(N4)
            CO(4,3)=CZ(N4)
            JJ=3*(J-1)

            DO K=1,L
*
                CALL NONSINGD(ZHELEM,ZGELEM,CO,CXI(K),CYI(K),CZI(K),
     $              ETAS(1:3,J),ZGE,ZCS,ZCP,DELTA,PI,FR,GI,OME,NPG)
*
                KK=3*(K-1)

                ZDSOL(KK+1:KK+3)=ZDSOL(KK+1:KK+3) +
     $              MATMUL(ZGELEM,ZDFI(JJ+1:JJ+3))-
     $              MATMUL(ZHELEM,ZFI(JJ+1:JJ+3))

            ENDDO
        ENDDO
!$OMP END PARALLEL DO

            t2 = OMP_GET_WTIME()
            PRINT *, "INTEREC1: Tempo na CPU: ", (t2-t1)
#endif

#ifdef TEST_CUDA
           ZDSOLP = ZDSOL
#endif

#ifdef USE_GPU
            t1 = OMP_GET_WTIME()

            CALL cuda_interec1(
     $          N,
     $          NBE,
     $          NPG,
     $          L,
     $          NP,
     $          CXI,
     $          CYI,
     $          CZI,
     $          ZGE,
     $          ZCS,
     $          ZCP,
     $          0.,
     $          0.,
     $          0.,
     $          0.,
     $          FR,
     $          ZDFI,
     $          ZFI,
     $          ZDSOL,
     $          STATS1
     $        )

            t2 = OMP_GET_WTIME()

            PRINT *, "INTEREC1: Tempo na GPU: ", (t2-t1)
#endif

#ifdef TEST_CUDA
            CALL ASSERT_ZDSOL(ZDSOL, ZDSOLP, L)
#endif
*
* ACRESCENTADO POSTERIORMANTE (APÓS O PROGRAMA ESTAR RODANDO ATÉ
* O CÁLCULO PARA OS DESLOCAMENTOS EM PONTOS INTERNOS).
*
* CÁLCULO DAS TENSÕES EM PONTOS INTERNOS
*

        ZSSOL = 0
        t1 = OMP_GET_WTIME()
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,K,CO,ZD,ZS)
!$OMP& REDUCTION(+:ZSSOL)
        DO K=1,L
            DO J=1,NBE
*
                N1=CONE(J,1)
                N2=CONE(J,2)
                N3=CONE(J,3)
                N4=CONE(J,4)
C                A=(CY(N2)-CY(N1))*(CZ(N3)-CZ(N1)) - 
C     $              (CZ(N2)-CZ(N1))*(CY(N3)-CY(N1))
C                B=(CZ(N2)-CZ(N1))*(CX(N3)-CX(N1)) - 
C     $              (CX(N2)-CX(N1))*(CZ(N3)-CZ(N1))
C                C=(CX(N2)-CX(N1))*(CY(N3)-CY(N1)) - 
C     $              (CY(N2)-CY(N1))*(CX(N3)-CX(N1))
C                
C                R=DSQRT(A*A+B*B+C*C)
C                ETA(1)=A/R
C                ETA(2)=B/R
C                ETA(3)=C/R
                CO(1,1)=CX(N1)
                CO(1,2)=CY(N1)
                CO(1,3)=CZ(N1)
                CO(2,1)=CX(N2)
                CO(2,2)=CY(N2)
                CO(2,3)=CZ(N2)
                CO(3,1)=CX(N3)
                CO(3,2)=CY(N3)
                CO(3,3)=CZ(N3)
                CO(4,1)=CX(N4)
                CO(4,2)=CY(N4)
                CO(4,3)=CZ(N4)
*
                CALL SIGMAEC(CO,CXI(K),CYI(K),CZI(K),ETAS(1:3,J),DELTA,
     $              PI,FR,ZGE,RHO,ZCS,ZCP,NPG,ZD,ZS,GI,OME)
*


                ZSSOL(9*K-8)=ZSSOL(9*K-8)+ZDFI(3*J-2)*ZD(1,1,1)+
     $                                    ZDFI(3*J-1)*ZD(2,1,1)+
     $                                    ZDFI(3*J)  *ZD(3,1,1)-
     $                                    ZFI (3*J-2)*ZS(1,1,1)-
     $                                    ZFI (3*J-1)*ZS(2,1,1)-
     $                                    ZFI (3*J)  *ZS(3,1,1)
                ZSSOL(9*K-7)=ZSSOL(9*K-7)+ZDFI(3*J-2)*ZD(1,2,1)+
     $                                    ZDFI(3*J-1)*ZD(2,2,1)+
     $                                    ZDFI(3*J)  *ZD(3,2,1)-
     $                                    ZFI (3*J-2)*ZS(1,2,1)-
     $                                    ZFI (3*J-1)*ZS(2,2,1)-
     $                                    ZFI (3*J)  *ZS(3,2,1)
                ZSSOL(9*K-6)=ZSSOL(9*K-6)+ZDFI(3*J-2)*ZD(1,3,1)+
     $                                    ZDFI(3*J-1)*ZD(2,3,1)+
     $                                    ZDFI(3*J)  *ZD(3,3,1)-
     $                                    ZFI (3*J-2)*ZS(1,3,1)-
     $                                    ZFI (3*J-1)*ZS(2,3,1)-
     $                                    ZFI (3*J)  *ZS(3,3,1)
                ZSSOL(9*K-5)=ZSSOL(9*K-5)+ZDFI(3*J-2)*ZD(1,1,2)+
     $                                    ZDFI(3*J-1)*ZD(2,1,2)+
     $                                    ZDFI(3*J)  *ZD(3,1,2)-
     $                                    ZFI (3*J-2)*ZS(1,1,2)-
     $                                    ZFI (3*J-1)*ZS(2,1,2)-
     $                                    ZFI (3*J)  *ZS(3,1,2)
                ZSSOL(9*K-4)=ZSSOL(9*K-4)+ZDFI(3*J-2)*ZD(1,2,2)+
     $                                    ZDFI(3*J-1)*ZD(2,2,2)+
     $                                    ZDFI(3*J)  *ZD(3,2,2)-
     $                                    ZFI (3*J-2)*ZS(1,2,2)-
     $                                    ZFI (3*J-1)*ZS(2,2,2)-
     $                                    ZFI (3*J)  *ZS(3,2,2)
                ZSSOL(9*K-3)=ZSSOL(9*K-3)+ZDFI(3*J-2)*ZD(1,3,2)+
     $                                    ZDFI(3*J-1)*ZD(2,3,2)+
     $                                    ZDFI(3*J)  *ZD(3,3,2)-
     $                                    ZFI (3*J-2)*ZS(1,3,2)-
     $                                    ZFI (3*J-1)*ZS(2,3,2)-
     $                                    ZFI (3*J)  *ZS(3,3,2)
                ZSSOL(9*K-2)=ZSSOL(9*K-2)+ZDFI(3*J-2)*ZD(1,1,3)+
     $                                    ZDFI(3*J-1)*ZD(2,1,3)+
     $                                    ZDFI(3*J)  *ZD(3,1,3)-
     $                                    ZFI (3*J-2)*ZS(1,1,3)-
     $                                    ZFI (3*J-1)*ZS(2,1,3)-
     $                                    ZFI (3*J)  *ZS(3,1,3)
                ZSSOL(9*K-1)=ZSSOL(9*K-1)+ZDFI(3*J-2)*ZD(1,2,3)+
     $                                    ZDFI(3*J-1)*ZD(2,2,3)+
     $                                    ZDFI(3*J)  *ZD(3,2,3)-
     $                                    ZFI (3*J-2)*ZS(1,2,3)-
     $                                    ZFI (3*J-1)*ZS(2,2,3)-
     $                                    ZFI (3*J)  *ZS(3,2,3)
                ZSSOL(9*K)=  ZSSOL(9*K)  +ZDFI(3*J-2)*ZD(1,3,3)+
     $                                    ZDFI(3*J-1)*ZD(2,3,3)+
     $                                    ZDFI(3*J)  *ZD(3,3,3)-
     $                                    ZFI (3*J-2)*ZS(1,3,3)-
     $                                    ZFI (3*J-1)*ZS(2,3,3)-
     $                                    ZFI (3*J)  *ZS(3,3,3)
            ENDDO
        ENDDO        
!$OMP END PARALLEL DO
        t2 = OMP_GET_WTIME()
        PRINT*, "SIGMAEC: Tempo na CPU: ", (t2-t1)
*
        RETURN
      END

      SUBROUTINE ASSERT_ZDSOL(ZDSOL, ZDSOLP, L)
        IMPLICIT NONE
        COMPLEX, DIMENSION(3*L), INTENT(IN) :: ZDSOL, ZDSOLP
        INTEGER :: L, i
        LOGICAL :: asserted = .TRUE. 
        REAL, PARAMETER :: eps = 1.0E-6
        REAL :: maxentry = 0

        DO i = 1, 3*L
            maxentry = MAX(maxentry, ABS(ZDSOL(i) - ZDSOLP(i)))
        ENDDO

        IF (maxentry > eps) THEN
            asserted = .FALSE.
            PRINT*, "||ZDSOL||_inf = ", maxentry
        ENDIF

 200    FORMAT (A,ES7.1)     
        WRITE(0,"(A)") "O vetor ZDSOL calculado em Interec1_cu é igual "
        WRITE(0,200) "calculado em Interec.for com um erro de ", eps
        IF (asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""
      END SUBROUTINE ASSERT_ZDSOL
