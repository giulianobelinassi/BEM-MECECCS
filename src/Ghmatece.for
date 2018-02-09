************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************
*
      SUBROUTINE GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HESTdiag,GESTdiag,NE,NX,
     $         NCOX,
     $    CONE,N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS,GI,OME)
*
        USE omp_lib        

        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)

#ifdef USE_GPU
        INCLUDE 'kernels/Ghmatece_cu.fd'
#endif

        COMMON INP,INQ,IPR,IPS,IPT
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
        INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
        REAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: HESTdiag
        REAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: GESTdiag
        REAL, DIMENSION(NPG), INTENT(IN) :: GI, OME 

        REAL, INTENT(IN) :: ETAS(3,n)
        INTEGER stats1, stats2
        DIMENSION DELTA(3,3)
        DOUBLE PRECISION :: t1, t2

#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif


#ifdef  TEST_CUDA
#undef  USE_GPU
#undef  USE_CPU
#define USE_GPU
#define USE_CPU
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: HdiagP, GdiagP
#endif

#ifdef USE_CPU
        DIMENSION HELEM(3,3),GELEM(3,3)
        DIMENSION CO(4,3)
        REAL, DIMENSION(:,:), ALLOCATABLE :: HEST
#endif
        PI=4.0*ATAN(1.0)


* CONSTANTES USADAS NAS SOLUÇÕES FUNDAMENTAIS ESTÁTICAS
*
        C1=1.0/(16.0*PI*RMU*(1.0-RNU))
        C2=3.0-(4.0*RNU)
        C3=-1.0/(8.0*PI*(1.0-RNU))
        C4=1.0-(2.0*RNU)
*
        DO 25 J=1,3
            DO 25 I=1,3
                IF (I == J) THEN
                    DELTA(I, J) = 1.0
                ELSE
                    DELTA(I, J) = 0.0
                ENDIF
  25    CONTINUE
*
* ZERANDO AS MATRIZES H E G
*
        ALLOCATE(GESTdiag(3,3,NBE), STAT = stats1)
        ALLOCATE(HESTdiag(3,3,NBE), STAT = stats2)
        
        IF (stats1 /= 0 .or. stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE!"
            STOP
        ENDIF

*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*

#ifdef USE_CPU
        t1 = OMP_GET_WTIME() 
        ALLOCATE(HEST(3*NBE, 3*N), STAT = stats1)        
        HEST = 0
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ,HELEM,GELEM)
        DO J=1,N
*
* CÁLCULO DAS COMPONENTES DO VETOR NORMAL
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
C            ETAS(1)=A/R
C            ETAS(2)=B/R
C            ETAS(3)=C/R
            
*
* ARMAZENA AS COORDENADAS DOS PONTOS EXTREMOS DO ELEMENTO NO VETOR CO
*
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
            
            JJ=3*(J-1) + 1
            DO I=1,NBE
*    
                IF (I == J) THEN
                    HEST(JJ:JJ+2, JJ:JJ+2) = 0
                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H NÃO SINGULAR
                    
                    CALL NONSINGE(HELEM,CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),N,NP,NPG,RNU,RMU,C3,C4,
     $                  DELTA,GI,OME)

                    II=3*(I-1) + 1
                    HEST(II:II+2, JJ:JJ+2) = HELEM
                ENDIF

*
            ENDDO
        ENDDO
!$OMP END PARALLEL DO

!Calcula Gest singular 

            CALL GEST_SINGULAR(
     $          GESTdiag,
     $          CXM,
     $          CYM, 
     $          CZM,
     $          ETAS,
     $          CX,
     $          CY,
     $          CZ,
     $          NCOX,
     $          N,
     $          NP,
     $          NPG,
     $          C1,
     $          C2,
     $          DELTA,
     $          GI,
     $          OME,
     $          CONE,
     $          NBE
     $          )



*
* ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H SINGULAR ATRAVÉS
* DA CONSIDERAÇÃO DO MOVIMENTO DO CORPO RÍGIDO
*
        
        CALL RIGID_BODY(NBE, N, HEST, HESTdiag)
        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECE: Tempo na CPU: ", (t2-t1)

        DEALLOCATE(HEST)
#endif

#ifdef TEST_CUDA
        ALLOCATE(HdiagP(3,3,NBE))
        ALLOCATE(GdiagP(3,3,NBE))
        HdiagP = HESTdiag
        GdiagP = GESTdiag

#endif

#ifdef USE_GPU
        HESTdiag = 0
        GESTdiag = 0
        t1 = OMP_GET_WTIME()
!$OMP PARALLEL NUM_THREADS(2)

        IF (OMP_GET_THREAD_NUM() == 0) THEN
            CALL cuda_ghmatece(
     $          NBE,
     $          NPG,
     $          N,
     $          NP,
     $          C3,
     $          C4,
     $          FR,
     $          HESTdiag,
     $          STATS1
     $          )
        ELSE

            CALL GEST_SINGULAR(
     $          GESTdiag,
     $          CXM,
     $          CYM, 
     $          CZM,
     $          ETAS,
     $          CX,
     $          CY,
     $          CZ,
     $          NCOX,
     $          N,
     $          NP,
     $          NPG,
     $          C1,
     $          C2,
     $          DELTA,
     $          GI,
     $          OME,
     $          CONE,
     $          NBE
     $          )

            CALL CUDA_SEND_GEST_DATA(NBE, GESTdiag) 
        ENDIF
!$OMP END PARALLEL
        
        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECE: Tempo na GPU: ", (t2-t1)
#endif


#ifdef TEST_CUDA
        CALL ASSERT_GHMATECE_H_G(HESTdiag, HdiagP,GESTdiag,GdiagP,NBE,N)
        DEALLOCATE(HdiagP)
        DEALLOCATE(GdiagP)
#endif

        RETURN
      END SUBROUTINE GHMATECE

      SUBROUTINE RIGID_BODY(NBE, N, HEST, HESTdiag)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NBE, N
        REAL, DIMENSION(3*NBE, 3*N), INTENT(INOUT) :: HEST
        REAL, DIMENSION(3,3,NBE), INTENT(OUT) :: HESTdiag
        INTEGER :: MA, MB, II, JJ

!$OMP PARALLEL DO PRIVATE(MA, II)
        DO MA=1, NBE
            II = 3*(MA-1)+1
            HESTdiag(1:3, 1:3, MA) = 0
        ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(MA, MB, II, JJ)
        DO MA=1,NBE
            II = 3*(MA-1)+1
            DO MB=1, N
                JJ = 3*(MB-1)+1
                HESTdiag(1:3, 1:3, MA) = HESTdiag(1:3, 1:3, MA) -
     $              HEST(II:II+2, JJ:JJ+2) 
            ENDDO
        ENDDO
!$OMP END PARALLEL DO
      END


!     Note que há cálculos muito diferentes dos demais no problema
!     singular, e portanto decidi (Giuliano) tratá-lo de maneira 
!     diferenciada.
      SUBROUTINE GEST_SINGULAR(GESTdiag, CXM, CYM, CZM, ETAS,
     $      CX, CY, CZ, NCOX, N, NP, NPG,
     $      C1, C2, DELTA, GI, OME, CONE, NBE)
       
        IMPLICIT NONE
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
        INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
        REAL, DIMENSION(3,3,NBE), INTENT(OUT) :: GESTdiag
        REAL, DIMENSION(NPG) :: GI, OME
        REAL, INTENT(IN) :: ETAS(3,n)
        INTEGER, INTENT(IN) :: NBE, N, NP, NPG, NCOX
        REAL, INTENT(IN) :: C1,C2
        REAL, INTENT(IN) :: DELTA(3,3)
        INTEGER :: J, JJ, N1, N2, N3, N4

        GESTdiag = 0

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,JJ)
        DO J=1,NBE
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
            
            JJ=3*(J-1) + 1
            CALL SINGGE(GESTdiag(1:3,1:3,J),
     $                  CXM(J),CYM(J),CZM(J),
     $                  ETAS(1:3,J),CX,
     $                  CY,CZ,N1,N2,N3,N4,NCOX,N,NP,NPG,C1,
     $                  C2,DELTA,GI,OME)
        ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GEST_SINGULAR

      SUBROUTINE ASSERT_GHMATECE_H_G(Hdiag, HdiagP, Gdiag, GdiagP,NBE,N)
        IMPLICIT NONE
        REAL :: eps = 1.0E-5
        REAL, INTENT(IN), DIMENSION(3,3,NBE) ::Hdiag,HdiagP,Gdiag,Gdiagp
        INTEGER, INTENT(IN) :: NBE, N  

        INTEGER :: i, j, k
        LOGICAL :: ghmatecd_asserted = .TRUE.
        REAL :: sum_norms = 0, local_sum, max_local_sum

        eps = 1.0E-5

        sum_norms = 0.0
        max_local_sum = 0.0

        DO k = 1, NBE
            local_sum = 0
            DO j = 1, 3
                DO i = 1, 3
                    local_sum = local_sum + ABS(HdiagP(i, j, k) - 
     $                       Hdiag(i, j, k))
                ENDDO
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
            if (max_local_sum > 1.0E10) then
                print*, "ERROU:", k
            endif

        ENDDO

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
        ENDIF

        PRINT*, "||[H]est|| = ", max_local_sum
        
        sum_norms = 0.0
        max_local_sum = 0.
        
        DO k = 1, NBE
            local_sum = 0
            DO j = 1, 3
                DO i = 1, 3
                    local_sum = local_sum + ABS(GdiagP(i, j, k) - 
     $                       Gdiag(i, j, k))
                ENDDO
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

        PRINT*, "||G||_inf = ", max_local_sum

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
        ENDIF

 200    FORMAT (A,ES7.1)     
       WRITE(0,"(A)") "As matrizes H e G em Ghmatecd_cu sao iguais as"
        WRITE(0,200) "calculadas em Ghmatece com um erro de ", eps
        IF (ghmatecd_asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""

      END SUBROUTINE ASSERT_GHMATECE_H_G
