************************************************************************
*                                                                      *
*         ESTA SUBROTINA CALCULA AS MATRIZES [H] E [G];                *
*                                                                      *
*        AINDA MONTA O SISTEMA DE EQUAÇÕES [A] {x} = {f}               *
*                                                                      *
************************************************************************



      SUBROUTINE GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HESTdiag,GESTdiag,ZH,
     $   KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,ZDFI,ZFI,
     $    L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS,GI,OME)
        
        USE omp_lib


        IMPLICIT NONE

#define USE_CPU
#ifdef USE_GPU
#undef USE_CPU
#endif

       
#ifdef USE_GPU
        INCLUDE 'kernels/Ghmatecd_cu.fd'
#endif
 
! O compilador reclama que eu estou passando parâmetros float para a 
! ZGEMV. Acontece que no caso de precisao simples esta chamada nunca
! é feita.
!#ifdef USE_CPU
!		INCLUDE 'Lapack.fd'
!#endif

 
        INTEGER INP,INQ,IPR,IPS,IPT
        COMMON  INP,INQ,IPR,IPS,IPT
        

        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N), INTENT(IN) :: CXM, CYM, CZM
        REAL, DIMENSION(3,3,NBE), INTENT(IN) :: HESTdiag, GESTdiag
        COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ZH
        INTEGER, INTENT(IN) :: KODE(3*NBE),NE,NX,NCOX,CONE(N,4)
        REAL, INTENT(IN) :: DELTA(3,3),PI
        INTEGER, INTENT(IN) :: L,N,NBE,NP,NPG
        REAL, INTENT(IN) :: GE,RNU,RMU,FR,DAM,RHO
        COMPLEX, INTENT(IN) :: ZGE,ZCS,ZCP
        REAL, INTENT(IN) :: C1,C2,C3,C4
        REAL, INTENT(IN) :: ETAS(3,N)
        COMPLEX, INTENT(IN) :: ZDFI(3*NBE)
        COMPLEX, INTENT(OUT), ALLOCATABLE :: ZFI(:)


        REAL, INTENT(IN) :: GI(NPG), OME(NPG)
        DOUBLE PRECISION :: t1, t2
        INTEGER :: stats1

#ifdef TEST_CUDA
#undef  USE_CPU
#undef  USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX, ALLOCATABLE :: ZHP(:,:), ZFIP(:)
#endif

#ifdef  USE_CPU
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ZG
        INTEGER I,J, stats2, II, JJ, NN
        INTEGER N1,N2,N3,N4
        REAL :: CO(4,3)
        INTEGER(kind=8) NNl
        INTEGER(kind=8) :: one = INT(1, 8)
#endif

#ifdef USE_GPU
        INTEGER RET
#endif

        ALLOCATE(ZFI(3*NBE) , STAT = stats1)
        IF (stats1 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
            STOP
        ENDIF
        
        ALLOCATE(ZH(3*NBE, 3*NBE), STAT = stats1)
        IF (stats1 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
            STOP
        ENDIF


#ifdef USE_CPU

*
* CÁLCULO DOS COEFICIENTES DAS MATRIZES H E G
*
*
        NN = 3*NBE

        ALLOCATE(ZG(3*NBE, 3*NBE), STAT = stats2)

        IF (stats2 /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE"
            STOP
        ENDIF

        t1 = OMP_GET_WTIME()


!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ)
        DO J=1,NBE
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
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
                II=3*(I-1) + 1

                IF (I == J) THEN
C                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G SINGULAR
C                   ATRAVÉS DA DIFERENÇA DINÂMICO - ESTÁTICO

                    CALL SING_DE (ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),
     $                  ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,GI,OME,NPG)


                    ZG(II:II+2, JJ:JJ+2) = ZG(II:II+2, JJ:JJ+2) +
     $                  GESTdiag(1:3, 1:3, I)
                    ZH(II:II+2, JJ:JJ+2) = ZH(II:II+2, JJ:JJ+2) + 
     $                  HESTdiag(1:3, 1:3, I)

                ELSE
*                   ACIONA ROTINA QUE CALCULA OS COEFICIENTES DE H E G NÃO SINGULAR
*
                    CALL NONSINGD (ZH(II:II+2, JJ:JJ+2),
     $                  ZG(II:II+2, JJ:JJ+2),
     $                  CO,CXM(I),CYM(I),CZM(I),
     $                  ETAS(1:3,J),
     $                  ZGE,ZCS,ZCP,DELTA,PI,FR,GI,OME,NPG)
                ENDIF
            ENDDO
        ENDDO
!$OMP END PARALLEL DO

        CALL BOUNDARY_CONDITIONS(ZH, ZG, KODE, ZGE, NN)

        NNl = NN

        IF (SIZEOF(1.0) == 8) THEN
            CALL ZGEMV('N',NNl,NNl,(1.0,0),ZG,NNl,ZDFI,one,(0,0),
     $          ZFI,one)
        ELSEIF (SIZEOF(1.0) == 4) THEN
            CALL CGEMV('N',NNl,NNl,(1.0,0),ZG,NNl,ZDFI,one,(0,0),
     $          ZFI,one)
        ELSE
            PRINT*, "ERRO FATAL: Precisão desconhecida"
            STOP
        ENDIF

        t2 = OMP_GET_WTIME()
        PRINT *, "GHMATECD: Tempo na CPU: ", (t2-t1)
#endif

#ifdef TEST_CUDA
!FAÇA UMA CÓPIA DAS MATRIZES ZG E ZH PARA COMPARAÇÃO COM O RESULTADO DA GPU.
        ALLOCATE(ZHP(3*NBE, 3*NBE))
        ALLOCATE(ZFIP(3*NBE))
        ZFIP = ZFI
        ZHP = ZH
#endif

#ifdef USE_GPU
        t1 = OMP_GET_WTIME()

        CALL cuda_ghmatecd(
     $        NBE,
     $        NPG,
     $        N,
     $        NP,
     $        ZGE,
     $        ZCS,
     $        ZCP,
     $        C1,
     $        C2,
     $        C3,
     $        C4,
     $        FR,
     $        HESTdiag,
     $        GESTdiag,
     $        ZFI,
     $        ZH,
     $        KODE,
     $        1,
     $        RET
     $        )

        IF (RET /= 0) THEN
            PRINT*, "GHMATECD: Erro: Matriz Singular."
        ENDIF

        t2 = OMP_GET_WTIME()
    
        PRINT *, "GHMATECD: Tempo na GPU: ", (t2 - t1)
#endif
        
#ifdef TEST_CUDA
        CALL ASSERT_GHMATECD_ZH_ZFI(ZH, ZHP, ZFI, ZFIP, NBE, NBE)
        DEALLOCATE(ZHP)
        DEALLOCATE(ZFIP)
#endif
#ifdef USE_CPU
        DEALLOCATE(ZG)
#endif

        RETURN
      END

      SUBROUTINE GHMATECD_SINGULAR(ZHdiag, ZGdiag, CX, CY, CZ, ZGE, ZCS,
     $         ZCP,
     $   C1, C2, C3, C4, DELTA, FR, GI, OME, NPG, N, NBE, NP, ETAS, 
     $   CXM, CYM, CZM, PI, GESTdiag, HESTdiag, CONE)
        
        IMPLICIT NONE 

        COMPLEX, DIMENSION(3,3,NBE),INTENT(OUT) :: ZHdiag, ZGdiag
        REAL, DIMENSION(3,3,NBE), INTENT(IN) :: HESTdiag, GESTdiag
        
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        COMPLEX,   INTENT(IN) :: ZGE,ZCS,ZCP
        REAL, INTENT(IN) :: C1,C2,C3,C4,FR,PI
        REAL, INTENT(IN) :: DELTA(3,3)
        REAL, INTENT(IN) :: GI(NPG), OME(NPG)
        INTEGER, INTENT(IN) :: N,NBE,NP,NPG
        REAL, INTENT(IN) :: ETAS(3,N)
        REAL, DIMENSION(N), INTENT(IN) :: CXM, CYM, CZM
        INTEGER, DIMENSION(N, 4) :: CONE

        INTEGER :: N1, N2, N3, N4, I, J, II, JJ
        REAL :: CO(4,3)

        ZGdiag = 0
        ZHdiag = 0

!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(N1,N2,N3,N4,J,I,CO,II,JJ)
        DO J=1,NBE
            N1=CONE(J,1)
            N2=CONE(J,2)
            N3=CONE(J,3)
            N4=CONE(J,4)
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
            I = J
            II=3*(I-1) + 1


            CALL SING_DE (ZHdiag(1:3, 1:3, J),
     $          ZGdiag(1:3, 1:3, J),
     $          CO,CXM(I),CYM(I),CZM(I),
     $          ETAS(1:3,J),
     $          ZGE,ZCS,ZCP,C1,C2,C3,C4,DELTA,PI,FR,GI,OME,NPG)


            ZGdiag(1:3, 1:3, J) = ZGdiag(1:3, 1:3, J) +
     $          GESTdiag(1:3, 1:3, J)
            ZHdiag(1:3, 1:3, J) = ZHdiag(1:3, 1:3, J) + 
     $          HESTdiag(1:3, 1:3, J)

        ENDDO
!$OMP END PARALLEL DO
      END SUBROUTINE GHMATECD_SINGULAR

*
* REORDENA AS COLUNAS DO SISTEMA DE EQUAÇÕES DE ACORDO COM AS
*
* CONDIÇÕES DE CONTORNO E FORMA A MATRIZ A QUE É ARMAZENADA EM G
*

      SUBROUTINE BOUNDARY_CONDITIONS(ZG, ZH, KODE, ZGE, NN)
        IMPLICIT NONE
        COMPLEX, INTENT(INOUT) :: ZH(NN, NN), ZG(NN, NN)
        COMPLEX, INTENT(IN) :: ZGE
        INTEGER, INTENT(IN) :: NN, KODE(NN)
        COMPLEX :: ZCH
        INTEGER :: I, J

        DO  J=1,NN
            IF (KODE(J) == 0) THEN
                DO I=1,NN
                    ZCH = ZG(I,J)*ZGE
                    ZG(I,J) = -ZH(I,J)
                    ZH(I,J) = -ZCH
                ENDDO
            ENDIF
        ENDDO
      END SUBROUTINE BOUNDARY_CONDITIONS


      SUBROUTINE ASSERT_GHMATECD_ZH_ZFI(ZH, ZHP, ZFI, ZFIP, NBE, N)
        IMPLICIT NONE
        COMPLEX, INTENT(IN), DIMENSION(3*NBE,3*N) :: ZH,ZHP
        COMPLEX, INTENT(IN), DIMENSION(3*NBE) :: ZFI, ZFIP
        INTEGER, INTENT(IN) :: NBE, N  

        INTEGER :: i, j
        INTEGER :: N3, NBE3
        LOGICAL :: ghmatecd_asserted = .TRUE.
        REAL :: sum_norms = 0, eps
        REAL :: local_sum = 0, max_local_sum = 0

        eps = 1.5E-6*nbe
        N3   = N*3
        NBE3 = NBE*3 

        sum_norms = 0.0
        

        DO j = 1, N3
            local_sum = 0
            DO i = 1, NBE3
                local_sum = local_sum + ABS(ZHP(i,j) - ZH(i,j))
            ENDDO
            sum_norms = sum_norms + local_sum
            max_local_sum = MAX(local_sum, max_local_sum)
        ENDDO

!        DO j = 1, N3
!            DO i = 1, NBE3
!                local_sum = ABS(ZHP(i,j) - ZH(i,j))
!                if (local_sum > max_local_sum) then
!                    max_local_sum = local_sum
!                endif
!            ENDDO
!        ENDDO
        PRINT*, max_local_sum
       

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
        ENDIF

        PRINT*, "||ZH||_inf = ", max_local_sum

        sum_norms = 0.0
        max_local_sum = 0
        
        DO j = 1, N3
            sum_norms = ABS(ZFI(j) - ZFIP(j))
            max_local_sum = MAX(sum_norms, max_local_sum)
        ENDDO

        PRINT*, "MAX ||ZFI|| = ", max_local_sum

        IF (max_local_sum > eps) THEN
            ghmatecd_asserted = .FALSE.
        ENDIF

 200    FORMAT (A,ES7.1)     
       WRITE(0,"(A)") "As matrizes ZH e ZG em Ghmatecd_cu sao iguais as"
        WRITE(0,200) "calculadas em Ghmatecd com um erro de ", eps
        IF (ghmatecd_asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""

      END SUBROUTINE ASSERT_GHMATECD_ZH_ZFI
