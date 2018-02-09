      SUBROUTINE LINSOLVE(NN, N, ZH, ZFI)
        USE omp_lib
    
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NN, N
        COMPLEX, INTENT(INOUT) ::  ZH(NN, NN)
        COMPLEX, INTENT(INOUT) :: ZFI(NN)
        DOUBLE PRECISION :: t1, t2
#define USE_CPU

#ifdef USE_GPU
#undef USE_CPU
#endif

#ifdef  TEST_CUDA
#undef  USE_CPU
#undef  USE_GPU
#define USE_CPU
#define USE_GPU
        COMPLEX :: ZFI_ORIG(NN), ZFIP(NN)
        COMPLEX, ALLOCATABLE :: ZH_ORIG(:,:), ZHP(:,:)
#endif


! FORMA O LADO DIREITO DO SISTEMA {VETOR f} QUE É ARMAZENADO EM ZFI

#ifdef TEST_CUDA
        ALLOCATE(ZH_ORIG(NN, NN))
        ALLOCATE(ZHP(NN, NN))
        ZH_ORIG = ZH
        ZFI_ORIG = ZFI
#endif

#ifdef USE_CPU

        t1 = OMP_GET_WTIME()
    
        CALL LINSOLVE_CPU(NN, ZH, ZFI)
     
        t2 = OMP_GET_WTIME()
        PRINT*, "LINSOLVE: Tempo na CPU: ", (t2-t1)

#endif

#ifdef  TEST_CUDA
        
        ZFIP = ZFI
        ZH = ZHP
        ZFI = ZFI_ORIG
        ZH = ZH_ORIG
#endif

#ifdef USE_GPU
        t1 = OMP_GET_WTIME()  
        CALL cuda_linsolve(NN, N, ZH, ZFI)
        t2 = OMP_GET_WTIME()
        PRINT*, "LINSOLVE: Tempo na GPU: ", (t2-t1)
#endif

#ifdef  TEST_CUDA
        CALL ASSERT_ZFI(ZFI, ZFIP, NN)
        DEALLOCATE(ZHP)
        DEALLOCATE(ZH_ORIG)
#endif

      END SUBROUTINE

      SUBROUTINE LINSOLVE_CPU(NN, ZH, ZFI)
        IMPLICIT NONE


! O compilador reclama que eu estou passando parâmetros float para a 
! ZGESV. Acontece que no caso de precisao simples esta chamada nunca
! é feita.
!		INCLUDE 'Lapack.fd'

        INTEGER, INTENT(IN) :: NN
        COMPLEX, INTENT(INOUT) ::  ZH(NN, NN)
        COMPLEX, INTENT(INOUT) :: ZFI(NN)
        INTEGER(kind=8), ALLOCATABLE :: PIV(:)
        INTEGER :: stats
        INTEGER(kind=8) :: stats8

        INTEGER(kind=8) :: NNl, one
        NNl = INT(NN, 8) 
        one = INT(1, 8) 

        ALLOCATE(PIV(NN), STAT = stats)
        IF (stats /= 0) THEN
            PRINT*, "MEMORIA INSUFICIENTE"
            STOP
        ENDIF

        IF (SIZEOF(1.0) == 8) THEN
            CALL ZGESV(NNl,one,ZH,NNl,PIV,ZFI,NNl,stats8)
        ELSEIF (SIZEOF(1.0) == 4) THEN
            CALL CGESV(NNl,one,ZH,NNl,PIV,ZFI,NNl,stats8)
        ELSE
            PRINT*, "ERRO FATAL: Precisão de float desconhecida"
        ENDIF

        IF (stats < 0) THEN
            PRINT *, "Erro em ZGESV :-("
        ELSE IF (stats > 0) THEN
            PRINT *, "Matriz Singular :-|"
        ENDIF
        DEALLOCATE(PIV)

      END SUBROUTINE LINSOLVE_CPU


      SUBROUTINE ASSERT_ZFI(ZFI, ZFIP, NN)
        IMPLICIT NONE
        COMPLEX, DIMENSION(NN), INTENT(IN) :: ZFI, ZFIP
        INTEGER :: NN, i
        LOGICAL :: asserted = .TRUE. 
        REAL, PARAMETER :: eps = 1.2E-6
        REAL :: maxentry = 0

        DO i = 1, NN
            maxentry = MAX(maxentry, ABS(ZFI(i) - ZFIP(i)))
!            PRINT*, ZFI(i), ZFIP(i)
        ENDDO

        IF (maxentry > eps) THEN
            asserted = .FALSE.
        ENDIF

        PRINT*, "||ZFI||_inf = ", maxentry


 200    FORMAT (A,ES7.1)     
        WRITE(0,"(A)") "O vetor ZFI calculado em Interec1_cu e igual ao"
        WRITE(0,200) "calculado em Interec.for com um erro de ", eps
        IF (asserted .EQV. .TRUE.) THEN
            WRITE(0,"(A)")"[OK]"
        ELSE
            WRITE(0,"(A)")"[FALHOU]"
        ENDIF
        WRITE(0,"(A)") ""
      END SUBROUTINE ASSERT_ZFI
