************************************************************************
************************************************************************
**                                                                    **
**     ESTE PROGRAMA RESOLVE PROBLEMAS VISCOELASTICOS, DINÂMICOS,     **
**                                                                    **
**   TRIDIMENSIONAIS ATRAVÉS DO MÉTODO DOS ELEMENTOS DE CONTORNO NUM  **
**                                                                    **
**                      SEMI ESPAÇO INFINITO                          **
**                                                                    **
************************************************************************
**                                                                    **
** PROGRAMA "MECECCS" = (MÉTODO DOS ELEMENTOS DE CONTORNO USANDO      **
**                                                                    **
**                     ELEMENTOS CONSTANTES COM SINGULARIDADE)        **
**                                                                    **
**                              DINÂMICO                              **
**                                                                    **
** RONALDO CARRION / EUCLIDES DE MESQUITA NETO                        **
**                                                                    **
** DATA: 18/07/01                                                     **
**                                                                    **
** MODIFICADO POR: GIULIANO AUGUSTO FAULIN BELINASSI                  ** 
** DATA: 12/2016 - 01/2018                                            **
************************************************************************
************************************************************************
*
      PROGRAM MAIN
*
        USE omp_lib

        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        
        INCLUDE 'Inputece.fd'
        INCLUDE 'Normvec.fd'
        INCLUDE 'Ghmatece.fd'
        INCLUDE 'Inputecd.fd'
        INCLUDE 'Ghmatecd.fd'
        INCLUDE 'Interec.fd'
        INCLUDE 'Outputec.fd'
        INCLUDE 'Linsolve.fd'

#ifdef USE_GPU
        INCLUDE 'kernels/shared.fd'
#endif

! Talvez estes parâmetros sejam desnecessários uma vez que nesta versão
! todas as matrizes são alocadas dinamicamente de acordo com a entrada
! do problema.
C       PARAMETER (NE=960,NX=3*NE,NCOX=962,NPIX=10,NFRX=7)
        PARAMETER (NE=2200,NX=3*NE,NCOX=2200,NPIX=10,NFRX=7)
        COMMON INP,INQ,IPR,IPS,IPT
        REAL, DIMENSION(:), ALLOCATABLE :: CX, CY, CZ, CXM, CYM, CZM
        INTEGER, ALLOCATABLE :: CONE(:,:), KODE(:)
        REAL, DIMENSION(:), ALLOCATABLE :: CXI,CYI, CZI
        REAL, DIMENSION(:), ALLOCATABLE :: BC, AFR
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: HESTdiag, GESTdiag
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ZH
        COMPLEX, DIMENSION(:), ALLOCATABLE :: ZDFI, ZFI, ZDSOL, ZSSOL
    
        REAL, DIMENSION(:), ALLOCATABLE :: GI, OME


        DIMENSION DELTA(3,3)
        REAL, DIMENSION(:,:), ALLOCATABLE :: ETAS
        INTEGER i
        CHARACTER(len=100) :: input_e, input_d, output_e, output_d
#ifdef USE_GPU
        DOUBLE PRECISION :: t1, t2
#endif
*
* NE = NÚMERO MÁXIMO DE ELEMENTOS DA MALHA (CONTORNO + ENCLOSING)
* NX = DIMENSÃO MÁXIMA DO SISTEMA DE EQUAÇÕES       
* NCOX = NÚMERO MÁXIMO DE PONTOS EXTREMOS DOS ELEMENTOS (NE+2)
* NPIX = NÚMERO MÁXIMO DE PONTOS INTERNOS
*
* PARÂMETROS ESPECIFICADORES DOS NÚMEROS DE UNIDADE (ARQUIVO) DE E/S
*
        INP=100
        INQ=101
        IPR=102
        IPS=103
        IPT=104
*
        IF (iargc() >= 4) THEN
            CALL getarg(1, input_e)
            CALL getarg(2, input_d)
            CALL getarg(3, output_e)
            CALL getarg(4, output_d)
        ELSE
!           Parâmetros estranhos, compute o caso abaixo:
            input_e  = 'ESOLO240E_-5+5.DAT'
            input_d  = 'ESOLO240D_-5+5.DAT'
            output_e = 'SSOLO240E_-5+5.DAT'
            output_d = 'SSOLO240D_-5+5.DAT'
        ENDIF

        OPEN(UNIT=INP,FILE=input_e,STATUS='OLD')
        OPEN(UNIT=INQ,FILE=input_d,STATUS='OLD')
        OPEN(UNIT=IPR,FILE=output_e,STATUS='UNKNOWN')
        OPEN(UNIT=IPS,FILE=output_d,STATUS='UNKNOWN')

*
*
* ACIONA ROTINA QUE LÊ OS DADOS DE ENTRADA PARA O PROBLEMA ESTÁTICO
*
      
        CALL INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $      N,NBE,NP,NPG,GE,RNU,RMU)

! ACIONA A SUBROTINA QUE CALCULA OS PONTOS E PESOS DE GAUSS.        
        ALLOCATE(GI(NPG))
        ALLOCATE(OME(NPG))
        CALL GAULEG(-1.0, 1.0, GI, OME, NPG)


! Aciona a rotina que calcula as normas originalmente usadas em
! Ghmatece.for Ghmatecd.for e Interec.for. Isto evita calculos
! reduntantes.
        ETAS = NORMVEC(CONE, CX, CY, CZ, N, NP)

! Aciona a rotina que envia dados que são usados em diversas subrotinas
! para a GPU
#ifdef TEST_CUDA
#define USE_GPU
#endif
#ifdef USE_GPU
        t1 = OMP_GET_WTIME()

        CALL send_shared_data_to_gpu(CX,CY,CZ,CXM,CYM,CZM,ETAS,GI,OME, 
     $      CONE,NP,NPG,N,NBE) 

        t2 = OMP_GET_WTIME()

      PRINT*, "Tempo gasto alocando os vetores compartilhados: ",(t2-t1)
#endif

*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA ESTÁTICO
*
        CALL GHMATECE (CX,CY,CZ,CXM,CYM,CZM,HESTdiag,GESTdiag,NE,NX,
     $      NCOX,CONE,
     $      N,NBE,NP,NPG,GE,RNU,RMU,DELTA,PI,C1,C2,C3,C4,ETAS,GI,OME)
*
* ACIONA ROTINA QUE LÊ OS DADOS DE ENTRADA PARA O PROBLEMA DINÂMICO
*
        CALL INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $     NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $     L,FR,DAM,RHO,ZGE,ZCS,ZCP)
       
        
        NN=3*NBE
        ALLOCATE(ZDFI(NN), STAT = ITER)
        IF (ITER /= 0) THEN
            PRINT*, "Memória insuficiente"
        ENDIF

        ZDFI = BC
        
#ifdef USE_GPU
        t1 = OMP_GET_WTIME()
        CALL send_dyn_data_gpu(NBE, ZDFI, KODE)
        t2 = OMP_GET_WTIME()
        PRINT*, "Vetores compartilhados do prob. dynamico: ", (t2-t1)
#endif

        DO 12 I=1,NFR
            FR=AFR(I)
!            PRINT *, 'rodando para frequencia ', i

*
* ACIONA ROTINA QUE CALCULA AS MATRIZES [H] E [G] DO PROBLEMA ESTÁTICO
*
            CALL GHMATECD (CX,CY,CZ,CXM,CYM,CZM,HESTdiag,GESTdiag,ZH,
     $         KODE,NE,NX,NCOX,CONE,DELTA,PI,N,NBE,NP,NPG,GE,RNU,RMU,
     $         ZDFI,ZFI,L,FR,DAM,RHO,ZGE,ZCS,ZCP,C1,C2,C3,C4,ETAS,GI,OME
     $      )
*
* ACIONA ROTINA QUE RESOLVE O SISTEMA DE EQUAÇÕES

            CALL LINSOLVE(NN, N, ZH, ZFI)
            DEALLOCATE(ZH)
*
* ACIONA ROTINA QUE CALCULA OS DESLOCAMENTOS EM PONTOS INTERNOS
*
            CALL INTEREC(ZFI,ZDFI,KODE,CX,CY,CZ,CXI,CYI,CZI,ZDSOL,ZSSOL,
     $          NE,NX,NCOX,NPIX,CONE,ZGE,ZCS,ZCP,DELTA,PI,FR,NPG,L,N,
     $          NBE,RHO,ETAS,GI,OME,NP
     $      )
     
*
* ACIONA ROTINA QUE IMPRIME AS VARIÁVEIS NO CONTORNO E PONTOS INTERNOS
*
            CALL OUTPUTEC(ZFI,ZDFI,ZDSOL,ZSSOL,NPIX,NX,N,NBE,L,FR)


            DEALLOCATE(ZFI)
            DEALLOCATE(ZDSOL)
            DEALLOCATE(ZSSOL)
 12     CONTINUE
*
        DEALLOCATE(ZDFI)

        DEALLOCATE(GI)
        DEALLOCATE(OME)

        DEALLOCATE(CX)
        DEALLOCATE(CY)
        DEALLOCATE(CZ)

        DEALLOCATE(CXM)
        DEALLOCATE(CYM)
        DEALLOCATE(CZM)
       
        DEALLOCATE(CXI)
        DEALLOCATE(CYI)
        DEALLOCATE(CZI)

        DEALLOCATE(BC)
        DEALLOCATE(KODE)
        DEALLOCATE(AFR)

        DEALLOCATE(CONE) 
        
        DEALLOCATE(ETAS)

        DEALLOCATE(HESTdiag)
        DEALLOCATE(GESTdiag)
    
#ifdef USE_GPU
        CALL deallocate_shared_gpu_data() 
#endif

        CLOSE (INP)          
        CLOSE (INQ)
        CLOSE (IPR)          
        CLOSE (IPS)
        STOP
      END
*
*#######################################################################
*#######################################################################
*                   FIM DO PROGRAMA PRINCIPAL                          #
*#######################################################################
*#######################################################################
