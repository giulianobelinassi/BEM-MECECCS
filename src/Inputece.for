************************************************************************
*                                                                      *
*      ESTA SUBROTINA LÊ E IMPRIME OS DADOS DO ARQUIVO DE ENTRADA      *
*                                                                      *
************************************************************************

      SUBROUTINE ALLOCATE_ASSERT(stats)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: stats
        IF (stats /= 0) THEN
            PRINT*, "MEMÓRIA INSUFICIENTE!"
            STOP
        ENDIF
      END
*
      SUBROUTINE INPUTECE (CX,CY,CZ,NE,NCOX,CONE,CXM,CYM,CZM,
     $    N,NBE,NP,NPG,GE,RNU,RMU)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        CHARACTER*75 TITULO
        REAL, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: CX, CY, CZ, CXM
        REAL, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: CYM, CZM
        INTEGER, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE :: CONE
        REAL cone_in(4)
        INTEGER stats
*
        WRITE(IPR,100)
 100    FORMAT(/' ',79('*')/)
*
* LEITURA DO TÍTULO DO TRABALHO
*
        READ(INP,'(A)')TITULO
        WRITE(IPR,'(A)')TITULO
*
* LEITURA DO NÚMERO DE NÓS E DAS PROPRIEDADES DO MATERIAL
*

        READ(INP,*)N,NBE,NP,NPG,GE,RNU 
        WRITE(IPR,200)N,NBE,NP,NPG,GE,RNU     
 200    FORMAT(//'DADOS'//2X,'NÚMERO DE ELEMENTOS DA MALHA= ',I4/
     $      2X,'NÚMERO DE ELEMENTOS DE CONTORNO= 'I4/
     $      2X,'NÚMERO DE PONTOS EXTREMOS DOS ELEMENTOS= 'I4/
     $      2X,'NÚMERO DE PONTOS DE GAUSS= ',I3/    
     $      2X,'MÓDULO DE CISALHAMENTO= 'D14.7/
     $      2X,'COEFICIENTE DE POISSON= ',D14.7)
*
        RMU=GE

        ALLOCATE(CX(NP), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)
        ALLOCATE(CY(NP), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)
        ALLOCATE(CZ(NP), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)


        ALLOCATE(CXM(N), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)
        ALLOCATE(CYM(N), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)
        ALLOCATE(CZM(N), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)

        ALLOCATE(CONE(N, 4), STAT = stats)
        CALL ALLOCATE_ASSERT(stats)

*
* LEITURA DAS COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CONTORNO
* NAS DIREÇÕES 'X','Y' E 'Z'
*
        WRITE(IPR,300)
 300    FORMAT(//2X,'COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CO
     $NTORNO'//4X,'PONTO',10X,'X',18X,'Y',18X,'Z')
        READ(INP,*)(CX(I),CY(I),CZ(I),I=1,NP)
        DO 10 I=1,NP
 10         WRITE(IPR,400)I,CX(I),CY(I),CZ(I)
 400    FORMAT(4X,I4,3(5X,D14.7))
*
* LEITURA DA CONECTIVIDADE ENTRE OS ELEMENTOS
*
        WRITE(IPR,500)
 500    FORMAT(//2X,'CONECTIVIDADE DOS ELEMENTOS DE CONTORNO'//
     $  1X,'ELEMENTO',3X,'EXTR I',3X,'EXTR J',3X,'EXTR K',
     $  3X,'EXTR L',5X,'XMED',8X,'YMED',8X,'ZMED')
        DO 18 I=1,N
            READ(INP,*) (cone_in(J),J=1,4)
            CONE(I, 1:4) = INT(cone_in(1:4))
*
* CALCULA AS COORDENADAS NODAIS E ARMAZENA NOS VETORES XM, YM E ZM
*
            N1=CONE(I,1)
            N2=CONE(I,2)
            N3=CONE(I,3)
            N4=CONE(I,4)
*
* ELEMENTOS QUADRILATERAIS DE QUATRO NÓS
*
            CXM(I)=(CX(N1)+CX(N2)+CX(N3)+CX(N4))*0.250
            CYM(I)=(CY(N1)+CY(N2)+CY(N3)+CY(N4))*0.250
            CZM(I)=(CZ(N1)+CZ(N2)+CZ(N3)+CZ(N4))*0.250
            WRITE(IPR,550) I,(CONE(I,J),J=1,4),CXM(I),CYM(I),CZM(I)
 550        FORMAT(I6,4I9,2X,3D12.4)
  18    CONTINUE
*
        RETURN
      END
