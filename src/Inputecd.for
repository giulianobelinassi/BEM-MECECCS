************************************************************************
*                                                                      *
*      ESTA SUBROTINA LÊ E IMPRIME OS DADOS DO ARQUIVO DE ENTRADA      *
*                                                                      *
************************************************************************
*
      SUBROUTINE INPUTECD (CX,CY,CZ,CXI,CYI,CZI,KODE,BC,NFR,AFR,NE,NX,
     $      NCOX,NPIX,NFRX,CONE,CXM,CYM,CZM,N,NBE,NP,NPG,GE,RNU,RMU,
     $      L,FR,DAM,RHO,ZGE,ZCS,ZCP)
*
        IMPLICIT REAL (A-H,O-Y)
        IMPLICIT COMPLEX (Z)
        COMMON INP,INQ,IPR,IPS,IPT
        CHARACTER*75 TITULO
        REAL, DIMENSION(NP), INTENT(IN) :: CX, CY, CZ
        REAL, DIMENSION(N) , INTENT(IN) :: CXM, CYM, CZM
        REAL, DIMENSION(:) , INTENT(OUT), ALLOCATABLE:: CXI, CYI, CZI 
        REAL, DIMENSION(:) , INTENT(OUT), ALLOCATABLE:: BC, AFR
        INTEGER, DIMENSION(:), INTENT(OUT), ALLOCATABLE:: KODE
        INTEGER, DIMENSION(N, 4), INTENT(IN) :: CONE
        INTEGER, DIMENSION(6) :: stats
*
        WRITE(IPS,100)
 100    FORMAT(/' ',79('*')/)
*
* LEITURA DO TÍTULO DO TRABALHO
*
        READ(INQ,'(A)')TITULO
        WRITE(IPS,'(A)')TITULO
*
* LEITURA DO NÚMERO DE FREQÜÊNCIAS E DAS PROPRIEDADES DO MATERIAL
*
        READ(INQ,*)L,NFR,DAM,RHO 
        WRITE(IPS,200)N,NBE,NP,L,NPG,NFR,GE,RNU,DAM,RHO     
 200    FORMAT(//'DADOS'//2X,'NÚMERO DE ELEMENTOS DA MALHA= ',I4/
     $      2X,'NÚMERO DE ELEMENTOS DE CONTORNO= 'I4/
     $      2X,'NÚMERO DE PONTOS EXTREMOS DOS ELEMENTOS= 'I4/
     $      2X,'NÚMERO DE PONTOS INTERNOS= ',I3/
     $      2X,'NÚMERO DE PONTOS DE GAUSS= ',I3/    
     $      2X,'NÚMERO DE FREQÜÊNCIAS= ',I3/    
     $      2X,'MÓDULO DE CISALHAMENTO (REAL)= 'D14.7/
     $      2X,'COEFICIENTE DE POISSON= ',D14.7/
     $      2X,'COEFICIENTE DE AMORTECIMENTO= ',D14.7/
     $      2X,'DENSIDADE DE MASSA= ',D14.7)
*
        ALLOCATE(CXI(L), STAT = stats(1))
        ALLOCATE(CYI(L), STAT = stats(2))
        ALLOCATE(CZI(L), STAT = stats(3))
        ALLOCATE(BC(3*NBE), STAT = stats(4))
        ALLOCATE(KODE(3*NBE), STAT = stats(5))
        ALLOCATE(AFR(NFR), STAT = stats(6))
        IF(stats(1) /= 0 .or. stats(2) /= 0 .or. stats(3) /= 0 .or.
     $    stats(4) /= 0 .or. stats(5) /= 0 .or. stats(6) /= 0) THEN
            PRINT*, "MEMORIA INSUFICIENTE!"
            STOP
        ENDIF
         
        ZGE=CMPLX(GE,GE*2.0*DAM)
        ZCS=CSQRT(ZGE/RHO)
        ZCP=ZCS*SQRT((2.0-2.0*RNU)/(1.0-2.0*RNU))
        WRITE(IPS,230)ZGE,ZCS,ZCP     
 230        FORMAT(//'DADOS DERIVADOS'//
     $          2X,'MÓDULO DE CISALHAMENTO (COMPLEXO)= '2(D14.7)/
     $          2X,'VELOCIDADE DA ONDA-S (COMPLEXO)= ',2(D14.7)/
     $          2X,'VELOCIDADE DA ONDA-P (COMPLEXO)= ',2(D14.7))

*
* LEITURA DAS FREQÜÊNCIAS
*
        READ(INQ,*) (AFR(I),I=1,NFR)
        WRITE(IPS,250) (AFR(I),I=1,NFR)
 250    FORMAT(/2X,'FREQUENCIAS'//5(2X,D14.7))
*
* IMPRIME AS COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DE CONTORNO
* NAS DIREÇÕES 'X','Y' E 'Z'
*
        WRITE(IPS,300)
 300      FORMAT(//2X,'COORDENADAS DOS PONTOS EXTREMOS DOS ELEMENTOS DA 
     $MALHA'//4X,'PONTO',10X,'X',18X,'Y',18X,'Z')
C      READ(INP,*)(CX(I),CY(I),CZ(I),I=1,NP)
        DO 10 I=1,NP
 10         WRITE(IPS,400)I,CX(I),CY(I),CZ(I)
 400            FORMAT(4X,I4,3(5X,D14.7))
*
* IMPRIME A CONECTIVIDADE ENTRE OS ELEMENTOS
*
        WRITE(IPS,500)
 500        FORMAT(//2X,'CONECTIVIDADE DOS ELEMENTOS DA MALHA'//
     $      1X,'ELEMENTO',3X,'EXTR I',3X,'EXTR J',3X,'EXTR K'
     $      ,3X,'EXTR L',5X,'XMED',8X,'YMED',8X,'ZMED')
        DO 18 I=1,N
C       READ(INP,*) (CONE(I,J),J=1,4)
*
* CALCULA AS COORDENADAS NODAIS E ARMAZENA NOS VETORES XM, YM E ZM
*
C           N1=CONE(I,1)
C           N2=CONE(I,2)
C           N3=CONE(I,3)
C           N4=CONE(I,4)
*
* ELEMENTOS QUADRILATERAIS DE QUATRO NÓS
*
C           CXM(I)=(CX(N1)+CX(N2)+CX(N3)+CX(N4))*0.250
C           CYM(I)=(CY(N1)+CY(N2)+CY(N3)+CY(N4))*0.250
C           CZM(I)=(CZ(N1)+CZ(N2)+CZ(N3)+CZ(N4))*0.250
            WRITE(IPS,550) I,(CONE(I,J),J=1,4),CXM(I),CYM(I),CZM(I)
 550        FORMAT(I6,4I9,2X,3D12.4)
  18    CONTINUE

*
* LEITURA DAS CONDIÇÕES DE CONTORNO ARMAZENADAS NO VETOR BC(I)
* SE KODE(I)=0; BC(I) CONTÉM O VALOR DO DESLOCAMENTO
* SE KODE(I)=1; BC(I) CONTÉM O VALOR DA "TRACTION"
*
        WRITE(IPS,600)
 600  FORMAT(//2X,'CONDIÇÕES DE CONTORNO'//8X,3('VALOR PRESCRITO',10X)
     $/1X,'ELEM',6X,'DIREÇÃO X',4X,'CÓDIGO',6X,'DIREÇÃO Y',4X,'CÓDIGO',
     $6X,'DIREÇÃO Z',4X,'CÓDIGO')

        DO 20 I=1,NBE
            READ(INQ,*) (KODE(3*(I-1)+J),BC(3*(I-1)+J),J=1,3)
  20        WRITE(IPS,650) I,(BC(3*(I-1)+J),KODE(3*(I-1)+J),J=1,3)
 650        FORMAT(I4,3(4X,D10.3,6X,I3,2X))
*
* LEITURA DAS COORDENADAS DOS PONTOS INTERNOS
*
        IF(L == 0) THEN
            RETURN
        ENDIF

        WRITE(IPS,700)
 700    FORMAT(//2X,'COORDENADAS DOS PONTOS INTERNOS ',
     $  //4X,'PONTO',10X,'X',18X,'Y',17X,' Z')
        READ(INQ,*)(CXI(I),CYI(I),CZI(I),I=1,L) 
        WRITE(IPS,750) (I,CXI(I),CYI(I),CZI(I),I=1,L)
 750    FORMAT(4X,I3,3D19.7)
*
        RETURN
      END
