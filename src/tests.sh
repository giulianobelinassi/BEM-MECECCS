#!/bin/bash
# Este arquivo especifica os testes automatizados.

# Especifica os caracteres de cores.
GREEN='\033[1;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color


function print_ok {
	echo -e "[${GREEN}OK${NC}]\n"
}

function print_failed {
	echo -e "[${RED}FALHOU${NC}]\n"
}

if [ ! -f dtester ]; then
    gcc -g -O3 -flto -std=c99 -o dtester dtester.c -lm
fi

#recompila tudo com as macros que habilitam testes que rodam internamente no programa.
make clean
make gpu EXTRA="-DTEST_CUDA" -j 4

#remove os arquivos antigos.
rm -f SSOLO240E_-5+5.DAT SSOLO240D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT

# Teste 1: O arquivo SSOLO240E_-5+5.DAT, que é a solução encontrada
# por este programa, deve ser igual ao SSOLO240E_-5+5.SOL, que é o
# arquivo calculado pelo programa compilado no compilador original
# (Compaq Visual Fortran 2003, no Windows) convertido para o formato
# Unix (sem o carriage return).

echo -e "Dado a execução do programa 'main' com 240 pontos"

 ./main ESOLO240E_-5+5.DAT ESOLO240D_-5+5.DAT SSOLO240E_-5+5.DAT SSOLO240D_-5+5.DAT

sleep 0.2s

echo -e "O arquivo gerado SSOLO240E_-5+5.DAT é igual ao SSOLO240E_-5+5.SOL ?"
if cmp -s "SSOLO240E_-5+5.DAT" "SSOLO240E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO240D_-5+5.DAT é compatível em precião ao SSOLO240D_-5+5.SOL ?"
./dtester "SSOLO240D_-5+5.DAT" "SSOLO240D_-5+5.SOL"
if [ $? -eq 0 ]; then
	print_ok
else 
	print_failed
fi

echo -e "Dado a execução do programa 'main' com 2160 pontos"

 ./main ESOLO2160E_-5+5.DAT ESOLO2160D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT

echo -e "O arquivo gerado SSOLO2160E_-5+5.DAT é igual ao SSOLO2160E_-5+5.SOL ?"
if cmp -s "SSOLO2160E_-5+5.DAT" "SSOLO2160E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO2160D_-5+5.DAT é compatível em precisão ao SSOLO2160D_-5+5.SOL ?"
./dtester "SSOLO2160D_-5+5.DAT" "SSOLO2160D_-5+5.SOL"
if [ $? -eq 0 ]; then
	print_ok
else 
	print_failed
fi

echo -e "Teste com precisao dupla"
make clean
FPREC=double make gpu EXTRA="-DTEST_CUDA" -j 4

#remove os arquivos antigos.
rm -f SSOLO240E_-5+5.DAT SSOLO240D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT

echo -e "Dado a execução do programa 'main' com 240 pontos e precisao dupla"

 ./main ESOLO240E_-5+5.DAT ESOLO240D_-5+5.DAT SSOLO240E_-5+5.DAT SSOLO240D_-5+5.DAT

sleep 0.2s

echo -e "O arquivo gerado SSOLO240E_-5+5.DAT é igual ao SSOLO240E_-5+5.SOL ?"
if cmp -s "SSOLO240E_-5+5.DAT" "SSOLO240E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO240D_-5+5.DAT é compatível em precião ao SSOLO240D_-5+5.SOL ?"
./dtester "SSOLO240D_-5+5.DAT" "SSOLO240D_-5+5.SOL"
if [ $? -eq 0 ]; then
	print_ok
else 
	print_failed
fi


echo -e "Dado a execução do programa 'main' com 2160 pontos e precisao dupla"

 ./main ESOLO2160E_-5+5.DAT ESOLO2160D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT

echo -e "O arquivo gerado SSOLO2160E_-5+5.DAT é igual ao SSOLO2160E_-5+5.SOL ?"
if cmp -s "SSOLO2160E_-5+5.DAT" "SSOLO2160E_-5+5.SOL"; then
	print_ok
else 
	print_failed
fi

echo -e "O arquivo gerado SSOLO2160D_-5+5.DAT é compatível em precisão ao SSOLO2160D_-5+5.SOL ?"
./dtester "SSOLO2160D_-5+5.DAT" "SSOLO2160D_-5+5.SOL"
if [ $? -eq 0 ]; then
	print_ok
else 
	print_failed
fi
echo -e "Este é o fim dos testes.\n\n"
