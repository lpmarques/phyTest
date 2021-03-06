# phyTest (versão 0.7 - alpha)

## Descrição

Este programa é capaz de realizar uma gama de testes topológicos, visando a seleção do modelo evolutivo (i.e. árvore filogenética) mais adequado para explicar o processo de geração de um alinhamento de sequências dado pelo usuário. 

O phyTest (v. 0.7.6) implementa os testes KH, SH, SOWH, ELW, BP e AU e possibilita sua realização por uma variedade de procedimentos, como o bootstrap paramétrico, não-paramétrico, RELL ou aproximação normal. O phyTest não é autônomo no cálculo de verossimilhança das árvores testadas ou na simulação de alinhamentos. Na necessidade destes procedimentos, recruta os programas IQTREE e Seq-Gen automaticamente.

Dada a variedade de combinações possíveis entre testes e procedimentos implementados, com o phyTest o usuário terá mais liberdade para avaliar o efeito de diferentes abordagens estatísticas sobre o cálculo da confiança (p-valor) de suas hipóteses evolutivas.


## Configurações de uso

### Tutorial rápido

Para rodar o programa chame:
```
perl phyTest.pl -s <arquivo_com_alinhamento> -m <modelo_de_substituição> -z <arquivo_com_filogenias> -n <número_de_réplicas> -t <testes_desejados>
```

Os formatos aceitos, sintaxe e parâmetros disponíveis seguem:
```
-s arquivo_com_alinhamento.[fasta|phy]
-m modelo_de_substituição [MODELO+I+(Gn|+Rn)]
-z arquivo_com_filogenias.tre [uma newick por linha]
-n número_de_réplicas [dígitos]
-t [BP/KH/SH/ELW/AU/SOWH]
```

### Detalhamento

`-s` (obrigatório; sem default)

Especifica o arquivo contendo o alinhamento que se deseja analisar.
Os formatos aceitos são fasta (.fas|.fasta|.fst) e phylip sequencial (.phy|.phylip).



`-m` (obrigatório; sem default)

Especifica o modelo de substituição de nucleotídeos a ser utilizado nas análises (otimização de parâmetros e cálculos de verossimilhança) tanto do alinhamento original (dado via -s) quanto de eventuais réplicas que venham a ser produzidas. 
Segue o mesmo formato aceito no programa IQTREE; matriz de substituição e demais parâmetros separados por '+'.

ex: `HKY+G6+I`


* [Matriz de substituição] - delimita a liberdade de variação das taxas de substituição de nucleotídeos e suas frequências de equilíbrio:

  * JC/JC69, F81, K2P/K80, HKY/HKY85, TN/TrN/TN93, TNe, K3P/K81, K81u, TPM2, TPM2u, TPM3, TPM3u, TIM, TIMe, TIM2, TIM2e, TIM3, TIM3e, TVM, TVMe, SYM, GTR ou codificação em 6 dígitos (ver http://www.iqtree.org/doc/Substitution-Models#dna-models para mais detalhes).


* [Taxas heterogêneas] - admite estimação de taxas de substituição heterogêneas ao longo dos sítios do alinhamento
  * G[n] - variação restrita a categorias de uma distribuição gama, com parâmetro alfa estimado por máxima verossimilhança (ou fixo em caso de otimização parcial) e beta = alfa
  * R[n] - variação restrita a categorias de uma distribuição livre, com média e probabilidade de cada categoria estimada por máxima verossimilhança (ou fixa em caso de otimização parcial)

Em ambos os casos, [n] é opcional e define o número de categorias da distribuição. Se 'G' ou 'R' for referido sem qualquer dígito, o número de categorias seguirá um default de 4.


* [Sítios invariáveis] - admite uma proporção de sítios invariáveis (com taxa de substituição de nucleotídeos = 0); é estimada por máxima verossimilhança (ou fixa em caso de otimização parcial)
  * I - categoria única de sítios invariáveis

**Atenção**: apesar da combinação de categorias de taxas heterogêneas e uma categoria de sítios invariáveis (e.g. GTR+G+I) ser viável, o IQTREE sofre significativa penalidade no tempo de estimação da proporção de sítios invariáveis e da distribuição de taxas simultaneamente. Para uma análise mais rápida, é recomendável aumentar o número de categorias de taxas heterogêneas, compensando assim a ausência de uma categoria de invariáveis.


`-z` (obrigatório; sem default):

Especifica o arquivo contendo o grupo de árvores filogenéticas a ser testado.
Neste aquivo, as árvores devem estar em formato newick (i.e. parentético) e cada uma deve ocupar exclusivamente uma linha.


`-n` (default `1000`)

Especifica o limite de réplicas (do alinhamento em -s) a serem produzidas. Utilizado apenas no caso de um ou mais testes requisitados envolverem procedimento de bootstrap paramétrico (simulações de alinhamentos), não-paramétrico (reamostras de sítios do alinhamento original) ou RELL (reamostras dos valores de verossimilhança por sítio calculados para o alinhamento original).


`-t` (default `BP:0/KH:0/SH:0`)

Especifica os tipos de testes desejados, separados por barra (/), e os parâmetros procedimentais, combinados aos respectivos testes por dois pontos (:). 
Se desejar repetir um mesmo tipo de teste via diferentes procedimentos numa só execução do programa, os parâmetros procedimentais vão separados entre si por vírgulas (,).

ex: `-t AU:1,2,3,4/BP:3,4/KH:-3,0` 

Tipos de testes:
* [BP] - Bootstrap Proportion
* [KH] - Kishino & Hasegawa (1989)
* [SH] - Shimodaira & Hasegawa (1999)
* [SOWH] - Teste Paramétrico (Swofford, Olsen, Waddel & Hillis, 1996)
* [ELW] - Expected Likelihood Weights (Strimmer & Rambaut, 2001) 
* [AU] - Approximatelly Unbiased (Shimodaira, 2002)
 
Procedimentos:
* [-3] - aproximação normal
* [-2] - bootstrap paramétrico (simulações) com otimização completa
* [-1] - bootstrap paramétrico com otimização parcial (default para o SOWH)
* [0] - RELL (default para os demais testes)
* [1] - bootstrap não-paramétrico (tradicional) com otimização parcial
* [2] - boot. não-paramétrico com otimização completa
* [3] - boot. não-paramétrico com busca pela árvore de ML e otimização parcial
* [4] - boot. não-paramétrico com busca pela árvore de ML e otimização completa
 
  * `-3`: válido apenas com `KH`
    * Assume que os valores esperados de diferença de verossimilhança (deltas) entre cada árvore e a árvore de maior verossimilhança têm distribuição normal com média zero.
    * Assume que o desvio padrão (DP) dessa distribuição é proporcional ao DP dos valores de delta por sítio calculados entre as árvores comparadas. Assim, dispensa procedimentos de replicação.
 
  * `-2`: válido apenas com `SOWH`
    * Utiliza cada uma das árvores, junto a seus comprimentos de ramo e parâmetros de substituição, para simular diferentes sets de alinhamentos-réplica.
    * Ao recalcular a verossimilhança de cada árvore, para as diferentes réplicas simuladas, reestima os parâmetros de substituição.
 
  * `-1`: válido apenas com `SOWH` (default)
    * ======//======
    * Ao recalcular a verossimilhança de cada árvore, para as diferentes réplicas simuladas, fixa os parâmetros de substituição já estimados para o alinhamento dado.
 
  * `0`: válido com `KH`, `SH`, `ELW`, `BP` e `AU` (default)
    * Sorteia diretamente dos valores de verossimilhança por sítio calculados para cada árvore com o alinhamento dado. 
    * Dessa forma, obtém verossimilhanças-réplica (soma de cada set de valores sorteado), dispensando alinhamentos-réplica e reestimação de parâmetros.
 
  * `1`: válido com `KH`, `SH`, `ELW`, `BP` e `AU`
    * Sorteia sítios do alinhamento original, com reposição, para formar alinhamentos-réplica.
    * Ao recalcular a verossimilhança de cada árvore, para as diferentes réplicas sorteadas, fixa os parâmetros de substituição já estimados para o alinhamento original.
 
  * `2`: válido com `KH`, `SH`, `ELW`, `BP` e `AU`
    * ======//======
    * Ao recalcular a verossimilhança de cada árvore, com as diferentes réplicas sorteadas, reestima os parâmetros de substituição.
 
  * `3`: válido com `BP` e `AU`
    * ======//======
    * Utiliza o algorítmo default do IQTREE para buscar a árvore de maxima verossimilhança (ML) para cada réplica sorteada, fixando os parâmetros de substituição já estimados para o alinhamento original.
 
  * `4`: válido com `BP` e `AU`
    * ======//======
    * Utiliza o algorítmo default do IQTREE para buscar a árvore e os parâmetros de substituição de máxima verossimilhança (ML) para cada réplica sorteada.


 
`-g` 

Define o seed (número inteiro) que alimentará o gerador aleatório de números do phyTest e o do Seq-Gen (se recrutado). Por default, o seed é definido com base na hora de partida do programa, mas utilizando a opção `-g` para forçar o mesmo seed de uma corrida anterior, o usuário pode garantir que as mesmas réplicas serão geradas.



`-nc` (default `2`)

Especifica o número de núcleos de processamento a ser utilizado sempre que o IQTREE é recrutado.



`-redo` (default `TRUE`)

Define se todas as análises para otimização de parâmetros e cálculo de verossimilhança devem ser repetidas ou se resultados anteriores devem ser reaproveitados. Com `-redo FALSE`, arquivos de saída do IQTREE gerado pelo phyTest anteriormente podem ser reaproveitados.

**Atenção**: se o grupo de árvores dado em `-t` for modificado entre uma e outra corrida do phyTest, garanta que `-redo TRUE` ou que o diretório de trabalho não contenha arquivos de saída relativos a análises anteriores.



`-verb` (default `FALSE`)

Modo "verbose". Com `-verb TRUE`, o programa se torna falador, entrando em detalhes mais específicos de cada passo percorrido.

