# Análises associando variantes genéticas e o desenvolvimento de mucosite oral

## Introdução

Esta consultoria foi prestada à uma pesquisadora de doutorado, cujo intuito era avaliar a relação entre variantes genéticas, protocolos de medicações utilizados e o desfecho clínico (mucosite bucal.

## Análises estatísticas

As análises englobadas neste projeto são divididas em duas etapas, sendo a primeira etapa composta dos seguintes passos:

1 - Testes qui-quadrado, para a seleção dos modelos de variantes genéticas. Cada variante genética pode possuir um ou mais modelo de representação, sendo ele aditivo, dominante ou recessivo, então para cada variante utilizamos os testes qui-quadrado, utilizando simulação de Monte Carlo para aumentar o tamanho da amostra, e teste exato de fisher para a seleção de variantes, caso a variante possua mais de um modelo que tenha um valor-p abaixo do nosso ponto de corte (0.05), foi dada a decisão à cliente.

2 - O processo anterior foi repetido para cada agrupamento de classificação da variável dependente (no caso, a classificação da mucosite bucal, variável PIORMB), possuíamos três agrupamentos de classificação, sendo eles: presença ou ausência de mucosite bucal, definido por PIORMB = 0 (ausência) ou 1,2,3,4 (presença); presença de ulcerações, definido por PIORMB = 0, 1 (sem ulcerações) ou 2, 3, 4 (com ulcerações); e por fim, severidade, definido por PIORMB = 0, 1, 2 (baixa severidade) ou 3, 4 (alta severidade).

3 - Executamos uma análise binomial univariada, para cada variante selecionada anteriormente (uma combinação de seu p-valor + seleção da cliente), também agrupadas por classificação da variável dependente, presença ou ausência, com ou sem ulceração e com severidade alta ou baixa.

4 - Uma análise binomial multivariada, filtrada pelo valor-p significativo (definimos como p < 0.05), do resultado da binomial univariada anterior. Foi utilizado o processo de step para a seleção do modelo, e a partir dele, removemos as variantes que possuíam um valor-p abaixo de 0.05, até termos apenas as variantes significativas, o processo usado foi de backward step regression.

E então, foi solicitada uma segunda análise, utilizando a mesma base de dados, porém agora agrupando por tratamento utilizado (a variável Protocolo), que indica quais medicações foram usadas nos pacientes.

O processo seguiu em mesma ordem e passando pelos mesmos passos que o anterior, apenas adicionando mais uma divisão, no caso agora do protocolo utilizado para o tratamento.

## Conclusão 

Ao final, entregamos então o resultado, indicando quais variantes genéticas que, segundo os padrões definidos, tiveram um impacto no desfecho de apresentar ou não mucosite bucal em pacientes que passaram por diversos protocolos de tratamentos médicos. 

Este projeto de consultoria foi desenvolvido pelos membros Mikael Marin Coletto e Matheus Borges
