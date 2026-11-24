# Perímetro da elipse via estatística: distribuição qui com graus de liberdade efetivos

## 1. A ponte com a probabilidade (identidade exata)

Sejam X, Y variáveis normais padrão independentes. Então

    P(a, b) = √(8π) · E[ √(a²X² + b²Y²) ].

Ou seja, o perímetro é proporcional ao **comprimento médio de um vetor gaussiano** com desvios-padrão a e b. Em telecomunicações, essa variável aleatória é conhecida como distribuição de Hoyt.

Demonstração: em coordenadas polares, X = R cos θ e Y = R sen θ, com R e θ independentes e E[R] = √(π/2). Logo E‖(aX, bY)‖ = √(π/2) · (1/2π) ∫ √(a²cos²θ + b²sen²θ) dθ = P/√(8π).

## 2. A ideia estatística

Q = a²X² + b²Y² é uma soma ponderada de qui-quadrados. Na estatística, a aproximação clássica de **Satterthwaite** troca Q por c·χ²_ν, um qui-quadrado escalado com ν "graus de liberdade efetivos". A média de √(c·χ²_ν) tem forma fechada, e o resultado é:

    P ≈ 4·√( π(a²+b²)/ν ) · Γ((ν+1)/2) / Γ(ν/2)

Esta estrutura é **exata automaticamente** nos dois extremos:
- círculo (a = b): Q é exatamente um χ² com 2 graus de liberdade, e ν = 2 dá 2πa;
- segmento (b = 0): Q é exatamente um χ² com 1 grau de liberdade, e ν = 1 dá 4a.

O ν de Satterthwaite, (a²+b²)²/(a⁴+b⁴), erra cerca de 2% no meio do intervalo. A contribuição nova é uma lei para os graus de liberdade efetivos.

## 3. A fórmula

    w = 2ab/(a²+b²)                    (vale 1 no círculo e 0 no segmento)
    W = w^p,   u = 1 − W
    ν = 1 + W / (1 + q·u + r·u²)

    P ≈ 4·√( π(a²+b²)/ν ) · Γ((ν+1)/2) / Γ(ν/2)

| constante | valor | origem |
|---|---|---|
| p | 1.613 | ajuste minimax |
| r | −0.107 | ajuste minimax |
| q | α/p − 1 = 0.3630872 | exato, com α = 1/(6 − 8 ln 2) |

O valor α = 1/(6 − 8 ln 2) = 2.19866 é a inclinação exata de ν perto do círculo (dν/dw em w = 1). Ele foi deduzido analiticamente: vem de ψ(3/2) − ψ(1) = 2 − 2 ln 2.

## 4. Precisão (comparada com a integral elíptica, 30 dígitos)

Erro relativo máximo: **1,08·10⁻⁴**. É cerca de 4× melhor que Ramanujan II (4,0·10⁻⁴) no pior caso.

| b/a | nova | Ramanujan II |
|---|---|---|
| 1 | 0 | 0 |
| 0.9 | 1.3e-7 | 3.8e-18 |
| 0.7 | 1.3e-5 | 7.0e-13 |
| 0.5 | 8.3e-5 | 4.6e-10 |
| 0.3 | 5.1e-5 | 7.4e-8 |
| 0.1 | 3.8e-5 | 1.2e-5 |
| 0.01 | 5.5e-5 | 2.4e-4 |
| 0.001 | 3.4e-6 | 3.8e-4 |
| 0 | 0 | 4.0e-4 |

Limitação: para elipses moderadas, Ramanujan II é bem mais preciso. A novidade aqui é a estrutura da fórmula, não a precisão recorde.

## 5. Literatura

Buscas (nov/2026) não encontraram uso de Satterthwaite, Nakagami ou da distribuição qui para aproximar o perímetro da elipse. A busca não foi exaustiva.

Verificação: `python3 verificacao_qui.py`.
