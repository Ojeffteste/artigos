# Aproximação do perímetro da elipse pela função Gama

Notação: semieixos `a ≥ b`, `h = ((a−b)/(a+b))²`, `P = π(a+b)·f(h)`.

## Fórmula

    P ≈ π(a+b) · (2/√π) · Γ(s + 3/2) / Γ(s + 1)
      = π(a+b) · (2/√π) · (s + ½)! / s!

    s = h·(c0 + c1·h + c2·h² + c3·h³ + c4·h⁴ + c5·h⁵ + c6·h⁶)

| constante | valor | origem |
|---|---|---|
| c0 | 0.4073614192 = 1/(8(1−ln 2)) | série exata (termo h) |
| c1 | 0.0705481897  | série exata (termo h²) |
| c2 | 0.01282088551 | série exata (termo h³) |
| c3 | 0.004100342805 | série exata (termo h⁴) |
| c4 | 0.012317 | ajuste minimax |
| c5 | −0.026923 | ajuste minimax |
| c6 | 0.019775163 = ½ − (c0+…+c5) | garante s(1) = ½ |

## De onde vem

O teorema de Gauss para a função hipergeométrica diz que

    Γ(2) / (Γ(3/2)·Γ(3/2)) = 4/π,

que é exatamente o valor de `f` para a elipse degenerada (b = 0). A função `Q(s) = (2/√π)·Γ(s+3/2)/Γ(s+1)` vale 1 em s = 0 (o círculo) e 4/π em s = ½. Ela cresce de forma monótona entre esses dois valores.

Assim, toda a aproximação se resume a escolher uma "variável interna" `s(h)` com s(0) = 0 e s(1) = ½. Essa função é suave e quase linear: `s(h)/h` vai de 0,407 a 0,5. Os quatro primeiros coeficientes reproduzem a série exata, e os outros foram otimizados.

## Propriedades

- Exata para o círculo (h = 0) e para a elipse degenerada (b = 0).
- A série de Taylor coincide com a exata até o termo h⁴, como em Ramanujan II.
- Erro relativo máximo: **1,84·10⁻⁵**, cerca de 22× melhor que Ramanujan II (4,0·10⁻⁴) no pior caso.
- Fica no mesmo nível da fórmula de Cantrell (1,45·10⁻⁵).

| b/a | nova | Ramanujan II | Cantrell |
|---|---|---|---|
| 0.7 | 1.7e-10 | 7.0e-13 | 7.0e-13 |
| 0.5 | 7.5e-8 | 4.6e-10 | 4.6e-10 |
| 0.3 | 4.3e-6 | 7.4e-8 | 7.4e-8 |
| 0.2 | 1.5e-5 | 8.5e-7 | 8.2e-7 |
| 0.1 | 5.8e-6 | 1.2e-5 | 8.0e-6 |
| 0.05 | 1.8e-5 | 5.2e-5 | 1.4e-5 |
| 0.01 | 1.7e-5 | 2.4e-4 | 1.2e-5 |
| 0 | 0 | 4.0e-4 | 0 |

Limitação: para elipses moderadas (b/a > 0.2), Ramanujan II é mais preciso. A vantagem está nas elipses alongadas.

## Revisão da literatura

Buscas por aproximações do perímetro com função Gama ou razão de fatoriais não encontraram nada desse tipo (nov/2026). Nenhuma das famílias conhecidas tem essa estrutura:
- racionais e Padé;
- com raízes (Ramanujan);
- exponenciais (arXiv:2510.16191);
- logarítmicas (ELIP, arXiv:2610.28689; família Underwood);
- médias de potências (Muir e outras).

A busca não foi exaustiva: os textos completos dos artigos do arXiv não puderam ser lidos.

Verificação: `python3 verificacao_gama.py` (mpmath, 30–60 dígitos).
