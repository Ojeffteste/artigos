# Fórmula de alta precisão (inclusive para elipses moderadas): razão de médias da distribuição qui

## Ideia

Seja χ_ν a distribuição qui com ν graus de liberdade, cuja média é E[χ_ν] = √2·Γ((ν+1)/2)/Γ(ν/2). Vale:

    E[χ_3] / E[χ_2] = 4/π

que é exatamente P/(π(a+b)) para a elipse degenerada (b = 0). E E[χ_2]/E[χ_2] = 1 corresponde ao círculo.

Assim, o perímetro é escrito como

    P = π(a+b) · E[χ_ν] / E[χ_2],     com ν = 2 + 2s,

e o grau de liberdade efetivo varia de 2 (círculo) a 3 (segmento). Em funções Gama:

    P ≈ π(a+b) · (2/√π) · Γ(s + 3/2) / Γ(s + 1)

## A fórmula

    h = ((a−b)/(a+b))²

    s = h · (n0 + n1·h + n2·h² + n3·h³) / (1 + d1·h + d2·h² + d3·h³)

    P ≈ π(a+b) · (2/√π) · Γ(s + 3/2) / Γ(s + 1)

| constante | valor |
|---|---|
| n0 | 1/(8(1 − ln 2)) = 0.4073614191588662 |
| n1 | −0.5017006478790753 |
| n2 | 0.1280463626575522 |
| n3 | 0.006629020704070353 |
| d1 | −1.404769353861877 |
| d2 | 0.5261406749235058 |
| d3 | −0.04069901177880171 |

Os coeficientes **não foram ajustados por otimização numérica**. Eles vêm de um aproximante de Padé [3/3] da série exata de s(h), com duas condições:
- a série de P coincide com a exata até o termo h⁶ (Ramanujan II só até h⁴);
- s(1) = ½ exatamente (segmento).

O denominador não se anula em 0 ≤ h ≤ 1: a primeira raiz está em h ≈ 1.205.

Use as constantes com 16 dígitos. Arredondá-las para 10 dígitos piora a precisão perto do círculo.

## Precisão (comparada com a integral elíptica, 40 dígitos)

| b/a | nova | Ramanujan II |
|---|---|---|
| 0.9 | 3e-22 | 3.8e-18 |
| 0.7 | 1.5e-16 | 7.0e-13 |
| 0.5 | 1.1e-12 | 4.6e-10 |
| 0.4 | 3.9e-11 | 6.3e-9 |
| 0.3 | 9.9e-10 | 7.4e-8 |
| 0.2 | 1.9e-8 | 8.5e-7 |
| 0.1 | 1.9e-7 | 1.2e-5 |
| 0.05 | 8.8e-7 | 5.2e-5 |
| 0.01 | 1.1e-5 | 2.4e-4 |
| 0 | 0 | 4.0e-4 |

- **Melhor que Ramanujan II em todos os pontos testados.** Entre a elipse moderada e o círculo, é de 40 a 4000 vezes mais precisa.
- Erro máximo para b/a ≥ 0.2: **1,9·10⁻⁸** (Ramanujan II: 8,5·10⁻⁷).
- Erro máximo geral: **1,2·10⁻⁵**, perto de b/a ≈ 0.008 (Ramanujan II: 4,0·10⁻⁴).

## Literatura

Buscas (nov/2026) não encontraram aproximações do perímetro baseadas na função Gama, em razões de fatoriais ou na distribuição qui. A busca não foi exaustiva: os textos completos do arXiv não puderam ser lidos.

Verificação: `python3 verificacao_moderada.py`.
