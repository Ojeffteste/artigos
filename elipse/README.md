# Nova aproximação para o perímetro da elipse (termo logarítmico)

Notação: semieixos `a ≥ b > 0`, `h = ((a−b)/(a+b))²`, `P = π(a+b)·f(h)`.
Valor exato: `f(h) = Σ C(1/2,n)² hⁿ = 1 + h/4 + h²/64 + h³/256 + …`, com `f(1) = 4/π` (elipse degenerada, b = 0).

## Ideia

As fórmulas clássicas (Ramanujan, Padé, Cantrell…) são racionais ou com raízes.
Mas a função exata tem uma **singularidade logarítmica** em `h = 1` (a integral elíptica E(k) contém `k'² ln k'`).
Nenhuma função racional reproduz isso, e é exatamente aí que elas erram mais.
As fórmulas abaixo incluem o termo `(1−h)·ln(1−h)`, que se anula em h = 1, mas imita o comportamento real perto do caso degenerado.

## Fórmula 1 (constantes exatas)

    f(h) ≈ 1 + h/4 + h²/64 + (4/π − 81/64)·h³ + (4/π − 325/256)·h²·(1−h)·ln(1−h)

- Exata em `h = 0` (círculo) e em `h = 1` (b = 0).
- A série de Taylor coincide com a exata até o termo h³.
- Erro relativo máximo: **3,7·10⁻⁵** (Ramanujan II: 4,0·10⁻⁴).

## Fórmula 2 (otimizada, minimax)

    f(h) ≈ 1 + h/4 + h²/64 + A·h³ + (4/π − 81/64 − A)·h⁴ + h²·(1−h)·ln(1−h)·[B + C·(1−h)]

    A = 0.000717     B = 0.001508     C = −0.005305

- Também exata em h = 0 e h = 1.
- Erro relativo máximo: **5,2·10⁻⁶**. É cerca de 80× melhor que Ramanujan II e 3× melhor que Cantrell (1,4·10⁻⁵) no pior caso.

Em h = 1, use `ln(1−h)·(1−h) = 0` (limite).

## Comparação (erro relativo)

| b/a   | Fórmula 1 | Fórmula 2 | Ramanujan II | Cantrell |
|-------|-----------|-----------|--------------|----------|
| 0.9   | 1.9e-14   | 1.3e-11   | 3.7e-18      | 3.7e-18  |
| 0.5   | 4.6e-8    | 5.6e-7    | 4.6e-10      | 4.6e-10  |
| 0.3   | 1.8e-6    | 4.2e-6    | 7.4e-8       | 7.4e-8   |
| 0.1   | 2.6e-5    | 4.3e-6    | 1.2e-5       | 8.0e-6   |
| 0.01  | 3.2e-5    | 2.0e-6    | 2.4e-4       | 1.2e-5   |
| 0.001 | 2.2e-5    | 4.8e-6    | 3.8e-4       | 5.4e-6   |
| 0     | 0         | 0         | 4.0e-4       | 0        |
| **máx** | **3.7e-5** | **5.2e-6** | **4.0e-4** | **1.4e-5** |

Ressalva: para elipses pouco achatadas (b/a ≳ 0.3), Ramanujan II continua muito mais preciso.
O ganho destas fórmulas está no **pior caso**, ou seja, em elipses muito alongadas.

Verificação numérica (mpmath, 30 dígitos, comparando com a integral elíptica completa E): `python3 verificacao.py`.
