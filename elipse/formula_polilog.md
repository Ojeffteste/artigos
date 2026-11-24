# Perímetro da elipse via números de Catalan e polilogaritmos

    h = ((a−b)/(a+b))²
    P ≈ π(a+b)·[ 1 + h/4 + h²/64 + Σ_{k=3..6} g_k·(Li_k(h) − h − h²/2^k) ]

    g3 = 0.07953808367625804   g4 = 0.06062630887318473
    g5 = 0.03505049186615928   g6 = 0.04933973528389241

- Exata no círculo e no segmento.
- Série exata até h⁵.
- Erro relativo máximo: 9,2·10⁻⁹. Para b/a ≥ 0,2: 7,3·10⁻¹¹.
- Supera Ramanujan II em todos os pontos testados.

A dedução completa está em `pdf_polilog/formula_polilog.pdf`. Verificação: `python3 verificacao_polilog.py`.
