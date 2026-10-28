## C1: `&GM?WUZeEie^O`  (29 arcs, hits line index 0)

Selection rule: fewest arcs (29); first such record in file order.

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 0 1 1 1 0 0 0
  1 | 0 0 0 0 0 1 1 0
  2 | 0 0 0 1 0 1 1 0
  3 | 0 1 1 0 1 1 1 0
  4 | 0 1 1 0 0 0 0 1
  5 | 1 0 1 0 1 0 1 0
  6 | 1 0 0 1 1 0 0 1
  7 | 1 1 1 1 0 1 0 0
```

digraph6 check: `M?WUZeEie^O` is the 64-bit row-major string `0011100000000110000101100110111001100001101010101001100111110100` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][1] = 0, i.e. the arc 0 -> 1 is missing (so the pair {0, 1} is not bidirected). Total arcs 29 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 1000, eps = 1/1001, delta = 124/123123, beta = 1 + delta/eps = 247/123, q = delta(1-eps)/eps = **124000/123123** (1.007123) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 1000/1001 if A[j][i] = 1 (j sends to i), else 1 + delta = 123247/123123.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  b  b  b  b  a  a  a ]
  [ b  1  b  a  a  b  b  a ]
  [ a  b  1  a  a  a  b  a ]
  [ a  b  a  1  b  b  a  a ]
  [ a  b  b  a  1  a  a  b ]
  [ b  a  a  a  b  1  b  a ]
  [ b  a  a  a  b  a  1  b ]
  [ b  b  b  b  a  b  a  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 22701631399731769485/721844766659949312607   (0.03144946)
  x[1] = 160748905929096509024/721844766659949312607   (0.22269179)
  x[2] = 15459459716899400202/721844766659949312607   (0.02141660)
  x[3] = 46756460986593023666/721844766659949312607   (0.06477357)
  x[4] = 120728804723040119606/721844766659949312607   (0.16725037)
  x[5] = 50237199897079734539/721844766659949312607   (0.06959557)
  x[6] = 186551151001384290121/721844766659949312607   (0.25843666)
  x[7] = 118670114621507932483/721844766659949312607   (0.16439839)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 27238594/15159273129   (0.0017968272)
  c3 = 23892094076/169677744132897   (0.0001408086)
  c4 = 4461892964849/229803561799621450641   (0.0000000194)
  c5 = 379773492147163/690100096084263216274923   (0.0000000006)
  c6 = 28012325843095549/497664994191070334153444446527   (0.0000000000)
  c7 = 66903739234526727850/142972683186170023087940594709402249   (0.0000000000)
  c8 = 721844766659949312607/17603225671930811752656509842405733103627   (0.0000000000)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 2415160567060/169677744132897
       sign +, approx 1.423381e-02
  D3 = 242593627268213180074408/316695905394317485370373738699
       sign +, approx 7.660144e-07
  D4 = 448699540213604817525225470936131427/61581255210513167111354590158670443358467358508973
       sign +, approx 7.286301e-15
  D5 = 5347823655331787970262859397878935724652111059837741/22651023663428569887697392621409201607588275605419397044274161117802694714107
       sign +, approx 2.360963e-25
  D6 = 590888196932209583762603859300633681091742424308272818324525715976033/78908350919173831390148644569471003084127063349108375603766372867596123321178821238556419969733079798394723
       sign +, approx 7.488285e-39
  D7 = 817040093421459200261037108199361522126769568227287606768154076871020193259027540429913/11281738656710168329693082121183098663204813463790395615346518800266340897608030293918976166849441832856531802732468842826019925763145085932027
       sign +, approx 7.242147e-56
  D8 = 589776115587636403816826532240652193810089521356794254587342839497640987433631926249297609289360985610813191/198594991545814666480049545553477446360366999633734969160207614451337934077165621432413580305982023968874032918321458931710685782584966400278168739647190032641240460717628930669161929
       sign +, approx 2.969743e-75
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=beta_min, beta = 2, q = 1
- hi: kind=root, beta in [2, 123/61] (root of a degree-21 integer polynomial, vanishing ['Delta']), q in [1, 62/61]

- witness beta = 247/123 (2.008130) lies strictly inside the stored interval:  lo: rational endpoint 2: witness > endpoint = True;  hi: isolating interval bisected exactly 1 more step(s) (sign of the integer minimal polynomial at the midpoint) to [245/122, 123/61]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -2.409e-06 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.

## C2: `&GNx|Ub]^j}^w`  (45 arcs, hits line index 9130)

Selection rule: most arcs (45); first such record in file order.

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 0 1 1 1 1 1 1
  1 | 1 0 0 1 1 1 1 1
  2 | 0 1 0 1 0 1 1 0
  3 | 1 0 0 0 1 1 0 1
  4 | 1 1 1 0 0 1 1 1
  5 | 1 1 1 0 1 0 1 1
  6 | 1 1 1 1 1 0 0 1
  7 | 1 1 1 1 1 1 1 0
```

digraph6 check: `Nx|Ub]^j}^w` is the 64-bit row-major string `0011111110011111010101101000110111100111111010111111100111111110` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][1] = 0, i.e. the arc 0 -> 1 is missing (so the pair {0, 1} is not bidirected). Total arcs 45 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 100, eps = 1/101, delta = 48/4747, beta = 1 + delta/eps = 95/47, q = delta(1-eps)/eps = **4800/4747** (1.011165) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 100/101 if A[j][i] = 1 (j sends to i), else 1 + delta = 4795/4747.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  a  b  a  a  a  a  a ]
  [ b  1  a  b  a  a  a  a ]
  [ a  b  1  b  a  a  a  a ]
  [ a  a  a  1  b  b  a  a ]
  [ a  a  b  a  1  a  a  a ]
  [ a  a  a  a  a  1  b  a ]
  [ a  a  a  b  a  a  1  a ]
  [ a  a  b  a  a  a  a  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 8482932819557/78316252311702   (0.10831638)
  x[1] = 1040884027349/156632504623404   (0.00664539)
  x[2] = 16614838323437/78316252311702   (0.21215058)
  x[3] = 24142486120949/156632504623404   (0.15413459)
  x[4] = 8482932819557/78316252311702   (0.10831638)
  x[5] = 12713272453589/156632504623404   (0.08116625)
  x[6] = 11777863723413/52210834874468   (0.22558275)
  x[7] = 8482932819557/78316252311702   (0.10831638)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 160016/479447   (0.3337511758)
  c3 = 1131084613/106968940723   (0.0105739536)
  c4 = 138060466920/507781561612081   (0.0002718895)
  c5 = 207029743886/51285937722820181   (0.0000040368)
  c6 = 494781418142037/11442354279400687762729   (0.0000000432)
  c7 = 409628236041657/1155677782219469464035629   (0.0000000004)
  c8 = 156632504623404/116723456004166415867598529   (0.0000000000)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 284476993339/106968940723
       sign +, approx 2.659435e+00
  D3 = 123029500047628819479/11442354279400687762729
       sign +, approx 1.075211e-02
  D4 = 1558225464448896701618693568367/5810216524512759026321981606467929049
       sign +, approx 2.681872e-07
  D5 = 1085343100325958613915416915512168042668093/14005172933116317525828793610195841007285248467986879843
       sign +, approx 7.749587e-14
  D6 = 95323254560138370475675008837442656830503887785454826986/160252150444990178055592913610379751360872495437997589085972435859817016771547
       sign +, approx 5.948329e-22
  D7 = 5962506629585496728703333754034739598803452364073520022198172796405478/185199849822167015555791670827491377236655675899238349583580334237561629386870769394933851031361448063
       sign +, approx 3.219499e-32
  D8 = 933922347225627317610728841702754441948812952482022498743924065244101932427272607112/21617166522695939057853705897849071403997873899848200339739931903803535572401168103094413363330517364165524742745729351768699327
       sign +, approx 4.320281e-44
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=beta_min, beta = 2, q = 1
- hi: kind=root, beta in [2, 47/23] (root of a degree-10 integer polynomial, vanishing ['Delta']), q in [1, 24/23]

- witness beta = 95/47 (2.021277) lies strictly inside the stored interval:  lo: rational endpoint 2: witness > endpoint = True;  hi: isolating interval bisected exactly 1 more step(s) (sign of the integer minimal polynomial at the midpoint) to [93/46, 47/23]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -5.834e-05 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.

## C3: `&GMx|\Aw^qz]w`  (40 arcs, hits line index 7613)

Selection rule: largest refined q upper end (named in the task: &GMx|\Aw^qz]w).

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 0 1 1 1 0 1 1
  1 | 1 0 0 1 1 1 1 1
  2 | 0 1 0 1 1 1 0 1
  3 | 0 0 0 0 1 0 1 1
  4 | 1 0 0 0 0 1 1 1
  5 | 1 1 1 1 0 0 1 0
  6 | 1 1 1 0 1 1 0 1
  7 | 1 1 1 0 1 1 1 0
```

digraph6 check: `Mx|\Aw^qz]w` is the 64-bit row-major string `0011101110011111010111010000101110000111111100101110110111101110` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][1] = 0, i.e. the arc 0 -> 1 is missing (so the pair {0, 1} is not bidirected). Total arcs 40 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 10, eps = 1/11, delta = 3/22, beta = 1 + delta/eps = 5/2, q = delta(1-eps)/eps = **15/11** (1.363636) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 10/11 if A[j][i] = 1 (j sends to i), else 1 + delta = 25/22.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  a  b  b  a  a  a  a ]
  [ b  1  a  b  b  a  a  a ]
  [ a  b  1  b  b  a  a  a ]
  [ a  a  a  1  b  a  b  b ]
  [ a  a  a  a  1  b  a  a ]
  [ b  a  a  b  a  1  a  a ]
  [ a  a  b  a  a  a  1  a ]
  [ a  a  a  a  a  b  a  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 176418/1074793   (0.16414137)
  x[1] = 150898/1074793   (0.14039727)
  x[2] = 214698/1074793   (0.19975753)
  x[3] = 6248/1074793   (0.00581321)
  x[4] = 48488/1074793   (0.04511380)
  x[5] = 272118/1074793   (0.25318178)
  x[6] = 192038/1074793   (0.17867441)
  x[7] = 48488/1074793   (0.04511380)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 188/121   (1.5537190083)
  c3 = 277/242   (1.1446280992)
  c4 = 2495/14641   (0.1704118571)
  c5 = 246017/5153632   (0.0477366254)
  c6 = 18291/3543122   (0.0051623963)
  c7 = 1398629/2494357888   (0.0005607171)
  c8 = 1074793/27437936768   (0.0000391718)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 2731/242
       sign +, approx 1.128512e+01
  D3 = 770727/322102
       sign +, approx 2.392804e+00
  D4 = 2477808049867/26559922791424
       sign +, approx 9.329124e-02
  D5 = 43585833016098343/136880068015412051968
       sign +, approx 3.184235e-04
  D6 = 6947425277938823065061/15519448971100888972574851072
       sign +, approx 4.476593e-07
  D7 = 527395082406829899274245891/38711059958478986452554295441868455936
       sign +, approx 1.362389e-11
  D8 = 566840542805283927930664563925563/1062151615363003135741913390420777313745762254848
       sign +, approx 5.336720e-16
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=root, beta in [2, 9/4] (root of a degree-4 integer polynomial, vanishing ['a3']), q in [1, 5/4]
- hi: kind=root, beta in [8/3, 11/4] (root of a degree-15 integer polynomial, vanishing ['Delta']), q in [5/3, 7/4]

- witness beta = 5/2 (2.500000) lies strictly inside the stored interval:  lo: isolating interval bisected exactly 0 more step(s) (sign of the integer minimal polynomial at the midpoint) to [2, 9/4]; witness is above it = True;  hi: isolating interval bisected exactly 0 more step(s) (sign of the integer minimal polynomial at the midpoint) to [8/3, 11/4]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -6.847e-03 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.

## C4: `&GNBIXBs\im[O`  (31 arcs, hits line index 103)

Selection rule: narrow-interval graph: first record in file order with witness_stored_legal = null.

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 0 1 1 1 1 0 0
  1 | 0 0 1 1 0 0 1 0
  2 | 1 0 0 1 1 0 0 1
  3 | 0 0 0 0 1 1 1 1
  4 | 0 1 0 0 0 1 1 1
  5 | 0 1 1 0 1 0 1 0
  6 | 1 0 1 1 1 0 0 1
  7 | 1 1 0 0 0 1 0 0
```

digraph6 check: `NBIXBs\im[O` is the 64-bit row-major string `0011110000110010100110010000111101000111011010101011100111000100` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][1] = 0, i.e. the arc 0 -> 1 is missing (so the pair {0, 1} is not bidirected). Total arcs 31 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 10000, eps = 1/10001, delta = 2146/21452145, beta = 1 + delta/eps = 4291/2145, q = delta(1-eps)/eps = **4292000/4290429** (1.000366) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 10000/10001 if A[j][i] = 1 (j sends to i), else 1 + delta = 21454291/21452145.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  b  a  b  b  b  a  a ]
  [ b  1  b  b  a  a  b  a ]
  [ a  a  1  b  b  a  a  b ]
  [ a  a  a  1  b  b  a  b ]
  [ a  b  a  a  1  a  a  b ]
  [ a  b  b  a  a  1  b  a ]
  [ b  a  b  a  a  a  1  b ]
  [ b  b  a  a  a  b  a  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 58306433717649172625549673706155/1516670872225445607406248437543188   (0.03844370)
  x[1] = 291674459633983470524712384564555/1516670872225445607406248437543188   (0.19231230)
  x[2] = 87506735607552475880627756298075/758335436112722803703124218771594   (0.11539318)
  x[3] = 175063676780724667786381335336675/1516670872225445607406248437543188   (0.11542628)
  x[4] = 58235309164461044278790986821105/1516670872225445607406248437543188   (0.03839680)
  x[5] = 43759644962003918585263774187325/379167718056361401851562109385797   (0.11540973)
  x[6] = 174971623487526370282408603484655/1516670872225445607406248437543188   (0.11536559)
  x[7] = 204186559891076169234784343730915/758335436112722803703124218771594   (0.26925626)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 275635666019/460194525101025   (0.0005989547)
  c3 = 1000589976172073/658143978711555196575   (0.0000015203)
  c4 = 7191390526830241151/70593000310985976262152016875   (0.0000000001)
  c5 = 66174322447148408905706/908622767193789753445345846826968125   (0.0000000000)
  c6 = 29224532569233407615816033/6497302450714140296807936227093303375959375   (0.0000000000)
  c7 = 1355427641719848934704962847509/2090716114223626367962003276415377088251050399140625   (0.0000000000)
  c8 = 1516670872225445607406248437543188/44850345236161795271344248776137749526839329564672562890625   (0.0000000000)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 3152586756483113663/658143978711555196575
       sign +, approx 4.790117e-03
  D3 = 991934548926443879742059909331637/1299460490142828059361587245418660675191875
       sign +, approx 7.633434e-10
  D4 = 170207962800360606040120849692910874976676537582793284/4127976665314499265105351608860819142136423233256378025850968898830078125
       sign +, approx 4.123278e-20
  D5 = 115462989339836250132058767360897319741223814499169127238369887197515116848/1704897082067933103490422213159507485598845002986554144307283826343408613049776771960363617055928704833984375
       sign +, approx 6.772432e-35
  D6 = 9179790180524201811061155152912227707485816704251765990581849149739082444033449877164172078733967936/234326061317094334245295041288053375930902097520836583773470911422336844002655410402115584693837959349563955337630713952555473265051186757564544677734375
       sign +, approx 3.917528e-53
  D7 = 531336017038896926878480669086309821635045239087428583753422847742036233476212941080850571542667977655639663245338566445750427136/19106461622749904296761030197008370255864304320014972197568404222086478233041039774037267669446536463307626001938937140507926473525784577114869460840440654080793793953554259509241237422861158847808837890625
       sign +, approx 2.780923e-77
  D8 = 805861860407178031105056401441303243293128186323049041574194698800934199347209193525962254272923981637935022753722719347800454927038357761500513098616669447149568/856931400021809334538373297487719550016576588635952704486867971424693659427886707730659939480505415038452697672044861127557486322461451225409407455505989221758444760227999639913449807200876426748845430312155767161358490700568672312881994308554567396640777587890625
       sign +, approx 9.404042e-103
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=beta_min, beta = 2, q = 1
- hi: kind=root, beta in [2, 537/268] (root of a degree-21 integer polynomial, vanishing ['Delta']), q in [1, 269/268]

- witness beta = 4291/2145 (2.000466) lies strictly inside the stored interval:  lo: rational endpoint 2: witness > endpoint = True;  hi: isolating interval bisected exactly 3 more step(s) (sign of the integer minimal polynomial at the midpoint) to [4289/2144, 2145/1072]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -6.261e-09 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.

## C5: `&GBw~UqQEX^NG`  (34 arcs, hits line index 980)

Selection rule: 34 arcs; first such record in file order (fallback: no closed-twin double blow-up of Geneson's &EZjtwWE is in the list).

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 0 0 0 1 1 1 1
  1 | 1 0 0 0 1 1 1 1
  2 | 1 1 0 1 0 1 1 0
  3 | 1 1 0 0 1 0 0 1
  4 | 0 0 1 0 0 0 0 1
  5 | 1 0 0 1 1 0 0 1
  6 | 0 1 1 1 1 1 0 0
  7 | 1 1 1 1 0 0 1 0
```

digraph6 check: `Bw~UqQEX^NG` is the 64-bit row-major string `0000111110001111110101101100100100100001100110010111110011110010` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][1] = 0, i.e. the arc 0 -> 1 is missing (so the pair {0, 1} is not bidirected). Total arcs 34 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 100, eps = 1/101, delta = 16/1515, beta = 1 + delta/eps = 31/15, q = delta(1-eps)/eps = **320/303** (1.056106) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 100/101 if A[j][i] = 1 (j sends to i), else 1 + delta = 1531/1515.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  a  a  a  b  a  b  a ]
  [ b  1  a  a  b  b  a  a ]
  [ b  b  1  b  a  b  a  a ]
  [ b  b  a  1  b  a  a  a ]
  [ a  a  b  a  1  a  a  b ]
  [ a  a  a  b  b  1  a  b ]
  [ a  a  a  b  b  b  1  a ]
  [ a  a  b  a  a  a  b  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 3897250985925/32176727667734   (0.12112018)
  x[1] = 1864851506805/32176727667734   (0.05795653)
  x[2] = 6810023667405/32176727667734   (0.21164438)
  x[3] = 4738478854725/64353455335468   (0.07363208)
  x[4] = 13253832346335/64353455335468   (0.20595370)
  x[5] = 4217820356085/64353455335468   (0.06554147)
  x[6] = 10045484122335/64353455335468   (0.15609860)
  x[7] = 7037657662335/64353455335468   (0.10935944)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 3244/30603   (0.1060026795)
  c3 = 16231373/1159088625   (0.0140035651)
  c4 = 41706224/351203853375   (0.0001187522)
  c5 = 53838915323/7981107567946875   (0.0000067458)
  c6 = 4686011167/149276271178265625   (0.0000000314)
  c7 = 2836875097004/3663687523528173234375   (0.0000000008)
  c8 = 64353455335468/27752432990725912250390625   (0.0000000000)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 966700627/1159088625
       sign +, approx 8.340179e-01
  D3 = 16658090832689533/4030459321813171875
       sign +, approx 4.133050e-03
  D4 = 6980783474344729903079831/63698078011138881944902822265625
       sign +, approx 1.095917e-07
  D5 = 46822929818973450918748398967767346/508381212478370958595180720705415631317138671875
       sign +, approx 9.210201e-14
  D6 = 7702678561078700218611289053209789939356068094/6147029390604399134194384532379091797123405656041862010955810546875
       sign +, approx 1.253073e-21
  D7 = 558044230899479667790678127376550331907639679292440486432/22520794885118326931893679994752365861621805421419899901593291423661716282367706298828125
       sign +, approx 2.477906e-32
  D8 = 35912074488405256375968999484014309874409941001952568616656247662370176/625006850947129237293933569580234469682164414524157699008696549943795994716117546765991489519365131855010986328125
       sign +, approx 5.745869e-44
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=beta_min, beta = 2, q = 1
- hi: kind=root, beta in [2, 15/7] (root of a degree-21 integer polynomial, vanishing ['Delta']), q in [1, 8/7]

- witness beta = 31/15 (2.066667) lies strictly inside the stored interval:  lo: rational endpoint 2: witness > endpoint = True;  hi: isolating interval bisected exactly 1 more step(s) (sign of the integer minimal polynomial at the midpoint) to [29/14, 15/7]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -1.153e-04 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.

## C6: `&GVruV@zLGxLw`  (36 arcs, hits line index 4567)

Selection rule: line index 4567 (0-based) of n8_limit_HITS.jsonl.

**Adjacency matrix A** (row i = senders; A[i][j] = 1 iff arc i -> j), decoded from the digraph6 string:

```
      0 1 2 3 4 5 6 7
  0 | 0 1 0 1 1 1 1 1
  1 | 0 0 1 1 1 1 0 1
  2 | 1 0 0 1 0 1 1 1
  3 | 0 0 0 0 0 1 1 1
  4 | 1 0 1 1 0 0 1 1
  5 | 0 1 0 0 1 0 0 0
  6 | 1 1 1 0 0 1 0 0
  7 | 1 1 0 1 1 1 1 0
```

digraph6 check: `VruV@zLGxLw` is the 64-bit row-major string `0101111100111101100101110000011110110011010010001110010011011110` padded to 66 bits, 6 bits per character, char = value + 63 (`&` = header, `G` = 71 = 8 + 63).

**Not a clique:** A[0][2] = 0, i.e. the arc 0 -> 2 is missing (so the pair {0, 2} is not bidirected). Total arcs 36 < 56.

**Legal point** (witness_in_first_stable_cell of `certificates/n8_candidates_legal_witnesses.json`):  t = 10, eps = 1/11, delta = 6/55, beta = 1 + delta/eps = 11/5, q = delta(1-eps)/eps = **12/11** (1.090909) > 1.  theta = 1.

**Matrix I - W = -delta I + (1+delta) J - (eps+delta) A^T.**  Entrywise: diagonal = 1; off-diagonal (i,j) = 1 - eps = 10/11 if A[j][i] = 1 (j sends to i), else 1 + delta = 61/55.  Writing `a` = 1 - eps, `b` = 1 + delta:

```
  [ 1  b  a  b  a  b  a  a ]
  [ a  1  b  b  b  a  a  a ]
  [ b  a  1  b  a  b  a  b ]
  [ a  a  a  1  a  b  b  a ]
  [ a  a  b  b  1  a  b  a ]
  [ a  a  a  a  b  1  a  a ]
  [ a  b  a  a  a  b  1  a ]
  [ a  a  a  a  a  b  b  1 ]
```

(Check: the pattern of `a` entries is the transpose of the `1` entries of A.)

**Fixed point x = (I - W)^{-1} 1** (exact, all entries positive):

```
  x[0] = 3186480/15106621   (0.21093268)
  x[1] = 11099605/90639726   (0.12245850)
  x[2] = 97350/15106621   (0.00644419)
  x[3] = 405625/90639726   (0.00447513)
  x[4] = 13718815/45319863   (0.30271087)
  x[5] = 6638500/45319863   (0.14648103)
  x[6] = 20011255/90639726   (0.22077797)
  x[7] = 405625/90639726   (0.00447513)
```

**Characteristic polynomial of -I + W = -(I - W)**, det(s I - (-I + W)) = s^8 + c1 s^7 + ... + c8:

```
  c1 = 8   (8.0000000000)
  c2 = 148/121   (1.2231404959)
  c3 = 161186/166375   (0.9688114200)
  c4 = 1078772/9150625   (0.1178905266)
  c5 = 11649348/503284375   (0.0231466514)
  c6 = 56124776/27680640625   (0.0020275823)
  c7 = 118900781/1522435234375   (0.0000780991)
  c8 = 90639726/16746787578125   (0.0000054124)
```

**Hurwitz minors** D1..D8 of the 8x8 Hurwitz matrix H[i][j] = c_{2j-i} (1-indexed, c_0 = 1, c_k = 0 outside 0..8):

```
  D1 = 8
       sign +, approx 8.000000e+00
  D2 = 1466814/166375
       sign +, approx 8.816313e+00
  D3 = 32705335324/27680640625
       sign +, approx 1.181524e+00
  D4 = 2700062277185938/50659032423828125
       sign +, approx 5.329873e-02
  D5 = 23456050914424438865864/127479497357655364990234375
       sign +, approx 1.839986e-04
  D6 = 910074737928555685885591059501/3528714153412895250897884368896484375
       sign +, approx 2.579055e-07
  D7 = 13023496727870623941928783548024696321/5372238759193540887448045018352568149566650390625
       sign +, approx 2.424222e-12
  D8 = 1180446174976089917545464852306266315768648046/89967741319184053676649209044662575938855297863483428955078125
       sign +, approx 1.312077e-17
```

My verdict at this point: positive = True, Routh-Hurwitz = **stable**, legal = True, complete = False.  => stable positive full-support fixed point of a non-complete digraph at a legal point: **counterexample confirmed at this point**.

**Stored certified limit-stable interval** (`cert.stable_intervals[0]`, from `n8/n8_limit_HITS.jsonl`; eps -> 0 ray, q -> beta - 1):

- lo: kind=beta_min, beta = 2, q = 1
- hi: kind=root, beta in [9/4, 16/7] (root of a degree-15 integer polynomial, vanishing ['Delta']), q in [5/4, 9/7]

- witness beta = 11/5 (2.200000) lies strictly inside the stored interval:  lo: rational endpoint 2: witness > endpoint = True;  hi: isolating interval bisected exactly 0 more step(s) (sign of the integer minimal polynomial at the midpoint) to [9/4, 16/7]; witness is below it = True.

- Non-decision float sanity (numpy eigvals of -I+W at this point): spectral abscissa approx -3.282e-03 < 0 (consistent with the exact Hurwitz verdict; not used for any decision).

- Is this graph a closed-twin blow-up of one of the 123 n = 7 limit counterexamples: False.
