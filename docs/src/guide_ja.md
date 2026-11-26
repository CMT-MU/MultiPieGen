---
orphan: true
---

# ユーザーガイド

英語版は[こちら](guide.md)です。

MultiPieGen が提供する量について、その定義（第 1 部）と、どう生成しているか（第 2 部）を説明します。

## 定義

この部では、MultiPieGen が提供する量を定義します。実際の関数形、$U$ 行列をデータから取り出す方法は、後の節で説明します。調和関数と SAMB には多極子のキー `PGMultipoleType` が付いていて、例ではその LaTeX 記号（例：$\mathbb{Q}_{1}(T_{1u})$）で書きます。どちらの説明もこの部の最後の節の[多極子のキー](#ja-multipole-key)にあります。

### 1. 点群調和関数

**球面調和関数。** ランク `l` の調和関数は

$$
O_{lm}(\boldsymbol{r}) = \sqrt{\frac{4\pi}{2l+1}}\, r^l\, Y_{lm}(\theta, \phi), \qquad m = l, l-1, \dots, -l
$$

から作ります。これは `l` 次の同次多項式です（例：$O_{10} = z$、$O_{1,\pm1} = \mp(x \pm iy)/\sqrt{2}$）。その実の組み合わせが方域調和関数です（$m > 0$）。

$$
C_{lm} = \frac{(-1)^m O_{lm} + O_{l,-m}}{\sqrt{2}}, \qquad
S_{lm} = \frac{(-1)^m O_{lm} - O_{l,-m}}{\sqrt{2}\, i}, \qquad C_{l0} = O_{l0}
$$

例：$C_{11} = x$、$S_{11} = y$、$C_{22} = \sqrt{3}(x^2 - y^2)/2$、$S_{22} = \sqrt{3}xy$。

**点群調和関数。** 点群の点群調和関数 `(X, l, Gamma, n)` は、既約表現 `Gamma` の基底となるランク `l` の関数の組です。成分 `gamma` は

$$
X_{l\Gamma n}^{(\gamma)}(\boldsymbol{r}) = \sum_{m} U_{m\gamma}\, O_{lm}(\boldsymbol{r})
$$

で、$U$（`(2l+1) x dim`、行は `m = l, ..., -l`）の列は正規直交です。

*例（Oh）。* 双極子 $\mathbb{Q}_{1}(T_{1u}) = (x, y, z)$ の $U$ は

$$
U = \begin{pmatrix} -\frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \\ 0 & 0 & 1 \\ \frac{1}{\sqrt{2}} & \frac{i}{\sqrt{2}} & 0 \end{pmatrix}
\quad (m = 1, 0, -1)
$$

で、例えば $x = (-O_{11} + O_{1,-1})/\sqrt{2}$ です。四極子は $\mathbb{Q}_{2}(E_{g}) = \left(z^2 - \frac{x^2+y^2}{2}, \frac{\sqrt{3}}{2}(x^2 - y^2)\right)$ です。

- **トップの群。** 調和関数は、2 つのトップの群 Oh（立方晶系）と D6h（三方晶系・六方晶系）について、ランク 11 まで一度だけ定義します。47 の点群（結晶点群 32、別の設定 `-1` が 5、複素の `-c` が 10）は、どれもそこから調和関数を受け取ります。各成分は、トップの調和関数の成分に符号を掛けたもの、$X_{\rm PG} = c\, X_{\rm top}$（$c = \pm 1$）です。
- 既約表現の表現行列は $U$ から決まり、同じ既約表現の調和関数はすべて同じ表現行列を与えます。
- 既約表現は指標表の順に並べます。複素の既約表現（`-c` の群）は名前が `a`/`b` で終わり、成分は 1 つです。
- *多重度。* 同じランクに同じ既約表現が複数回現れるときは、`n = 1, 2, ...` で区別します。例：C4h で $\mathbb{Q}_{2}(B_{g},1) = \frac{\sqrt{3}}{2}(x^2 - y^2)$、$\mathbb{Q}_{2}(B_{g},2) = \sqrt{3}xy$。
- *複素の群。* 例：C4-c で $\mathbb{Q}_{1}(E^{(a)}) = (x + iy)/\sqrt{2}$、$\mathbb{Q}_{1}(E^{(b)}) = (x - iy)/\sqrt{2}$。
- キー：`PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q")`（[多極子のキー](#ja-multipole-key)を参照）。

**極性（`Q`）と軸性（`G`）の調和関数。** 軸性の調和関数は、非固有な操作 `R`（鏡映、反転、回反）で `det R` の因子が余分に付きます。どちらも同じトップの調和関数から作りますが、点群では、同じ関数から作った極性と軸性の調和関数が別の既約表現に属することがあり、同じ既約表現でも成分の順や符号が違うことがあります（係数 `c` と割り当てる成分が `Q` と `G` で異なるため）。そのため、同じ `(l, Gamma, n)` でも、`Q` と `G` の関数形は一般に異なります。

- トップの群 Oh と D6h では、軸性の調和関数は、極性の調和関数の既約表現の `g` と `u` を入れ替えたもので、関数・成分・順序は同じです（$c = 1$）。例：Oh で $\mathbb{Q}_{1}(T_{1u}) = \mathbb{G}_{1}(T_{1g}) = (x, y, z)$。
- 回転だけからなる群（C1、C2、D2、C4、D4、C3、D3、D3-1、C6、D6、T、O）では、非固有な操作がないので、`Q` と `G` は一致します。例：C4 で $\mathbb{Q}_{1}(E) = \mathbb{G}_{1}(E) = (x, y)$。
- それ以外の群では異なります。例：

| 群 | 極性 | 軸性 | 違い |
|---|---|---|---|
| C4v | $\mathbb{Q}_{1}(A_{1}) = z$ | $\mathbb{G}_{1}(A_{2}) = z$ | 関数は同じで、既約表現が別 |
| C4v | $\mathbb{Q}_{1}(E) = (x, y)$ | $\mathbb{G}_{1}(E) = (-y, x)$ | 成分の入れ替えと符号 |
| C4h | $\mathbb{Q}_{1}(E_{u}) = (x, y)$ | $\mathbb{G}_{1}(E_{g}) = (y, x)$ | 成分の入れ替え |
| D2d | $\mathbb{Q}_{1}(E) = (x, y)$ | $\mathbb{G}_{1}(E) = (x, -y)$ | 成分の符号 |
| C4v | $\mathbb{Q}_{2}(B_{1}) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | $\mathbb{G}_{2}(B_{1}) = \sqrt{3}xy$ | 同じ既約表現で関数が別 |
| D2d | $\mathbb{Q}_{2}(A_{1}) = z^2 - \frac{x^2+y^2}{2}$ | $\mathbb{G}_{2}(A_{1}) = \frac{\sqrt{3}}{2}(x^2 - y^2)$ | 同じ既約表現で関数が別 |

(ja-multipolar-harmonics)=
### 2. 点群多極子調和関数

点群多極子調和関数は、点群調和関数と内部自由度を組み合わせたものです。

**内部自由度。** まず、内部自由度をそのランク `s`（`s = 1`：ベクトル、`s = 2, 3`：より高いランクのテンソル）とタイプ `x` で固定します。タイプは極性（`x = q`。例：変位、電気分極）か軸性（`x = g`。例：スピン、回転）です。その基底 $e_\nu$ は、そのタイプのランク `s` の軌道の組で、立方晶系（Oh）か六方晶系（D6h）のものを使い、記号 $Q_\nu$（`x = q`）または $G_\nu$（`x = g`）で書きます。例えば `s = 1` では $(Q_x, Q_y, Q_z)$ または $(G_x, G_y, G_z)$ です。

**多極子球面調和関数。** 全ランク `l` の多極子球面調和関数は、ランク `l + k` の軌道部分と内部基底を結合したものです。

$$
O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} \sqrt{2l+1} \sum_{\nu}
\begin{pmatrix} l+k & l & s \\ m-\nu & -m & \nu \end{pmatrix} O_{l+k, m-\nu}\, e_\nu,
\qquad k = -s, \dots, s
$$

**多極子のタイプ。** 内部タイプ `x` を固定すると、多極子のタイプ `X` は `s + k` の偶奇で決まります。

| 内部タイプ | `s + k` 偶 | `s + k` 奇 |
|---|---|---|
| `x = q`（極性） | `Q` | `G` |
| `x = g`（軸性） | `G` | `Q` |

つまり、`X` の極性・軸性は、`s + k` が偶なら内部基底と同じ、奇なら逆です。点群の点群多極子調和関数は、点群調和関数 `(X, l, Gamma, n)` の $U$ を使った $\sum_m U_{m\gamma} O_{lm}^{(s)}(k)$ です。

- 範囲：`s = 1, 2, 3`、`l <= 3`、47 の点群すべて。キー：`PGMultipoleType(X, l, Gamma, n, p, s, k, x)`（[多極子のキー](#ja-multipole-key)を参照）。
- 値は、内部基底の各関数の係数（x, y, z の多項式）と、内部基底の記号で書いた式です。

*例（Oh、`s = 1`）。* 内部基底を極性ベクトル $\boldsymbol{Q} = (Q_x, Q_y, Q_z)$（`x = q`）または軸性ベクトル $\boldsymbol{G} = (G_x, G_y, G_z)$（`x = g`）とすると：

| 記号 | 式 | 意味 |
|---|---|---|
| $\vec{\mathbb{Q}}_{0}^{(1,1)}[q](A_{1g})$ | $(xQ_x + yQ_y + zQ_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{Q}$、電気単極子 |
| $\vec{\mathbb{Q}}_{1}^{(1,-1)}[q](T_{1u})$ | $(Q_x, Q_y, Q_z)$ | 極性ベクトルそのもの |
| $\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$ | $(zQ_y - yQ_z, xQ_z - zQ_x, yQ_x - xQ_y)/\sqrt{2}$ | $\boldsymbol{Q} \times \boldsymbol{r}$、電気トロイダル双極子 |
| $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | $(xG_x + yG_y + zG_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{G}$、電気トロイダル単極子 |

### 3. 原子軌道

ランク `L = 0, 1, 2, 3`（s, p, d, f）の原子軌道が、原子 SAMB の基底です。基底は 3 種類あります。

**`lg`（実の軌道、スピンなし）。** 方域調和関数を $C_{LL}, S_{LL}, C_{L,L-1}, S_{L,L-1}, \dots, C_{L0}$ の順に並べたものです。

| `L` | 軌道 | 関数 |
|---|---|---|
| 0 | `s` | $1$ |
| 1 | `px`, `py`, `pz` | $x$, $y$, $z$ |
| 2 | `dv`, `dxy`, `dxz`, `dyz`, `du` | $\frac{\sqrt{3}}{2}(x^2-y^2)$, $\sqrt{3}xy$, $\sqrt{3}xz$, $\sqrt{3}yz$, $z^2 - \frac{x^2+y^2}{2}$ |
| 3 | `fc`, `fs`, `fbz`, `f3`, `fx`, `fy`, `faz` | $\frac{\sqrt{10}}{4}x(x^2-3y^2)$, $\frac{\sqrt{10}}{4}y(3x^2-y^2)$, $\frac{\sqrt{15}}{2}z(x^2-y^2)$, $\sqrt{15}xyz$, $\frac{\sqrt{6}}{4}x(4z^2-x^2-y^2)$, $\frac{\sqrt{6}}{4}y(4z^2-x^2-y^2)$, $\frac{1}{2}z(2z^2-3x^2-3y^2)$ |

**`lgs`（実の軌道、スピンあり）。** `lg` にスピン `u`（上向き）と `d`（下向き）を掛けたもので、スピンを先に回します：`(px, u), (px, d), (py, u), ...`。

**`jm`（全角運動量）。** $|J, M\rangle$（$J = L \mp 1/2$）を、`J` の昇順、`M` の降順に並べたものです。例えば `L = 1` では $|1/2, 1/2\rangle, |1/2, -1/2\rangle, |3/2, 3/2\rangle, \dots, |3/2, -3/2\rangle$ です（`L = 0` は $J = 1/2$ だけ）。データでは、この基底を `jml` と呼んでいます。

**ブロック。** これらの基底での行列は、bra のランク `L1` と ket のランク `L2` のブロックごとに与えます。ブロックは `(L1, L2) = (0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)`（`L1 <= L2`）です。
例えばブロック `(1, 2)`（p-d）は、`lg` で `3 x 5`、`lgs` と `jml` で `6 x 10` です。`lgs` での p 軌道の並びは `(px, u), (px, d), (py, u), (py, d), (pz, u), (pz, d)` です。

### 4. 原子 SAMB

原子 SAMB（対称性適合多極子基底）は、原子多極子演算子の軌道基底での行列を、点群に適合させたものです。

1. **原子多極子。** `X = Q, G, T, M` の球面多極子演算子 $X_{lm}^{(s)}(k)$ です（`Q` と `G` は時間反転で偶、`T` と `M` は奇）。`s = 0`（軌道）と `s = 1`（スピンに依存）があり、`k = -s, ..., s` です。ランク `L1` と `L2` の軌道の間の行列要素は、Wigner–Eckart の定理から求めます。
   - `s = 1` では、内部自由度を軸性ベクトルのスピンに固定します。これは[点群多極子調和関数](#ja-multipolar-harmonics)の `x = g` の場合です。$X_{lm}^{(1)}(k)$ はランク `l + k` の軌道多極子とスピンを結合したものです。点群多極子調和関数では軌道部分 $O_{l+k,m}$ はいつも極性の関数ですが、ここでの軌道多極子は極性（`Q`、`T`）のことも軸性（`G`、`M`。例えば軌道角運動量）のこともあります。そのため、そこでの `x = g` の規則を軌道多極子に対して当てはめます。つまり、`X` の極性・軸性は、`s + k` が偶なら軌道多極子と逆、奇なら同じです。さらに、スピンは時間反転で奇なので、`X` の時間反転の偶奇は軌道多極子と逆になります。
   - 例えば、スピンの磁気双極子 $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ は、軌道の単極子 `Q`（ランク 0、極性、偶）とスピンの結合です。`s + k = 0` が偶なので軸性、スピンで奇になり、`M` です。スピン軌道相互作用 $\boldsymbol{l} \cdot \boldsymbol{s}$ は、p 軌道で $\sqrt{3}\, \mathbb{Q}_{0}^{(1,1;a)}(A_{1g})$ です。軌道角運動量 `M`（ランク 1、軸性、奇）とスピンの結合で、`s + k = 2` が偶なので極性、スピンで偶になり、`Q` です。
2. **対称性への適合。** 多極子全体のタイプの点群調和関数（`Q, T` は極性、`G, M` は軸性）の $U$ を使って、$X_{l\Gamma n}^{(s)}(k)[\gamma] = \sum_m U_{m\gamma} X_{lm}^{(s)}(k)$ とします。
3. **規格化。** 各成分を、対角ブロック（`L1 = L2`）ではフロベニウスノルム 1 に、非対角ブロックでは `1/sqrt(2)` に規格化します（エルミート共役のブロック `(L2, L1)` と合わせてノルム 1）。

- 基底：`lg`（`s = 0` だけ）、`lgs`、`jml`。
- キー：`PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")`（[多極子のキー](#ja-multipole-key)を参照）。値：`(dim, n1, n2)`。ランク `L1` の `n1` 個の軌道とランク `L2` の `n2` 個の軌道の間の、各成分の行列です。
- 並び：`Q, G, T, M`、次に `s`、`k`、`l`、`Gamma`、`n`。47 の点群すべて。

*例（Oh）。* p 軌道（`lg`、ブロック `(1, 1)`、基底 `px, py, pz`）では：

- $\mathbb{Q}_{0}^{(a)}(A_{1g}) = \frac{1}{\sqrt{3}}\,\mathrm{diag}(1, 1, 1)$：電荷。
- $\mathbb{Q}_{2}^{(a)}(E_{g}) = \left(\frac{1}{\sqrt{6}}\,\mathrm{diag}(-1, -1, 2),\ \frac{1}{\sqrt{2}}\,\mathrm{diag}(1, -1, 0)\right)$：軌道の四極子。
- $\mathbb{M}_{1}^{(a)}(T_{1g})$：$(l_x, l_y, l_z)/\sqrt{2}$、軌道角運動量。例えば `x` 成分は $l_x/\sqrt{2}$ で、その行列要素は $\langle p_y | l_x | p_z \rangle / \sqrt{2} = -i/\sqrt{2}$ です。

スピンを含む s 軌道（`lgs`、ブロック `(0, 0)`、基底 `(s, u), (s, d)`）では、$\mathbb{Q}_{0}^{(a)}(A_{1g})$ は $\sigma_0/\sqrt{2}$、$\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ は $(\sigma_x, \sigma_y, \sigma_z)/\sqrt{2}$、つまりスピンです。

### 5. クラスター SAMB

Wyckoff 位置のクラスター SAMB は、その `N` 個の等価なサイト（またはボンド）上の関数で、群の既約表現の基底になります。実の点群 37 と空間群 230 について与えます。

- **`site`。** サイトの Wyckoff 位置（例：`4b`）の `N` 個のサイト上の関数です。値は実数で `X = Q` です。完全な正規直交系（`N` 個）をなします。
- **`bond`。** ボンドの Wyckoff 位置の `N` 本のボンド上の関数です。ラベル `2b@1o` は、中心がサイトの Wyckoff 位置 `1o` にある、2 本からなるボンドの Wyckoff 位置 `2b` を表します。次の 2 種類をまとめています。
  - 対称ボンド（ボンドの 2 つの向きで値が同じ。実数の飛び移り）：`X = Q`、値は実数、完全な正規直交系（`N` 個）。
  - 反対称ボンド（2 つの向きで符号が逆。虚数の飛び移り）：反対称な関数に `i` を掛けたもので、`X = T`（極性）または `M`（軸性）、値は純虚数。

  Q と T/M は、それぞれの中で正規直交です。両者は、エルミートなボンド行列として互いに直交します（時間反転の偶奇が異なるため）。
- キー：`PGMultipoleType(X, l, Gamma, n, p, 0, 0, "q")`（[多極子のキー](#ja-multipole-key)を参照）。`(X, l, Gamma, n)` はその関数に対応する点群調和関数を表し（`T` と `M` は `Q` と `G` のもの）、`p` はそれが同じ関数どうしを区別します。
- 値：`(dim, N)`。各成分の `N` 個のサイト（ボンド）上の値で、各成分のノルムは 1 です。
- 並び：`Gamma`（指標表の順）。`bond` では、各 `Gamma` の中で `Q` を `T/M` より先に並べます。

*例（D4h）。*

- サイト `2a`：$(0, 0, \pm z)$ の 2 つのサイト。$\mathbb{Q}_{0}^{(s)}(A_{1g}) = (1, 1)/\sqrt{2}$（一様）、$\mathbb{Q}_{1}^{(s)}(A_{2u}) = (1, -1)/\sqrt{2}$（`z` 方向の双極子）。
- ボンド `2b@1o`：原点を通る、`x` 方向と `y` 方向の 2 本のボンド。$\mathbb{Q}_{0}^{(b)}(A_{1g}) = (1, 1)/\sqrt{2}$、$\mathbb{Q}_{2}^{(b)}(B_{1g}) = (-1, 1)/\sqrt{2}$（実数の飛び移り）、$\mathbb{T}_{1}^{(b)}(E_{u},a) = \left((i, 0), (0, i)\right)$（虚数の飛び移り）。`x` 成分は `x` 方向のボンドに、`y` 成分は `y` 方向のボンドに値を持ちます。

### 6. 応答テンソル

ランク 0〜4 の応答テンソルは、2 つの物理量（例えば電流と電場）を結び付けます。MultiPieGen は、157 の磁気点群それぞれについて、対称性で許されるテンソルを与えます。

- タイプ `X` は、空間反転と時間反転での偶奇で決まります：`Q`（極性、時間反転で偶）、`G`（軸性、偶）、`T`（極性、奇）、`M`（軸性、奇）。
- `opt` は、添字の対称性でテンソルの部分を指定します。添字の組について `s` は対称、`a` は反対称で、関係は下の表のとおりです（極性で時間反転に偶なテンソル `X = Q` の場合。ほかのタイプも同じ関係です）。`;` の前の添字が行列の行の組（または 3 つ組）、後ろが列です。添字の集合は、`c = [x, y, z]`、`v = [xx, yy, zz, yz, zx, xy]`（対称な組）、`a = [yz, zx, xy]`（反対称な組）、`u = [xx, yy, zz, yz, zx, xy, zy, xz, yx]`（すべての組）、`t = [xxx, yyy, zzz, yyz, zzx, xxy, yzz, zxx, xyy, xyz]`（対称な 3 つ組）です。「独立な数」は、その部分の独立なパラメータの数です。

| ランク | `opt` | 成分間の関係 | 独立な数 | 行列（行 x 列） |
|---|---|---|---|---|
| 0 | （なし） | - | 1 | 1 x 1 |
| 1 | （なし） | - | 3 | 1 x 3（`c`） |
| 2 | `s` | $X_{ab} = X_{ba}$ | 6 | 3 x 3（`c` x `c`） |
| 2 | `a` | $X_{ab} = -X_{ba}$ | 3 | 3 x 3（`c` x `c`） |
| 3 | `ss` | $X_{abc}$ は $a, b, c$ について完全対称 | 10 | 6 x 3（`v` x `c`） |
| 3 | `sa` | $X_{ab;c} = X_{ba;c} = -(X_{ac;b} + X_{bc;a})$ | 8 | 6 x 3（`v` x `c`） |
| 3 | `a` | $X_{ab;c} = -X_{ba;c}$ | 9 | 3 x 3（`a` x `c`） |
| 4 | `sss` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = X_{cd;ab}$ | 21 | 6 x 6（`v` x `v`） |
| 4 | `ssa` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = -X_{cd;ab}$ | 15 | 6 x 6（`v` x `v`） |
| 4 | `aas` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = X_{cd;ab}$ | 6 | 3 x 3（`a` x `a`） |
| 4 | `aaa` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = -X_{cd;ab}$ | 3 | 3 x 3（`a` x `a`） |
| 4 | `sa` | $X_{ab;cd} = X_{ba;cd} = -X_{ab;dc}$ | 18 | 6 x 3（`v` x `a`） |
| 4 | `as` | $X_{ab;cd} = -X_{ba;cd} = X_{ab;dc}$ | 18 | 3 x 6（`a` x `v`） |
| 4 | `s` | $X_{ab;c;d} = X_{ba;c;d}$（`sss` + `ssa` + `sa`） | 54 | 6 x 9（`v` x `u`） |
| 4 | `a` | $X_{ab;c;d} = -X_{ba;c;d}$（`aas` + `aaa` + `as`） | 27 | 3 x 9（`a` x `u`） |
| 4 | `t` | $X_{abc;d}$ は $a, b, c$ について完全対称 | 30 | 10 x 3（`t` x `c`） |

- テンソルの各成分を多極子（既約なデカルト成分）に分解します。0 でない値をとれるのは、群の active multipole（時間反転を含むすべての操作で不変な多極子）だけです。これから、0 でない成分とその間の関係が決まります。
- `(X, rank, opt)` の値：成分の行列（$Q_{xyz}$ のような記号。0 になる成分は `0`。行と列は上の表のとおり）と、各成分を active multipole で書いた式です。

*例。*

- 4/mmm（`MPG:15.1.53`）の `(Q, 2, s)`（例えば電気伝導度）：$\mathrm{diag}(Q_{xx}, Q_{xx}, Q_{zz})$ で、$Q_{xx} = \mathbb{Q}_0^{(1)} - \mathbb{Q}_u^{(1)}$、$Q_{zz} = \mathbb{Q}_0^{(1)} + 2\mathbb{Q}_u^{(1)}$。
- 4/mm'm'（`MPG:15.6.58`、`z` 方向の強磁性体）の `(T, 2, a)`（例えば異常ホール効果）：0 でない成分は $T_{xy} = -T_{yx} = \mathbb{M}_z^{(1)}$ だけで、active な磁気双極子 $M_z$ で与えられます。

### 7. 親の群での剰余類

親 `P`（`data_parent_pg/sg`）を持つ群 `H`（実の点群または空間群）について、`P` を左剰余類 `g H` に分けます。代表
`g` は `H` の構造を別の構造に写し、その群は `g H g^-1` です。`P` で不変な構造（親構造）では、`g` は対称操作です。
MultiPieGen は、親を持つすべての群（複素の点群を除く）について次を与えます。

- **親と剰余類。** トップの群（Oh、D6h、SG 221–230、191–194）までの親の列（すぐ上の親が先）と、各親での左剰余類
  （各剰余類にある親の対称操作。SG では並進を `H` の座標で、`H` の格子を法として書きます）。`H` 自身と並進だけの
  剰余類は除きます。
- **最初の剰余類。** すぐ上の親での、`H` 以外の最初の剰余類、その代表 `g`（最初の操作）、`g` の位数 `n`
  （`g^n` が `H` に入る最小の `n`）、`g` の種類（回転、または空間反転。SG では並進を伴うこともあります）。空間反転は
  `H` がカイラルなときだけに現れます（77 群：点群 12、空間群 65）。
  - 回転：像の構造は、`H` の構造を `g` で動かした座標系で書いたものです。すべての SAMB で `g Z_j = +Z'_j` です
    （`Z'_j` は像の同じラベルの SAMB）。
  - 空間反転：像の群 `H' = I H I^-1`（鏡像異性体、または `H` 自身）、`H` の座標から `H'` の座標への写像 `M`、
    各 Wyckoff のサイトとボンドの像（Wyckoff 位置、サイトやボンドの順、ボンドの反転）、クラスター SAMB の像。
    クラスター SAMB はラベルごとに `I Z_j = eps_j Z'_j`（`eps = +-1`）で写ります。ただし自分自身に写る Wyckoff
    位置の一部（38 個。例：SG 210 の 48f）では、既約表現の中で `I Z_a = sum_b U_ab Z'_b` です。
- **親構造。** 自由なパラメータを持つ Wyckoff 位置ごとに、`H` の軌道がサイトを重ねずにすぐ上の親で不変になる
  パラメータの値（その最大の集合）、そこで占める親の Wyckoff 位置、それを保つ剰余類。

*例（Te、P3_121 = No.152）。* すぐ上の親は P-3m1（No.164）、最初の剰余類の代表は `g = -1`（`n = 2`）で、その像は
`x -> -x` で P3_221（No.154）です。`3a` のサイトとボンドの SAMB はすべて `eps = +1` で写ります。`3a` の親構造は
`x = 0`（P-3m1 の Wyckoff 位置 `1a`。その単位胞は `H` の 3 分の 1）です。

(ja-multipole-key)=
### 8. 多極子のキー（`PGMultipoleType`）

この節では、調和関数や SAMB を使うときに必要な、多極子のキーを説明します。

点群の多極子（点群調和関数、点群多極子調和関数、原子 SAMB、クラスター SAMB）は、すべて `multipie` パッケージの名前付きタプル `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` で表します。

| フィールド | 意味 | 値 |
|---|---|---|
| `X` | タイプ。空間反転（$\mathcal{P}$）と時間反転（$\mathcal{T}$）での偶奇で決まる | `Q` 電気（$\mathcal{P} = (-1)^l$、$\mathcal{T}$ 偶）、`G` 電気トロイダル（$(-1)^{l+1}$、偶）、`T` 磁気トロイダル（$(-1)^l$、奇）、`M` 磁気（$(-1)^{l+1}$、奇） |
| `l` | ランク | 0, 1, 2, ... |
| `Gamma` | 点群の既約表現 | 指標表での名前。例：`A1g`、`Eu`、`T2g` |
| `n` | 点群調和関数 `(X, l, Gamma)` の多重度。同じ `(X, l, Gamma)` の調和関数のうちどれか | 1, 2, ...。1 つだけなら `-1` |
| `p` | 追加の多重度。ほかのフィールドがすべて同じものに番号を付ける | 1, 2, ...。1 つだけなら `-1` |
| `s` | 内部ランク（内部自由度のランク） | 0（スカラー）、1、2、3 |
| `k` | 内部成分。軌道部分のランクは `l + k` | `-s, ..., s` |
| `x` | 内部タイプ。内部基底の偶奇 | `q`（極性）、`g`（軸性） |

- `(X, l, Gamma, n)` は、常に 1 つの点群調和関数（その xyz の多項式と $U$）を指します。`T` と `M` は、`Q` と `G` の調和関数を使います。`n` は振り直さないので、別のデータでも `(X, l, Gamma, n)` が同じなら同じ調和関数です。
- `p` が使われるのは、`(X, l, Gamma, n, s, k, x)` が同じものが複数あるときだけです（例：ある Wyckoff 位置で、同じ調和関数に対応するクラスター SAMB が 2 つある場合）。
- スカラーの量（点群調和関数、クラスター SAMB）は `s = 0, k = 0, x = "q"` です。
- キーは多重項全体を表します。値には、すべての成分を既約表現の成分の順に（配列の行として）まとめて格納しています。

#### LaTeX 記号

ドキュメントと PDF では、多極子をキーから作った LaTeX 記号で書きます（`multipie.util.util_tag.TagMultipole.latex` に、下の tag を付けたもの）。

$$
\mathbb{X}_{l,\gamma}^{(s,k;\,\mathrm{tag})}[x](\Gamma, n\,p)
$$

| 部分 | キーとの対応 | 書く場合 |
|---|---|---|
| $\mathbb{X}$ | `X` | $\mathbb{Q}, \mathbb{G}, \mathbb{T}, \mathbb{M}$ |
| $l$ | `l` | 常に |
| $\gamma$ | 成分（1, 2, ...） | 多重項の 1 つの成分を指すときだけ。多重項全体や 1 次元の既約表現では省く |
| $(s, k)$ | `s`, `k` | `s > 0` のときだけ |
| tag | データの種類 | 調和関数ではなし。原子 SAMB では `a`。クラスター SAMB では、`site` が `s`、`bond` が `b` |
| $[x]$, $\vec{\ }$ | `x`, `s` | 点群多極子調和関数（`s > 0`）で、$\mathbb{X}$ に矢印を付け、内部タイプ $[q]$ または $[g]$ を書く |
| $\Gamma$ | `Gamma` | 常に。例：$T_{1u}$。複素の既約表現 `Ea` は $E^{(a)}$ と書く |
| $n$ | `n` | `n != -1` のときだけ |
| $p$ | `p` | `p != -1` のときだけ。$n$ の後ろに英字で書く（`1 -> a, 2 -> b, ...`） |

つまり、上付きの括弧と最後の括弧から `s, k` と `n, p` の値が読み取れます。書かれていないフィールドは既定値（`n = p = -1`、`s = k = 0`、`x = q`）です。

例：

| キー | 記号 | 意味 |
|---|---|---|
| `("Q", 1, "T1u", -1, -1, 0, 0, "q")` | $\mathbb{Q}_{1}(T_{1u})$ | Oh の電気双極子。成分は $(x, y, z)$。その `x` 成分は $\mathbb{Q}_{1,1}(T_{1u})$ |
| `("Q", 2, "Bg", 2, -1, 0, 0, "q")` | $\mathbb{Q}_{2}(B_{g},2)$ | C4h の $B_g$ の 2 番目（$n = 2$）の四極子、$\sqrt{3}xy$ |
| `("G", 0, "A1u", -1, -1, 1, 1, "g")` | $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | 点群多極子調和関数：ランク 1 の軌道部分と軸性ベクトルから作る電気トロイダル単極子 |
| `("M", 1, "T1g", -1, -1, 1, -1, "q")` | $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ | 原子 SAMB：スピンの磁気双極子 |
| `("T", 1, "Eu", -1, 1, 0, 0, "q")` | $\mathbb{T}_{1}^{(b)}(E_{u},a)$ | ボンドのクラスター SAMB：D4h の反対称ボンドの関数の 1 つ。$(x, y)$ の対称性を持つ（$p = 1$） |
| `("Q", 2, "E", 1, 2, 0, 0, "q")` | $\mathbb{Q}_{2}(E,1b)$ | $n = 1$、$p = 2$ |

## 生成規則：アルゴリズム、キー、軌道基底

この部では、MultiPieGen が多極子のキーを使うデータをどう作っているか、そしてキーと軌道基底の規則を、1 か所にまとめます。このパッケージのアルゴリズムだけを説明し、各節に処理を実装しているコードの場所を書いています。MultiPie Ver.2 からの変更点は、「詳細」の章の「MultiPie Ver.2 からの変更点」に別にまとめています。

目次：

1. 共通の規則（キー、既約表現、トップの群、軌道名）
2. 点群の調和関数
3. multipolar harmonics（`multipolar_harmonics.pkl`）
4. トップの群の cluster 基底（`cluster_scalar_basis.pkl`、`cluster_vector_basis.pkl`、`cluster_bond_basis.pkl`）
5. 各 Wyckoff のサイトとボンドの cluster SAMB（`cluster_samb.pkl`）
6. 原子多極子（`atomic_multipole.pkl`）と点群の SAMB（`atomic_multipole_group.pkl`）
7. 厳密計算と式の簡約
8. MultiPie 向けの binary データ（`binary_data`）

---

### 1. 共通の規則

#### 1.1 多極子のキー

多極子のキーは、すべて `multipie` パッケージの `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` です。対称性に合わせる前の原子多極子
には `SphericalMultipoleType(X, l, s, k, x)` を使います。

| 項目 | 意味 |
|---|---|
| `X` | 種類：`Q`（電気）、`G`（電気トロイダル）、`T`（磁気トロイダル）、`M`（磁気）。調和関数と cluster のデータは `Q` と `G` だけ。 |
| `l` | ランク |
| `Gamma` | 点群の既約表現（`data_character_table` の名前） |
| `n` | その点群の調和関数 `(X, l, Gamma)` の多重度。**1 から数え**、ただ 1 つなら `-1`。常に調和関数の xyz の多項式を特定する（振り直さない）。 |
| `p` | 追加の多重度。`(X, l, Gamma, n)` とほかの項目が同じものを番号付けするときだけ使う。ただ 1 つなら `-1` |
| `s` | 内部ランク（スカラー関数は 0） |
| `k` | 内部成分、`-s <= k <= s` |
| `x` | 内部の種類。`q`（極性の内部基底）か `g`（軸性の内部基底） |

多重項のすべての成分は、成分の順に行として、1 つの値にまとめて入れます。

そのほかのキー：

- 調和関数のキー `(l, irrep, n, comp)`：`n` は 1 から（ただ 1 つなら `-1`）、`comp` は 0 から（1 成分なら `-1`）。
- tesseral のキー `(l, m)`：`m >= 0` は cos 型の `C_{l,m}`、`m < 0` は sin 型の `S_{l,|m|}`。
- mapping の中の対称操作の番号は 1 から数える。

#### 1.2 既約表現

- 並べ替えでは、`data_character_table` の既約表現の並びを正とします。`g` が `u` より先、`'` のまとまりが `''` の
  まとまりより先です（C3h、D3h、D3h-1、C3h-c もこの順です）。
- 複素の既約表現（`-c` の点群）は名前が `a`/`b` で終わります（`Ega`、`T1ua` など）。実の点群では、複素型の実 2 次元
  既約表現（C3、C4、C6、S4、C3i、C4h、C6h、T、Th の `E` など）を使います。

#### 1.3 トップの群

- 調和関数は、2 つのトップの群 Oh と D6h について一度だけ定義し、47 点群（結晶点群 32、別設定 `-1` が 5、複素
  `-c` が 10）のそれぞれに対応させます。
- 立方晶の点群は Oh を、三方晶と六方晶の点群は D6h を使います。それより対称性の低い群は、親の群をたどって着く
  トップの群を使います（`data_harmonics_mapping[pg]["top"]`）。

#### 1.4 軌道名

軌道名は `data_basis` が正です（`data_harmonics` には古い名前の表があるので、名前には使いません）。

軌道名は、それぞれ調和関数のキーを指します（`data_basis.harmonics_key`）。tesseral のキーか、Oh/D6h のトップの
キーです。式と u-matrix は `info["harmonics"]["basis_function"][name] = (式, u-matrix)` で、u-matrix は
`m = l, l-1, ..., -l` の列 `<m|orbital>` です。

| `data_basis` のリスト | 軌道（ランク 0..3） | 用途 |
|---|---|---|
| `atomic_orbital_basis`（tesseral の並び） | s; px, py, pz; dv, dxy, dxz, dyz, du; fc, fs, fbz, f3, fx, fy, faz | 原子の基底（Lg/Lgs） |
| `cubic_orbital_basis` | s; px, py, pz; du, dv, dyz, dxz, dxy; f3, fax, fay, faz, fbx, fby, fbz | 立方晶の群（Oh 系）の内部基底 |
| `hexagonal_orbital_basis` | s; px, py, pz; du, dxz, dyz, dxy, dv; faz, fs, fc, fx, fy, f3, fbz | 三方・六方晶の群（D6h 系）の内部基底 |

LaTeX の名前：原子の基底には `atomic_latex_name`（`p_{x}`、`f_{z}^{\alpha}` など）を使います。multipolar harmonics
の内部基底の記号には `orbital_latex_name`（`_{x}`、`_{3c}` など）に Q/G を付けたものを使います
（`info["harmonics"]["internal_basis"]`）。

スピンを含む基底（`data_basis`）：スピンは `["u", "d"]`。LM は `M = L..-L`。JM は `J` が昇順、`M` が降順。LMs と Lgs は
`(m または軌道, スピン)` で、スピンがいちばん内側で回ります。

---

### 2. 点群の調和関数

データ：`data_top_harmonics`（トップの調和関数）、`data_compatibility_relation`。
コード：`util_top_harmonics`（トップの調和関数の式と u-matrix）、`util_harmonics_mapping`（各点群への対応
`data_harmonics_mapping`）、`util_harmonics_multipole.pg_harmonics`。

1. **トップの調和関数。** Oh と D6h について、ランク 11 までのすべての調和関数を、デカルト座標の多項式
   （`top_harmonics_symmetry[top][top_key]`）と、球面調和関数への u-matrix（`top_harmonics_u_matrix[top][top_key]`、
   `m = l..-l`）で与えます。持っているのは極性の調和関数だけです。同じ関数の軸性のものは、既約表現の `g` と `u`
   を入れ替えたものになります。
   - 実と複素の組み合わせが混在しています。複素の既約表現は名前が `a`/`b` で終わり、1 成分です。
   - 2 次元既約表現の基底は、いちばん単純な `(x, y)` の基底で `C_n` の表現行列が回転行列になるように選んでいます。
2. **対応。** `data_harmonics_mapping[pg] = {"top": top, "Q": {pg_key: (c, top_key)}, "G": {...}}`、`c = ±1`。
   点群の調和関数はトップの調和関数の `c` 倍、`X(pg_key) = c T(top_key)` で、u-matrix は `c U_T` です。
   - 適合関係（`data_compatibility_relation`）から作ります。トップの調和関数の各成分 `(irrep, comp)` は、部分群の
     極性の既約表現・成分と、軸性の既約表現・成分に、係数つきで対応します。部分群の既約表現が複素のときは、
     代わりに複素のトップの調和関数（`irrep + "a"/"b"`）を参照します。
   - 部分群の多重度 `n` は、`top_harmonics_tag` の順に振ります（ただ 1 つなら `-1`）。
   - ある群の極性 `Q` と軸性 `G` は、同じトップの調和関数から来ます。極性の調和関数に対応する軸性のものは、
     `top_key` が一致するもので、符号は `cQ cG` です（`util_harmonics_symbol.polar_to_axial`）。部分群では適合関係
     だけでは足りません。対応がトップの調和関数のパリティによるからです。
3. **表現行列** は u-matrix から求めます（`util_representation_matrix`）。1 つの既約表現のすべての調和関数が同じ
   行列を与えます（ランク 11 までテスト済み）。
4. **binary（`harmonics.pkl`）。** `{PG_id: {PGMultipoleType(X, l, Gamma, n, -1, 0, 0, "q"): (式 (dim), u-matrix
   (2l+1, dim))}}`。47 点群、ランク 11 まで、並びは Gamma、l、Q/G、n です（`util_harmonics`、
   `create_binary/create_harmonics.py`）。`PGMultipoleType` をキーに持つほかの binary は、調和関数の式も u-matrix も
   持ちません。`util_harmonics.harmonics_expression(pg, key)` と `harmonics_u_matrix(pg, key)` で、どのキーからでも
   その `(X, l, Gamma, n)` で引けます。T は Q、M は G の調和関数を引きます（空間群はその点群で引く）。

---

### 3. multipolar harmonics（`multipolar_harmonics.pkl`）

コード：`util_harmonics_multipole.create_multipolar_harmonics`、`util_multipolar_harmonics`（Olm）、
`create_binary/create_multipolar_harmonics.py`。範囲：s = 1, 2, 3、l <= 3、47 点群。

1. **多極子球面調和関数。** `O_{lm}^{(s)}(k)`（`util_multipolar_harmonics.Olm`）は、軌道ランク `l+k` の球面調和関数と、
   ランク `s` の内部基底を結合したものです：
   `O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} sqrt(2l+1) sum_n (l+k, l, s; m-n, -m, n) O_{l+k, m-n} e_n`。
   ここで `O_{lm} = sqrt(4 pi/(2l+1)) r^l Y_lm`、`e_n` は内部基底のベクトルで、選んだ内部基底（`InternalBasis`）で
   表します。
2. **内部基底。** `data_basis.cubic_orbital_basis[s]`（トップが Oh）または `hexagonal_orbital_basis[s]`（トップが D6h）
   の軌道で、u-matrix は `info["harmonics"]["basis_function"]` から取ります。
3. **射影。** 群の各調和関数 `(X, l, Gamma, n)`（u-matrix `U`、2 節）と各 k について、内部係数
   `U^T O_l^{(s)}(k)` を作ります。x, y, z の多項式の `(dim, 2s+1)` 配列（各要素は因数分解）で、内部基底の各関数に掛かる
   係数です。0 になるもの（`w3j = 0`）は除きます。
4. **内部の種類。** (`X = Q` かつ `s+k` が偶数) または (`X = G` かつ `s+k` が奇数) なら `x = q`、それ以外は `x = g`。
5. **式。** 多極子の式 `= 内部係数 @ symbols`。内部の記号は `Q/G + orbital_latex_name`
   （`info["harmonics"]["internal_basis"][crystal][(s, x)]`）。
6. **キー。** `PGMultipoleType(X, l, Gamma, n, p, s, k, x)`：`n` は調和関数の多重度（振り直さない）。`p` は
   `(x, l, X, Gamma, n, k)` が同じものを番号付けします（ただ 1 つなら `-1`）。
7. **並び。** `x`（q, g）、`s`、`Gamma`、`l`、`X`（Q, G）、`k`、`n`、`p`。
8. **値。** `(internal_coefficient (dim, 2s+1), multipole_expression (dim))`：内部係数と多極子の式。調和関数
   `(X, l, Gamma, n)` そのものの式は `harmonics.pkl` にあります。

---

### 4. トップの群の cluster 基底

コード：`util_cluster_basis`、`util_cluster_vector_basis`、`util_bond_cluster_basis` と、対応する `create_binary` の
モジュール。データ：`data_root_cluster`、`info["root_cluster"]`。

#### 4.1 クラスター

- Oh のクラスターは 48 点、D6h のクラスターは 24 点で、デカルト座標に `x, y, z` を入れる形で書いてあります
  （`data_root_cluster`）。点 `i` は `g_i site0` で、`g_i` は Oh または D6h の `i` 番目の対称操作です（`group.pkl` の順）。
- Oh の `site0` は `(1, 5, 6)`（正規化していない、長さ `sqrt(62)`）。D6h の `(2 sqrt(2), sqrt(6), 4 sqrt(3))` は、Oh
  のものを主軸が `[111]` に向くよう回したものです。`vector0`（ボンドの向き）は `(6, 1, 5)` の単位ベクトルです（D6h
  では同じように回したもの）。`(5, 1, 6)` は使いません。`x + y = z` の面について `site0` の鏡像になっていて、いくつかの
  ボンド関数が 0 になったからです。
- 一般点は Oh と D6h のどちらでも一般の位置にあるので、クラスターの点の数は群の位数と同じで、クラスター上の関数は
  正則表現になります。

#### 4.2 スカラーの cluster 基底（`cluster_scalar_basis.pkl`）

1. 実のトップの調和関数 `T`（ランク 0..11）のそれぞれについて、クラスターの点での値を並べたベクトル。
2. トップの群の (既約表現, 成分) ごとに、低いランク、次に小さい多重度の順に Gram–Schmidt します。一次従属なものは
   多重項ごと除きます（`"dropped"` に記録）。
3. 結果：48 本（Oh）または 24 本（D6h）の正規直交基底で、正則表現です。
4. キー：`PGMultipoleType("Q", l, Gamma, n, -1, 0, 0, "q")`。`(l, Gamma, n)` はトップの調和関数のもの。

#### 4.3 ベクトルの cluster 基底（`cluster_vector_basis.pkl`、s = 1）

ボンドの cluster 基底を作るための中間の基底です。

1. 実のトップの調和関数 `T` と `k = -1, 0, 1` のそれぞれについて、内部基底 `[px, py, pz]` のベクトル調和関数
   `sum_m U_T[m] O_{lm}^{(1)}(k)` を、クラスターの点で評価したもの（`(N, 3)`）。
2. 種類：`x = q` では、`s + k` が偶数なら `X = Q`、奇数なら `G`。`x = g` ではその逆。極性のものは `T` の既約表現・
   成分のまま。軸性のものは `data_harmonics_mapping[top]["G"]` の既約表現・成分と係数 `c` を使います。
3. (既約表現, 成分) ごとに、ランク `l`、次に `k = -1, 0, 1`、次に多重度の順に Gram–Schmidt します。一次従属なものは
   除きます（数値で判定し、残したものを厳密に正規直交化）。
4. 結果：各 `x` について `3N` 本、つまり 144 本（Oh）または 72 本（D6h）。
5. キー：`PGMultipoleType(X, l, Gamma, n, -1, 1, k, x)`（`p = -1`。多重項は `k` で区別される）。

#### 4.4 ボンドの cluster 基底（`cluster_bond_basis.pkl`）

1. ボンドベクトル：`b_i` はクラスターの点 `i` の `x, y, z` を `vector0` に置き換えたもの、つまり `g_i vector0`。
2. ベクトルの cluster 基底 `v_key`（`(N, 3)`）のそれぞれについて、ボンド関数 `f_key(i) = v_key(i) . b_i`（同じキー）。
   `N` 次元の空間の中の `3N` 本のベクトルになります。
3. (既約表現, 成分) ごとに、ランク、`k = -1, 0, 1`、多重度の順に Gram–Schmidt します。一次従属なものと 0 になる
   ものは除きます。結果：各 `x` について `N` 本。
4. `f_key` はクラスター上のスカラー関数として変換します（`x = g` では符号 `det R` が付く）。
5. キー：`PGMultipoleType(X, l, Gamma, n, p, 0, 0, x)`。`(X, l, Gamma, n)` は元の
   ベクトル基底のキーのもの。`(X, l, Gamma, n)` が同じ多重項には `p = 1, 2, ...`（`k` の順）、ただ 1 つなら `-1`。
   `x` は元のもの。`"source"` に、各キーの元のベクトルのキー（`k` を含む）を入れています。

---

### 5. 各 Wyckoff のサイトとボンドの cluster SAMB（`cluster_samb.pkl`）

コード：`util_cluster_samb`、`create_binary/create_cluster_samb.py`。範囲：実の点群 37 と空間群 230。`"site"`、
`"bond"`。ボンドの SAMB は、対称ボンドの SAMB（以下 `"bond_s"`）と反対称ボンドの SAMB（以下 `"bond_a"`）から
作ります（7. を参照）。複素の点群は含みません。

1. **一般点。** 群 `H` の一般点のクラスターは、トップのクラスターの一部です。`H` の点 `i` は、
   `util_op_mapping.op_mapping_d6h_oh_cluster[PG]` の coset 0 の `i` 番目にあるトップの点です。
   - トップの基底をこれらの点に制限したものが、`H` の一般点の基底です。`"site"` と `"bond_s"` には
     `cluster_scalar_basis` を、`"bond_a"` には `cluster_bond_basis`（`x = q` だけ。それだけで完全な基底）を使います。
     一般のボンドでは、対称ボンドと反対称ボンドは同じ空間を張ります。
   - トップの群のキーは引き継ぎますが、既約表現と成分はトップの調和関数を通して `H` のものに付け替えます。
     トップの群のキー `= c_top T`、`H` のキー `= c T` なので、係数は `c c_top` です。
2. **空間群。** 空間群の対称操作（plus set を除く）は、その点群と同じ点操作を同じ順で持ちます。各操作は点群を
   通して対応させます（`op_mapping_parent_sg` は軸を変えて親の空間群に対応させるもので、使いません）。
3. **特殊点。** Wyckoff のサイト（ボンド）の mapping `M_j`（代表をサイトまたはボンド `j` に移す操作、`group.pkl`
   にある）を使い、`f_j = sum_{i in M_j} f_i` とします。`"bond_a"` では、ボンドの向きを反転させる操作（mapping で
   負の番号）は `-f_i` として足します。
4. **正規直交化。** `H` の既約表現ごとに Gram–Schmidt します（成分ごとに順に。複素型の実 2 次元既約表現では、
   別の多重項の成分どうしが直交しないため）。順番は、`"site"` と `"bond_s"` が `[Gamma, l, n, p]`、
   `"bond_a"` が `[Gamma, l, k, Q/G, n, p]`（`k` は元のベクトル基底のもの）です。一次従属なものと 0 になるものは除き
   ます（数値で選び、残したものを厳密に正規直交化）。各成分は `N` 個のサイト（ボンド）上でノルム 1 です。
5. **キー。** `"site"` と `"bond_s"`：`PGMultipoleType("Q", l, Gamma, n, p, 0, 0, "q")`。`"bond_a"`：
   `cluster_bond_basis` のキー（`X = Q/G`、`x = q`）。`n` は `H` の調和関数の多重度、`p` はトップのキーのもの。
6. **値。** `samb (dim, N)`。点群の調和関数の式は `harmonics.pkl` にあります。
7. **ボンドの SAMB。** `"bond" = "bond_s" + i "bond_a"` です。反対称ボンドの SAMB に `i` を掛け（虚数の飛び移り）、
   タイプを `Q` から `T` に、`G` から `M` に付け替えます（時間反転に奇）。`"bond_s"`（`Q`）と `"bond_a"`（`T/M`）
   のキーは決して重ならないので、1 つにまとめます。時間反転の偶奇が異なる SAMB は（エルミートなボンド行列として、
   実対称と虚反対称で）直交するので、`Q` と `T/M` は別々に正規直交化します（4.）。順番は `Gamma` の順で、各
   `Gamma` の中では `Q`（`[l, n, p]`）を `T/M`（`[l, k, T/M, n, p]`）より先に並べます。値は `Q` で実数、`T, M` で
   純虚数です。`T`（`M`）の調和関数は `Q`（`G`）のものを使います（`util_harmonics.harmonics_type`）。
8. **空間反転の剰余類のための基準サイト。** 特殊点の SAMB は、平均をとる元の代表（0 行目）によって変わります。
   最初の剰余類が空間反転の群（9. を参照）では、Wyckoff 位置 `W` の代表の像が、像 `W'` の代表になるとは限りません。
   そのとき、同じ既約表現の SAMB が空間反転で混ざることがあります。そこで、`W != W'`（または `H != H'`）の組で
   SAMB がまだラベルごとに写らないものは、（群、Wyckoff の順）が後の方を、もう片方の代表の像から平均して作ります
   （`util_cluster_samb_image.cluster_samb_reference`、`create_cluster_samb(..., reference=...)`）。対象は SG 98、
   208、211、213、214 の 20 個の Wyckoff 位置です。自分自身に写る Wyckoff 位置はデータを変えません（その行列 `U` は
   `coset.pkl` にあります）。

---

### 6. 原子多極子

コード：`util_atomic_multipole`、`create_binary/create_atomic_multipole.py`、
`create_binary/create_atomic_multipole_group.py`。

#### 6.1 原子多極子の行列（`atomic_multipole.pkl`）

1. 基底（`data_basis`）：`L = 0..3` のスピンを含む LMs（`atomic_lms_basis`、`(m, spin)`）と JML（`atomic_jm_basis`、
   `(J, M)`）。スピンを含まない LM は上向きスピンの部分です。
2. 行列要素 `<b1|X_{lm}^{(s)}(k)|b2>`（`X = Q, G, T, M`、`s = 0, 1`、`l <= 7`）を Wigner–Eckart の定理で求めます
   （`multipie.util.util_atomic_multipole`）。Jz の保存で許される要素（`m1 + s1 = m2 + s2 + m`、`M1 = M2 + m`）だけを
   計算し、パリティや三角条件で禁止されるブロックは飛ばします。
3. `data_basis.atomic_rank_block`（上三角）の各 bra-ket ブロック `(L1, L2)` を、0 でなければ残します。
4. キー：`SphericalMultipoleType(X, l, s, k, "q")`。並びは `Q, G, T, M`、`s`、`k`、`l`。値：`m = l..-l` の
   `(2l+1, n1, n2)`。`"lm"` は `"lms"` の `s = 0`、上向きスピンの部分です。

#### 6.2 点群の SAMB（`atomic_multipole_group.pkl`）

1. **トップの SAMB。** トップの調和関数 `T`（u-matrix `U_T`、`m = l..-l`）の各成分と、パリティが合う各原子多極子
   （極性の調和関数には `Q, T`、軸性には `G, M`）について、`sum_m U_T[m] X_{lm}^{(s)}(k)` を作り、Frobenius ノルムが
   1（対角ブロック）または `1/sqrt(2)`（非対角ブロック）になるよう正規化します。トップの群ごとに一度だけ計算します。
2. 成分どうしは直交しています（球テンソルが直交し、`U_T` の列が正規直交なため）。そのため、この正規化は (`X`, 既約
   表現) ごとの Gram–Schmidt による正規直交化と同じです。非対角ブロックでは、同じ既約表現の Q と T（G と M）は内積の
   虚部だけが重なります（エルミート共役のブロックと合わせると直交します）。
3. **点群。** 群の調和関数 `(X, l, Gamma, n)` の SAMB は、そのトップのキーのトップの SAMB の `c` 倍で
   （`data_harmonics_mapping`）、成分を重ねたものです。
4. **基底。** `"jml"` は JML の行列から作ります。`"lgs"` は LMs の行列 `M` から `U_bra^dagger M U_ket` で作ります。`U`
   は `data_basis.atomic_lgs_basis` の軌道（`atomic_orbital_basis`、基底関数の u-matrix をスピンに広げたもの）です。
   `"lg"` は `"lgs"` の `s = 0`、上向きスピンの部分に `sqrt(2)` を掛けたものです。
5. **キー。** `PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")`。`n` は調和関数の多重度（振り直さない）。
   並びは `Q, G, T, M`、`s`、`k`、`l`、`Gamma`、`n`。
6. **値。** ブロックごとに `行列 (dim, n1, n2)`（調和関数の式は `harmonics.pkl` にある）。47 点群すべて。

---

### 7. 厳密計算と式の簡約

- 保存する値はすべて sympy です。保存する前に `util_symmetry_operation.simplify_value`
  （`expand(radsimp(expand(x)))`、有理数はそのまま）を通し、分母を有理化します（`1/sqrt(3)` ではなく `sqrt(3)/3`）。
- 一次従属かどうかは数値で判定し（相対許容誤差 `1e-8`）、残したベクトルを厳密に正規直交化します。
- pickle を読み込むと、sympy が一部の積を評価し直します（`(a+b)/2` → `a/2 + b/2`）。そのため「binary is up to date」
  のテストは、pickle に一度書き出して読み直した後で比べます。

---

### 8. MultiPie 向けの binary データ（`binary_data`）

コード：`util_binary_data`、`create_binary/create_binary_data.py`。

1. **構成。** MultiPie の `multipie/binary_data` と同じく、`info.pkl`（`binary/info.pkl` と同じデータ。ヘッダに構成の説明を追加）と、群ごとのファイル `{PG,SG,MPG,MSG}/{id}.pkl`（47、230、157、1651 個。ファイル名では id の `:` を `_` に置き換えます。例：`PG/PG_32.pkl`。Windows の git は `:` を含むファイル名を checkout できないためです）です。`_opt` のファイルには分けません。PG のファイルは最大でも約 260 KB です。
2. **中身。** `multipiegen/binary` の binary を群ごとに並べ直したもので、計算し直すのは PG の表現行列（`util_representation_matrix`、1 群 1 KB 未満）だけです。
   - PG：`info`、`symmetry_operation`、`wyckoff`、`character`、`representation_matrix`、`harmonics`、`multipolar_harmonics`、`atomic_samb`、`cluster_samb`（実の点群だけ）、`coset`（親を持つ実の点群だけ。`coset.pkl` の `"first"`）。
   - SG：`info`、`symmetry_operation`、`wyckoff`、`cluster_samb`、`coset`（親を持つ空間群だけ）。
   - `info.pkl` には `"parent"` と `"coset"`（`coset.pkl` の `"parent"` と `"coset"`）も入れます。
   - MPG：`info`、`symmetry_operation`、`wyckoff`、`active_multipole`、`response_tensor`。
   - MSG：`info`、`symmetry_operation`、`wyckoff`。
3. **入れないもの。** 中間データのクラスター基底（`cluster_scalar_basis`、`cluster_vector_basis`、`cluster_bond_basis`）と、球面の原子多極子（`atomic_multipole`）。
4. **調和関数。** 重複させません。SAMB と点群多極子調和関数の調和関数は、PG のファイルの `"harmonics"` から `(X, l, Gamma, n)` で引きます（T → Q、M → G）。SG と MPG では、`info.PG` の PG のファイルを使います。
5. **ヘッダ。** 各ファイルに、元の binary のヘッダから作ったヘッダ（`str()`）を付けます。そのため、`BinaryManager` だけでファイルの意味が分かります。

---

### 9. 剰余類（`coset.pkl`）

コード：`util_coset`、`util_cluster_samb_image`、`create_binary/create_coset.py`。範囲：親を持つ実の点群 35 と空間群
216。群と対称操作はこのリポジトリの表（`data_parent_pg/sg`、`data_wyckoff_site`、`group_symmetry_operation`）だけ
から取り、表は保存しません。

1. **剰余類。** 親 `P` を左剰余類 `g H` に分けます（SG では `H` の格子を法として。`H` の格子の並進でない `P` の
   並進は別の剰余類になります）。順番は、回転部分が恒等のものを先に、次に並進の長さの順です。親の列は
   `data_representative_parent_pg/sg` をたどり、変換 `x_P = R x_H` を掛け合わせます。
2. **像。** すぐ上の親の各剰余類について、像 `g H g^-1` を `P` の子のうち軸が同じものから探します（SG では原点の
   移動 `k` を許し、その中で自分自身に写る Wyckoff 位置が最も多いもの、次に最も短いものを選びます）。写像
   `M = R'^-1 k^-1 g R` で、`H` の各 Wyckoff 位置を像の Wyckoff 位置に、パラメータの対応とともに写します。
3. **最初の剰余類。** `data_parent_pg/sg[H][P]` の 2 番目の剰余類（`H` 以外の最初のもの）と、その代表（最初の操作）、
   位数 `n`。種類は `"rotation"`（det > 0）か `"inversion"` です。
4. **サイト、ボンド、クラスター SAMB の像（空間反転）。** `H` の各 Wyckoff 位置 `W` を、`M` でその軌道が写る `H'` の
   Wyckoff 位置 `W'` と組にします（軌道の大きさと自由パラメータの数が同じもの。格子と面心・体心の並進とボンドの
   反転を許します）。`W` の `i` 行目は `W'` の `perm[i]` 行目に写ります。こうして動かしたクラスター SAMB（サイト `i`
   の値をサイト `perm[i]` に移し、`T/M` ではボンドの反転の符号を付ける）を、既約表現ごとに（ボンドでは `Q` と
   `T/M` を分けて）`W'` のクラスター SAMB で展開します。対角なら `"eps"`（成分ごとの `+-1`）、そうでなければ `"U"`
   （できるところは厳密な `+-sqrt(有理数)`）です。77 群（1847 個の Wyckoff 位置）の結果：ラベルごとに写るものが
   1809（そのうち 20 は基準サイトを替えた後。5.8）、`U` が 38。
5. **親構造。** 自由なパラメータを持つ Wyckoff 位置ごとに、右剰余類 `H c` のすべての代表 `c` で `c p` が `H` の軌道に
   載る（このとき軌道は `P` で不変）パラメータの最大のアフィン部分空間を、軌道の大きさが `W` の多重度のもの
   （サイトが重ならないもの）に限って求めます。最初の探索範囲で見つかった極限を先に並べ、それ以外は、`H` の
   軌道が既出のものと同じでなければ後に加えます。
6. **ヘッダとキー。** `coset.pkl` のヘッダ（`str(load_coset())`）を参照。`binary_data` では、`"first"` を各群の
   ファイルの `coset` に、`"parent"` と `"coset"` を `info.pkl` に入れます。
