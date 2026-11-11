# 多極子データの生成：アルゴリズム、キー、軌道基底

MultiPieGen が多極子のキーを使うデータをどう作っているか、そしてキーと軌道基底の規則を、1 か所にまとめた確認用の
メモです（今のところパッケージのドキュメントではありません）。英語版 `generation_rules.md` と同じ内容です。各節に、
処理を実装しているコードの場所を書いています。「MultiPie」はすべて MultiPie Ver.2 を指します。MultiPie Ver.2 との
違いを数字つきでまとめたものは `multipie_v2_comparison/README.md` にあります。

目次：

1. 共通の規則（キー、既約表現、トップの群、軌道名）
2. 点群の調和関数
3. multipolar harmonics（`multipolar_harmonics.pkl`）
4. トップの群の cluster 基底（`cluster_scalar_basis.pkl`、`cluster_vector_basis.pkl`、`cluster_bond_basis.pkl`）
5. 各 Wyckoff のサイトとボンドの cluster SAMB（`cluster_samb.pkl`）
6. 原子多極子（`atomic_multipole.pkl`）と点群の SAMB（`atomic_multipole_group.pkl`）
7. 厳密計算と式の簡約

---

## 1. 共通の規則

### 1.1 多極子のキー

多極子のキーは、すべて MultiPie の `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` です。対称性に合わせる前の原子多極子
には `SphericalMultipoleType(X, l, s, k, x)` を使います。

| 項目 | 意味 |
|---|---|
| `X` | 種類：`Q`（電気）、`G`（電気トロイダル）、`T`（磁気トロイダル）、`M`（磁気）。調和関数と cluster のデータは `Q` と `G` だけ。 |
| `l` | ランク |
| `Gamma` | 点群の既約表現（`data_character_table` の名前） |
| `n` | その点群の調和関数 `(X, l, Gamma)` の多重度。**1 から数え**、ただ 1 つなら `-1`。常に調和関数の xyz の多項式を特定する（MultiPie と違い、振り直さない）。 |
| `p` | 追加の多重度。`(X, l, Gamma, n)` とほかの項目が同じものを番号付けするときだけ使う。ただ 1 つなら `-1` |
| `s` | 内部ランク（スカラー関数は 0） |
| `k` | 内部成分、`-s <= k <= s` |
| `x` | 内部の種類。`q`（極性の内部基底）か `g`（軸性の内部基底） |

多重項のすべての成分は、成分の順に行として、1 つの値にまとめて入れます。

そのほかのキー：

- 調和関数のキー `(l, irrep, n, comp)`：`n` は 1 から（ただ 1 つなら `-1`）、`comp` は 0 から（1 成分なら `-1`）。
- tesseral のキー `(l, m)`：`m >= 0` は cos 型の `C_{l,m}`、`m < 0` は sin 型の `S_{l,|m|}`。
- mapping の中の対称操作の番号は 1 から数える。

### 1.2 既約表現

- 並べ替えでは、`data_character_table` の既約表現の並びを正とします。`g` が `u` より先、`'` のまとまりが `''` の
  まとまりより先です（C3h、D3h、D3h-1、C3h-c をこの順に並べ替えました）。MultiPie は既約表現を名前の文字列で
  並べています。
- 複素の既約表現（`-c` の点群）は名前が `a`/`b` で終わります（`Ega`、`T1ua` など）。実の点群では、複素型の実 2 次元
  既約表現（C3、C4、C6、S4、C3i、C4h、C6h、T、Th の `E` など）を使います。

### 1.3 トップの群

- 調和関数は、2 つのトップの群 Oh と D6h について一度だけ定義し、47 点群（結晶点群 32、別設定 `-1` が 5、複素
  `-c` が 10）のそれぞれに対応させます。
- 立方晶の点群は Oh を、三方晶と六方晶の点群は D6h を使います。それより対称性の低い群は、親の群をたどって着く
  トップの群を使います（`data_harmonics_mapping[pg]["top"]`）。

### 1.4 軌道名

軌道名は MultiPie から変えています：f2→fc、f1→fs、f3x→fx、f3y→fy、ga→gc、gb→gs、gc→gc4、i→i3、l→l3、i0→ixy、
j→j3（`gc` は旧 `ga` を指すようになった点に注意）。名前は `data_basis` が正です。`data_harmonics` には古い名前の表が
残っています。

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

## 2. 点群の調和関数

データ：`data_top_harmonics_info`（トップの調和関数）、`data_harmonics_mapping`（各点群への対応）。
コード：`util_harmonics_multipole.pg_harmonics`。表の生成スクリプトは `others/creation/create_top_harm.py`、
`create_harmonics_mapping.py`。

1. **トップの調和関数。** Oh と D6h について、ランク 11 までのすべての調和関数を、デカルト座標の多項式
   （`top_harmonics_symmetry[top][top_key]`）と、球面調和関数への u-matrix（`top_harmonics_u_matrix[top][top_key]`、
   `m = l..-l`）で与えます。持っているのは極性の調和関数だけです。同じ関数の軸性のものは、既約表現の `g` と `u`
   を入れ替えたものになります。
   - 実と複素の組み合わせが混在しています。複素の既約表現は名前が `a`/`b` で終わり、1 成分です。
   - 2 次元既約表現の基底は、いちばん単純な `(x, y)` の基底で `C_n` の表現行列が回転行列になるように選んでいます。
   - MultiPie からの意図的な符号の変更：D6h の rank 4 `E2g` 多重度 2（2 成分とも）を −1 倍しています。これは D6h 系の
     部分群（rank 4 の E の多重項）と、調和関数から作るデータにそのまま引き継がれます。
   - 軌道 `gv`（Oh の rank 4 `Eg` の成分 1）も、MultiPie の `gv` に対して −1 倍しています（`data_response_tensor`）。
     これは軌道（と応答テンソルの多極子）の定義の変更で、トップの調和関数は変わりません。Oh 系の調和関数は
     MultiPie と一致します。
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

---

## 3. multipolar harmonics（`multipolar_harmonics.pkl`）

コード：`util_harmonics_multipole.create_multipolar_harmonics`、`util_multipolar_harmonics`（Olm）、
`create_binary/create_multipolar_harmonics.py`。範囲：s = 1, 2, 3、l <= 3、47 点群。

1. **多極子球面調和関数。** `O_{lm}^{(s)}(k)`（`util_multipolar_harmonics.Olm`）は、軌道ランク `l+k` の球面調和関数と、
   ランク `s` の内部基底を結合したものです：
   `O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} sqrt(2l+1) sum_n (l+k, l, s; m-n, -m, n) O_{l+k, m-n} e_n`。
   ここで `O_{lm} = sqrt(4 pi/(2l+1)) r^l Y_lm`、`e_n` は内部基底のベクトルで、選んだ内部基底（`InternalBasis`）で
   表します。
2. **内部基底。** `data_basis.cubic_orbital_basis[s]`（トップが Oh）または `hexagonal_orbital_basis[s]`（トップが D6h）
   の軌道で、u-matrix は `info["harmonics"]["basis_function"]` から取ります（MultiPie はどの群でも tesseral の基底を
   使います）。
3. **射影。** 群の各調和関数 `(X, l, Gamma, n)`（u-matrix `U`、2 節）と各 k について、`mharm = U^T O_l^{(s)}(k)`。
   多項式の `(dim, 2s+1)` 配列で、各要素は因数分解します。0 になるもの（`w3j = 0`）は除きます。
4. **内部の種類。** (`X = Q` かつ `s+k` が偶数) または (`X = G` かつ `s+k` が奇数) なら `x = q`、それ以外は `x = g`。
5. **式。** `m_ex = mharm @ symbols`。内部の記号は `Q/G + orbital_latex_name`
   （`info["harmonics"]["internal_basis"][crystal][(s, x)]`）。
6. **キー。** `PGMultipoleType(X, l, Gamma, n, p, s, k, x)`：`n` は調和関数の多重度（振り直さない）。`p` は
   `(x, l, X, Gamma, n, k)` が同じものを番号付けします（ただ 1 つなら `-1`）。
7. **並び。** `x`（q, g）、`s`、`Gamma`、`l`、`X`（Q, G）、`k`、`n`、`p`。
8. **値。** `(ex, mharm, m_ex)`：調和関数のデカルト座標の式、内部基底の係数、多極子の式。

---

## 4. トップの群の cluster 基底

コード：`util_cluster_basis`、`util_cluster_vector_basis`、`util_bond_cluster_basis` と、対応する `create_binary` の
モジュール。データ：`data_root_cluster`、`info["root_cluster"]`。

### 4.1 クラスター

- Oh のクラスターは 48 点、D6h のクラスターは 24 点で、デカルト座標に `x, y, z` を入れる形で書いてあります
  （`data_root_cluster`）。点 `i` は `g_i site0` で、`g_i` は Oh または D6h の `i` 番目の対称操作です（`group.pkl` の順）。
- Oh の `site0` は `(1, 5, 6)`（正規化していない、長さ `sqrt(62)`）。D6h の `(2 sqrt(2), sqrt(6), 4 sqrt(3))` は、Oh
  のものを主軸が `[111]` に向くよう回したものです。`vector0`（ボンドの向き）は `(6, 1, 5)` の単位ベクトルです（D6h
  では同じように回したもの）。`(5, 1, 6)` は使いません。`x + y = z` の面について `site0` の鏡像になっていて、いくつかの
  ボンド関数が 0 になったからです。
- 一般点は Oh と D6h のどちらでも一般の位置にあるので、クラスターの点の数は群の位数と同じで、クラスター上の関数は
  正則表現になります。

### 4.2 スカラーの cluster 基底（`cluster_scalar_basis.pkl`）

1. 実のトップの調和関数 `T`（ランク 0..11）のそれぞれについて、クラスターの点での値を並べたベクトル。
2. トップの群の (既約表現, 成分) ごとに、低いランク、次に小さい多重度の順に Gram–Schmidt します。一次従属なものは
   多重項ごと除きます（`"dropped"` に記録）。
3. 結果：48 本（Oh）または 24 本（D6h）の正規直交基底で、正則表現です。
4. キー：`PGMultipoleType("Q", l, Gamma, n, -1, 0, 0, "q")`。`(l, Gamma, n)` はトップの調和関数のもの。

### 4.3 ベクトルの cluster 基底（`cluster_vector_basis.pkl`、s = 1）

ボンドの cluster 基底を作るための中間の基底です。

1. 実のトップの調和関数 `T` と `k = -1, 0, 1` のそれぞれについて、内部基底 `[px, py, pz]` のベクトル調和関数
   `sum_m U_T[m] O_{lm}^{(1)}(k)` を、クラスターの点で評価したもの（`(N, 3)`）。
2. 種類：`x = q` では、`s + k` が偶数なら `X = Q`、奇数なら `G`。`x = g` ではその逆。極性のものは `T` の既約表現・
   成分のまま。軸性のものは `data_harmonics_mapping[top]["G"]` の既約表現・成分と係数 `c` を使います。
3. (既約表現, 成分) ごとに、ランク `l`、次に `k = -1, 0, 1`、次に多重度の順に Gram–Schmidt します。一次従属なものは
   除きます（数値で判定し、残したものを厳密に正規直交化）。
4. 結果：各 `x` について `3N` 本、つまり 144 本（Oh）または 72 本（D6h）。
5. キー：`PGMultipoleType(X, l, Gamma, n, -1, 1, k, x)`（`p = -1`。多重項は `k` で区別される）。

### 4.4 ボンドの cluster 基底（`cluster_bond_basis.pkl`）

1. ボンドベクトル：`b_i` はクラスターの点 `i` の `x, y, z` を `vector0` に置き換えたもの、つまり `g_i vector0`。
2. ベクトルの cluster 基底 `v_key`（`(N, 3)`）のそれぞれについて、ボンド関数 `f_key(i) = v_key(i) . b_i`（同じキー）。
   `N` 次元の空間の中の `3N` 本のベクトルになります。
3. (既約表現, 成分) ごとに、ランク、`k = -1, 0, 1`、多重度の順に Gram–Schmidt します。一次従属なものと 0 になる
   ものは除きます。結果：各 `x` について `N` 本。
4. `f_key` はクラスター上のスカラー関数として変換します（`x = g` では符号 `det R` が付く）。
5. キー（MultiPie の `bond_a` と同じ形）：`PGMultipoleType(X, l, Gamma, n, p, 0, 0, x)`。`(X, l, Gamma, n)` は元の
   ベクトル基底のキーのもの。`(X, l, Gamma, n)` が同じ多重項には `p = 1, 2, ...`（`k` の順）、ただ 1 つなら `-1`。
   `x` は元のもの。`"source"` に、各キーの元のベクトルのキー（`k` を含む）を入れています。
6. この作り方（内積を取ってから、この順で Gram–Schmidt）は MultiPieGen 独自のものです。MultiPie は、各ボンド関数に
   最初に重なる vector SAMB のキーを付けますが、その規則は使いません。

---

## 5. 各 Wyckoff のサイトとボンドの cluster SAMB（`cluster_samb.pkl`）

コード：`util_cluster_samb`、`create_binary/create_cluster_samb.py`。範囲：実の点群 37 と空間群 230。`"site"`、
`"bond_s"`（対称ボンド）、`"bond_a"`（反対称ボンド）。複素の点群と、MultiPie の vector SAMB は作りません。

1. **一般点。** 群 `H` の一般点のクラスターは、トップのクラスターの一部です。`H` の点 `i` は、
   `data_op_mapping.op_mapping_d6h_oh_cluster[PG]` の coset 0 の `i` 番目にあるトップの点です。
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
   別の多重項の成分どうしが直交しないため）。順番は、`"site"` と `"bond_s"` が `[Gamma, l, n, p]`（MultiPie と同じ）、
   `"bond_a"` が `[Gamma, l, k, Q/G, n, p]`（`k` は元のベクトル基底のもの）です。一次従属なものと 0 になるものは除き
   ます（数値で選び、残したものを厳密に正規直交化）。各成分は `N` 個のサイト（ボンド）上でノルム 1 です。
5. **キー。** `"site"` と `"bond_s"`：`PGMultipoleType("Q", l, Gamma, n, p, 0, 0, "q")`。`"bond_a"`：
   `cluster_bond_basis` のキー（`X = Q/G`、`x = q`）。`n` は `H` の調和関数の多重度、`p` はトップのキーのもの。
6. **値。** `(samb (dim, N), ex (dim))`。`ex` は点群の調和関数のデカルト座標の式です。

---

## 6. 原子多極子

コード：`util_atomic_multipole`、`create_binary/create_atomic_multipole.py`、
`create_binary/create_atomic_multipole_group.py`。

### 6.1 原子多極子の行列（`atomic_multipole.pkl`）

1. 基底（`data_basis`）：`L = 0..3` のスピンを含む LMs（`atomic_lms_basis`、`(m, spin)`）と JML（`atomic_jm_basis`、
   `(J, M)`）。スピンを含まない LM は上向きスピンの部分です。
2. 行列要素 `<b1|X_{lm}^{(s)}(k)|b2>`（`X = Q, G, T, M`、`s = 0, 1`、`l <= 7`）を Wigner–Eckart の定理で求めます
   （`multipie.util.util_atomic_multipole`）。Jz の保存で許される要素（`m1 + s1 = m2 + s2 + m`、`M1 = M2 + m`）だけを
   計算し、パリティや三角条件で禁止されるブロックは飛ばします。
3. `data_basis.atomic_rank_block`（上三角）の各 bra-ket ブロック `(L1, L2)` を、0 でなければ残します。
4. キー：`SphericalMultipoleType(X, l, s, k, "q")`。並びは `Q, G, T, M`、`s`、`k`、`l`。値：`m = l..-l` の
   `(2l+1, n1, n2)`。`"lm"` は `"lms"` の `s = 0`、上向きスピンの部分です。

### 6.2 点群の SAMB（`atomic_multipole_group.pkl`）

1. **トップの SAMB。** トップの調和関数 `T`（u-matrix `U_T`、`m = l..-l`）の各成分と、パリティが合う各原子多極子
   （極性の調和関数には `Q, T`、軸性には `G, M`）について、`sum_m U_T[m] X_{lm}^{(s)}(k)` を作り、Frobenius ノルムが
   1（対角ブロック）または `1/sqrt(2)`（非対角ブロック）になるよう正規化します。トップの群ごとに一度だけ計算します。
2. 成分どうしは直交しています（球テンソルが直交し、`U_T` の列が正規直交なため）。そのため、この正規化は MultiPie の
   Gram–Schmidt と同じ結果になります。非対角ブロックでは、同じ既約表現の Q と T（G と M）は内積の虚部だけが重なり
   ます（エルミート共役のブロックと合わせると直交します）。MultiPie も `X` ごとに直交化しています。
3. **点群。** 群の調和関数 `(X, l, Gamma, n)` の SAMB は、そのトップのキーのトップの SAMB の `c` 倍で
   （`data_harmonics_mapping`）、成分を重ねたものです。
4. **基底。** `"jml"` は JML の行列から作ります。`"lgs"` は LMs の行列 `M` から `U_bra^dagger M U_ket` で作ります。`U`
   は `data_basis.atomic_lgs_basis` の軌道（`atomic_orbital_basis`、基底関数の u-matrix をスピンに広げたもの）です。
   `"lg"` は `"lgs"` の `s = 0`、上向きスピンの部分に `sqrt(2)` を掛けたものです。
5. **キー。** `PGMultipoleType(X, l, Gamma, n, -1, s, k, "q")`。`n` は調和関数の多重度（MultiPie と違い振り直さない）。
   並びは `Q, G, T, M`、`s`、`k`、`l`、`Gamma`、`n`。
6. **値。** ブロックごとに `(行列 (dim, n1, n2), ex (dim))`。47 点群すべて（MultiPie は実の 37 群）。

---

## 7. 厳密計算と式の簡約

- 保存する値はすべて sympy です。保存する前に `util_symmetry_operation.simplify_value`
  （`expand(radsimp(expand(x)))`、有理数はそのまま）を通し、分母を有理化します（`1/sqrt(3)` ではなく `sqrt(3)/3`）。
- 一次従属かどうかは数値で判定し（相対許容誤差 `1e-8`）、残したベクトルを厳密に正規直交化します。
- pickle を読み込むと、sympy が一部の積を評価し直します（`(a+b)/2` → `a/2 + b/2`）。そのため「binary is up to date」
  のテストは、pickle に一度書き出して読み直した後で比べます。
