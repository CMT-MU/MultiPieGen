# MultiPie Ver.2 との違い（多極子のキーを使うデータ）

確認用のメモです（パッケージのドキュメントではありません）。MultiPieGen が作る、多極子のキー
（`PGMultipoleType` / `SphericalMultipoleType`）を使うデータについて、MultiPie Ver.2 の対応するデータとの違いを、
キーと内容に分けて網羅的にまとめます。数字はすべて `compare.py` で実際に比べた結果です（2026-11-11、develop
369c8bf 時点）。生成の方法やキーの規則そのものは `../generation_rules.md` にまとめています。

## まとめ

| MultiPieGen | MultiPie Ver.2 | キー | 内容 | 違いの原因 |
|---|---|---|---|---|
| 点群の調和関数（`data_harmonics_mapping`、rank ≤ 11、47 群） | `harmonics.pkl` | 同じ（10,296 多重項） | 32 多重項で符号だけ逆 | D6h の rank 4 `E2g` 多重度 2 を意図的に −1 倍 |
| `multipolar_harmonics.pkl`（s = 1..3、l ≤ 3、47 群） | `harmonics_multipole.pkl`（l ≤ 11） | 同じ（14,302 件、n と p も同じ） | 同じ（内部基底を変換すれば） | 内部基底の軌道と並び、記号が違う |
| `cluster_scalar_basis.pkl` など cluster 基底 3 本 | `harmonics_root_cluster.pkl` | キーの形式は同じ | 別物（直交化した基底 vs 生の調和関数） | クラスターと作り方が違う |
| `cluster_samb.pkl`（実の点群 37 + 空間群 230） | `cluster_samb.pkl` | site と bond_s は同じ、bond_a は 6663 中 5990 で違う | どの既約表現も同じ空間を張る、値は違う | 代表点、bond_a のキーの付け方 |
| `atomic_multipole.pkl` | `atomic_multipole.pkl` | 同じ | 同じ | なし |
| `atomic_multipole_group.pkl`（47 群） | `atomic_multipole_group.pkl`（37 群） | 同じ（37 群） | 896 行列で符号だけ逆 | 調和関数の符号の変更（D6h の rank 4 `E2g`） |
| 表現行列（`util_representation_matrix`） | `representation_matrix.pkl` | 同じ | 同じ（241 既約表現） | なし |
| MPG の active multipole（`group.pkl`） | MPG ごとの binary | 51 群で名前が違う | 名前を対応させれば同じ | 軌道名の変更 |
| 応答テンソル（`response_tensor.pkl`） | `Group.response_tensor_all` | 多極子の名前が違う | 名前以外は同じ（作成時に確認） | 軌道名の変更 |

共通する原因は次の 4 つです。

1. **調和関数の符号の変更。** D6h の rank 4 `E2g` 多重度 2（2 成分とも）を −1 倍しています（`data_top_harmonics.py`
   の注記）。D6h 系の部分群に引き継がれ、調和関数を使うデータ（atomic multipole group、応答テンソルなど）にも
   そのまま現れます。立方晶の `gv` も −1 倍していますが、これは軌道（基底関数）の定義の変更で、点群の調和関数
   には現れません。
2. **軌道名の変更。** f2→fc、f1→fs、f3x→fx、f3y→fy、ga→gc、gb→gs、gc→gc4、i→i3、l→l3、i0→ixy、j→j3
   （`gc` は旧 `ga` を指すようになった点に注意）。
3. **多重度 n の扱い。** MultiPieGen では、n は常に「その点群の調和関数 (X, l, Γ) の多重度」で、xyz の多項式を
   特定します。振り直しが必要なときは p を使います。今回比べた範囲では、multipolar harmonics と atomic multipole
   group で、MultiPie Ver.2 の振り直した n と結果的に一致していました。
4. **既約表現の並び。** `data_character_table` の並び（g が先、`'` が先）を正とします。MultiPie Ver.2 は名前の
   文字列で並べています。キーの集合は同じでも、辞書の並び順が違うことがあります。

---

## 1. 点群の調和関数

- **対象**：47 点群すべて、rank ≤ 11。こちらは `pg_harmonics(pg, 11)`（`data_harmonics_mapping` の `c × T`）、
  MultiPie Ver.2 は `harmonics.pkl` の `"harmonics"`。
- **キー**：(X, l, Γ, n) の集合は 47 群すべてで同じ（10,296 多重項）。複素の点群（`-c`）も含みます。
- **内容**：xyz の式を 3 点で数値にして比べました。10,264 多重項が同じで、**32 多重項が符号だけ逆**です。
  - 16 群（C3、C3h、C3i、C3v、C3v-1、C6、C6h、C6v、D3、D3-1、D3d、D3d-1、D3h、D3h-1、D6、D6h）の rank 4 の E 系列。
  - 例：D6h の Q₄(E2g) n=2、G₄(E2u) n=2。C3 の Q₄(E) n=3、G₄(E) n=3。D3h の Q₄(E') n=2、G₄(E'') n=2。
  - 原因は、D6h の rank 4 `E2g` 多重度 2 の符号の変更です。部分群ではこのトップの調和関数に対応する多重度
    （2 または 3）に現れます。
- **立方晶**：Oh 系の群では違いはありません。

## 2. multipolar harmonics

- **対象**：47 点群、s = 1, 2, 3、l ≤ 3。こちらは `multipolar_harmonics.pkl`、MultiPie Ver.2 は
  `harmonics_multipole.pkl`（l ≤ 11 まであるので l ≤ 3 だけ比較）。
- **キー**：(X, l, Γ, n, p, s, k, x) は 14,302 件すべてで同じ。各 (X, l, Γ, s, k, x) の件数も同じ。
- **内容**：係数 `mharm`（内部基底の係数）は、そのままでは比べられません。内部基底が違うからです。
  - MultiPie Ver.2：どの群でも tesseral の並び（s=2 は v, xy, xz, yz, u、s=3 は旧名 2, 1, bz, 3, 3x, 3y, az）。
    記号は `Q_v`、`Q_2` など。
  - MultiPieGen：立方晶系は `cubic_orbital_basis[s]`（du, dv, dyz, dxz, dxy、f3, fax, …）、三方・六方晶系は
    `hexagonal_orbital_basis[s]`（du, dxz, dyz, dxy, dv、faz, fs, fc, …）。記号は `Q_{u}`、`Q_{x}^{\alpha}` など
    （`orbital_latex_name`）。
  - こちらの係数を、軌道の u-matrix を使って tesseral の基底の係数に変換して比べると、**14,302 件すべてが符号まで
    一致**しました。つまり、違いは内部基底の選び方（軌道、並び、記号）だけです。
- **多極子の式 `m_ex`**：内部基底の記号が違うので、式の見た目は違います（中身は同じ関数）。
- **範囲**：MultiPie Ver.2 は l ≤ 11、こちらは l ≤ 3。

## 3. cluster 基底（cluster_scalar_basis、cluster_vector_basis、cluster_bond_basis）

MultiPie Ver.2 に同じものはありません。いちばん近いのは `harmonics_root_cluster.pkl`（各点群の調和関数を root
cluster の上で評価したもの）です。

| | MultiPie Ver.2 `harmonics_root_cluster` | MultiPieGen の cluster 基底 |
|---|---|---|
| クラスター | Oh と D6h に共通の 64 点の root cluster | Oh は 48 点、D6h は 24 点で別々 |
| 代表点 | MultiPie Ver.2 の代表点 | site0 = (1, 5, 6)、vector0 = (6, 1, 5) の単位ベクトル |
| 中身 | rank ≤ 11 の調和関数の値そのもの（直交化していない、一次従属なものも含む） | トップの群の (既約表現, 成分) ごとに低い rank から Gram–Schmidt し、一次従属なものを除いた正規直交基底 |
| 多重項の数（Oh） | monopole 120、dipole 178 | scalar 20（48 本）、vector x=q 60 と x=g 60（各 144 本）、bond x=q 20 と x=g 20（各 48 本） |
| 多重項の数（D6h） | monopole 192、dipole 286 | scalar 16（24 本）、vector 各 48（各 72 本）、bond 各 16（各 24 本） |
| 内部基底 | x = q だけ | vector と bond は x = q と x = g の両方 |

- **キー**：形式は MultiPie Ver.2 の site SAMB（scalar）、vector SAMB（vector）、bond_a（bond）に合わせています。
  scalar と vector（x = q）のキーは、すべて MultiPie Ver.2 の root cluster のキーの中にあります。
- **bond 基底**：MultiPie Ver.2 は「最初に重なる vector SAMB」からキーを決めますが、こちらはベクトル基底と
  ボンドベクトルの内積を、rank → k → 多重度の順に Gram–Schmidt して作ります（ユーザー指定の規則）。
- **値**：クラスターも作り方も違うので、値は比べられません。

## 4. cluster_samb

- **対象**：実の点群 37 と空間群 230 の、各 Wyckoff の site（1889 個）、bond_s と bond_a（各 6663 個）。
  Wyckoff の一覧は MultiPie Ver.2 と同じです。MultiPie Ver.2 の `"vector"` は作っていません。
- **張る空間**：site、bond_s、bond_a のすべてで、**既約表現ごとに張る空間が MultiPie Ver.2 と同じ**です。
- **site と bond_s**：キー (X, l, Γ, n, p, x) が全件同じです。
- **bond_a**：6663 個のうち 673 個はキーが同じです。違う 5990 個の内訳は次のとおりです。
  - **p だけ違う（4866 個）**：こちらはトップの群のキーの p を引き継ぎます（例：Oh の `3a@1o` は、こちらが
    Q₁(T1u) p=1、MultiPie Ver.2 は p=−1）。MultiPie Ver.2 は Wyckoff ごとに p を振り直します。
  - **n が違う（93 個）**：同じ (X, l, Γ) の中で、どの調和関数の多重度をキーにするかが違います（例：C4 の
    `4c@1a` は、こちらが Q₂(B) n=1、MultiPie Ver.2 は n=2）。
  - **(X, l, Γ) が違う（1031 個、168 群、うち点群 15 群）**：同じ空間に付けるキーが違います。こちらのほうが低い
    rank のキーになる傾向があります（例：C4v の `4f@1a` は、こちらが G₁(E)、MultiPie Ver.2 は Q₁(E)）。
  - 原因は、bond_a のキーの付け方です。こちらは `cluster_bond_basis`（x = q）を一般点に制限し、[Γ, l, k, Q/G, n, p]
    の順に直交化して、トップのキーをそのまま使います。MultiPie Ver.2 は monopole の順に直交化してから、最初に
    重なる vector SAMB のキーに付け替えます。
- **値**：一般点の代表点（site0）が違うので、各キーの値は MultiPie Ver.2 と一致しません（張る空間は同じ）。
- **空間群の対応**：こちらは空間群の対称操作を点群の対応（`op_mapping_d6h_oh_cluster[PG]`）でクラスターの点に
  対応させます。MultiPie Ver.2 も点操作のタグで対応させているので、考え方は同じです。

## 5. atomic multipole

- **atomic_multipole.pkl**：キー、並び、値のすべてが MultiPie Ver.2 と同じです（380 件）。
- **atomic_multipole_group.pkl**（実の点群 37、jml/lgs/lg）：
  - キー (X, l, Γ, n, p, s, k, x) は同じです。n はすべて一致しました（61,653 行列）。
  - **896 行列で符号だけ逆**です。どれも六方・三方晶系の群の rank 4 の E 系列で、調和関数の式の符号も同時に逆に
    なっています（1 と同じ原因）。それ以外の 60,757 行列は同じです。
  - こちらは複素の点群（`-c`）10 群も作っています（MultiPie Ver.2 は 37 群）。
  - 正規化は MultiPie Ver.2 と同じです（成分ごとに Frobenius ノルム 1、非対角ブロックは 1/√2）。

## 6. 表現行列

- 47 点群の 241 既約表現すべてで、対称操作の並びも行列も MultiPie Ver.2 と同じです。

## 7. MPG の active multipole と応答テンソル

- **active multipole**（`group.pkl` の MPG）：157 群のうち 106 群は同じです。51 群は軌道名が違うだけで（例：
  `Qf1` → `Qfs`、`Qga` → `Qgc`）、名前を対応させれば同じです。
- **応答テンソル**（`response_tensor.pkl`）：MultiPie Ver.2 の `Group.response_tensor_all` と、多極子の名前以外は
  同じです。これは `response_tensor.pkl` を作ったときに確かめたもので（PR の検証）、`compare.py` では比べていません。

---

## 再現方法

```bash
python others/claude/multipie_v2_comparison/compare.py              # 全項目（約 15 分）
python others/claude/multipie_v2_comparison/compare.py harmonics    # 1 項目だけ
```

- MultiPie Ver.2 は `../MultiPie`（または環境変数 `MULTIPIE_DIR`）から読みます。結果は `result/<項目>.json` に
  書き出します。
- **メモリに注意。** MultiPie Ver.2 の `harmonics_multipole.pkl` は、読み込むだけで約 8 GB、約 4 分かかります。
  `compare.py` は項目ごとに別のプロセスで、1 つずつ順番に実行します（並列にしない）。`mp_cache` で一度だけ読み込み、
  l ≤ 3 の部分を数値にして `cache/`（git に入れない）に書き出し、`multipolar` はそれを読みます。ほかの項目は
  1 GB 以下です。
- 項目：`harmonics`、`mp_cache`、`multipolar`、`cluster_basis`、`cluster_samb`、`atomic`、`repmat`、`active`。
