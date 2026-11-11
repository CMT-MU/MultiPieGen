# MultiPie Ver.2 からの変更点（多極子のデータ）

多極子のキーを使うデータについて、MultiPieGen が MultiPie Ver.2 と違うことをしている点（規則、アルゴリズム、
その結果の違い）をまとめたものです。MultiPieGen のアルゴリズム
そのものは [多極子データの生成](generation_rules_ja.md) にあります。下の数字は 2026-11-11 に MultiPie Ver.2 と
比べた結果です（develop 369c8bf の時点）。[英語版](changes_from_multipie_v2.md) もあります。

目次：

1. 規則（キー、既約表現、軌道名、符号）
2. 点群の調和関数
3. multipolar harmonics
4. cluster 基底（新規）
5. cluster SAMB
6. 原子多極子
7. cluster SAMB で使う対称操作の対応
8. そのほかのデータ（表現行列、active multipole、応答テンソル）

---

## 1. 規則

| 項目 | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| 多重度 `n` | データごとに振り直す（multipolar harmonics、atomic multipole group、cluster SAMB） | 常にその点群の調和関数 `(X, l, Gamma)` の多重度で、xyz の多項式を特定する。重複は `p` で番号付け |
| 既約表現の並び | 名前の文字列で並べる | `data_character_table` の並び（`g` が `u` より先、`'` が `''` より先） |
| 軌道名 | f2, f1, f3x, f3y, ga, gb, gc, i, l, i0, j | fc, fs, fx, fy, gc, gs, gc4, i3, l3, ixy, j3（`gc` は旧 `ga` を指す） |
| 軌道の定義 | `data_definition.py` | `data_basis.py`。各軌道は調和関数のキーを指し、式はトップの調和関数から取る |
| D6h の rank 4 `E2g` 多重度 2 の符号 | MultiPie のまま | −1 倍（2 成分とも） |
| 軌道 `gv` の符号 | MultiPie のまま | −1 倍（`data_response_tensor`） |
| 値の中の調和関数の式 | 多極子のキーを持つ binary が調和関数のデカルト座標の式 `ex` を持つ | 持たない。キーの `(X, l, Gamma, n)` で `harmonics.pkl` から引く（値は、multipolar harmonics が `(internal_coefficient, multipole_expression)`、atomic multipole group が `matrix`、cluster SAMB が `samb`） |

- `E2g` の符号の変更はトップの調和関数の変更です。D6h 系の部分群（六方・三方晶系 16 群の rank 4 の E の多重項）と、
  調和関数から作るすべてのデータに引き継がれます。
- `gv` の符号の変更は軌道（と応答テンソルの多極子）の変更です。Oh のトップの調和関数は変わらず、Oh 系の調和関数は
  MultiPie Ver.2 と一致します。
- これまで比べた範囲（multipolar harmonics、atomic multipole group）では、調和関数の多重度 `n` は MultiPie Ver.2 の
  振り直した `n` と結果的に一致していました。規則の違いが効くのは、cluster SAMB の bond_a です。

## 2. 点群の調和関数

- Oh と D6h について一度だけ定義し（`data_top_harmonics_info`）、`data_harmonics_mapping` で 47 点群に対応させます
  （トップの調和関数の `c = ±1` 倍）。MultiPie Ver.2 は各点群の調和関数を `harmonics.pkl` に持っています。
- 結果：キー `(X, l, Gamma, n)` は 47 群すべてで同じです（rank ≤ 11）。32 多重項で符号が逆で、すべて `E2g` の符号の
  変更（1 節）によるものです。

## 3. multipolar harmonics

| 項目 | MultiPie Ver.2（`harmonics_multipole.pkl`） | MultiPieGen（`multipolar_harmonics.pkl`） |
|---|---|---|
| 内部基底 | どの群でも tesseral の基底（`Q_v`、`Q_xy`、…、f は旧名） | Oh 系は立方晶の軌道（`cubic_orbital_basis[s]`）、D6h 系は六方晶の軌道（`hexagonal_orbital_basis[s]`） |
| 内部の記号 | `Q_v`、`Q_2` など | `Q + orbital_latex_name`（`Q_{u}`、`Q_{x}^{\alpha}` など） |
| `n` | `(x, l, X, Gamma, k)` ごとに振り直す | 調和関数の多重度。重複は `p` |
| ランク | l ≤ 11 | l ≤ 3 |

- 結果：キーは同じです（14,302 件、l ≤ 3）。内部基底の係数を tesseral の基底に変換すると、すべての件が一致します。
  違いは内部基底の選び方（軌道、並び、記号）だけです。

## 4. cluster 基底（新規）

`cluster_scalar_basis.pkl`、`cluster_vector_basis.pkl`、`cluster_bond_basis.pkl` に当たるものは MultiPie Ver.2 には
ありません。いちばん近いのは `harmonics_root_cluster.pkl` です。

| 項目 | MultiPie Ver.2 `harmonics_root_cluster` | MultiPieGen の cluster 基底 |
|---|---|---|
| クラスター | Oh と D6h に共通の 64 点の root cluster | Oh（48 点）と D6h（24 点）を別々に持つ |
| 代表点 | MultiPie Ver.2 の点 | site0 = (1, 5, 6)、vector0 = (6, 1, 5) の単位ベクトル |
| 中身 | 各点群の rank 11 までの調和関数の値（直交化していない、一次従属なものも含む） | トップの群の (既約表現, 成分) ごとに低いランクから Gram–Schmidt し、一次従属なものを除いた正規直交基底（正則表現） |
| 内部の種類 | x = q | x = q と x = g（vector と bond の基底） |

- キーの形式は、MultiPie Ver.2 の site SAMB（scalar）、vector SAMB（vector）、bond_a（bond）に合わせています。
- bond の基底は、ベクトル基底とボンドベクトルの内積を取り、ランク、k、多重度の順に Gram–Schmidt して作ります。
  MultiPie Ver.2 は反対称ボンドの関数ごとに、最初に重なる vector SAMB のキーを付けますが、その規則は使いません。

## 5. cluster SAMB

| 項目 | MultiPie Ver.2（`cluster_samb.pkl`） | MultiPieGen（`cluster_samb.pkl`） |
|---|---|---|
| 元にするもの | 各点群の調和関数を 64 点の root cluster で評価した値 | Oh/D6h の cluster 基底を、その群の操作（`op_mapping_d6h_oh_cluster` の coset 0）に制限したもの |
| 特殊点 | Wyckoff の mapping で足し合わせてから直交化 | 同じ（一般点の SAMB を Wyckoff の mapping で平均。bond_a は符号つき） |
| site、bond_s | [Gamma, l, Q/G, n, p] の順に Gram–Schmidt | 既約表現ごとに [Gamma, l, n, p] の順に Gram–Schmidt。キーは同じ |
| bond_a | monopole の順に直交化し、最初に重なる vector SAMB のキーを付ける。`p` は Wyckoff ごとに振り直す | `cluster_bond_basis`（x = q）を制限し、[Gamma, l, k, Q/G, n, p] の順に Gram–Schmidt。トップの群のキー（`p` も）を引き継ぐ |
| vector SAMB | 作る | 作らない |
| 群 | 実の点群と空間群 | 実の点群と空間群（複素の点群は含まない） |

- 結果：すべての Wyckoff のサイトとボンドで、既約表現ごとに張る空間は MultiPie Ver.2 と同じです（サイト 1889、
  ボンド 6663）。site と bond_s はキーも同じです。
- bond_a は、6663 個のボンド Wyckoff のうち 5990 個でキーが違います。
  - p だけ違う（4866）：トップのキーの `p` を引き継ぐため（例：Oh の `3a@1o` は、こちらが Q₁(T1u) p = 1、MultiPie
    Ver.2 は p = −1）。
  - n が違う（93）：同じ (X, l, Gamma) の別の多重度を使う（例：C4 の `4c@1a` は、こちらが Q₂(B) n = 1、MultiPie Ver.2
    は n = 2）。
  - (X, l, Gamma) が違う（1031、168 群、うち点群 15 群）：同じ空間に別のキーを付ける。こちらのほうが低いランクに
    なることが多い（例：C4v の `4f@1a` は、こちらが G₁(E)、MultiPie Ver.2 は Q₁(E)）。
- 値は一致しません。一般点の代表点（site0）が違うためです。

## 6. 原子多極子

- `atomic_multipole.pkl`：MultiPie Ver.2 と同じです（キー、並び、値）。違うのは計算の仕方だけで、Jz の保存で許される
  行列要素だけを計算します（約 10 秒）。
- `atomic_multipole_group.pkl`：

| 項目 | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| 作り方 | 点群ごとに、その調和関数で原子多極子を射影し、(X, 既約表現) ごとに Gram–Schmidt して、sqrt(dim) と非対角ブロックでは 1/sqrt(2) を掛ける | トップの調和関数の成分ごとに正規化した SAMB を、トップの群ごとに一度だけ計算し、各点群ではその `c` 倍を使う（正規化は同じ） |
| `n` | 振り直す | 調和関数の多重度 |
| 群 | 実の点群 37 | 47 点群（複素の `-c` を含む） |
| Lg の基底 | MultiPie Ver.2 の lg の軌道（旧名） | `data_basis.atomic_lgs_basis` |

- 結果（実の点群 37）：キーは同じです。896 行列で符号が逆（`E2g` の符号の変更）で、残りの 60,757 行列は同じです。

## 7. cluster SAMB で使う対称操作の対応

| 項目 | MultiPie Ver.2 | MultiPieGen |
|---|---|---|
| 点群の `symmetry_operation["mapping"]` | 積表から作る右剰余類 `H g` | 親の群をたどる左剰余類 `g H`（`data_op_mapping.op_mapping_d6h_oh_cluster`） |
| 空間群 | なし | その点群の対応（空間群の操作は、点群と同じ点操作を同じ順で持つ） |
| 空間群から親の空間群へ | なし | `op_mapping_parent_sg`（軸の変換を含む。クラスターには使わない） |

## 8. そのほかのデータ

- 表現行列：47 点群の 241 既約表現すべてで同じです。
- MPG の active multipole（`group.pkl`）：軌道名を対応させれば同じです（157 群のうち 51 群で軌道名が変わっている。
  例：`Qf1` → `Qfs`）。
- 応答テンソル（`response_tensor.pkl`）：MultiPie Ver.2 の `Group.response_tensor_all` と、多極子の名前以外は同じです
  （binary を作ったときに確認したもので、2026-11-11 の比較には含めていません）。
