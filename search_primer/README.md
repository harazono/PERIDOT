# Search Primer

DNA配列から高頻度で出現するLR-tuple（Left-Right tuple）を効率的に検出するRustプログラム。

## ビルド方法

### 必要な環境

- Rust 1.70以上
- Cargo

### ビルド手順

```bash
# プロジェクトディレクトリに移動
cd search_primer

# 全バイナリをビルド
cargo build --release

# 特定のバイナリのみビルド
cargo build --release --bin search_primer_pipeline
```

### 利用可能なバイナリ

- `search_primer` - オリジナル版（互換性維持）
- `search_primer_pipeline` - 統合パイプライン版（推奨）
- `build_cbf` - CBF構築専用
- `extract_candidates` - 候補抽出専用
- `remove_false_positives` - 偽陽性除去専用
- `primer3_caller` - Primer3統合ツール（NEW）

## 使い方

### 1. 設定ファイル作成（オプション）

```toml
# config.toml
[sequence]
l_len = 30
r_len = 30
chunk_max = 200
margin_size = 50

[bloom_filter]
bloomfilter_table_size = 1073741824  # 2^30
hashset_size = 536870912             # 2^29
```

### 2. 基本的な使用方法

#### 統合パイプライン（推奨）

```bash
# デフォルト設定で実行
./target/release/search_primer_pipeline input.fasta

# カスタム設定で実行
./target/release/search_primer_pipeline input.fasta -c config.toml -t 16 -a 500

# CBFを保存して再利用
./target/release/search_primer_pipeline input.fasta --save-cbf -o results.txt

# 既存CBFを使用
./target/release/search_primer_pipeline input.fasta --cbf input.fasta.cbf -a 1000
```

#### 段階的実行

```bash
# 1. CBF構築
./target/release/build_cbf input.fasta -o input.cbf -c config.toml -t 8

# 2. 候補抽出
./target/release/extract_candidates input.cbf input.fasta -a 1000 -t 8

# 3. 偽陽性除去
./target/release/remove_false_positives input.cbf input.fasta -a 1000 -t 8 -o final_results.txt
```

#### Primer3統合ツール

```bash
# LR-tupleからプライマー設計
./target/release/primer3_caller lr_tuples.bin -c primer3_config.txt -o primers.txt

# 設定ファイル指定
./target/release/primer3_caller lr_tuples.bin -c primer3_config.txt --search-config config.toml -t 8

# 一時ファイルディレクトリ指定
./target/release/primer3_caller lr_tuples.bin -c primer3_config.txt -m /tmp/primer3_work
```

### 3. コマンドラインオプション

#### 共通オプション

- `-c, --config CONFIG` - 設定ファイルパス
- `-t, --threads THREADS` - スレッド数（デフォルト: 8）
- `-a, --threshold THRESHOLD` - 閾値（デフォルト: 1000）
- `-o, --output OUTPUT` - 出力ファイルパス
- `-h, --help` - ヘルプ表示

#### 出力形式オプション

- `-b, --binary` - バイナリ形式で出力
- `-r, --only-num` - カウントのみ出力

#### パイプライン専用オプション

- `--cbf CBF_FILE` - 既存CBFファイルを使用
- `--save-cbf` - CBFファイルを保存

#### Primer3ツール専用オプション

- `-c, --config CONFIG` - Primer3設定ファイル（必須）
- `--search-config CONFIG` - Search Primer設定ファイル
- `-m, --tmpfile PREFIX` - 一時ファイル名プレフィックス

### 4. 出力形式

#### テキスト形式（デフォルト）

```
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
...
```

#### カウントのみ（-r オプション）

```
lr_tuple count: 1234 threshold: 1000 input file: input.fasta
```

#### バイナリ形式（-b オプション）

16バイトのu128値が連続して出力されます。

### 5. パフォーマンス最適化

#### メモリ使用量

- CBFファイル: 約2GB
- 作業メモリ: 入力ファイルサイズの2-3倍

#### 推奨設定

```bash
# 大規模データセット用
./target/release/search_primer_pipeline large_input.fasta -t 32 -a 500 --save-cbf

# 高精度検索用
./target/release/search_primer_pipeline input.fasta -a 2000 -c high_precision_config.toml
```

### 6. トラブルシューティング

#### CBF互換性エラー

```
Error: CBF file incompatible with current config
```

→ 設定が変更された場合、CBFを再構築してください

#### メモリ不足

```
Error: failed to allocate memory
```

→ スレッド数を減らすか、より小さな閾値を使用してください

#### ファイルが見つからない

```
Error: No such file or directory
```

→ ファイルパスを確認してください

## Primer3統合機能

### 機能概要
- LR-tupleからPrimer3を使用したプライマー設計
- バージョン情報の自動記録（`primer3_info.log`）
- エラーハンドリングとログ出力の改善
- マルチスレッド対応

### 出力ファイル
- **プライマー結果**: 指定した出力ファイル
- **バージョンログ**: `primer3_info.log`（実行ディレクトリ）
- **エラーログ**: 標準エラー出力

### Primer3設定ファイル例
```
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=30
PRIMER_MIN_SIZE=16
PRIMER_MAX_SIZE=32
PRIMER_PRODUCT_SIZE_RANGE=101-200 201-301
PRIMER_OPT_TM=66.0
PRIMER_MAX_TM=72.0
PRIMER_EXPLAIN_FLAG=1
```

## アルゴリズム概要

このプログラムは、DNA配列データから種特異的なプライマー設計に有用な高頻度LR-tupleを検出するために、Counting Bloom Filter（CBF）とハッシュテーブルを組み合わせた3段階のアルゴリズムを実装している。

### 1. データ構造と定数

- **L_LEN**: 30 (左側セグメントの長さ)
- **R_LEN**: 30 (右側セグメントの長さ)
- **CHUNK_MAX**: 200 (LR-tuple間の最大距離)
- **MARGIN_SIZE**: 50 (L-Rセグメント間のマージン)
- **BLOOMFILTER_TABLE_SIZE**: 2^30 (Bloom Filterのサイズ)
- **HASHSET_SIZE**: 2^29 (HashSetの容量)

### 2. 3段階処理アルゴリズム

#### 第1段階: Counting Bloom Filter構築

```
for each DNA sequence:
    for each L-window position:
        if L-segment has repeats: skip
        for each R-window position (with margin):
            if R-segment has repeats: skip
            if distance > CHUNK_MAX: break
            
            lr_tuple = encode(L-segment + R-segment) as u128
            hash_indices = SHA256_hash(lr_tuple) % table_size
            for each hash_index:
                CBF[hash_index] += 1
```

#### 第2段階: 高頻度LR-tuple候補抽出

```
for each DNA sequence:
    for each L-window position:
        if L-segment has repeats: skip
        for each R-window position (with margin):
            if R-segment has repeats: skip
            
            lr_tuple = encode(L-segment + R-segment) as u128
            hash_indices = SHA256_hash(lr_tuple)
            estimated_count = min(CBF[hash_indices])
            
            if estimated_count >= threshold:
                candidate_set.insert(lr_tuple)
```

#### 第3段階: 正確なカウントと偽陽性除去

```
for each DNA sequence:
    for each L-window position:
        if L-segment has repeats: skip
        for each R-window position (with margin):
            if R-segment has repeats: skip
            
            lr_tuple = encode(L-segment + R-segment) as u128
            if lr_tuple in candidate_set:
                exact_count[lr_tuple] += 1

final_results = filter(exact_count, count >= threshold)
```

### 3. 主要な最適化技術

#### マルチスレッド処理

- 各段階でDNA配列をチャンクに分割し、並列処理
- スレッド間でのCBFマージとHashSet統合

#### メモリ効率化

- DNA配列のu128エンコーディング（2bit/base）
- Counting Bloom Filterによる偽陽性を許容した高速フィルタリング
- 段階的な候補絞り込みによるメモリ使用量削減

#### リピート配列の除外

- 各セグメントでリピート配列を検出し、スキップ
- プライマー設計に不適切な配列の事前除去

### 4. ハッシュ関数

SHA256を使用して8つの独立したハッシュ値を生成：

```rust
hash_indices[i] = SHA256(lr_tuple)[i*4:(i+1)*4] % table_size
```

### 5. CBFファイル形式

```
[Header: 32 bytes]
- Magic number: 8 bytes ("CBFV0001")
- Config hash: 8 bytes (設定変更検出用)
- Table size: 8 bytes
- Reserved: 8 bytes

[Data: table_size * 2 bytes]
- CBF array (u16 values)
```

### 6. 計算量

- **時間計算量**: O(n × L × R × t) （n: 配列数, L,R: ウィンドウ数, t: スレッド数）
- **空間計算量**: O(2^30) （Bloom Filter） + O(候補数) （HashSet/HashMap）

このアルゴリズムにより、大規模なDNA配列データセットから効率的に高頻度LR-tupleを抽出し、種特異的プライマー設計の基盤データを提供する。
