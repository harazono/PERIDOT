# Search Primer

DNA配列から高頻度で出現するLR-tuple（Left-Right tuple）を効率的に検出するRustプログラム。

## アルゴリズム概要

このプログラムは、DNA配列データから種特異的なプライマー設計に有用な高頻度LR-tupleを検出するために、Counting Bloom Filter（CBF）とハッシュテーブルを組み合わせた3段階のアルゴリズムを実装している。

### 1. データ構造と定数

- **L_LEN**: 32 (左側セグメントの長さ)
- **R_LEN**: 32 (右側セグメントの長さ)
- **CHUNK_MAX**: 200 (LR-tuple間の最大距離)
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

### 5. 出力形式

- **テキスト形式**: デコードされたDNA配列
- **バイナリ形式**: u128値の16バイト表現
- **カウントのみ**: 閾値を超えるLR-tupleの総数

### 6. 使用方法

```bash
./search_primer input.fasta -o output.txt -t 8 -a 1000 -m 0
```

- `-t`: スレッド数（デフォルト: 8）
- `-a`: 閾値（デフォルト: 1000）
- `-m`: L-Rセグメント間のマージン（デフォルト: 0）
- `-b`: バイナリ出力
- `-r`: カウントのみ出力

### 7. 計算量

- **時間計算量**: O(n × L × R × t) （n: 配列数, L,R: ウィンドウ数, t: スレッド数）
- **空間計算量**: O(2^30) （Bloom Filter） + O(候補数) （HashSet/HashMap）

このアルゴリズムにより、大規模なDNA配列データセットから効率的に高頻度LR-tupleを抽出し、種特異的プライマー設計の基盤データを提供する。
