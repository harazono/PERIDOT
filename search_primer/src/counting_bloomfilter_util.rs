use crate::config::Config;
use crate::sequence_encoder_util::DnaSequence;
use sha2::Digest;
use sha2::Sha256;
use std::collections::HashMap;
use std::collections::HashSet;
use std::time::Instant;

// 後方互換性のため定数を再エクスポート
pub const L_LEN: usize = 30;
pub const R_LEN: usize = 30;
pub const BLOOMFILTER_TABLE_SIZE: usize = 1 << 30;
pub const HASHSET_SIZE: usize = 1 << 29;

//全てのL, Rと、hash値を出力する
//部分配列のdecoderを書き、テストする
pub fn build_counting_bloom_filter(
    sequences: &Vec<DnaSequence>,
    start_idx: usize,
    end_idx: usize,
    config: &Config,
    thread_id: usize,
) -> Vec<u16> {
    let mut l_window_start_idx: usize;
    let mut l_window_end_idx: usize;
    let mut r_window_start_idx: usize;
    let mut r_window_end_idx: usize;

    let mut loop_cnt: usize = 0;
    eprintln!(
        "Allocating Vec<u16> where BLOOMFILTER_TABLE_SIZE = {}",
        config.bloomfilter_table_size
    );
    let mut ret_array: Vec<u16> = vec![0u16; config.bloomfilter_table_size];
    eprintln!("Filling Vec<u16; {}> with 0", config.bloomfilter_table_size);
    eprintln!("finish allocating");

    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();

    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize = 0;
        loop_cnt += 1;
        l_window_start_idx = 0;
        if current_sequence.len() < config.l_len || current_sequence.len() < config.r_len {
            continue 'each_read;
        }
        'each_l_window: loop {
            l_window_end_idx = l_window_start_idx + config.l_len;
            if l_window_end_idx >= current_sequence.len() + 1 {
                let end = start_time.elapsed();
                eprintln!("1st loop[{:02}]({:05}-{:05},length is {})\t{:05?}({:.4}%)\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",
                    thread_id,
                    start_idx,
                    end_idx,
                    end_idx - start_idx,
                    loop_cnt,
                    loop_cnt as f64 / (end_idx - start_idx) as f64 * 100f64,
                    current_sequence.len(),
                    end.as_secs() - previous_time.as_secs(),
                    end.subsec_millis() - previous_time.subsec_millis(),
                    add_bloom_filter_cnt,
                    l_window_cnt
                );
                previous_time = end;
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_repeat_bool, l_has_repeat_offset) =
                current_sequence.has_repeat(l_window_start_idx, l_window_end_idx);
            if l_has_repeat_bool {
                l_window_start_idx += l_has_repeat_offset + 1;
                continue 'each_l_window;
            }
            r_window_start_idx = l_window_end_idx + config.margin_size;
            'each_r_window: loop {
                r_window_end_idx = r_window_start_idx + config.r_len;
                if r_window_end_idx > current_sequence.len() {
                    break 'each_r_window;
                }
                if r_window_end_idx - l_window_start_idx > config.chunk_max - config.r_len + config.margin_size {
                    break 'each_r_window;
                }
                let (r_has_repeat_bool, r_has_repeat_offset) =
                    current_sequence.has_repeat(r_window_start_idx, r_window_end_idx);
                if r_has_repeat_bool {
                    r_window_start_idx += r_has_repeat_offset + 1;
                    continue 'each_r_window;
                }
                add_bloom_filter_cnt += 1;
                let lmr_string: u128 = current_sequence.subsequence_as_u128(vec![
                    [l_window_start_idx, l_window_end_idx],
                    [r_window_start_idx, r_window_end_idx],
                ]);
                let table_indice: [u32; 8] = hash_from_u128(lmr_string, config.bloomfilter_table_size);
                for &i in table_indice.iter() {
                    let idx: usize = i as usize;
                    if ret_array[idx] == u16::MAX {
                        eprintln!("index {} reaches u16::MAX", idx);
                    } else {
                        ret_array[idx] += 1;
                    }
                }
                r_window_start_idx += 1;
            }
            l_window_start_idx += 1;
        }
        let end: std::time::Duration = start_time.elapsed();
        eprintln!("1st loop[{:02}]({:05}-{:05},length is {})\t{:05?}({:.4}%)\tlength: {}\tsec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}",
            thread_id,
            start_idx,
            end_idx,
            end_idx - start_idx,
            loop_cnt,
            loop_cnt as f64 / (end_idx - start_idx) as f64 * 100f64,
            current_sequence.len(),
            end.as_secs() - previous_time.as_secs(),
            end.subsec_millis() - previous_time.subsec_millis(),
            add_bloom_filter_cnt,
            l_window_cnt
        );
        previous_time = end;
    }
    return ret_array;
}

pub fn hash_from_u128(source: u128, table_size: usize) -> [u32; 8] {
    let mut ret_val: [u32; 8] = [0; 8];
    let mut hasher = Sha256::new();
    let u8_array: [u8; 16] = source.to_le_bytes();
    hasher.update(u8_array);
    let result = hasher.finalize();
    for i in 0..8 {
        ret_val[i] = ((result[i * 4] as u32) << 24)
            | ((result[i * 4 + 1] as u32) << 16)
            | ((result[i * 4 + 2] as u32) << 8)
            | (result[i * 4 + 3] as u32);
        ret_val[i] %= table_size as u32;
    }
    ret_val
}

pub fn count_occurence_from_counting_bloomfilter_table(
    counting_bloomfilter_table: &Vec<u16>,
    indice: [u32; 8],
) -> u16 {
    let mut retval: u16 = u16::MAX;
    for index in indice {
        if counting_bloomfilter_table[index as usize] < retval {
            retval = counting_bloomfilter_table[index as usize];
        }
    }
    return retval;
}

pub fn number_of_high_occurence_lr_tuple(
    source_table: &Vec<u16>,
    sequences: &Vec<DnaSequence>,
    start_idx: usize,
    end_idx: usize,
    threshold: u16,
    config: &Config,
    thread_id: usize,
) -> HashSet<u128> {
    let mut ret_table: HashSet<u128> = HashSet::with_capacity(config.hashset_size);
    let mut l_window_start_idx: usize;
    let mut l_window_end_idx: usize;
    let mut r_window_start_idx: usize;
    let mut r_window_end_idx: usize;

    let start: Instant = Instant::now();
    let mut previous_time: std::time::Duration = start.elapsed();
    let mut loop_cnt: usize = 0;
    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize = 0;
        loop_cnt += 1;
        l_window_start_idx = 0;
        if current_sequence.len() <= config.l_len || current_sequence.len() <= config.r_len {
            continue 'each_read;
        }
        'each_l_window: loop {
            l_window_end_idx = l_window_start_idx + config.l_len;
            if l_window_end_idx > current_sequence.len() {
                let end: std::time::Duration = start.elapsed();
                eprintln!("2nd loop[{:02}]({:05}-{:05},length is {})\t{:05?}({:.4}%)\tlength: {}\tsec: {}.{:03}\t high occurence LR-tuple: {}\tl_window_cnt: {}",
                    thread_id,
                    start_idx,
                    end_idx,
                    end_idx - start_idx,
                    loop_cnt,
                    loop_cnt as f64 / (end_idx - start_idx) as f64 * 100f64,
                    current_sequence.len(),
                    end.as_secs() - previous_time.as_secs(),
                    end.subsec_millis() - previous_time.subsec_millis(),
                    add_bloom_filter_cnt,
                    l_window_cnt
                );
                previous_time = end;
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let (l_has_repeat_bool, l_has_repeat_offset) =
                current_sequence.has_repeat(l_window_start_idx, l_window_end_idx);
            if l_has_repeat_bool {
                l_window_start_idx += l_has_repeat_offset + 1;
                continue 'each_l_window;
            }
            r_window_start_idx = l_window_end_idx + config.margin_size;
            'each_r_window: loop {
                r_window_end_idx = r_window_start_idx + config.r_len;
                if r_window_end_idx >= current_sequence.len() + 1 {
                    break 'each_r_window;
                }
                if r_window_end_idx - l_window_start_idx > config.chunk_max - config.r_len + config.margin_size {
                    break 'each_r_window;
                }
                let (r_has_repeat_bool, r_has_repeat_offset) =
                    current_sequence.has_repeat(r_window_start_idx, r_window_end_idx);
                if r_has_repeat_bool {
                    r_window_start_idx += r_has_repeat_offset + 1;
                    continue 'each_r_window;
                }
                add_bloom_filter_cnt += 1;
                let lmr_string: u128 = current_sequence.subsequence_as_u128(vec![
                    [l_window_start_idx, l_window_end_idx],
                    [r_window_start_idx, r_window_end_idx],
                ]);
                let table_indice: [u32; 8] = hash_from_u128(lmr_string, config.bloomfilter_table_size);
                let occurence: u16 =
                    count_occurence_from_counting_bloomfilter_table(source_table, table_indice);
                if occurence >= threshold {
                    if ret_table.len() >= (ret_table.capacity() as f64 * 0.9) as usize {
                        eprintln!("2nd loop[{:02}] attempts to reallocate HashSet.", thread_id);
                        break 'each_read;
                    }
                    ret_table.insert(lmr_string);
                }

                r_window_start_idx += 1;
            }
            l_window_start_idx += 1;
        }
        let end: std::time::Duration = start.elapsed();
        eprintln!("2nd loop[{:02}]({:05}-{:05},length is {})\t{:05?}({:.4}%)\tlength: {}\tsec: {}.{:03}\t high occurence LR-tuple: {}\tl_window_cnt: {}",
            thread_id,
            start_idx,
            end_idx,
            end_idx - start_idx,
            loop_cnt,
            loop_cnt as f64 / (end_idx - start_idx) as f64 * 100f64,
            current_sequence.len(),
            end.as_secs() - previous_time.as_secs(),
            end.subsec_millis() - previous_time.subsec_millis(),
            add_bloom_filter_cnt,
            l_window_cnt
        );
        previous_time = end;
    }
    return ret_table;
}

pub fn count_lr_tuple_with_hashtable(
    sequences: &Vec<DnaSequence>,
    start_idx: usize,
    end_idx: usize,
    high_occurence_lr_tuple: &HashSet<u128>,
    config: &Config,
    thread_id: usize,
) -> HashMap<u128, u16> {
    let hash_size_to_allocate: usize = (high_occurence_lr_tuple.len() as f64 * 1.2) as usize;
    eprintln!(
        "thread [{:02}] Allocating HashMap<u128, u16> where hash_size_to_allocate = {}",
        thread_id, hash_size_to_allocate
    );
    let mut lr_tuple_hashmap: HashMap<u128, u16> = HashMap::with_capacity(hash_size_to_allocate);
    eprintln!(
        "thread [{:02}] finish Allocating HashMap<u128, u16> where hash_size_to_allocate = {}",
        thread_id, hash_size_to_allocate
    );
    let mut l_window_start_idx: usize;
    let mut l_window_end_idx: usize;
    let mut r_window_start_idx: usize;
    let mut r_window_end_idx: usize;

    let start_time: Instant = Instant::now();
    let mut previous_time: std::time::Duration = start_time.elapsed();
    let mut loop_cnt: usize = 0;

    'each_read: for current_sequence in sequences[start_idx..end_idx].iter() {
        let mut add_hashmap_cnt: usize = 0;
        let mut l_window_cnt: usize = 0;
        loop_cnt += 1;
        l_window_start_idx = 0;
        if current_sequence.len() < config.l_len || current_sequence.len() < config.r_len {
            continue 'each_read;
        }
        'each_l_window: loop {
            l_window_end_idx = l_window_start_idx + config.l_len;
            if l_window_end_idx >= current_sequence.len() + 1 {
                let end: std::time::Duration = start_time.elapsed();
                eprintln!("hs loop[{:02}]({:05}-{:05},length is {})\t{:05?}({:.4}%)\tlength: {}\tsec: {}.{:03}\tadd_hashmap_cnt: {}\tl_window_cnt: {}, lr_tuple_hashmap.len():{}",
                    thread_id,
                    start_idx,
                    end_idx,
                    end_idx - start_idx,
                    loop_cnt,
                    loop_cnt as f64 / (end_idx - start_idx) as f64 * 100f64,
                    current_sequence.len(),
                    end.as_secs() - previous_time.as_secs(),
                    end.subsec_millis() - previous_time.subsec_millis(),
                    add_hashmap_cnt,
                    l_window_cnt,
                    lr_tuple_hashmap.len()
                );
                previous_time = end;
                continue 'each_read;
            }
            l_window_cnt += 1;
            let (l_has_repeat_bool, l_has_repeat_offset) =
                current_sequence.has_repeat(l_window_start_idx, l_window_end_idx);
            if l_has_repeat_bool {
                l_window_start_idx += l_has_repeat_offset + 1;
                continue 'each_l_window;
            }
            r_window_start_idx = l_window_end_idx + config.margin_size;
            'each_r_window: loop {
                r_window_end_idx = r_window_start_idx + config.r_len;
                if r_window_end_idx > current_sequence.len() {
                    l_window_start_idx += 1;
                    continue 'each_l_window;
                }
                if r_window_end_idx - l_window_start_idx > config.chunk_max - config.r_len + config.margin_size {
                    break 'each_r_window;
                }
                let (r_has_repeat_bool, r_has_repeat_offset) =
                    current_sequence.has_repeat(r_window_start_idx, r_window_end_idx);
                if r_has_repeat_bool {
                    r_window_start_idx += r_has_repeat_offset + 1;
                    continue 'each_r_window;
                }
                add_hashmap_cnt += 1;
                let lmr_string: u128 = current_sequence.subsequence_as_u128(vec![
                    [l_window_start_idx, l_window_end_idx],
                    [r_window_start_idx, r_window_end_idx],
                ]);
                if lr_tuple_hashmap.len() > high_occurence_lr_tuple.len() {
                    break 'each_read;
                }
                if high_occurence_lr_tuple.contains(&lmr_string) {
                    *lr_tuple_hashmap.entry(lmr_string).or_insert(0) += 1;
                }
                r_window_start_idx += 1;
            }
            l_window_start_idx += 1;
        }
    }
    return lr_tuple_hashmap;
}

pub fn aggregate_length_between_lr_tuple(
    sequences: &Vec<DnaSequence>,
    thread_id: usize,
    primer: &Vec<(Vec<u8>, DnaSequence, DnaSequence)>,
    product_size_max: usize,
) -> Vec<u8> {
    let mut l_window_start: usize;
    let mut l_window_end: usize;
    let mut r_window_start: usize;
    let mut r_window_end: usize;
    let mut primer_l_seq: u128;
    let mut primer_l_size: usize;
    let mut primer_r_seq: u128;
    let mut primer_r_size: usize;
    let mut mask_l: u128;
    let mut mask_r: u128;
    let mut loop_cnt: usize = 0;
    let mut ret_array: Vec<u8> = Vec::with_capacity(4_000_000_000);
    let mut lr_hit_counter: usize = 0;
    let mut l_hit_counter: usize = 0;
    eprintln!("[{}]primer pairs: {}", thread_id, primer.len());

    let start_time = Instant::now();
    let mut previous_time = start_time.elapsed();
    '_each_primer: for current_primer in primer.iter() {
        primer_l_size = current_primer.1.len();
        primer_l_seq = current_primer
            .1
            .subsequence_as_u128(vec![[0, primer_l_size]]);
        primer_r_size = current_primer.2.len();
        primer_r_seq = current_primer
            .2
            .subsequence_as_u128(vec![[0, primer_r_size]]);
        mask_l = u128::MAX >> (64 - primer_l_size) * 2;
        mask_r = u128::MAX >> (64 - primer_r_size) * 2;
        loop_cnt += 1;

        'each_read: for current_sequence in sequences.iter() {
            l_window_start = 0;
            'each_l_window: loop {
                l_window_end = l_window_start + primer_l_size;
                if l_window_end >= current_sequence.len() + 1 {
                    continue 'each_read;
                }
                let l_window_as_u128: u128 =
                    current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end]]);
                if (l_window_as_u128 & mask_l) != primer_l_seq {
                    l_window_start += 1;
                    continue 'each_l_window;
                }

                r_window_start = l_window_end;
                l_hit_counter += 1;
                'each_r_window: loop {
                    r_window_end = r_window_start + primer_r_size;
                    if r_window_end >= current_sequence.len() + 1 {
                        l_window_start += 1;
                        continue 'each_l_window;
                    }
                    if r_window_end - l_window_start > product_size_max {
                        l_window_start += 1;
                        continue 'each_l_window;
                    }
                    let r_window_as_u128: u128 =
                        current_sequence.subsequence_as_u128(vec![[r_window_start, r_window_end]]);
                    if (r_window_as_u128 & mask_r) != primer_r_seq {
                        r_window_start += 1;
                        continue 'each_r_window;
                    }

                    let primer_id = &current_primer.0;
                    let sequence_slice = current_sequence.decode(l_window_start, r_window_end);
                    let length = r_window_end - l_window_start;

                    ret_array.push(b'>');
                    ret_array.extend(primer_id);
                    ret_array.push(b'_');

                    for digit in length.to_string().as_bytes() {
                        ret_array.push(*digit);
                    }

                    ret_array.push(b'\n');
                    ret_array.extend(sequence_slice);
                    ret_array.push(b'\n');
                    lr_hit_counter += 1;
                    r_window_start += 1;
                }
            }
        }
        let end = start_time.elapsed();
        eprintln!(
            "loop[{:02?}]: {:06?}\t{:09?}\t{}\t{}\tsec: {}.{:03}",
            thread_id,
            primer.len(),
            loop_cnt,
            lr_hit_counter,
            l_hit_counter,
            end.as_secs() - previous_time.as_secs(),
            end.subsec_millis() - previous_time.subsec_millis()
        );
        previous_time = end;
        lr_hit_counter = 0;
        l_hit_counter = 0;
    }
    return ret_array;
}
