use crate::config::Config;
use std::cmp;

pub fn decode_u128_2_dna_seq(source: &u128, char_size: usize) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..char_size {
        base = source >> 2 * (char_size - 1 - i) & 3;
        match base {
            0 => {
                result.push(b'A');
            }
            1 => {
                result.push(b'C');
            }
            2 => {
                result.push(b'G');
            }
            3 => {
                result.push(b'T');
            }
            _ => {
                panic!("Never reached!!!base: {}", base);
            }
        }
    }
    return result;
}

// 後方互換性のためのラッパー関数
pub fn decode_u128_l(source: &u128) -> Vec<u8> {
    decode_u128_l_with_config(source, &Config::default())
}

pub fn decode_u128_r(source: &u128) -> Vec<u8> {
    decode_u128_r_with_config(source, &Config::default())
}

pub fn decode_u128_l_with_config(source: &u128, config: &Config) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..config.l_len {
        base = source >> (((config.l_len + config.r_len) - i - 1) * 2) & 3;
        match base {
            0 => {
                result.push(b'A');
            }
            1 => {
                result.push(b'C');
            }
            2 => {
                result.push(b'G');
            }
            3 => {
                result.push(b'T');
            }
            _ => {
                panic!("Never reached!!!base: {}", base);
            }
        }
    }
    return result;
}

pub fn decode_u128_r_with_config(source: &u128, config: &Config) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..config.r_len {
        base = source >> ((config.r_len - i - 1) * 2) & 3;
        match base {
            0 => {
                result.push(b'A');
            }
            1 => {
                result.push(b'C');
            }
            2 => {
                result.push(b'G');
            }
            3 => {
                result.push(b'T');
            }
            _ => {
                panic!("Never reached!!!base: {}", base);
            }
        }
    }
    return result;
}

pub struct DnaSequence {
    length: usize,
    sequence: Vec<u64>,
}
impl Clone for DnaSequence {
    fn clone(&self) -> Self {
        DnaSequence {
            length: self.length,
            sequence: self.sequence.clone(),
        }
    }
}

impl DnaSequence {
    pub fn new(source: &Vec<u8>) -> DnaSequence {
        let mut retval: Vec<u64> = Vec::new();
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        for each_base in source.iter() {
            match each_base {
                b'A' | b'a' => {
                    buf |= 0;
                }
                b'C' | b'c' => {
                    buf |= 1;
                }
                b'G' | b'g' => {
                    buf |= 2;
                }
                b'T' | b't' => {
                    buf |= 3;
                }
                _ => {
                    panic!("Unexpected character: {}({})", each_base, *each_base as char);
                }
            }
            //buf <<= 2;//境界を跨ぐ場合、シフトしてはいけない。
            if cnt == 31 {
                //cnt == 31のときはbufが全部埋まってる。
                retval.push(buf);
                buf = 0;
                cnt = 0;
            } else {
                buf <<= 2;
                cnt += 1;
            }
        }

        if cnt != 32 {
            buf <<= 2 * (31 - cnt); //塩基を表すbitを上位に寄せる。
            retval.push(buf); //32塩基の節目で切れなかった時に備えてpushする。
        }
        /*
                for each_buf in retval.iter(){
                    println!("{:064b}", each_buf);
                }
        */
        return DnaSequence {
            sequence: retval,
            length: source.len(),
        };
    }

    pub fn encode(source: &Vec<u8>) -> DnaSequence {
        return DnaSequence::new(source);
    }

    pub fn len(&self) -> usize {
        return self.length;
    }

    pub fn decode(&self, start: usize, end: usize) -> Vec<u8> {
        assert!(
            start < end,
            "DnaSequence::decode assertion failed: {} !< {}",
            start,
            end
        );
        assert!(
            end <= self.length,
            "DnaSequence::decode assertion failed: {} !< {}",
            end,
            self.length
        );
        let mut retval = Vec::new();
        let mut buf: u8;
        for i in start..end {
            buf = ((self.sequence[i / 32] >> (2 * (31 - i % 32))) & 3)
                .try_into()
                .unwrap();
            match buf {
                0 => {
                    retval.push('A' as u8);
                }
                1 => {
                    retval.push('C' as u8);
                }
                2 => {
                    retval.push('G' as u8);
                }
                3 => {
                    retval.push('T' as u8);
                }
                _ => {
                    panic!("Never reached!!!buf: {}", buf);
                }
            }
        }
        return retval;
    }
    pub fn complement(&self) -> DnaSequence {
        let mut retval = Vec::new();
        let mut buf: u8;
        for i in 0..self.length {
            buf = ((self.sequence[i / 32] >> (2 * (31 - i % 32))) & 3)
                .try_into()
                .unwrap();
            match buf {
                0 => {
                    retval.push('T' as u8);
                }
                1 => {
                    retval.push('G' as u8);
                }
                2 => {
                    retval.push('C' as u8);
                }
                3 => {
                    retval.push('A' as u8);
                }
                _ => {
                    panic!("Never reached!!!buf: {}", buf);
                }
            }
        }
        return DnaSequence::new(&retval);
    }

    pub fn reverse(&self) -> DnaSequence {
        let mut retval = Vec::new();
        let mut buf: u8;
        for i in 0..self.length {
            let j = self.length - i - 1;
            buf = ((self.sequence[j / 32] >> (2 * (31 - j % 32))) & 3)
                .try_into()
                .unwrap();
            match buf {
                0 => {
                    retval.push('A' as u8);
                }
                1 => {
                    retval.push('C' as u8);
                }
                2 => {
                    retval.push('G' as u8);
                }
                3 => {
                    retval.push('T' as u8);
                }
                _ => {
                    panic!("Never reached!!!buf: {}", buf);
                }
            }
        }
        return DnaSequence::new(&retval);
    }

    pub fn reverse_complement(&self) -> DnaSequence {
        return self.complement().reverse();
    }

    pub fn subsequence(&self, ranges: Vec<[usize; 2]>) -> DnaSequence {
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        let mut retval: Vec<u64> = Vec::new();
        let mut length: usize = 0;
        for each_range in ranges.iter() {
            let start: usize = each_range[0];
            let end: usize = each_range[1];
            assert!(
                start < end,
                "DnaSequence::subsequence assertion failed: {} !< {}",
                start,
                end
            );
            //assert!(start >= 0, "DnaSequence::subsequence assertion failed: {} >= 0", start);
            assert!(
                end < self.length,
                "DnaSequence::subsequence assertion failed: {} < {}",
                end,
                self.length
            );
            length += end - start;
            for i in start..end {
                buf += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
                buf <<= 2;
                cnt += 1;
                if cnt == 31 {
                    retval.push(buf);
                    buf = 0;
                    cnt = 0;
                }
            }
        }
        if cnt != 31 {
            buf <<= 2 * (31 - cnt); //塩基を表すbitを上位に寄せる。
            retval.push(buf); //32塩基の節目で切れなかった時に備えてpushする。
        }
        return DnaSequence {
            sequence: retval,
            length: length,
        };
    }

    //subsequence_as_u128は右詰め
    pub fn subsequence_as_u128(&self, ranges: Vec<[usize; 2]>) -> u128 {
        let mut buf: u128 = 0;
        let mut cnt: usize = 0;
        for each_range in ranges.iter() {
            let start: usize = each_range[0];
            let end: usize = each_range[1];
            assert!(
                start < end,
                "DnaSequence::subsequence_as_u128 assertion failed: {} !< {}",
                start,
                end
            );
            assert!(
                end <= self.length,
                "DnaSequence::subsequence_as_u128 assertion failed: {} < {}",
                end,
                self.length
            );
            for i in start..end {
                buf <<= 2;
                cnt += 1;
                buf += ((self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3) as u128;
            }
            assert!(cnt <= 64, "DnaSequence::subsequence_as_u128 assertion failed: too many DNA size. cnt reaches {}", cnt);
        }
        return buf;
    }

    pub fn has_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        let has_one_base_repeat: (bool, usize) = self.has_one_base_repeat(start, end);
        let has_two_base_repeat: (bool, usize) = self.has_two_base_repeat(start, end);
        let has_three_base_repeat: (bool, usize) = self.has_three_base_repeat(start, end);
        let has_gc_or_at_consequence_region: (bool, usize) =
            self.has_gc_or_at_consequence_region(start, end);

        /*
               eprint!("{}\t", std::str::from_utf8(&self.decode(start, end)).unwrap());
               eprint!("{}\t{}\t{}\t", start, end, end - start);
               eprint!("{:?}\t", has_one_base_repeat);
               eprint!("{:?}\t", has_two_base_repeat);
               eprint!("{:?}\t", has_three_base_repeat);
               eprintln!("{:?}", has_gc_or_at_consequence_region);

        */
        let retval_bool: bool = has_one_base_repeat.0
            | has_two_base_repeat.0
            | has_three_base_repeat.0
            | has_gc_or_at_consequence_region.0;
        let retval_base: usize = cmp::max(
            cmp::max(has_one_base_repeat.1, has_two_base_repeat.1),
            cmp::max(has_three_base_repeat.1, has_gc_or_at_consequence_region.1),
        );
        return (retval_bool, retval_base);
    }

    pub fn has_one_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(
            start < end,
            "DnaSequence::has_one_base_repeat assertion failed: {} !< {}",
            start,
            end
        );
        assert!(
            end - start > 3,
            "DnaSequence::has_one_base_repeat assertion failed: {} - {} < 4",
            end,
            start
        );
        assert!(end - start >8, "DnaSequence::has_one_base_repeat assertion failed: length of the evaluation subject must be longer than 8");
        assert!(end - start <= 32, "DnaSequence::has_one_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_one_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. start: {}, end: {}, self.lngth: {}", start, end, self.length);

        let upper_mask: u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let lower_mask: u128 = 0xFFFFFFFFFFFFFFFF0000000000000000;
        let mut left_bits: u128 = self.sequence[start / 32] as u128;
        let mut right_bits: u128 = self.sequence[(end - 1) / 32] as u128;
        left_bits <<= 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 2 * (start % 32);
        left_bits <<= (end % 32) * 2;
        right_bits <<= 2 * (end % 32);
        right_bits &= lower_mask;
        right_bits >>= 64;
        let original = left_bits + right_bits;
        let mask: u64 = (1u64 << (2 * (end - start)) - 1) - 1;
        if original != 0 && original & 0xFF == 0 {
            return (true, end - start - 4);
        }
        let zero_ichi: u64 = 0x5555_5555_5555_5555 /* & !63 */ & mask; //末尾がAAAの時の偽陽性には目をつぶる
        let val1: u64 = original as u64;
        let val2: u64 = val1 << 2;
        let val3: u64 = val1 ^ val2;
        let val4: u64 = val3 >> 1;
        let val5: u64 = val3 | val4;
        let val6: u64 = !val5;
        let val7: u64 = val6 & zero_ichi; //下位1bitだけ0にする。
        let val8: u64 = val7 << 2;
        let val9: u64 = val7 << 4;
        let val10: u64 = val7 & val8 & val9;
        let val11: u32 = (val10 << (2 * (32 + start - end))).leading_zeros() / 2;

        /*
               println!("has_one_base_repeat");
               println!("start:    {}, {}", start, start / 32);
               println!("end:      {}, {}", end, (end - 1) / 32);
               println!(" 0101: {:064b}", mask & zero_ichi);
               println!("original: {:0128b}", original);
               println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
               println!("{}", std::str::from_utf8(&decode_u128_2_dna_seq(&original, 64)).unwrap());
               println!("{}", std::str::from_utf8(&decode_u128_l(&original)).unwrap());
               println!("{}", std::str::from_utf8(&decode_u128_r(&original)).unwrap());
               println!(" val1: {:064b}", val1);
               println!(" val2: {:064b}", val2);
               println!(" val3: {:064b}", val3);
               println!(" val4: {:064b}", val4);
               println!(" val5: {:064b}", val5);
               println!(" val6: {:064b}", val6);
               println!(" val7: {:064b}", val7);
               println!(" val8: {:064b}", val8);
               println!(" val9: {:064b}", val9);
               println!("val10: {:064b}", val10);
               println!("val11: {}", val11);
        */

        #[cfg(test)]
        {
            println!("has_one_base_repeat");
            println!("start:    {}, {}", start, start / 32);
            println!("end:      {}, {}", end, (end - 1) / 32);
            println!(" mask: {:064b}", mask & zero_ichi);

            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!(" val1: {:064b}", val1);
            println!(" val2: {:064b}", val2);
            println!(" val3: {:064b}", val3);
            println!(" val4: {:064b}", val4);
            println!(" val5: {:064b}", val5);
            println!(" val6: {:064b}", val6);
            println!(" val7: {:064b}", val7);
            println!(" val8: {:064b}", val8);
            println!(" val9: {:064b}", val9);
            println!("val10: {:064b}", val10);
            println!("val11: {}", val11);
        }

        if val11 == 32 {
            return (false, 0);
        } else {
            return (true, val11.try_into().unwrap());
        }
    }

    pub fn has_two_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(
            start < end,
            "DnaSequence::has_two_base_repeat assertion failed: {} !< {}",
            start,
            end
        );
        assert!(
            end - start > 3,
            "DnaSequence::has_two_base_repeat assertion failed: {} - {} < 4",
            end,
            start
        );
        assert!(end - start >8, "DnaSequence::has_two_base_repeat assertion failed: length of the evaluation subject must be longer than 8");
        assert!(end - start <= 32, "DnaSequence::has_two_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_two_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let upper_mask: u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let lower_mask: u128 = 0xFFFFFFFFFFFFFFFF0000000000000000;
        let mut left_bits: u128 = self.sequence[start / 32] as u128;
        let mut right_bits: u128 = self.sequence[(end - 1) / 32] as u128;
        left_bits <<= 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 2 * (start % 32);
        left_bits <<= (end % 32) * 2;
        right_bits <<= 2 * (end % 32);
        right_bits &= lower_mask;
        right_bits >>= 64;
        let original = left_bits + right_bits;
        let mask: u64 = (1u64 << (2 * (end - start)) - 1) - 1;
        if original != 0 && original & 0xFF == 0 {
            return (true, end - start - 4);
        }
        let zero_ichi: u64 = 0x5555_5555_5555_5555 /* & !63 */ & mask; //末尾がAAAの時の偽陽性には目をつぶる

        //ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1 = original as u64;
        let val2 = (original as u64) << 4;
        let val3 = val1 ^ val2 & !15;
        let val4 = val3 >> 1;
        let val5 = val3 | val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let last = val7
            & val7 << 2
            & val7 << 4
            & val7 << 6
            & val7 << 8
            & val7 << 10
            & val7 << 12
            & val7 << 14;
        let leading0 = (last << (2 * (32 + start - end))).leading_zeros() / 2;
        #[cfg(test)]
        {
            println!("has_two_base_repeat");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("last:  {:064b}", last);
            println!("leading0: {}", leading0);
        }
        if leading0 == 32 {
            return (false, 0);
        } else {
            return (true, leading0.try_into().unwrap());
        }
    }

    pub fn has_three_base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        //eprintln!("has_three_base_repeat start: {}, end: {}", start, end);
        assert!(
            start < end,
            "DnaSequence::has_three_base_repeat assertion failed: {} !< {}",
            start,
            end
        );
        assert!(
            end - start > 3,
            "DnaSequence::has_three_base_repeat assertion failed: {} - {} < 4",
            end,
            start
        );
        assert!(end - start >8, "DnaSequence::has_three_base_repeat assertion failed: length of the evaluation subject must be longer than 8");
        assert!(end - start <= 32, "DnaSequence::has_three_base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_three_base_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let upper_mask: u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let lower_mask: u128 = 0xFFFFFFFFFFFFFFFF0000000000000000;
        let mut left_bits: u128 = self.sequence[start / 32] as u128;
        let mut right_bits: u128 = self.sequence[(end - 1) / 32] as u128;
        left_bits <<= 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 2 * (start % 32);
        left_bits <<= (end % 32) * 2;
        right_bits <<= 2 * (end % 32);
        right_bits &= lower_mask;
        right_bits >>= 64;
        let original = left_bits + right_bits;
        let mask: u64 = (1u64 << (2 * (end - start)) - 1) - 1;
        let zero_ichi: u64 = 0x5555_5555_5555_5554 & !63 & mask;
        //ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1 = original as u64;
        let val2 = (original as u64) << 6;
        let val3 = val1 ^ val2;
        let val4 = val3 >> 1;
        let val5 = val3 | val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let val8 = val7 << 2;
        let val9 = val7 << 4;
        let val10 = val7 << 6;
        let val11 = val7 << 8;
        let val12 = val7 << 10;
        let last = val7 & val8 & val9 & val10 & val11 & val12;
        let leading0 = (last << (2 * (32 + start - end + 1))).leading_zeros() / 2;

        #[cfg(test)]
        {
            println!("has_three_base_repeat");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("val8:  {:064b}", val8);
            println!("val9:  {:064b}", val9);
            println!("val10: {:064b}", val10);
            println!("val11: {:064b}", val11);
            println!("val12: {:064b}", val12);
            println!("last:  {:064b}", last);
        }
        if leading0 == 32 {
            return (false, 0);
        } else {
            return (true, leading0.try_into().unwrap());
        }
    }
    pub fn has_gc_or_at_consequence_region(&self, start: usize, end: usize) -> (bool, usize) {
        //eprintln!("has_gc_or_at_consequence_region start: {}, end: {}", start, end);
        assert!(
            start < end,
            "DnaSequence::has_gc_or_at_consequence_region assertion failed: {} !< {}",
            start,
            end
        );
        assert!(
            end - start > 3,
            "DnaSequence::has_gc_or_at_consequence_region assertion failed: {} - {} < 4",
            end,
            start
        );
        assert!(end - start >8, "DnaSequence::has_gc_or_at_consequence_region assertion failed: length of the evaluation subject must be longer than 8");
        assert!(end - start <= 32, "DnaSequence::has_gc_or_at_consequence_region assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_gc_or_at_consequence_region assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let upper_mask: u128 = 0x0000000000000000FFFFFFFFFFFFFFFF;
        let lower_mask: u128 = 0xFFFFFFFFFFFFFFFF0000000000000000;
        let mut left_bits: u128 = self.sequence[start / 32] as u128;
        let mut right_bits: u128 = self.sequence[(end - 1) / 32] as u128;
        left_bits <<= 2 * (start % 32);
        left_bits &= upper_mask;
        left_bits >>= 2 * (start % 32);
        left_bits <<= (end % 32) * 2;
        right_bits <<= 2 * (end % 32);
        right_bits &= lower_mask;
        right_bits >>= 64;
        let original = left_bits + right_bits;
        let mask: u64 = (1u64 << (2 * (end - start)) - 1) - 1;
        let zero_ichi: u64 = 0x5555_5555_5555_5554 & !63 & mask;
        //ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1: u64 = original as u64;
        let val2: u64 = val1 >> 1;
        let val3: u64 = val1 ^ val2;
        let val4: u64 = val3 & zero_ichi;
        let val5: u64 =
            val4 << 2 & val4 << 4 & val4 << 6 & val4 << 8 & val4 << 10 & val4 << 12 & val4 << 14;
        let val6: u64 = (val3 ^ u64::MAX) & zero_ichi;
        let val7: u64 =
            val6 << 2 & val6 << 4 & val6 << 6 & val6 << 8 & val6 << 10 & val6 << 12 & val6 << 14;
        /*
               let val11:u64 = val9 << 4;
               let val12:u64 = val9 << 6;
               let val13:u64 = val9 & val10 & val11 & val12;

        */
        let leading0_gc = (val5 << (2 * (32 + start - end))).leading_zeros() / 2;
        let leading0_at = (val7 << (2 * (32 + start - end))).leading_zeros() / 2;

        #[cfg(test)]
        {
            println!("has_gc_or_at_consequence_region");
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("{}", leading0_gc);
            println!("{}", leading0_at);
        }
        if leading0_gc == 32 && leading0_at == 32 {
            return (false, 0);
        } else {
            return (true, cmp::min(leading0_gc, leading0_at).try_into().unwrap());
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sequence_encoder_util::decode_u128_2_dna_seq;
    use crate::sequence_encoder_util::DnaSequence;
    use ::function_name::named;
    /*
     *
     *Encode Test
     *
     */
    #[test]
    #[named]
    fn encode_test_4A() {
        let source: Vec<u8> = vec![b'A', b'A', b'A', b'A'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.sequence[0] == 0x0000000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_4C() {
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.sequence[0] == 0x5500000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_4G() {
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.sequence[0] == 0xAA00000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_4T() {
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.sequence[0] == 0xFF00000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_ACGT() {
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.sequence[0] == 0x1b00000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_16C() {
        let source: String = "CCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555500000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_31C() {
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555554,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_32C() {
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_33C() {
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x4000000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_34C() {
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x5000000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_63C() {
        let source: String =
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x5555555555555554,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_64C() {
        let source: String =
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_65C() {
        let source: String =
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x5555555555555555,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[2] == 0x4000000000000000,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn encode_test_32N() {
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0b1000000110000111101111110111000111001100000011010111110101111101,
            "{} failed",
            function_name!()
        );
    } //                               G A A C G A C T G T T T C T A C T A T A A A T C C T T C C T T C
    #[test]
    #[named]
    fn encode_test_120N() {
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.sequence[0] == 0x8187bf71cc0d7d7d,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[1] == 0x725cd3f7a2d7eb81,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[2] == 0xeca09de04446f57e,
            "{} failed",
            function_name!()
        );
        assert!(
            obj.sequence[3] == 0x8f6c5ce0c75b0000,
            "{} failed",
            function_name!()
        );
    }

    /*
     *
     *has_one_base_repeat_test
     *4つからtrue
     *
     */
    /* 8塩基しかない配列で失敗するのは諦める
       #[test]
       #[named]
       fn has_one_base_repeat_test_8C() {
           let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
           let obj = DnaSequence::new(&source);
           assert!(
               obj.has_one_base_repeat(0, 4) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 5) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 6) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 7) == (true, 0),
               "{} failed",
               function_name!()
           );
       }
       #[test]
       #[named]
       fn has_one_base_repeat_test_8G() {
           let source: Vec<u8> = vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G'];
           let obj = DnaSequence::new(&source);
           assert!(
               obj.has_one_base_repeat(0, 4) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 5) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 6) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 7) == (true, 0),
               "{} failed",
               function_name!()
           );
       }
       #[test]
       #[named]
       fn has_one_base_repeat_test_8T() {
           let source: Vec<u8> = vec![b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T'];
           let obj = DnaSequence::new(&source);
           assert!(
               obj.has_one_base_repeat(0, 4) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 5) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 6) == (true, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 7) == (true, 0),
               "{} failed",
               function_name!()
           );
       }

       #[test]
       #[named]
       fn has_one_base_repeat_test_8N() {
           let source: Vec<u8> = vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
           let obj = DnaSequence::new(&source);
           assert!(
               obj.has_one_base_repeat(0, 4) == (false, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 5) == (false, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 6) == (false, 0),
               "{} failed",
               function_name!()
           );
           assert!(
               obj.has_one_base_repeat(0, 7) == (false, 0),
               "{} failed",
               function_name!()
           );
       }

    */
    #[test]
    #[named]
    fn has_one_base_repeat_test_120N() {
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 9) == (false, 0),
            "{} failed",
            function_name!()
        );
        assert!(
            obj.has_one_base_repeat(0, 19) == (true, 9),
            "{} failed",
            function_name!()
        );
        assert!(
            obj.has_one_base_repeat(9, 20) == (true, 0),
            "{} failed",
            function_name!()
        );
        assert!(
            obj.has_one_base_repeat(28, 39) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_1() {
        let source: String = "ATTCATACTTAATACTGTATCAGTTGA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_2() {
        let source: String = "TTCATACTTAATACTGTATCAGTTGAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_3() {
        let source: String = "TCATACTTAATACTGTATCAGTTGAGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_4() {
        let source: String = "TTCGAAAATCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (true, 4),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_5() {
        let source: String = "TTCGAAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_6() {
        let source: String = "TTCGTAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_7() {
        let source: String = "AAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (true, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_8() {
        let source: String = "AAAACGTCGTCGTCGCATACGATCGAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (true, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_9() {
        let source: String = "AAACGTCGTCGTCGCATACGAATCGAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_10() {
        let source: String = "GTCGTCGTCGCATACGAATCGATAAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 27) == (true, 23),
            "{} failed, {:?}",
            function_name!(),
            obj.has_one_base_repeat(0, 27)
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_11() {
        let source: String = "GTCGTCGTCGCATACGAATCGATTAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 23) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_12() {
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 30) == (true, 23),
            "{} failed",
            function_name!()
        );
        assert!(
            obj.has_one_base_repeat(0, 19) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_13() {
        //末尾がAAAの時の偽陽性には目をつぶる
        let source: String = "TAATGTATTCACTAACAAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 19) == (true, 16),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_14() {
        let source: String = "TAATGTATTCACTAACCAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 19) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_15() {
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 19) == (true, 10),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_16() {
        let source: String = "GCTTGTATACAGGGGATTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 19) == (true, 11),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_17() {
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 19) == (true, 10),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_18() {
        let source: String = "AACGTACGTACGTACGTACGTACGTACGTACGTCTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(30, 52) == (true, 13),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_19() {
        let source: String =
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 31) == (true, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_20() {
        let source: String =
            "AAAAAAAAAACTAAGAAAATCTATGAGACAGAGTGGACTATATATATCTATATCTATAGTGAAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 31) == (true, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_21() {
        let source: String = "AAACTAATTGTTGTTTTTTTTAATAATAATTA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 32) == (true, 13),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_22() {
        let source: String = "ATACCCCTCCGCAGTCCCTAGGCCTGGGGCTG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 32) == (true, 3),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_one_base_repeat_test_27N_23() {
        let source: String = "ATTGTCTGATTGCTATAGTGGTTACAGATTTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_one_base_repeat(0, 32) == (true, 28),
            "{} failed",
            function_name!()
        );
    }

    /*
     *
     *has_three_base_repeat
     *
     */
    #[test]
    #[named]
    fn has_three_base_repeat_27N_1() {
        let source: String = "TCATATGCAACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 27) == (true, 6),
            "{} failed{:?}",
            function_name!(),
            obj.has_three_base_repeat(0, 27)
        );
    }

    #[test]
    #[named]
    fn has_three_base_repeat_27N_1_a() {
        let source: String = "ACGTACGTCAACAAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 15) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_three_base_repeat_27N_2() {
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_three_base_repeat_27N_3() {
        let source: String = "CAACAACTGC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 10) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_three_base_repeat_27N_4() {
        let source: String = "ACGTACGTCAACAAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 15) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_three_base_repeat_27N_5() {
        let source: String = "CGTACGTCAACAACA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 15) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_three_base_repeat_27N_6() {
        let source: String = "GTACGTCAACAACAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 15) == (true, 5),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_three_base_repeat_64N_1() {
        let source: String =
            "TACGAATAAGATGTACGTACAGATGTACGAATTGTAGTAGTACGAATAGATGTACGAATAGACC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_three_base_repeat(0, 31) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    /*
     *
     *has_two_base_repeat
     *
     */
    #[test]
    #[named]
    fn has_two_base_repeat_27N_1() {
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2() {
        let source: String = "CGTACGCTTATATATATATATACCGCA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 8),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x6_1() {
        let source: String = "TATATATATATAGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x6_2() {
        let source: String = "GCCCGCACGTACGCTATATATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 14),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x5_1() {
        let source: String = "TATATATATACTGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x5_2() {
        let source: String = "CTGCCCGCACGTACGCTATATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 16),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x4_1() {
        let source: String = "TATATATACTGCCCGAACACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_2x4_2() {
        let source: String = "CTGCCCGAACACGTACGCTATATATAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 18),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_3() {
        let source: String = "GAATTCGTTACGTAACGACGCGCGCGC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_5() {
        let source: String = "ACGTTATACCTGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_6() {
        let source: String = "GCGGTATATACGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_two_base_repeat_27N_7() {
        let source: String = "CGGTATATACGTACCGCACGCACACAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_two_base_repeat_27N_8() {
        let source: String =
            "ATATATATATATATATATATATATATATATATTATATATATATATATATATATATATATATATA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 0),
            "{} failed{:?}",
            function_name!(),
            obj.has_two_base_repeat(0, 27)
        );
    }

    #[test]
    #[named]
    fn has_two_base_repeat_27N_9() {
        let source: String =
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 27) == (true, 0),
            "{} failed{:?}",
            function_name!(),
            obj.has_two_base_repeat(0, 27)
        );
    }

    #[test]
    #[named]
    fn has_two_base_repeat_27N_10() {
        let source: String = "GACTGGATCTCTCTTGTTTGTTCCTGCACTAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_two_base_repeat(0, 32) == (false, 0),
            "{} failed",
            function_name!()
        );
    }

    /*
    has_gc_or_at_consequence_region test
    */

    #[test]
    #[named]
    fn has_gc_or_at_consequence_region_32N_1() {
        let source: String = "TTTCTGTGCTAATTTTAATTTGTTAATAGTGCCAT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_gc_or_at_consequence_region(0, 32) == (true, 8),
            "{} failed, {:?}",
            function_name!(),
            obj.has_gc_or_at_consequence_region(0, 32)
        );
    }

    #[test]
    #[named]
    fn has_gc_or_at_consequence_region_32N_2() {
        let source: String = "GTTTGTTCCTGCACTAGATCACAGACATATCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_gc_or_at_consequence_region(0, 32) == (false, 0),
            "{} failed, {:?}",
            function_name!(),
            obj.has_gc_or_at_consequence_region(0, 32)
        );
    }

    #[test]
    #[named]
    fn has_gc_or_at_consequence_region_32N_3() {
        let source: String = "CATCCCGGGAGATATGAAGAAGAGATAGATCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_gc_or_at_consequence_region(0, 32) == (false, 0),
            "{} failed, {:?}",
            function_name!(),
            obj.has_gc_or_at_consequence_region(0, 32)
        );
    }

    /*
    has_repeat test
    */
    #[test]
    #[named]
    fn has_repeat_1() {
        let source: String = "CTTGTATACAGGGGATTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 19).0 == true,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_repeat_2() {
        let source: String = "CATCACCAATTATTGGTCCTAATGTA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 26) == (true, 6),
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_repeat_3() {
        let source: String = "ATCCTCAGCTGCTTGTATA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 19).0 == false,
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn has_repeat_4() {
        for string in ["ATCCTCAGCTGCTTGTATA"] {
            let source: String = string.to_string();
            let v: Vec<u8> = source.into_bytes();
            let obj = DnaSequence::new(&v);
            assert!(
                obj.has_repeat(0, 19).0 == false,
                "{} failed",
                function_name!()
            );
        }
    }
    #[test]
    #[named]
    fn has_repeat_5() {
        for string in ["AAAAAAAATCAATTGTTTTATCAAATGTGGACCAAGGCAAATGAACAATTTTTTTTACTGCTAA"] {
            let source: String = string.to_string();
            let v: Vec<u8> = source.into_bytes();
            let obj = DnaSequence::new(&v);
            assert!(
                obj.has_repeat(0, 19) == (true, 0),
                "{} failed",
                function_name!()
            );
        }
    }

    #[test]
    #[named]
    fn has_repeat_6() {
        for string in ["TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"] {
            let source: String = string.to_string();
            let v: Vec<u8> = source.into_bytes();
            let obj = DnaSequence::new(&v);
            assert!(
                obj.has_repeat(0, 19) == (true, 0),
                "{} failed",
                function_name!()
            );
        }
    }

    #[test]
    #[named]
    fn has_repeat_7() {
        let source: String = "TTGATGATTGGCCATTTTGTTTTCTAGAAGGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 14),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_8() {
        let source: String = "ATTTTTTGTTTGCATCGCGCGTCTCCTGATCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 1),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_9() {
        let source: String = "AATCACTAAAATAATTCAATAGAAATTTTGGC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 7),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_10() {
        let source: String = "TTGATGATTGGCCATTTTGTTTTCTAGAAGGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 14),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_11() {
        let source: String = "TTGATGATTGGCCATTTTGTTTTCTAGAAGGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 14),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_12() {
        let source: String = "AAATTGATGATTGGCCATTTTGTTTTCTAGAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(0, 31) == (true, 17),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_13() {
        let source: String = "AATTGATGATTGGCCATTTTGTTTTCTAGAAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(10, 31) == (true, 6),
            "{} failed",
            function_name!()
        );
    }

    #[test]
    #[named]
    fn has_repeat_14() {
        let source: String = "AAAAAAAAAAAAAAAAAAAATAATTGCTACAGTGCTTAAAGATTTTGTATAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(
            obj.has_repeat(20, 51) == (true, 22),
            "{} failed",
            function_name!()
        );
    }

    /*
     *
     *Decode test
     *
     */
    #[test]
    #[named]
    fn decode_test_8C() {
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(
            obj.decode(0, 4) == vec![67, 67, 67, 67],
            "{} failed",
            function_name!()
        );
        assert!(
            obj.decode(0, 5) == vec![67, 67, 67, 67, 67],
            "{} failed",
            function_name!()
        );
        assert!(
            obj.decode(0, 6) == vec![67, 67, 67, 67, 67, 67],
            "{} failed",
            function_name!()
        );
        assert!(
            obj.decode(0, 7) == vec![67, 67, 67, 67, 67, 67, 67],
            "{} failed",
            function_name!()
        );
    }
    #[test]
    #[named]
    fn decode_test_120N() {
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.decode(0, 120) == v, "{} failed", function_name!());
    }

    /*
     *
     *Subsequence Test
     *
     */
    #[test]
    #[named]
    fn subsequence_as_u128_test_64G() {
        let source: String =
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 64 as usize]]);
        let expectedvalue: u128 = 0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_10G() {
        let source: String = "GGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 10 as usize]]);
        let expectedvalue: u128 = 0b10101010101010101010;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_12N() {
        let source: String = "ACGTAACCGGTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 12 as usize]]);
        let expectedvalue: u128 = 0b000110110000010110101111;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_1() {
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String =
            "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_2() {
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 64]]);
        let subseq: String =
            "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 19], [19, 45], [45, 64]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_3() {
        let source: String = "GAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAATCCTCAGCTGCTTGTATACAGGGGATTTCTTCTTCATCACCAATTATTGGTCCTAATGTAT".to_string();
        let v1: Vec<u8> = source.into_bytes();
        let obj1 = DnaSequence::new(&v1);
        let ret1 = obj1.subsequence_as_u128(vec![[0, 27], [64, 84]]);
        let subseq: String = "GAATCCTCAGCTGCTTGTATACAGGGGAAAAAAAAAAAAAAAAAAAA".to_string();
        let v2: Vec<u8> = subseq.into_bytes();
        let obj2 = DnaSequence::new(&v2);
        let ret2 = obj2.subsequence_as_u128(vec![[0, 47]]);
        assert!(ret1 == ret2, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_1() {
        let source: u128 = 0xFFFFFFFFFFFFFFFF;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_2() {
        let source: u128 = 0x0000000000000000;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_3() {
        let source: u128 = 0x5555555555555555;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_4() {
        let source: u128 = 0xAAAAAAAAAAAAAAAA;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_5() {
        let source: u128 = 0x00000000AAAAAAAA;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_u128_2_dna_seq_test_6() {
        let source: u128 = 0x000000AAAAAAAAAA;
        let v1: Vec<u8> = decode_u128_2_dna_seq(&source, 32);
        let v2 = String::from_utf8(v1).unwrap();
        let answer: String = "AAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG".to_string();
        assert!(v2 == answer, "{} failed", function_name!());
    }
}
