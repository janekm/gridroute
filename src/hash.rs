//! Fast 128-bit content hash of a window of a dense array (for memoising searches on unchanged board regions).
//! Not cryptographic: two independent 64-bit multiply-xorshift lanes over 8-byte words.

const K1: u64 = 0x9E37_79B9_7F4A_7C15;
const K2: u64 = 0xC2B2_AE3D_27D4_EB4F;

#[inline(always)]
fn mix(h: u64, v: u64, k: u64) -> u64 {
    let x = (h ^ v).wrapping_mul(k);
    x ^ (x >> 29)
}

pub struct Hasher {
    a: u64,
    b: u64,
}

impl Hasher {
    pub fn new(seed: u64) -> Self {
        Hasher { a: seed ^ K1, b: seed.rotate_left(32) ^ K2 }
    }
    pub fn bytes(&mut self, s: &[u8]) {
        let mut ch = s.chunks_exact(8);
        for c in &mut ch {
            let v = u64::from_le_bytes(c.try_into().unwrap());
            self.a = mix(self.a, v, K1);
            self.b = mix(self.b, v.rotate_left(17), K2);
        }
        let r = ch.remainder();
        let mut t = [0u8; 8];
        t[..r.len()].copy_from_slice(r);
        let v = u64::from_le_bytes(t) ^ ((r.len() as u64) << 56);
        self.a = mix(self.a, v, K1);
        self.b = mix(self.b, v.rotate_left(17), K2);
    }
    pub fn finish(&self) -> (u64, u64) {
        (mix(self.a, self.b, K2), mix(self.b, self.a, K1))
    }
}

/// Hash rows i0..=i1, columns j0..=j1 of every plane of a [nl][ny][nx] array of `esize`-byte elements.
pub unsafe fn hash_window(p: *const u8, nl: usize, ny: usize, nx: usize, esize: usize, win: (usize, usize, usize, usize), h: &mut Hasher) {
    let (i0, j0, i1, j1) = win;
    let rowlen = (j1 - j0 + 1) * esize;
    for l in 0..nl {
        for i in i0..=i1 {
            let row = p.add(((l * ny + i) * nx + j0) * esize);
            h.bytes(std::slice::from_raw_parts(row, rowlen));
        }
    }
}
