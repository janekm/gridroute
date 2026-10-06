//! Backends agree bit for bit with the naive reference on random grids.
use gridroute::dilate::{dilate_cpu, dilate_naive, Element, Grid3, Window};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn disc(r: f64) -> Vec<i32> {
    let rc = r.ceil() as i32;
    (-rc..=rc)
        .map(|di| {
            let d = (di * di) as f64;
            if d >= r * r - 1e-9 { -1 } else { ((r * r - d).sqrt() - 1e-9).floor() as i32 }
        })
        .collect()
}

#[test]
fn dilate_backends_match() {
    let mut rng = Rng(0x9e3779b97f4a7c15);
    for case in 0..60 {
        let (nl, ny, nx) = (1 + rng.below(4) as usize, 20 + rng.below(200) as usize, 20 + rng.below(300) as usize);
        let mut src = vec![0i16; nl * ny * nx];
        let density = 1 + rng.below(40);
        for v in src.iter_mut() {
            if rng.below(1000) < density {
                *v = 1 + rng.below(5) as i16;
            }
        }
        let extra: Vec<u8> = (0..ny * nx).map(|_| (rng.below(1000) < 3) as u8).collect();
        let g = Grid3 { ptr: src.as_ptr(), nl, ny, nx };
        let eg = Grid3 { ptr: extra.as_ptr(), nl: 1, ny, nx };
        let i0 = rng.below(ny as u64 / 2) as usize;
        let j0 = rng.below(nx as u64 / 2) as usize;
        let win = Window { i0, j0, h: 1 + rng.below((ny - i0) as u64) as usize, w: 1 + rng.below((nx - j0) as u64) as usize };
        let span = disc(0.5 + rng.below(120) as f64 / 10.0);
        let el = Element { span: &span };
        let layers: Vec<usize> = (0..nl).filter(|_| rng.below(3) > 0).collect();
        let layers = if layers.is_empty() { vec![0] } else { layers };
        let excl = [1 + rng.below(5) as i16];
        let on: Vec<bool> = (0..nl).map(|l| l == 0).collect();
        let use_extra = case % 2 == 0;
        for reduce in [false, true] {
            let planes = if reduce { 1 } else { layers.len() };
            let n = planes * win.h * win.w;
            let (mut a, mut b) = (vec![0u8; n], vec![0u8; n]);
            let ex = if use_extra { Some(&eg) } else { None };
            dilate_naive(&g, &layers, &excl, ex, &on, &win, &el, reduce, &mut a);
            dilate_cpu(&g, &layers, &excl, ex, &on, &win, &el, reduce, &mut b);
            assert_eq!(a, b, "cpu case {case} reduce {reduce}");
            #[cfg(all(target_os = "macos", feature = "metal"))]
            if gridroute::gpu::available() {
                let mut c = vec![7u8; n];
                assert!(gridroute::gpu::dilate_metal(&g, &layers, &excl, ex, &on, &win, &el, reduce, &mut c));
                assert_eq!(a, c, "metal case {case} reduce {reduce}");
            }
        }
    }
}

fn random_problem(rng: &mut Rng, nl: usize, h: usize, w: usize) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<i32>) {
    let n = nl * h * w;
    let blk: Vec<u8> = (0..n).map(|_| (rng.below(100) < 25) as u8).collect();
    let vok: Vec<u8> = (0..h * w).map(|_| (rng.below(100) < 30) as u8).collect();
    let mut tgt = vec![0u8; n];
    for _ in 0..3 {
        let k = rng.below(n as u64) as usize;
        tgt[k] = 1 - blk[k];
    }
    let src: Vec<i32> = (0..3).map(|_| rng.below(n as u64) as i32).collect();
    (blk, vok, tgt, src)
}

fn problem<'a>(nl: usize, h: usize, w: usize, b: &'a [u8], v: &'a [u8], t: &'a [u8], s: &'a [i32], mc: &'a [f32],
               lok: &'a [u8]) -> gridroute::astar::Problem<'a> {
    let (mut tx0, mut ty0, mut tx1, mut ty1) = (i32::MAX, i32::MAX, -1, -1);
    for (k, &x) in t.iter().enumerate() {
        if x != 0 {
            let c = k % (h * w);
            let (i, j) = ((c / w) as i32, (c % w) as i32);
            tx0 = tx0.min(j); tx1 = tx1.max(j); ty0 = ty0.min(i); ty1 = ty1.max(i);
        }
    }
    gridroute::astar::Problem { nl, h, w, blk: b, vok: v, tgt: t, src: s, mcost: mc, lay_ok: lok, vcost: 7.0, turn: 2.0,
        hmul: 1.2, tbox: (tx0.max(0), ty0.max(0), tx1.max(0), ty1.max(0)), max_exp: 1 << 30, touched: std::ptr::null_mut(), cost: None }
}

#[test]
fn astar_matches_reference_and_fields_agree() {
    use gridroute::astar::{astar, astar_ref, Outcome};
    let mut rng = Rng(0x1234_5678_9abc_def1);
    let mc: Vec<f32> = (0..4 * 8).map(|k| if k % 8 < 4 { 1.0 + (k / 8) as f32 * 0.3 } else { 1.5 }).collect();
    let lok = [1u8, 1, 1, 1];
    for _ in 0..40 {
        let (nl, h, w) = (1 + rng.below(4) as usize, 20 + rng.below(80) as usize, 20 + rng.below(120) as usize);
        let (b, v, t, s) = random_problem(&mut rng, nl, h, w);
        if !t.iter().any(|&x| x != 0) {
            continue;
        }
        let p = problem(nl, h, w, &b, &v, &t, &s, &mc, &lok[..nl]);
        let (a, _) = astar(&p);
        let (r, _) = astar_ref(&p);
        match (a, r) {
            (Outcome::Found(x), Outcome::Found(y)) => assert_eq!(x, y),
            (Outcome::NotFound, Outcome::NotFound) => {}
            _ => panic!("optimised and reference search disagree"),
        }
        let mut fc = vec![0f32; nl * h * w];
        gridroute::field::field_cpu(&p, &mut fc);
        #[cfg(all(target_os = "macos", feature = "metal"))]
        {
            let mut fg = vec![0f32; nl * h * w];
            if gridroute::gpu::field_metal(&p, &mut fg).is_some() {
                for (x, y) in fc.iter().zip(&fg) {
                    assert!((x.is_infinite() && y.is_infinite()) || (x - y).abs() <= 1e-3 * x.max(1.0), "{x} vs {y}");
                }
            }
        }
    }
}
