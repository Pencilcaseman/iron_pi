use crate::{
    fact::{PrimeFactorSieve, PrimeFactors},
    util,
};

const PAR_THRESHOLD: u64 = 8192;

#[derive(Debug)]
pub struct NumFac {
    pub num: flint3_sys::fmpz_t,
    pub fac: PrimeFactors,
}

pub fn binary_split_work_stealing(
    a: u64,
    b: u64,
    sieve: &PrimeFactorSieve,
    thread_pool: &rayon::ThreadPool,
) -> (NumFac, NumFac, NumFac) {
    if b - a == 1 {
        // Base case

        // P(a, b) = (6 * a - 1) * (6 * a - 5) * (2 * a - 1)
        let mut p_fac = PrimeFactors::new(sieve, 6 * a - 1);
        p_fac = p_fac * PrimeFactors::new(sieve, 6 * a - 5);
        p_fac = p_fac * PrimeFactors::new(sieve, 2 * a - 1);
        p_fac.neg = true;

        let p_num = p_fac.to_int();
        let factored_p = NumFac { num: p_num, fac: p_fac.clone() };

        // Q(a, b) = (640320^3) / 24 * a^3
        let mut q_fac = PrimeFactors::new_with_pow(sieve, a, 3);
        q_fac = q_fac * PrimeFactors::new_with_pow(sieve, 640320, 2);
        q_fac = q_fac * PrimeFactors::new(sieve, 640320 / 24);

        let q_num = q_fac.to_int();
        let factored_q = NumFac { num: q_num, fac: q_fac };

        unsafe {
            let mut r_num = util::new_fmpz_with(545_140_134);

            flint3_sys::fmpz_mul_ui(&mut r_num[0], &r_num[0], a);
            flint3_sys::fmpz_add_ui(&mut r_num[0], &r_num[0], 13_591_409);
            flint3_sys::fmpz_mul(&mut r_num[0], &r_num[0], &p_num[0]);

            let factored_r = NumFac { num: r_num, fac: p_fac };

            (factored_p, factored_q, factored_r)
        }
    } else {
        // Recursive case

        let mid = (a + b) / 2;

        let ((mut p1, mut q1, mut r1), (mut p2, mut q2, mut r2)) = thread_pool
            .join(
                || binary_split_work_stealing(a, mid, sieve, thread_pool),
                || binary_split_work_stealing(mid, b, sieve, thread_pool),
            );

        unsafe {
            let mut p_fac = &p1.fac * &p2.fac;
            let mut q_fac = &q1.fac * &p2.fac;

            let mut p_num = util::new_fmpz();
            flint3_sys::fmpz_mul(&mut p_num[0], &p1.num[0], &p2.num[0]);

            let mut q_num = util::new_fmpz();
            flint3_sys::fmpz_mul(&mut q_num[0], &q1.num[0], &q2.num[0]);

            let mut r_num = util::new_fmpz();

            flint3_sys::fmpz_fmma(
                &mut r_num[0],
                &q2.num[0],
                &r1.num[0],
                &p1.num[0],
                &r2.num[0],
            );

            let gcd_fac = q1.fac.gcd(&r1.fac);

            let gcd_num = gcd_fac.to_int();

            flint3_sys::fmpz_divexact(&mut p_num[0], &p_num[0], &gcd_num[0]);
            p_fac.div_exact_mut(&gcd_fac);

            flint3_sys::fmpz_divexact(&mut q_num[0], &q_num[0], &gcd_num[0]);
            q_fac.div_exact_mut(&gcd_fac);

            flint3_sys::fmpz_divexact(&mut r_num[0], &r_num[0], &gcd_num[0]);

            let r_fac = p1.fac.gcd(&r2.fac);

            let p = NumFac { num: p_num, fac: p_fac };
            let q = NumFac { num: q_num, fac: q_fac };
            let r = NumFac { num: r_num, fac: r_fac };

            flint3_sys::fmpz_clear(&mut p1.num[0]);
            flint3_sys::fmpz_clear(&mut p2.num[0]);
            flint3_sys::fmpz_clear(&mut q1.num[0]);
            flint3_sys::fmpz_clear(&mut q2.num[0]);
            flint3_sys::fmpz_clear(&mut r1.num[0]);
            flint3_sys::fmpz_clear(&mut r2.num[0]);

            (p, q, r)
        }
    }
}

pub fn binary_split_2(
    a: u64,
    b: u64,
    q_partial: &flint3_sys::fmpz_t,
    thread_pool: &rayon::ThreadPool,
) -> (flint3_sys::fmpz_t, flint3_sys::fmpz_t, flint3_sys::fmpz_t) {
    unsafe {
        if b - a == 1 {
            // P(a, a + 1) = -(6a - 1)(6a - 5)(2a - 1)
            let mut p = util::new_fmpz_with(6 * a - 1);
            flint3_sys::fmpz_mul_ui(&mut p[0], &p[0], 6 * a - 5);
            flint3_sys::fmpz_mul_ui(&mut p[0], &p[0], 2 * a - 1);
            flint3_sys::fmpz_neg(&mut p[0], &p[0]);

            // Q(a, a + 1) = a^3 * 640320^3 / 24
            let mut q = util::new_fmpz_with(a * a); // Safe for digits < 2e+20
            flint3_sys::fmpz_mul(&mut q[0], &q[0], &q_partial[0]);
            flint3_sys::fmpz_mul_ui(&mut q[0], &q[0], a);

            // R(a, a + 1) = P(a, a + 1) * (545140134 * a + 13591409)
            let mut r = util::new_fmpz_with(545_140_134);
            flint3_sys::fmpz_mul_ui(&mut r[0], &r[0], a);
            flint3_sys::fmpz_add_ui(&mut r[0], &r[0], 13_591_409);
            flint3_sys::fmpz_mul(&mut r[0], &r[0], &p[0]);

            (p, q, r)
        } else {
            let mid = (a + b) / 2;

            let ((mut p1, mut q1, mut r1), (mut p2, mut q2, mut r2)) =
                if b - a < PAR_THRESHOLD {
                    (
                        binary_split_2(a, mid, q_partial, thread_pool),
                        binary_split_2(mid, b, q_partial, thread_pool),
                    )
                } else {
                    thread_pool.join(
                        || binary_split_2(a, mid, q_partial, thread_pool),
                        || binary_split_2(mid, b, q_partial, thread_pool),
                    )
                };

            flint3_sys::fmpz_mul(&mut p2[0], &p1[0], &p2[0]);
            flint3_sys::fmpz_mul(&mut q1[0], &q1[0], &q2[0]);
            flint3_sys::fmpz_fmma(&mut r2[0], &q2[0], &r1[0], &p1[0], &r2[0]);

            flint3_sys::fmpz_clear(&mut p1[0]);
            flint3_sys::fmpz_clear(&mut q2[0]);
            flint3_sys::fmpz_clear(&mut r1[0]);

            (p2, q1, r2)
        }
    }
}

pub fn binary_split(
    a: u64,
    b: u64,
    thread_pool: &rayon::ThreadPool,
) -> (flint3_sys::fmpz_t, flint3_sys::fmpz_t, flint3_sys::fmpz_t) {
    unsafe {
        let mut q_partial = util::new_fmpz_with(640320);
        flint3_sys::fmpz_mul_ui(&mut q_partial[0], &q_partial[0], 640320);
        flint3_sys::fmpz_mul_ui(&mut q_partial[0], &q_partial[0], 640320 / 24);

        let res = binary_split_2(a, b, &q_partial, thread_pool);

        flint3_sys::fmpz_clear(&mut q_partial[0]);

        res
    }
}
