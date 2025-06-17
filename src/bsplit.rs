use crate::{
    fact::{PrimeFactorSieve, PrimeFactors},
    util,
};

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

        // let mut r_num = Integer::from(545_140_134);
        unsafe {
            // (545_140_134 * a + 13_591_409) * p_num

            // let mut r_num: flint3_sys::fmpz_t =
            //     std::mem::MaybeUninit::uninit().assume_init();
            // flint3_sys::fmpz_init_set_ui(&mut r_num[0], 545_140_134);

            let mut r_num = util::new_fmpz_with(545_140_134);

            // r_num.mul_assign(a);
            // r_num.add_assign(13_591_409);
            // r_num.mul_assign(&p_num);

            // r = (r * a + 13_591_409) * p_num

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

        // let mut p_fac = &p1.fac * &p2.fac;
        // let mut p_num = (&p1.num * &p2.num).complete();
        //
        // let mut q_fac = &q1.fac * &p2.fac;
        // let mut q_num = (&q1.num * &q2.num).complete();

        // let mut r_num =
        //     (&q2.num * &r1.num).complete() + (&p1.num * &r2.num).complete();

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

            // p_num.div_exact_mut(&gcd_num);
            // p_fac.div_exact_mut(&gcd_fac);
            //
            // q_num.div_exact_mut(&gcd_num);
            // q_fac.div_exact_mut(&gcd_fac);
            //
            // r_num.div_exact_mut(&gcd_num);

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
