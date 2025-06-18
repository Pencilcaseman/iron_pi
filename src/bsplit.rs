use crate::util;

fn binary_split_2(
    a: u64,
    b: u64,
    q_partial: &flint3_sys::fmpz_t,
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

            let (mut p1, mut q1, mut r1) = binary_split_2(a, mid, q_partial);
            let (mut p2, mut q2, mut r2) = binary_split_2(mid, b, q_partial);

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

pub struct SharedArb(flint3_sys::arb_t);

unsafe impl Send for SharedArb {}
unsafe impl Sync for SharedArb {}

fn binary_split_2_arb(
    a: u64,
    b: u64,
    q_partial_fmpz: &flint3_sys::fmpz_t,
    thread_pool: &rayon::ThreadPool,
    depth: usize,
    max_par_depth: usize,
    prec: i64,
) -> (SharedArb, SharedArb, SharedArb) {
    unsafe {
        if b - a < 4 {
            let mut p = util::new_arb();
            let mut q = util::new_arb();
            let mut r = util::new_arb();

            let (mut p_tmp, mut q_tmp, mut r_tmp) =
                binary_split_2(a, b, q_partial_fmpz);

            flint3_sys::arb_set_fmpz(&mut p[0], &p_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut q[0], &q_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut r[0], &r_tmp[0]);

            flint3_sys::fmpz_clear(&mut p_tmp[0]);
            flint3_sys::fmpz_clear(&mut q_tmp[0]);
            flint3_sys::fmpz_clear(&mut r_tmp[0]);

            (SharedArb(p), SharedArb(q), SharedArb(r))
        } else {
            let mid = (a + b) / 2;

            let ((mut p1, mut q1, mut r1), (mut p2, mut q2, mut r2)) =
                if depth <= max_par_depth {
                    thread_pool.join(
                        || {
                            binary_split_2_arb(
                                a,
                                mid,
                                q_partial_fmpz,
                                thread_pool,
                                depth + 1,
                                max_par_depth,
                                prec,
                            )
                        },
                        || {
                            binary_split_2_arb(
                                mid,
                                b,
                                q_partial_fmpz,
                                thread_pool,
                                depth + 1,
                                max_par_depth,
                                prec,
                            )
                        },
                    )
                } else {
                    (
                        binary_split_2_arb(
                            a,
                            mid,
                            q_partial_fmpz,
                            thread_pool,
                            depth + 1,
                            max_par_depth,
                            prec,
                        ),
                        binary_split_2_arb(
                            mid,
                            b,
                            q_partial_fmpz,
                            thread_pool,
                            depth + 1,
                            max_par_depth,
                            prec,
                        ),
                    )
                };

            flint3_sys::arb_mul(&mut p2.0[0], &p1.0[0], &p2.0[0], prec);
            flint3_sys::arb_mul(&mut q1.0[0], &q1.0[0], &q2.0[0], prec);

            let mut temp1 = util::new_arb();
            let mut temp2 = util::new_arb();

            flint3_sys::arb_init(&mut temp1[0]);
            flint3_sys::arb_init(&mut temp2[0]);

            flint3_sys::arb_mul(&mut temp1[0], &q2.0[0], &r1.0[0], prec);
            flint3_sys::arb_mul(&mut temp2[0], &p1.0[0], &r2.0[0], prec);
            flint3_sys::arb_add(&mut r2.0[0], &temp1[0], &temp2[0], prec);

            flint3_sys::arb_clear(&mut temp1[0]);
            flint3_sys::arb_clear(&mut temp2[0]);

            flint3_sys::arb_clear(&mut p1.0[0]);
            flint3_sys::arb_clear(&mut q2.0[0]);
            flint3_sys::arb_clear(&mut r1.0[0]);

            (p2, q1, r2)
        }
    }
}

pub fn binary_split_arb(
    a: u64,
    b: u64,
    thread_pool: &rayon::ThreadPool,
    max_par_depth: usize,
    prec: i64,
) -> (flint3_sys::arb_t, flint3_sys::arb_t, flint3_sys::arb_t) {
    unsafe {
        let mut q_partial = util::new_fmpz_with(640320);
        flint3_sys::fmpz_mul_ui(&mut q_partial[0], &q_partial[0], 640320);
        flint3_sys::fmpz_mul_ui(&mut q_partial[0], &q_partial[0], 640320 / 24);

        let (SharedArb(p), SharedArb(q), SharedArb(r)) = binary_split_2_arb(
            a,
            b,
            &q_partial,
            thread_pool,
            0,
            max_par_depth,
            prec,
        );

        flint3_sys::fmpz_clear(&mut q_partial[0]);

        (p, q, r)
    }
}
