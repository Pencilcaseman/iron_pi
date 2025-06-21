use crate::util;

pub struct SharedArb(flint3_sys::arb_t);
pub struct SharedPoly(flint3_sys::fmpz_poly_t);

unsafe impl Send for SharedArb {}
unsafe impl Sync for SharedArb {}

unsafe impl Send for SharedPoly {}
unsafe impl Sync for SharedPoly {}

#[allow(clippy::too_many_arguments)]
fn binary_split_fmpz(
    a: u64,
    b: u64,
    p: &mut flint3_sys::fmpz_t,
    q: &mut flint3_sys::fmpz_t,
    r: &mut flint3_sys::fmpz_t,
    poly_p: &SharedPoly,
    poly_q: &SharedPoly,
    poly_r: &SharedPoly,
) {
    unsafe {
        if b - a == 1 {
            let mut tmp_a = util::new_fmpz_with(a);

            // P(a, a + 1) = -(6a - 1)(6a - 5)(2a - 1)
            flint3_sys::fmpz_poly_evaluate_divconquer_fmpz(
                &mut p[0],
                &poly_p.0[0],
                &tmp_a[0],
            );

            // Q(a, a + 1) = a^3 * 640320^3 / 24
            flint3_sys::fmpz_poly_evaluate_divconquer_fmpz(
                &mut q[0],
                &poly_q.0[0],
                &tmp_a[0],
            );

            // R(a, a + 1) = P(a, a + 1) * (545140134 * a + 13591409)
            flint3_sys::fmpz_poly_evaluate_divconquer_fmpz(
                &mut r[0],
                &poly_r.0[0],
                &tmp_a[0],
            );
            flint3_sys::fmpz_mul(&mut r[0], &r[0], &p[0]);

            flint3_sys::fmpz_clear(&mut tmp_a[0]);
        } else {
            let mid = (a + b) / 2;

            let mut p2 = util::new_fmpz();
            let mut q2 = util::new_fmpz();
            let mut r2 = util::new_fmpz();

            binary_split_fmpz(a, mid, p, q, r, poly_p, poly_q, poly_r);
            binary_split_fmpz(
                mid, b, &mut p2, &mut q2, &mut r2, poly_p, poly_q, poly_r,
            );

            // R(a, b) = R(a, m) * Q(m, b) + P(a, m) * R(m, b)
            flint3_sys::fmpz_fmma(&mut r[0], &q2[0], &r[0], &p[0], &r2[0]);

            // Q(a, b) = Q(a, m) * Q(m, b)
            flint3_sys::fmpz_mul(&mut q[0], &q[0], &q2[0]);

            // P(a, b) = P(a, m) * P(m, b)
            flint3_sys::fmpz_mul(&mut p[0], &p[0], &p2[0]);

            flint3_sys::fmpz_clear(&mut p2[0]);
            flint3_sys::fmpz_clear(&mut q2[0]);
            flint3_sys::fmpz_clear(&mut r2[0]);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn binary_split_arb_single(
    a: u64,
    b: u64,
    p: &mut SharedArb,
    q: &mut SharedArb,
    r: &mut SharedArb,
    poly_p: &SharedPoly,
    poly_q: &SharedPoly,
    poly_r: &SharedPoly,
    prec: i64,
) {
    unsafe {
        if b - a < 4 {
            let mut p_tmp = util::new_fmpz();
            let mut q_tmp = util::new_fmpz();
            let mut r_tmp = util::new_fmpz();

            binary_split_fmpz(
                a, b, &mut p_tmp, &mut q_tmp, &mut r_tmp, poly_p, poly_q,
                poly_r,
            );

            flint3_sys::arb_set_fmpz(&mut p.0[0], &p_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut q.0[0], &q_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut r.0[0], &r_tmp[0]);

            flint3_sys::fmpz_clear(&mut p_tmp[0]);
            flint3_sys::fmpz_clear(&mut q_tmp[0]);
            flint3_sys::fmpz_clear(&mut r_tmp[0]);
        } else {
            let mid = (a + b) / 2;

            let mut p2 = SharedArb(util::new_arb());
            let mut q2 = SharedArb(util::new_arb());
            let mut r2 = SharedArb(util::new_arb());

            binary_split_arb_single(
                a, mid, p, q, r, poly_p, poly_q, poly_r, prec,
            );
            binary_split_arb_single(
                mid, b, &mut p2, &mut q2, &mut r2, poly_p, poly_q, poly_r, prec,
            );

            // R(a, b) = Q(m, b) * R(a, m) + P(a, m) * R(m, b)
            flint3_sys::arb_mul(&mut r.0[0], &r.0[0], &q2.0[0], prec);
            flint3_sys::arb_addmul(&mut r.0[0], &p.0[0], &r2.0[0], prec);

            // Q(a, b) = Q(a, m) * Q(m, b)
            flint3_sys::arb_mul(&mut q.0[0], &q.0[0], &q2.0[0], prec);

            // if depth > 0 {
            // P(a, b) = P(a, m) * P(m, b)
            flint3_sys::arb_mul(&mut p.0[0], &p.0[0], &p2.0[0], prec);
            // }

            // Cleanup
            flint3_sys::arb_clear(&mut p2.0[0]);
            flint3_sys::arb_clear(&mut q2.0[0]);
            flint3_sys::arb_clear(&mut r2.0[0]);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn binary_split_arb(
    a: u64,
    b: u64,
    p: &mut SharedArb,
    q: &mut SharedArb,
    r: &mut SharedArb,
    poly_p: &SharedPoly,
    poly_q: &SharedPoly,
    poly_r: &SharedPoly,
    thread_pool: &rayon::ThreadPool,
    depth: usize,
    max_par_depth: usize,
    prec: i64,
) {
    unsafe {
        if b - a < 4 {
            let mut p_tmp = util::new_fmpz();
            let mut q_tmp = util::new_fmpz();
            let mut r_tmp = util::new_fmpz();

            binary_split_fmpz(
                a, b, &mut p_tmp, &mut q_tmp, &mut r_tmp, poly_p, poly_q,
                poly_r,
            );

            flint3_sys::arb_set_fmpz(&mut p.0[0], &p_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut q.0[0], &q_tmp[0]);
            flint3_sys::arb_set_fmpz(&mut r.0[0], &r_tmp[0]);

            flint3_sys::fmpz_clear(&mut p_tmp[0]);
            flint3_sys::fmpz_clear(&mut q_tmp[0]);
            flint3_sys::fmpz_clear(&mut r_tmp[0]);
        } else {
            let mid = (a + b) / 2;

            let mut p2 = SharedArb(util::new_arb());
            let mut q2 = SharedArb(util::new_arb());
            let mut r2 = SharedArb(util::new_arb());

            if thread_pool.current_num_threads() > 1 && depth <= max_par_depth {
                thread_pool.join(
                    || {
                        binary_split_arb(
                            a,
                            mid,
                            p,
                            q,
                            r,
                            poly_p,
                            poly_q,
                            poly_r,
                            thread_pool,
                            depth + 1,
                            max_par_depth,
                            prec,
                        )
                    },
                    || {
                        binary_split_arb(
                            mid,
                            b,
                            &mut p2,
                            &mut q2,
                            &mut r2,
                            poly_p,
                            poly_q,
                            poly_r,
                            thread_pool,
                            depth + 1,
                            max_par_depth,
                            prec,
                        )
                    },
                )
            } else {
                (
                    binary_split_arb_single(
                        a, mid, p, q, r, poly_p, poly_q, poly_r, prec,
                    ),
                    binary_split_arb_single(
                        mid, b, &mut p2, &mut q2, &mut r2, poly_p, poly_q,
                        poly_r, prec,
                    ),
                )
            };

            // R(a, b) = Q(m, b) * R(a, m) + P(a, m) * R(m, b)
            flint3_sys::arb_mul(&mut r.0[0], &r.0[0], &q2.0[0], prec);
            flint3_sys::arb_addmul(&mut r.0[0], &p.0[0], &r2.0[0], prec);

            // Q(a, b) = Q(a, m) * Q(m, b)
            flint3_sys::arb_mul(&mut q.0[0], &q.0[0], &q2.0[0], prec);

            if depth > 0 {
                // P(a, b) = P(a, m) * P(m, b)
                flint3_sys::arb_mul(&mut p.0[0], &p.0[0], &p2.0[0], prec);
            }

            // Cleanup
            flint3_sys::arb_clear(&mut p2.0[0]);
            flint3_sys::arb_clear(&mut q2.0[0]);
            flint3_sys::arb_clear(&mut r2.0[0]);
        }
    }
}

pub fn binary_split(
    a: u64,
    b: u64,
    thread_pool: &rayon::ThreadPool,
    max_par_depth: usize,
    prec: i64,
) -> (flint3_sys::arb_t, flint3_sys::arb_t, flint3_sys::arb_t) {
    unsafe {
        let mut poly_p = SharedPoly(util::new_poly());
        let mut poly_q = SharedPoly(util::new_poly());
        let mut poly_r = SharedPoly(util::new_poly());

        let poly_p_str = std::ffi::CString::new("4  5 -46 108 -72").unwrap();

        // 10939058860032000 = (640320^3)/24
        let poly_q_str =
            std::ffi::CString::new("4  0 0 0 10939058860032000").unwrap();

        let poly_r_str =
            std::ffi::CString::new("2  13591409 545140134").unwrap();

        flint3_sys::fmpz_poly_set_str(&mut poly_p.0[0], poly_p_str.as_ptr());
        flint3_sys::fmpz_poly_set_str(&mut poly_q.0[0], poly_q_str.as_ptr());
        flint3_sys::fmpz_poly_set_str(&mut poly_r.0[0], poly_r_str.as_ptr());

        let mut p = SharedArb(util::new_arb());
        let mut q = SharedArb(util::new_arb());
        let mut r = SharedArb(util::new_arb());

        binary_split_arb(
            a,
            b,
            &mut p,
            &mut q,
            &mut r,
            &poly_p,
            &poly_q,
            &poly_r,
            thread_pool,
            0,
            max_par_depth,
            prec,
        );

        flint3_sys::fmpz_poly_clear(&mut poly_p.0[0]);
        flint3_sys::fmpz_poly_clear(&mut poly_q.0[0]);
        flint3_sys::fmpz_poly_clear(&mut poly_r.0[0]);

        (p.0, q.0, r.0)
    }
}
