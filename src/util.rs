pub fn new_fmpz() -> flint3_sys::fmpz_t {
    unsafe {
        let mut fmpz = std::mem::MaybeUninit::<flint3_sys::fmpz_t>::uninit();
        flint3_sys::fmpz_init(&mut (*fmpz.as_mut_ptr())[0]);
        fmpz.assume_init()
    }
}

pub fn new_fmpz_with(value: u64) -> flint3_sys::fmpz_t {
    unsafe {
        let mut fmpz = std::mem::MaybeUninit::<flint3_sys::fmpz_t>::uninit();
        flint3_sys::fmpz_init_set_ui(&mut (*fmpz.as_mut_ptr())[0], value);
        fmpz.assume_init()
    }
}

pub fn new_arb() -> flint3_sys::arb_t {
    unsafe {
        let mut fmpz = std::mem::MaybeUninit::<flint3_sys::arb_t>::uninit();
        flint3_sys::arb_init(&mut (*fmpz.as_mut_ptr())[0]);
        fmpz.assume_init()
    }
}

pub fn new_arb_with(value: f64) -> flint3_sys::arb_t {
    let mut arb = new_arb();

    unsafe {
        flint3_sys::arb_set_d(&mut arb[0], value);
    }

    arb
}

pub fn new_arb_with_i64(value: i64) -> flint3_sys::arb_t {
    let mut arb = new_arb();

    unsafe {
        flint3_sys::arb_set_si(&mut arb[0], value);
    }

    arb
}

pub fn new_arb_with_u64(value: u64) -> flint3_sys::arb_t {
    let mut arb = new_arb();

    unsafe {
        flint3_sys::arb_set_ui(&mut arb[0], value);
    }

    arb
}

pub fn new_arb_with_fmpz(value: &flint3_sys::fmpz_t) -> flint3_sys::arb_t {
    let mut arb = new_arb();

    unsafe {
        flint3_sys::arb_set_fmpz(&mut arb[0], &value[0]);
    }

    arb
}

pub fn new_poly() -> flint3_sys::fmpz_poly_t {
    unsafe {
        let mut poly =
            std::mem::MaybeUninit::<flint3_sys::fmpz_poly_t>::uninit();
        flint3_sys::fmpz_poly_init(&mut (*poly.as_mut_ptr())[0]);
        poly.assume_init()
    }
}

pub fn new_hypgeom() -> flint3_sys::hypgeom_t {
    unsafe {
        let mut hypgeom =
            std::mem::MaybeUninit::<flint3_sys::hypgeom_t>::uninit();
        flint3_sys::hypgeom_init(&mut (*hypgeom.as_mut_ptr())[0]);
        hypgeom.assume_init()
    }
}
