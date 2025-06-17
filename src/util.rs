pub fn new_fmpz() -> flint3_sys::fmpz_t {
    unsafe {
        let mut p_num = std::mem::MaybeUninit::<flint3_sys::fmpz_t>::uninit();
        flint3_sys::fmpz_init(&mut (*p_num.as_mut_ptr())[0]);
        p_num.assume_init()
    }
}

pub fn new_fmpz_with(value: u64) -> flint3_sys::fmpz_t {
    unsafe {
        let mut p_num = std::mem::MaybeUninit::<flint3_sys::fmpz_t>::uninit();
        flint3_sys::fmpz_init_set_ui(&mut (*p_num.as_mut_ptr())[0], value);
        p_num.assume_init()
    }
}

pub fn new_arb() -> flint3_sys::arb_t {
    unsafe {
        let mut p_num = std::mem::MaybeUninit::<flint3_sys::arb_t>::uninit();
        flint3_sys::arb_init(&mut (*p_num.as_mut_ptr())[0]);
        p_num.assume_init()
    }
}

pub fn new_arb_with(value: f64) -> flint3_sys::arb_t {
    let mut arb = new_arb();

    unsafe {
        flint3_sys::arb_set_d(&mut arb[0], value);
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
