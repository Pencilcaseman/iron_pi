use rug::{ops::NegAssign, Integer};

const INITIAL_CAPACITY: usize = 512;

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct FactorSieveElement {
    pub base: usize,
    pub exponent: usize,
    pub next: usize,
}

#[derive(Debug)]
pub struct PrimeFactorSieve {
    pub max: usize,
    pub sieve: Vec<FactorSieveElement>,
}

impl PrimeFactorSieve {
    pub fn new(max: usize) -> Self {
        let sqrt_max = max.isqrt();

        let mut sieve = vec![FactorSieveElement::default(); max / 2];

        sieve[0].base = 2;
        sieve[0].exponent = 1;

        // Step by all odd numbers
        for i in (3..max).step_by(2) {
            if sieve[i / 2].base == 0 {
                // Value is prime
                sieve[i / 2].base = i;
                sieve[i / 2].exponent = 1;
                sieve[i / 2].next = 0;

                if i < sqrt_max {
                    // Fill in powers of i
                    let mut j = i * i;
                    let mut k = i / 2;

                    while j < max {
                        if sieve[j / 2].base == 0 {
                            sieve[j / 2].base = i;
                            if sieve[k].base == i {
                                sieve[j / 2].exponent = sieve[k].exponent + 1;
                                sieve[j / 2].next = sieve[k].next;
                            } else {
                                sieve[j / 2].exponent = 1;
                                sieve[j / 2].next = k;
                            }
                        }

                        j += i * 2; // Skip even multiples
                        k += 1; // Counts for powers of i
                    }
                }
            }
        }

        Self { max, sieve }
    }
}

#[derive(Debug, Clone)]
pub struct PrimeFactors {
    pub neg: bool,
    pub factors: Vec<(usize, usize)>,
}

impl PrimeFactors {
    pub fn empty() -> Self {
        Self { neg: false, factors: Vec::with_capacity(INITIAL_CAPACITY) }
    }

    pub fn new(sieve: &PrimeFactorSieve, mut value: usize) -> Self {
        #[cold]
        fn out_of_range() {
            panic!("Value is out of range for the sieve");
        }

        let mut factors = Vec::with_capacity(INITIAL_CAPACITY);

        // Handle factors of two
        let pow_2 = value.trailing_zeros() as usize;
        value >>= pow_2;

        if pow_2 > 0 {
            factors.push((2, pow_2));
        }

        let mut base = value / 2;
        while base > 0 {
            if base > sieve.max {
                out_of_range();
            }

            let element = sieve.sieve[base];
            factors.push((element.base, element.exponent));
            base = element.next;
        }

        Self { neg: false, factors }
    }

    pub fn new_with_pow(
        sieve: &PrimeFactorSieve,
        value: usize,
        pow: usize,
    ) -> Self {
        let mut tmp = Self::new(sieve, value);
        for p in tmp.factors.iter_mut() {
            p.1 *= pow;
        }

        tmp
    }

    pub fn to_int(&self) -> Integer {
        fn split_mul(
            factors: &[(usize, usize)],
            a: usize,
            b: usize,
        ) -> rug::Integer {
            if b - a < 32 {
                // Repeated multiplication
                let mut result = rug::Integer::from(1);

                factors[a..b].iter().for_each(|(base, pow)| {
                    for _ in 0..*pow {
                        result *= base;
                    }
                });

                result
            } else {
                let mid = (a + b) / 2;
                let lhs = split_mul(factors, a, mid);
                let rhs = split_mul(factors, mid, b);
                lhs * rhs
            }
        }

        let mut res = split_mul(&self.factors, 0, self.factors.len());

        if self.neg {
            res.neg_assign();
        }

        res
    }

    pub fn gcd(&self, rhs: &Self) -> Self {
        let mut result = Vec::with_capacity(self.factors.len());
        let mut i = 0;
        let mut j = 0;

        while i < self.factors.len() && j < rhs.factors.len() {
            match self.factors[i].0 {
                base if base == rhs.factors[j].0 => {
                    result
                        .push((base, self.factors[i].1.min(rhs.factors[j].1)));
                    i += 1;
                    j += 1;
                }
                base if base < rhs.factors[j].0 => {
                    i += 1;
                }
                _ => {
                    j += 1;
                }
            }
        }

        // Remove zero exponents
        result.retain(|term| term.1 != 0);

        Self {
            neg: false, // self.neg && rhs.neg,
            factors: result,
        }
    }

    pub fn remove_gcd(
        sieve: &PrimeFactorSieve,
        (lhs_factors, lhs_int): (&mut Self, &mut Integer),
        (rhs_factors, rhs_int): (&mut Self, &mut Integer),
    ) {
        let mut i = 0;
        let mut j = 0;
        let lhs = &mut lhs_factors.factors;
        let rhs = &mut rhs_factors.factors;

        let mut to_remove = PrimeFactors::new(sieve, 1);

        while i < lhs.len() && j < rhs.len() {
            match lhs[i].0 {
                base if base == rhs[j].0 => {
                    let gcd_pow = lhs[i].1.min(rhs[j].1);
                    lhs[i].1 -= gcd_pow;
                    rhs[j].1 -= gcd_pow;

                    to_remove = to_remove
                        * PrimeFactors::new_with_pow(sieve, base, gcd_pow);

                    i += 1;
                    j += 1;
                }
                base if base < rhs[j].0 => {
                    i += 1;
                }
                _ => {
                    j += 1;
                }
            }
        }

        // Remove zero exponents
        lhs.retain(|term| term.1 != 0);
        rhs.retain(|term| term.1 != 0);

        let to_remove = to_remove.to_int();
        println!("Removing GCD: {to_remove}");

        lhs_int.div_exact_mut(&to_remove);
        rhs_int.div_exact_mut(&to_remove);

        if lhs_factors.neg && rhs_factors.neg {
            lhs_factors.neg = false;
            rhs_factors.neg = false;
            lhs_int.neg_assign();
            rhs_int.neg_assign();
        }
    }

    /// # Safety
    ///
    /// rhs must be a factor of self
    pub unsafe fn div_exact_mut(&mut self, rhs: &Self) {
        let mut i = 0;
        let mut j = 0;

        while i < self.factors.len() && j < rhs.factors.len() {
            match self.factors[i].0 {
                base if base == rhs.factors[j].0 => {
                    self.factors[i].1 -= rhs.factors[j].1;
                    i += 1;
                    j += 1;
                }
                base if base < rhs.factors[j].0 => {
                    i += 1;
                }
                _ => {
                    j += 1;
                }
            }
        }
    }

    pub fn mul_into(out: &mut Self, lhs: &Self, rhs: &Self) {
        out.factors.clear(); // Clear out the existing factors. Pushes are now free

        let mut i = 0;
        let mut j = 0;

        while i < lhs.factors.len() && j < rhs.factors.len() {
            match lhs.factors[i].0 {
                base if base == rhs.factors[j].0 => {
                    out.factors
                        .push((base, lhs.factors[i].1 + rhs.factors[j].1));
                    i += 1;
                    j += 1;
                }
                base if base < rhs.factors[j].0 => {
                    out.factors.push(lhs.factors[i]);
                    i += 1;
                }
                _ => {
                    out.factors.push(rhs.factors[j]);
                    j += 1;
                }
            }
        }

        while i < lhs.factors.len() {
            out.factors.push(lhs.factors[i]);
            i += 1;
        }

        while j < rhs.factors.len() {
            out.factors.push(rhs.factors[j]);
            j += 1;
        }
    }
}

impl std::ops::Mul<PrimeFactors> for PrimeFactors {
    type Output = PrimeFactors;

    fn mul(self, rhs: Self) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        let mut result = Self {
            neg: self.neg ^ rhs.neg,
            factors: Vec::with_capacity(
                self.factors.len().max(rhs.factors.len()),
            ),
        };

        Self::mul_into(&mut result, &self, &rhs);
        result
    }
}

impl std::ops::Mul<&PrimeFactors> for PrimeFactors {
    type Output = PrimeFactors;

    fn mul(self, rhs: &Self) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        let mut result = Self {
            neg: self.neg ^ rhs.neg,
            factors: Vec::with_capacity(
                self.factors.len().max(rhs.factors.len()),
            ),
        };

        Self::mul_into(&mut result, &self, rhs);
        result
    }
}

impl<'a> std::ops::Mul<PrimeFactors> for &'a PrimeFactors {
    type Output = PrimeFactors;

    fn mul(self, rhs: PrimeFactors) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        let mut result = PrimeFactors {
            neg: self.neg ^ rhs.neg,
            factors: Vec::with_capacity(
                self.factors.len().max(rhs.factors.len()),
            ),
        };

        PrimeFactors::mul_into(&mut result, self, &rhs);
        result
    }
}

impl<'a> std::ops::Mul<&'a PrimeFactors> for &'a PrimeFactors {
    type Output = PrimeFactors;

    fn mul(self, rhs: &'a PrimeFactors) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        let mut result = PrimeFactors {
            neg: self.neg ^ rhs.neg,
            factors: Vec::with_capacity(
                self.factors.len().max(rhs.factors.len()),
            ),
        };

        PrimeFactors::mul_into(&mut result, self, rhs);
        result
    }
}
