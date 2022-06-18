//! A `FieldElement` represents an element of the finite field
//! modulo `2^252 + 27742317777372353535851937790883648493`.
//!
//! The `FieldElement` type is an alias for one of the backend
//! implementations.
//!
//! `ConstantTimeEq` and `PartialEq` traits have been implemented
//! here since they will be the samme across all of the different backends.
//!
//! # Examples
//! ```rust
//! use zerocaf::field::FieldElement;
//! use zerocaf::traits::ops::*;
//! use zerocaf::constants::EDWARDS_D;
//!
//! use subtle::Choice;
//! use rand::rngs::OsRng;
//!
//! // You can create a FieldElement from a byte-array as follows:
//! let a = FieldElement::from_bytes(&[0u8;32]);
//!
//! // You ca also create a FieldElement from an uint type as follows:
//! let b = FieldElement::from(126296u128);
//! let c = FieldElement::from(126297u64);
//!
//! // You can create random FieldElements by calling:
//! let rand = FieldElement::random(&mut OsRng);
//!
//! // The last way of creating a FieldElement it by calling the
//! // constructor. THIS IS NOT RECOMMENDED since NO checks about
//! // the correctness of the input will be done at all.
//! // It can be done as follows:
//! let d: FieldElement = FieldElement([0, 1, 0, 0, 0]); // d = 2^52.
//! assert!(d == FieldElement::two_pow_k(52u64));
//!
//! // All of the basuc modular operations are implemented
//! // for FieldElement type:
//! let mut res = &a + &b; // Performs a + b (mod l).
//! res = a - b; // Performs a - b (mod l).
//! res = a * b; // Performs a * b (mod l).
//! res = a.square(); // Performs a^2 (mod l).
//! res = -&a; // Performs Negation over the modulo l.
//! res = a.pow(&b); // Performs Modular exponentiation.
//! res = a.mod_sqrt(Choice::from(1u8)).unwrap(); //Performs
//! // modular sqrt.
//! // Returs `None` if the input is not a QR on the field.
//! // Returns Some(result) if everything is correct.
//!
//! // Division has been also implemented. Remember that when we write
//! // a/b (mod l), we are indeed performing a * inverse_mod(b, l) (mod l).
//! assert!((-b / c) == EDWARDS_D);
//!
//! // Dividing by two even FieldElements is recommended through the `Half`
//! // trait implmementation since it's much faster.
//! if a.is_even() {
//!     let half_a = &a.half(); // This will panic if a isn't even.
//! };
//!
//! // We can finally perform inversion modulo l for a FieldElement:
//! let inv_a = &c.inverse(); // Performs a^-1 (mod l).
//!
//! // You can export your `FieldElement` as an slice of 32 bytes in Little
//! // Endian encoding by:
//! let c_bytes: [u8; 32] = c.to_bytes();
//! ```
//!
//! `PartialOrd`, `Ord`, `PartialEq` and `Eq` are also implemented for
//! `FieldElement` type.
//!
//! All `std::core::ops traits -> (Add, Sub, Mul, Div)` are implemented
//! for both, `&FieldElement` and `FieldElement`.

use core::cmp::PartialEq;
use std::ops::{MulAssign, Mul, SubAssign, Neg, Sub, AddAssign, Add};

use ff::{PrimeField, Field, PrimeFieldBits, FieldBits};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};


use rand::{CryptoRng, Rng};

use curve25519_dalek::scalar::Scalar;

use crate::backend;

#[cfg(feature = "u64_backend")]
pub use backend::u64::field::*;
/// A `FieldElement` represents an element of the field
/// `2^252 + 27742317777372353535851937790883648493`
///
/// The `FieldElement` type is an alias for one of the platform-specific
/// implementations.
#[cfg(feature = "u64_backend")]
pub type FieldElement = backend::u64::field::FieldElement;

impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl ConstantTimeEq for FieldElement {
    /// Test equality between two `FieldElement`s.  Since the
    /// internal representation is not canonical, the field elements
    /// are normalized to wire format before comparison.
    fn ct_eq(&self, other: &FieldElement) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

impl ConditionallySelectable for FieldElement {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        FieldElement([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
            u64::conditional_select(&a.0[4], &b.0[4], choice),
        ])
    }
}

impl From<&FieldElement> for Scalar {
    fn from(val: &FieldElement) -> Self {
        Scalar::from_bytes_mod_order(val.to_bytes())
    }
}

// impl Into<Scalar> for &FieldElement {
//     fn into(self) -> Scalar {
//         Scalar::from_bytes_mod_order(self.to_bytes())
//     }
// }

impl FieldElement {
    /// Generate a valid FieldElement choosen uniformly using user-
    /// provided rng.
    ///
    /// By `rng` we mean any Rng that implements: `Rng` + `CryptoRng`.
    pub fn random<T>(rand: &mut T) -> FieldElement
    where
        T: Rng + CryptoRng,
    {
        let mut bytes = [0u8; 32];
        rand.fill_bytes(&mut bytes);
        // Ensure that the value is lower than `FIELD_L`.
        bytes[31] &= 0b0000_0111;
        FieldElement::from_bytes(&bytes)
    }
}

impl PrimeField for FieldElement {
    type Repr = [u8;32];

    fn from_repr(repr: Self::Repr) -> subtle::CtOption<Self> {
        CtOption::new(
            Self::from_bytes(&repr),
            Choice::from(1u8),
        )
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        todo!()
    }

    const NUM_BITS: u32 = 256;

    const CAPACITY: u32 = 256;

    const S: u32 = 4;

    fn multiplicative_generator() -> Self {
        todo!()
    }

    fn root_of_unity() -> Self {
        todo!()
    }

    fn from_str_vartime(s: &str) -> Option<Self> {
        if s.is_empty() {
            return None;
        }

        if s == "0" {
            return Some(Self::zero());
        }

        let mut res = Self::zero();

        let ten = Self::from(10u32);

        let mut first_digit = true;

        for c in s.chars() {
            match c.to_digit(10) {
                Some(c) => {
                    if first_digit {
                        if c == 0 {
                            return None;
                        }

                        first_digit = false;
                    }

                    res.mul_assign(&ten);
                    res.add_assign(&Self::from(u64::from(c)));
                }
                None => {
                    return None;
                }
            }
        }

        Some(res)
    }

    fn from_repr_vartime(repr: Self::Repr) -> Option<Self> {
        Self::from_repr(repr).into()
    }

    fn is_even(&self) -> Choice {
        !self.is_odd()
    }
}

impl Field for FieldElement {
    fn random(rng: impl ff::derive::rand_core::RngCore) -> Self {
        todo!()
    }

    fn zero() -> Self {
        FieldElement::zero()
    }

    fn one() -> Self {
        FieldElement::one()
    }

    fn square(&self) -> Self {
        self * self
    }

    fn double(&self) -> Self {
        self + self
    }

    fn invert(&self) -> CtOption<Self> {
        todo!()
    }

    fn sqrt(&self) -> CtOption<Self> {
        todo!()
    }
}

impl<'b> Add<&'b FieldElement> for FieldElement {
  type Output = Self;

  #[inline]
  fn add(self, rhs: &'b FieldElement) -> Self::Output {
    self + *rhs
  }
}

impl<'a> Add<FieldElement> for &'a FieldElement {
  type Output = FieldElement;

  #[inline]
  fn add(self, rhs: FieldElement) -> Self::Output {
    *self + rhs
  }
}

impl AddAssign for FieldElement {
  #[inline]
  fn add_assign(&mut self, rhs: Self) {
    *self = *self + rhs;
  }
}

impl<'b> AddAssign<&'b FieldElement> for FieldElement {
  #[inline]
  fn add_assign(&mut self, rhs: &'b FieldElement) {
    *self = *self + rhs;
  }
}

impl<'b> Sub<&'b FieldElement> for FieldElement {
  type Output = Self;

  #[inline]
  fn sub(self, rhs: &'b FieldElement) -> Self::Output {
    self - *rhs
  }
}

impl<'a> Sub<FieldElement> for &'a FieldElement {
  type Output = FieldElement;

  #[inline]
  fn sub(self, rhs: FieldElement) -> Self::Output {
    *self - rhs
  }
}

impl SubAssign for FieldElement {
  #[inline]
  fn sub_assign(&mut self, rhs: Self) {
    *self = *self - rhs;
  }
}

impl<'b> SubAssign<&'b FieldElement> for FieldElement {
  #[inline]
  fn sub_assign(&mut self, rhs: &'b FieldElement) {
    *self = *self - rhs;
  }
}

impl<'b> Mul<&'b FieldElement> for FieldElement {
  type Output = Self;

  #[inline]
  fn mul(self, rhs: &'b FieldElement) -> Self::Output {
    self * *rhs
  }
}

impl<'a> Mul<FieldElement> for &'a FieldElement {
  type Output = FieldElement;

  #[inline]
  fn mul(self, rhs: FieldElement) -> Self::Output {
    *self * rhs
  }
}

impl MulAssign for FieldElement {
  #[inline]
  fn mul_assign(&mut self, rhs: Self) {
    *self = *self * rhs;
  }
}

impl<'b> MulAssign<&'b FieldElement> for FieldElement {
  #[inline]
  fn mul_assign(&mut self, rhs: &'b FieldElement) {
    *self = *self * rhs;
  }
}

impl PrimeFieldBits for FieldElement {
    type ReprBits = [u8; 32];

    fn to_le_bits(&self) -> FieldBits<Self::ReprBits> {
        unimplemented!();
    }

    fn char_le_bits() -> FieldBits<Self::ReprBits> {
        unimplemented!();
    }
}

#[cfg(test)]
pub mod test_field_element {
    use super::*;

    #[test]
    fn test_scalar_mul() {
        let f = FieldElement::from_str_vartime("340282366920938463463374607431768211455").unwrap();
        println!("f1*f2: {:?}", (f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f * f).to_bytes());
        let d = Scalar::from(u128::MAX);
        println!("d1*d2: {:?}", d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d * d);
    }
}
