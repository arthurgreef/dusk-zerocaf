//! A `Scalar` represents an element of the finite field 
//! modulo `2^249 - 15145038707218910765482344729778085401`.
//! 
//! The `Scalar` type is an alias for one of the platform-specific
//! implementations. 
//! 
//! The most commonly used on all architectures will be the `u64` backend
//! which is the one that has been implemented at this moment.
use crate::backend;



/// A `Scalar` represents an element of the field GF(l), optimized for speed.
///
/// This is a type alias for one of the Scalar types in the `backend`
/// module.
#[cfg(feature = "u64_backend")]
pub type Scalar = backend::u64::scalar::Scalar;

/// This is a type alias for the Scalar type in the `curve25519-dalek` lib.
pub type Ristretto255Scalar = curve25519_dalek::scalar::Scalar;
