//! Module for Public Trait implementations.


/// Gives the Identity element for the
/// type which it has been implemented on.
/// 
/// This trait is implemented following the rules that
/// mandate over the Type that is being implemented. 
pub trait Identity {

    #[must_use]
    /// Returns the identity element for the implemented
    /// type, following it's mathematical rules.
    fn identity() -> Self;
}

pub mod ops {
    use subtle::Choice;

    /// Trait that represents the `^2` operation for any
    /// kind of element on the library. 
    /// 
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented. 
    pub trait Square {
        type Output;

        #[must_use]
        /// Returns the square of the input: `x^2`.
        fn square(self) -> Self::Output;
    }

    /// Trait that represents the Point doubling operation
    /// for any type of Point that is used on the lib. 
    /// 
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Double {
        type Output;

        #[must_use]
        /// Performs the point-doubling operation over the
        /// coordinates which this trait has been implemented
        /// for.
        fn double(self) -> Self::Output;
    }

    /// Trait that represents the `/2` operation for any
    /// kind of element on the library.
    /// 
    /// This is a more performant way of performing the
    /// division by 2 that dividing by 2 with the `Div`
    /// trait implementation.
    /// 
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Half {
        type Output;

        #[must_use]
        /// Returns the half of the input: `x/2`.
        fn half(self) -> Self::Output;
    }

    /// Trait that represents the modular exponentiation 
    /// operation, ie.`a ^ b (mod l)`, for any
    /// kind of element on the library (except points).
    ///
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Pow<T> {
        type Output;

        #[must_use]
        /// Returns  `a^b (mod l)`. 
        fn pow(self, exp: T) -> Self::Output;
    }

    pub trait ModSqrt {
        type Output;

        #[must_use]
        /// Performs the modular Square Root operation over a finite
        /// field ie. `sqrt(x) (mod l)`. 
        /// 
        /// With the given `Choice`, the impl is able to provide the 
        /// result that corresponds to the positive or negative sign choosen.
        /// 
        /// # Returns 
        /// 
        /// `Some(symb_choosen_result)` if the input is a QR for the prime modulo.
        /// Otherways it returns `None`
        fn mod_sqrt(self, choice: Choice) -> Self::Output;
    }

    pub trait InvSqrt {
        type Output;

        #[must_use]
        /// Performs the Inverse Square root as it's defined on the
        /// Ristretto paper.
        /// 
        /// This operation returns always the positive result of the
        /// modular sqrt, understanding positive as the definition that
        /// appears on the Decaf paper: 0 < result < (P - 1)/2.  
        fn inv_sqrt(self) -> Self::Output;
    }
}