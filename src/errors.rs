use thiserror::Error;

pub type Result<T, E = Error> = core::result::Result<T, E>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Reference bases must not be empty")]
    RefBasesEmptyError(),

    #[error("Alternate bases must not be empty")]
    AltBasesEmptyError(),

    #[error("Reference bases contains non-ACGT characters: {0}")]
    RefBasesInvalidSymbolError(String),

    #[error("Alternate bases contains non-ACGT characters: {0}")]
    AltBasesInvalidSymbolError(String),
}
