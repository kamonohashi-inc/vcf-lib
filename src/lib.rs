pub mod errors;
pub mod record;

#[derive(Debug, PartialEq, Eq)]
pub enum VariantType {
    SNV,
    Deletion,
    Insertion,
    Indel,
    MNV,
}
