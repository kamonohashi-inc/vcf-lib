use crate::errors::{Error, Result};
use crate::VariantType;
use once_cell::sync::Lazy;
use regex::Regex;

static REGEX_ALLELES: Lazy<Regex> = Lazy::new(|| Regex::new(r"\A[ACGTURYKMSWBDHVN]+\z").unwrap());

pub fn normalize<'a>(
    position: u64,
    reference: &'a str,
    alternate: &'a str,
) -> Result<(u64, &'a str, &'a str)> {
    if reference.len() == 0 {
        Err(Error::RefBasesEmptyError())?
    }

    if alternate.len() == 0 {
        Err(Error::AltBasesEmptyError())?
    }

    if !REGEX_ALLELES.is_match(reference) {
        Err(Error::RefBasesInvalidSymbolError(reference.to_string()))?
    }

    if !REGEX_ALLELES.is_match(alternate) {
        Err(Error::AltBasesInvalidSymbolError(alternate.to_string()))?
    }

    let (r, a) = trim_trailing_shared_bases(reference, alternate);

    Ok(trim_leading_shared_bases(position, r, a))
}

pub fn variant_type(reference: &str, alternate: &str) -> Option<VariantType> {
    match (reference, alternate) {
        (r, a) if r.len() == 1 && a.len() == 1 && a != r => Some(VariantType::SNV),
        (r, a) if r.len() == a.len() && a != r => Some(VariantType::MNV),
        (r, a) if r.len() == 1 && a.len() > 1 && r.get(0..1) == a.get(0..1) => {
            Some(VariantType::Insertion)
        }
        (r, a) if a.len() == 1 && r.len() > 1 && r.get(0..1) == a.get(0..1) => {
            Some(VariantType::Deletion)
        }
        (r, a) if r.len() != a.len() => Some(VariantType::Indel),
        _ => None,
    }
}

fn trim_trailing_shared_bases<'b>(reference: &'b str, alternate: &'b str) -> (&'b str, &'b str) {
    let mut itr_r = reference.chars().rev();
    let mut itr_a = alternate.chars().rev();
    let i = count_shared(&mut itr_r, &mut itr_a);

    let mut p1 = reference.len() - i;
    let mut p2 = alternate.len() - i;

    if p1 == 0 || p2 == 0 {
        p1 += 1;
        p2 += 1;
    }

    (&reference[0..p1], &alternate[0..p2])
}

fn trim_leading_shared_bases<'b>(
    position: u64,
    reference: &'b str,
    alternate: &'b str,
) -> (u64, &'b str, &'b str) {
    let mut itr_r = reference.chars();
    let mut itr_a = alternate.chars();
    let mut i = count_shared(&mut itr_r, &mut itr_a);

    if i == reference.len() || i == alternate.len() {
        i -= 1;
    }

    (position + i as u64, &reference[i..], &alternate[i..])
}

fn count_shared<T: Iterator<Item = char>>(itr_r: &mut T, itr_a: &mut T) -> usize {
    let mut i = 0;

    while let (Some(c1), Some(c2)) = (itr_r.next(), itr_a.next()) {
        if c1 == c2 {
            i += 1
        } else {
            break;
        }
    }

    i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_1() {
        let (p, r, a) = normalize(1000, "A", "T").unwrap();

        assert_eq!(p, 1000);
        assert_eq!(r, "A");
        assert_eq!(a, "T");
        assert_eq!(variant_type(r, a), Some(VariantType::SNV));
    }

    #[test]
    fn test_normalize_2() {
        let (p, r, a) = normalize(1000, "A", "AT").unwrap();

        assert_eq!(p, 1000);
        assert_eq!(r, "A");
        assert_eq!(a, "AT");
        assert_eq!(variant_type(r, a), Some(VariantType::Insertion));
    }

    #[test]
    fn test_normalize_3() {
        let (p, r, a) = normalize(1000, "AT", "A").unwrap();

        assert_eq!(p, 1000);
        assert_eq!(r, "AT");
        assert_eq!(a, "A");
        assert_eq!(variant_type(r, a), Some(VariantType::Deletion));
    }

    #[test]
    fn test_normalize_4() {
        let (p, r, a) = normalize(1000, "AT", "ATA").unwrap();

        assert_eq!(p, 1001);
        assert_eq!(r, "T");
        assert_eq!(a, "TA");
        assert_eq!(variant_type(r, a), Some(VariantType::Insertion));
    }

    #[test]
    fn test_normalize_5() {
        let (p, r, a) = normalize(1000, "ATCC", "ATACC").unwrap();

        assert_eq!(p, 1001);
        assert_eq!(r, "T");
        assert_eq!(a, "TA");
        assert_eq!(variant_type(r, a), Some(VariantType::Insertion));
    }

    #[test]
    fn test_normalize_6() {
        let (p, r, a) = normalize(1000, "ACTCC", "AGTTCC").unwrap();

        assert_eq!(p, 1001);
        assert_eq!(r, "C");
        assert_eq!(a, "GT");
        assert_eq!(variant_type(r, a), Some(VariantType::Indel));
    }

    #[test]
    fn test_normalize_7() {
        let (p, r, a) = normalize(1000, "A", "A").unwrap();

        assert_eq!(p, 1000);
        assert_eq!(r, "A");
        assert_eq!(a, "A");
        assert_eq!(variant_type(r, a), None);
    }

    #[test]
    fn test_normalize_err_1() {
        assert!(normalize(1000, "", "A").is_err());
    }

    #[test]
    fn test_normalize_err_2() {
        assert!(normalize(1000, "A", "").is_err());
    }

    #[test]
    fn test_normalize_err_3() {
        assert!(normalize(1000, ".", "A").is_err());
    }

    #[test]
    fn test_normalize_err_4() {
        assert!(normalize(1000, "A", ".").is_err());
    }
}
