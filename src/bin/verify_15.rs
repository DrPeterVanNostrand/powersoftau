extern crate pairing;
extern crate powersoftau;

use std::fs::OpenOptions;
use std::io::{Read, Write, BufWriter, BufReader};

use pairing::{
    bls12_381::{G1Affine, G2Affine},
    CurveAffine,
};
use powersoftau::{
    power_pairs, same_ratio, Accumulator, CheckForCorrectness, HashReader, UseCompression,
    ACCUMULATOR_BYTE_SIZE,
};

// Taken from phase1 participant attestations:
// https://github.com/arielgabizon/perpetualpowersoftau#ceremony-progress.
const RESPONSE_14_DIGEST: &str = "357738325a328e2c9fc5c49a8665baf9e5bd30d0475e2d842a7541e90930b9b8c8c74380f82f503255c2dffb35ef64232b3cd05681da013e06f1b6430a0f2589";
const RESPONSE_15_DIGEST: &str= "c47f56c5ab6aa5234fcbfbaa9ecfa5b450672e7d65a276bfc31f3b7e43ee65ca9f7c2ad15cc3a5c22040e11aaa749e1fecb0e6d7dd8f6976865d03ffa3a4bd85";

// Splits `points` (of length `n`) into two windows of length `n - 1`, then samples `n - 1` random
// scalars and returns each window's random linear combination (RLC) with the sampled scalars.
//
// The two RLCs will differ by a common factor `tau` whp. iff. each point in `points` was multiplied
// by a consecutive power of `tau`, i.e. `[tau]RLC_1 == RLC_2`.
#[inline]
fn rlc_ratio<G: CurveAffine>(points: &[G]) -> (G, G) {
    power_pairs(points)
}

fn main() {
    // Try to load `challenge_14` and `challenge_15` from disk.
    let challenge_14_reader = OpenOptions::new()
        .read(true)
        .open("challenge_14")
        .expect("unable to open file `challenge_14` in this directory");

    let challenge_15_reader = OpenOptions::new()
        .read(true)
        .open("challenge_15")
        .expect("unable to open file `challenge_15` in this directory");

    // Check challenge files' sizes.
    let challenge_14_size = challenge_14_reader
        .metadata()
        .expect("unable to get filesystem metadata for `challenge_14`")
        .len();
    let challenge_15_size = challenge_15_reader
        .metadata()
        .expect("unable to get filesystem metadata for `challenge_15`")
        .len();
    assert_eq!(challenge_14_size, ACCUMULATOR_BYTE_SIZE as u64);
    assert_eq!(challenge_15_size, ACCUMULATOR_BYTE_SIZE as u64);

    let challenge_14_reader = BufReader::new(challenge_14_reader);
    let mut challenge_14_reader = HashReader::new(challenge_14_reader);

    let challenge_15_reader = BufReader::new(challenge_15_reader);
    let mut challenge_15_reader = HashReader::new(challenge_15_reader);

    // Check participant 14 and 15's response digests.
    let response_14_digest: String = {
        let mut digest = [0u8; 64];
        challenge_14_reader
            .read_exact(&mut digest)
            .expect("unable to read `response_14` digest from `challenge_14`");
        digest.iter().map(|byte| format!("{:02x}", byte)).collect()
    };
    assert_eq!(response_14_digest, RESPONSE_14_DIGEST);
    let response_15_digest: String = {
        let mut digest = [0u8; 64];
        challenge_15_reader
            .read_exact(&mut digest)
            .expect("unable to read `response_15` digest from `challenge_15`");
        digest.iter().map(|byte| format!("{:02x}", byte)).collect()
    };
    assert_eq!(response_15_digest, RESPONSE_15_DIGEST);

    // Load into memory the accumulators output by participants 14 and 15.
    let before = Accumulator::deserialize(
        &mut challenge_14_reader,
        UseCompression::No,
        // This challenge file's points will have already been checked for valid subgroups and
        // for points at infinity.
        CheckForCorrectness::No,
    )
    .expect("unable to deserialize `challenge_14` accumulator");
    // Note that `challenge_15` does not contain participant 15's public key.
    let after = Accumulator::deserialize(
        &mut challenge_15_reader,
        UseCompression::No,
        CheckForCorrectness::Yes,
    )
    .expect("unable to deserialize `challenge_15` accumulator");

    // Get the digest of the accumulators output by participants 14 and 15.
    let acc_digest_before = challenge_14_reader.into_hash();
    let acc_digest_after = challenge_15_reader.into_hash();

    // Check that participant 15's first powers of tau are generators. Particpant 14's tau powers
    // have already passed this check.
    assert_eq!(after.tau_powers_g1[0], G1Affine::one());
    assert_eq!(after.tau_powers_g2[0], G2Affine::one());

    // Products of the participants' taus before and after the 15th participant in G1 and G2:
    // `[tau_1 * ... * tau_14]_1`
    let taus_14_g1 = before.tau_powers_g1[1];
    // `[tau_1 * ... * tau_15]_1`
    let taus_15_g1 = after.tau_powers_g1[1];
    // `[tau_1 * ... * tau_14]_2`
    let taus_14_g2 = before.tau_powers_g2[1];
    // `[tau_1 * ... * tau_15]_2`
    let taus_15_g2 = after.tau_powers_g2[1];

    // Check that the participant multiplied the powers of tau in G2 and G2 by the same value (i.e.
    // `tau_15`):
    // `[tau_1...tau_14]_1 / [tau_1...tau_15]_1 == [tau_1...tau_14]_2 / [tau_1...tau_15]_2`.
    //
    // For an honest participant this is equivalent to checking the ratio `1 / tau_15` in G1 and G2:
    // `[1]_1 / [tau_15]_1 == [1]_2 / [tau_15]_2`.
    let taus_ratio_g1 = (taus_14_g1, taus_15_g1);
    let taus_ratio_g2 = (taus_14_g2, taus_15_g2);
    assert!(same_ratio(taus_ratio_g1, taus_ratio_g2));

    // Check that the alpha powers of tau were updated by participant 15.
    // `[alpha_1 * ... * alpha_14]_1`
    let alphas_14_g1 = before.alpha_tau_powers_g1[0];
    // `[alpha_1 * ... * alpha_15]_1`
    let alphas_15_g1 = after.alpha_tau_powers_g1[0];
    assert_ne!(alphas_14_g1, alphas_15_g1);

    // Products of the participants' betas befefore and after the 15th participant in G1 and G2:
    // `[beta_1 * ... * beta_14]_1`
    let betas_14_g1 = before.beta_tau_powers_g1[0];
    // `[beta_1 * ... * beta_15]_1`
    let betas_15_g1 = after.beta_tau_powers_g1[0];
    // `[beta_1 * ... * beta_14]_2`
    let betas_14_g2 = before.beta_g2;
    // `[beta_1 * ... * beta_15]_2`
    let betas_15_g2 = after.beta_g2;

    // Check that the participant multiplied the beta powers of tau in G1 and the product of betas
    // in G2 by the same value (i.e. `beta_15`).
    // `[beta_1...beta_14]_1 / [beta_1...beta_15]_1 == [beta_1...beta_14]_2 / [beta_1...beta_15]_2`.
    //
    // For an honest participant this is equivalent to checking the ratio `1 / beta_15` in G1 and G2:
    // `[1]_1 / [beta_15]_1 == [1]_2 / [beta_15]_2`.
    //
    // The ratio `[beta_1 * ... * beta_14]_1 / [beta_1 * ... * beta_15]_1` in G1.
    let betas_ratio_g1 = (betas_14_g1, betas_15_g1);
    // The ratio `[beta_1 * ... * beta_14]_2 / [beta_1 * ... * beta_15]_2` in G2.
    let betas_ratio_g2 = (betas_14_g2, betas_15_g2);
    assert!(same_ratio(betas_ratio_g1, betas_ratio_g2));

    // Ratios in G1 and G2 of `1 / [tau_1 * ... * tau_15]`.
    let one_over_taus_g1 = (G1Affine::one(), taus_15_g1);
    let one_over_taus_g2 = (G2Affine::one(), taus_15_g2);

    // Check that the participant multiplied the powers of tau in G1 and G2 by consecutive powers
    // of their own tau (i.e. `tau_15`).
    let tau_powers_ratio_g1 = rlc_ratio(&after.tau_powers_g1);
    let tau_powers_ratio_g2 = rlc_ratio(&after.tau_powers_g2);
    assert!(same_ratio(tau_powers_ratio_g1, one_over_taus_g2));
    assert!(same_ratio(tau_powers_ratio_g2, one_over_taus_g1));

    // Check that the participant multiplied the alpha and beta powers of tau in G1 by consecutive
    // powers of their own tau (i.e. `tau_15`).
    let alpha_tau_powers_ratio_g1 = rlc_ratio(&after.alpha_tau_powers_g1);
    let beta_tau_powers_ratio_g1 = rlc_ratio(&after.beta_tau_powers_g1);
    assert!(same_ratio(alpha_tau_powers_ratio_g1, one_over_taus_g2));
    assert!(same_ratio(beta_tau_powers_ratio_g1, one_over_taus_g2));
}
