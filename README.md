# Deciphering Beale Cipher 1

A statistical and geographical exploration of Beale Cipher 1 by **William Duckworth**.

**[Read the paper (PDF)](deciphering-beale-cipher-1.pdf)** · 7 pages · February 2025

---

## Status

**Exploratory hobbyist analysis — not a peer-reviewed result.**

This paper is presented as a thought experiment, not as a solved cipher. The proposed coordinate recovery depends on:

- A chosen two-digit segmentation of the cipher's number sequence
- A choice of which halves of the cipher are read as latitude vs. longitude

Both of these introduce researcher degrees of freedom, which means the result is vulnerable to multiple-comparison / p-hacking concerns that the paper itself acknowledges but does not fully defeat. The Beale Papers' historical authenticity also remains disputed in the cryptography community. Read accordingly.

## Summary

The paper proposes that two halves of Beale Cipher 1 encode digit sequences that, when interpreted as Degrees–Minutes–Seconds (DMS), yield geographic coordinates pointing into Bedford County, Virginia:

- **Latitude digits:** 371221
- **Longitude digits:** 792316

These two segments do **not** appear as one continuous 12-digit block; they are recovered from separate halves of the cipher and combined.

## What the paper actually does

- Reviews historical context and authenticity debates around the Beale Papers (1885).
- Justifies a two-digit segmentation as historically plausible for 19th-century surveying notation.
- Estimates, using a simplified statistical model (with brief treatment of digit correlations), the probability of coordinate-shaped digit strings emerging by chance.
- Cross-references the already-solved Cipher 2, which states: *"Paper number '1' describes the exact locality of the vault so that no difficulty will be had in finding it."*
- Discusses counterarguments — coincidence, data-dredging, single-coordinate scope, and authenticity — without claiming to have refuted them.

## Citing

> Duckworth, W. *Deciphering Beale Cipher 1: A Comprehensive Statistical and Geographical Analysis.* 2025.

## License

Released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
