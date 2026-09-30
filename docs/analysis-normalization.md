# Original automatic normalization

The original EXCAVATOR2 normalization and ratio method is restored. There is no
calibration file, exposure input, half-read pseudocount or raw-count-only analysis
requirement.

## Preparation

1. Divide counts by window length (`end - start`).
2. Apply the original size correction to IN windows.
3. Apply mappability correction separately to IN and OUT.
4. Apply GC correction separately to IN and OUT.
5. Replace zero corrected values with the smallest nonzero corrected value in
   the same sample and IN/OUT class.

Corrections use the original five-unit covariate bins and median rescaling.
The original treatment of partial size bins and nonpositive bin medians is
retained. All-zero classes fail explicitly. Preparation saves the weighted,
size, MAP, GC and final normalized checkpoints alongside observed counts.

## Analysis

Analysis reads the normalized column of the prepared seven-column matrices.
Original converted prepared artifacts do not need raw-count arrays or a new
normalization declaration.

Paired analysis uses the matching control. Pooling averages the normalized
control values at each window, retaining the original 15-significant-digit text
roundtrip. Analysis computes `log2(test / reference)`, then subtracts the IN median
from IN windows and the OUT median from OUT windows. No pseudocount is added.
Nonfinite ratios fail explicitly.

```sh
excavator2 analyze --samples design.yaml --input prepared/ --target target/ \
  --output results/ --experiment pooling
```

This restores the original normalization assumptions and limitations, including
possible attenuation of covariate-associated or broad alterations, replacement
of zero observations, and control-depth influence on pooling. P13, P14, A17 and
P15 are not claimed resolved through a new normalization model.

The separately implemented MAPQ, counting, target geometry, feature alignment,
segmentation, FastCall and reporting fixes remain in place. Consequently, whole
pipeline output is not expected to equal the original for every input.

## Verification

Tests compare all normalization checkpoints and paired/pooling ratios against
saved original-R fixtures. CLI tests use original normalized prepared matrices
without calibration files or raw-count arrays. Additional cases cover class-wise
zero replacement and arithmetic pooling before separate IN/OUT centering.
