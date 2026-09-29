# Remaining Python plotting scope

Plotting remains part of the first-port roadmap. It is not implemented or
qualified by the current scientific-output comparisons. This inventory separates
the original outputs and their inputs so they can be ported and checked without
claiming pixel-identical rendering.

| Stage | Original outputs | Inputs / behavior |
| --- | --- | --- |
| Preparation | `InSizeBias.pdf`, `InMAPBias.pdf`, `OutMAPBias.pdf`, `InGCBias.pdf`, `OutGCBias.pdf` | Original `EXCAVATORNormalizationExome.R`: before/after correction panels using bin medians and 10th/90th quantiles; size applies to IN only |
| Analysis, chromosome | `Plots/<sample>/PlotResults_<chromosome>.pdf` | Original `EXCAVATORPlotsExome.R`: two panels, IN/OUT log2 ratios and segment line above, ratio line and colored call intervals below; fixed log2 range -3 to 3 |
| Analysis, genome | `Plots/<sample>/PlotResults.pdf` | Same script: gtrellis ideograms and call intervals, assembly selected by name; depends on assembly annotations supplied by the R plotting stack |

Chromosome plots use the result tables directly and do not require a UCSC lookup.
The genome-wide layout is the part implicated in the offline synthetic-assembly
failure. A failure of that later plot must not be described as proof that no
chromosome PDFs were generated. The existing end-to-end gate does not inspect
their validity or appearance.

## Compatibility decisions and acceptance gates

1. Preserve output names, scientific data, ordering, axes and plotted interval
   coordinates; compare the arrays driving each plot independently of rendering.
2. Reproduce original normalization bin membership and correction checkpoints
   before drawing preparation diagnostics. Final normalized counts alone cannot
   reconstruct all intermediate panels.
3. Render small original/current chromosome examples with IN/OUT points, losses,
   gains and a no-call result; inspect PDF pages for missing data, clipping and
   incorrect class/interval associations. Pixel equality is not required.
4. Qualify genome-wide annotations separately. The target artifact currently
   carries centromeres but is not a complete ideogram source. Define and record
   local chromosome lengths/cytobands before claiming offline genome-wide parity;
   do not silently replace ideograms with another visualization.
5. Keep rendering failures visible and test the installed optional `plots` extra.
   Define publication behavior before integrating plotting into stage outputs.

## Original inconsistency to preserve or explicitly version

`LibraryFastCall.R:LabelAss` assigns negative labels to deletions and positive
labels to gains. In `EXCAVATORPlotsExome.R`, the genome-wide `sv` mapping follows
that convention, but chromosome rectangles use deletion colors for positive
labels and amplification colors for negative labels. Their vertical direction
does follow the sign. This is a source-confirmed original plotting inconsistency,
recorded in the local legacy backlog. It must not be silently repaired while
claiming strict visual compatibility. Any corrected presentation should be
documented as a deliberate difference; scientific tables/calls remain unchanged.

## M6 scope decision

The user explicitly deferred plotting on 2026-09-26. M6 closure covers scientific
arrays, decisions and tabular/VCF artifacts only. All preparation and analysis
plots above remain deferred work; M6 completion makes no visual-parity claim.
