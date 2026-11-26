# KASM comparison with the supplied plots

Checked 2026-11-02 against all 72 images in the supplied shares:

- [Disturbed interval, 8–15 May 2024](https://cloud.irf.se/s/WWcmCyBQwoGFXHr)
- [Quiet interval, 21–30 December 2023](https://cloud.irf.se/s/xqtr7ZiJdsnZxxz)

The images cover BEL, HYB, KIR and WIC, with eight three-hour UTC indices per
image. The `KAsm` rows were visually transcribed using image contact sheets;
the final reference values do not rely on OCR. Calculations use the repository's
IAF minute data, both adjacent-day wings, smoothing 1, and K9 limits BEL=450,
HYB=300, KIR=1500 and WIC=500 nT.

## K indices

| Station | Matching intervals, default coverage | Total |
| --- | ---: | ---: |
| BEL | 143 | 144 |
| HYB | 144 | 144 |
| KIR | 144 | 144 |
| WIC | 143 | 144 |
| Total | 574 | 576 |

Two intervals differ:

- **BEL, 2024-05-13, 18–21 UTC:** image K=3; implementation K=2.
  The implementation's filtered component indices are X=3 and Y=2. Its final
  combined index is capped at the unfiltered combined index, as required by
  the author's `wwsk` routine. The raw normalized ranges are X=19.8888889 and
  Y=9.5555556; K=3 starts at 20. The final result is therefore K=2.
  A fresh build of the complete author's KASM 1.11 distribution, run against
  the same local IAF input, independently returns `4 5 3 4 4 4 2 2` for this
  day, agreeing with our implementation. The image shows `4 5 3 4 4 4 3 2`.
  Rounding the input to whole nT reproduces that day's image K sequence, but
  this does not establish how the image was generated. Its exact source data,
  input conversion and program version would be needed to resolve the difference.
- **WIC, 2023-12-22, 21–24 UTC:** image K=0; implementation gives no valid K
  under the default 90% observed-coverage cutoff. Each component has 159/180
  valid observed minutes (88.33%). With `minimum_coverage=.01`, KASM's own gap
  checks accept the interval and return K=0, matching the plot.

With the additional application coverage restriction relaxed, **575/576**
indices match. There is no solver change justified by these results.

## Component ranges and curves

The plots state that signals are normalized by 500/K9. The comparison CSV
therefore includes computed component ranges in those same normalized units.
The decimal ranges are not all identical: for BEL on 13 May at 18–21 UTC,
the image prints X=22.3 and Y=11.0, whereas the local calculation gives
X=22.2236097 and Y=11.1446731 (22.2 and 11.1 at one decimal place).
Thus the plots cannot be claimed to be exact numerical reproductions of this
implementation with the local inputs. The full plotted curves were not digitized;
this check establishes the K-index agreement and documents an example of range
disagreement, rather than asserting pointwise curve equality.

## Reproduction and checks

Run `python validation/check_kasm_images.py` to repeat all 576 K comparisons.
The reference JSON records input and image SHA-256 hashes and transcribed K
values. The CSV contains every interval's reference K, computed K, relaxed
coverage K, observed coverage and computed normalized ranges. Downloads and OCR
are not needed to reproduce the K comparison.

Existing tests passed:

- `python -m unittest discover -s tests -p 'test_kasm_k_index.py'` — 8 tests,
  including 1,240 independent standalone-author reference indices.
- `python -m unittest discover -s tests -p 'test_python_ports.py'` — 2 tests,
  including native diagnostic component ranges and gap flags.

The additional author-program check used the distribution at
https://magneto.igf.edu.pl/soft/kasm/KAsm.zip, built with its supplied Makefile,
and the command below (run with copies of `bel24apr.bin` and `bel24may.bin`
in `INPUT_DIRECTORY/`, without modifying the originals):

```
KAsm BEL:13MAY2024:1 450:1 xy OUT INPUT_DIRECTORY/ -IAF
```

No production calculation code was changed.
