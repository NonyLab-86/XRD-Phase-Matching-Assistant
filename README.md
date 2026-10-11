# XRD Phase-Matching Assistant

A Streamlit application for screening up to **5 experimental XRD patterns** against up to **20 reference files**.

## What it does

- Reads two-column experimental patterns (`2θ`, intensity).
- Optionally applies asymmetric least-squares (ALS) baseline correction.
- Detects local maxima using configurable prominence and peak-separation settings.
- Matches each experimental peak to its nearest reference reflection within a configurable `2θ` tolerance.
- Produces an Excel workbook with the requested summary columns:
  1. Chemical formula / reference name
  2. Crystal structure and space group
  3. Phase / classification
  4. PDF reference number
  5. Matched peak count for each uploaded experimental file
- Includes per-reference match-detail sheets, processed experimental data and parsed reference peaks.
- Infers the output filename from the uploaded experimental filenames.

## Project structure

```text
xrd-phase-matching/
├── app.py
├── xrd_matcher.py
├── requirements.txt
├── requirements-dev.txt
├── README.md
├── tests/
│   ├── test_peak_matching.py
│   ├── test_reference_parser.py
│   ├── test_output_summary.py
│   └── data/
│       ├── synthetic_sample.xy
│       └── synthetic_reference.txt
└── .github/
    └── workflows/
        └── tests.yml
```

## Run locally

Python 3.10+ recommended.

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
streamlit run app.py
```

## Run automated tests

Install development dependencies and run the tests:

```bash
pip install -r requirements-dev.txt
pytest -v
```

The tests cover known-answer peak matching, tolerance boundaries, empty reference handling, reference metadata parsing, reference peak-list parsing, and inferred output names.

GitHub Actions runs the test suite when code is pushed or a pull request is opened.

## Supported files

**Experimental:** `.xy`, `.txt`, `.dat`, `.csv`, `.xlsx`, `.xls`  
**Reference:** `.rtf`, `.txt`, `.dat`, `.xy`, `.csv`, `.xlsx`, `.xls`, `.pdf`

The most reliable reference input is a text/RTF export with metadata and a peak list. Two-column reference patterns are also supported through local-maximum detection.

**PDF limitation:** PDF uploads are accepted by the interface, but this version does not include robust PDF table extraction or OCR. If a PDF does not parse, export its reference peak table to TXT or RTF first. Scanned PDFs are not parsed automatically.

## Matching definition

For each detected experimental peak, the program finds the nearest reference peak. It counts the experimental peak as matched when the absolute difference in `2θ` is less than or equal to the chosen tolerance. Each experimental peak counts at most once per reference, but can match more than one reference if different references contain nearby reflections.

## Scientific limitations

- Peak-position matching is screening evidence, not definitive phase identification.
- Phase labels are heuristic if the reference metadata does not explicitly state the phase.
- A reference with no parsed peaks is not meaningfully tested; inspect the warnings and `Reference peaks parsed` column.
- Baseline correction can affect peak detection. Assess broad low-angle features using blank/sample-holder scans and instrument context.
- Overlapping peaks can be consistent with multiple candidate phases.
- For publication-grade assignments, inspect reference quality, multiple diagnostic reflections, chemistry/EDS, and consider profile fitting or Rietveld refinement.
- Reference database content may be subject to licence restrictions. Ensure you have permission to redistribute reference files.

## Reproducibility

The files under `tests/data/` are synthetic and created only to test software behaviour. They are not real diffraction standards and must not be used to identify physical phases.
