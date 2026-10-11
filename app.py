from pathlib import Path
import io
import pandas as pd
import streamlit as st

from xrd_matcher import (
    decode_text, detect_experimental_peaks, extract_reference_metadata,
    infer_output_name, match_peaks, parse_reference_peak_table, read_experimental_xy
)

EXP_TYPES = ["xy", "txt", "dat", "csv", "xlsx", "xls"]
REF_TYPES = ["rtf", "txt", "dat", "xy", "csv", "xlsx", "xls", "pdf"]

def load_reference(uploaded_file):
    uploaded_file.seek(0)
    data = uploaded_file.read()
    suffix = Path(uploaded_file.name).suffix.lower()
    text = decode_text(data, suffix)
    meta = extract_reference_metadata(text, uploaded_file.name)
    peaks = parse_reference_peak_table(text)
    if peaks.empty and suffix in (".xy", ".txt", ".dat", ".csv", ".xlsx", ".xls"):
        try:
            frame = read_experimental_xy(data, uploaded_file.name)
            detected, _ = detect_experimental_peaks(frame, lower_limit=0, baseline_correction=False)
            peaks = detected[["peak_2theta", "peak_intensity"]].rename(
                columns={"peak_2theta": "two_theta", "peak_intensity": "reference_intensity"}
            )
        except Exception:
            pass
    return meta, peaks

st.set_page_config(page_title="XRD Phase-Matching Assistant", layout="wide")
st.title("XRD Phase-Matching Assistant")
st.write("Upload up to **5 experimental patterns** and **20 reference files**. This is a screening tool: positional matches do not prove a phase assignment.")

with st.sidebar:
    tolerance = st.number_input("Peak matching tolerance (° 2θ)", 0.05, 2.0, 0.60, 0.05)
    lower_limit = st.number_input("Start analysis at (° 2θ)", 0.0, 150.0, 20.0, 1.0)
    baseline_on = st.checkbox("Apply ALS baseline correction", value=True)
    prominence = st.number_input("Minimum peak prominence (0 = automatic)", min_value=0.0, value=0.0, step=10.0)
    separation = st.number_input("Minimum peak separation (° 2θ)", 0.01, 5.0, 0.15, 0.05)

experimental_files = st.file_uploader("1. Experimental XRD data (up to 5)", type=EXP_TYPES, accept_multiple_files=True)
reference_files = st.file_uploader("2. Standard reference files (up to 20)", type=REF_TYPES, accept_multiple_files=True)

if len(experimental_files) > 5:
    st.error("Upload no more than 5 experimental files.")
if len(reference_files) > 20:
    st.error("Upload no more than 20 reference files.")

if st.button("Run peak matching", type="primary",
              disabled=not experimental_files or not reference_files or len(experimental_files) > 5 or len(reference_files) > 20):
    try:
        exp_peaks, processed_data = {}, {}
        for uploaded in experimental_files:
            uploaded.seek(0)
            frame = read_experimental_xy(uploaded.read(), uploaded.name)
            peaks, processed = detect_experimental_peaks(
                frame, lower_limit=lower_limit, baseline_correction=baseline_on,
                prominence=prominence, min_separation_deg=separation
            )
            exp_peaks[uploaded.name] = peaks
            processed_data[uploaded.name] = processed

        refs, reference_data, warnings = {}, {}, []
        for uploaded in reference_files:
            meta, peaks = load_reference(uploaded)
            refs[uploaded.name] = (meta, peaks)
            reference_data[uploaded.name] = peaks
            if peaks.empty:
                warnings.append(f"{uploaded.name}: no reference peaks parsed. Export its peak list to TXT/RTF or provide a two-column pattern.")

        summary_rows, match_details = [], {}
        for reference_name, (meta, ref_peaks) in refs.items():
            row = {
                "Chemical formula / reference name": meta["chemical_formula_or_name"],
                "Crystal structure and space group": (
                    f'{meta["crystal_system"]}, {meta["space_group"]}'
                    if meta["space_group"] != "Not stated" else meta["crystal_system"]
                ),
                "Phase / classification": meta["phase_classification"],
                "PDF reference number": meta["pdf_reference_number"],
                "Reference file": reference_name,
                "Reference peaks parsed": len(ref_peaks),
            }
            for experimental_name, peaks in exp_peaks.items():
                count, detail = match_peaks(peaks, ref_peaks, tolerance=tolerance)
                row[f"Matched peaks — {Path(experimental_name).stem}"] = count
                match_details[f"{Path(reference_name).stem}_{Path(experimental_name).stem}"] = detail
            summary_rows.append(row)

        summary = pd.DataFrame(summary_rows)
        required = ["Chemical formula / reference name", "Crystal structure and space group",
                    "Phase / classification", "PDF reference number"]
        count_columns = [c for c in summary.columns if c.startswith("Matched peaks —")]
        other = [c for c in summary.columns if c not in required + count_columns]
        summary = summary[required + count_columns + other]
        output_stem = infer_output_name([f.name for f in experimental_files])

        output = io.BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            summary.to_excel(writer, sheet_name="Phase summary", index=False)
            for name, frame in match_details.items():
                frame.to_excel(writer, sheet_name=("Match_" + name)[:31], index=False)
            for name, frame in processed_data.items():
                frame.to_excel(writer, sheet_name=("Data_" + Path(name).stem)[:31], index=False)
            for name, frame in reference_data.items():
                frame.to_excel(writer, sheet_name=("Ref_" + Path(name).stem)[:31], index=False)

        st.success(f"Analysis completed: {output_stem}.xlsx")
        st.subheader("Phase-matching summary")
        st.dataframe(summary, use_container_width=True)
        st.download_button("Download Excel results", output.getvalue(), f"{output_stem}.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        st.download_button("Download summary CSV", summary.to_csv(index=False).encode("utf-8-sig"),
                           f"{output_stem}_summary.csv", "text/csv")
        if warnings:
            st.warning("\n\n".join(warnings))
        st.info("Check reference quality, multiple diagnostic reflections, chemistry/EDS and peak overlap before confirming phase assignments.")
    except Exception as error:
        st.error(f"Analysis failed: {error}")
        st.exception(error)
