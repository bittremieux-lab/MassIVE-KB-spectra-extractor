# !/usr/bin/env python3
import os
import sys
from ftplib import FTP, error_perm
from pathlib import Path

import pandas as pd
from pyteomics import mzml, mzxml, mgf

FTP_REPLACEMENTS = {"MSV000083508/ccms_peak_centroided": "MSV000083508/ccms_peak"}


def fix_MSV000080620(mzml_file):
    if "MSV000080620/ccms_peak/RAW/" in mzml_file:
        f = mzml_file.split("/")[3]
        day, a, b, c, d, e, ext = f.split("_")
        d0 = d[:3]
        f = f.replace(".mzML", ".mzXML")
        mzml_file = f"MSV000080620/ccms_peak/{d0}/{d}/{f}"
    return mzml_file


def _download(mskb_version, mzml_file, ftp_host):
    mzml_file = fix_MSV000080620(mzml_file)
    for k, v in FTP_REPLACEMENTS.items():
        if k in mzml_file:
            mzml_file = mzml_file.replace(k, v)
            break
    remote_path = f"{mskb_version}/{mzml_file}"
    local_file = Path(mzml_file).name
    with open(local_file, "wb") as out_f:
        ftp_host.retrbinary(f"RETR {remote_path}", out_f.write)
    return mzml_file, local_file


def download_ftp(host="massive-ftp.ucsd.edu", mzml_file=None):
    with FTP(host) as ftp:
        ftp.login()  # Anonymous login
        mskb_version = "z01" if "ccms_peak" in mzml_file else "v01"
        try:
            return _download(mskb_version, mzml_file, ftp)
        except error_perm as e:
            return _download("x01", mzml_file, ftp)


def mzml_spectrum_to_mgf(spectrum, mzml_file_name, modified_peptide, scan, charge):
    mgf_spectrum = {
        "m/z array": spectrum["m/z array"],
        "intensity array": spectrum["intensity array"],
    }

    if mzml_file_name.endswith(".mzML"):
        precursor = spectrum["precursorList"]["precursor"][0]
        precursor_ion = precursor["selectedIonList"]["selectedIon"][0]
        precursor_mz = precursor_ion["selected ion m/z"]
        if "charge state" in precursor_ion:
            precursor_charge = int(precursor_ion["charge state"])
        elif "possible charge state" in precursor_ion:
            precursor_charge = int(precursor_ion["possible charge state"])
        else:
            precursor_charge = charge
        retention_time = spectrum["scanList"]["scan"][0].get("scan start time", -1)

    elif mzml_file_name.endswith(".mzXML"):
        retention_time = spectrum.get("retentionTime", -1)

        precursor_mz = spectrum["precursorMz"][0]["precursorMz"]
        if "precursorCharge" in spectrum["precursorMz"][0]:
            precursor_charge = spectrum["precursorMz"][0]["precursorCharge"]
        else:
            precursor_charge = charge
    else:
        raise ValueError(f"Unsupported file type: {mzml_file_name}")

    mgf_spectrum["params"] = {
        "TITLE": f"{mzml_file_name}:scan:{scan}",
        "PEPMASS": precursor_mz,
        "RTINSECONDS": retention_time,
        "CHARGE": precursor_charge,
        "SEQ": modified_peptide,
    }
    return mgf_spectrum


def write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count):
    """Write failure information as a single CSV row"""
    # Escape quotes in error message for CSV format
    error_msg_escaped = error_msg.replace('"', '""')
    csv_row = f'"{mzml_file}","{error_msg_escaped}",{spectra_count}\n'
    failed_file_path.write_text(csv_row)


def process_mzml_group(tsv_file_path):
    df = pd.read_csv(tsv_file_path, sep="\t")
    mzml_file = df.loc[0, "filename"]
    spectra_count = len(df)  # Get number of spectra from dataframe length

    # Create persistent failed logs directory in the main pipeline directory
    # Get pipeline directory from environment variable set by Nextflow
    pipeline_dir = Path(os.environ.get("PIPELINE_DIR", "."))
    task_id = os.environ.get("TASK_ID", "unknown")
    failed_logs_dir = pipeline_dir / f"failed_logs_{task_id}"
    failed_logs_dir.mkdir(exist_ok=True)

    # Create unique failed file name based on input file (now .csv extension)
    input_basename = Path(tsv_file_path).name
    failed_file_path = failed_logs_dir / f"{input_basename}.csv"

    try:
        mzml_file, local_file = download_ftp(mzml_file=mzml_file)
    except error_perm as e:
        mskb_version = "z01" if "ccms_peak" in mzml_file else "v01"
        error_msg = f"Tried getting {mzml_file} from {mskb_version} and x01 failed. Error: {str(e)}"
        write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count)
        print(f"Wrote {failed_file_path}")
        return False  # Return False instead of raising exception

    if mzml_file.endswith(".mzML"):
        parser = mzml.MzML
        id_fmt_l = ["controllerType=0 controllerNumber=1 scan=%i", "scan=%i"]
    elif mzml_file.endswith(".mzXML"):
        parser = mzxml.MzXML
        id_fmt_l = ["%i"]
    else:
        error_msg = f"Unsupported file type: {mzml_file}"
        write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count)
        print(f"Wrote {failed_file_path}")
        return False

    try:
        mzml_reader = parser(local_file)
        for i, id_fmt in enumerate(id_fmt_l):
            formatted_scan_numbers = df["scan"].map(lambda x: id_fmt % x)
            try:
                spectra = mzml_reader.get_by_ids(formatted_scan_numbers)
                break  # on success, stop trying and don't execute the else statement
            except KeyError:
                continue
        else:
            # This runs only if the loop never break'ed
            error_msg = f"Tried to get scans with {id_fmt_l} from {mzml_file} resulted in a KeyError"
            write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count)
            print(f"Wrote {failed_file_path}")
            return False

        mgf.write(
            [
                mzml_spectrum_to_mgf(s, local_file, pep, scan)
                for s, pep, scan, charge in zip(
                    spectra, df["annotation"], df["scan"], df["charge"]
                )
            ],
            f"{tsv_file_path}.mgf",
            fragment_format="%.5f %.1f",
            use_numpy=True,
        )
    except Exception as e:
        error_msg = f"Error processing {mzml_file}: {str(e)}"
        write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count)
        print(f"Wrote {failed_file_path}")
        return False
    finally:
        # Always clean up the downloaded file
        if Path(local_file).exists():
            Path(local_file).unlink()

    return True  # Success


if __name__ == "__main__":
    try:
        tsv_file_path = sys.argv[1]
        success = process_mzml_group(tsv_file_path)
        if not success:
            sys.exit(1)  # Exit with error code 1 for failures
    except Exception as e:
        # General exception handler to catch any unexpected errors
        try:
            # Try to get basic information for error logging
            tsv_file_path = sys.argv[1] if len(sys.argv) > 1 else "unknown_file"

            # Try to read the TSV file to get mzml_file and spectra_count
            try:
                df = pd.read_csv(tsv_file_path, sep="\t")
                mzml_file = df.loc[0, "filename"]
                spectra_count = len(df)
            except:
                # If we can't read the file, use fallback values
                mzml_file = "unknown_mzml_file"
                spectra_count = 0

            # Set up failed logs directory
            pipeline_dir = Path(os.environ.get("PIPELINE_DIR", "."))
            task_id = os.environ.get("TASK_ID", "unknown")
            failed_logs_dir = pipeline_dir / f"failed_logs_{task_id}"
            failed_logs_dir.mkdir(exist_ok=True)

            # Create failed file path
            input_basename = Path(tsv_file_path).name
            failed_file_path = failed_logs_dir / f"{input_basename}.csv"

            # Write the exception details to failed logs
            error_msg = f"Unexpected exception in mzml_group_to_mgf: {type(e).__name__}: {str(e)}"
            write_failure_csv(failed_file_path, mzml_file, error_msg, spectra_count)
            print(f"Unexpected exception caught and logged to {failed_file_path}")
            print(f"Exception: {type(e).__name__}: {str(e)}")

        except Exception as logging_error:
            # If even the error logging fails, print to stderr
            print(
                f"Critical error: Failed to log exception. Original error: {type(e).__name__}: {str(e)}",
                file=sys.stderr,
            )
            print(
                f"Logging error: {type(logging_error).__name__}: {str(logging_error)}",
                file=sys.stderr,
            )

        # Always exit with code 1 for any exception
        sys.exit(1)
