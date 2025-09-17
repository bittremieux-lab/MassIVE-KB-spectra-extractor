# !/usr/bin/env python3
import sys
from ftplib import FTP, error_perm
from pathlib import Path

import pandas as pd
from pyteomics import mzml, mzxml, mgf


def _download(mskb_version, mzml_file, ftp_host, out_f):
    remote_path = f"{mskb_version}/{mzml_file}"
    ftp_host.retrbinary(f"RETR {remote_path}", out_f.write)


def download_ftp(host="massive-ftp.ucsd.edu", mzml_file=None, local_path=None):
    with FTP(host) as ftp:
        ftp.login()  # Anonymous login
        with open(local_path, "wb") as f:
            mskb_version = "z01" if "ccms_peak" in mzml_file else "v01"
            try:
                _download(mskb_version, mzml_file, ftp, f)
            except error_perm as e:
                _download("x01", mzml_file, ftp, f)


def mzml_spectrum_to_mgf(spectrum, mzml_file_name, modified_peptide, scan):
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
            precursor_charge = None
        retention_time = spectrum["scanList"]["scan"][0].get("scan start time", -1)

    elif mzml_file_name.endswith(".mzXML"):
        retention_time = spectrum.get("retentionTime", -1)

        precursor_mz = spectrum["precursorMz"][0]["precursorMz"]
        if "precursorCharge" in spectrum["precursorMz"][0]:
            precursor_charge = spectrum["precursorMz"][0]["precursorCharge"]
        else:
            precursor_charge = None
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
    local_file = Path(mzml_file).name
    spectra_count = len(df)  # Get number of spectra from dataframe length

    # Create persistent failed logs directory in the main pipeline directory
    # Get pipeline directory from environment variable set by Nextflow
    import os

    pipeline_dir = Path(os.environ.get("PIPELINE_DIR", "."))
    failed_logs_dir = pipeline_dir / "failed_logs"
    failed_logs_dir.mkdir(exist_ok=True)

    # Create unique failed file name based on input file (now .csv extension)
    input_basename = Path(tsv_file_path).name
    failed_file_path = failed_logs_dir / f"{input_basename}.csv"

    try:
        download_ftp(mzml_file=mzml_file, local_path=local_file)
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
                for s, pep, scan in zip(spectra, df["annotation"], df["scan"])
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
    tsv_file_path = sys.argv[1]
    success = process_mzml_group(tsv_file_path)
    if not success:
        sys.exit(1)  # Exit with error code 1 for failures
