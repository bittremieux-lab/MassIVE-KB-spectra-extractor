#!/usr/bin/env python3
import sys
from ftplib import FTP
from pathlib import Path

import pandas as pd
from pyteomics import mzml, mzxml, mgf


def download_ftp(host="massive-ftp.ucsd.edu", remote_path=None, local_path=None):
    with FTP(host) as ftp:
        ftp.login()  # Anonymous login
        with open(local_path, "wb") as f:
            ftp.retrbinary(f"RETR {remote_path}", f.write)


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


def process_mzml_group(tsv_file_path):
    df = pd.read_csv(tsv_file_path, sep="\t")
    mzml_file = df.loc[0, "filename"]
    local_file = Path(mzml_file).name

    if "ccms_peak" in mzml_file:
        massivekb_version = "z01"
    else:
        massivekb_version = "v01"
    download_ftp(remote_path=f"{massivekb_version}/{mzml_file}", local_path=local_file)

    if mzml_file.endswith(".mzML"):
        parser = mzml.MzML
        id_fmt = "controllerType=0 controllerNumber=1 scan=%i"
    elif mzml_file.endswith(".mzXML"):
        parser = mzxml.MzXML
        id_fmt = "%i"
    else:
        raise ValueError(f"Unsupported file type: {mzml_file}")

    mzml_reader = parser(local_file)
    formatted_scan_numbers = df["scan"].map(lambda x: id_fmt % x)
    try:
        spectra = mzml_reader.get_by_ids(formatted_scan_numbers)
    except KeyError:
        print(
            f"Trying to get {formatted_scan_numbers} from {mzml_file} resulted in a KeyError"
        )
        raise

    mgf.write(
        [
            mzml_spectrum_to_mgf(s, local_file, pep, scan)
            for s, pep, scan in zip(spectra, df["annotation"], df["scan"])
        ],
        f"{local_file}.mgf",
        fragment_format="%.5f %.1f",
        use_numpy=True,
    )


tsv_file_path = sys.argv[1]
process_mzml_group(tsv_file_path)
