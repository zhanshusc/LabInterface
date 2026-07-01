"""
How to use:
1. Edit the section below.
2. Run this script
3. In LifeFit UI, choose:
       Fileformat: time_intensity
   Then upload the generated decay .txt and IRF .txt files.

Output format: .txt in two columns 
    time_ns    intensity
"""

from pathlib import Path

import numpy as np
from sdtfile import SdtFile


# Edit these for change 

# Put your SDT file paths here.
DECAY_SDT_PATH = r"../tcspc/tryptamine_0p067odat261ex_350em_200uw_50ns_20mins_1938_may20_2025.sdt"
IRF_SDT_PATH = r"../tcspc/irf_261ex_261emt261ex_100uw_2143_nodblexc_05202025_correctcuvetteposition.sdt"

DECAY_BLOCK = 0
IRF_BLOCK = 0

# Folder where the LifeFit .txt files should be written.
OUTPUT_FOLDER = r"lifefit_txt"

# Auto name the outputs with none else name them 
DECAY_TXT_OUTPUT = None
IRF_TXT_OUTPUT = None



# Conversion taken from earlier code 


def load_sdt_decay(path, block=0):

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SDT file not found: {path}")

    with SdtFile(str(path)) as sdt:
        data = np.asarray(sdt.data[block], dtype=float)
        time = np.asarray(sdt.times[block], dtype=float)

    if data.ndim == 1:
        decay = data
    else:
        decay = data.reshape(-1, data.shape[-1]).sum(axis=0)

    time = np.squeeze(time)

    # sdtfile often returns seconds. A 50 ns window appears as 5e-8 seconds.
    if np.nanmax(time) < 1e-3:
        time_ns = time * 1e9
    else:
        time_ns = time

    if time_ns.shape[0] != decay.shape[0]:
        raise ValueError(
            f"Time and intensity length mismatch for {path}: "
            f"time={time_ns.shape[0]}, intensity={decay.shape[0]}"
        )

    return time_ns, decay


def default_output_name(input_sdt_path, label, output_folder):
    stem = Path(input_sdt_path).stem
    return Path(output_folder) / f"{stem}_{label}_lifefit_time_intensity.txt"


def write_lifefit_time_intensity_txt(input_sdt_path, output_txt_path, block=0):
    """
    Write a two-column LifeFit text file:
    """
    time_ns, intensity = load_sdt_decay(input_sdt_path, block=block)

    output_txt_path = Path(output_txt_path)
    output_txt_path.parent.mkdir(parents=True, exist_ok=True)

    table = np.column_stack([time_ns, intensity])
    np.savetxt(
        output_txt_path,
        table,
        fmt=["%.9f", "%.0f"],
        delimiter="\t",
        header="time_ns\tintensity",
        comments="",
    )

    return output_txt_path


def main():
    output_folder = Path(OUTPUT_FOLDER)

    decay_output = (
        Path(DECAY_TXT_OUTPUT)
        if DECAY_TXT_OUTPUT is not None
        else default_output_name(DECAY_SDT_PATH, "decay", output_folder)
    )
    irf_output = (
        Path(IRF_TXT_OUTPUT)
        if IRF_TXT_OUTPUT is not None
        else default_output_name(IRF_SDT_PATH, "irf", output_folder)
    )

    decay_written = write_lifefit_time_intensity_txt(
        DECAY_SDT_PATH,
        decay_output,
        block=DECAY_BLOCK,
    )
    irf_written = write_lifefit_time_intensity_txt(
        IRF_SDT_PATH,
        irf_output,
        block=IRF_BLOCK,
    )

    print("Wrote LifeFit time_intensity files:")
    print(f"  Decay: {decay_written}")
    print(f"  IRF:   {irf_written}")
    print()
    print("LifeFit UI settings:")
    print("  Mode: Analyze your own data")
    print("  Fileformat: time_intensity")
    print("  IRF: experimental IRF")


if __name__ == "__main__":
    main()
