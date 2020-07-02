# MITMAPS

An open-source software package developed for fast analysis of rheological data obtained using the multi-tone MAPS experimental protocol. See the following paper for details about the MAPS experimental protocol:

K. R. Lennon, M. Geri, G. H. McKinley, and J. W. Swan, ["Medium Amplitude Parallel Superposition (MAPS) Rheology, Part 2: Experimental Protocols and Data Analysis"](https://arxiv.org/abs/2006.09465), arXiv preprint (2020).

## Compatibility

MITMAPS has been tested using [Python](https://www.python.org/downloads/) 3.7.6 with [NumPy](https://numpy.org/install/) 1.18.1, [Matplotlib](https://matplotlib.org/downloads.html) 3.1.3, and [SciPy](https://www.scipy.org/install.html) 1.4.1.

The data used to test this software was obtained from the TRIOS v5.0.0 software supplied by TA Instruments. File I/O details are confined to the source file `mapsio.py` and depend on the specific file format exported by the TRIOS software. The remainder of the software functions independently from the specifics of the input data file format.

## Installation

To install the MITMAPS software package, download all of the files in the `src` directory and place these files in the same local directory.

## Contributing
Inquiries about the MITMAPS software can be directed to krlennon[at]mit.edu.

## License
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/)
