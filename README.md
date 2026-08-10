# DNA Short Reads Mapping

This project maps short DNA reads to reference sequences using a binary k-mer index and Hamming distance. It generates read-mapping positions, reference-sequence coverage statistics, and match-distribution plots.

## Project Structure

- `scripts/`: Main program and analysis modules.
- `tests/`: Unit tests.
- `data/references/`: Reference-sequence FASTA files.
- `data/reads/`: Short-read FASTA files.
- `results/output/`: Generated results. The directory is retained, but its generated files are excluded from Git.
- `docs/`: Project report, answers to the project questions, and their PDF versions.
- `launchers/`: macOS shortcuts for opening a project terminal and running the tests.

## Installation

Creating a virtual environment in the project root is recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Usage

Place reference FASTA files in `data/references/` and read FASTA files in `data/reads/`. Then run the following command from the project root:

```bash
python3 scripts/main.py --k 15 refs0.fa reads0.fa Refs0_Reads0
```

All generated files are written to `results/output/`. The output prefix should not contain path characters such as `/` or `.`. The supported k-mer length is 5–32; values between 10 and 30 are recommended.

## Running the Tests

After activating an environment with the required dependencies, run:

```bash
PYTHONPATH=scripts python3 -m unittest discover -s tests -v
```

On macOS, you can also double-click `launchers/RunTest.command`. If a `.venv` directory exists in the project root, the launcher uses it automatically. `launchers/OpenTerminal.command` opens a terminal in the project root and displays an example command.

## Documentation

- Main report: `docs/Report_Zhaokun_Wang.pdf`
- Figure explanations and answers to the project questions: `docs/Answer_to_Question_Zhaokun_Wang.pdf`

## GitHub Repository Notes

This directory contains only the project files suitable for version control. Generated results, historical result snapshots, the old ZIP archive, PyCharm settings, Python caches, and macOS metadata are not included.

The `.gitignore` covers macOS AppleDouble files (`._*`), `.DS_Store`, Python virtual environments, test caches, IDE settings, local environment files, and generated files under `results/output/`.
