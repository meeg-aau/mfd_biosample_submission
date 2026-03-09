# BioSample Structured Data Submission

Internal scripts to submit structured data and curation links to BioSamples.

## Prerequisites

- Python 3.6+
- `requests`
- `pandas`

Install dependencies:

```bash
pip install requests pandas
```

## Setup

The script requires ENA (Webin) credentials to authenticate with the BioSamples API. Set them as environment variables:

```bash
export ENA_USERNAME="Webin-XXXX"
export ENA_PASSWORD="your_password"
```

## Usage

The script processes a list of BioSample accessions and metadata from a CSV file to generate and/or submit structured data.

### Basic Command

```bash
python make_biosample_structureddata_json.py --data-csv {metadata_csv} --accessions-file {accessions_txt} --json --submit
```

### Options

| Argument | Description |
| :--- | :--- |
| `--data-csv` | Path to the CSV file containing sample metadata. |
| `--accessions-file` | Path to a line-separated text file with BioSample accessions (e.g., SAMEA12345). |
| `--json` | Write the generated structured data JSON to local files (named `structureddata_{accession}.json`). |
| `--submit` | Submit the structured data and a "project" curation link to BioSamples. |
| `--is-genome` | Indicate that the samples relate to a genome (uses different metadata parsing logic). |

### Output

- **JSON Files:** If `--json` is used, individual JSON files are created for each accession.
- **Results Log:** Execution results are appended to `submission_results.tsv`, including status for JSON generation, submission, and curation.

## Example

```bash
python make_biosample_structureddata_json.py \
  --data-csv tests/fixtures/mfd_metadata.csv \
  --accessions-file tests/fixtures/mfd_biosamples.txt \
  --json \
  --submit
```
