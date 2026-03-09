#!/usr/bin/env python3

import json
import os
import requests
import pandas as pd
import argparse
import logging
import time
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed



SAMPLE_URL = "https://www.ebi.ac.uk/biosamples"
AUTH_URL = "https://www.ebi.ac.uk/ena/submit/webin/auth/token"

ENA_USERNAME = os.getenv("ENA_USERNAME")
ENA_PASSWORD = os.getenv("ENA_PASSWORD")

PROJECT_NAME = "Microflora Danica"
MAX_WORKERS = 4
SLEEP_SECONDS = 0.1

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)


def clean_value(v):
    if v is None:
        return None
    if pd.isna(v):
        return None
    if isinstance(v, str):
        v = v.strip()
        if not v or v.lower() == "not provided":
            return None
    return v


def auth_token() -> str:
    credentials = {
        "authRealms": ["ENA"],
        "username": ENA_USERNAME,
        "password": ENA_PASSWORD,
    }
    headers = {"Accept": "*/*", "Content-Type": "application/json"}
    r = requests.post(AUTH_URL, json=credentials, headers=headers, timeout=30)
    r.raise_for_status()
    return r.text.strip()


def fetch_biosample_json(accession: str, session: requests.Session) -> dict:
    r = session.get(f"{SAMPLE_URL}/samples/{accession}.json", timeout=30)
    r.raise_for_status()
    return r.json()


def submit_structureddata(accession: str, token: str, data: dict, session: requests.Session):
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {token}",
    }
    url = f"{SAMPLE_URL}/structureddata/{accession}"
    r = session.put(url, json=data, headers=headers, timeout=60)
    logging.info(f"StructuredData PUT {accession}: {r.status_code}")
    r.raise_for_status()


def has_project_curation(accession: str, token: str, session: requests.Session) -> bool:
    headers = {
        "Accept": "application/json",
        "Authorization": f"Bearer {token}",
    }
    url = f"{SAMPLE_URL}/samples/{accession}/curationlinks"
    r = session.get(url, headers=headers, timeout=30)
    r.raise_for_status()
    data = r.json()

    for item in data.get("items", []):
        for attr in item.get("curation", {}).get("attributesPost", []):
            if attr.get("type") == "project" and attr.get("value") == PROJECT_NAME:
                return True
    return False


def submit_project_curation(accession: str, token: str, session: requests.Session):
    curation = {
        "sample": accession,
        "curation": {
            "attributesPre": [],
            "attributesPost": [
                {"type": "project", "value": PROJECT_NAME}
            ],
            "externalReferencesPre": [],
            "externalReferencesPost": [],
        }
    }
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "Authorization": f"Bearer {token}",
    }
    url = f"{SAMPLE_URL}/samples/{accession}/curationlinks"
    r = session.post(url, json=curation, headers=headers, timeout=30)
    logging.info(f"Curation POST {accession}: {r.status_code}")
    r.raise_for_status()


# -------------------- METADATA --------------------

def parse_metadata(csv_path: str) -> dict:
    df = pd.read_csv(csv_path)
    meta = {}

    for _, row in df.iterrows():
        acc = str(row["BioSample"]).strip()

        meta[acc] = {
            "experimental_data": {
                "project_identifier": clean_value(row.get("project_id")),
                "extraction_method": clean_value(row.get("extraction_method")),
                "library_method": clean_value(row.get("library_method")),
            },
            "marker_gene_operons": {
                "eukaryote_rRNA_operon": clean_value(row.get("EUK_operon")),
                "bacteria_rRNA_operon": clean_value(row.get("BAC_operon")),
                "UMI_16SrRNA": clean_value(row.get("UMI_16SrRNA")),
            },
            "mfd_location_data": {
                "ontology": "Microflora Danica curated environmental ontology",
                "01_mfd_sampletype": clean_value(row.get("mfd_sampletype")),
                "02_mfd_areatype": clean_value(row.get("mfd_areatype")),
                "03_mfd_hab1": clean_value(row.get("mfd_hab1")),
                "04_mfd_hab2": clean_value(row.get("mfd_hab2")),
                "05_mfd_hab3": clean_value(row.get("mfd_hab3")),
                "coordinate_reliability": clean_value(row.get("coords_reliable")) or "masked",
            },
        }

    return meta


def parse_genome_metadata(csv_path: str) -> dict:
    df = pd.read_csv(csv_path)
    meta = {}

    for _, row in df.iterrows():
        acc = str(row["BioSample"]).strip()
        meta[acc] = {
            "genome_quality":{
                "completeness score": clean_value(row.get("Completeness")),
                "contamination score": clean_value(row.get("Contamination")),
                "completeness software": "CheckM2",
            }

        }

    return meta

def build_structured_block(block_type, webin_id, fields, iri=False):
    row = {}
    for k, v in fields.items():
        v = clean_value(v)
        if v is None:
            continue
        row[k] = {
            "value": v,
            "iri": f"https://www.ebi.ac.uk/ena/browser/view/{v}" if iri else None,
        }

    if not row:
        return None

    return {
        "domain": None,
        "schema": None,
        "type": block_type,
        "webinSubmissionAccountId": webin_id,
        "content": [row],
    }


def build_structureddata(accession, sample_json, metadata, is_genome=False):
    meta = metadata[accession]
    blocks = []

    if is_genome:
        blocks.append(build_structured_block("Genome quality", ENA_USERNAME, meta["genome_quality"]))
    else:
        for block in (
            build_structured_block("Experiment", ENA_USERNAME, meta["experimental_data"]),
            build_structured_block("Marker genes", ENA_USERNAME, meta["marker_gene_operons"], iri=True),
            build_structured_block("Location ontology", ENA_USERNAME, meta["mfd_location_data"]),
        ):
            if block:
                blocks.append(block)

    return {
        "accession": accession,
        "create": sample_json.get("create"),
        "update": sample_json.get("update"),
        "data": blocks,
    }


# -------------------- WORKER --------------------

def process_accession(acc, metadata, token, write_json, submit, is_genome=False):
    session = requests.Session()
    result = {
        "accession": acc,
        "structureddata": False,
        "curation": False,
        "error": None,
    }

    try:
        sample = fetch_biosample_json(acc, session)
        data = build_structureddata(acc, sample, metadata, is_genome=is_genome)

        if write_json:
            with open(f"structureddata_{acc}.json", "w") as f:
                json.dump(data, f, indent=2)

        if submit:
            # ---- structured data ----
            time.sleep(SLEEP_SECONDS)
            submit_structureddata(acc, token, data, session)
            result["structureddata"] = True

            # ---- curation ----
            time.sleep(SLEEP_SECONDS)
            if not has_project_curation(acc, token, session):
                submit_project_curation(acc, token, session)
            result["curation"] = True

    except Exception as e:
        result["error"] = str(e)

    finally:
        session.close()

    return (
        acc,
        write_json,
        result["structureddata"],
        result["curation"],
        result["error"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-csv", required=True)
    parser.add_argument("--is-genome", action="store_true", help="metadata is genome", default=False)
    parser.add_argument("--accessions-file", required=True)
    parser.add_argument("--json", action="store_true", help="write structured data JSON")
    parser.add_argument("--submit", action="store_true", help="submit to BioSamples")
    args = parser.parse_args()

    if not args.json and not args.submit:
        parser.error("At least one of --json or --submit is required")

    if not ENA_USERNAME or not ENA_PASSWORD:
        raise EnvironmentError("ENA_USERNAME and ENA_PASSWORD must be set")

    if args.is_genome:
        metadata = parse_genome_metadata(args.data_csv)
    else:
        metadata = parse_metadata(args.data_csv)
    with open(args.accessions_file) as f:
        accessions = [l.strip() for l in f if l.strip()]

    token = auth_token()

    logging.info(f"Processing {len(accessions)} samples with {MAX_WORKERS} workers")

    results_fh = open("submission_results.tsv", "a", newline="")
    results_writer = csv.writer(results_fh, delimiter="\t")

    # write header if file is empty
    if results_fh.tell() == 0:
        results_writer.writerow(["accession", "json", "submit", "curate", "error"])


    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as exe:
        futures = [
            exe.submit(
                process_accession,
                acc,
                metadata,
                token,
                args.json,
                args.submit,
                is_genome=args.is_genome,
            )
            for acc in accessions
        ]

        for future in as_completed(futures):
            acc, json_ok, submit_ok, curate_ok, error = future.result()

            results_writer.writerow([
                acc,
                json_ok,
                submit_ok,
                curate_ok,
                error or "",
            ])
            results_fh.flush()

            if error:
                logging.error(f"[FAIL] {acc}: {error}")
            else:
                logging.info(f"[OK] {acc}")

    logging.info("All done.")




if __name__ == "__main__":
    main()
