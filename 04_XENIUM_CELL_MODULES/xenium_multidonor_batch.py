# PUBLIC RELEASE CANDIDATE COPY
# Source: tka_inflammation_redesign/STAGE7G_17DONOR_ECM_SPATIAL/12_batch_process_stage7g.py
# Source SHA256: 0f57aa8d525083d7d4a44b29d6d2e36847ca902d50b0d792a308f7f85bfa1c34
# Purpose: Checkpointed public Xenium donor processing
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env python3
"""Checkpointed Stage 7G donor processing using the verified Stage 7F pipeline."""
from __future__ import annotations
import argparse, csv, json, os, shutil, subprocess, sys, time
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--coverage', type=Path, required=True)
    ap.add_argument('--donors', nargs='+', required=True)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--chunk-mb', type=int, default=32)
    ap.add_argument('--tmp-dir', type=Path, default=Path('.xenium_tmp'))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    py = os.environ.get('XENIUM_PYTHON', sys.executable)
    analysis_root = Path(__file__).resolve().parents[1]
    downloader = analysis_root / 'STAGE7F_MULTIDONOR_EXECUTION' / 'download_range.py'
    extractor = analysis_root / 'STAGE7F_MULTIDONOR_EXECUTION' / 'XENIUM_MINIMAL_EXTRACTION_PIPELINE.py'
    processor = analysis_root / 'STAGE7F_MULTIDONOR_EXECUTION' / '13_process_donor.py'
    log = args.out / 'STAGE7G_BATCH_EXECUTION.log'
    rows = {r['donor_id']: r for r in csv.DictReader(args.manifest.open())}
    with log.open('a', encoding='utf-8') as lf:
        for donor in args.donors:
            if donor not in rows:
                lf.write(f'{donor}\tNOT_IN_MANIFEST\n'); continue
            r = rows[donor]
            result = args.out / f'{donor}_results'
            summary = result / f'{donor}_summary.json'
            if summary.exists() and summary.stat().st_size > 0:
                lf.write(f'{donor}\tSKIP_ALREADY_COMPLETE\n'); continue
            temp = args.tmp_dir / f'stage7g_{donor}'
            minimal = args.tmp_dir / f'stage7g_{donor}_minimal'
            temp.mkdir(parents=True, exist_ok=True)
            archive = temp / f'{donor}_full.tar.gz'
            dmanifest = args.out / f'{donor}_download_manifest.json'
            emanifest = args.out / f'{donor}_extraction_manifest.json'
            lf.write(f'{donor}\tSTART\t{time.strftime("%Y-%m-%dT%H:%M:%S%z")}\n'); lf.flush()
            cmd = [py, str(downloader), '--url', r['source_url'], '--output', str(archive), '--size', r['archive_size_bytes'], '--workers', str(args.workers), '--chunk-mb', str(args.chunk_mb), '--expected-md5', r['md5'], '--manifest', str(dmanifest)]
            if subprocess.run(cmd).returncode != 0:
                lf.write(f'{donor}\tDOWNLOAD_FAILED\n'); lf.flush(); continue
            cmd = [py, str(extractor), '--archive', str(archive), '--expected-md5', r['md5'], '--out', str(minimal), '--manifest', str(emanifest)]
            if subprocess.run(cmd).returncode != 0:
                lf.write(f'{donor}\tEXTRACTION_FAILED\n'); lf.flush(); continue
            cmd = [py, str(processor), '--input', str(minimal), '--donor', donor, '--coverage', str(args.coverage), '--out', str(result)]
            if subprocess.run(cmd).returncode != 0:
                lf.write(f'{donor}\tPROCESSING_FAILED\n'); lf.flush(); continue
            archive.unlink(missing_ok=True)
            for child in minimal.glob('*'):
                if child.is_file() or child.is_symlink(): child.unlink()
            minimal.rmdir()
            lf.write(f'{donor}\tCOMPLETE\n'); lf.flush()

if __name__ == '__main__':
    main()
