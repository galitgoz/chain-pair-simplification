from pathlib import Path
import urllib.request, hashlib, json, datetime
from concurrent.futures import ThreadPoolExecutor
ROOT = Path(__file__).resolve().parent
CACHE = ROOT / 'data' / 'pdb_cache'
CACHE.mkdir(exist_ok=True)
IDS = '1o7j 1hfj 1qd1 1toh 4eca 1d9q 3ntx 1wls 2eq5 2zsk 1zq1 3jq0 2fep 4cea 107j'.split()
def fetch(pid):
    url = f'https://files.rcsb.org/download/{pid.upper()}.cif'
    path = CACHE / f'{pid}.cif'
    record = dict(pdb_id=pid, url=url, retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    try:
        if not path.exists():
            with urllib.request.urlopen(url, timeout=90) as response:
                path.write_bytes(response.read())
        record.update(status='available', sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    except Exception as exc:
        record.update(status='unavailable', error=str(exc))
    return record
if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=6) as pool:
        records = list(pool.map(fetch, IDS))
    (CACHE / 'sources.json').write_text(json.dumps(records, indent=2))
    print(json.dumps(records, indent=2))
    path = ROOT / 'papers' / 'chain_pair_simplification.pdf'
    if not path.exists():
        with urllib.request.urlopen('https://arxiv.org/pdf/1409.2457', timeout=90) as response:
            path.write_bytes(response.read())
