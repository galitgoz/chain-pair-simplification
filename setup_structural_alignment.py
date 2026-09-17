"""Cache the official US-align Windows executable with source provenance."""
from pathlib import Path
import urllib.request,zipfile,hashlib,json,datetime,subprocess
root=Path(__file__).resolve().parent
target=root/'tools'/'usalign'; target.mkdir(parents=True,exist_ok=True)
url='https://seq2fun.dcmb.med.umich.edu/US-align/bin/module/USalignWin64.zip'
archive=target/'USalignWin64.zip'
if not archive.exists():
    with urllib.request.urlopen(url,timeout=90) as response: archive.write_bytes(response.read())
with zipfile.ZipFile(archive) as z:
    names=[s for s in z.namelist() if s.lower().endswith('.exe')]
    print(names)
    for name in names:
        (target/Path(name).name).write_bytes(z.read(name))
exe=next(target.glob('*.exe'))
version=subprocess.run([str(exe),'-v'],capture_output=True,text=True).stdout
(target/'source.json').write_text(json.dumps(dict(url=url,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),version=version),indent=2))
print(version)
