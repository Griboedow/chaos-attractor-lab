"""Cache retrieved primary-source HTML outside the repository."""
from pathlib import Path
import datetime,hashlib,json,urllib.request
folder=Path('sources');folder.mkdir(exist_ok=True)
sources=[('lehmann2021','https://arxiv.org/html/2012.12562v2'),
         ('thibault2021','https://arxiv.org/html/1711.01851v2'),
         ('soma2026','https://arxiv.org/html/2410.14104v2'),
         ('eisenmann2026','https://arxiv.org/html/2603.12955v1')]
records=[]
for name,url in sources:
    data=urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'chaos-attractor-lab reproducibility audit'}),timeout=60).read()
    path=folder/(name+'.html');path.write_bytes(data)
    records.append(dict(name=name,url=url,local_path=str(path.resolve()),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
Path('sources.json').write_text(json.dumps(dict(sources=records),indent=2)+'\n')
print(json.dumps(records))
