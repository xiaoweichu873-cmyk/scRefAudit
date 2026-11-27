"""Download a public input and verify the retained source checksum when present."""
import sys,json,hashlib,urllib.request
from pathlib import Path

def main(url,destination):
    dest=Path(destination);dest.parent.mkdir(parents=True,exist_ok=True)
    source=Path(str(dest)+'.source.json')
    known=json.loads(source.read_text()) if source.exists() else None
    if known and known['url']!=url:raise ValueError('URL differs from retained provenance')
    temp=Path(str(dest)+'.part')
    with urllib.request.urlopen(url,timeout=90) as response,temp.open('wb') as out:
        while chunk:=response.read(4*1024*1024):out.write(chunk)
    digest=hashlib.sha256(temp.read_bytes()).hexdigest()
    if known and digest!=known['sha256']:raise ValueError('Downloaded bytes fail retained SHA-256')
    temp.replace(dest)
    if not known:source.write_text(json.dumps({'url':url,'path':str(dest),'bytes':dest.stat().st_size,'sha256':digest},indent=2))
    print(dest,dest.stat().st_size,digest)

if __name__=='__main__':main(*sys.argv[1:])
