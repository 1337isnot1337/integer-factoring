#!/usr/bin/env python3
"""Reconstruct one exact certificate archive from checked repository segments."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path


def assemble(manifest_path, output=None):
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 'sha256-segmented-certificate-v1':
        raise ValueError('Unsupported certificate segment manifest.')
    name = manifest['filename']
    if Path(name).name != name or name in ('', '.', '..'):
        raise ValueError('The certificate filename must be a single path component.')
    output = Path(output) if output is not None else manifest_path.parent/name
    if output.exists():
        digest = sha256()
        with output.open('rb') as stream:
            for chunk in iter(lambda:stream.read(1048576),b''):
                digest.update(chunk)
        if output.stat().st_size==manifest['bytes'] and digest.hexdigest()==manifest['sha256']:
            return {'output':str(output),'bytes':manifest['bytes'],'sha256':digest.hexdigest(),'already_complete':True}
        raise FileExistsError('The output exists with a different length or checksum.')
    temporary = output.with_name(output.name+'.assembling')
    total, digest, created = 0, sha256(), False
    try:
        with temporary.open('xb') as destination:
            created = True
            for part in manifest['segments']:
                part_name = part['filename']
                if Path(part_name).name != part_name or part_name in ('', '.', '..'):
                    raise ValueError('A segment filename must be a single path component.')
                if part['offset'] != total:
                    raise ValueError('Certificate segments are not contiguous and ordered.')
                part_hash, size = sha256(), 0
                with (manifest_path.parent/part_name).open('rb') as source:
                    for chunk in iter(lambda:source.read(1048576),b''):
                        destination.write(chunk)
                        digest.update(chunk)
                        part_hash.update(chunk)
                        size += len(chunk)
                if size != part['bytes'] or part_hash.hexdigest()!=part['sha256']:
                    raise ValueError('A certificate segment has the wrong length or checksum.')
                total += size
        if total != manifest['bytes'] or digest.hexdigest()!=manifest['sha256']:
            raise ValueError('The reconstructed certificate has the wrong length or checksum.')
        os.link(temporary,output)
        temporary.unlink()
    except BaseException:
        if created and temporary.exists():
            temporary.unlink()
        raise
    return {'output':str(output),'bytes':total,'sha256':digest.hexdigest(),'already_complete':False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest',type=Path)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    print(json.dumps(assemble(args.manifest,args.output),indent=2))


if __name__=='__main__':
    main()
