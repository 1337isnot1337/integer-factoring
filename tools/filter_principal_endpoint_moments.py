#!/usr/bin/env python3
"""Exact first/second endpoint moments of every recorded principal profile."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import lzma
from pathlib import Path
import time


def moments(values):
    assert values and all(int(v)==v for v in values)
    values = list(map(int,values))
    n,total = len(values),sum(values)
    return Fraction(total,n),Fraction(n*sum(v*v for v in values)-total*total,n*n)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root,started = Path(__file__).resolve().parents[1],time.monotonic()
    manifest = json.loads(args.manifest.read_text())
    path = root/manifest['input_archive']
    assert sha256(path.read_bytes()).hexdigest()==manifest['input_archive_sha256']
    packet = json.loads(lzma.decompress(path.read_bytes()))
    assert packet['summary']['status']=='complete_registered_screen'
    records = [json.loads(line) for line in packet['profiles_jsonl'].splitlines()]
    assert len(records)==45292
    gap_signs,multiplicities,variances = Counter(),Counter(),Counter()
    minimum_gap,minimum_variance,candidates = None,None,[]
    count,saturated_simple = 0,0
    for index,record in enumerate(records):
        ex,ey = record['physical_energies']
        mx,vx = moments(ex)
        my,vy = moments(ey)
        for strong in (0,1):
            r,s = record['orders'][strong],record['orders'][1-strong]
            mr,vr = moments(r)
            gap,variance = mx+my-mr-min(s),vr-vx-vy
            multiplicity = s.count(min(s))
            detail = {'record_index':index,'source_case_index':record['source_case_index'],
                      'alphabet_x':record['alphabet_x'],'alphabet_y':record['alphabet_y'],
                      'strong_place':strong,'first_order_gap':str(gap),
                      'weak_minimum_multiplicity':multiplicity,'variance_difference':str(variance)}
            gap_signs['negative' if gap<0 else 'positive' if gap>0 else 'zero'] += 1
            if minimum_gap is None or gap < minimum_gap[0]:
                minimum_gap = gap,detail
            if gap==0:
                multiplicities[multiplicity] += 1
                if multiplicity==1:
                    saturated_simple += 1
                    variances[str(variance)] += 1
                    if minimum_variance is None or variance < minimum_variance[0]:
                        minimum_variance = variance,detail
            if gap<0 or (gap==0 and multiplicity==1 and variance<=0):
                candidates.append(detail)
            count += 1
    assert count==manifest['registered_endpoint_states']==90584
    report = {'status':'candidate_requires_exact_profile_replay' if candidates else 'complete_recorded_profile_filter',
              'registration':manifest,'recorded_profiles':len(records),'endpoint_states':count,
              'first_order_gap_signs':dict(gap_signs),'minimum_first_order_gap':minimum_gap[1],
              'zero_gap_weak_multiplicities':dict(multiplicities),'saturated_simple_states':saturated_simple,
              'saturated_simple_variance_histogram':dict(sorted(variances.items())),
              'minimum_saturated_simple_variance':minimum_variance[1] if minimum_variance else None,
              'candidates':candidates,'seconds':time.monotonic()-started,
              'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':'Exact arithmetic on every recorded inherited profile and physical energy '
                      'vector. This does not reconstruct all native profiles, optimize free '
                      'heights or prove a global endpoint/factoring bound.'}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','endpoint_states','first_order_gap_signs',
        'saturated_simple_states','minimum_saturated_simple_variance')}))


if __name__=='__main__':
    main()
