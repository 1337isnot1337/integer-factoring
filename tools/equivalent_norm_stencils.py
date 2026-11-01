"""Exact extra norms from proper quadratic variable changes.

Supplied-prefix experiment: all norms, duplicates, and unusable values are
charged. No uniform improvement to relation supply or factoring is claimed.
"""
from math import gcd, prod
from time import perf_counter, process_time
from arithmetic_two_group import power
from batched_coprime_bucket_basis import BufferedBucketOrientedBasis
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import RationalRowBasis, verify_partial_identity
from oriented_norm_basis_lift import verify_lifted_oriented_basis
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from small_pivot_integer_row_basis import SmallPivotIntegerRowBasis
from sparse_oriented_power_search import transport_sparse_witnesses
from sparse_projected_signed_graph import eliminate_sparse_projected_signed_graph, verify_sparse_projected_signed_graph
from unit_presentation_trace import eliminate_unit_residual_columns_compact, verify_unit_presentation_trace

STENCILS = {
    'leading': (((1,0),(0,1)),),
    'leading_and_trailing': (((1,0),(0,1)),((0,-1),(1,0))),
    'four_norms': (((1,0),(0,1)),((0,-1),(1,0)),((1,0),(1,1)),((1,0),(-1,1))),
}


def proper_transform(form, matrix):
    if (len(form)!=3 or len(matrix)!=2 or any(len(row)!=2 for row in matrix)
            or any(not isinstance(v,int) for row in matrix for v in row)):
        raise ValueError('Supply three coefficients and a two-by-two integer matrix.')
    a,b,c = form
    (x,u),(y,v) = matrix
    if x*v-y*u != 1:
        raise ValueError('A proper variable change must have determinant one.')
    return (a*x*x+b*x*y+c*y*y,2*a*x*u+b*(x*v+y*u)+2*c*y*v,a*u*u+b*u*v+c*v*v)


def partial_represented_word(integer, represented, primes, records, prime_product=None):
    """Decode a primitive, possibly unreduced, norm using only known primes."""
    delta = -integer if integer%4==3 else -4*integer
    a,b,c = represented
    if a<1 or b*b-4*a*c!=delta or gcd(gcd(a,b),c)!=1 or gcd(a,integer)!=1:
        raise ValueError('Require a primitive positive form with leading norm coprime to the input.')
    if prime_product is None: prime_product = prod(primes)
    smooth = gcd(a,pow(prime_product,a.bit_length(),a))
    residual,remaining,parts,word,divisions = a//smooth,smooth,[],{},0
    for p in primes:
        exponent = 0
        while True:
            divisions += 1
            if remaining%p: break
            remaining //= p
            exponent += 1
        if exponent:
            parts.append((p,exponent))
            if p==2 and integer%4==1:
                if exponent!=1 or b%4!=2: raise ValueError('Inconsistent ramified-two branch.')
                word[p] = 1
            else:
                root,residue = records[p][0],b%(2*p)
                if residue not in (root,(-root)%(2*p)): raise ValueError('Inconsistent known-prime root.')
                word[p] = exponent if residue==root else -exponent
        if remaining==1: break
    assert remaining==1 and gcd(smooth,residual)==1
    key,orientation,residual_form = None,1,None
    if residual>1:
        if gcd(residual,2*delta)!=1: raise ValueError('Residual norm is not odd and unramified.')
        residue = b%(2*residual)
        canonical = min(residue,(-residue)%(2*residual))
        key,orientation = (residual,canonical),1 if residue==canonical else -1
        residual_form = reduce_form((residual,canonical,(canonical*canonical-delta)//(4*residual)))
    return {'smooth_part':smooth,'residual':residual,'factors':parts,'word':word,
            'key':key,'orientation':orientation,'residual_form':residual_form,'trial_divisions':divisions}


def analyze_power_prefix(reference, prefix, *, stencil='four_norms', max_columns=512):
    """Analyze all declared norms for a supplied prefix; never stop at a factor."""
    if stencil not in STENCILS or not 0<prefix<=len(reference['witnesses']):
        raise ValueError('Supply a known stencil and a nonempty available power prefix.')
    wall,cpu = perf_counter(),process_time()
    integer,delta = reference['integer'],reference['discriminant']
    identity = (1,delta%2,(delta%2-delta)//4)
    primes = list(reference['active_primes'])
    records = {int(p):(int(v[0]),tuple(v[1])) for p,v in reference['prime_records'].items()}
    columns,prime_product = {p:i for i,p in enumerate(primes)},prod(primes)
    state = BufferedBucketOrientedBasis(delta)
    indices,witnesses,rejected,direct = {},[],[],[]
    old_basis,old_keys = SmallPivotIntegerRowBasis(len(primes)),{}
    trials = norm_queries = maximum_bits = 0
    for power_index,original in enumerate(reference['witnesses'][:prefix]):
        p,exponent,target = original['prime'],original['exponent'],tuple(original['target'])
        for stencil_index,matrix in enumerate(STENCILS[stencil]):
            represented = proper_transform(target,matrix)
            norm_queries += 1
            maximum_bits = max(maximum_bits,represented[0].bit_length())
            divisor = gcd(represented[0],integer)
            if divisor!=1:
                event = {'power_index':power_index,'stencil_index':stencil_index,'norm':represented[0],'gcd':divisor}
                rejected.append(event)
                if divisor<integer: direct.append(event)
                continue
            split = partial_represented_word(integer,represented,primes,records,prime_product)
            trials += split['trial_divisions']
            key = split['key']
            if key is not None and key not in indices: indices[key] = state.add(*key)
            word = {columns[p]:exponent}
            for q,value in split['word'].items():
                j = columns[q]
                word[j] = word.get(j,0)-value
                if not word[j]: del word[j]
            witnesses.append({'prime':p,'exponent':exponent,'target':target,
                'power_index':power_index,'stencil_index':stencil_index,'represented_form':represented,
                'split':split,'base_word':sorted(word.items()),'key_index':indices[key] if key is not None else None})
            candidate = word if key is None else None
            if key is not None:
                signed = {j:split['orientation']*value for j,value in word.items()}
                if key not in old_keys: old_keys[key] = signed
                else:
                    candidate = dict(signed)
                    for j,value in old_keys[key].items(): candidate[j] = candidate.get(j,0)-value
            if candidate is not None: old_basis.add([candidate.get(j,0) for j in range(len(primes))])
    refinement = state.result(verify=False)
    rows = transport_sparse_witnesses(witnesses,len(primes),refinement['sparse_words'])
    graph = eliminate_sparse_projected_signed_graph(rows,len(primes),
        generator_count=len(primes)+len(refinement['atoms']),allow_extra_rows=True,verify=True)
    units = eliminate_unit_residual_columns_compact(rows,graph,verify=True,verify_graph=False)
    width = len(units['torsion_generator_indices'])
    factor = top = extraction = independent = None
    if width<=max_columns:
        selector = SmallPivotIntegerRowBasis(width)
        for row in units['torsion_rows']: selector.add(row)
        independent = selector.rows
        generators = [records[p][1] for p in primes]+[tuple(a['form']) for a in refinement['atoms']]
        selected = [generators[i] for i in units['torsion_generator_indices']]
        top,factor,extraction = first_useful_relation_top(selected,independent,identity,
            compose_forms_hnf,inverse_form,lambda f:ambiguous_factor(f,integer),verify_relations=True)
    return {'integer':integer,'discriminant':delta,'stencil':stencil,'power_prefix':prefix,
        'norm_queries':norm_queries,'trial_divisions':trials,'maximum_norm_bits':maximum_bits,
        'active_primes':primes,'prime_records':records,'keys':state.keys,'witnesses':witnesses,
        'rejected_norms':rejected,'direct_factors':direct,'refinement':refinement,
        'transported_rows':rows,'graph_reduction':graph,'unit_reduction':units,
        'old_base_rank':len(old_basis.rows),'core_columns':width,'independent_rows':independent,
        'independent_rank':len(independent) if independent is not None else None,
        'factor':direct[0]['gcd'] if direct else factor,'relation_factor':factor,
        'ambiguous_form':top,'extraction':extraction,
        'analysis_wall_seconds':perf_counter()-wall,'analysis_cpu_seconds':process_time()-cpu,
        'factoring_oracle_calls':0,'residual_primality_tests':0,
        'scope':'All stencil norms on a supplied power prefix; power generation is separately charged.'}


def verify_prefix_analysis(reference, result):
    """Replay powers, independent coordinate evaluations, ideals, and maps."""
    integer,delta = result['integer'],result['discriminant']
    assert integer==reference['integer'] and delta==reference['discriminant']
    identity = (1,delta%2,(delta%2-delta)//4)
    normalize = lambda r:{int(p):(int(v[0]),tuple(v[1])) for p,v in r.items()}
    records = normalize(result['prime_records'])
    assert records==normalize(reference['prime_records'])
    primes = result['active_primes']; columns = {p:i for i,p in enumerate(primes)}
    assert primes==reference['active_primes']
    prefix = result['power_prefix']; matrices = STENCILS[result['stencil']]
    assert 0<prefix<=len(reference['witnesses']) and result['norm_queries']==prefix*len(matrices)
    expected,rejected = {},[]
    for i,original in enumerate(reference['witnesses'][:prefix]):
        p,exponent,form = original['prime'],original['exponent'],tuple(original['target'])
        assert power(records[p][1],exponent,identity,compose_forms_hnf)==form
        for j,matrix in enumerate(matrices):
            (x,u),(y,v) = matrix
            assert x*v-y*u==1
            evaluate = lambda s,t:form[0]*s*s+form[1]*s*t+form[2]*t*t
            a,c = evaluate(x,y),evaluate(u,v)
            b = evaluate(x+u,y+v)-a-c
            d = gcd(a,integer)
            if d!=1: rejected.append({'power_index':i,'stencil_index':j,'norm':a,'gcd':d})
            else: expected[i,j] = ((a,b,c),p,exponent,form)
    assert result['rejected_norms']==rejected
    assert result['direct_factors']==[e for e in rejected if e['gcd']<integer]
    rebuilt,old_keys = [],{}
    old_basis = RationalRowBasis(len(primes))
    for witness in result['witnesses']:
        represented,p,exponent,target = expected.pop((witness['power_index'],witness['stencil_index']))
        assert tuple(witness['represented_form'])==represented and tuple(witness['target'])==target
        assert witness['prime']==p and witness['exponent']==exponent
        split = witness['split']; word = {int(q):v for q,v in split['word'].items()}
        smooth = prod(q**abs(v) for q,v in word.items())
        residual = represented[0]//smooth
        assert represented[0]==smooth*residual and split['smooth_part']==smooth
        assert split['residual']==residual and gcd(residual,prod(primes))==1
        for q,value in word.items():
            assert q in records and isinstance(value,int) and value
            if q==2 and integer%4==1: assert value==1 and represented[1]%4==2
            else: assert represented[1]%(2*q)==(records[q][0] if value>0 else -records[q][0])%(2*q)
        key_index = witness['key_index']
        if residual==1: assert key_index is None and split['key'] is None
        else:
            root = represented[1]%(2*residual); canonical = min(root,(-root)%(2*residual))
            assert tuple(split['key'])==(residual,canonical)
            assert tuple(result['keys'][key_index])==(residual,canonical)
            assert split['orientation']==(1 if root==canonical else -1)
        assert verify_partial_identity(integer,target,dict(split,word=word),records)
        row = {columns[p]:exponent}
        for q,value in word.items(): row[columns[q]] = row.get(columns[q],0)-value
        row = {j:v for j,v in row.items() if v}
        assert sorted(row.items())==[tuple(e) for e in witness['base_word']]
        candidate = row if key_index is None else None
        if key_index is not None:
            key = tuple(result['keys'][key_index])
            signed = {j:split['orientation']*v for j,v in row.items()}
            if key not in old_keys: old_keys[key] = signed
            else:
                candidate = dict(signed)
                for j,value in old_keys[key].items(): candidate[j] = candidate.get(j,0)-value
        if candidate is not None: old_basis.add([candidate.get(j,0) for j in range(len(primes))])
        if key_index is not None:
            for j,value in result['refinement']['sparse_words'][key_index]: row[len(primes)+j] = -split['orientation']*value
        rebuilt.append(row)
    assert not expected
    stored = [{int(j):v for j,v in row.items()} for row in result['transported_rows']]
    assert rebuilt==stored
    assert result['old_base_rank']==len(old_basis.rows)
    assert result['refinement']['discriminant']==delta
    assert [tuple(k) for k in result['refinement']['keys']]==[tuple(k) for k in result['keys']]
    ideals = verify_lifted_oriented_basis(result['refinement'])
    graph = verify_sparse_projected_signed_graph(stored,result['graph_reduction'])
    units = verify_unit_presentation_trace(stored,result['graph_reduction'],result['unit_reduction'],verify_graph=False)
    assert result['core_columns']==len(result['unit_reduction']['torsion_generator_indices'])
    if result['independent_rows'] is not None:
        selector = RationalRowBasis(result['core_columns'])
        for row in result['unit_reduction']['torsion_rows']: selector.add(row)
        assert selector.rows==result['independent_rows'] and len(selector.rows)==result['independent_rank']
    assert result['factor']==(result['direct_factors'][0]['gcd'] if result['direct_factors'] else result['relation_factor'])
    if result['factor'] is not None: assert 1<result['factor']<integer and integer%result['factor']==0
    return {'powers_replayed':prefix,'coordinate_norms_replayed':result['norm_queries'],
        'class_identities':len(rebuilt),'ideals':ideals,'graph':graph,'units':units,
        'complete_declared_stencil_replayed':True}
