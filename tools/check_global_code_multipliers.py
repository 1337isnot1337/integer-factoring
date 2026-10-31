#!/usr/bin/env python3
"""Regenerate every global code family and replay the complete multiplier cost."""

import argparse
import json
from pathlib import Path

from boolean_multiplier_profiles import complete_multiplier_automaton
from check_symmetric_filter_catalogues import physical_entropy
from prime_state_filter_pressure import gaussian_pressure
from search_global_code_multipliers import (
    CODE_MATRICES, GlobalCodePressure, build_families, digest,
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, default=root / 'evidence' / 'global-code-multiplier-searches.json')
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'global-code-multiplier-replay.json')
    args = parser.parse_args()
    record = json.loads(args.record.read_text())
    assert record['status'].startswith('All requested')
    assert json.loads(json.dumps(CODE_MATRICES)) == record['source_code_matrices']
    for rank, stored in record['automata'].items():
        automaton = complete_multiplier_automaton(int(rank))
        assert digest(automaton) == stored['sha256']
        assert len(automaton['continuous_profiles']) == stored['continuous_profiles']
        assert sum(len(tail['cycle']) for tail in automaton['singleton_tails']) == stored['negative_gap_periodic_rows']
    rows = []
    supports = 0
    for stored in record['optimized_runs']:
        n, rank = stored['n'], stored['rank']
        families = build_families(n, rank, record['seed'])
        assert digest(families) == stored['support_families_sha256']
        assert json.loads(json.dumps(families)) == stored['support_families']
        for family in families:
            all_energies = [value for group in family['actual_mode_energy_groups'] for value in group]
            all_energies += list(family['spectator_mode_energies'])
            assert sorted(all_energies) == list(range(1, 2 * n, 2))
            assert sum(map(len, family['actual_mode_energy_groups'])) == family['active_modes']
        evaluator = GlobalCodePressure(n, rank, families)
        relative, selected, pieces = evaluator.relative_cost(stored['tau'], stored['scaled_charge'])
        assert abs(relative - stored['relative_log_cost']) < 1e-9
        assert selected.tolist() == stored['selected_catalogue_indices']
        assert json.loads(json.dumps(evaluator.chosen_detail)) == stored['selected_multiplier']
        residual = gaussian_pressure(2, n, stored['charge'], stored['scale'], stored['rho']) + relative
        critical = gaussian_pressure(2, n, stored['charge'], stored['scale']) + relative
        assert abs(residual) < 1e-7 and abs(critical - stored['full_critical_pressure']) < 1e-8
        entropy = physical_entropy(2, n, stored['charge'], stored['scale'], stored['rho'])
        assert entropy > 1e-8
        supports += len(families) * (1 << rank)
        rows.append({
            'n': n, 'rank': rank, 'code_supports': len(families) * (1 << rank),
            'rho': stored['rho'], 'root_residual': residual,
            'full_critical_pressure': critical, 'physical_entropy': entropy,
            'selected_kind': evaluator.chosen_detail['kind'],
        })
        evaluator.release()
    result = {
        'status': 'All explicit global code supports regenerated and complete multiplier-envelope roots replayed.',
        'scope': 'Exact support/multiplier reductions; the outer scale and charge searches remain exploratory.',
        'total_code_supports_in_completed_runs': supports,
        'complete_root_replays': rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
