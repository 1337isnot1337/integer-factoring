#!/usr/bin/env python3
"""Different divided-power alphabets with complete profiles and asymmetric cost."""
import argparse
from hashlib import sha256
import json
from math import cos, exp, log, pi, sin
from pathlib import Path

from scipy.optimize import minimize

from coherent_divided_powers import paired_top
from coherent_inverse_systems import inverse_system
from coherent_unequal_compounds import profiles
from screen_unequal_coherent_hot_codes import exact_point, support


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, required=True)
    parser.add_argument('--left-power', type=int, default=1)
    parser.add_argument('--right-power', type=int, default=2)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 2 <= args.modes <= 16 or args.modes % 2:
        parser.error('Use a positive even mode count through sixteen.')
    codes = [inverse_system([sum(1 << word for word in paired_top(args.modes, r))], args.modes)
             for r in (args.left_power, args.right_power)]
    x, y = codes
    print(json.dumps(dict(q_x=len(x['embedding']), q_y=len(y['embedding']))), flush=True)
    result = profiles(x['embedding'], x['retraction'], y['embedding'], y['retraction'], args.modes)
    sides = [result['infinity_orders_t'], result['zero_orders_t']]
    ex, ey = ([2*h for h in code['heights']] for code in codes)
    objective = lambda values: support(sides, ex, ey, exp(values[0]), cos(values[1]), sin(values[1]), 1)
    trials = []
    for angle in (pi/6, pi/4, pi/3):
        fit = minimize(objective, [log(.3), angle], method='L-BFGS-B',
                       bounds=((log(.005), log(2)), (.05, pi/2-.05)))
        trials.append(dict(support=float(fit.fun), scale=exp(float(fit.x[0])),
                           angle=float(fit.x[1]), completed=bool(fit.success)))
    best = min(trials, key=lambda trial: trial['support'])
    # Every retained rational test is independently meaningful, even when
    # the floating-point optimizer has not found the global best branch.
    rational = []
    for p, q, h in ((3, 4, 5), (4, 3, 5), (5, 12, 13), (12, 5, 13), (20, 21, 29), (21, 20, 29)):
        rational.append(exact_point(sides, ex, ey, best['scale'], p, q, h))
    names = (Path(__file__).name, 'coherent_divided_powers.py', 'coherent_inverse_systems.py',
             'coherent_unequal_compounds.py', 'coherent_hot_compounds.py',
             'coherent_hot_compounds.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'screen_unequal_coherent_hot_codes.py',
             'check_binary_coefficient_entropy.py', 'certify_hot_charge_pressure.py')
    report = dict(modes=args.modes, powers=[args.left_power, args.right_power], codes=codes,
                  profile=result, numerical_trials=trials, best_numerical=best,
                  exact_rational_tests=rational,
                  source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in names},
                  scope='One finite pair and six exact rational points, not a global asymmetric optimum.')
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps(dict(best=best, exact_below_one=any(point['strict_support_below_one'] for point in rational))))


if __name__ == '__main__':
    main()
