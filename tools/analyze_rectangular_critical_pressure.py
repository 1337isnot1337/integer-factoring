#!/usr/bin/env python3
"""Analyze saved rectangular profiles without making new profile-oracle calls.

All pressure values and starts are floating research data, not certificates.
The tool uses only the standard library. It verifies each saved Gibbs root
and the physical-energy identities of its transported starts.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path


CRITICAL_R = math.sqrt(2)


def softplus(value):
    return max(0.0, value) + math.log1p(math.exp(-abs(value)))


def logsumexp(values):
    largest = max(values)
    return largest + math.log(sum(math.exp(value - largest) for value in values))


def probabilities(energies, R):
    normalizer = logsumexp([-value / R for value in energies])
    return [math.exp(-value / R - normalizer) for value in energies]


def entropy(values):
    return -sum(value * math.log(value) for value in values if value)


def physical_energies(M, theta):
    stride = (len(theta) - 5) // 2
    if stride not in (7, 10) or len(theta) != 2 * stride + 5:
        raise ValueError("Expected 19 or 25 parameters.")
    infinity_charge, zero_charge = sum(theta[:2]), sum(theta[2:4])
    gx = [(4 * a - 2 * M + 2) * theta[-1] + infinity_charge + zero_charge
          for a in range(M)]
    gy = [(4 * a - 2 * M + 2) * theta[-1] + infinity_charge - zero_charge
          for a in range(M)]
    vx = [sum(theta[4 + side * stride + j] for side in (0, 1)) for j in range(3)]
    if stride == 10:
        vx = [value + sum(theta[4 + side * stride + 7 + j] for side in (0, 1))
              for j, value in enumerate(vx)]
    uy = [sum(theta[4 + side * stride + leg + j]
              for side in (0, 1) for leg in (3, 5)) for j in range(2)]
    cx = [0., gx[0] + vx[0], gx[-1] + vx[1], gx[0] + gx[-1] + vx[2]]
    cy = [0., gy[0] + uy[0], gy[-1] + uy[1]]
    return gx, gy, cx, cy


def partition_and_entropy(energies, R):
    gx, gy, cx, cy = energies
    logZ = sum(softplus(-g / R) for g in gx[1:-1] + gy[1:-1])
    logZ += logsumexp([-g / R for g in cx]) + logsumexp([-g / R for g in cy])
    px, py = probabilities(cx, R), probabilities(cy, R)
    H = entropy(px) + entropy(py)
    for g in gx[1:-1] + gy[1:-1]:
        H += entropy(probabilities([0., g], R))
    return logZ, H, px, py


def gaussian_control(M, theta, R=CRITICAL_R):
    ell = theta[-1]
    ci, cz = sum(theta[:2]), sum(theta[2:4])
    logB = sum(softplus(-((2 * a + 1) * ell + sign * cz))
               for a in range(M) for sign in (-1, 1))
    logB += sum(softplus(-((2 * j + 1 - 2 * M) * ell + ci)) for j in range(2 * M))
    gx, gy, _, _ = physical_energies(M, theta)
    logZ = sum(softplus(-g / R) for g in gx + gy)
    return {"log_profile_product": logB, "log_physical_partition": logZ,
            "critical_free_energy": logB - R * logZ}


def analyze(record):
    M, theta, logB = record['M'], record['parameters'], record['log_profile_product']
    energies = physical_energies(M, theta)
    Z, H, px, py = partition_and_entropy(energies, CRITICAL_R)
    F = logB - CRITICAL_R * Z
    given_R = 2 * record['rho']
    root_Z, _, _, _ = partition_and_entropy(energies, given_R)
    residual = logB - given_R * root_Z
    assert abs(residual) < 1e-8 * max(1., abs(logB))
    gaussian = gaussian_control(M, theta)
    return {"M": M, "source_rho": record['rho'], "source_parameters": theta,
            "critical_R": CRITICAL_R, "critical_free_energy": F,
            "critical_entropy": H, "linearized_rho": CRITICAL_R / 2 + F / (2 * H),
            "critical_type_x_core": px, "critical_type_y_core": py,
            "physical_core_energies_x": energies[2], "physical_core_energies_y": energies[3],
            "gaussian_control": gaussian,
            "critical_total_defect_excess": F - gaussian['critical_free_energy'],
            "saved_root_residual": residual}


def critical_gaussian_scale(M):
    """Golden-section search of the zero-charge Gaussian critical free energy."""
    def objective(ell):
        return gaussian_control(M, [0.] * 18 + [ell])['critical_free_energy']
    left, right = .05 / M, 8. / M
    ratio = (math.sqrt(5) - 1) / 2
    a, b = right - ratio * (right - left), left + ratio * (right - left)
    fa, fb = objective(a), objective(b)
    for _ in range(100):
        if fa < fb:
            right, b, fb = b, a, fa
            a = right - ratio * (right - left)
            fa = objective(a)
        else:
            left, a, fa = a, b, fb
            b = left + ratio * (right - left)
            fb = objective(b)
    return (left + right) / 2


def change_core_energies(theta, dx=None, dy=None):
    out = theta[:]
    stride = (len(theta) - 5) // 2
    if dx is not None:
        legs = (0, 7) if stride == 10 else (0,)
        for side in (0, 1):
            for leg in legs:
                for j, value in enumerate(dx):
                    out[4 + side * stride + leg + j] += value / (2 * len(legs))
    if dy is not None:
        for side in (0, 1):
            for leg in (3, 5):
                for j, value in enumerate(dy):
                    out[4 + side * stride + leg + j] += value / 4
    return out


def transported_starts(record, target_M):
    M, original = record['M'], record['parameters']
    source_energies = physical_energies(M, original)
    starts, metadata = [], []
    def save(name, theta):
        starts.append(theta)
        _, _, px, py = partition_and_entropy(physical_energies(target_M, theta), CRITICAL_R)
        metadata.append({"name": name, "log_s": theta[-1],
                         "largest_normalized_height": max(abs(x / theta[-1]) for x in theta[:-1]),
                         "critical_type_x_core": px, "critical_type_y_core": py})
    theta = original[:]
    theta[-1] *= (M - 1) / (target_M - 1)
    assert all(abs(a - b) < 1e-10 for old, new in zip(source_energies[2:], physical_energies(target_M, theta)[2:])
               for a, b in zip(old, new))
    save("preserve_extreme_energies", theta)
    theta = original[:]
    theta[-1] = critical_gaussian_scale(target_M)
    delta = (target_M - 1) * theta[-1] - (M - 1) * original[-1]
    theta = change_core_energies(theta, [2 * delta, -2 * delta, 0.], [2 * delta, -2 * delta])
    assert all(abs(a - b) < 1e-10 for old, new in zip(source_energies[2:], physical_energies(target_M, theta)[2:])
               for a, b in zip(old, new))
    save("critical_Gaussian_scale_preserve_core_energies", theta)
    def target(name, px=None, py=None):
        _, _, cx, cy = physical_energies(target_M, theta)
        dx = [-CRITICAL_R * math.log(p / px[0]) - e for p, e in zip(px[1:], cx[1:])] if px else None
        dy = [-CRITICAL_R * math.log(p / py[0]) - e for p, e in zip(py[1:], cy[1:])] if py else None
        candidate = change_core_energies(theta, dx, dy)
        _, _, actual_x, actual_y = partition_and_entropy(physical_energies(target_M, candidate), CRITICAL_R)
        if px:
            assert max(abs(a - b) for a, b in zip(px, actual_x)) < 1e-12
        if py:
            assert max(abs(a - b) for a, b in zip(py, actual_y)) < 1e-12
        save(name, candidate)
    target("activate_Y_vacuum", py=[.12, .879, .001])
    target("activate_Y_top", py=[.05, .475, .475])
    target("uniform_Y_core", py=[1/3] * 3)
    target("activate_both_core_top_states", px=[.05, .45, .45, .05], py=[.05, .475, .475])
    target("uniform_both_cores", px=[.25] * 4, py=[1/3] * 3)
    return starts, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/rectangular-critical-pressure-inputs.json'))
    parser.add_argument('--output', type=Path, default=Path('evidence/rectangular-critical-pressure.json'))
    parser.add_argument('--starts-directory', type=Path)
    args = parser.parse_args()
    records = json.loads(args.input.read_text())['records']
    analyses, transports = [], []
    for record in records:
        analyses.append(analyze(record))
        target_M = 2 * record['M']
        starts, metadata = transported_starts(record, target_M)
        transports.append({'source_M': record['M'], 'target_M': target_M,
                           'starts': starts, 'metadata': metadata})
        if args.starts_directory:
            args.starts_directory.mkdir(parents=True, exist_ok=True)
            name = args.starts_directory / f'defect-M{target_M}-transported-starts.json'
            name.write_text(json.dumps(starts, indent=2) + '\n')
    report = {'status': 'Floating discriminator and starts only; no exact pressure bound or new profile query.',
              'critical_R': CRITICAL_R, 'analyses': analyses, 'transports': transports,
              'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'analyses': analyses, 'transported_start_counts': [len(t['starts']) for t in transports]}, indent=2))


if __name__ == '__main__':
    main()
