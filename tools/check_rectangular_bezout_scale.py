#!/usr/bin/env python3
"""Independent integer controls for binary rectangular scales and support windows."""
import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from rectangular_bezout_scale import rectangular_bezout_scale, rectangular_newton_window


def valuation(value):
    return (value & -value).bit_length() - 1


def controls():
    maxima = windows = 0
    for rows in range(1, 65):
        for columns in range(1, 65):
            literal = [columns + k - valuation(comb(columns + k - 1, k))
                       for k in range(rows)]
            data = rectangular_bezout_scale(rows, columns)
            assert data["exponent"] == max(literal)
            assert literal[data["maximizing_index"]] == data["exponent"]
            assert data["exponent"] == rectangular_bezout_scale(columns, rows)["exponent"]
            assert data["examined_indices"] <= (rows + columns - 1).bit_length() + 1
            if rows == columns:
                assert data["exponent"] == 2 * rows - 1 - (rows - 1).bit_count()
            maxima += 1
            for precision in (1, 3, 8):
                result = rectangular_newton_window(rows, columns, precision)
                expected = [{"index": k, "valuation_gap": max(literal) - value}
                            for k, value in enumerate(literal) if max(literal) - value < precision]
                assert result["retained_newton_terms"] == expected
                assert result["examined_newton_indices"] <= data["defect"] + precision
                assert result["reverse_v_depth"] < data["defect"] + precision
                assert result["reverse_u_depth"] < data["defect"] + precision
                windows += 1
    nonmonotone = rectangular_newton_window(4, 6, 1)
    assert nonmonotone["retained_newton_terms"] == [{"index": 2, "valuation_gap": 0}]
    # The top coefficient vanishes at this precision, but the next one survives.
    assert nonmonotone["exponent"] - 6 - 3 + valuation(comb(8, 3)) == 2
    huge = []
    for rows, columns, precision in [
            ((1 << 200) + 13, (1 << 399) + 17, 9),
            ((1 << 1000) + 1, (1 << 2000) - 1, 4),
            (1 << 200, 1 << 200, 5),
            (3, 1 << 200, 7)]:
        result = rectangular_newton_window(rows, columns, precision)
        assert result["exponent"] == rectangular_bezout_scale(columns, rows)["exponent"]
        assert max(rows, columns) <= result["exponent"] <= rows + columns - 1
        huge.append({"rows": rows, "columns": columns, **result})
    invalid = 0
    for args in ((0, 2), (2, 0), (-1, 1), (2.5, 1)):
        try:
            rectangular_bezout_scale(*args)
        except ValueError:
            invalid += 1
        else:
            raise AssertionError("Invalid count accepted")
    return {"literal_binomial_maximum_cases": maxima, "complete_newton_windows": windows,
            "nonmonotone_gap_control": nonmonotone, "huge_binary_controls": huge,
            "invalid_count_rejections": invalid}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    result = {"status": "passed", "controls": controls(),
              "source_sha256": {name: sha256((root / 'tools' / name).read_bytes()).hexdigest()
                  for name in ('rectangular_bezout_scale.py', 'check_rectangular_bezout_scale.py')},
              "scope": "Binary computation of the proved rectangular scale and full terminal windows; no factor search."}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result['controls'].items() if isinstance(v, int)}, indent=2))


if __name__ == '__main__':
    main()
