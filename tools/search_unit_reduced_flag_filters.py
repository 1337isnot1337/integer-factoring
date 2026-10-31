#!/usr/bin/env python3
"""Run the F4 flag search after exact zero-unit source normalization."""
from unit_reduced_lattice_basis import annotate_output, install


if __name__ == '__main__':
    install()
    from search_weyl_flag_filters import main
    main()
    annotate_output()
