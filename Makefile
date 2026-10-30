.PHONY: verify verify-docs verify-frontier verify-post-half verify-further verify-energy-full verify-search-native verify-witt-full verify-formal verify-k6 verify-half-appendix help

verify: verify-docs verify-frontier verify-further
	python3 -m unittest discover -s tests -v
	python3 -m unittest tools.check_crossing -v
	python3 tools/check_coefficient_basis.py
	python3 tools/parameters.py
	python3 tools/check_k6_sharp_rate.py
	python3 tools/check_half_rate.py
	python3 tools/check_half_piece_ranks.py --no-numpy
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-docs:
	python3 tools/check_docs.py

verify-further:
	python3 tools/check_vanishing_slack_transfer.py
	python3 tools/check_optimized_divisor_bound.py
	python3 tools/check_gap_free_seed_transfer.py
	python3 tools/check_seed_field_selection.py
	python3 tools/check_modular_binary_word_barrier.py
	python3 tools/check_all_field_binary_word_barrier.py
	python3 tools/check_hankel_digit_box.py
	python3 tools/check_padic_digit_carry.py
	python3 tools/check_small_space_berkowitz.py
	python3 tools/check_polylog_space_rank.py
	python3 tools/check_truncated_ring_factorization.py
	python3 tools/check_border_nonunit_transfer.py
	python3 tools/check_commuting_halves_matrix_seeds.py
	python3 tools/check_binary_energy_symmetry.py
	python3 tools/check_weighted_energy.py
	python3 tools/check_witt_interfaces.py
	python3 tools/check_weighted_weyl_interfaces.py
	python3 tools/check_weighted_weyl_matching.py
	python3 tools/check_weighted_weyl_filtration.py

verify-energy-full:
	python3 tools/check_weighted_energy.py --full

verify-search-native:
	python3 tools/check_markov_pairing_transfer.py
	python3 tools/check_square_zero_rank_reduction.py
	python3 tools/check_two_state_run_profile.py
	python3 tools/check_parity_matrix_crossings.py

verify-witt-full:
	python3 tools/witt_artin_hasse.py --h 32 --gram

verify-post-half:
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-frontier:
	python3 tools/check_piece_transfer.py
	python3 tools/check_boson_truncation.py
	python3 tools/check_graded_space.py
	python3 tools/check_rectangular_profile.py
	python3 tools/check_partition_channels.py --max-weight 6
	python3 tools/check_homogeneous_duals.py
	python3 tools/check_binary_prefix_channel.py
	python3 tools/check_near_ones_channels.py
	python3 tools/check_primary_central_gram.py
	python3 tools/check_modular_crossing.py
	python3 tools/check_sparse_divisor_nonzeroness.py

verify-formal:
	bash tools/check_formal.sh

verify-k6:
	bash tools/check_k6.sh

verify-half-appendix:
	bash tools/check_half_appendix.sh

help:
	@printf '%s\n' \
	  'make verify               - unit tests, exact certificates, and local documentation links' \
	  'make verify-frontier      - focused frontier certificates' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-further       - new transfers, arithmetic bounds, and word certificates' \
	  'make verify-energy-full   - larger optional weighted-energy certificate' \
	  'make verify-search-native - exact finite-state search engines (C++17 required)' \
	  'make verify-witt-full     - larger optional Witt Gram and vacuum calculation' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
