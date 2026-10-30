.PHONY: verify verify-docs verify-frontier verify-code-width5 verify-post-half verify-formal verify-k6 verify-half-appendix help

verify: verify-docs verify-frontier
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

verify-post-half:
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-frontier:
	python3 tools/check_weighted_cut_rank_budget.py
	python3 tools/check_weighted_cut_residual_budget.py
	python3 tools/check_weighted_cut_dual_certificate.py
	python3 tools/check_two_product_complement.py
	python3 tools/check_two_product_seed_cutoff.py
	python3 tools/check_weighted_cut_relaxation_limit.py
	python3 tools/check_weighted_cut_border_stability.py
	python3 tools/check_quadratic_crossing.py
	python3 tools/check_word_restriction_barrier.py
	python3 tools/check_code_sections.py
	python3 tools/check_nonlinear_layer_probe.py
	python3 tools/check_growing_number_field_transfer.py
	python3 tools/check_characteristic_independent_recovery.py
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

verify-code-width5:
	g++ -O3 -std=c++20 -fopenmp tools/check_code_kernel_width5.cpp -o /tmp/integer-factoring-code-width5
	/tmp/integer-factoring-code-width5 0 64 4 /tmp/integer-factoring-code-width5
	python3 tools/check_code_sections.py --full-classifier --native-prefix /tmp/integer-factoring-code-width5

verify-k6:
	bash tools/check_k6.sh

verify-half-appendix:
	bash tools/check_half_appendix.sh

help:
	@printf '%s\n' \
	  'make verify               - unit tests, exact certificates, and local documentation links' \
	  'make verify-frontier      - focused frontier certificates' \
	  'make verify-code-width5   - regenerate every width-five coordinate-section class' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
