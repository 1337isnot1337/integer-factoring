.PHONY: verify verify-docs verify-frontier verify-parity verify-recovery verify-constructive verify-native verify-exterior verify-entangler-density verify-post-half verify-formal verify-k6 verify-half-appendix help

verify: verify-docs verify-frontier verify-recovery verify-constructive
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

verify-parity:
	python3 tools/check_parity_power_divisor_recovery.py
	python3 tools/check_weighted_network_ports.py

verify-recovery: verify-parity
	python3 tools/check_fixed_characteristic_recovery.py
	python3 tools/check_root_threshold_counterexample.py
	python3 tools/check_unitary_orbit_recovery.py
	python3 tools/check_quadratic_congruence_recovery.py
	python3 tools/check_high_order_characteristic.py
	python3 tools/check_fixed_precision_count_lifting.py

verify-constructive:
	python3 tools/check_explicit_entanglers.py
	python3 tools/check_explicit_entanglers.py --control evidence/binary-entangler-extension-control.json
	python3 tools/check_smith_profile_transfer.py
	python3 tools/check_degree_sensitive_smith_paths.py
	python3 tools/check_trace_shear_general.py
	python3 tools/check_auxiliary_gauge_smith.py
	python3 tools/check_auxiliary_lattice_gauges.py
	python3 tools/check_balanced_physical_gauge.py
	python3 tools/check_sum_preserving_smith_barrier.py

verify-native: verify-exterior verify-entangler-density

verify-exterior:
	python3 tools/check_coherent_exterior_boundaries.py

verify-entangler-density:
	python3 tools/check_density16_certificate.py

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
	  'make verify-parity        - parity recovery and variable-output multiplication ports' \
	  'make verify-recovery      - factor recovery in every fixed prime characteristic' \
	  'make verify-constructive  - exact entangler, Smith-profile, and gauge certificates' \
	  'make verify-native        - optional coherent-boundary and complete density screens (C++)' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
