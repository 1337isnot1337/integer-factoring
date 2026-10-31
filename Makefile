.PHONY: verify verify-docs verify-frontier verify-constructions verify-native-crossings verify-post-half verify-formal verify-k6 verify-half-appendix help

verify: verify-docs verify-frontier verify-constructions
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

verify-constructions:
	python3 tools/check_interior_jacobian_budget.py --output /tmp/integer-factoring-interior-jacobian-budget.json
	python3 tools/check_rational_cut_mesh_budget.py --output /tmp/integer-factoring-rational-cut-mesh-budget.json
	python3 tools/check_compatible_compression_budget.py --cores evidence/compatible-compression-cores.json --output /tmp/integer-factoring-compatible-compression-budget.json
	python3 tools/check_carry_swap_circuits.py --witnesses evidence/carry-swap-finite-minors.json --output /tmp/integer-factoring-carry-swap-circuit-audit.json
	python3 tools/check_low_choi_cycle_searches.py --output /tmp/integer-factoring-low-choi-cycle-searches.json
	python3 tools/check_modified_transpose_roots.py --output /tmp/integer-factoring-modified-transpose-roots.json
	python3 tools/check_rank_two_channel_normal_form.py --output /tmp/integer-factoring-rank-two-channel-normal-form.json
	python3 tools/check_rank_one_fixed_pairs.py --output /tmp/integer-factoring-rank-one-fixed-pairs.json
	python3 tools/check_partial_transpose_rank.py --output /tmp/integer-factoring-partial-transpose-rank.json
	python3 tools/check_cnot_crossing_modules.py --output /tmp/integer-factoring-cnot-crossing-modules.json
	python3 tools/cnot_full_span.py --output /tmp/integer-factoring-cnot-full-span-audit.json
	python3 tools/check_cnot_cut_budgets.py --output /tmp/integer-factoring-cnot-cut-budgets.json
	python3 tools/check_geometric_transpose_crossings.py --polynomial-screen --output /tmp/integer-factoring-geometric-transpose-audit.json
	python3 tools/check_geometric_row_inverses.py --output /tmp/integer-factoring-geometric-row-inverse-audit.json
	python3 tools/check_nicolas_threshold_centers.py --output /tmp/integer-factoring-nicolas-threshold-centers.json
	python3 tools/check_adaptive_divisor_preprocessing.py --output /tmp/integer-factoring-adaptive-divisor-preprocessing.json
	python3 tools/check_certified_divisor_oracle.py --output /tmp/integer-factoring-certified-divisor-oracle.json
	python3 tools/check_large_divisor_checkpoint.py --source-dir evidence --output-dir /tmp/integer-factoring-large-divisor/evidence
	python3 tools/check_quantum_crossing.py --output /tmp/integer-factoring-quantum-crossing-audit.json
	python3 tools/check_nilpotent_flag_crossings.py --output /tmp/integer-factoring-nilpotent-flag-crossings.json
	python3 tools/check_cauchy_two_adic_smith.py --output /tmp/integer-factoring-cauchy-two-adic-smith.json
	python3 tools/check_tapered_hermite_precision.py --output /tmp/integer-factoring-tapered-hermite-precision.json
	python3 tools/check_exterior_precision_witness.py --output /tmp/integer-factoring-exterior-precision-witness.json
	python3 tools/check_arithmetic_crossing.py --output /tmp/integer-factoring-arithmetic-crossing-audit.json
	python3 tools/check_primorial_weighted_preprocessing.py
	python3 tools/check_exchange_divisor_preprocessing.py
	python3 tools/check_simultaneous_seed_endpoint.py
	python3 tools/check_homogeneous_central_symbols.py --central 16
	python3 tools/check_angle_pfaffian.py
	python3 tools/check_homogeneous_central_selection.py
	python3 tools/check_noncommutative_crossing_orbit.py
	python3 tools/check_bosonic_hermite_profile.py
	python3 tools/check_noncommutative_crossing_coded.py
	python3 tools/check_homogeneous_central_eighth_jet.py
	python3 tools/check_binary_quadratic_charge_shells.py
	python3 tools/check_noncommutative_crossing_hopf.py
	python3 tools/check_residue_folded_channels.py
	python3 tools/check_noncommutative_crossing_pfaffian_boundaries.py
	python3 tools/check_noncommutative_crossing_pfaffian_isotropic.py
	python3 tools/check_noncommutative_crossing_optimal_retracts.py
	python3 tools/check_noncommutative_crossing_finite_gates.py
	python3 tools/check_noncommutative_crossing_abelian_midpoints.py
	python3 tools/check_noncommutative_crossing_involution_support.py
	python3 tools/check_noncommutative_crossing_symplectic_radicals.py
	python3 tools/check_noncommutative_crossing_radical_ward.py
	python3 tools/check_noncommutative_crossing_ward_rank_descent.py
	python3 tools/check_noncommutative_crossing_ward_rectangles.py
	python3 tools/check_noncommutative_crossing_ward_euler.py
	python3 tools/check_pfaffian_charge_gauges.py
	python3 tools/check_ward_charge_filtration.py
	python3 tools/check_ward_dual_incidence.py
	python3 tools/check_ward_code_charge_entropy.py
	python3 tools/check_ward_rectangular_envelope.py
	python3 tools/check_ward_stable_operators.py --output /tmp/integer-factoring-ward-stable-operators.json
	python3 tools/check_ward_form_operator_transport.py --output /tmp/integer-factoring-ward-form-operator-transport.json
	python3 tools/check_ward_disjoint_transport.py --output /tmp/integer-factoring-ward-disjoint-transport.json
	python3 tools/check_ward_interval_transport.py --output /tmp/integer-factoring-ward-interval-transport.json
	python3 tools/check_ward_fan_family.py --output /tmp/integer-factoring-ward-fan-family.json
	python3 tools/check_ward_pfaffian_queries.py --output /tmp/integer-factoring-ward-pfaffian-queries.json
	python3 tools/check_ward_operator_types.py --output /tmp/integer-factoring-ward-operator-types-check.json
	python3 tools/check_overlapping_symplectic_profiles.py
	python3 tools/check_cubic_radical_output_maps.py
	python3 tools/check_cubic_cyclic_fourier_bound.py --output /tmp/integer-factoring-cubic-cyclic-fourier-bound.json
	python3 tools/check_cubic_jordan_filtration.py --output /tmp/integer-factoring-cubic-jordan-filtration.json
	python3 tools/check_noncommutative_crossing_seed_algebra_audit.py
	python3 tools/check_symplectic_seed_tensor_modules.py
	python3 tools/symplectic_seed_weighted_modules.py --output /tmp/integer-factoring-weighted-modules.json
	python3 tools/check_noncommutative_crossing_character_compression.py --output /tmp/integer-factoring-character-compression.json
	python3 tools/check_symplectic_seed_coefficient_barrier.py --output /tmp/integer-factoring-seed-coefficient-barrier.json
	python3 tools/check_double_endpoint_retracts.py --cases full:3:3
	python3 tools/check_noncommutative_crossing_group_degeneration.py

verify-native-crossings:
	mkdir -p build/native
	$(CXX) -O3 -std=c++17 -Wall -Wextra -Werror tools/explore_group_correction_rank.cpp -o build/native/group-correction-screen
	python3 tools/check_noncommutative_crossing_native.py --binary $(CURDIR)/build/native/group-correction-screen --lookup-summary evidence/noncommutative-crossing-lookup-complete.json
	python3 tools/check_homogeneous_central_profiles.py
	python3 tools/check_complement_shell_crossings.py
	$(CXX) -O3 -std=c++17 -Wall -Wextra -Werror tools/explore_graded_hopf_crossings.cpp -o build/native/graded-hopf-screen
	python3 tools/check_graded_hopf_crossings.py --binary $(CURDIR)/build/native/graded-hopf-screen
	$(CXX) -O3 -std=c++17 -Wall -Wextra -Werror tools/explore_graded_hopf_gauges.cpp -o build/native/graded-hopf-gauges
	python3 tools/check_graded_hopf_gauges.py --binary $(CURDIR)/build/native/graded-hopf-gauges
	$(CXX) -O3 -std=c++17 -Wall -Wextra -Werror tools/screen_pfaffian_span_f4.cpp -o build/native/pfaffian-span-screen
	python3 tools/export_pfaffian_span_native.py evidence/pfaffian-basis-all_P2_L2.json --output build/native/pfaffian-span-input.txt
	build/native/pfaffian-span-screen < build/native/pfaffian-span-input.txt > build/native/pfaffian-span-F4.jsonl
	python3 tools/check_pfaffian_span_f4.py --basis evidence/pfaffian-basis-all_P2_L2.json --screen build/native/pfaffian-span-F4.jsonl --output build/native/pfaffian-span-audit.json
	$(CXX) -O3 -std=c++20 -Wall -Wextra -Werror tools/screen_weighted_pfaffian_span.cpp -o build/native/weighted-pfaffian-screen
	python3 tools/check_uniform_f4_two_even_sections.py --native $(CURDIR)/build/native/weighted-pfaffian-screen --screen-output build/native/uniform-f4-section-lift.jsonl --output build/native/uniform-f4-section.json
	python3 tools/check_symplectic_seed_tensor_fields.py --types evidence/symplectic-seed-tensor-modules-certificate.json --power 2 --native $(CURDIR)/build/native/weighted-pfaffian-screen --compare-binary evidence/symplectic-seed-tensor-modules.json --export-dir build/native/seed-power-2 --output build/native/seed-power-2-audit.json
	python3 tools/check_symplectic_mixed_grade_codes.py --native $(CURDIR)/build/native/pfaffian-span-screen --field4 --export-dir build/native/mixed-grade --output build/native/mixed-grade-audit.json
	$(CXX) -O3 -std=c++20 -Wall -Wextra -Werror tools/screen_cubic_gf8.cpp -o build/native/cubic-gf8-screen
	python3 tools/explore_cubic_complement_gf8.py --replay evidence/cubic-complement-F8.json --output build/native/cubic-complement-F8-replay.json
	python3 tools/check_cubic_complement_gf8.py --input evidence/cubic-complement-F8.json --native $(CURDIR)/build/native/cubic-gf8-screen --output build/native/cubic-complement-F8-audit.json
	python3 tools/explore_cubic_mixed_gf8.py --replay evidence/cubic-mixed-F8.json --output build/native/cubic-mixed-F8-replay.json
	python3 tools/check_cubic_mixed_gf8.py --native $(CURDIR)/build/native/cubic-gf8-screen --output build/native/cubic-mixed-F8-audit.json

verify-post-half:
	python3 tools/check_common_character_compressions.py --output /tmp/integer-factoring-common-compressions.json
	python3 tools/check_two_even_site_sections.py --scope graphs --work-dir /tmp/integer-factoring-two-even-sections --output /tmp/integer-factoring-two-even-sections.json
	python3 tools/assemble_two_even_site_cores.py --pairs 6 --screen evidence/two-even-site-core-P6-F2.jsonl --output /tmp/integer-factoring-two-even-P6.json
	python3 tools/check_two_even_site_compressions.py --output /tmp/integer-factoring-two-even-P6-compressions.json
	python3 tools/check_two_even_site_jordan_compressions.py --output /tmp/integer-factoring-two-even-P7-compressions.json
	python3 tools/check_two_even_site_P8_compressions.py --output /tmp/integer-factoring-two-even-P8-compressions.json
	python3 tools/check_two_even_site_row_corners.py --output /tmp/integer-factoring-two-even-row-corners.json
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py
	python3 tools/check_complement_shell_ranks.py
	python3 tools/check_noncommutative_crossing_pbw.py
	python3 tools/check_folded_gram_channels.py
	python3 tools/check_folded_gram_aggregate.py
	python3 tools/check_pfaffian_p6_boundary_orbit.py
	python3 tools/check_unequal_grade_gram.py
	python3 tools/check_pfaffian_particle_hole.py
	python3 tools/check_stable_pfaffian_bounds.py
	python3 tools/check_stable_boundary_envelope.py --maximum-pairs 128
	python3 tools/check_stable_pfaffian_coordinates.py
	python3 tools/check_additional_crossing_screens.py
	python3 tools/explore_symplectic_boundary_codes.py --cases 4:4 --optimize
	python3 tools/explore_cubic_off_syndrome_outputs.py --output /tmp/integer-factoring-cubic-output-corrections.json
	python3 tools/check_cubic_output_identity.py --multirow evidence/cubic-multirow-outputs.json --output /tmp/integer-factoring-cubic-output-identity.json
	python3 tools/check_cubic_complement_outputs.py --binary evidence/cubic-complement-outputs.json --field4 evidence/cubic-complement-fields.json --output /tmp/integer-factoring-cubic-complement-audit.json
	python3 tools/check_cubic_affordable_pool.py --input evidence/cubic-multirow-outputs.json --output /tmp/integer-factoring-cubic-affordable-pool.json

verify-frontier:
	python3 tools/check_strict_mixed_cut_gate.py --output /tmp/integer-factoring-strict-mixed-cut-gate.json
	python3 tools/check_cnot_sixteen_span.py --output /tmp/integer-factoring-cnot-sixteen-span.json
	python3 tools/check_fixed_prime_factor_transfer.py --output /tmp/integer-factoring-fixed-prime-factor-transfer.json
	python3 tools/check_mixed_cut_rank_budget.py --output /tmp/integer-factoring-mixed-cut-rank-budget.json
	python3 tools/check_parity_factor_transfer.py --output /tmp/integer-factoring-parity-factor-transfer.json
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
	python3 tools/check_rough_sparse_transfer.py --output /tmp/integer-factoring-rough-sparse-transfer.json

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
	  'make verify-constructions - new crossing and resource-interface certificates' \
	  'make verify-native-crossings - compile and independently audit the native group screen' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
