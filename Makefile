.PHONY: verify verify-docs verify-frontier verify-post-half verify-schur-fermion verify-schur-central verify-schur-gaussian verify-schur-sparse verify-fixed-modes verify-cancellation-modules verify-minimal-modules verify-exact-modules verify-cancellation-crossing verify-block-smith verify-linear-smith verify-type-pressure verify-crossing-interpolation verify-smith-oracle verify-controlled-smith verify-sparse-holonomy verify-charge-pruning verify-cubic-full-images verify-formal verify-k6 verify-half-appendix help

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

.PHONY: verify-mixed-place verify-network-codecs
verify-mixed-place:
	python3 tools/check_weighted_multipoint_hermite.py
	python3 tools/check_mixed_place_hermite.py
	python3 tools/check_ternary_basis_mixed.py

.PHONY: verify-adaptive-mixed-place
verify-adaptive-mixed-place:
	python3 tools/check_adaptive_mixed_place.py --certificate evidence/fused-q4-adaptive-certificate.json

verify-network-codecs:
	python3 tools/check_factoring_network_structure.py
	python3 tools/check_factoring_primitive_recipients.py
	python3 tools/check_factoring_recipient_rank_envelope.py

.PHONY: verify-adaptive-transform
verify-adaptive-transform:
	python3 tools/check_adaptive_transform_constructor.py

.PHONY: verify-operand-length
verify-operand-length:
	python3 tools/check_operand_length_transform.py

.PHONY: verify-fused-places
verify-fused-places:
	python3 tools/check_fused_multipoint_profiles.py
	python3 tools/check_fused_closed_profiles.py

.PHONY: verify-fused-weights
verify-fused-weights:
	python3 tools/check_fused_weighted_places.py

.PHONY: verify-compressed-fusion
verify-compressed-fusion:
	python3 tools/check_compressed_fusion.py --output /tmp/integer-factoring-compressed-fusion-replay.json
	python3 tools/check_f16_weighted_local_oracle.py --output /tmp/integer-factoring-f16-local-oracle-replay.json

.PHONY: verify-multilevel-packet
verify-multilevel-packet:
	python3 tools/check_multilevel_thue_morse_packet.py

.PHONY: verify-read-once-filters
verify-read-once-filters:
	python3 tools/check_read_once_boolean_filter_floor.py

.PHONY: verify-all-digit-gaussian
verify-all-digit-gaussian:
	python3 tools/check_all_digit_charged_gaussian_floor.py

.PHONY: verify-boolean-charge-bands
verify-boolean-charge-bands:
	python3 tools/check_homogeneous_boolean_charge_bands.py

.PHONY: verify-polynomial-flags
verify-polynomial-flags:
	python3 tools/check_binary_polynomial_flags.py

.PHONY: verify-unrestricted-gaussian
verify-unrestricted-gaussian:
	python3 tools/check_reflection_unrestricted_gaussian.py --output /tmp/integer-factoring-unrestricted-gaussian-replay.json

.PHONY: verify-global-homogeneous-filters
verify-global-homogeneous-filters:
	python3 tools/check_global_homogeneous_boolean_filter_floor.py --output /tmp/integer-factoring-global-homogeneous-replay.json

.PHONY: verify-multipoint-hermite
verify-multipoint-hermite:
	python3 tools/check_multipoint_hermite_crossings.py --output /tmp/integer-factoring-multipoint-hermite-replay.json

.PHONY: verify-gaussian-fixed-types
verify-gaussian-fixed-types:
	python3 tools/check_reflection_gaussian_fixed_type_duals.py --output /tmp/integer-factoring-gaussian-fixed-types-replay.json

.PHONY: verify-symmetric-binary-filters
verify-symmetric-binary-filters:
	python3 tools/check_symmetric_binary_filter_synthesis.py --output /tmp/integer-factoring-symmetric-binary-replay.json

.PHONY: verify-binary-smith-floor
verify-binary-smith-floor:
	python3 tools/check_binary_fixed_ray_floor.py --output /tmp/integer-factoring-binary-fixed-ray-replay.json
	python3 tools/check_binary_all_gauge_floor.py --output /tmp/integer-factoring-binary-all-gauge-replay.json

.PHONY: verify-symmetric-prime-filters verify-symmetric-filter-catalogues
verify-symmetric-prime-filters:
	python3 tools/check_symmetric_prime_filter_synthesis.py --output /tmp/integer-factoring-symmetric-prime-replay.json

verify-symmetric-filter-catalogues:
	python3 tools/check_symmetric_filter_catalogues.py --output /tmp/integer-factoring-symmetric-catalogue-replay.json

.PHONY: verify-gaussian-physical-frames verify-gaussian-coupled-gates
verify-gaussian-physical-frames:
	python3 tools/check_reflection_gaussian_physical_frames.py --output /tmp/integer-factoring-gaussian-physical-frames-replay.json

verify-gaussian-coupled-gates:
	python3 tools/check_reflection_gaussian_coupled_gates.py --output /tmp/integer-factoring-gaussian-coupled-gates-replay.json

.PHONY: verify-physical-splits verify-reflection-code-alphabets
verify-physical-splits:
	python3 tools/check_physical_split_crossings.py --output /tmp/integer-factoring-physical-split-replay.json

verify-reflection-code-alphabets:
	python3 tools/check_reflection_code_alphabets.py --output /tmp/integer-factoring-reflection-code-alphabets-replay.json

.PHONY: verify-nonresonant-filters verify-affine-filter-tiles
verify-nonresonant-filters:
	python3 tools/check_nonresonant_boolean_filter_floor.py --output /tmp/integer-factoring-nonresonant-filter-replay.json

verify-affine-filter-tiles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_affine_symmetric_tiles.py --output /tmp/integer-factoring-affine-filter-tile-replay.json

.PHONY: verify-frontier-extensions
verify-frontier-extensions:
	python3 tools/check_all_scale_balanced_obstruction.py
	python3 tools/check_truncated_symmetric_interval_frontier.py

.PHONY: verify-cauchy-factorial
verify-cauchy-factorial:
	python3 tools/check_arithmetic_cauchy_factorial.py

verify-schur-fermion:
	python3 tools/check_schur_fermion.py

verify-schur-central:
	python3 tools/check_schur_homogeneous_central.py

verify-schur-gaussian:
	python3 tools/check_schur_gaussian_reduction.py

verify-schur-sparse:
	python3 tools/check_schur_downset_crossing.py
	python3 tools/check_schur_offset_interval.py
	python3 tools/check_schur_two_block_holonomy.py

verify-fixed-modes:
	python3 tools/check_fixed_mode_crossings.py
	python3 tools/check_three_mode_pressure_certificate.py --certificate evidence/three-mode-pressure-refinement.json
	python3 tools/check_schur_one_species_interval.py
	python3 tools/check_large_interval_pressure.py
	python3 tools/check_schur_pascal_complement.py
	python3 tools/check_constant_flag_crossings.py

verify-cancellation-modules:
	python3 tools/check_cancellation_module_automata.py

verify-minimal-modules:
	python3 tools/check_minimal_boundary_modules.py

verify-exact-modules:
	python3 tools/check_exact_module_orbits.py

verify-cancellation-crossing:
	python3 tools/check_cancellation_module_crossing.py

verify-block-smith:
	python3 tools/check_block_smith_crossings.py

verify-linear-smith:
	python3 tools/check_linear_smith_barrier.py

verify-type-pressure:
	python3 tools/check_type_optimized_pressure.py

verify-crossing-interpolation:
	python3 tools/check_polynomial_crossing_interpolation.py
	python3 tools/check_tensorized_smith_sections.py
	python3 tools/check_positive_energy_quartic_endpoints.py
	python3 tools/check_localized_quartic_profiles.py
	python3 tools/check_residual_quartic_codes.py
	python3 tools/check_coded_quartic_mode_families.py
	python3 tools/check_conjugate_residual_phase_codes.py
	python3 tools/check_residual_three_letter_obstruction.py
	python3 tools/check_signed_residual_codes.py
	python3 tools/check_fixed_field_descent.py
	python3 tools/check_involutive_core_crossings.py
	python3 tools/check_coordinate_core_crossings.py
	python3 tools/check_even_core_postselection.py
	python3 tools/check_paired_type_crossing_pressure.py
	python3 tools/check_two_parameter_even_core.py
	python3 tools/check_rectangular_residual_repair.py
	python3 tools/check_rectangular_f4_repair.py
	python3 tools/check_four_parameter_even_core.py

verify-smith-oracle:
	python3 tools/check_weighted_smith_oracle.py
	python3 tools/check_sparse_bulk_species_floor.py

.PHONY: verify-native-smith
verify-native-smith:
	python3 tools/check_native_weighted_smith_f4.py

.PHONY: verify-native-rectangular
verify-native-rectangular:
	python3 tools/check_native_rectangular_smith.py
	python3 tools/check_native_rectangular_smith_f4.py

.PHONY: verify-frontier-reductions verify-native-rectangular-search
verify-frontier-reductions:
	python3 tools/check_rectangular_f4_repair.py
	python3 tools/check_rectangular_defect_support.py
	python3 tools/check_one_species_dyadic_affine.py
	python3 tools/check_pascal_complement_convolution.py
	python3 tools/check_prime_power_vacuum.py
	python3 tools/check_cauchy_packet_polynomials.py
	python3 tools/check_canonical_alternating_packet.py
	python3 tools/check_cauchy_packet_lattice.py

.PHONY: verify-bounded-defect
verify-bounded-defect:
	python3 tools/check_rectangular_defect_core_basis.py
	python3 tools/check_rectangular_defect_fock.py --compare-active
	python3 tools/check_rectangular_defect_coordinates.py
	python3 tools/check_rectangular_quartic_perturbations.py
	python3 tools/check_rectangular_core_pressure.py
	python3 tools/check_rectangular_refined_repairs.py

verify-native-rectangular-search:
	python3 tools/check_native_rectangular_smith_f4.py
	python3 tools/check_native_matrix_smith_f4.py
	python3 tools/check_matrix_f4_entry_shifts.py
	python3 tools/check_rectangular_f4_optimizer.py
	python3 tools/check_native_matrix_precision.py

.PHONY: verify-full-core verify-canonical-packets verify-all-field-crossings
verify-full-core:
	python3 tools/check_rectangular_full_core_support.py
	python3 tools/check_rectangular_defect_core_basis.py --full-core
	python3 tools/check_rectangular_full_core_coordinates.py
	python3 tools/check_rectangular_core_pressure.py --full-core
	python3 tools/check_rectangular_input_repairs.py

verify-canonical-packets:
	python3 tools/check_canonical_alternating_packet.py
	python3 tools/check_canonical_packet_exhaustive.py
	python3 tools/check_canonical_packet_decoder.py

verify-all-field-crossings:
	python3 tools/check_all_field_pascal_crossings.py
	python3 tools/check_graded_even_power.py
	python3 tools/check_truncated_symmetric_crossings.py

.PHONY: verify-core-transports verify-boolean-cauchy
verify-core-transports:
	python3 tools/check_rectangular_identity_transports.py
	python3 tools/check_rectangular_input_coordinates.py
	python3 tools/check_rectangular_coordinate_transforms.py
	python3 tools/check_rectangular_scalar_cache.py

verify-boolean-cauchy:
	python3 tools/check_boolean_cauchy_smith.py
	python3 tools/check_grouped_cauchy_inverse.py

.PHONY: verify-full-alphabet verify-core-polynomials verify-leading-core
verify-full-alphabet:
	python3 tools/check_full_alphabet_quartic_repair.py --full
	python3 tools/check_full_alphabet_quartic_profiles.py

verify-core-polynomials:
	python3 tools/check_rectangular_polynomial_shears.py --dense
	python3 tools/check_rectangular_extra_quartics.py
	python3 tools/check_rectangular_extra_coordinates.py

verify-leading-core:
	python3 tools/check_rectangular_low_background.py
	python3 tools/check_rectangular_leading_face_lp.py
	python3 tools/check_rectangular_leading_face_lp.py --input evidence/rectangular-full-leading-face-lp.json

verify-uniform-core:
	python3 tools/check_rectangular_uniform_entry_bounds.py
	python3 tools/check_rectangular_balanced_dual_gaps.py
	python3 tools/check_full_alphabet_leading_obstruction.py

verify-core-monomials:
	python3 tools/check_rectangular_core_monomials.py
	python3 tools/check_rectangular_monomial_coordinates.py
	python3 tools/check_rectangular_monomial_charges.py

verify-critical-cubics:
	python3 tools/check_full_alphabet_monomial_residual.py
	python3 tools/check_full_alphabet_monomial_profiles.py

verify-accelerated-interval:
	python3 tools/check_accelerated_interval_pressure.py

verify-walsh-permanents:
	python3 tools/check_walsh_cauchy_permanents.py
	python3 tools/check_walsh_twirl_group.py

verify-controlled-smith:
	python3 tools/check_controlled_permutation_smith.py
	python3 tools/check_qutrit_weyl_crossing.py

verify-sparse-holonomy:
	python3 tools/check_sparse_holonomy_singletons.py

verify-charge-pruning:
	python3 tools/check_one_sided_charge_pruning.py
	python3 tools/check_typed_gaussian_generic_blocks.py
	python3 tools/check_typed_gaussian_charge_blocks.py

verify-cubic-full-images:
	python3 tools/check_cubic_full_image_span.py --output /tmp/cubic-full-image-span.json

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
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-schur-fermion - exact full Schur entry and crossing queries' \
	  'make verify-schur-central - hook identities and exact Schur central ranks' \
	  'make verify-schur-gaussian - exact Schur digit reduction and Gaussian seed' \
	  'make verify-schur-sparse - exact down-set and offset interval reductions' \
	  'make verify-fixed-modes - arbitrary mode sets and independent occupations' \
	  'make verify-cancellation-modules - exact finite module and word certificates' \
	  'make verify-minimal-modules - exact boundary reduction and complete word rates' \
	  'make verify-exact-modules - complete rational lattice orbits and pressures' \
	  'make verify-cancellation-crossing - strict module improvement on a complete crossing' \
	  'make verify-block-smith - ordered block profiles and submultiplicativity' \
	  'make verify-linear-smith - uniform degree-one obstruction and binary catalogue' \
	  'make verify-type-pressure - rational type and physical-gauge certificates' \
	  'make verify-crossing-interpolation - exact polynomial synthesis from local flags' \
	  'make verify-smith-oracle - exact weighted Smith minima and sparse entropy dual' \
	  'make verify-native-smith - optional compiled F4 search controls (C++17)' \
	  'make verify-native-rectangular - optional rectangular paired-type controls (C++17)' \
	  'make verify-controlled-smith - controlled monomial and qutrit crossings' \
	  'make verify-sparse-holonomy - sparse Smith profiles and Gaussian identities' \
	  'make verify-charge-pruning - one-sided charges and typed Gaussian block bounds' \
	  'make verify-cubic-full-images - complete varying-image cubic spaces and ranks' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'

.PHONY: verify-weighted-pascal
verify-weighted-pascal:
	python3 tools/check_lucas_weighted_pascal.py

.PHONY: verify-frozen-core
verify-frozen-core:
	python3 tools/check_rectangular_frozen_core_rows.py
	python3 tools/check_frozen_core_row_profiles.py

.PHONY: verify-retained-profile-precision
verify-retained-profile-precision:
	python3 tools/check_retained_profile_precision.py

.PHONY: verify-cauchy-powers
verify-cauchy-powers:
	python3 tools/check_arithmetic_cauchy_powers.py

.PHONY: verify-strict-gaussian
verify-strict-gaussian:
	python3 tools/check_strict_gaussian_crossings.py

.PHONY: verify-slow-ramification
verify-slow-ramification:
	python3 tools/check_slow_ramification_crossings.py

.PHONY: verify-state-filters verify-state-filters-extended
verify-state-filters:
	python3 tools/check_truncated_symmetric_identity_filters.py
	python3 tools/check_prime_state_filter_profiles.py

verify-state-filters-extended:
	python3 tools/check_prime_state_filter_profiles.py --extended

.PHONY: verify-polynomial-state-filters
verify-polynomial-state-filters:
	python3 tools/check_polynomial_state_filter_synthesis.py

.PHONY: verify-coupled-filters
verify-coupled-filters:
	python3 tools/check_coupled_identity_mode_filters.py
	python3 tools/check_coupled_identity_extreme_filters.py
	python3 tools/check_binary_coupled_filter_profiles.py
	python3 tools/check_coupled_filter_pressure.py

.PHONY: verify-prime-filter-pressure
verify-prime-filter-pressure:
	python3 tools/check_prime_state_filter_pressure.py

.PHONY: verify-core-charge-projections
verify-core-charge-projections:
	python3 tools/check_strict_core_charge_profiles.py
	python3 tools/check_strict_core_charge_gauges.py
	python3 tools/check_raw_extreme_charge_projections.py

.PHONY: verify-frozen-core-complete
verify-frozen-core-complete:
	python3 tools/check_frozen_core_complete_search.py

.PHONY: verify-reflection-compression
verify-reflection-compression:
	python3 tools/check_reflection_dual_rail_compression.py

.PHONY: verify-coupled-synthesis
verify-coupled-synthesis:
	python3 tools/check_coupled_filter_linear_synthesis.py

.PHONY: verify-coupled-operators
verify-coupled-operators:
	python3 tools/check_coupled_identity_operator_filters.py --extended

.PHONY: verify-reflection-hot-cores
verify-reflection-hot-cores:
	python3 tools/check_reflection_hot_core_profiles.py --extended

.PHONY: verify-coupled-weyl
verify-coupled-weyl:
	python3 tools/check_coupled_identity_weyl_profiles.py --extended --screen evidence/coupled-identity-operator-screens.json

.PHONY: verify-repaired-state-filters
verify-repaired-state-filters:
	python3 tools/check_repaired_endpoint_state_filters.py
	python3 tools/check_repaired_endpoint_polynomial_synthesis.py
	python3 tools/check_f4_repaired_identity_filters.py
	python3 tools/check_f4_weyl_filter_compiler.py

.PHONY: verify-raw-core-masks
verify-raw-core-masks:
	python3 tools/check_raw_core_mask_profiles.py

.PHONY: verify-f4-filter-catalogue
verify-f4-filter-catalogue:
	python3 tools/check_f4_filter_catalogue.py

.PHONY: verify-reflection-hot-grid verify-reflection-hot-grid-exhaustive
verify-reflection-hot-grid:
	python3 tools/check_reflection_hot_core_grid.py
	python3 tools/reflection_hot_core_uniform.py

verify-reflection-hot-grid-exhaustive:
	python3 tools/check_reflection_hot_core_grid.py --exhaustive
	python3 tools/reflection_hot_core_uniform.py

.PHONY: verify-split-identity-modes verify-f4-identity-operators
verify-split-identity-modes:
	python3 tools/check_split_identity_mode_filters.py --extended

verify-f4-identity-operators:
	python3 tools/check_f4_identity_operator_filters.py
	python3 tools/check_compiled_f4_identity_operators.py

.PHONY: verify-weyl-lattice-synthesis
verify-weyl-lattice-synthesis:
	python3 tools/check_weyl_lattice_filter_synthesis.py
	python3 tools/check_weyl_lattice_endpoint_cover.py

.PHONY: verify-reflection-hot-state-gauges
verify-reflection-hot-state-gauges:
	python3 tools/check_reflection_hot_state_gauges.py --extended

.PHONY: verify-operator-lattice-synthesis
verify-operator-lattice-synthesis:
	python3 tools/check_weyl_operator_lattice_synthesis.py

.PHONY: verify-f4-operator-catalogues
verify-f4-operator-catalogues:
	python3 tools/check_f4_operator_catalogues.py

.PHONY: verify-field-lattice-lifts
verify-field-lattice-lifts:
	python3 tools/check_f4_lattice_linear.py

.PHONY: verify-f16-identity-operators
verify-f16-identity-operators:
	python3 tools/check_f16_identity_operator_profiles.py --extended

.PHONY: verify-reflection-compensated-gauges
verify-reflection-compensated-gauges:
	python3 tools/check_reflection_hot_compensated_gauges.py --extended

.PHONY: verify-prime-power-filters verify-prime-power-filters-extended
verify-prime-power-filters:
	python3 tools/check_prime_power_identity_filters.py
	python3 tools/check_joint_identity_filter_pressure.py
	python3 tools/check_selected_identity_operators.py --extended

verify-prime-power-filters-extended:
	python3 tools/check_prime_power_identity_filters.py --include-qutrit
	python3 tools/check_joint_identity_filter_pressure.py
	python3 tools/check_selected_identity_operators.py --extended

.PHONY: verify-operator-field-flags
verify-operator-field-flags:
	python3 tools/check_f4_operator_lattice_flags.py
	python3 tools/check_identity_operator_minimum_tradeoff.py

.PHONY: verify-resonant-operators
verify-resonant-operators:
	python3 tools/check_resonant_identity_operators.py
	python3 tools/check_resonant_laurent_catalogue.py

.PHONY: verify-reflection-leg-bulks
verify-reflection-leg-bulks:
	python3 tools/check_reflection_hot_leg_bulk_gauges.py --extended

.PHONY: verify-reflection-all-bulks
verify-reflection-all-bulks:
	python3 tools/check_reflection_hot_all_bulk_gauges.py --extended
	python3 tools/check_reflection_hot_state_gauges.py --extended

.PHONY: verify-near-power-operators
verify-near-power-operators:
	python3 tools/check_near_power_terminal_operators.py --extended

.PHONY: verify-joint-polynomial-synthesis verify-weyl-filter-flags
verify-joint-polynomial-synthesis:
	python3 tools/check_joint_binary_filter_synthesis.py
	python3 tools/check_joint_prime_filter_synthesis.py
	python3 tools/check_joint_filter_catalogues.py

verify-weyl-filter-flags:
	python3 tools/check_weyl_flag_filter_synthesis.py
	python3 tools/check_unit_reduced_lattice_basis.py

.PHONY: verify-binary-two-cut-lattices
verify-binary-two-cut-lattices:
	python3 tools/check_binary_two_cut_lattices.py

.PHONY: verify-deficit-identity-geometry
verify-deficit-identity-geometry:
	python3 tools/check_deficit_identity_mode_geometry.py

.PHONY: verify-above-power-operators
verify-above-power-operators:
	python3 tools/check_above_power_terminal_operators.py --extended

.PHONY: verify-homogeneous-boolean-filters
verify-homogeneous-boolean-filters:
	python3 tools/check_homogeneous_boolean_filter_bounds.py

.PHONY: verify-prime-deficit-operators
verify-prime-deficit-operators:
	python3 tools/check_prime_deficit_identity_operators.py --extended

.PHONY: verify-selected-deficit-geometry verify-surplus-identity-geometry
verify-selected-deficit-geometry:
	python3 tools/check_selected_one_deficit_geometry.py

verify-surplus-identity-geometry:
	python3 tools/check_surplus_identity_mode_geometry.py

.PHONY: verify-bulk-search-records
verify-bulk-search-records:
	python3 tools/check_reflection_hot_bulk_search_records.py

.PHONY: verify-charged-gaussian-floor verify-one-mode-laurent-filters
verify-charged-gaussian-floor:
	python3 tools/check_charged_gaussian_floor.py

verify-one-mode-laurent-filters:
	python3 tools/check_one_mode_binary_laurent_filters.py

.PHONY: verify-selected-deficit-operators verify-deficit-energy-conjugation
verify-selected-deficit-operators:
	python3 tools/check_selected_deficit_operators.py --extended

verify-deficit-energy-conjugation:
	python3 tools/check_deficit_energy_conjugation.py

.PHONY: verify-exterior-gaussian-floor verify-charged-one-mode-filters
verify-exterior-gaussian-floor:
	python3 tools/check_exterior_charged_gaussian_floor.py

verify-charged-one-mode-filters:
	python3 tools/check_charged_one_mode_filters.py

.PHONY: verify-folded-gaussian-floor verify-small-gaussian-bases
verify-folded-gaussian-floor:
	python3 tools/check_folded_charged_gaussian_floor.py

verify-small-gaussian-bases:
	python3 tools/check_small_charged_gaussian_bases.py

.PHONY: verify-global-gaussian-floor verify-global-one-mode-filters
verify-global-gaussian-floor: verify-charged-gaussian-floor verify-exterior-gaussian-floor verify-folded-gaussian-floor verify-small-gaussian-bases
	python3 tools/check_global_charged_gaussian_floor.py

verify-global-one-mode-filters: verify-global-gaussian-floor verify-one-mode-laurent-filters verify-charged-one-mode-filters

.PHONY: verify-prime-deficit-affine-synthesis
verify-prime-deficit-affine-synthesis:
	python3 tools/check_prime_deficit_operator_synthesis.py

.PHONY: verify-binary-target-box
verify-binary-target-box:
	python3 tools/check_binary_target_box.py --output /tmp/integer-factoring-binary-target-box-replay.json
.PHONY: verify-full-alphabet-affine
verify-full-alphabet-affine:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_full_alphabet_two_cut_lattices.py --output /tmp/integer-factoring-full-alphabet-affine-replay.json --search-records evidence/full-alphabet-q4-natural-affine-search.json evidence/full-alphabet-q4-degenerate-affine-search.json

.PHONY: verify-multiplier-automata verify-global-code-multipliers
verify-multiplier-automata:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_boolean_multiplier_automata.py --output /tmp/integer-factoring-multiplier-automata-replay.json

verify-global-code-multipliers:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_global_code_multipliers.py --output /tmp/integer-factoring-global-code-multipliers-replay.json

.PHONY: verify-quartic-tiles
verify-quartic-tiles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_quartic_tiles.py --output /tmp/integer-factoring-quartic-tiles-replay.json

.PHONY: verify-shared-minor-targets
verify-shared-minor-targets:
	python3 tools/check_shared_minor_profile_targets.py --output /tmp/integer-factoring-shared-minor-targets-replay.json

.PHONY: verify-regular-frame-cells verify-finite-smith-frames
verify-regular-frame-cells:
	python3 tools/check_formal_row_flags.py --output /tmp/integer-factoring-formal-row-flags-replay.json

verify-finite-smith-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_precision_smith_frames.py --output /tmp/integer-factoring-finite-smith-frames-replay.json

.PHONY: verify-odd-fiber-codes
verify-odd-fiber-codes:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_odd_fiber_crossings.py --output /tmp/integer-factoring-odd-fiber-codes-replay.json

.PHONY: verify-joint-syndrome-filters
verify-joint-syndrome-filters:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_joint_syndrome_filter_synthesis.py --output /tmp/integer-factoring-joint-syndrome-filters-replay.json

.PHONY: verify-typed-charge-lattices
verify-typed-charge-lattices:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_typed_charge_lattice_synthesis.py --output /tmp/integer-factoring-typed-charge-lattice-replay.json

.PHONY: verify-joint-multiplier-separation
verify-joint-multiplier-separation:
	python3 tools/check_joint_boolean_multiplier_separation.py --output /tmp/integer-factoring-joint-multiplier-separation-replay.json

.PHONY: verify-colored-row-fibers
verify-colored-row-fibers:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_colored_row_fibers.py --output /tmp/integer-factoring-colored-row-fibers-replay.json

.PHONY: verify-fused-six-vertex
verify-fused-six-vertex:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_fused_six_vertex_crossings.py --output /tmp/integer-factoring-fused-six-vertex-replay.json --search-records evidence/fused-six-vertex-gauges-m2-theta2.json evidence/fused-six-vertex-gauges-m3-theta4.json

.PHONY: verify-nonlinear-symbol-charges
verify-nonlinear-symbol-charges:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_nonlinear_symbol_charges.py --output /tmp/integer-factoring-nonlinear-symbol-charges-replay.json

.PHONY: verify-paired-quartic-poles
verify-paired-quartic-poles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_paired_quartic_words.py --output /tmp/integer-factoring-paired-quartic-poles-replay.json

.PHONY: verify-cube-charge-vertices
verify-cube-charge-vertices:
	python3 tools/check_cube_charge_vertices.py --output /tmp/integer-factoring-cube-charge-vertices-replay.json

.PHONY: verify-coefficient-height-walks
verify-coefficient-height-walks:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_full_alphabet_walk_records.py --records evidence/full-alphabet-alternating-walk-0.json evidence/full-alphabet-alternating-walk-1.json --output /tmp/integer-factoring-coefficient-height-walks-replay.json

.PHONY: verify-mixed-quartic-poles
verify-mixed-quartic-poles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_mixed_paired_words.py --output /tmp/integer-factoring-mixed-quartic-poles-replay.json

.PHONY: verify-binary-typed-endpoints
verify-binary-typed-endpoints:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_binary_typed_endpoint_gauges.py --output /tmp/integer-factoring-binary-typed-endpoints-replay.json

.PHONY: verify-joint-residue-retiming
verify-joint-residue-retiming:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_joint_syndrome_residue_retiming.py --output /tmp/integer-factoring-joint-residue-retiming-replay.json

.PHONY: verify-explicit-fusion verify-carry-fourier
.PHONY: verify-asymmetric-crossings
.PHONY: verify-weighted-adaptive verify-basis-searches
.PHONY: verify-frobenius-recovery
.PHONY: verify-physical-frames verify-ternary-majorization
.PHONY: verify-weighted-carries
.PHONY: verify-square-zero-frames verify-dense-coupled-frames
.PHONY: verify-constant-carry-weights
.PHONY: verify-uniform-reflection-jets verify-coupled-searches
.PHONY: verify-triangular-series
.PHONY: verify-paired-trace verify-companion-search verify-monoid-blocks
.PHONY: verify-unsplit-monoid
.PHONY: verify-degree16-crossing verify-reflection-truncations
.PHONY: verify-observable-constructors verify-prime-constructors
.PHONY: verify-frobenius-stars
.PHONY: verify-finite-crossing-certificates
.PHONY: verify-stable-hot-codes verify-hot64-certificate
verify-stable-hot-codes:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_implicit_reflection_codes.py --output /tmp/integer-factoring-implicit-hot.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_natural_cut_duality.py --output /tmp/integer-factoring-natural-duality.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_hot_code_charge_profiles.py --output /tmp/integer-factoring-hot-charge.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_weighted_hot_compounds.py --output /tmp/integer-factoring-weighted-compound.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_weighted_hot_code_oracle.py --output /tmp/integer-factoring-weighted-hot-oracle.json

verify-hot64-certificate:
	OPENBLAS_NUM_THREADS=1 python3 tools/certify_hot_charge_pressure.py --size 64 --modes 12 --base 459 460 --support 64 45 --output /tmp/integer-factoring-hot64-certificate.json

verify-finite-crossing-certificates:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_block_amplification.py --output /tmp/integer-factoring-finite-amplification.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_typed_adaptive_certificate.py --output /tmp/integer-factoring-finite-typed.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_rational_typed_adaptive_bounds.py --output /tmp/integer-factoring-rational-finite.json

verify-frobenius-stars:
	python3 tools/check_frobenius_rational_trace.py --output /tmp/integer-factoring-rational-star.json
	python3 tools/check_fermat_torus_binomials.py --output /tmp/integer-factoring-torus-star.json
	python3 tools/check_frobenius_trace_infinity.py --output /tmp/integer-factoring-star-infinity.json
	python3 tools/check_relative_frobenius_trace.py --output /tmp/integer-factoring-relative-trace.json
	python3 tools/check_semisimple_cyclic_box_contraction.py --output /tmp/integer-factoring-cyclic-box.json

verify-observable-constructors:
	python3 tools/check_observable_boundary_constructor.py --output /tmp/integer-factoring-observable-boundary.json
	python3 tools/check_lattice_boundary_decoder.py --output /tmp/integer-factoring-lattice-boundary.json

verify-prime-constructors:
	python3 tools/check_power_two_progression_dyadic_prime.py --output /tmp/integer-factoring-dyadic-prime.json
	python3 tools/check_prime_transform_boundary_constructor.py --output /tmp/integer-factoring-prime-transform.json
	python3 tools/check_width_sensitive_prime_constructor.py --output /tmp/integer-factoring-width-sensitive.json
	python3 tools/check_polynomial_entry_prime_constructor.py --output /tmp/integer-factoring-polynomial-entry.json
	python3 tools/check_rough_observable_factor_constructor.py --output /tmp/integer-factoring-rough-observable.json
	python3 tools/check_prime_root_field_construction.py --output /tmp/integer-factoring-prime-root-field.json

verify-degree16-crossing:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_degree16_reflection_crossing.py --output /tmp/integer-factoring-degree16-replay.json

verify-reflection-truncations:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_truncated_reflection_crossings.py --output /tmp/integer-factoring-truncation-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_truncated_reflection4_search.py --output /tmp/integer-factoring-truncation-pressure-replay.json

verify-unsplit-monoid:
	python3 tools/check_unsplit_monoid_contraction.py --output /tmp/integer-factoring-unsplit-monoid-replay.json

verify-paired-trace:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_paired_trace.py --output /tmp/integer-factoring-paired-trace-replay.json

verify-companion-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_coupled_companions.py --output /tmp/integer-factoring-companion-search-replay.json

verify-monoid-blocks:
	python3 tools/check_multiplicative_monoid_blocks.py --output /tmp/integer-factoring-monoid-blocks-replay.json

verify-triangular-series:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_triangular_reflection_jet_families.py --output /tmp/integer-factoring-triangular-jets-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_laurent_series_coupled_frames.py --output /tmp/integer-factoring-laurent-series-replay.json

verify-uniform-reflection-jets:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_low_jet_frames.py --include-interval16 --output /tmp/integer-factoring-uniform-reflection-jets-replay.json

verify-coupled-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection4_coupled_search.py --output /tmp/integer-factoring-coordinate-search-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_rank_one_pressure.py --output /tmp/integer-factoring-dense-pressure-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_affine_coupled_frames.py --output /tmp/integer-factoring-affine-frames-replay.json

verify-constant-carry-weights:
	python3 tools/check_constant_field_carry_weights.py --output /tmp/integer-factoring-constant-carry-weights-replay.json
	python3 tools/check_first_jet_carry_weights.py --output /tmp/integer-factoring-first-jet-carry-weights-replay.json

verify-square-zero-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_square_zero_coupled_frames.py --output /tmp/integer-factoring-square-zero-first-jet-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_square_zero_second_jet.py --output /tmp/integer-factoring-square-zero-second-jet-replay.json

verify-dense-coupled-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_coupled_reflection_frames.py --catalogue evidence/dense-rank-one-coupled-frames.json --output /tmp/integer-factoring-dense-coupled-frames-replay.json

verify-weighted-carries:
	python3 tools/check_weighted_dyadic_carry.py --output /tmp/integer-factoring-weighted-carries-replay.json

verify-physical-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_physical_polynomial_frames.py --output /tmp/integer-factoring-polynomial-frames-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_laurent_frames.py --output /tmp/integer-factoring-laurent-frames-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_coupled_frames.py --output /tmp/integer-factoring-coupled-frames-replay.json

verify-ternary-majorization:
	python3 tools/check_ternary_gauge_majorization.py --output /tmp/integer-factoring-ternary-majorization-replay.json

verify-frobenius-recovery:
	python3 tools/check_frobenius_sparse_recovery.py --output /tmp/integer-factoring-frobenius-recovery-replay.json

verify-weighted-adaptive:
	python3 tools/check_weighted_adaptive_catalogues.py --output /tmp/integer-factoring-weighted-adaptive-replay.json

verify-basis-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_ternary_basis_search.py --output /tmp/integer-factoring-ternary-basis-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_fused_weighted_winners.py evidence/fused-mixed-weighted-summary.json --output /tmp/integer-factoring-fused-winners-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_fused_physical_bases.py --output /tmp/integer-factoring-fused-bases-replay.json

verify-asymmetric-crossings:
	python3 tools/check_asymmetric_support_line.py --output /tmp/integer-factoring-asymmetric-support-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_asymmetric_pressure_oracle.py --output /tmp/integer-factoring-asymmetric-oracle-replay.json

verify-explicit-fusion:
	python3 tools/check_explicit_fused_row.py --output /tmp/integer-factoring-explicit-fused-row-replay.json
	python3 tools/check_fused_frobenius_spectators.py --output /tmp/integer-factoring-fused-spectators-replay.json

verify-carry-fourier:
	python3 tools/check_carry_fourier_elimination.py --output /tmp/integer-factoring-carry-fourier-replay.json

.PHONY: verify-paired-affine-space verify-bell-jets verify-finite-fibers verify-closed-place-oracles

verify-paired-affine-space:
	OPENBLAS_NUM_THREADS=1 python3 tools/saturate_paired_trace_affine_space.py evidence/independent-paired-affine-degree-basis.json.gz --output /tmp/integer-factoring-paired-affine-saturation.json --basis-output /tmp/integer-factoring-paired-affine-basis.json.gz

verify-bell-jets:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_entrywise_reflection_constraints.py --output /tmp/integer-factoring-entrywise-jets.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_bell_jet_crt_crossings.py --output /tmp/integer-factoring-bell-crt.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_bell_jet_transversality.py --output /tmp/integer-factoring-bell-transversality.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_ramified_carrier_minors.py --output /tmp/integer-factoring-ramified-carrier-minors.json

verify-finite-fibers:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_regular_smith_frames.py --output /tmp/integer-factoring-f4-regular-frames.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_jet_smith_fibers.py --output /tmp/integer-factoring-finite-jet-fibers.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_two_cut_smith_fibers.py --output /tmp/integer-factoring-f4-two-cut-fibers.json

verify-closed-place-oracles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_closed_place_oracle.py --output /tmp/integer-factoring-f4-closed-oracle.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_paired_mixed_place_pressure.py --output /tmp/integer-factoring-paired-mixed-pressure.json

.PHONY: verify-filtered-affine-search verify-laurent-search

verify-filtered-affine-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_filtered_reflection16_affine_search.py --output /tmp/integer-factoring-filtered-affine-search.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_filtered_reflection16_affine_combinations.py --output /tmp/integer-factoring-filtered-affine-combinations.json

verify-laurent-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_laurent_series_search.py --output /tmp/integer-factoring-laurent-series-search.json

.PHONY: verify-coherent-boundary

verify-coherent-boundary:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_coherent_hot_boundary.py --output /tmp/integer-factoring-coherent-hot8.json
	OPENBLAS_NUM_THREADS=1 python3 tools/checkpoint_coherent_hot_boundaries.py --output /tmp/integer-factoring-coherent-boundary-search.json

.PHONY: verify-root-worker-searches

verify-root-worker-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/checkpoint_root_worker_searches.py --output /tmp/integer-factoring-root-worker-searches.json

