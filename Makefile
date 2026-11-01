.PHONY: verify verify-docs verify-frontier verify-arithmetic verify-arithmetic-framework verify-post-half verify-formal verify-k6 verify-half-appendix verify-bell help verify-schur-fermion verify-schur-central verify-schur-gaussian verify-schur-sparse verify-fixed-modes verify-cancellation-modules verify-minimal-modules verify-exact-modules verify-cancellation-crossing verify-block-smith verify-linear-smith verify-type-pressure verify-crossing-interpolation verify-smith-oracle verify-controlled-smith verify-sparse-holonomy verify-charge-pruning verify-cubic-full-images verify-arithmetic-sieve verify-polynomial-crossings verify-polynomial-native verify-nine-letter-all

verify: verify-docs verify-frontier verify-arithmetic
	python3 -m unittest discover -s tests -v
	python3 -m unittest tools.check_crossing -v
	python3 tools/check_coefficient_basis.py
	python3 tools/parameters.py
	python3 tools/check_k6_sharp_rate.py
	python3 tools/check_half_rate.py
	python3 tools/check_half_piece_ranks.py --no-numpy
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py
	python3 tools/check_weighted_prefix.py
	python3 tools/check_fixed_characteristic_transfer.py
	python3 tools/check_modular_pieces.py
	python3 tools/check_modular_pfaffian.py
	python3 tools/check_fermi_sea_schur.py
	python3 tools/check_sparse_grid_geometry.py

verify-docs:
	python3 tools/check_docs.py

.PHONY: verify-modular-descent
verify-modular-descent:
	python3 tools/check_modular_parameter_descent.py

.PHONY: verify-smooth-norms verify-neutral-duals verify-f4-recentered verify-odd-cyclic
verify-smooth-norms:
	python3 tools/check_smooth_norm_gap.py
	python3 tools/check_smooth_norm_relation.py
	python3 tools/check_factored_norm_search.py

.PHONY: verify-class-square-conic
verify-class-square-conic:
	python3 tools/check_class_square_conic.py

.PHONY: verify-divisor-moments
verify-divisor-moments:
	python3 tools/check_divisor_moment_recovery.py

.PHONY: verify-sparse-divisor-residues
verify-sparse-divisor-residues:
	python3 tools/check_sparse_divisor_residues.py

.PHONY: verify-odd-character-residues
verify-odd-character-residues:
	python3 tools/check_odd_character_residues.py

verify-neutral-duals:
	python3 tools/check_neutral_reflected_crossings.py
	python3 tools/check_homogeneous_dual_completion.py

verify-f4-recentered:
	python3 tools/check_f4_recentered_lattices.py

verify-odd-cyclic:
	python3 tools/check_odd_cyclic_crossing.py

verify-arithmetic:
	python3 tools/check_dyadic_hyperbola.py
	python3 tools/check_dyadic_orbit_aggregation.py
	python3 tools/check_interval_orbit_gcd.py
	python3 tools/check_interval_unit_jets.py
	python3 tools/check_interval_packet_circuits.py
	python3 tools/check_normalized_interval_jets.py
	python3 tools/check_interval_pole_trace.py
	python3 tools/check_pole_cauchy_slope.py

.PHONY: verify-interval-pole-trace
verify-interval-pole-trace:
	python3 tools/check_interval_pole_trace.py
	python3 tools/check_pole_cauchy_slope.py

verify-bell:
	python3 tools/check_bell_crossing.py
	python3 tools/check_bell_pencils.py
	python3 tools/check_bell_profiles.py evidence/bell-synthesis-d4h4.json evidence/bell-synthesis-d8h3.json evidence/bell-folded-adaptive-d8.json evidence/bell-folded-adaptive-d16.json
	python3 tools/check_folded_bell.py
	python3 tools/check_code_folded_bell.py
	python3 tools/check_bell_projection.py
	python3 tools/check_bell_cubic_floor.py
	python3 tools/check_bell_quartic_rigidity.py


verify-arithmetic-framework:
	python3 tools/check_arithmetic_two_group.py
	python3 tools/check_positive_relation_kernel.py
	python3 tools/check_arithmetic_bias_seed.py
	python3 tools/check_arithmetic_class_groups.py
	python3 tools/check_arithmetic_density_lift.py
	python3 tools/check_bounded_smooth_norm.py

.PHONY: verify-cyclic-algebra
verify-cyclic-algebra:
	python3 tools/check_cyclic_algebra_crossings.py

.PHONY: verify-shift-register
verify-shift-register:
	python3 tools/check_shift_register_crossings.py

.PHONY: verify-register-barrier verify-quadratic-register
verify-register-barrier:
	python3 tools/check_shift_register_copy_wire_barrier.py

verify-quadratic-register:
	python3 tools/check_quadratic_register_crossings.py

.PHONY: verify-quantum-pascal
verify-quantum-pascal:
	python3 tools/check_quantum_pascal_crossings.py

.PHONY: verify-matrix-vacuum
verify-matrix-vacuum:
	python3 tools/check_matrix_vacuum_barrier.py

.PHONY: verify-steiner-cubic
verify-steiner-cubic:
	python3 tools/check_steiner_cubic_crossings.py

.PHONY: verify-steiner-coefficients
verify-steiner-coefficients:
	python3 tools/check_steiner_coefficient_barrier.py

.PHONY: verify-ordered-steiner-gauges
verify-ordered-steiner-gauges:
	python3 tools/check_ordered_steiner_gauges.py

.PHONY: verify-weyl-ordered-boundaries
verify-weyl-ordered-boundaries:
	python3 tools/check_weyl_ordered_boundaries.py

.PHONY: verify-weyl-middle-coefficients
verify-weyl-middle-coefficients:
	python3 tools/check_weyl_middle_coefficients.py

.PHONY: verify-subgroup-return
verify-subgroup-return:
	python3 tools/check_subgroup_return_crossings.py

.PHONY: verify-group-shears
verify-group-shears:
	python3 tools/check_group_shear_crossings.py

.PHONY: verify-smith-capacity
verify-smith-capacity:
	python3 tools/check_smith_pressure_convexity.py

.PHONY: verify-multipoint-rank-one
verify-multipoint-rank-one:
	python3 tools/check_multipoint_rank_one_seeds.py

.PHONY: verify-multipoint-rank-two
verify-multipoint-rank-two:
	python3 tools/check_multipoint_rank_two_seeds.py

.PHONY: verify-multipoint-smith-floor
verify-multipoint-smith-floor:
	python3 tools/check_multipoint_smith_floor.py

.PHONY: verify-transpose-intertwiner
verify-transpose-intertwiner:
	python3 tools/check_transpose_intertwiner_barrier.py

.PHONY: verify-joint-determinant
verify-joint-determinant:
	python3 tools/check_joint_determinant_barrier.py

.PHONY: verify-bounded-spin
verify-bounded-spin:
	python3 tools/check_bounded_spin_crossing.py

.PHONY: verify-bounded-spin-algebra
verify-bounded-spin-algebra:
	python3 tools/check_bounded_spin_rank_algebra.py

.PHONY: verify-bounded-spin-six-vertex
verify-bounded-spin-six-vertex:
	python3 tools/check_bounded_spin_six_vertex.py

.PHONY: verify-downset-crossings verify-regular-code-crossings
verify-downset-crossings:
	python3 tools/check_downset_crossing_normal_form.py

verify-regular-code-crossings:
	python3 tools/search_regular_code_crossings.py --verify evidence/regular-code-crossings.json

.PHONY: verify-auxiliary-coefficient-span
verify-auxiliary-coefficient-span:
	python3 tools/search_auxiliary_coefficient_span.py --verify evidence/auxiliary-coefficient-span.json

.PHONY: verify-auxiliary-extension-span
verify-auxiliary-extension-span:
	python3 tools/search_auxiliary_extension_span.py --verify evidence/auxiliary-extension-span.json

.PHONY: verify-binary-site-symmetry verify-matching-curve-controls verify-matching-curve-jets
verify-binary-site-symmetry:
	python3 tools/check_binary_site_symmetry_barrier.py

verify-matching-curve-controls:
	python3 tools/check_degree_five_matching_candidate.py
	python3 tools/check_frobenius_cofactor_span.py

verify-matching-curve-jets:
	python3 tools/check_matching_curve_jets.py --prime 5 --output /tmp/integer-factoring-matching-five-jets.json
	python3 tools/check_matching_curve_jets.py --prime 11 --output /tmp/integer-factoring-matching-eleven-jets.json

.PHONY: verify-common-cofactor
verify-common-cofactor:
	python3 tools/check_common_cofactor_barrier.py

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

.PHONY: verify-polynomial-crossings verify-polynomial-native verify-nine-letter-all
verify-polynomial-crossings:
	python3 tools/check_bare_seven_to_twelve_letters.py
	python3 tools/check_quadratic_inverse_degree_caps.py
	python3 tools/check_quadratic_inverse_paths.py --alphabet 9
	python3 tools/check_leading_rank_entropy_floor.py
	python3 tools/check_quadratic_inverse_paths.py
	python3 tools/check_shared_minor_smith_floor.py
	python3 tools/check_all_degree_smith_floor.py
	python3 tools/check_quadratic_crossings.py
	python3 tools/check_extended_binary_crossings.py

verify-polynomial-native:
	python3 tools/check_polynomial_local_smith.py --nine-letter

verify-nine-letter-all:
	python3 tools/check_nine_letter_level_catalogue.py

verify-post-half:
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-arithmetic-sieve:
	python3 tools/check_arithmetic_sieve_certificate.py

verify-frontier: verify-bell verify-class-square-conic verify-divisor-moments verify-sparse-divisor-residues verify-odd-character-residues
	python3 tools/check_bare_seven_to_twelve_letters.py > /dev/null
	python3 tools/check_quadratic_inverse_degree_caps.py > /dev/null
	python3 tools/check_quadratic_inverse_paths.py --alphabet 9 > /dev/null
	python3 tools/check_leading_rank_entropy_floor.py > /dev/null
	python3 tools/check_quadratic_inverse_paths.py > /dev/null
	python3 tools/check_shared_minor_smith_floor.py > /dev/null
	python3 tools/check_all_degree_smith_floor.py > /dev/null
	python3 tools/check_quadratic_crossings.py > /dev/null
	python3 tools/check_arithmetic_two_group.py
	python3 tools/check_positive_relation_kernel.py
	python3 tools/check_arithmetic_bias_seed.py
	python3 tools/check_arithmetic_class_groups.py
	python3 tools/check_arithmetic_density_lift.py
	python3 tools/check_bounded_smooth_norm.py
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
  'make verify-polynomial-crossings - focused polynomial crossing certificates' \
  'make verify-polynomial-native - optional native local Smith profile checks' \
  'make verify-nine-letter-all - complete 31-path nine-letter profile catalogue' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-arithmetic-sieve - replay the 60/80-bit arithmetic certificates and finite search limits' \
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

.PHONY: verify-cell-minor-caps
verify-cell-minor-caps:
	python3 tools/check_cell_minor_caps.py --output /tmp/integer-factoring-cell-minor-caps-replay.json

.PHONY: verify-auxiliary-binary-floor
verify-auxiliary-binary-floor:
	python3 tools/check_auxiliary_binary_gauge_floor.py --output /tmp/integer-factoring-auxiliary-binary-floor-replay.json

.PHONY: verify-bare-seven-to-twelve-letters
verify-bare-seven-to-twelve-letters:
	python3 tools/check_bare_seven_to_twelve_letters.py

.PHONY: verify-modular
verify-modular:
	python3 tools/check_modular_pieces.py --full
	python3 tools/check_modular_pfaffian.py --full
	python3 tools/check_fermi_sea_schur.py
	python3 tools/check_weighted_prefix.py
	python3 tools/check_fixed_characteristic_transfer.py
	python3 tools/check_sparse_grid_geometry.py
