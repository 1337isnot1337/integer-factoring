.PHONY: verify verify-docs verify-frontier verify-arithmetic verify-arithmetic-framework verify-post-half verify-formal verify-k6 verify-half-appendix verify-bell help verify-schur-fermion verify-schur-central verify-schur-gaussian verify-schur-sparse verify-fixed-modes verify-cancellation-modules verify-minimal-modules verify-exact-modules verify-cancellation-crossing verify-block-smith verify-linear-smith verify-type-pressure verify-crossing-interpolation verify-smith-oracle verify-controlled-smith verify-sparse-holonomy verify-charge-pruning verify-cubic-full-images verify-arithmetic-sieve verify-polynomial-crossings verify-polynomial-native verify-nine-letter-all

verify: verify-docs verify-frontier verify-arithmetic verify-further verify-constructions verify-recovery verify-constructive
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

.PHONY: verify-moving-boundary
verify-moving-boundary:
	python3 tools/check_moving_boundary_crossings.py --output /tmp/moving-boundary-controls.json
	python3 tools/check_moving_boundary_core_floor.py --output /tmp/moving-boundary-core-floor.json
	python3 tools/check_moving_interior_crossings.py --output /tmp/moving-interior-controls.json
	python3 tools/check_moving_physical_frames.py --output /tmp/moving-physical-frame-controls.json
	python3 tools/check_moving_five_core_floor.py --output /tmp/moving-five-core-floor.json
	python3 tools/check_moving_frame_prefix.py --output /tmp/moving-frame-height1-prefix16384-floor.json

.PHONY: verify-pole-smith verify-joint-top-boundary
verify-pole-smith:
	python3 tools/check_pole_smith_stabilization.py

verify-joint-top-boundary:
	python3 tools/check_joint_top_boundary_rank.py

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

.PHONY: verify-odd-character-residues verify-cyclotomic-decoder
verify-cyclotomic-decoder:
	python3 tools/check_cyclotomic_character_decoder.py

.PHONY: verify-lambert-theta
verify-lambert-theta:
	python3 tools/check_lambert_theta_sections.py

verify-frontier: verify-lambert-theta

.PHONY: verify-tate-coordinate
verify-tate-coordinate:
	python3 tools/check_tate_coordinate_compiler.py

verify-frontier: verify-tate-coordinate

.PHONY: verify-theta-programs
verify-theta-programs:
	python3 tools/check_theta_coordinate_program.py

verify-frontier: verify-theta-programs

.PHONY: verify-q7-parameter
verify-q7-parameter:
	python3 tools/check_q7_cubic_parameter.py

verify-frontier: verify-q7-parameter

.PHONY: verify-q5-quadratic
verify-q5-quadratic:
	python3 tools/check_q5_trace_zero.py

verify-frontier: verify-q5-quadratic

verify-odd-character-residues: verify-cyclotomic-decoder
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

.PHONY: verify-ordered-charge-area
verify-ordered-charge-area:
	python3 tools/check_ordered_charge_area.py

.PHONY: verify-joint-energy-conjugations
verify-joint-energy-conjugations:
	python3 tools/check_joint_energy_conjugations.py

.PHONY: verify-joint-leading-tangent
verify-joint-leading-tangent:
	python3 tools/check_joint_leading_tangent.py

.PHONY: verify-joint-energy-dressings
verify-joint-energy-dressings:
	python3 tools/check_joint_energy_dressings.py

.PHONY: verify-joint-energy-riccati
verify-joint-energy-riccati:
	python3 tools/check_joint_energy_riccati.py

.PHONY: verify-joint-polynomial-tangent
verify-joint-polynomial-tangent:
	python3 tools/check_joint_polynomial_tangent.py

.PHONY: verify-joint-monomial-lattices
verify-joint-monomial-lattices:
	python3 tools/check_joint_monomial_lattices.py

.PHONY: verify-joint-projector-lattices
verify-joint-projector-lattices:
	python3 tools/check_joint_projector_lattices.py

.PHONY: verify-joint-evaluation-algebra
verify-joint-evaluation-algebra:
	python3 tools/check_joint_evaluation_algebra.py

.PHONY: verify-joint-endpoint-loops
verify-joint-endpoint-loops:
	python3 tools/check_joint_endpoint_loops.py

.PHONY: verify-joint-nilpotent-depths
verify-joint-nilpotent-depths:
	python3 tools/check_joint_nilpotent_depths.py

.PHONY: verify-joint-filtered-gauges
verify-joint-filtered-gauges:
	python3 tools/check_joint_filtered_gauges.py

.PHONY: verify-fused-scalar-projection
verify-fused-scalar-projection:
	python3 tools/check_fused_scalar_projection.py

.PHONY: verify-fused-ordered-projection verify-fused-ordered-projection-full
verify-fused-ordered-projection:
	python3 tools/check_fused_ordered_projection.py

verify-fused-ordered-projection-full:
	python3 tools/check_fused_ordered_projection.py --full-replay --jobs 4

.PHONY: verify-fused-primary-code
verify-fused-primary-code:
	python3 tools/check_fused_primary_code.py

.PHONY: verify-joint-graded-newton
verify-joint-graded-newton:
	python3 tools/check_joint_graded_newton.py

.PHONY: verify-joint-zero-prefix
verify-joint-zero-prefix:
	python3 tools/check_joint_zero_prefix.py

.PHONY: verify-joint-two-place-obstruction
verify-joint-two-place-obstruction:
	python3 tools/check_joint_two_place_obstruction.py

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

.PHONY: verify-bounded-spin-second
verify-bounded-spin-second:
	python3 tools/check_bounded_spin_second.py

.PHONY: verify-bounded-spin-algebra
verify-bounded-spin-algebra:
	python3 tools/check_bounded_spin_rank_algebra.py

.PHONY: verify-bounded-spin-six-vertex
verify-bounded-spin-six-vertex:
	python3 tools/check_bounded_spin_six_vertex.py

.PHONY: verify-bounded-spin-projector
verify-bounded-spin-projector:
	python3 tools/check_bounded_spin_projector.py

.PHONY: verify-bounded-spin-support
verify-bounded-spin-support:
	python3 tools/check_bounded_spin_support.py

.PHONY: verify-bounded-spin-growth
verify-bounded-spin-growth:
	python3 tools/check_bounded_spin_growth.py

.PHONY: verify-paired-paths
verify-paired-paths:
	python3 tools/check_paired_path_crossing.py

.PHONY: verify-paired-height
verify-paired-height:
	python3 tools/check_paired_path_height.py

.PHONY: verify-quantum-braid verify-quantum-braid-profiles
verify-quantum-braid:
	python3 tools/check_quantum_block_braid.py

verify-quantum-braid-profiles:
	python3 tools/check_quantum_block_braid_profile.py

.PHONY: verify-braid-path verify-braid-path-profiles
verify-braid-path:
	python3 tools/check_quantum_braid_path.py

verify-braid-path-profiles:
	python3 tools/check_quantum_braid_path_profiles.py

.PHONY: verify-braid-reversal
verify-braid-reversal:
	python3 tools/check_quantum_braid_reversal.py

.PHONY: verify-braid-matrix-aux
verify-braid-matrix-aux:
	python3 tools/check_quantum_braid_matrix_aux.py

.PHONY: verify-braid-rank-one
verify-braid-rank-one:
	python3 tools/check_quantum_braid_rank_one.py

.PHONY: verify-braid-nilpotent-seam
verify-braid-nilpotent-seam:
	python3 tools/check_quantum_braid_nilpotent_seam.py

.PHONY: verify-braid-rank-envelope
verify-braid-rank-envelope:
	python3 tools/check_nilpotent_seam_rank_bounds.py

.PHONY: verify-braid-fusion-seam
verify-braid-fusion-seam:
	python3 tools/check_quantum_braid_fusion_seam.py

.PHONY: verify-braid-graded-seam
verify-braid-graded-seam:
	python3 tools/check_quantum_braid_graded_seam.py

.PHONY: verify-braid-spectral-fusion
verify-braid-spectral-fusion:
	python3 tools/check_quantum_braid_spectral_fusion.py

.PHONY: verify-braid-fusion-occupancy
verify-braid-fusion-occupancy:
	python3 tools/check_quantum_braid_fusion_occupancy.py

.PHONY: verify-braid-leading-automaton
verify-braid-leading-automaton:
	python3 tools/check_nilpotent_seam_top_rank.py

.PHONY: verify-catalan-cup-crossing
verify-catalan-cup-crossing:
	python3 tools/check_catalan_cup_crossing.py

.PHONY: verify-catalan-cup-jets
verify-catalan-cup-jets:
	python3 tools/check_catalan_cup_jets.py

.PHONY: verify-catalan-gram-smith
verify-catalan-gram-smith:
	python3 tools/check_catalan_gram_smith.py --output evidence/catalan-gram-dual-smith.json

.PHONY: verify-catalan-spherical-ranks
verify-catalan-spherical-ranks:
	python3 tools/check_catalan_spherical_ranks.py --output evidence/catalan-spherical-ranks.json

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
	python3 tools/check_native_polynomial_certificates.py
	python3 tools/check_polynomial_local_smith.py --nine-letter

verify-nine-letter-all:
	python3 tools/check_nine_letter_level_catalogue.py

verify-post-half: verify-pr29-extended
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-arithmetic-sieve:
	python3 tools/check_arithmetic_sieve_certificate.py

verify-frontier: verify-bell verify-class-square-conic verify-divisor-moments verify-sparse-divisor-residues verify-odd-character-residues verify-pr31-frontier verify-f4-hyperoval verify-hot-boundary-floor
	python3 tools/check_direct_sublinear_counts.py
	python3 tools/check_rational_simultaneous_transfer.py
	python3 tools/check_minimum_cut_rank_budget.py
	python3 tools/check_two_piece_rank_sum.py
	python3 tools/check_permutation_coordinate_restrictions.py
	python3 tools/check_swap_product_rank_sum.py
	python3 tools/check_conserved_convolution.py
	python3 tools/check_ideal_hankel_endpoint.py
	python3 tools/check_seed_subspace_recombination.py
	python3 tools/check_permutation_span_cycles.py
	python3 tools/check_distinguished_rank_budget.py
	python3 tools/check_simultaneous_time_space_transfer.py
	python3 tools/check_fast_matrix_clique.py
	python3 tools/check_sublinear_space_recovery.py
	python3 tools/check_streamed_top_space.py
	python3 tools/check_blocked_half_family.py
	python3 tools/check_near_linear_factor_space.py
	python3 tools/check_linear_space_recovery.py
	python3 tools/check_streamed_seed_transfer.py
	python3 tools/check_six_label_rank_sum.py
	python3 tools/check_cyclic_orbit_crossing.py
	python3 tools/check_nilpotent_cycle_labels.py
	python3 tools/check_split_permutation_seed.py
	python3 tools/check_balanced_partition_crossings.py
	python3 tools/check_schur_fusion_rank.py
	python3 tools/check_regular_fusion_rank.py
	python3 tools/check_permutation_homogeneous_labels.py
	python3 tools/check_intersection_matching_lift.py
	python3 tools/check_reflected_support_reduction.py
	python3 tools/check_paired_mode_johnson.py
	python3 tools/check_sos_partial_transpose.py
	python3 tools/check_finite_seed_amplification.py
	python3 tools/check_fixed_defect_invariants.py
	python3 tools/check_balanced_three_label_barrier.py
	python3 tools/check_one_defect_seed.py
	python3 tools/check_isotropic_three_label_barrier.py
	python3 tools/check_integral_counting_transfer.py
	python3 tools/check_partial_matching_minors.py
	python3 tools/check_prime_power_gram_involution.py
	python3 tools/check_schur_fusion_gram.py
	python3 tools/check_dual_gauge_filters.py
	python3 tools/check_linear_carry_blocks.py
	python3 tools/check_matching_code_lift.py
	python3 tools/check_matching_gram_unipotence.py
	python3 tools/check_logarithmic_wire_network.py
	python3 tools/check_fourier_odd_characteristic.py
	python3 tools/check_polynomial_top_queries.py
	python3 tools/check_bit_tagged_witness_recovery.py
	python3 tools/check_modular_factor_space_odd.py
	python3 tools/check_binary_coefficient_entropy.py
	python3 tools/check_symmetric_three_matching.py
	python3 tools/check_fourier_multiplication_gates.py
	python3 tools/check_integral_schur_crossing.py
	python3 tools/check_prime_power_capacities.py
	python3 tools/check_rough_input_transfer.py
	python3 tools/check_modular_factor_space.py
	python3 tools/check_fermionic_piece_polynomial.py
	python3 tools/check_finite_order_crossing_transfer.py
	python3 tools/check_coordinate_support_rank.py
	python3 tools/check_rational_height_transfer.py
	python3 tools/check_subexponential_piece_transfer.py
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
	python3 tools/check_weighted_cut_rank_budget.py
	python3 tools/check_weighted_cut_residual_budget.py
	python3 tools/check_weighted_cut_dual_certificate.py
	python3 tools/check_two_product_complement.py
	python3 tools/check_two_product_seed_cutoff.py
	python3 tools/check_flagged_cut_rank_budget.py
	python3 tools/check_geometric_subspace_compression.py
	python3 tools/check_geometric_cut_rank_budget.py
	python3 tools/check_weighted_cut_relaxation_limit.py
	python3 tools/check_geometric_cut_relaxation_limit.py
	python3 tools/check_weighted_cut_border_stability.py
	python3 tools/check_geometric_border_cutoff.py
	python3 tools/check_rank_one_formal_cutoff.py
	python3 tools/check_quadratic_crossing.py
	python3 tools/check_word_restriction_barrier.py
	python3 tools/check_code_sections.py
	python3 tools/check_frame_filtration.py
	python3 tools/check_feistel_filtration.py
	python3 tools/check_diagonal_word_rank_sum.py
	python3 tools/check_incidence_kernel_rank_sum.py
	python3 tools/check_full_conservation_rank_sum.py
	python3 tools/check_p_group_convolution_degree.py
	python3 tools/check_bi_exterior_crossing_algebra.py
	python3 tools/check_bi_exterior_degree_budget.py
	python3 tools/check_nonlinear_layer_probe.py
	python3 tools/check_pauli_walsh_dictionary.py
	python3 tools/check_pauli_walsh_pairs.py
	python3 tools/check_growing_number_field_transfer.py
	python3 tools/check_rank_sum_global_transfer.py
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
	python3 tools/check_modular_recombination.py
	python3 tools/check_modular_particle_hole.py
	python3 tools/check_hamming_pfaffian.py
	python3 tools/check_hamming_tensor_power.py
	python3 tools/check_sparse_width_interpolation.py
	python3 tools/check_hamming_sectors.py
	python3 tools/check_hamming_complement.py
	python3 tools/check_hamming_joint.py
	python3 tools/check_hamming_joint_lower.py
	python3 tools/check_extension_descent.py
	python3 tools/check_hamming_f4_probe.py
	python3 tools/check_frozen_core_leak.py
	python3 tools/check_word_symmetry_barrier.py
	python3 tools/check_pairing_word_barrier.py
	python3 tools/check_pairing_small_seeds.py
	python3 tools/check_affine_code_barrier.py
	python3 tools/check_hamming_relabeling.py
	python3 tools/check_hamming_relabeling_optimum.py
	python3 tools/check_sl2_boolean_rank.py
	python3 tools/check_all_field_pairing.py
	python3 tools/check_modular_central_rank.py
	python3 tools/check_modular_binary_central_rank.py
	python3 tools/check_three_label_span.py
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

.PHONY: verify-pr23
verify-pr23:
	python3 tools/check_direct_sublinear_counts.py
	python3 tools/check_rational_simultaneous_transfer.py
	python3 tools/check_minimum_cut_rank_budget.py
	python3 tools/check_two_piece_rank_sum.py
	python3 tools/check_permutation_coordinate_restrictions.py
	python3 tools/check_swap_product_rank_sum.py
	python3 tools/check_conserved_convolution.py
	python3 tools/check_ideal_hankel_endpoint.py
	python3 tools/check_seed_subspace_recombination.py
	python3 tools/check_permutation_span_cycles.py
	python3 tools/check_distinguished_rank_budget.py
	python3 tools/check_simultaneous_time_space_transfer.py
	python3 tools/check_fast_matrix_clique.py
	python3 tools/check_sublinear_space_recovery.py
	python3 tools/check_streamed_top_space.py
	python3 tools/check_blocked_half_family.py
	python3 tools/check_near_linear_factor_space.py
	python3 tools/check_linear_space_recovery.py
	python3 tools/check_streamed_seed_transfer.py
	python3 tools/check_six_label_rank_sum.py
	python3 tools/check_cyclic_orbit_crossing.py
	python3 tools/check_nilpotent_cycle_labels.py
	python3 tools/check_split_permutation_seed.py
	python3 tools/check_balanced_partition_crossings.py
	python3 tools/check_schur_fusion_rank.py
	python3 tools/check_regular_fusion_rank.py
	python3 tools/check_permutation_homogeneous_labels.py
	python3 tools/check_intersection_matching_lift.py
	python3 tools/check_reflected_support_reduction.py
	python3 tools/check_paired_mode_johnson.py
	python3 tools/check_sos_partial_transpose.py
	python3 tools/check_finite_seed_amplification.py
	python3 tools/check_fixed_defect_invariants.py
	python3 tools/check_balanced_three_label_barrier.py
	python3 tools/check_one_defect_seed.py
	python3 tools/check_isotropic_three_label_barrier.py
	python3 tools/check_integral_counting_transfer.py
	python3 tools/check_partial_matching_minors.py
	python3 tools/check_prime_power_gram_involution.py
	python3 tools/check_schur_fusion_gram.py
	python3 tools/check_dual_gauge_filters.py
	python3 tools/check_linear_carry_blocks.py
	python3 tools/check_matching_code_lift.py
	python3 tools/check_matching_gram_unipotence.py
	python3 tools/check_logarithmic_wire_network.py
	python3 tools/check_fourier_odd_characteristic.py
	python3 tools/check_polynomial_top_queries.py
	python3 tools/check_bit_tagged_witness_recovery.py
	python3 tools/check_modular_factor_space_odd.py
	python3 tools/check_binary_coefficient_entropy.py
	python3 tools/check_symmetric_three_matching.py
	python3 tools/check_fourier_multiplication_gates.py
	python3 tools/check_integral_schur_crossing.py
	python3 tools/check_prime_power_capacities.py
	python3 tools/check_rough_input_transfer.py
	python3 tools/check_modular_factor_space.py
	python3 tools/check_fermionic_piece_polynomial.py
	python3 tools/check_finite_order_crossing_transfer.py
	python3 tools/check_coordinate_support_rank.py
	python3 tools/check_rational_height_transfer.py
	python3 tools/check_subexponential_piece_transfer.py

.PHONY: verify-pr24
verify-pr24:
	python3 tools/check_modular_recombination.py
	python3 tools/check_modular_particle_hole.py
	python3 tools/check_hamming_pfaffian.py
	python3 tools/check_hamming_tensor_power.py
	python3 tools/check_sparse_width_interpolation.py
	python3 tools/check_hamming_sectors.py
	python3 tools/check_hamming_complement.py
	python3 tools/check_hamming_joint.py
	python3 tools/check_hamming_joint_lower.py
	python3 tools/check_extension_descent.py
	python3 tools/check_hamming_f4_probe.py
	python3 tools/check_frozen_core_leak.py
	python3 tools/check_word_symmetry_barrier.py
	python3 tools/check_pairing_word_barrier.py
	python3 tools/check_pairing_small_seeds.py
	python3 tools/check_affine_code_barrier.py
	python3 tools/check_hamming_relabeling.py
	python3 tools/check_hamming_relabeling_optimum.py
	python3 tools/check_sl2_boolean_rank.py
	python3 tools/check_all_field_pairing.py
	python3 tools/check_modular_central_rank.py
	python3 tools/check_modular_binary_central_rank.py
	python3 tools/check_three_label_span.py

.PHONY: verify-further verify-energy-full verify-search-native verify-weyl-orbit verify-witt-full
verify-further:
	python3 tools/check_vanishing_slack_transfer.py
	python3 tools/check_optimized_divisor_bound.py
	python3 tools/check_exact_divisor_maximum.py
	python3 tools/check_half_second_order.py
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
	python3 tools/check_boundary_stable_catalytic.py
	python3 tools/check_reduction_map_fierz.py
	python3 tools/check_binary_two_piece_f4_boundary.py
	python3 tools/check_three_piece_binary_boundary.py
	python3 tools/check_four_piece_f4_boundary.py
	python3 tools/check_observable_transition_realization.py
	python3 tools/check_observable_rank_detector.py
	python3 tools/check_local_alphabet_spectral_barriers.py
	python3 tools/check_finite_functional_detectors.py
	python3 tools/check_overlapping_subalphabet_detectors.py
	python3 tools/check_twisted_pairing_catalogue.py
	python3 tools/check_twisted_witnesses_dense.py
	python3 tools/check_mismatched_commutator.py
	python3 tools/check_twisted_channel_walks.py
	python3 tools/check_channel_flow_duals.py
	python3 tools/check_component_channel_walks.py
	python3 tools/check_exact_channel_growth.py
	python3 tools/check_channel_support_feasibility.py
	python3 tools/check_scalarized_channel_automata.py
	python3 tools/check_generic_factor_nonscalarization.py
	python3 tools/check_binary_energy_symmetry.py
	python3 tools/check_weighted_energy.py
	python3 tools/check_witt_interfaces.py
	python3 tools/check_weighted_weyl_interfaces.py
	python3 tools/check_weighted_weyl_matching.py
	python3 tools/check_weighted_weyl_filtration.py
	python3 tools/check_weighted_weyl_additive.py
	python3 tools/check_weyl_alternate_algebra.py
	python3 tools/check_eight_feature_monoid.py
	python3 tools/check_trace_bilinear_kernel.py
	python3 tools/check_two_patch_trace.py --output /tmp/integer-factoring-two-patch-trace.json
	python3 tools/check_two_patch_span.py --output /tmp/integer-factoring-two-patch-span.json
	python3 tools/check_tagged_two_patch_span.py --time-limit 300 --output /tmp/integer-factoring-tagged-two-patch-span.json
	python3 tools/check_bivariate_two_patch_span.py --time-limit 300 --output /tmp/integer-factoring-bivariate-two-patch-span.json
	python3 tools/check_polynomial_shear_probe.py --output /tmp/integer-factoring-polynomial-shear-verification.json

verify-energy-full:
	python3 tools/check_weighted_energy.py --full

verify-search-native:
	python3 tools/check_markov_pairing_transfer.py
	python3 tools/check_square_zero_rank_reduction.py
	python3 tools/check_two_state_run_profile.py
	python3 tools/check_parity_matrix_crossings.py
	python3 tools/check_site_weight_crossings.py
	python3 tools/check_eight_feature_crossings.py
	python3 tools/check_matrix_unit_monodromy.py
	python3 tools/check_exterior_monodromy.py
	python3 tools/check_weyl_rank_one_ideal.py --energy 8
	python3 tools/check_odd_weyl_spectrum.py
	python3 tools/check_mismatched_weyl_spectrum.py
	python3 tools/check_trace_kernel_support.py
	python3 tools/check_trace_weight_restriction.py --output /tmp/integer-factoring-trace-restriction.json
	python3 tools/check_partial_mode_sweep.py

verify-weyl-orbit:
	python3 tools/check_weyl_rank_one_ideal.py

verify-witt-full:
	python3 tools/witt_artin_hasse.py --h 32 --gram

.PHONY: verify-pr27 verify-code-width5
verify-pr27:
	python3 tools/check_weighted_cut_rank_budget.py
	python3 tools/check_weighted_cut_residual_budget.py
	python3 tools/check_weighted_cut_dual_certificate.py
	python3 tools/check_two_product_complement.py
	python3 tools/check_two_product_seed_cutoff.py
	python3 tools/check_flagged_cut_rank_budget.py
	python3 tools/check_geometric_subspace_compression.py
	python3 tools/check_geometric_cut_rank_budget.py
	python3 tools/check_weighted_cut_relaxation_limit.py
	python3 tools/check_geometric_cut_relaxation_limit.py
	python3 tools/check_weighted_cut_border_stability.py
	python3 tools/check_geometric_border_cutoff.py
	python3 tools/check_rank_one_formal_cutoff.py
	python3 tools/check_quadratic_crossing.py
	python3 tools/check_word_restriction_barrier.py
	python3 tools/check_code_sections.py
	python3 tools/check_frame_filtration.py
	python3 tools/check_feistel_filtration.py
	python3 tools/check_diagonal_word_rank_sum.py
	python3 tools/check_incidence_kernel_rank_sum.py
	python3 tools/check_full_conservation_rank_sum.py
	python3 tools/check_p_group_convolution_degree.py
	python3 tools/check_bi_exterior_crossing_algebra.py
	python3 tools/check_bi_exterior_degree_budget.py
	python3 tools/check_nonlinear_layer_probe.py
	python3 tools/check_pauli_walsh_dictionary.py
	python3 tools/check_pauli_walsh_pairs.py
	python3 tools/check_growing_number_field_transfer.py
	python3 tools/check_rank_sum_global_transfer.py
	python3 tools/check_characteristic_independent_recovery.py

verify-code-width5:
	g++ -O3 -std=c++20 -fopenmp tools/check_code_kernel_width5.cpp -o /tmp/integer-factoring-code-width5
	/tmp/integer-factoring-code-width5 0 64 4 /tmp/integer-factoring-code-width5
	python3 tools/check_code_sections.py --full-classifier --native-prefix /tmp/integer-factoring-code-width5

# beyond-half-checks:start
.PHONY: verify-beyond-half verify-norm-processing verify-adaptive-channels verify-crossing-pressure-certificates verify-protected-scalar-pressure verify-reflection-codes-full

# Core exact replays; exhaustive coordinate-code classification has its own target.
verify-beyond-half: verify-norm-processing verify-adaptive-channels verify-crossing-pressure-certificates verify-protected-scalar-pressure

verify-norm-processing:
	python3 tools/check_prime_form_bias.py
	python3 tools/check_prime_form_cutoffs.py
	python3 tools/check_direct_discriminant_short_norms.py
	python3 tools/check_early_ambiguous_factor.py
	python3 tools/check_factored_norm_representation.py
	python3 tools/check_prime_form_factor_base.py
	python3 tools/check_smooth_band_density.py
	python3 tools/check_powered_prime_norm_parameters.py
	python3 tools/check_powered_prime_norm_witnesses.py
	python3 tools/check_quasipolynomial_norm_parameters.py
	python3 tools/check_relation_two_saturation.py
	python3 tools/check_weighted_norm_relations.py
	python3 tools/check_arbitrary_norm_representation.py
	python3 tools/check_selected_norm_factoring.py
	python3 tools/check_selected_norm_factoring_124_bit.py
	python3 tools/check_arbitrary_selected_norm_factoring.py
	python3 tools/check_dynamic_leading_norm_relations.py
	python3 tools/check_leading_power_smooth_relations.py
	python3 tools/check_fresh_leading_power_control.py
	python3 tools/check_odd_pivot_two_primary_elimination.py
	python3 tools/check_adaptive_leading_power_windows.py

verify-adaptive-channels:
	python3 tools/check_adaptive_smith_lattices.py
	python3 tools/check_matroid_smith_channels.py
	python3 tools/check_dvr_smith_channels.py
	python3 tools/check_subspace_smith_channels.py
	python3 tools/check_channel_pressure_separation.py
	python3 tools/check_block_subspace_channels.py

verify-crossing-pressure-certificates:
	python3 tools/check_direct_type_elementary_floor.py
	python3 tools/check_arbitrary_binary_code_pressure.py
	python3 tools/check_integer_shared_minor_floor.py
	python3 tools/check_adjugate_seven_letter_floor.py
	python3 tools/check_bounded_order_smith_floor.py
	python3 tools/check_clifford_compression.py
	python3 tools/check_clifford_threshold_sharpness.py
	python3 tools/check_homogeneous_restriction_floor.py
	python3 tools/check_joint_height_conjugations.py
	python3 tools/check_eight_letter_reflection_retraction.py
	python3 tools/check_eight_letter_pressure_factorization.py
	python3 tools/check_polynomial_joint_height_rigidity.py

verify-protected-scalar-pressure:
	python3 tools/check_triangular_auxiliary_row_fiber.py
	python3 tools/check_protected_subtensor_pressure.py
	python3 tools/check_scalar_terminal_projections.py

verify-reflection-codes-full:
	python3 tools/check_reflection_coordinate_codes.py --full --progress --output evidence/reflection-coordinate-codes.json

.PHONY: verify-adaptive-power-full
verify-adaptive-power-full:
	python3 tools/check_adaptive_leading_power_windows.py --full

# beyond-half-checks:end

.PHONY: verify-constructions
verify-constructions: verify-pr31-constructions
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
	python3 tools/check_ward_fan_entry_queries.py --output /tmp/integer-factoring-ward-fan-entry-queries.json
	python3 tools/check_ward_fan_leading_block.py --output /tmp/integer-factoring-ward-fan-leading-block.json
	python3 tools/check_ward_fan_jet_probe.py --output /tmp/integer-factoring-ward-fan-jet-probe.json
	python3 tools/check_ward_operator_types.py --output /tmp/integer-factoring-ward-operator-types-check.json
	python3 tools/check_overlapping_symplectic_profiles.py
	python3 tools/check_cubic_radical_output_maps.py
	python3 tools/check_cubic_cyclic_fourier_bound.py --output /tmp/integer-factoring-cubic-cyclic-fourier-bound.json
	python3 tools/check_cubic_jordan_filtration.py --output /tmp/integer-factoring-cubic-jordan-filtration.json
	python3 tools/check_cubic_fixed_image_rank_one.py --output /tmp/integer-factoring-cubic-fixed-image-rank-one.json
	python3 tools/check_cubic_coordinate_image_mixtures.py --output /tmp/integer-factoring-cubic-coordinate-image-mixtures.json
	python3 tools/check_cubic_varying_image_escape.py --output /tmp/integer-factoring-cubic-varying-image-escape.json
	python3 tools/check_cubic_balanced_peels.py --output /tmp/integer-factoring-cubic-balanced-peels.json
	python3 tools/check_mixed_common_functional_barrier.py --output /tmp/integer-factoring-mixed-common-functional-barrier.json
	python3 tools/check_noncommutative_crossing_seed_algebra_audit.py
	python3 tools/check_symplectic_seed_tensor_modules.py
	python3 tools/symplectic_seed_weighted_modules.py --output /tmp/integer-factoring-weighted-modules.json
	python3 tools/check_noncommutative_crossing_character_compression.py --output /tmp/integer-factoring-character-compression.json
	python3 tools/check_symplectic_seed_coefficient_barrier.py --output /tmp/integer-factoring-seed-coefficient-barrier.json
	python3 tools/check_double_endpoint_retracts.py --cases full:3:3
	python3 tools/check_noncommutative_crossing_group_degeneration.py


.PHONY: verify-native-crossings
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
	$(CXX) -O3 -std=c++20 -Wall -Wextra -Werror tools/search_cubic_pairing_regroupings.cpp -o build/native/cubic-pairing-regroupings
	python3 tools/explore_cubic_pairing_regroupings.py --native $(CURDIR)/build/native/cubic-pairing-regroupings --replay evidence/cubic-pairing-regroupings.json --output build/native/cubic-pairing-regroupings-audit.json


.PHONY: verify-pr29-extended
verify-pr29-extended:
	python3 tools/check_common_character_compressions.py --output /tmp/integer-factoring-common-compressions.json
	python3 tools/check_two_even_site_sections.py --scope graphs --work-dir /tmp/integer-factoring-two-even-sections --output /tmp/integer-factoring-two-even-sections.json
	python3 tools/assemble_two_even_site_cores.py --pairs 6 --screen evidence/two-even-site-core-P6-F2.jsonl --output /tmp/integer-factoring-two-even-P6.json
	python3 tools/check_two_even_site_compressions.py --output /tmp/integer-factoring-two-even-P6-compressions.json
	python3 tools/check_two_even_site_jordan_compressions.py --output /tmp/integer-factoring-two-even-P7-compressions.json
	python3 tools/check_two_even_site_P8_compressions.py --output /tmp/integer-factoring-two-even-P8-compressions.json
	python3 tools/check_two_even_site_row_corners.py --output /tmp/integer-factoring-two-even-row-corners.json
	python3 tools/check_two_even_site_frobenius_probe.py --output /tmp/integer-factoring-two-even-P10-frobenius.json
	python3 tools/check_ward_fan_boundary.py --output /tmp/integer-factoring-ward-fan-boundary.json
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

.PHONY: verify-parity
verify-parity:
	python3 tools/check_parity_power_divisor_recovery.py
	python3 tools/check_weighted_network_ports.py

.PHONY: verify-recovery
verify-recovery: verify-parity
	python3 tools/check_fixed_characteristic_recovery.py
	python3 tools/check_root_threshold_counterexample.py
	python3 tools/check_unitary_orbit_recovery.py
	python3 tools/check_quadratic_congruence_recovery.py
	python3 tools/check_high_order_characteristic.py
	python3 tools/check_fixed_precision_count_lifting.py

.PHONY: verify-constructive
verify-constructive:
	python3 tools/check_explicit_entanglers.py
	python3 tools/check_explicit_entanglers.py --control evidence/binary-entangler-extension-control.json
	python3 tools/check_smith_profile_transfer.py
	python3 tools/check_degree_sensitive_smith_paths.py
	python3 tools/check_tropical_smith_channel_words.py
	python3 tools/check_integral_schur_modes.py
	python3 tools/check_schur_construction.py
	python3 tools/check_ternary_resolvent_probe.py
	python3 tools/check_trace_shear_general.py
	python3 tools/check_auxiliary_gauge_smith.py
	python3 tools/check_auxiliary_lattice_gauges.py
	python3 tools/check_balanced_physical_gauge.py
	python3 tools/check_sum_preserving_smith_barrier.py

.PHONY: verify-native
verify-native: verify-exterior verify-entangler-density verify-polynomial-native

.PHONY: verify-schur-ranks
verify-schur-ranks:
	python3 tools/full_schur_crossing.py --weights 1 2 3 4 5 --prime 2 --output /tmp/full-schur-f2.json
	python3 tools/homogeneous_full_schur.py --weights 2 4 8 --output /tmp/homogeneous-full-schur-f2.json

.PHONY: verify-exterior
verify-exterior:
	python3 tools/check_coherent_exterior_boundaries.py

.PHONY: verify-entangler-density
verify-entangler-density:
	python3 tools/check_density16_certificate.py

.PHONY: verify-pole-arithmetic
verify-pole-arithmetic:
	python3 tools/check_pole_smith_flow.py
	python3 tools/check_pole_smith_balanced.py
	python3 tools/check_pole_cauchy_conditioning.py
	python3 tools/check_pole_cauchy_polynomial_solve.py
	python3 tools/check_modular_interval_jets.py

.PHONY: verify-cofactor-collisions
verify-cofactor-collisions:
	python3 tools/check_cofactor_collision_relations.py
	python3 tools/check_cofactor_collision_quota.py

.PHONY: verify-coherent-retractions verify-coherent-q32-catalogue
verify-coherent-retractions:
	python3 tools/check_coherent_hot_compounds.py --output /tmp/coherent-hot-compound-controls.json
	python3 tools/check_coherent_hot_fast.py --output /tmp/coherent-hot-fast-controls.json
	python3 tools/check_coherent_hot_mixed.py --output /tmp/coherent-hot-mixed-controls.json
	python3 tools/check_coherent_q32_pressure.py --output /tmp/coherent-q32-pressure.json

verify-coherent-q32-catalogue:
	python3 tools/check_coherent_q32_catalogue.py --output /tmp/coherent-q32-catalogue-optimum.json

.PHONY: verify-cube-transfer
verify-cube-transfer:
	python3 tools/check_cube_crossing_transfer.py --output /tmp/cube-crossing-transfer.json

.PHONY: verify-cube-gaussian
verify-cube-gaussian:
	python3 tools/check_cube_homogeneous_gaussian.py --output /tmp/cube-homogeneous-gaussian.json

.PHONY: verify-two-axis-cube-regrouping
verify-two-axis-cube-regrouping:
	PYTHONPATH=tools python3 tools/check_two_axis_cube_regrouping.py

.PHONY: verify-two-axis-global-pairing
verify-two-axis-global-pairing:
	PYTHONPATH=tools python3 tools/check_two_axis_global_pairing.py

.PHONY: verify-pr31-constructions
verify-pr31-constructions:
	python3 tools/check_interior_jacobian_budget.py --output /tmp/integer-factoring-interior-jacobian-budget.json
	python3 tools/check_rational_cut_mesh_budget.py --output /tmp/integer-factoring-rational-cut-mesh-budget.json
	python3 tools/check_compatible_compression_budget.py --cores evidence/compatible-compression-cores.json --output /tmp/integer-factoring-compatible-compression-budget.json
	python3 tools/check_compatible_formal_cutoff.py --cores evidence/compatible-compression-cores.json --output /tmp/integer-factoring-compatible-formal-cutoff.json
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

.PHONY: verify-pr31-frontier
verify-pr31-frontier:
	python3 tools/check_strict_mixed_cut_gate.py --output /tmp/integer-factoring-strict-mixed-cut-gate.json
	python3 tools/check_cnot_sixteen_span.py --output /tmp/integer-factoring-cnot-sixteen-span.json
	python3 tools/check_fixed_prime_factor_transfer.py --output /tmp/integer-factoring-fixed-prime-factor-transfer.json
	python3 tools/check_mixed_cut_rank_budget.py --output /tmp/integer-factoring-mixed-cut-rank-budget.json
	python3 tools/check_parity_factor_transfer.py --output /tmp/integer-factoring-parity-factor-transfer.json
	python3 tools/check_rough_sparse_transfer.py --output /tmp/integer-factoring-rough-sparse-transfer.json

.PHONY: verify-cofactor-ideal-closure
verify-cofactor-ideal-closure:
	python3 tools/check_cofactor_ideal_closure.py --output /tmp/cofactor-ideal-closure.json
	python3 tools/check_cofactor_ideal_closure.py --replay evidence/cofactor-ideal-closure-92.json
	python3 tools/check_cofactor_ideal_closure.py --replay evidence/cofactor-ideal-closure-124.json
	python3 tools/check_cofactor_ideal_closure.py --replay evidence/cofactor-ideal-closure-144.json

.PHONY: verify-mixed-place
verify-mixed-place:
	python3 tools/check_weighted_multipoint_hermite.py
	python3 tools/check_mixed_place_hermite.py
	python3 tools/check_ternary_basis_mixed.py

.PHONY: verify-adaptive-mixed-place
verify-adaptive-mixed-place:
	python3 tools/check_adaptive_mixed_place.py --certificate evidence/fused-q4-adaptive-certificate.json

.PHONY: verify-network-codecs
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

.PHONY: verify-stable-hot-codes
verify-stable-hot-codes:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_implicit_reflection_codes.py --output /tmp/integer-factoring-implicit-hot.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_natural_cut_duality.py --output /tmp/integer-factoring-natural-duality.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_hot_code_charge_profiles.py --output /tmp/integer-factoring-hot-charge.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_weighted_hot_compounds.py --output /tmp/integer-factoring-weighted-compound.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_weighted_hot_code_oracle.py --output /tmp/integer-factoring-weighted-hot-oracle.json

.PHONY: verify-hot64-certificate
verify-hot64-certificate:
	OPENBLAS_NUM_THREADS=1 python3 tools/certify_hot_charge_pressure.py --size 64 --modes 12 --base 459 460 --support 64 45 --output /tmp/integer-factoring-hot64-certificate.json

.PHONY: verify-finite-crossing-certificates
verify-finite-crossing-certificates:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_block_amplification.py --output /tmp/integer-factoring-finite-amplification.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_typed_adaptive_certificate.py --output /tmp/integer-factoring-finite-typed.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_rational_typed_adaptive_bounds.py --output /tmp/integer-factoring-rational-finite.json

.PHONY: verify-frobenius-stars
verify-frobenius-stars:
	python3 tools/check_frobenius_rational_trace.py --output /tmp/integer-factoring-rational-star.json
	python3 tools/check_fermat_torus_binomials.py --output /tmp/integer-factoring-torus-star.json
	python3 tools/check_frobenius_trace_infinity.py --output /tmp/integer-factoring-star-infinity.json
	python3 tools/check_relative_frobenius_trace.py --output /tmp/integer-factoring-relative-trace.json
	python3 tools/check_semisimple_cyclic_box_contraction.py --output /tmp/integer-factoring-cyclic-box.json

.PHONY: verify-observable-constructors
verify-observable-constructors:
	python3 tools/check_observable_boundary_constructor.py --output /tmp/integer-factoring-observable-boundary.json
	python3 tools/check_lattice_boundary_decoder.py --output /tmp/integer-factoring-lattice-boundary.json

.PHONY: verify-prime-constructors
verify-prime-constructors:
	python3 tools/check_power_two_progression_dyadic_prime.py --output /tmp/integer-factoring-dyadic-prime.json
	python3 tools/check_prime_transform_boundary_constructor.py --output /tmp/integer-factoring-prime-transform.json
	python3 tools/check_width_sensitive_prime_constructor.py --output /tmp/integer-factoring-width-sensitive.json
	python3 tools/check_polynomial_entry_prime_constructor.py --output /tmp/integer-factoring-polynomial-entry.json
	python3 tools/check_rough_observable_factor_constructor.py --output /tmp/integer-factoring-rough-observable.json
	python3 tools/check_prime_root_field_construction.py --output /tmp/integer-factoring-prime-root-field.json

.PHONY: verify-degree16-crossing
verify-degree16-crossing:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_degree16_reflection_crossing.py --output /tmp/integer-factoring-degree16-replay.json

.PHONY: verify-reflection-truncations
verify-reflection-truncations:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_truncated_reflection_crossings.py --output /tmp/integer-factoring-truncation-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_truncated_reflection4_search.py --output /tmp/integer-factoring-truncation-pressure-replay.json

.PHONY: verify-unsplit-monoid
verify-unsplit-monoid:
	python3 tools/check_unsplit_monoid_contraction.py --output /tmp/integer-factoring-unsplit-monoid-replay.json

.PHONY: verify-paired-trace
verify-paired-trace:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_paired_trace.py --output /tmp/integer-factoring-paired-trace-replay.json

.PHONY: verify-companion-search
verify-companion-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_coupled_companions.py --output /tmp/integer-factoring-companion-search-replay.json

.PHONY: verify-monoid-blocks
verify-monoid-blocks:
	python3 tools/check_multiplicative_monoid_blocks.py --output /tmp/integer-factoring-monoid-blocks-replay.json

.PHONY: verify-triangular-series
verify-triangular-series:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_triangular_reflection_jet_families.py --output /tmp/integer-factoring-triangular-jets-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_laurent_series_coupled_frames.py --output /tmp/integer-factoring-laurent-series-replay.json

.PHONY: verify-uniform-reflection-jets
verify-uniform-reflection-jets:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_low_jet_frames.py --include-interval16 --output /tmp/integer-factoring-uniform-reflection-jets-replay.json

.PHONY: verify-coupled-searches
verify-coupled-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection4_coupled_search.py --output /tmp/integer-factoring-coordinate-search-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_rank_one_pressure.py --output /tmp/integer-factoring-dense-pressure-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_affine_coupled_frames.py --output /tmp/integer-factoring-affine-frames-replay.json

.PHONY: verify-constant-carry-weights
verify-constant-carry-weights:
	python3 tools/check_constant_field_carry_weights.py --output /tmp/integer-factoring-constant-carry-weights-replay.json
	python3 tools/check_first_jet_carry_weights.py --output /tmp/integer-factoring-first-jet-carry-weights-replay.json

.PHONY: verify-square-zero-frames
verify-square-zero-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_square_zero_coupled_frames.py --output /tmp/integer-factoring-square-zero-first-jet-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_square_zero_second_jet.py --output /tmp/integer-factoring-square-zero-second-jet-replay.json

.PHONY: verify-dense-coupled-frames
verify-dense-coupled-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_dense_coupled_reflection_frames.py --catalogue evidence/dense-rank-one-coupled-frames.json --output /tmp/integer-factoring-dense-coupled-frames-replay.json

.PHONY: verify-weighted-carries
verify-weighted-carries:
	python3 tools/check_weighted_dyadic_carry.py --output /tmp/integer-factoring-weighted-carries-replay.json

.PHONY: verify-physical-frames
verify-physical-frames:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_physical_polynomial_frames.py --output /tmp/integer-factoring-polynomial-frames-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_laurent_frames.py --output /tmp/integer-factoring-laurent-frames-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection_coupled_frames.py --output /tmp/integer-factoring-coupled-frames-replay.json

.PHONY: verify-ternary-majorization
verify-ternary-majorization:
	python3 tools/check_ternary_gauge_majorization.py --output /tmp/integer-factoring-ternary-majorization-replay.json

.PHONY: verify-frobenius-recovery
verify-frobenius-recovery:
	python3 tools/check_frobenius_sparse_recovery.py --output /tmp/integer-factoring-frobenius-recovery-replay.json

.PHONY: verify-weighted-adaptive
verify-weighted-adaptive:
	python3 tools/check_weighted_adaptive_catalogues.py --output /tmp/integer-factoring-weighted-adaptive-replay.json

.PHONY: verify-basis-searches
verify-basis-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_ternary_basis_search.py --output /tmp/integer-factoring-ternary-basis-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_fused_weighted_winners.py evidence/fused-mixed-weighted-summary.json --output /tmp/integer-factoring-fused-winners-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_fused_physical_bases.py --output /tmp/integer-factoring-fused-bases-replay.json

.PHONY: verify-asymmetric-crossings
verify-asymmetric-crossings:
	python3 tools/check_asymmetric_support_line.py --output /tmp/integer-factoring-asymmetric-support-replay.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_asymmetric_pressure_oracle.py --output /tmp/integer-factoring-asymmetric-oracle-replay.json

.PHONY: verify-asymmetric-scaling-floor
verify-asymmetric-scaling-floor:
	python3 tools/check_asymmetric_scaling_floor.py --output /tmp/integer-factoring-asymmetric-scaling-floor.json

.PHONY: verify-explicit-fusion
verify-explicit-fusion:
	python3 tools/check_explicit_fused_row.py --output /tmp/integer-factoring-explicit-fused-row-replay.json
	python3 tools/check_fused_frobenius_spectators.py --output /tmp/integer-factoring-fused-spectators-replay.json

.PHONY: verify-carry-fourier
verify-carry-fourier:
	python3 tools/check_carry_fourier_elimination.py --output /tmp/integer-factoring-carry-fourier-replay.json

.PHONY: verify-paired-affine-space
verify-paired-affine-space:
	OPENBLAS_NUM_THREADS=1 python3 tools/saturate_paired_trace_affine_space.py evidence/independent-paired-affine-degree-basis.json.gz --output /tmp/integer-factoring-paired-affine-saturation.json --basis-output /tmp/integer-factoring-paired-affine-basis.json.gz

.PHONY: verify-bell-jets
verify-bell-jets:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_entrywise_reflection_constraints.py --output /tmp/integer-factoring-entrywise-jets.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_bell_jet_crt_crossings.py --output /tmp/integer-factoring-bell-crt.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_bell_jet_transversality.py --output /tmp/integer-factoring-bell-transversality.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_ramified_carrier_minors.py --output /tmp/integer-factoring-ramified-carrier-minors.json

.PHONY: verify-finite-fibers
verify-finite-fibers:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_regular_smith_frames.py --output /tmp/integer-factoring-f4-regular-frames.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_finite_jet_smith_fibers.py --output /tmp/integer-factoring-finite-jet-fibers.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_two_cut_smith_fibers.py --output /tmp/integer-factoring-f4-two-cut-fibers.json

.PHONY: verify-closed-place-oracles
verify-closed-place-oracles:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_closed_place_oracle.py --output /tmp/integer-factoring-f4-closed-oracle.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_paired_mixed_place_pressure.py --output /tmp/integer-factoring-paired-mixed-pressure.json

.PHONY: verify-filtered-affine-search
verify-filtered-affine-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_filtered_reflection16_affine_search.py --output /tmp/integer-factoring-filtered-affine-search.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_filtered_reflection16_affine_combinations.py --output /tmp/integer-factoring-filtered-affine-combinations.json

.PHONY: verify-laurent-search
verify-laurent-search:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_laurent_series_search.py --output /tmp/integer-factoring-laurent-series-search.json

.PHONY: verify-coherent-boundary
verify-coherent-boundary:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_coherent_hot_boundary.py --output /tmp/integer-factoring-coherent-hot8.json
	OPENBLAS_NUM_THREADS=1 python3 tools/checkpoint_coherent_hot_boundaries.py --output /tmp/integer-factoring-coherent-boundary-search.json

.PHONY: verify-root-worker-searches
verify-root-worker-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/checkpoint_root_worker_searches.py --output /tmp/integer-factoring-root-worker-searches.json

.PHONY: verify-hot96-pressure
verify-hot96-pressure:
	python3 tools/check_natural_pressure_intervals.py --output /tmp/integer-factoring-pressure-interval-controls.json
	python3 tools/certify_natural_pressure_interval.py evidence/hot96-charge-profile.json --rho 710005874881191 1000000000000000 --base 5571156672673465 5571156672673466 --output /tmp/integer-factoring-hot96-interval-pressure.json

.PHONY: verify-hot96-profile
verify-hot96-profile:
	OPENBLAS_NUM_THREADS=1 python3 tools/certify_hot_charge_pressure.py --size 96 --modes 12 --base 768 769 --support 98 69 --output /tmp/integer-factoring-hot96-full-replay.json

.PHONY: verify-final-searches
verify-final-searches:
	OPENBLAS_NUM_THREADS=1 python3 tools/checkpoint_algorithm_final_searches.py --output /tmp/integer-factoring-final-searches.json
	OPENBLAS_NUM_THREADS=1 python3 tools/check_unequal_coherent_hot_codes.py --output /tmp/integer-factoring-unequal-coherent-codes.json

.PHONY: verify-bilinear-cauchy
verify-bilinear-cauchy:
	python3 tools/check_bilinear_cauchy_closure.py --output /tmp/integer-factoring-bilinear-cauchy-closure.json
	python3 tools/check_cauchy_bilinear_phase_counterexample.py --output /tmp/integer-factoring-cauchy-phase-counterexample.json

.PHONY: verify-hot96-type
verify-hot96-type:
	python3 tools/certify_hot_rational_type_rounding.py --output /tmp/integer-factoring-hot96-small-type.json

.PHONY: verify-bilinear-pole-packets
verify-bilinear-pole-packets:
	python3 tools/check_bilinear_pole_packets.py

.PHONY: verify-cube-profile-pressure
verify-cube-profile-pressure:
	python3 tools/check_cube_gaussian_profile_budget.py --output /tmp/integer-factoring-cube-gaussian-profile-budget.json
	python3 tools/check_cube_gaussian_tight_pressure.py --output /tmp/integer-factoring-cube-gaussian-tight-pressure.json

.PHONY: verify-pole-cauchy-analytic
verify-pole-cauchy-analytic:
	python3 tools/check_pole_cauchy_analytic_weights.py --output /tmp/integer-factoring-pole-cauchy-analytic-weights.json
	python3 tools/check_pole_cauchy_analytic_kernels.py --output /tmp/integer-factoring-pole-cauchy-analytic-kernels.json
	python3 tools/check_regular_analytic_selector_precision.py --output /tmp/integer-factoring-regular-analytic-selector-precision.json

.PHONY: verify-pole-cauchy-domains
verify-pole-cauchy-domains:
	python3 tools/check_pole_cauchy_domain_precision.py --output /tmp/integer-factoring-pole-cauchy-domain-precision.json

.PHONY: verify-pole-cauchy-circuit-moments
verify-pole-cauchy-circuit-moments:
	python3 tools/check_pole_cauchy_circuit_moments.py

.PHONY: verify-pole-cauchy-resolvents
verify-pole-cauchy-resolvents:
	python3 tools/check_pole_cauchy_rational_resolvents.py

.PHONY: verify-prime-cofactor-saturation
verify-prime-cofactor-saturation:
	python3 tools/check_prime_cofactor_saturation.py --output /tmp/integer-factoring-prime-cofactor-saturation.json
	python3 tools/check_prime_cofactor_saturation.py --replay evidence/prime-cofactor-saturation-92.json
	python3 tools/check_prime_cofactor_saturation.py --replay evidence/prime-cofactor-saturation-124.json
	python3 tools/check_prime_cofactor_saturation.py --replay evidence/prime-cofactor-saturation-144.json
	python3 tools/check_prime_cofactor_saturation.py --replay evidence/prime-cofactor-saturation-selective-144.json

.PHONY: verify-cofactor-orbit-moments verify-entry-moment-cycles
verify-cofactor-orbit-moments:
	python3 tools/check_cofactor_orbit_moments.py --output /tmp/integer-factoring-cofactor-orbit-moments.json

verify-entry-moment-cycles:
	python3 tools/check_entry_moment_cycles.py --output /tmp/integer-factoring-entry-moment-cycles.json
	python3 tools/check_entry_moment_cycles.py --replay evidence/entry-moment-cycles.json

.PHONY: verify-output-space-transfer
verify-output-space-transfer:
	python3 tools/check_output_sensitive_space_transfer.py

.PHONY: verify-small-word-fourier
verify-small-word-fourier:
	python3 tools/check_small_word_fourier_constructor.py


.PHONY: verify-bell-walsh
verify-bell-walsh:
	python3 tools/check_bell_walsh_constructions.py evidence/bell-walsh-energy-search-d8x1.json evidence/bell-walsh-energy-search-d8x3.json evidence/bell-walsh-energy-search-d16x1.json evidence/bell-walsh-energy-search-d16x3.json

.PHONY: verify-isotonic-retiming
verify-isotonic-retiming:
	python3 tools/check_isotonic_coefficient_retiming.py


.PHONY: verify-resolvent-determinant-budget verify-register-toeplitz-budget verify-binary-rank-two-channels
verify-resolvent-determinant-budget:
	python3 tools/check_regular_resolvent_determinant_budget.py
verify-register-toeplitz-budget:
	python3 tools/check_shift_register_toeplitz_budget.py --output /tmp/integer-factoring-shift-register-toeplitz-budget.json
verify-binary-rank-two-channels:
	python3 tools/check_binary_rank_two_channels.py

.PHONY: verify-coherent-graded-quotients verify-coherent-shared-gauges
verify-coherent-graded-quotients:
	python3 tools/check_coherent_graded_quotients.py --output /tmp/integer-factoring-coherent-graded-quotients.json
verify-coherent-shared-gauges:
	python3 tools/check_coherent_shared_gauges.py --output /tmp/integer-factoring-coherent-shared-gauges.json

.PHONY: verify-f4-hyperoval
verify-f4-hyperoval:
	python3 tools/check_coherent_hot_f4_hyperoval_entries.py
	python3 tools/certify_coherent_hot_f4_hyperoval.py

.PHONY: verify-four-axis-floor
verify-four-axis-floor:
	python3 tools/check_four_axis_gaussian_profile_floor.py

.PHONY: verify-weighted-four-axis-floor
verify-weighted-four-axis-floor:
	python3 tools/check_weighted_four_axis_gaussian_profile_floor.py

.PHONY: verify-multiaxis-floor
verify-multiaxis-floor:
	python3 tools/check_multiaxis_gaussian_profile_floor.py


.PHONY: verify-pole-cauchy-rational-blocks
verify-pole-cauchy-rational-blocks:
	python3 tools/check_pole_cauchy_rational_poles.py --small-only --output /tmp/integer-factoring-pole-cauchy-rational-poles.json
	python3 tools/check_pole_cauchy_repeated_poles.py --output /tmp/integer-factoring-pole-cauchy-repeated-poles.json
	python3 tools/check_pole_cauchy_quadratic_selectors.py --output /tmp/integer-factoring-pole-cauchy-quadratic-selectors.json
	python3 tools/check_pole_cauchy_residue_blocks.py --output /tmp/integer-factoring-pole-cauchy-residue-blocks.json

.PHONY: verify-coherent-boundary-floor
verify-coherent-boundary-floor:
	python3 tools/check_coherent_boundary_protected_floor.py --recompute-profiles

.PHONY: verify-coherent-f4-boundary
verify-coherent-f4-boundary:
	python3 tools/check_coherent_f4_boundary.py --catalogue evidence/coherent-f4-boundary-orbits.json


.PHONY: verify-coherent-q96-pressure verify-coherent-mixed-catalogue
verify-coherent-q96-pressure:
	python3 tools/check_coherent_q96_pressure.py --output /tmp/integer-factoring-coherent-q96-pressure.json
verify-coherent-mixed-catalogue:
	python3 tools/check_coherent_q32_mixed_catalogue.py --output /tmp/integer-factoring-coherent-q32-mixed-boundary-optimum.json


.PHONY: verify-pole-cauchy-batched-poles
verify-pole-cauchy-batched-poles:
	python3 tools/check_pole_cauchy_pole_polynomials.py --output /tmp/integer-factoring-pole-cauchy-pole-polynomials.json

.PHONY: verify-rough-length-fourier
verify-rough-length-fourier:
	python3 tools/check_rough_length_fourier_constructor.py --output /tmp/integer-factoring-rough-length-fourier-constructor.json

.PHONY: verify-pole-cauchy-localization
verify-pole-cauchy-localization:
	python3 tools/check_pole_cauchy_localization.py --output /tmp/integer-factoring-pole-cauchy-localization.json

.PHONY: verify-coherent-protected-height25
verify-coherent-protected-height25:
	python3 tools/check_coherent_protected_height25.py --recompute-profiles --output /tmp/integer-factoring-coherent-protected-height25.json


.PHONY: verify-prime-cofactor-fresh-160
verify-prime-cofactor-fresh-160:
	python3 tools/check_fresh_prime_cofactor_control.py --bits 160 --replay evidence/prime-cofactor-fresh-160-fixed.json.gz
	python3 tools/check_fresh_prime_cofactor_control.py --bits 160 --replay evidence/prime-cofactor-fresh-160-selective.json.gz


.PHONY: verify-pole-cauchy-general-weights
verify-pole-cauchy-general-weights:
	python3 tools/check_pole_cauchy_rational_weights.py --output /tmp/integer-factoring-pole-cauchy-rational-weights.json
	python3 tools/check_pole_cauchy_restricted_kernels.py --output /tmp/integer-factoring-pole-cauchy-restricted-kernels.json
	python3 tools/check_rational_interval_selector_bound.py --output /tmp/integer-factoring-rational-interval-selector-bound.json


.PHONY: verify-cube-shared-budgets
verify-cube-shared-budgets:
	python3 tools/check_cube_homogeneous_profile_budget.py --output /tmp/integer-factoring-cube-homogeneous-profile-budget.json
	python3 tools/check_cube_graded_minor_budget.py --output /tmp/integer-factoring-cube-graded-minor-budget.json
	python3 tools/check_cube_scalar_color_clone_floor.py --output /tmp/integer-factoring-cube-scalar-color-clone-floor.json
	python3 tools/check_cube_equal_height_fock_floor.py --output /tmp/integer-factoring-cube-equal-height-fock-floor.json

.PHONY: verify-coherent-full-tail
verify-coherent-full-tail:
	python3 tools/check_coherent_full_tail_core.py --evidence evidence/coherent-full-tail-core.json --output /tmp/integer-factoring-coherent-full-tail-controls.json


.PHONY: verify-coherent-free-tail verify-coherent-asymmetric-gauges
verify-coherent-free-tail:
	python3 tools/check_coherent_free_tail.py --output /tmp/integer-factoring-coherent-free-tail-controls.json
verify-coherent-asymmetric-gauges:
	python3 tools/check_coherent_asymmetric_gauges.py --output /tmp/integer-factoring-coherent-asymmetric-gauge-controls.json


.PHONY: verify-cube-generic-smith
verify-cube-generic-smith:
	python3 tools/check_cube_generic_smith128.py --output /tmp/integer-factoring-cube-generic-smith128.json

.PHONY: verify-f4-global-restarts
verify-f4-global-restarts:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_recentered_walks.py evidence/f4-global-restart-controls-search.json --output /tmp/integer-factoring-f4-global-restart-controls.json


.PHONY: verify-pole-cauchy-spectral-trace
verify-pole-cauchy-spectral-trace:
	python3 tools/check_pole_cauchy_spectral_trace.py --output /tmp/integer-factoring-pole-cauchy-spectral-trace.json


.PHONY: verify-prime-cofactor-fresh-176
verify-prime-cofactor-fresh-176:
	python3 tools/check_fresh_prime_cofactor_control.py --bits 176 --replay evidence/prime-cofactor-fresh-176-selective.json.xz

.PHONY: verify-fixed-rational-frames
verify-fixed-rational-frames:
	python3 tools/check_fixed_rational_frame_transfer.py --output /tmp/integer-factoring-fixed-rational-frame-transfer.json

.PHONY: verify-hot-boundary-floor
verify-hot-boundary-floor:
	python3 tools/check_hot_boundary_core_barrier.py


.PHONY: verify-coherent-physical-tail
verify-coherent-physical-tail:
	python3 tools/check_coherent_physical_tail_core.py --output /tmp/integer-factoring-coherent-physical-tail-core.json
	python3 tools/check_coherent_physical_tail_low_band.py --evidence evidence/coherent-physical-tail-low-band.json --output /tmp/integer-factoring-coherent-physical-tail-low-band.json


.PHONY: verify-cube-pair-specialization
verify-cube-pair-specialization:
	python3 tools/check_cube_pair_encoding_specialization.py --output /tmp/integer-factoring-cube-pair-specialization.json


.PHONY: verify-coherent-tail-injectivity
verify-coherent-tail-injectivity:
	python3 tools/check_coherent_tail_injectivity.py --output /tmp/integer-factoring-coherent-tail-injectivity.json

.PHONY: verify-joint-grid-compressions
verify-joint-grid-compressions:
	python3 tools/check_joint_grid_compressions.py

.PHONY: verify-joint-split-jet-guards
verify-joint-split-jet-guards:
	python3 tools/check_joint_split_jet_guards.py

.PHONY: verify-joint-split-flags
verify-joint-split-flags:
	python3 tools/check_joint_split_flag_queries.py

.PHONY: verify-joint-moving-retractions
verify-joint-moving-retractions:
	python3 tools/check_joint_moving_retractions.py

.PHONY: verify-joint-small-rational
verify-joint-small-rational:
	python3 tools/check_joint_small_rational_retractions.py

.PHONY: verify-joint-six-to-five
verify-joint-six-to-five:
	python3 tools/check_joint_six_to_five_retractions.py

.PHONY: verify-zero-energy-parity-core
verify-zero-energy-parity-core:
	python3 tools/check_zero_energy_parity_core.py

.PHONY: verify-zero-energy-hereditary
verify-zero-energy-hereditary:
	python3 tools/check_zero_energy_hereditary_adjugate.py

.PHONY: verify-hereditary-adjugate
verify-hereditary-adjugate:
	python3 tools/check_hereditary_adjugate_caps.py

.PHONY: verify-hereditary-q6-rejection
verify-hereditary-q6-rejection:
	python3 tools/check_hereditary_q6_rectangular_rejection.py

.PHONY: verify-multichannel-cofactors
verify-multichannel-cofactors:
	python3 tools/check_multichannel_cofactor_caps.py

.PHONY: verify-multichannel-bare-entropy
verify-multichannel-bare-entropy:
	python3 tools/check_multichannel_bare_entropy.py

.PHONY: verify-inherited-second-jet
verify-inherited-second-jet:
	python3 tools/check_inherited_second_jet_caps.py

.PHONY: verify-q12-second-jet-target
verify-q12-second-jet-target:
	python3 tools/check_q12_second_jet_target.py

.PHONY: verify-higher-jet-words
verify-higher-jet-words:
	python3 tools/check_higher_jet_word_caps.py

.PHONY: verify-q12-integer-word-target
verify-q12-integer-word-target:
	python3 tools/check_q12_integer_word_target.py

.PHONY: verify-q13-box-target-exclusion
verify-q13-box-target-exclusion:
	python3 tools/check_q13_box_target_exclusion.py

.PHONY: verify-q13-joint-word-target
verify-q13-joint-word-target:
	PYTHONPATH=tools python3 -c 'from joint_bare_profile_search import majorant_controls; print(majorant_controls())'
	python3 tools/check_q13_joint_word_target.py

.PHONY: verify-negation-odd-qutrit
verify-negation-odd-qutrit:
	python3 tools/check_negation_odd_qutrit.py

.PHONY: verify-nonclifford-css
verify-nonclifford-css:
	python3 tools/check_nonclifford_css_crossing.py

.PHONY: verify-box-adjugate-words
verify-box-adjugate-words:
	python3 tools/check_box_adjugate_word_caps.py

.PHONY: verify-graded-box-minors
verify-graded-box-minors:
	python3 tools/check_graded_box_minor_caps.py

.PHONY: verify-high-axis-frobenius-residues
verify-high-axis-frobenius-residues:
	python3 tools/check_high_axis_frobenius_residues.py

.PHONY: verify-global-axis-alignment
verify-global-axis-alignment:
	python3 tools/check_global_axis_alignment.py

.PHONY: verify-multiaxis-cube
verify-multiaxis-cube:
	python3 tools/check_multiaxis_cube_transfer.py


.PHONY: verify-cube-coefficient-fibers verify-cube-coefficient-fibers-full
verify-cube-coefficient-fibers:
	python3 tools/check_cube_coefficient_fibers.py --output /tmp/integer-factoring-cube-coefficient-fibers.json
verify-cube-coefficient-fibers-full:
	python3 tools/check_cube_coefficient_fibers.py --recompute-census --census /tmp/integer-factoring-cube-coefficient-census.json --output /tmp/integer-factoring-cube-coefficient-fibers-full.json

.PHONY: verify-multiaxis-types
verify-multiaxis-types:
	python3 tools/check_multiaxis_type_amplification.py

.PHONY: verify-product-power-weights
verify-product-power-weights:
	python3 tools/check_pole_cauchy_product_power_weights.py


.PHONY: verify-routed-rectangular-fourier
verify-routed-rectangular-fourier:
	python3 tools/check_routed_rectangular_fourier.py --output /tmp/integer-factoring-routed-rectangular-fourier.json

.PHONY: verify-symmetric-exterior
verify-symmetric-exterior:
	python3 tools/check_cube_symmetric_exterior_blocks.py --output /tmp/integer-factoring-cube-symmetric-exterior-blocks.json

.PHONY: verify-interval-bezout verify-balanced-bezout
verify-interval-bezout:
	python3 tools/check_interval_bezout_precision.py --output /tmp/integer-factoring-interval-bezout-precision.json
verify-balanced-bezout:
	python3 tools/check_balanced_bezout_denominator.py --output /tmp/integer-factoring-balanced-bezout-denominator.json


.PHONY: verify-cube-flat-parity-floor
verify-cube-flat-parity-floor:
	python3 tools/check_cube_flat_parity_floor.py

.PHONY: verify-balanced-bezout-modular
verify-balanced-bezout-modular:
	python3 tools/check_balanced_bezout_modular.py --output /tmp/integer-factoring-balanced-bezout-modular.json


.PHONY: verify-coherent-moving-pair
verify-coherent-moving-pair:
	python3 tools/check_coherent_moving_pair_core.py --evidence evidence/coherent-moving-pair-core.json --output /tmp/integer-factoring-coherent-moving-pair.json


.PHONY: verify-cube-f4-generic-frames
verify-cube-f4-generic-frames:
	python3 tools/check_cube_f4_bivariate_smith.py --output /tmp/integer-factoring-cube-f4-generic-frames.json


.PHONY: verify-cube-flat-low-order
verify-cube-flat-low-order:
	python3 tools/cube_flat_seeds.py
	python3 tools/check_cube_flat_low_order_floors.py

.PHONY: verify-matched-progression-bezout
verify-matched-progression-bezout:
	python3 tools/check_matched_progression_bezout.py --output /tmp/matched-progression-bezout.json

.PHONY: verify-tensor-pascal-boxes
verify-tensor-pascal-boxes:
	python3 tools/check_tensor_pascal_box_crossings.py --output /tmp/tensor-pascal-box-crossings.json


.PHONY: verify-cube-flat-coefficient-fibers
verify-cube-flat-coefficient-fibers:
	python3 tools/check_cube_flat_coefficient_fibers.py --compare evidence/cube-flat-coefficient-fibers.json --output /tmp/integer-factoring-cube-flat-coefficient-fibers.json


.PHONY: verify-local-parity-gauge
verify-local-parity-gauge:
	python3 tools/check_local_parity_gauge.py --output /tmp/integer-factoring-local-parity-gauge.json

.PHONY: verify-rectangular-bezout-denominator
verify-rectangular-bezout-denominator:
	python3 tools/check_rectangular_bezout_denominator.py --output /tmp/rectangular-bezout-denominator.json
	python3 tools/check_rectangular_bezout_scale.py --output /tmp/rectangular-bezout-scale.json

.PHONY: verify-rectangular-bezout-modular
verify-rectangular-bezout-modular:
	python3 tools/check_rectangular_bezout_modular.py --output /tmp/rectangular-bezout-modular.json

.PHONY: verify-weighted-pascal-flags
verify-weighted-pascal-flags:
	python3 tools/check_weighted_pascal_flags.py --output /tmp/weighted-pascal-flags.json

.PHONY: verify-pole-cauchy-hadamard-rank
verify-pole-cauchy-hadamard-rank:
	python3 tools/check_pole_cauchy_hadamard_rank.py --output /tmp/pole-cauchy-hadamard-rank.json

.PHONY: verify-quadratic-coalescence-bezout
verify-quadratic-coalescence-bezout:
	python3 tools/check_quadratic_coalescence_bezout.py --output /tmp/quadratic-coalescence-bezout.json


.PHONY: verify-canonical-affine-codes
verify-canonical-affine-codes:
	python3 tools/check_canonical_affine_codes.py --full --output /tmp/canonical-affine-codes.json


.PHONY: verify-regularized-bezout-factor-moments verify-displaced-bezout-extraction verify-leading-bezout-factoring
verify-regularized-bezout-factor-moments:
	python3 tools/check_regularized_bezout_factor_moments.py --output /tmp/regularized-bezout-factor-moments.json

verify-displaced-bezout-extraction:
	python3 tools/check_displaced_bezout_extraction.py --output /tmp/displaced-bezout-extraction.json

verify-leading-bezout-factoring:
	python3 tools/check_leading_bezout_factoring.py --output /tmp/leading-bezout-factoring.json

.PHONY: verify-joint-retraction-module
verify-joint-retraction-module:
	python3 tools/check_joint_hot_retraction_module.py


.PHONY: verify-filtered-pascal-exponents
verify-filtered-pascal-exponents:
	python3 tools/check_filtered_pascal_exponents.py --output /tmp/filtered-pascal-exponents.json


.PHONY: verify-prime-bezout-counting
verify-prime-bezout-counting:
	python3 tools/check_prime_bezout_counting.py --output /tmp/prime-bezout-counting.json

.PHONY: verify-cube-optimized-type
verify-cube-optimized-type:
	python3 tools/certify_cube_optimized_type.py


.PHONY: verify-one-prime-bezout-recovery
verify-one-prime-bezout-recovery:
	python3 tools/check_one_prime_bezout_recovery.py --output /tmp/one-prime-bezout-recovery.json


.PHONY: verify-filtered-unit-crossings
verify-filtered-unit-crossings:
	python3 tools/check_filtered_unit_crossings.py --output /tmp/filtered-unit-crossings.json


.PHONY: verify-recursive-symmetric-exterior-blocks
verify-recursive-symmetric-exterior-blocks:
	python3 tools/check_recursive_symmetric_exterior_blocks.py --large --output /tmp/recursive-symmetric-exterior-blocks.json


.PHONY: verify-bezout-count-values
verify-bezout-count-values:
	python3 tools/check_bezout_count_values.py --output /tmp/bezout-count-values.json

.PHONY: verify-orthogonal-differential-pairs
verify-orthogonal-differential-pairs:
	python3 tools/check_orthogonal_differential_pairs_f4.py


.PHONY: verify-linear-pascal-matroid
verify-linear-pascal-matroid:
	python3 tools/check_linear_pascal_matroid.py --output /tmp/linear-pascal-matroid-energies.json


.PHONY: verify-permuted-pascal-kernels
verify-permuted-pascal-kernels:
	python3 tools/check_permuted_pascal_kernels.py --output /tmp/permuted-pascal-zero-spectra.json

.PHONY: verify-joint-hot-jet-rankone
verify-joint-hot-jet-rankone:
	python3 tools/check_joint_hot_jet_rankone.py


.PHONY: verify-permuted-pascal-dominance
verify-permuted-pascal-dominance:
	python3 tools/check_permuted_pascal_dominance.py --output /tmp/permuted-pascal-cone-dominance.json

.PHONY: verify-permuted-pascal-free-gauges
verify-permuted-pascal-free-gauges:
	python3 tools/check_permuted_pascal_free_gauges.py --output /tmp/permuted-pascal-free-gauges.json

.PHONY: verify-pascal-letter-gauges
verify-pascal-letter-gauges:
	python3 tools/check_pascal_letter_gauges.py --output /tmp/pascal-letter-gauge-controls.json

.PHONY: verify-pascal-coordinate-sectors
verify-pascal-coordinate-sectors:
	python3 tools/check_pascal_coordinate_sectors.py --output /tmp/pascal-coordinate-sectors-controls.json

.PHONY: verify-f4-pascal-flags
verify-f4-pascal-flags:
	python3 tools/check_f4_pascal_flags.py --output /tmp/f4-pascal-flag-controls.json

.PHONY: verify-two-place-determinant verify-two-place-f4-census
verify-two-place-determinant:
	python3 tools/check_two_place_determinant.py --output /tmp/two-place-determinant.json

verify-two-place-f4-census:
	$(CXX) -O3 -std=c++17 tools/scan_two_place_f4_pencil.cpp -o /tmp/integer-factoring-two-place-f4-census-check
	python3 tools/check_two_place_f4_census.py --binary /tmp/integer-factoring-two-place-f4-census-check --output /tmp/two-place-f4-census-controls.json

verify-frontier: verify-two-place-determinant

.PHONY: verify-graded-zero-lattice verify-two-place-f4-census-complete
verify-graded-zero-lattice:
	python3 tools/check_graded_zero_lattice.py --output /tmp/graded-zero-lattice.json

verify-two-place-f4-census-complete:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_two_place_census_complete.py --output /tmp/two-place-census-complete-controls.json

verify-frontier: verify-graded-zero-lattice

.PHONY: verify-fused-f4-envelope verify-fused-f4-finite-controls
verify-fused-f4-envelope:
	$(CXX) -O3 -std=c++17 tools/search_fused_f4_envelope.cpp -o /tmp/integer-factoring-fused-f4-search-check
	$(CXX) -O3 -std=c++17 tools/scan_fused_f4_subspace.cpp -o /tmp/integer-factoring-fused-f4-scan-check
	python3 tools/check_fused_f4_envelope.py --binary /tmp/integer-factoring-fused-f4-search-check --scanner /tmp/integer-factoring-fused-f4-scan-check --output /tmp/fused-f4-envelope-controls.json

verify-fused-f4-finite-controls:
	python3 tools/check_fused_f4_finite_controls.py --output /tmp/fused-f4-finite-controls.json

.PHONY: verify-natural-five-letter-floor
verify-natural-five-letter-floor:
	python3 tools/check_natural_five_letter_floor.py --replay evidence/natural-five-letter-floor.json.xz

verify-frontier: verify-natural-five-letter-floor

.PHONY: verify-graded-adjugate-caps
verify-graded-adjugate-caps:
	python3 tools/check_graded_adjugate_caps.py --output /tmp/graded-adjugate-caps.json

verify-frontier: verify-graded-adjugate-caps

.PHONY: verify-windowed-smith-frames
verify-windowed-smith-frames:
	python3 tools/check_windowed_binary_smith.py --output /tmp/windowed-binary-smith.json

.PHONY: verify-fused-f4-zero-census
verify-fused-f4-zero-census:
	python3 tools/check_fused_f4_zero_census.py --output /tmp/fused-f4-zero-census-controls.json

.PHONY: verify-windowed-f4-affine
verify-windowed-f4-affine:
	python3 tools/check_windowed_f4_affine.py --output /tmp/windowed-f4-affine.json

.PHONY: verify-subalphabet-adjugate-caps
verify-subalphabet-adjugate-caps:
	python3 tools/check_subalphabet_adjugate_caps.py --output /tmp/subalphabet-adjugate-caps.json

verify-frontier: verify-subalphabet-adjugate-caps

.PHONY: verify-subalphabet-bare-entropy
verify-subalphabet-bare-entropy:
	python3 tools/check_subalphabet_bare_entropy.py --output /tmp/subalphabet-bare-entropy.json

verify-frontier: verify-subalphabet-bare-entropy

.PHONY: verify-f4-pascal-coordinate-sectors
verify-f4-pascal-coordinate-sectors:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_f4_pascal_coordinate_sectors.py --output /tmp/f4-pascal-coordinate-sectors.json

.PHONY: verify-rectangular-adjugate-jets
verify-rectangular-adjugate-jets:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_rectangular_adjugate_jets.py --output /tmp/rectangular-adjugate-jets.json

verify-frontier: verify-rectangular-adjugate-jets

.PHONY: verify-scaled-asymmetric-pressure
verify-scaled-asymmetric-pressure:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_scaled_asymmetric_pressure.py --output /tmp/scaled-asymmetric-pressure.json

.PHONY: verify-pencil-uniform-transport
verify-pencil-uniform-transport:
	python3 tools/check_pencil_uniform_transport.py --output /tmp/pencil-uniform-transport.json

verify-frontier: verify-pencil-uniform-transport

.PHONY: verify-polynomial-crossing-axis
verify-polynomial-crossing-axis:
	python3 tools/check_polynomial_crossing_axis.py --recorded-principal

verify-frontier: verify-polynomial-crossing-axis

.PHONY: verify-fused-f4-obstruction
verify-fused-f4-obstruction:
	python3 tools/check_fused_f4_decomposition_obstruction.py

.PHONY: verify-smith-chambers-f4
verify-smith-chambers-f4:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_smith_chambers_f4.py --output /tmp/smith-chambers-f4.json

verify-frontier: verify-smith-chambers-f4

.PHONY: verify-coherent-f4-partial
verify-coherent-f4-partial:
	python3 tools/check_coherent_hot_f4_partial.py

.PHONY: verify-orthogonal-pair-core-floor
verify-orthogonal-pair-core-floor:
	python3 tools/check_orthogonal_pair_core_floor.py

.PHONY: verify-circuit-box-products
verify-circuit-box-products:
	python3 tools/check_circuit_box_products.py

.PHONY: verify-batch-smooth
verify-batch-smooth:
	python3 tools/batch_smooth_check.py

.PHONY: verify-circuit-box-normalization verify-circuit-box-bezout verify-circuit-box-averages verify-circuit-box-bezout-averages verify-circuit-box-pole-trace
verify-circuit-box-normalization: verify-circuit-box-bezout verify-circuit-box-averages verify-circuit-box-bezout-averages verify-circuit-box-pole-trace

verify-circuit-box-bezout:
	python3 tools/check_circuit_box_bezout.py

verify-circuit-box-averages:
	python3 tools/check_circuit_box_averages.py

verify-circuit-box-bezout-averages:
	python3 tools/check_circuit_box_bezout_averages.py

verify-circuit-box-pole-trace:
	python3 tools/check_circuit_box_pole_trace.py

.PHONY: verify-tensor-orthogonal verify-tensor-orthogonal-scaling verify-tensor-orthogonal-capacity verify-tensor-f4-catalogues verify-tensor-orthogonal-comparisons

verify-tensor-orthogonal: verify-tensor-orthogonal-scaling verify-tensor-orthogonal-capacity verify-tensor-f4-catalogues verify-tensor-orthogonal-comparisons

verify-tensor-orthogonal-scaling:
	python3 tools/check_tensor_orthogonal_scaling.py --output /tmp/integer-factoring-tensor-scaling-check.json

verify-tensor-orthogonal-capacity:
	python3 tools/check_tensor_orthogonal_capacity.py --output /tmp/integer-factoring-tensor-capacity-check.json

verify-tensor-f4-catalogues:
	python3 tools/check_tensor_f4_rank_two.py --output /tmp/integer-factoring-tensor-rank-two-check.json
	python3 tools/tensor_f4_rank_three_catalogue.py --output /tmp/integer-factoring-tensor-rank-three-check.json
	python3 tools/check_tensor_f4_hyperoval_obstruction.py --output /tmp/integer-factoring-tensor-hyperoval-check.json

verify-tensor-orthogonal-comparisons:
	python3 -c 'import pathlib, shutil, subprocess, tempfile; d=tempfile.TemporaryDirectory(prefix="factoring-tensor-"); p=pathlib.Path(d.name)/"certificate.json"; shutil.copyfile("evidence/tensor-orthogonal-pressure-dominance.json",p); subprocess.run(["python3","tools/check_tensor_orthogonal_pressure_dominance.py","--output",str(p)],check=True)'
	python3 tools/check_tensor_f4_five_pair_dominance.py --output /tmp/integer-factoring-tensor-five-pair-check.json
	python3 tools/check_tensor_coupled_quadratic_dominance.py --output /tmp/integer-factoring-tensor-coupled-check.json

.PHONY: verify-oriented-coprime-basis

verify-oriented-coprime-basis:
	python3 tools/check_oriented_coprime_basis.py

.PHONY: benchmark-batch-cofactor-policies

benchmark-batch-cofactor-policies:
	python3 tools/batch_smooth_policy_check.py --output /tmp/integer-factoring-batch-cofactor-policies.json

.PHONY: verify-paired-mds-low-core
verify-paired-mds-low-core:
	python3 tools/check_paired_mds_low_core_floor.py

.PHONY: verify-tensor-f16

verify-tensor-orthogonal: verify-tensor-f16

verify-tensor-f16:
	python3 tools/check_tensor_hot_compounds_f16.py --output /tmp/integer-factoring-tensor-f16-check.json
	python3 tools/check_tensor_coupled_quadratic_core_floor.py --recompute-profiles --output /tmp/integer-factoring-tensor-field-floor-check.json

.PHONY: verify-oriented-refinement verify-oriented-cofactor-prefixes verify-signed-residual-graph

verify-oriented-refinement: verify-oriented-coprime-basis verify-oriented-cofactor-prefixes verify-signed-residual-graph

verify-oriented-cofactor-prefixes:
	python3 tools/check_oriented_cofactor_prefixes.py

verify-signed-residual-graph:
	python3 tools/check_signed_residual_graph.py --prefix evidence/oriented-cofactor-prefixes.json

.PHONY: verify-circuit-box-inverse-windows

verify-circuit-box-normalization: verify-circuit-box-inverse-windows

verify-circuit-box-inverse-windows:
	python3 tools/check_circuit_box_inverse_windows.py --output /tmp/integer-factoring-inverse-window-check.json

.PHONY: verify-incremental-oriented-basis

verify-oriented-refinement: verify-incremental-oriented-basis

verify-incremental-oriented-basis:
	python3 tools/check_incremental_oriented_coprime_basis.py

.PHONY: verify-tensor-f16-fast

verify-tensor-f16: verify-tensor-f16-fast

verify-tensor-f16-fast:
	python3 tools/check_tensor_hot_compounds_f16_fast.py --output /tmp/integer-factoring-tensor-f16-fast-check.json

.PHONY: verify-sparse-oriented-power

verify-oriented-refinement: verify-sparse-oriented-power

verify-sparse-oriented-power:
	python3 tools/check_sparse_oriented_power_search.py --replay evidence/sparse-oriented-power-64.json.xz
	python3 tools/check_sparse_oriented_power_search.py --replay evidence/sparse-oriented-power-80.json.xz

.PHONY: verify-binary-extension-smith
verify-binary-extension-smith:
	python3 tools/check_binary_extension_smith.py

.PHONY: verify-circuit-box-kloosterman-blocks

verify-circuit-box-normalization: verify-circuit-box-kloosterman-blocks

verify-circuit-box-kloosterman-blocks:
	python3 tools/check_circuit_box_kloosterman_blocks.py --output /tmp/integer-factoring-kloosterman-blocks.json

.PHONY: verify-tensor-cubic-graph

verify-tensor-f16: verify-tensor-cubic-graph

verify-tensor-cubic-graph:
	python3 tools/check_tensor_hot_compounds_f16_mixed.py --output /tmp/integer-factoring-tensor-cubic-graph.json

.PHONY: verify-sparse-unit-presentation

verify-oriented-refinement: verify-sparse-unit-presentation

verify-sparse-unit-presentation:
	python3 tools/check_sparse_unit_presentation.py

.PHONY: verify-quotient-update-oriented-basis

verify-oriented-refinement: verify-quotient-update-oriented-basis

verify-quotient-update-oriented-basis:
	python3 tools/check_quotient_update_oriented_basis.py --repetitions 1

.PHONY: verify-circuit-box-ratio-lookahead

verify-circuit-box-normalization: verify-circuit-box-ratio-lookahead

verify-circuit-box-ratio-lookahead:
	python3 tools/check_circuit_box_ratio_lookahead.py --output /tmp/integer-factoring-ratio-lookahead.json

.PHONY: verify-unit-presentation-trace

verify-oriented-refinement: verify-unit-presentation-trace

verify-unit-presentation-trace:
	python3 tools/check_unit_presentation_trace.py

.PHONY: verify-optimized-oriented-power

verify-oriented-refinement: verify-optimized-oriented-power

verify-optimized-oriented-power:
	python3 tools/check_optimized_oriented_power_search.py --replay evidence/optimized-oriented-power-search.json
	python3 tools/check_optimized_oriented_fresh_run.py --replay evidence/optimized-oriented-power-96.json.xz

.PHONY: verify-tensor-literal-hot

verify-tensor-f16: verify-tensor-literal-hot

verify-tensor-literal-hot:
	python3 tools/check_tensor_literal_hot_pfaffians_f16.py --output /tmp/integer-factoring-tensor-literal-hot.json

.PHONY: verify-circuit-box-ratio-lattice

verify-circuit-box-normalization: verify-circuit-box-ratio-lattice

verify-circuit-box-ratio-lattice:
	python3 tools/check_circuit_box_ratio_lattice.py --output /tmp/integer-factoring-ratio-lattice.json

.PHONY: verify-primitive-integer-row-basis

verify-oriented-refinement: verify-primitive-integer-row-basis

verify-primitive-integer-row-basis:
	python3 tools/check_primitive_integer_row_basis.py --repetitions 1

.PHONY: verify-batched-coprime-bucket-basis
verify-batched-coprime-bucket-basis:
	python3 tools/check_batched_coprime_bucket_basis.py --output /tmp/batched-coprime-bucket-basis-verify.json

verify-oriented-coprime-basis: verify-batched-coprime-bucket-basis

.PHONY: verify-compact-oriented-power-search
verify-compact-oriented-power-search:
	python3 tools/check_compact_oriented_power_search.py --replay evidence/compact-oriented-power-search.json

verify-oriented-coprime-basis: verify-compact-oriented-power-search

.PHONY: verify-tensor-literal-height15-gate verify-tensor-literal-rank-one-cells
verify-tensor-literal-height15-gate:
	python3 tools/check_tensor_literal_height15_pressure_gate.py --recompute-profiles --output /tmp/tensor-literal-height15-gate-verify.json

verify-tensor-literal-rank-one-cells:
	python3 tools/check_tensor_literal_rank_one_cells.py --output /tmp/tensor-literal-rank-one-cells-verify.json

.PHONY: verify-small-pivot-integer-rows verify-prepared-oriented-power-search
verify-small-pivot-integer-rows:
	python3 tools/check_small_pivot_integer_row_basis.py --controls-only --output /tmp/small-pivot-integer-row-verify.json

verify-prepared-oriented-power-search:
	python3 tools/check_prepared_oriented_power_search.py --replay evidence/prepared-oriented-power-search.json

verify-oriented-coprime-basis: verify-small-pivot-integer-rows verify-prepared-oriented-power-search

.PHONY: verify-oriented-basis-linear-charges
verify-oriented-basis-linear-charges:
	python3 tools/check_oriented_norm_basis_lift.py --output /tmp/oriented-norm-basis-lift-verify.json

verify-oriented-coprime-basis: verify-oriented-basis-linear-charges

.PHONY: verify-three-row-hecke verify-three-row-hecke-full
verify-three-row-hecke:
	python3 tools/check_three_row_hecke_crossing.py --alphabet 5
verify-three-row-hecke-full:
	python3 tools/check_three_row_hecke_crossing.py --alphabet 9

.PHONY: verify-oriented-power-112
verify-oriented-power-112:
	python3 tools/check_optimized_oriented_fresh_run.py --replay evidence/optimized-oriented-power-112.json.xz

.PHONY: verify-circuit-box-ratio-small-roots
verify-circuit-box-ratio-small-roots:
	python3 tools/check_circuit_box_ratio_small_roots.py --output /tmp/circuit-box-ratio-small-roots-verify.json

.PHONY: verify-batch-smooth-fresh-192
verify-batch-smooth-fresh-192:
	python3 tools/assemble_factoring_certificate.py evidence/batch-smooth-fresh-192-selective.json.xz.manifest.json
	python3 tools/check_fresh_prime_cofactor_control.py --bits 192 --replay evidence/batch-smooth-fresh-192-selective.json.xz

.PHONY: verify-classical-oriented-basis
verify-classical-oriented-basis:
	python3 tools/check_classical_coprime_basis.py

.PHONY: verify-tensor-literal-full-tail
verify-tensor-literal-full-tail:
	python3 tools/check_tensor_literal_full_tail_core.py --output /tmp/tensor-full-tail-replay.json
	python3 tools/check_tensor_literal_full_tail_rank_one.py --output /tmp/tensor-full-tail-rank-one-replay.json

.PHONY: verify-sparse-projected-graph
verify-sparse-projected-graph:
	python3 tools/check_sparse_projected_signed_graph.py
	python3 tools/check_projected_oriented_power_search.py --replay evidence/projected-oriented-power-search.json

.PHONY: verify-two-pass-oriented-basis
verify-two-pass-oriented-basis:
	python3 tools/check_circuit_box_oriented_batch.py --output /tmp/two-pass-oriented-replay.json

.PHONY: verify-tensor-literal-moving-frames
verify-tensor-literal-moving-frames:
	python3 tools/check_tensor_literal_moving_frames.py --output /tmp/tensor-moving-frame-replay.json

.PHONY: verify-tensor-affine-full-tail
verify-tensor-affine-full-tail:
	python3 tools/check_tensor_literal_full_tail_affine_frames.py --output /tmp/tensor-affine-full-tail-replay.json

.PHONY: verify-equivalent-norm-stencils
verify-equivalent-norm-stencils:
	python3 tools/check_equivalent_norm_stencils.py --output /tmp/equivalent-norm-stencils-replay.json

.PHONY: verify-compact-power-witnesses
verify-compact-power-witnesses:
	python3 tools/check_compact_witness_oriented_power_search.py --replay evidence/compact-witness-oriented-power-search.json

.PHONY: verify-compact-oriented-basis
verify-compact-oriented-basis:
	python3 tools/check_compact_oriented_basis.py --output /tmp/compact-oriented-basis-replay.json

.PHONY: verify-tensor-affine-clifford
verify-tensor-affine-clifford:
	python3 tools/check_tensor_affine_clifford_profiles.py --output /tmp/tensor-affine-clifford-replay.json

.PHONY: check-stencil-overlap
check-stencil-overlap:
	python3 tools/check_circuit_box_stencil_overlap.py

.PHONY: check-direct-unit-presentation check-direct-unit-search
check-direct-unit-presentation:
	python3 tools/check_direct_unit_presentation.py --replay evidence/direct-unit-presentation.json
check-direct-unit-search:
	python3 tools/check_direct_unit_oriented_power_search.py --replay evidence/direct-unit-oriented-power-search.json.xz

.PHONY: check-centered-windows
check-centered-windows:
	python3 tools/check_circuit_box_centered_windows.py

.PHONY: check-cached-lll
check-cached-lll:
	python3 tools/check_circuit_box_cached_lll.py

.PHONY: check-singular-affine-clifford
check-singular-affine-clifford:
	python3 tools/check_tensor_singular_affine_clifford.py

.PHONY: verify-orthogonal-joint-global
verify-orthogonal-joint-global:
	python3 tools/check_orthogonal_joint_global_rigidity.py
	python3 tools/check_orthogonal_joint_finite_profiles.py

.PHONY: verify-orthogonal-joint-symbolic
verify-orthogonal-joint-symbolic:
	python3 tools/check_orthogonal_joint_symbolic_profiles.py

.PHONY: verify-height36-core-floor verify-height36-core-replay
verify-height36-core-floor:
	python3 tools/check_height36_core_floor.py
verify-height36-core-replay:
	python3 tools/check_height36_core_floor.py --replay-core

.PHONY: verify-coherent-f8
verify-coherent-f8:
	python3 tools/check_coherent_hot_f8.py

.PHONY: verify-f8-hyperoval verify-f8-hyperoval-full
verify-f8-hyperoval:
	python3 tools/replay_f8_hyperoval_q96.py
	python3 tools/certify_f8_hyperoval_q96.py
verify-f8-hyperoval-full:
	python3 tools/replay_f8_hyperoval_q96.py --full
	python3 tools/certify_f8_hyperoval_q96.py

.PHONY: verify-hecke-all-places verify-hecke-all-places-full
verify-hecke-all-places:
	python3 tools/check_hecke_all_algebraic_places.py --alphabet 5
verify-hecke-all-places-full:
	python3 tools/check_hecke_all_algebraic_places.py --alphabet 9

.PHONY: verify-pencil-endpoint-moments
verify-pencil-endpoint-moments:
	python3 tools/check_pencil_endpoint_moments.py

.PHONY: verify-pencil-principal-alphabets
verify-pencil-principal-alphabets:
	python3 tools/check_pencil_principal_alphabets.py

.PHONY: verify-pencil-principal-screen
verify-pencil-principal-screen:
	python3 tools/check_pencil_principal_screen.py

.PHONY: verify-pascal-coordinate-bundle
verify-pascal-coordinate-bundle:
	python3 tools/check_pascal_coordinate_bundle.py

.PHONY: verify-rectangular-pascal-coordinates
verify-rectangular-pascal-coordinates:
	python3 tools/check_rectangular_pascal_sectors.py

.PHONY: verify-jensen-core-floors
verify-jensen-core-floors:
	python3 tools/check_jensen_tail_core_floor.py

.PHONY: verify-coherent-inverse-f4 verify-coherent-gram-f4
verify-coherent-inverse-f4:
	python3 tools/check_coherent_inverse_systems_f4.py

verify-coherent-gram-f4:
	python3 tools/check_coherent_gram_f4.py

.PHONY: verify-coherent-divided-powers
verify-coherent-divided-powers:
	python3 tools/check_coherent_divided_powers.py

.PHONY: verify-coherent-divided-probes
verify-coherent-divided-probes:
	python3 tools/check_coherent_divided_probes.py

.PHONY: verify-hecke-q5-determinants
verify-hecke-q5-determinants:
	python3 tools/check_hecke_q5_determinants.py

.PHONY: verify-f8-retraction-refinement verify-f8-retraction-refinement-full verify-f8-retraction-census-full
verify-f8-retraction-refinement:
	python3 tools/check_f8_hyperoval_retraction_entries.py evidence/f8-hyperoval-retraction-678.json --output evidence/f8-hyperoval-retraction-678-entry-controls.json
	python3 tools/check_f8_hyperoval_retraction_entries.py evidence/f8-hyperoval-retraction-789.json --output evidence/f8-hyperoval-retraction-789-entry-controls.json
	python3 tools/check_f8_hyperoval_retraction_census.py
	python3 tools/certify_f8_hyperoval_refinement.py
verify-f8-retraction-refinement-full:
	python3 tools/replay_f8_hyperoval_refinement.py
	$(MAKE) verify-f8-retraction-refinement
verify-f8-retraction-census-full:
	python3 tools/search_f8_hyperoval_retractions.py
	python3 tools/check_f8_hyperoval_retraction_census.py

.PHONY: verify-mixed-graph-retractions
verify-mixed-graph-retractions:
	python3 tools/check_mixed_graph_retractions.py

.PHONY: verify-asymmetric-profile-bundle
verify-asymmetric-profile-bundle:
	python3 tools/check_asymmetric_profile_bundle.py

.PHONY: verify-q96-asymmetric-scale
verify-q96-asymmetric-scale:
	python3 tools/check_q96_asymmetric_scales.py --output /tmp/q96-asymmetric-scale-check.json

.PHONY: verify-q12-reflection-graphs verify-q12-reflection-graphs-independent
verify-q12-reflection-graphs:
	python3 tools/check_q12_reflection_graphs.py
verify-q12-reflection-graphs-independent:
	python3 tools/check_q12_reflection_graphs.py --full-independent

.PHONY: verify-normalized-support-gap
verify-normalized-support-gap:
	python3 tools/check_normalized_support_gap.py

.PHONY: verify-rectangular-coordinate-support
verify-rectangular-coordinate-support:
	python3 tools/check_rectangular_coordinate_support.py

.PHONY: verify-regularized-transpose-words
verify-regularized-transpose-words:
	python3 tools/check_regularized_transpose_words.py

.PHONY: verify-cyclic-pauli-crossings
verify-cyclic-pauli-crossings:
	python3 tools/check_cyclic_pauli_crossings.py
.PHONY: verify-joint-frobenius
verify-joint-frobenius:
	python3 tools/check_joint_frobenius_carries.py
	python3 tools/check_joint_geometric_box_network.py

.PHONY: verify-trace-constraints
verify-trace-constraints:
	python3 tools/check_trace_constraint_contraction.py
.PHONY: verify-fused-q12-endpoints
verify-fused-q12-endpoints:
	python3 tools/check_fused_q12_endpoint_block_cover.py

.PHONY: verify-frobenius-prefix
verify-frobenius-prefix:
	python3 tools/check_frobenius_prefix_contraction.py

.PHONY: verify-square-zero-traces
verify-square-zero-traces:
	python3 tools/check_square_zero_trace_contraction.py

.PHONY: verify-collective-css-pauli
verify-collective-css-pauli:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_collective_css_pauli.py

.PHONY: verify-integral-spin-annihilators
verify-integral-spin-annihilators:
	python3 tools/check_integral_spin_annihilators.py

verify-constructions: verify-integral-spin-annihilators

.PHONY: verify-high-axis-cayley
verify-high-axis-cayley:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_high_axis_cayley.py

.PHONY: verify-high-axis-plucker
verify-high-axis-plucker:
	OPENBLAS_NUM_THREADS=1 python3 tools/check_high_axis_plucker.py

.PHONY: verify-cold-gap
verify-cold-gap:
	python3 tools/check_cold_gap_bound.py

.PHONY: verify-uniform-axis-window
verify-uniform-axis-window:
	python3 tools/check_uniform_axis_window.py

.PHONY: verify-approximate-strictification
verify-approximate-strictification:
	python3 tools/check_approximate_strictification.py

.PHONY: verify-elliptic-word-costs
verify-elliptic-word-costs:
	python3 tools/check_elliptic_word_costs.py

verify-arithmetic: verify-elliptic-word-costs

.PHONY: verify-teichmuller-carries
verify-teichmuller-carries:
	python3 tools/check_teichmuller_cyclic_carries.py

verify-arithmetic: verify-teichmuller-carries

.PHONY: verify-teichmuller-inverses
verify-teichmuller-inverses:
	python3 tools/check_teichmuller_unit_inverses.py

verify-arithmetic: verify-teichmuller-inverses

.PHONY: verify-teichmuller-packed-systems
verify-teichmuller-packed-systems:
	python3 tools/check_teichmuller_packed_systems.py

verify-arithmetic: verify-teichmuller-packed-systems

.PHONY: verify-multiaxis-face-rank
verify-multiaxis-face-rank:
	python3 tools/check_multiaxis_face_rank_floor.py

.PHONY: verify-precision-indicators
verify-precision-indicators:
	python3 tools/check_precision_residue_indicators.py

.PHONY: verify-precision-rational-carries
verify-precision-rational-carries:
	python3 tools/check_precision_rational_carries.py

.PHONY: verify-lazy-potential-tail verify-routing-graphs verify-window-graph-data
verify-lazy-potential-tail:
	python3 tools/check_lazy_potential_tail.py

GRAPH_PYTHON ?= python3
verify-routing-graphs:
	$(GRAPH_PYTHON) tools/routing_rank_nogo/check_graph_minrank.py

verify-window-graph-data:
	cd tools/window_lemma/graph/data && sha256sum --check SHA256SUMS
	python3 -m py_compile tools/window_lemma/graph/cayley_small.py tools/window_lemma/graph/summarize.py

.PHONY: verify-multiplier-treewidth-controls verify-multiplier-treewidth-brackets
verify-multiplier-treewidth-controls:
	python3 tools/multiplier_treewidth/check_shift_lemma.py
	python3 tools/multiplier_treewidth/check_cayley_relations.py

verify-multiplier-treewidth-brackets:
	python3 tools/multiplier_treewidth/tw_multipliers.py

.PHONY: verify-witt-trace-words
verify-witt-trace-words:
	python3 tools/check_witt_trace_word_circuits.py

verify-arithmetic: verify-witt-trace-words
.PHONY: verify-divided-power-traces
verify-divided-power-traces:
	python3 tools/check_divided_power_trace_carries.py

verify-arithmetic: verify-divided-power-traces

.PHONY: verify-artin-hasse-traces
verify-artin-hasse-traces:
	python3 tools/check_artin_hasse_trace_carries.py

verify-arithmetic: verify-artin-hasse-traces

.PHONY: verify-normal-basis-and
verify-normal-basis-and:
	python3 tools/check_normal_basis_and_tradeoff.py

.PHONY: verify-fixed-cube-projection
verify-fixed-cube-projection:
	python3 tools/check_fixed_cube_projection_floor.py

.PHONY: verify-block-gaussian-floor
verify-block-gaussian-floor:
	python3 tools/check_block_gaussian_profile_floor.py

.PHONY: verify-trace-selection-weil
verify-trace-selection-weil:
	python3 tools/check_trace_selection_weil.py

.PHONY: verify-mixed-direction-carries
verify-mixed-direction-carries:
	python3 tools/check_mixed_direction_carries.py

.PHONY: verify-quartic-interactions
verify-quartic-interactions:
	python3 tools/check_quartic_interaction_panel.py

.PHONY: verify-quartic-weighted-run
verify-quartic-weighted-run:
	python3 tools/check_quartic_weighted_run.py

.PHONY: verify-cube-favorable-reweighting
verify-cube-favorable-reweighting:
	python3 tools/check_cube_favorable_reweighting.py

.PHONY: verify-recombination-before-orientation
verify-recombination-before-orientation:
	PYTHONPATH=tools python3 tools/check_recombination_before_orientation.py

.PHONY: verify-cube-global-functional-profile
verify-cube-global-functional-profile:
	python3 tools/check_cube_global_functional_profile.py

.PHONY: verify-determinant-inverse-codes
verify-determinant-inverse-codes:
	python3 tools/check_determinant_inverse_codes.py

.PHONY: verify-determinant-inverse-profiles
verify-determinant-inverse-profiles:
	python3 tools/check_determinant_inverse_profiles.py

.PHONY: verify-pfaffian-derivative-codes
verify-pfaffian-derivative-codes:
	python3 tools/check_pfaffian_derivative_codes.py

.PHONY: verify-window-orbit-duality verify-window-duality-analysis
verify-window-orbit-duality:
	python3 tools/check_window_orbit_duality.py

verify-window-duality-analysis:
	python3 tools/check_window_duality_analysis.py

.PHONY: verify-trace-anf-counterexample
verify-trace-anf-counterexample:
	python3 tools/check_trace_anf_counterexample.py

verify-arithmetic: verify-trace-anf-counterexample

.PHONY: verify-clustered-power-relations
verify-clustered-power-relations:
	python3 tools/check_clustered_power_relations.py

verify-arithmetic: verify-clustered-power-relations

.PHONY: verify-power-character-lifts
verify-power-character-lifts: verify-clustered-power-relations
	python3 tools/check_power_character_lifts.py

verify-arithmetic: verify-power-character-lifts

.PHONY: verify-relation-change-events
verify-relation-change-events: verify-clustered-power-relations
	python3 tools/check_relation_change_events.py

verify-arithmetic: verify-relation-change-events
