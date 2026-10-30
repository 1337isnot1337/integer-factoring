.PHONY: verify verify-frontier verify-post-half verify-formal verify-k6

verify: verify-frontier
	python3 -m unittest discover -s tests -v
	python3 -m unittest tools.check_crossing -v
	python3 tools/check_coefficient_basis.py
	python3 tools/parameters.py
	python3 tools/check_k6_sharp_rate.py
	python3 tools/check_half_rate.py
	python3 tools/check_half_piece_ranks.py --no-numpy
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-post-half:
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-frontier:
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
