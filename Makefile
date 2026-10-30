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
	python3 tools/check_modular_central_rank.py
	python3 tools/check_modular_binary_central_rank.py
	python3 tools/check_three_label_span.py
	python3 tools/check_sparse_divisor_nonzeroness.py

verify-formal:
	bash tools/check_formal.sh

verify-k6:
	bash tools/check_k6.sh
