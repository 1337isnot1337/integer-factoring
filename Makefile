.PHONY: verify verify-docs verify-frontier verify-post-half verify-schur-fermion verify-schur-central verify-schur-gaussian verify-schur-sparse verify-cancellation-modules verify-minimal-modules verify-exact-modules verify-cancellation-crossing verify-block-smith verify-linear-smith verify-type-pressure verify-controlled-smith verify-sparse-holonomy verify-cubic-full-images verify-formal verify-k6 verify-half-appendix help

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

verify-controlled-smith:
	python3 tools/check_controlled_permutation_smith.py
	python3 tools/check_qutrit_weyl_crossing.py

verify-sparse-holonomy:
	python3 tools/check_sparse_holonomy_singletons.py

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
	  'make verify-cancellation-modules - exact finite module and word certificates' \
	  'make verify-minimal-modules - exact boundary reduction and complete word rates' \
	  'make verify-exact-modules - complete rational lattice orbits and pressures' \
	  'make verify-cancellation-crossing - strict module improvement on a complete crossing' \
	  'make verify-block-smith - ordered block profiles and submultiplicativity' \
	  'make verify-linear-smith - uniform degree-one obstruction and binary catalogue' \
	  'make verify-type-pressure - rational type and physical-gauge certificates' \
	  'make verify-controlled-smith - controlled monomial and qutrit crossings' \
	  'make verify-sparse-holonomy - sparse Smith profiles and Gaussian identities' \
	  'make verify-cubic-full-images - complete varying-image cubic spaces and ranks' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
