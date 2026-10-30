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
	  'make verify-constructions - new crossing and resource-interface certificates' \
	  'make verify-native-crossings - compile and independently audit the native group screen' \
	  'make verify-post-half     - larger optional matrix-rank checks' \
	  'make verify-docs          - check local Markdown file links' \
	  'make verify-formal        - rebuild original Lean lemmas (Lean 4.35.0)' \
	  'make verify-k6            - rebuild the certified 1/30 machine (Lean 4.35.0)' \
	  'make verify-half-appendix - rebuild the certified 1/2 machine (Lean 4.33.1)'
