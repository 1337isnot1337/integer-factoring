.PHONY: verify verify-post-half verify-modular verify-formal verify-k6

verify:
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

verify-modular:
	python3 tools/check_modular_pieces.py --full
	python3 tools/check_weighted_prefix.py
	python3 tools/check_fixed_characteristic_transfer.py

verify-post-half:
	python3 tools/check_half_piece_ranks.py --full
	python3 tools/check_half_projection.py
	python3 tools/check_half_restricted_labels.py

verify-formal:
	bash tools/check_formal.sh

verify-k6:
	bash tools/check_k6.sh
