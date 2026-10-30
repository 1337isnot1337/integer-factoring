.PHONY: verify verify-formal verify-k6

verify:
	python3 -m unittest discover -s tests -v
	python3 -m unittest tools.check_crossing -v
	python3 tools/check_coefficient_basis.py
	python3 tools/parameters.py

verify-formal:
	bash tools/check_formal.sh

verify-k6:
	bash tools/check_k6.sh
