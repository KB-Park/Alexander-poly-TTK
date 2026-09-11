PYTHON ?= python3

.PHONY: test checks smoke-volume

test:
	$(PYTHON) -m unittest discover -s code -p 'test_*.py'

checks: test
	$(PYTHON) code/check_alexander_degree.py
	$(PYTHON) code/check_non_l_space_family.py --k-max 100

smoke-volume:
	$(PYTHON) code/generalized_ttk_volume_experiment.py \
		--p-max 8 --s-max 2 --region all --volume-scope collisions \
		--isometry-signatures --workers 2 \
		--output-dir output/smoke-volume
