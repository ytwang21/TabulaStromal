PYTHON ?= python

.PHONY: check test-benchmark test-cellage smoke-cellage figure
check:
	$(PYTHON) scripts/validate_repository.py --verify-checksums

test-benchmark:
	$(PYTHON) -m unittest discover -s workflows/benchmark/tests -p 'test_*.py' -v

test-cellage:
	cd workflows/cellage && $(PYTHON) -m pytest -q

smoke-cellage:
	$(PYTHON) scripts/run_workflow.py cellage --phase smoke-test --output-root outputs/cellage-smoke

figure:
	$(PYTHON) scripts/plot_summary.py
