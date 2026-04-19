.PHONY: help install stage1 task1.1 task1.2 task1.3 task1.4 task1.5 clean

PY := python
CONFIG := configs/stage1.yaml

help:
	@echo "Targets:"
	@echo "  install        Install package + dev deps + pre-commit hooks"
	@echo "  stage1         Run all Stage 1 tasks (1.1 -> 1.5)"
	@echo "  task1.1        Species inventory"
	@echo "  task1.2        Threshold sensitivity"
	@echo "  task1.3        Environmental coverage"
	@echo "  task1.4        Basin overlap matrix"
	@echo "  task1.5        Invasion contamination audit"
	@echo "  clean          Remove generated outputs in data/processed and reports/stage1"

install:
	pip install -e ".[dev]"
	pre-commit install

stage1: task1.1 task1.2 task1.3 task1.4 task1.5
	@echo ""
	@echo "Stage 1 outputs in data/processed/ and reports/stage1/"
	@echo "Next: write reports/stage1/etapa1_raport.md and convene Decision Point P1."

task1.1:
	$(PY) scripts/01_species_inventory.py --config $(CONFIG)

task1.2: task1.1
	$(PY) scripts/02_threshold_sensitivity.py --config $(CONFIG)

task1.3: task1.1
	$(PY) scripts/03_env_coverage.py --config $(CONFIG)

task1.4: task1.1
	$(PY) scripts/04_basin_overlap.py --config $(CONFIG)

task1.5: task1.1
	$(PY) scripts/05_invasion_audit.py --config $(CONFIG)

clean:
	rm -f data/processed/*.csv data/processed/*.geojson
	rm -f reports/stage1/*.png
