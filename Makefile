.PHONY: help demo test lint fix fetch analyse clean
.DEFAULT_GOAL := help

PY ?= python
PER_YEAR ?= 400

help:            ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	  | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-9s\033[0m %s\n", $$1, $$2}'

demo:            ## Show the reference check and the finding (no network)
	$(PY) demo.py

test:            ## Run the test suite (no network)
	$(PY) -m pytest -q

lint:            ## Check formatting and lint
	$(PY) -m ruff check src tests scripts demo.py
	$(PY) -m ruff format --check src tests scripts demo.py

fix:             ## Apply formatting and autofixes
	$(PY) -m ruff format src tests scripts demo.py
	$(PY) -m ruff check --fix src tests scripts demo.py

fetch:           ## Download genomes from GenBank (PER_YEAR=400 by default)
	$(PY) scripts/fetch.py --per-year $(PER_YEAR)

analyse:         ## Score every assay against the downloaded genomes
	$(PY) scripts/analyse.py

clean:           ## Remove caches (NOT the downloaded sequences)
	rm -rf .pytest_cache .ruff_cache src/assaydrift/__pycache__ tests/__pycache__
