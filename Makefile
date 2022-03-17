define USAGE
superrec2 maintenance tasks.
Please refer to the README for general usage guidance.

Available commands:

    test         Run the test suite.
    fixtures     Regenerate the expected test outputs for the
                 rendering module and show differences.
    lint         Check for common errors and correct typing.
    format       Check that the source code follows formatting rules.
    format-fix   Automatically fix formatting errors.
endef
export USAGE

help:
	@echo "$$USAGE"

test:
	python -m unittest discover --buffer --verbose

.ONESHELL:
SHELL=/usr/bin/bash
fixtures:
	for input in superrec2/render/fixtures/*/input.json; do
	    for orient in horizontal vertical; do
	    	echo "draw $$input $$orient"
		./draw.py --orientation "$$orient" < "$$input" \
	            > "$${input%input.json}output-$$orient.tex"
	    done
	done
	git diff superrec2/render/fixtures

lint:
	pylint --ignore-patterns "test_" *.py superrec2
	mypy .

format:
	black --line-length 80 --check --diff .

format-fix:
	black --line-length 80 .

.PHONY: help test lint format format-fix
