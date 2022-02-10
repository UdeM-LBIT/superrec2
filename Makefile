
test:
	python -m unittest discover --buffer --verbose

.ONESHELL:
SHELL=/usr/bin/bash
fixtures:
	for input in superrec2/render/fixtures/*/input.json; do
	    for orient in horizontal vertical; do
		./draw.py --orientation "$$orient" < "$$input" \
	            > "$${input%input.json}output-$$orient.tex"
	    done
	done

lint:
	pylint --ignore-patterns "test_" *.py superrec2
	mypy .

format:
	black --line-length 80 --check --diff .

format-fix:
	black --line-length 80 .

.PHONY: test lint format format-fix
