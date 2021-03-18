test:
	python -m unittest discover --buffer

lint:
	pylint --ignore-patterns "test_" *.py superrec2
	mypy .

format:
	black --line-length 80 --check --diff .

format-fix:
	black --line-length 80 .

.PHONY: test lint format format-fix
