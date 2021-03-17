test:
	python -m unittest discover --buffer

lint:
	mypy .

format:
	black --line-length 80 --check --diff .

format-fix:
	black --line-length 80 .

.PHONY: test lint format format-fix
