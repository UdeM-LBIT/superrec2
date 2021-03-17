test:
	python -m unittest discover --buffer

lint:
	mypy .

.PHONY: test lint
