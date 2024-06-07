.PHONY: help install lint format test docs download-suggested-structures

help:
	@echo "Available commands:"
	@echo "  install                			Install dependencies using Poetry"
	@echo "  lint                   			Run linters using pre-commit"
	@echo "  format                 			Run formatters using pre-commit"
	@echo "  test                   			Run tests using pytest"
	@echo "  docs                   			Generate documentation using pdoc"
	@echo "  download-structures				Download suggested structures from the web"

install:
	poetry install

lint:
	pre-commit run ruff --all-files

format:
	pre-commit run ruff-format --all-files

test:
	poetry run pytest

docs:
	poetry run lazydocs --overview-file="README.md" supra_molecular_simulator

download-structures:
	poetry run download-suggested-structures

.DEFAULT_GOAL := help