.PHONY: help install lint format pytest mypy mypy-types docs download-structures

help:
	@echo "Available commands:"
	@echo "  install                    Install dependencies using Poetry"
	@echo "  lint                       Run linters using pre-commit"
	@echo "  format                     Run formatters using pre-commit"
	@echo "  pytest                     Run tests using pytest"
	@echo "  mypy                       Run type-checking using mypy"
	@echo "  mypy-types                 Install missing types using mypy"
	@echo "  docs                       Generate documentation using lazydocs"
	@echo "  download-structures        Run supra_molecular_simulator.download:download_suggested_structures"

install:
	poetry install

lint:
	pre-commit run ruff --all-files

format:
	pre-commit run ruff-format --all-files

pytest:
	poetry run pytest

mypy:
	poetry run mypy --ignore-missing-imports supra_molecular_simulator

mypy-types:
	poetry run mypy --install-types

docs:
	rm -rf docs
	poetry run lazydocs --remove-package-prefix --overview-file="README.md" src

download-structures:
	poetry run download-structures

.DEFAULT_GOAL := help