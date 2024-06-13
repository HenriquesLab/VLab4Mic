from datetime import date
import toml


def main():
    pyproject = toml.load("pyproject.toml")
    today = date.today().strftime("%Y.%m.%d")
    if pyproject["tool"]["poetry"]["version"] != today:
        pyproject["tool"]["poetry"]["version"] = today
        toml.dump(pyproject, open("pyproject.toml", "w"))


if __name__ == "__main__":
    main()
