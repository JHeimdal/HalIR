import argparse
import json
import sys

from pathlib import Path
from typing import Optional

from jsonschema import Draft202012Validator


def load_json(path: Path) -> dict:
    with path.open() as infile:
        return json.load(infile)


def resolve_schema_path(input_path: Path, schema_arg: Optional[str]) -> Path:
    if schema_arg is not None:
        return Path(schema_arg).expanduser().resolve()

    input_doc = load_json(input_path)
    schema_ref = input_doc.get("$schema")
    if not isinstance(schema_ref, str):
        raise ValueError("Input JSON does not define a $schema field and no --schema was provided")

    schema_path = (input_path.parent / schema_ref).resolve()
    return schema_path


def validate_input(input_path: Path, schema_path: Path) -> int:
    input_doc = load_json(input_path)
    schema_doc = load_json(schema_path)

    Draft202012Validator.check_schema(schema_doc)
    validator = Draft202012Validator(schema_doc)
    errors = sorted(validator.iter_errors(input_doc), key=lambda err: list(err.path))

    if errors:
        for error in errors:
            json_path = ".".join(str(part) for part in error.absolute_path)
            if json_path:
                print(f"Validation error at {json_path}: {error.message}", file=sys.stderr)
            else:
                print(f"Validation error: {error.message}", file=sys.stderr)
        return 1

    print(f"Validated {input_path} against {schema_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate a HalIR input JSON file against its schema")
    parser.add_argument(
        "input",
        nargs="?",
        default=str(Path(__file__).with_name("input.json")),
        help="Input JSON file to validate",
    )
    parser.add_argument(
        "--schema",
        "-s",
        help="Optional explicit schema file path; otherwise the input file's $schema field is used",
    )

    args = parser.parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists() or not input_path.is_file():
        print(f"Input file not found: {input_path}", file=sys.stderr)
        return 2

    try:
        schema_path = resolve_schema_path(input_path, args.schema)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    if not schema_path.exists() or not schema_path.is_file():
        print(f"Schema file not found: {schema_path}", file=sys.stderr)
        return 2

    return validate_input(input_path, schema_path)


if __name__ == "__main__":
    raise SystemExit(main())