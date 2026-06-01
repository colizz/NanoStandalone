#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def load_dataset_yaml(path):
    """Load a simple mapping of sample names to dataset-path lists."""
    try:
        import yaml

        with open(path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        return {name: list(datasets or []) for name, datasets in data.items()}
    except ImportError:
        pass

    data = {}
    current_name = None
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if not raw_line.startswith((" ", "\t")) and line.endswith(":"):
                current_name = line[:-1].strip()
                data[current_name] = []
                continue
            if current_name is None or not line.startswith("-"):
                raise ValueError(f"Unsupported YAML line in {path}: {raw_line.rstrip()}")
            data[current_name].append(line[1:].strip())
    return data


def das_files(dataset):
    """Return all files for one DAS dataset path."""
    query = f"file dataset={dataset}"
    result = subprocess.run(
        ["dasgoclient", f"-query={query}"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    files = [line.strip() for line in result.stdout.splitlines() if line.strip().startswith("/store/")]
    if not files:
        message = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"DAS query returned no /store files for {dataset}:\n{message}")
    return files


def main():
    parser = argparse.ArgumentParser(description="Build a DAS file list from a YAML sample definition.")
    parser.add_argument("-i", "--input", required=True, help="Input YAML path, e.g. samples/bkg.yaml")
    parser.add_argument("-o", "--output", required=True, help="Output txt path, e.g. samples/bkg.txt")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parents[1]
    input_path = Path(args.input)
    output_path = Path(args.output)
    if not input_path.is_absolute():
        input_path = base_dir / input_path
    if not output_path.is_absolute():
        output_path = base_dir / output_path

    samples = load_dataset_yaml(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_lines = 0
    with open(output_path, "w", encoding="utf-8") as output:
        for sample_name, datasets in samples.items():
            sample_index = 0
            for dataset in datasets:
                files = das_files(dataset)
                print(f"{sample_name}: {dataset} -> {len(files)} files")
                for file_path in files:
                    output.write(f"{file_path} {sample_name}_file{sample_index}\n")
                    sample_index += 1
                    n_lines += 1

    print(f"Wrote {n_lines} lines to {output_path}")


if __name__ == "__main__":
    main()
