from __future__ import annotations

import argparse

from seacofs_eddy_dataset.config import load_config
from seacofs_eddy_dataset.pipeline import STAGE_BY_NAME, run_all, run_stage


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the modular SEACOFS eddy dataset pipeline.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_stage_parser = subparsers.add_parser("run-stage", help="Run one pipeline stage.")
    run_stage_parser.add_argument("stage", choices=sorted(STAGE_BY_NAME))
    run_stage_parser.add_argument("--config", required=True, help="Path to a YAML config file.")

    run_all_parser = subparsers.add_parser("run-all", help="Run all pipeline stages in order.")
    run_all_parser.add_argument("--config", required=True, help="Path to a YAML config file.")

    list_parser = subparsers.add_parser("list-stages", help="List available pipeline stages.")
    list_parser.add_argument("--config", required=False, help="Accepted for command symmetry; ignored.")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "list-stages":
        for name, stage in STAGE_BY_NAME.items():
            print(f"{name}: {stage.description}")
        return

    config = load_config(args.config)

    if args.command == "run-stage":
        run_stage(args.stage, config)
    elif args.command == "run-all":
        run_all(config)
    else:
        parser.error(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
