from __future__ import annotations

import argparse

from .config import load_config
from .pipeline import STAGE_BY_NAME, run_all, run_stage


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the AVISO surface eddy dataset pipeline.")
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("run-stage")
    one.add_argument("stage", choices=sorted(STAGE_BY_NAME))
    one.add_argument("--config", required=True)
    all_stages = commands.add_parser("run-all")
    all_stages.add_argument("--config", required=True)
    commands.add_parser("list-stages")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "list-stages":
        for stage in STAGE_BY_NAME.values():
            print(f"{stage.name}: {stage.description}")
        return
    config = load_config(args.config)
    if args.command == "run-stage":
        run_stage(args.stage, config)
    else:
        run_all(config)


if __name__ == "__main__":
    main()
