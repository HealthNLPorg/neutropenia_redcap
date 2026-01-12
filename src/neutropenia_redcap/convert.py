import argparse
import logging

from .formats import Formats, valid_format_choices

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data_location", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument(
    "--input_format",
    type=str,
    choices=valid_format_choices,
    default=Formats.RAW_TSV.name,
)
parser.add_argument(
    "--output_format",
    type=str,
    choices=valid_format_choices,
    default=Formats.REDCAP.name,
)
logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


def raw_output_to_redcap(data_location: str, output_dir: str) -> None:
    pass


def convert(
    data_location: str, output_dir: str, input_format: Formats, output_format: Formats
) -> None:
    if input_format == output_format:
        logger.error("Input and output formats are both %s, exiting", input_format.name)
    match input_format, output_format:
        case Formats.RAW_TSV, Formats.REDCAP:
            raw_output_to_redcap(data_location, output_dir)
        case _:
            raise ValueError(
                "%s to %s not currently supported",
                input_format.name,
                output_format.name,
            )


def main() -> None:
    args = parser.parse_args()
    input_format = Formats(args.input_format)
    output_format = Formats(args.output_format)
    convert(args.data_location, args.output_dir, input_format, output_format)


if __name__ == "__main__":
    main()
