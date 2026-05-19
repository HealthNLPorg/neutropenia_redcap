import argparse
import logging
import os
from collections.abc import Iterable
from itertools import chain
from operator import itemgetter, methodcaller

import polars as pl
from more_itertools import bucket

from .redcap.forms.generic import REDCapForm
from .utils.dataframe import mrn_cluster_to_forms
from .utils.filename import get_mrn
from .utils.formats import Formats

VALID_FORMAT_CHOICES = [_format.name for _format in Formats]
parser = argparse.ArgumentParser(description="")
parser.add_argument("--data_location", type=str)
parser.add_argument("--output_dir", type=str)
parser.add_argument(
    "--input_format",
    type=str,
    help="Format of gene findings",
    choices=VALID_FORMAT_CHOICES,
    default=Formats.RAW_TSV.name,
)
parser.add_argument(
    "--output_format",
    type=str,
    help="REDCap input format",
    choices=VALID_FORMAT_CHOICES,
    default=Formats.REDCAP.name,
)

parser.add_argument(
    "--corpus",
    type=str,
    help="The corpus source of notes we are parsing (SDS or SCNIR), eventually would like to automatically infer this from the notes themselves",
    choices=("scnir", "sds"),
    default="scnir",
)
parser.add_argument("--debug", action="store_true", help="Debug options for REDCap")
logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


def cluster_forms(forms: Iterable[REDCapForm]) -> pl.LazyFrame:
    return pl.concat(map(methodcaller("to_data_frame"), forms))


def tsv_to_redcap(
    data_location: str, output_dir: str, corpus: str, debug: bool
) -> None:
    raw_output_frame = pl.read_csv(data_location, separator="\t").with_columns(
        pl.col(pl.String).replace("__UNK__", None)
    )
    if debug:

        def mrn_fn(fn: str) -> str:
            return f"UPLOAD_TEST_{get_mrn(fn)}"
    else:
        mrn_fn = get_mrn
    forms = chain.from_iterable(
        mrn_cluster_to_forms(mrn, corpus, mrn_cluster_df).items()
        for (mrn,), mrn_cluster_df in raw_output_frame.with_columns(
            MRN=raw_output_frame["Filename"].map_elements(mrn_fn)
        ).group_by("MRN")
    )
    form_bucket = bucket(forms, key=itemgetter(0))
    extant_form_types = list(form_bucket)
    for form_type in extant_form_types:
        frame = cluster_forms(map(itemgetter(1), form_bucket[form_type]))
        frame.sink_csv(
            os.path.join(output_dir, f"{form_type.name.lower()}_redcap_upload.csv")
        )


def convert(
    data_location: str,
    output_dir: str,
    input_format: Formats,
    output_format: Formats,
    corpus: str,
    debug: bool,
) -> None:
    if input_format == output_format:
        logger.error("Input and output formats are both %s, exiting", input_format.name)
    match input_format, output_format:
        case Formats.RAW_TSV, Formats.REDCAP:
            tsv_to_redcap(
                data_location=data_location,
                output_dir=output_dir,
                corpus=corpus,
                debug=debug,
            )
        case _:
            raise ValueError(
                f"{input_format.name} to {output_format.name} not currently supported",
            )


def main() -> None:
    args = parser.parse_args()
    input_format = Formats(args.input_format)
    output_format = Formats(args.output_format)
    convert(
        data_location=args.data_location,
        output_dir=args.output_dir,
        input_format=input_format,
        output_format=output_format,
        corpus=args.corpus,
        debug=args.debug,
    )


if __name__ == "__main__":
    main()
