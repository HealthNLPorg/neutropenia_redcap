import argparse
import logging
import os
from collections.abc import Collection, Iterable
from itertools import chain
from operator import attrgetter, methodcaller

import polars as pl
from more_itertools import map_reduce
from variants.sources import TextSource

from .redcap.forms.generic import REDCapForm
from .utils.filename import get_mrn
from .utils.formats import Formats
from .utils.interpretation import attributes_to_text_source

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
parser.add_argument("--smoke_test", action="store_true")
logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


def parse_text_sources(clustered_attribute_df: pl.DataFrame) -> Collection[TextSource]:
    return {
        attributes_to_text_source(sentence=sentence, section=section, filename=filename)
        for sentence, section, filename in zip(
            clustered_attribute_df["Sentence"],
            clustered_attribute_df["Section"],
            clustered_attribute_df["Filename"],
        )
    }


def parse_forms(
    mrn: int, corpus: str, mrn_cluster_df: pl.DataFrame
) -> Iterable[REDCapForm]:
    raise NotImplementedError


def tsv_to_redcap(
    data_location: str, output_dir: str, corpus: str, smoke_test: bool
) -> None:
    raw_output_frame = pl.read_csv(data_location, separator="\t").with_columns(
        pl.col(pl.String).replace("__UNK__", None)
    )
    if smoke_test:

        def mrn_fn(fn: str) -> str:
            return f"UPLOAD_TEST_{get_mrn(fn)}"
    else:
        mrn_fn = get_mrn
    forms = chain.from_iterable(
        parse_forms(mrn, corpus, mrn_cluster_df)
        for (mrn,), mrn_cluster_df in raw_output_frame.with_columns(
            MRN=raw_output_frame["Filename"].map_elements(mrn_fn)
        ).group_by("MRN")
    )
    for form_type_name, form_frame in map_reduce(
        iterable=forms,
        keyfunc=attrgetter("form_name"),
        valuefunc=methodcaller("to_data_frame"),
        reducefunc=pl.concat,
    ).items():
        if smoke_test:
            form_frame = form_frame.head(2)
        form_frame.write_csv(
            os.path.join(output_dir, f"{form_type_name}_redcap_upoad.csv")
        )


def convert(
    data_location: str,
    output_dir: str,
    input_format: Formats,
    output_format: Formats,
    corpus: str,
    smoke_test: bool,
) -> None:
    if input_format == output_format:
        logger.error("Input and output formats are both %s, exiting", input_format.name)
    match input_format, output_format:
        case Formats.RAW_TSV, Formats.REDCAP:
            tsv_to_redcap(
                data_location=data_location,
                output_dir=output_dir,
                corpus=corpus,
                smoke_test=smoke_test,
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
        smoke_test=args.smoke_test,
    )


if __name__ == "__main__":
    main()
