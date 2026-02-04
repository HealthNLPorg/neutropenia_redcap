import argparse
import datetime
import logging
import os
from collections.abc import Collection
from functools import cache

import polars as pl

from .filename_utils import get_mrn
from .formats import Formats, valid_format_choices
from .redcap.scnir import SCNIRForm, SCNIRGeneMention, SCNIRVariant, TextSource

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
parser.add_argument("--smoke_test", action="store_true")
logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


@cache
def parse_date_from_filename(filename: str) -> datetime.date | None:
    raw_date = (
        filename.split("-")[-2] if filename.endswith("lab") else filename.split("-")[-1]
    )
    try:
        month, day, year = (int(attr) for attr in raw_date.split("_"))
        return datetime.date(year=year, month=month, day=day)
    except Exception:
        logger.error("Bad filename, no date information: %s", filename)
        return None


def attributes_to_text_source(sentence: str, section: str, filename: str) -> TextSource:
    return TextSource(
        filename=filename,
        section=section,
        sentence=sentence,
        file_date=parse_date_from_filename(filename),
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


def get_variant(
    gene: str,
    clustered_attributes: tuple[str, ...],
    clustered_attribute_df: pl.DataFrame,
) -> SCNIRVariant:
    syntax_n, syntax_p, variant_type, vaf = clustered_attributes
    heterozygous = (
        True if vaf is not None and vaf.strip().lower() == "heterozygous" else None
    )
    return SCNIRVariant(
        gene=gene,
        syntax_p=syntax_p,
        syntax_n=syntax_n,
        variant_type=variant_type,
        vaf=vaf,
        heterozygous=heterozygous,
        text_sources=parse_text_sources(clustered_attribute_df=clustered_attribute_df),
        specimen_collection_dates={
            specimen_collection_date
            for specimen_collection_date in clustered_attribute_df[
                "Specimen_Collection_Date"
            ]
            if specimen_collection_date is not None
        },
        sample_sources={
            sample_source
            for sample_source in clustered_attribute_df["Sample_Source"]
            if sample_source is not None
        },
    )


def get_variants(
    gene: str,
    gene_cluster_df: pl.DataFrame,
    attributes: tuple[str, ...] = ("Syntax_N", "Syntax_P", "Type", "Vaf"),
) -> Collection[SCNIRVariant]:
    return {
        get_variant(
            gene=gene,
            clustered_attributes=clustered_attributes,
            clustered_attribute_df=clustered_attribute_df,
        )
        for clustered_attributes, clustered_attribute_df in gene_cluster_df.group_by(
            *attributes
        )
    }


def get_gene_mention(gene: str, gene_cluster_df: pl.DataFrame) -> SCNIRGeneMention:
    return SCNIRGeneMention(
        gene=gene, variants=get_variants(gene=gene, gene_cluster_df=gene_cluster_df)
    )


def get_gene_mentions(mrn_cluster_df: pl.DataFrame) -> Collection[SCNIRGeneMention]:
    return {
        get_gene_mention(gene=gene, gene_cluster_df=gene_cluster_df)
        for (gene,), gene_cluster_df in mrn_cluster_df.group_by("Gene")
    }


def mrn_cluster_to_form(mrn: int, mrn_cluster_df: pl.DataFrame) -> SCNIRForm:
    return SCNIRForm(mrn=mrn, gene_mentions=get_gene_mentions(mrn_cluster_df))


def raw_output_to_redcap(data_location: str, output_dir: str, smoke_test: bool) -> None:
    raw_output_frame = pl.read_csv(data_location, separator="\t").with_columns(
        pl.col(pl.String).replace("__UNK__", None)
    )
    if smoke_test:

        def mrn_fn(fn: str) -> str:
            return f"UPLOAD_TEST_{get_mrn(fn)}"
    else:
        mrn_fn = get_mrn
    final_frame = pl.concat(
        mrn_cluster_to_form(mrn, mrn_cluster_df).to_data_frame()
        for (mrn,), mrn_cluster_df in raw_output_frame.with_columns(
            MRN=raw_output_frame["Filename"].map_elements(mrn_fn)
        ).group_by("MRN")
    )
    if smoke_test:
        final_frame = final_frame.head(2)
    final_frame.write_csv(os.path.join(output_dir, "redcap_upoad.csv"))


def convert(
    data_location: str,
    output_dir: str,
    input_format: Formats,
    output_format: Formats,
    smoke_test: bool,
) -> None:
    if input_format == output_format:
        logger.error("Input and output formats are both %s, exiting", input_format.name)
    match input_format, output_format:
        case Formats.RAW_TSV, Formats.REDCAP:
            raw_output_to_redcap(data_location, output_dir, smoke_test)
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
    convert(
        args.data_location,
        args.output_dir,
        input_format,
        output_format,
        args.smoke_test,
    )


if __name__ == "__main__":
    main()
