import argparse
import logging
import os

import polars as pl

from .filename_utils import get_mrn
from .formats import Formats, valid_format_choices
from .redcap.scnir import SCNIRForm, SCNIRGeneMention, SCNIRVariant

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


def get_variant(
    clustered_attributes: tuple[str, ...], clustered_attribute_df: pl.DataFrame
) -> SCNIRVariant:
    syntax_n, syntax_p, variant_type, vaf, section = clustered_attributes
    known_heterozygous = vaf.strip().lower() == "heterozygous"
    return SCNIRVariant(
        syntax_p=syntax_p,
        syntax_n=syntax_n,
        variant_type=variant_type,
        vaf=vaf,
        known_heterozygous=known_heterozygous,
        specimen_collection_date=next(
            (
                specimen_collection_date
                for specimen_collection_date in clustered_attribute_df[
                    "Specimen_Collection_Date"
                ]
                if specimen_collection_date != "__UNK__"
            ),
            None,
        ),
        sample_source=next(
            (
                sample_source
                for sample_source in clustered_attribute_df["Sample_Source"]
                if sample_source != "__UNK__"
            ),
            None,
        ),
        source_filenames=[
            fn for fn in clustered_attribute_df["Filename"] if fn != "__UNK__"
        ],
    )


def get_variants(
    gene_cluster_df: pl.DataFrame,
    attributes: tuple[str, ...] = ("Syntax_N", "Syntax_P", "Type", "Vaf", "Section"),
) -> list[SCNIRVariant]:
    return [
        get_variant(clustered_attributes, clustered_attribute_df)
        for clustered_attributes, clustered_attribute_df in gene_cluster_df.group_by(
            *attributes
        )
    ]


def get_gene_mention(gene: str, gene_cluster_df: pl.DataFrame) -> SCNIRGeneMention:
    return SCNIRGeneMention(variants=get_variants(gene_cluster_df))


def get_gene_mentions(mrn_cluster_df: pl.DataFrame) -> list[SCNIRGeneMention]:
    return [
        get_gene_mention(gene, gene_cluster_df)
        for (gene,), gene_cluster_df in mrn_cluster_df.group_by("Gene")
    ]


def mrn_cluster_to_form(mrn: int, mrn_cluster_df: pl.DataFrame) -> SCNIRForm:
    return SCNIRForm(mrn=mrn, gene_mentions=get_gene_mentions(mrn_cluster_df))


def raw_output_to_redcap(data_location: str, output_dir: str) -> None:
    raw_output_frame = pl.read_csv(data_location, separator="\t")
    final_frame = pl.concat(
        mrn_cluster_to_form(mrn, mrn_cluster_df).to_data_frame()
        for (mrn,), mrn_cluster_df in raw_output_frame.with_columns(
            MRN=raw_output_frame["Filename"].map_elements(get_mrn)
        ).group_by("MRN")
    )
    final_frame.write_csv(os.path.join(output_dir, "redcap_upoad.csv"))


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
