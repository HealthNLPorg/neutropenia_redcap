import logging
from collections import Counter
from collections.abc import Collection

import polars as pl
from variants import SupportedVariantTypes, Variant
from variants.germline.generic import GermlineVariant
from variants.somatic import SomaticVariant
from variants.sources import TextSource

from .interpretation import attributes_to_text_source, attributes_to_variant_type

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


def parse_variant_types(
    corpus: str, clustered_attribute_df: pl.DataFrame
) -> Counter[SupportedVariantTypes]:
    return Counter(
        attributes_to_variant_type(corpus=corpus, section=section, filename=filename)
        for section, filename in zip(
            clustered_attribute_df["Section"],
            clustered_attribute_df["Filename"],
        )
    )


def get_variant(
    gene: str,
    corpus: str,
    clustered_attributes: tuple[str, ...],
    clustered_attribute_df: pl.DataFrame,
) -> Variant:
    syntax_n, syntax_p, variant_type, vaf = clustered_attributes
    heterozygous = (
        True if vaf is not None and vaf.strip().lower() == "heterozygous" else None
    )
    text_sources = parse_text_sources(clustered_attribute_df=clustered_attribute_df)
    variant_types = parse_variant_types(
        corpus=corpus, clustered_attribute_df=clustered_attribute_df
    )
    total_distinct_variant_types = len(variant_types)
    match total_distinct_variant_types:
        case 0:
            raise ValueError(
                f"Issue with parsing cluster for variant for gene {gene} at protein {syntax_p} nucleotide {syntax_n}"
            )
        case 1:
            selected_variant_type = next(iter(variant_types.keys()))
            if selected_variant_type == SupportedVariantTypes.NA:
                raise ValueError(
                    f"No known variant type found for gene {gene} at protein {syntax_p} nucleotide {syntax_n}"
                )
        case _:
            logger.warning(
                "Variant for gene %s at protein %s nucleotide %s has %d distinct variant type, selecting most common variant type",
                gene,
                syntax_p,
                syntax_n,
                total_distinct_variant_types,
            )
            selected_variant_type, _ = next(
                filter(
                    lambda t: t[0] != SupportedVariantTypes.NA,
                    variant_types.most_common(),
                )
            )
    match selected_variant_type:
        case SupportedVariantTypes.SCNIR_GERMLINE | SupportedVariantTypes.SDS_GERMLINE:
            return GermlineVariant(
                gene=gene,
                syntax_p=syntax_p,
                syntax_n=syntax_n,
                variant_type=variant_type,
                vaf=vaf,
                heterozygous=heterozygous,
                text_sources=text_sources,
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

        case SupportedVariantTypes.SOMATIC:
            return SomaticVariant(
                gene=gene,
                syntax_p=syntax_p,
                syntax_n=syntax_n,
                variant_type=variant_type,
                vaf=vaf,
                heterozygous=heterozygous,
                text_sources=text_sources,
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
        case _:
            raise ValueError(f"Unsupported variant type {_}")
