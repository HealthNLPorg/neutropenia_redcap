import logging
from collections import Counter
from collections.abc import Collection, Iterable, Mapping, Set
from itertools import chain
from operator import itemgetter
from typing import cast

import polars as pl
from more_itertools import map_reduce

from neutropenia_redcap.redcap.forms.generic import REDCapForm, SupportedFormTypes
from neutropenia_redcap.redcap.forms.germline.scnir import SCNIRGermlineForm
from neutropenia_redcap.redcap.forms.germline.sds import SDSGermlineForm
from neutropenia_redcap.redcap.forms.somatic import SomaticForm
from neutropenia_redcap.variants import (
    GeneMention,
    SupportedMentionTypes,
    SupportedVariantTypes,
    Variant,
)
from neutropenia_redcap.variants.germline.generic import (
    GermlineGeneMention,
    GermlineVariant,
)
from neutropenia_redcap.variants.germline.scnir import SCNIRGermlineGeneMention
from neutropenia_redcap.variants.germline.sds import SDSGermlineGeneMention
from neutropenia_redcap.variants.somatic import SomaticGeneMention, SomaticVariant
from neutropenia_redcap.variants.sources import TextSource

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
    clustered_attribute_df: pl.DataFrame,
) -> Counter[SupportedVariantTypes]:
    return Counter(
        attributes_to_variant_type(section=section, filename=filename)
        for section, filename in zip(
            clustered_attribute_df["Section"],
            clustered_attribute_df["Filename"],
        )
    )


def get_variant(
    gene: str,
    clustered_attributes: tuple[str, ...],
    clustered_attribute_df: pl.DataFrame,
) -> Variant:
    syntax_n, syntax_p, variant_type, vaf = clustered_attributes
    heterozygous = (
        True if vaf is not None and vaf.strip().lower() == "heterozygous" else None
    )
    text_sources = parse_text_sources(clustered_attribute_df=clustered_attribute_df)
    variant_types = parse_variant_types(clustered_attribute_df=clustered_attribute_df)
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
        case SupportedVariantTypes.GERMLINE:
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


def get_variants(
    gene: str,
    gene_cluster_df: pl.DataFrame,
    attributes: tuple[str, ...] = ("Syntax_N", "Syntax_P", "Type", "Vaf"),
) -> Set[Variant]:
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


def get_mentions_for_gene(
    gene: str, gene_cluster_df: pl.DataFrame
) -> Mapping[SupportedMentionTypes, GeneMention]:
    def variant_to_mention_type(variant: Variant) -> SupportedMentionTypes:
        if isinstance(variant, SomaticVariant):
            return SupportedMentionTypes.SOMATIC
        elif isinstance(variant, GermlineVariant):
            return SupportedMentionTypes.GERMLINE
        else:
            raise ValueError(f"Unsupported variant type for variant {variant}")

    def variants_to_mention(variants: Iterable[Variant]) -> GeneMention:
        variants = set(variants)
        variant = next(iter(variants))
        mention_type = variant_to_mention_type(variant)
        match mention_type:
            case SupportedMentionTypes.SOMATIC:
                return SomaticGeneMention(
                    gene=gene, variants=cast(Set[SomaticVariant], variants)
                )
            case SupportedMentionTypes.GERMLINE:
                return GermlineGeneMention(
                    gene=gene, variants=cast(Set[GermlineVariant], variants)
                )
            case _:
                raise ValueError(f"Unsupported variant type for variant {variant}")

    return map_reduce(
        iterable=get_variants(gene, gene_cluster_df=gene_cluster_df),
        keyfunc=variant_to_mention_type,
        reducefunc=variants_to_mention,
    )


def collect_gene_mentions_for_mrn(
    mrn_cluster_df: pl.DataFrame,
) -> Mapping[SupportedMentionTypes, Set[GeneMention]]:
    return map_reduce(
        iterable=chain.from_iterable(
            get_mentions_for_gene(gene=gene, gene_cluster_df=gene_cluster_df).items()
            for (gene,), gene_cluster_df in mrn_cluster_df.group_by("Gene")
        ),
        keyfunc=itemgetter(0),
        valuefunc=itemgetter(1),
        reducefunc=set,
    )


def mrn_cluster_to_forms(
    mrn: int, corpus: str, mrn_cluster_df: pl.DataFrame
) -> Mapping[SupportedFormTypes, REDCapForm]:
    def mention_type_to_form_type(
        mention_type: SupportedMentionTypes,
    ) -> SupportedFormTypes:
        match corpus, mention_type:
            case "sds", SupportedMentionTypes.GERMLINE:
                return SupportedFormTypes.SDS_GERMLINE
            case "scnir", SupportedMentionTypes.GERMLINE:
                return SupportedFormTypes.SCNIR_GERMLINE
            case _, SupportedMentionTypes.SOMATIC:
                return SupportedFormTypes.SOMATIC
            case _:
                raise ValueError(
                    f"Unsupported corpus and mention type combination {corpus} {mention_type}"
                )

    result = {}
    for mention_type, mentions in collect_gene_mentions_for_mrn(
        mrn_cluster_df=mrn_cluster_df
    ).items():
        form_type = mention_type_to_form_type(mention_type=mention_type)
        match form_type:
            case SupportedFormTypes.SDS_GERMLINE:
                result[form_type] = SDSGermlineForm(
                    mrn=mrn,
                    gene_mentions={
                        SDSGermlineGeneMention(
                            gene=mention.gene,
                            variants={
                                GermlineVariant(
                                    gene=variant.gene,
                                    syntax_n=variant.syntax_n,
                                    syntax_p=variant.syntax_p,
                                    specimen_collection_dates=variant.specimen_collection_dates,
                                    variant_type=variant.variant_type,
                                    vaf=variant.vaf,
                                    sample_sources=variant.sample_sources,
                                    heterozygous=variant.heterozygous,
                                    text_sources=variant.text_sources,
                                )
                                for variant in mention.variants
                            },
                        )
                        for mention in mentions
                    },
                )

            case SupportedFormTypes.SCNIR_GERMLINE:
                result[form_type] = SCNIRGermlineForm(
                    mrn=mrn,
                    gene_mentions={
                        SCNIRGermlineGeneMention(
                            gene=mention.gene,
                            variants={
                                GermlineVariant(
                                    gene=variant.gene,
                                    syntax_n=variant.syntax_n,
                                    syntax_p=variant.syntax_p,
                                    specimen_collection_dates=variant.specimen_collection_dates,
                                    variant_type=variant.variant_type,
                                    vaf=variant.vaf,
                                    sample_sources=variant.sample_sources,
                                    heterozygous=variant.heterozygous,
                                    text_sources=variant.text_sources,
                                )
                                for variant in mention.variants
                            },
                        )
                        for mention in mentions
                    },
                )

            case SupportedFormTypes.SOMATIC:
                result[form_type] = SomaticForm(
                    mrn=mrn,
                    gene_mentions={
                        SomaticGeneMention(
                            gene=mention.gene,
                            variants={
                                SomaticVariant(
                                    gene=variant.gene,
                                    syntax_n=variant.syntax_n,
                                    syntax_p=variant.syntax_p,
                                    specimen_collection_dates=variant.specimen_collection_dates,
                                    variant_type=variant.variant_type,
                                    vaf=variant.vaf,
                                    sample_sources=variant.sample_sources,
                                    heterozygous=variant.heterozygous,
                                    text_sources=variant.text_sources,
                                )
                                for variant in mention.variants
                            },
                        )
                        for mention in mentions
                    },
                )
            case _:
                raise ValueError(f"Unsupported form type {form_type}")
    return result
