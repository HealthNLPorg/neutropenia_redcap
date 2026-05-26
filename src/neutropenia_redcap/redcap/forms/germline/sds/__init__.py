from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum
from itertools import chain

import polars as pl

from neutropenia_redcap.utils.iter import up_to_n
from neutropenia_redcap.variants.germline.sds import (
    MAXIMUM_SDS_ALLELES,
    MAXIMUM_SDS_GERMLINE_VARIANTS,
    MAXIMUM_SDS_GERMLINES,
    MINIMUM_SDS_ALLELES,
    MINIMUM_SDS_GERMLINE_VARIANTS,
    MINIMUM_SDS_GERMLINES,
    SDSGermlineGeneMention,
)

from ..generic import (
    GermlineForm,
)


# gotta love ridiculous database internal shorthand
class Parent(StrEnum):
    father = "fath"
    mother = "mot"


INTERVENING_COLUMNS = ("vus_gen_sds", "comments_gene_gen_sds", "gen_family")

FINAL_COLUMNS = (
    "sdbs_upload_gen_sds",
    "parental_upload_gen_sds",
    "msc_upload_gen_sds",
    "genetic_comments",
)


def parent_and_allele_indexed_columns(
    parent: Parent, allele_index: int
) -> Sequence[str]:
    if not (1 <= allele_index <= 2):
        raise ValueError(f"Unsupported parental allele index: {allele_index}")
    return (
        f"{parent.value}_allele{allele_index}_gen_sds",
        f"{parent.value}_allele{allele_index}_oth_gen_sds",
    )


def parent_indexed_columns(parent: Parent) -> Iterable[str]:
    return chain(
        chain.from_iterable(
            parent_and_allele_indexed_columns(parent=parent, allele_index=allele_index)
            for allele_index in range(MINIMUM_SDS_ALLELES, MAXIMUM_SDS_ALLELES + 1)
        ),
        (f"{parent.value}_vus",),
    )


def gene_and_variant_indexed_columns(
    gene_index: int, variant_index: int
) -> Sequence[str]:
    gene_verbal_ordinals = ["gene_", "second", "third"]
    suffix_ordinals = ["", "_2", "_3"]
    if not (1 <= variant_index <= 4):
        raise ValueError(f"Unsupported variant index: {variant_index}")
    if not (1 <= gene_index <= 3):
        raise ValueError(f"Unsupported gene index: {gene_index}")

    def alleles_nonsense(_gene_index: int) -> str:
        match _gene_index:
            case 1:
                return ""
            case 2 | 3:
                return f"_gene{_gene_index}"
            case other:
                raise ValueError(f"Unsupported index {other}")

    return (
        f"{gene_verbal_ordinals[gene_index - 1]}mut_gen_sds",
        f"gene_mut_gen_sds_{gene_index}",  # What gene has mutations and/or was tested
        f"gene_mut_oth{suffix_ordinals[gene_index - 1]}",
        f"alleles{alleles_nonsense(gene_index)}_gen_sds",  # Number of variants in gene
        f"allele{variant_index}{alleles_nonsense(gene_index)}_gen_sds",
        f"allele{variant_index}{alleles_nonsense(gene_index)}_oth_gen_sds",
        f"protein{variant_index}{alleles_nonsense(gene_index)}_gen_sds",
        f"protein{variant_index}{alleles_nonsense(gene_index)}_oth_gen_sds",
    )


def gene_indexed_columns(gene_index: int) -> Iterable[str]:
    if not (1 <= gene_index <= 3):
        raise ValueError(f"Unsupported gene index: {gene_index}")
    return chain.from_iterable(
        gene_and_variant_indexed_columns(
            gene_index=gene_index, variant_index=variant_index
        )
        for variant_index in range(
            MINIMUM_SDS_GERMLINE_VARIANTS, MAXIMUM_SDS_GERMLINE_VARIANTS + 1
        )
    )


_SDS_GERMLINE_COLUMNS = chain(
    (
        "patient_id",
        "pat_done_gen_sds",  # Was genetic testing done for this participant
        "gene_pos_gen_sds",  # Gene positive or negative NB: How is this a global question when it would apply to each gene?
    ),
    chain.from_iterable(
        gene_indexed_columns(gene_index=germline_index)
        for germline_index in range(MINIMUM_SDS_GERMLINES, MAXIMUM_SDS_GERMLINES + 1)
    ),
    INTERVENING_COLUMNS,
    chain.from_iterable(
        parent_indexed_columns(parent=parent) for parent in sorted(Parent, reverse=True)
    ),
    FINAL_COLUMNS,
)


@dataclass
class SDSGermlineForm(GermlineForm):
    gene_mentions: Collection[SDSGermlineGeneMention]
    schema = [(column_name, pl.String) for column_name in _SDS_GERMLINE_COLUMNS]

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        yield min(len(self.gene_mentions), MAXIMUM_SDS_GERMLINES)
        for germline in up_to_n(
            self.gene_mentions, n=MAXIMUM_SDS_GERMLINES, fillvalue=None
        ):
            yield from (
                SDSGermlineGeneMention.blank_row_fragment()
                if germline is None
                else germline.to_row_fragment()
            )
