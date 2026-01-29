import logging
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from datetime import date
from functools import lru_cache
from itertools import chain, islice, repeat
from typing import cast

import polars as pl

from .redcap_import import (
    MAXIMUM_GERMLINES,
    MAXIMUM_VARIANTS,
    MINIMUM_VARIANTS,
    SCNIR_COLUMNS,
)

logger = logging.getLogger(__name__)


def pad_iterable[T, S](iterable: Iterable[T], n: int, pad: S) -> Iterable[T | S]:
    return islice(chain(iterable, cast(Iterable[S], repeat(pad))), n)


logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)

SCNIR_SCHEMA = [(column_name, pl.String) for column_name in SCNIR_COLUMNS]


@lru_cache
def map_variant_type(variant_type: str | None) -> int | None:
    if variant_type is None:
        return None
    normalized_variant_type = " ".join(variant_type.lower().split())
    is_likely = "likely" in normalized_variant_type
    is_benign = "benign" in normalized_variant_type
    is_pathogenic = "patho" in normalized_variant_type
    is_uncertain = (
        "uncertain" in normalized_variant_type
        or "vus" in normalized_variant_type
        or "unknown significance" in normalized_variant_type
    )
    match is_likely, is_benign, is_pathogenic, is_uncertain:
        # Something weird
        # case False, False, False, False:
        #     # Don't know - leave to comment or further tweaking
        #     logger.warning("Not sure how to map: %s", variant_type)
        #     return None
        # Straightforward pathogenic
        case False, False, True, False:
            return 1
        # Likely pathogenic
        case True, False, True, False:
            return 2
        # Straightforward benign
        case False, True, False, False:
            return 3
        # Likely benign
        case True, True, False, False:
            return 4
        # Straightforward uncertain
        case False, False, False, True:
            return 5
    logger.warning("Cannot currently map: %s", variant_type)
    return None


@dataclass
class TextSource:
    filename: str
    section: str
    sentence: str


@dataclass
class SCNIRVariant:
    syntax_p: str | None
    syntax_n: str | None
    variant_type: str | None
    vaf: str | None
    heterozygous: (
        bool | None
    )  # True for is heterozygous, False for definitely isn't, None for unknown
    specimen_collection_dates: Collection[date] = set()
    sample_sources: Collection[str] = set()
    text_sources: Collection[TextSource] = set()

    # Another weird thing is I can't find the field where the specimen collection date
    # would go
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | None | bool]:
        if blank:
            yield from SCNIRVariant.blank_row_fragment()
        # sum_germ_var{variant_index}_cdna_{germline_index}
        yield self.syntax_n
        # sum_germ_var{variant_index}_pro_{germline_index}
        yield self.syntax_p
        # sum_germ_var{variant_index}_acmg_{germline_index}
        yield map_variant_type(self.variant_type)
        # sum_germ_var{variant_index}_comment_{germline_index}
        yield self.build_comment()

    # Heterozygosity in "Clinical and Research Sequencing Summary Form Comments"
    # or variant level comments
    def build_comment(self) -> str:
        return ""

    @staticmethod
    def blank_row_fragment() -> Iterable[None]:
        for _ in range(4):
            yield None


@dataclass
class SCNIRGeneMention:
    gene: str
    variants: Collection[SCNIRVariant] = set()

    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SCNIRGeneMention.blank_row_fragment()
        # sum_germ_gene_{germline_index}
        yield self.gene
        # sum_germ_num_var_{germline_index}
        yield min(len(self.variants), MAXIMUM_VARIANTS)
        for variant in pad_iterable(self.variants, n=MAXIMUM_VARIANTS, pad=None):
            yield from (
                variant.to_row_fragment()
                if variant is not None
                else SCNIRVariant.blank_row_fragment()
            )

    @staticmethod
    def blank_row_fragment() -> Iterable[None]:
        yield None
        yield None
        for i in range(MINIMUM_VARIANTS, MAXIMUM_VARIANTS + 1):
            yield from SCNIRVariant.blank_row_fragment()


@dataclass
class SCNIRForm:
    mrn: int
    gene_mentions: list[SCNIRGeneMention]

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        yield min(len(self.gene_mentions), MAXIMUM_GERMLINES)
        for germline in pad_iterable(self.gene_mentions, n=MAXIMUM_GERMLINES, pad=None):
            yield from (
                germline.to_row_fragment()
                if germline is not None
                else SCNIRGeneMention.blank_row_fragment()
            )

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCNIR_SCHEMA, orient="row")
