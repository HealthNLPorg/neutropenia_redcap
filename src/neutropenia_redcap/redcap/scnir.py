from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date

import polars as pl

from .redcap_import import (
    MAXIMUM_GERMLINES,
    MAXIMUM_VARIANTS,
    MINIMUM_GERMLINES,
    MINIMUM_VARIANTS,
)


@dataclass
class SCNIRVariant:
    syntax_p: str | None
    syntax_n: str | None
    variant_type: str | None
    vaf: str | None
    known_heterozygous: bool
    specimen_collection_date: date | None
    sample_source: str | None
    source_filenames: list[str]

    def to_row_fragment(self) -> Iterable[str | None | bool]:
        return []


@dataclass
class SCNIRGeneMention:
    gene: str
    variants: list[SCNIRVariant]

    def __post_init__(self) -> None:
        if (
            len(self.variants) < MINIMUM_VARIANTS
            or len(self.variants) > MAXIMUM_VARIANTS
        ):
            raise ValueError(
                f"Too many or too few variant mentions {len(self.variants)}"
            )

        def to_row_fragment(self) -> Iterable[str | bool | None]:
            # Hard-coded weirdness - hopefully only for now
            opening_cells = 18
            closing_cells = 2
            for _ in range(opening_cells):
                yield None

            for _ in range(closing_cells):
                yield None


@dataclass
class SCNIRForm:
    mrn: int
    gene_mentions: list[SCNIRGeneMention]

    def __post_init__(self) -> None:
        if (
            len(self.gene_mentions) < MINIMUM_GERMLINES
            or len(self.gene_mentions) > MAXIMUM_GERMLINES
        ):
            raise ValueError(
                f"Too many or too few gene/germline mentions {len(self.gene_mentions)}"
            )

    def to_row(self) -> Iterable[str | bool | None]:
        return []

    def to_data_frame(self) -> pl.DataFrame:
        return pl.DataFrame()
