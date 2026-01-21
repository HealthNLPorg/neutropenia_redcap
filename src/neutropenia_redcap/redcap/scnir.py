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

    # Another weird thing is I can't find the field where the specimen collection date
    # would go
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | None | bool]:
        if blank:
            yield from SCNIRVariant.blank_row_fragment()
        # f"Germline variant {variant} (genomic or mitochondrial/mtDNA) [{n}]:",
        yield None  # even the clinicians don't really fill this out
        # f"Germline variant {variant} (cDNA) [{n}]:",
        yield self.syntax_n
        # f"Germline variant {variant} (protein) [{n}]:",
        yield self.syntax_p
        # f"Was germline variant {variant_id} confirmed by parental genetic testing? [{n}]",
        yield None  # clinician needs to assess this
        # type of... they'll need to check this
        for _ in range(12):
            yield "Unchecked"
        # Specify... Germline... Mitochondrial
        yield None
        yield None
        yield None
        # ACMG
        yield self.variant_type
        yield self.build_comment()

    # Heterozygosity in "Clinical and Research Sequencing Summary Form Comments"
    # or variant level comments
    def build_comment(self) -> str:
        return ""

    @staticmethod
    def blank_row_fragment() -> Iterable[None]:
        for _ in range(21):
            yield None


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

        # Corresponds to nth_germline_testing_information
        def to_row_fragment(self) -> Iterable[str | bool | None]:
            # Hard-coded weirdness - hopefully only for now
            total_opening_cells = 18
            total_closing_cells = 2
            # this is where we would have to do some logic with self.sample_source
            for _ in range(total_opening_cells):
                yield None
            # nth_germline_information
            # total_open_internal_cells = 3
            yield self.gene
            for _ in range(2):
                yield None
            for i in range(MINIMUM_VARIANTS, MAXIMUM_VARIANTS):
                if i < len(self.variants):
                    yield from self.variants[i].get_row_fragment()
                else:
                    yield from SCNIRVariant.blank_row_fragment()
            for _ in range(total_closing_cells):
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
