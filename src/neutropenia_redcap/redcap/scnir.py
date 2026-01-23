from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date

import polars as pl

from .redcap_import import (
    MAXIMUM_GERMLINES,
    MAXIMUM_VARIANTS,
    MINIMUM_GERMLINES,
    MINIMUM_VARIANTS,
    build_scnir_columns,
)

SCHEMA = [(column, pl.String) for column in build_scnir_columns()]


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

    # This is getting handled later
    # def __post_init__(self) -> None:
    #     if (
    #         len(self.variants) < MINIMUM_VARIANTS
    #         or len(self.variants) > MAXIMUM_VARIANTS
    #     ):
    #         raise ValueError(
    #             f"Too many or too few variant mentions {len(self.variants)}"
    #         )

    # Corresponds to nth_germline_testing_information
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SCNIRGeneMention.blank_row_fragment()
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
        for i in range(MINIMUM_VARIANTS - 1, MAXIMUM_VARIANTS):
            if i < len(self.variants):
                yield from self.variants[i].to_row_fragment()
            else:
                yield from SCNIRVariant.blank_row_fragment()
        for _ in range(total_closing_cells):
            yield None

    @staticmethod
    def blank_row_fragment() -> Iterable[None]:
        # Hard-coded weirdness - hopefully only for now
        total_opening_cells = 18
        total_closing_cells = 2
        # this is where we would have to do some logic with self.sample_source
        for _ in range(total_opening_cells):
            yield None
        # nth_germline_information
        # total_open_internal_cells = 3
        for _ in range(3):
            yield None
        for i in range(MINIMUM_VARIANTS - 1, MAXIMUM_VARIANTS):
            yield from SCNIRVariant.blank_row_fragment()
        for _ in range(total_closing_cells):
            yield None


@dataclass
class SCNIRForm:
    mrn: int
    gene_mentions: list[SCNIRGeneMention]

    # This is handled downstream
    # def __post_init__(self) -> None:
    #     if (
    #         len(self.gene_mentions) < MINIMUM_GERMLINES
    #         or len(self.gene_mentions) > MAXIMUM_GERMLINES
    #     ):
    #         raise ValueError(
    #             f"Too many or too few gene/germline mentions {len(self.gene_mentions)}"
    #         )

    def to_row(self) -> Iterable[str | bool | None]:
        # Testing information
        # Subject ID
        yield self.mrn
        # Repeat... Repeat... Date last updated
        for _ in range(3):
            yield None
        # Were any potentially disease causing variants found
        yield "Yes" if len(self.gene_mentions) > 0 else None
        # How many GENES with germline variants...
        yield (
            len(self.gene_mentions)
            if len(self.gene_mentions) <= MAXIMUM_GERMLINES
            else f"More than {MAXIMUM_GERMLINES}"
        )
        # Other...
        # this never has anything
        yield None
        for i in range(MINIMUM_GERMLINES - 1, MAXIMUM_GERMLINES):
            if i < len(self.gene_mentions):
                yield from self.gene_mentions[i].to_row_fragment()
            else:
                yield from SCNIRGeneMention.blank_row_fragment()
        # Form level information (deal with this later)
        yield None
        yield None

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCHEMA, orient="row")
