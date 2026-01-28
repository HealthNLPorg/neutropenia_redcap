from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date

import polars as pl

from .redcap_import import (
    MAXIMUM_GERMLINES,
    MAXIMUM_VARIANTS,
    MINIMUM_GERMLINES,
    MINIMUM_VARIANTS,
    SCNIR_COLUMNS,
)

SCNIR_SCHEMA = [(column_name, pl.String) for column_name in SCNIR_COLUMNS]


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
        # sum_germ_var{variant_index}_cdna_{germline_index}
        yield self.syntax_n
        # sum_germ_var{variant_index}_pro_{germline_index}
        yield self.syntax_p
        # sum_germ_var{variant_index}_acmg_{germline_index}
        yield self.variant_type
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
    variants: list[SCNIRVariant]

    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SCNIRGeneMention.blank_row_fragment()
        # sum_germ_gene_{germline_index}
        yield self.gene
        # sum_germ_num_var_{germline_index}
        yield min(len(self.variants), MAXIMUM_VARIANTS)
        for i in range(MINIMUM_VARIANTS, MAXIMUM_VARIANTS + 1):
            if i <= len(self.variants):
                yield from self.variants[i - 1].to_row_fragment()
            else:
                yield from SCNIRVariant.blank_row_fragment()

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
        for i in range(MINIMUM_GERMLINES, MAXIMUM_GERMLINES + 1):
            if i <= len(self.gene_mentions):
                yield from self.gene_mentions[i - 1].to_row_fragment()
            else:
                yield from SCNIRGeneMention.blank_row_fragment()

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCNIR_SCHEMA, orient="row")
