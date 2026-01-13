from dataclasses import dataclass
from datetime import date

import polars as pl


@dataclass
class SCNIRVariant:
    syntax_p: str | None
    syntax_n: str | None
    variant_type: str | None
    vaf: str | None
    known_heterozygous: bool
    specimen_collection_date: date | None
    sample_source: str | None


@dataclass
class SCNIRGeneMention:
    variants: list[SCNIRVariant]


@dataclass
class SCNIRForm:
    mrn: int
    gene_mentions: list[SCNIRGeneMention]

    def to_data_frame(self) -> pl.DataFrame:
        return pl.DataFrame()
