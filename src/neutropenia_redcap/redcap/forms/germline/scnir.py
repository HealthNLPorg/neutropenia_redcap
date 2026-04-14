from __future__ import annotations

from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from functools import partial
from itertools import chain

import polars as pl
from more_itertools import padded
from variants.germline.scnir import SCNIRGermlineGeneMention

from .generic import GermlineForm

MINIMUM_SCNIR_GERMLINES = 1
MAXIMUM_SCNIR_GERMLINES = 3

MINIMUM_SCNIR_GERMLINE_VARIANTS = 1
MAXIMUM_SCNIR_GERMLINE_VARIANTS = 4


def germline_and_variant_index_to_columns(
    germline_index: int, variant_index: int
) -> Sequence[str]:
    return [
        f"sum_germ_var{variant_index}_cdna_{germline_index}",
        f"sum_germ_var{variant_index}_pro_{germline_index}",
        f"sum_germ_var{variant_index}_acmg_{germline_index}",
        f"sum_germ_var{variant_index}_comment_{germline_index}",
    ]


def germline_index_to_columns(germline_index: int) -> Sequence[str]:
    # Other way around from data labels format
    variant_index_to_columns = partial(
        germline_and_variant_index_to_columns, germline_index
    )
    return [
        f"sum_germ_gene_{germline_index}",
        f"sum_germ_num_var_{germline_index}",
        *chain.from_iterable(
            map(
                variant_index_to_columns,
                range(
                    MINIMUM_SCNIR_GERMLINE_VARIANTS, MAXIMUM_SCNIR_GERMLINE_VARIANTS + 1
                ),
            )
        ),
    ]


SCNIR_GERMLINE_COLUMNS = list(
    chain(
        (
            "patient_id",
            "sum_germ",
            "sum_germ_num_gen",
        ),
        chain.from_iterable(
            germline_index_to_columns(germline_index)
            for germline_index in range(
                MINIMUM_SCNIR_GERMLINES, MAXIMUM_SCNIR_GERMLINES + 1
            )
        ),
    )
)

SCNIR_GERMLINE_SCHEMA = [
    (column_name, pl.String) for column_name in SCNIR_GERMLINE_COLUMNS
]


@dataclass
class SCNIRGermlineForm(GermlineForm):
    gene_mentions: Collection[SCNIRGermlineGeneMention]

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        yield min(len(self.gene_mentions), MAXIMUM_SCNIR_GERMLINES)
        for germline in padded(
            self.gene_mentions, n=MAXIMUM_SCNIR_GERMLINES, fillvalue=None
        ):
            yield from (
                SCNIRGermlineGeneMention.blank_row_fragment()
                if germline is None
                else germline.to_row_fragment()
            )

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCNIR_GERMLINE_SCHEMA, orient="row")
