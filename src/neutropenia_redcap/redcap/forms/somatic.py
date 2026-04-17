from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from functools import partial
from itertools import chain

import polars as pl
from more_itertools import padded

from neutropenia_redcap.variants.somatic import (
    MAXIMUM_SOMATIC_VARIANTS,
    MAXIMUM_SOMATICS,
    MINIMUM_SOMATIC_VARIANTS,
    MINIMUM_SOMATICS,
    SomaticGeneMention,
)

from .generic import REDCapForm


def somatic_and_variant_index_to_columns(
    germline_index: int, variant_index: int
) -> Sequence[str]:
    return [
        f"sum_germ_var{variant_index}_cdna_{germline_index}",
        f"sum_germ_var{variant_index}_pro_{germline_index}",
        f"sum_germ_var{variant_index}_acmg_{germline_index}",
        f"sum_germ_var{variant_index}_comment_{germline_index}",
    ]


def somatic_index_to_columns(somatic_index: int) -> Sequence[str]:
    # Other way around from data labels format
    variant_index_to_columns = partial(
        somatic_and_variant_index_to_columns, somatic_index
    )
    return [
        f"sum_germ_gene_{somatic_index}",
        f"sum_germ_num_var_{somatic_index}",
        *chain.from_iterable(
            map(
                variant_index_to_columns,
                range(
                    MINIMUM_SOMATIC_VARIANTS,
                    MAXIMUM_SOMATIC_VARIANTS + 1,
                ),
            )
        ),
    ]


COLUMNS = list(
    chain(
        (
            "patient_id",
            "sum_germ",
            "sum_germ_num_gen",
        ),
        chain.from_iterable(
            somatic_index_to_columns(somatic_index)
            for somatic_index in range(MINIMUM_SOMATICS, MAXIMUM_SOMATICS + 1)
        ),
    )
)

SCHEMA = [(column_name, pl.String) for column_name in COLUMNS]


@dataclass
class SomaticForm(REDCapForm):
    gene_mentions: Collection[SomaticGeneMention]

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        yield min(len(self.gene_mentions), MAXIMUM_SOMATICS)
        for somatic in padded(self.gene_mentions, n=MAXIMUM_SOMATICS, fillvalue=None):
            yield from (
                SomaticGeneMention.blank_row_fragment()
                if somatic is None
                else somatic.to_row_fragment()
            )

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCHEMA, orient="row")
