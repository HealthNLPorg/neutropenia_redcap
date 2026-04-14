from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from functools import partial
from itertools import chain

import polars as pl
from more_itertools import padded
from variants.germline.sds import SDSGermlineGeneMention

from .generic import GermlineForm

MINIMUM_SDS_GERMLINES = 1
MAXIMUM_SDS_GERMLINES = 3

MINIMUM_SDS_GERMLINE_VARIANTS = 1
MAXIMUM_SDS_GERMLINE_VARIANTS = 4


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
                    MINIMUM_SDS_GERMLINE_VARIANTS, MAXIMUM_SDS_GERMLINE_VARIANTS + 1
                ),
            )
        ),
    ]


SDS_GERMLINE_COLUMNS = list(
    chain(
        (
            "patient_id",
            "sum_germ",
            "sum_germ_num_gen",
        ),
        chain.from_iterable(
            germline_index_to_columns(germline_index)
            for germline_index in range(
                MINIMUM_SDS_GERMLINES, MAXIMUM_SDS_GERMLINES + 1
            )
        ),
    )
)

SDS_GERMLINE_SCHEMA = [
    (column_name, pl.String) for column_name in SDS_GERMLINE_COLUMNS
]


@dataclass
class SDSGermlineForm(GermlineForm):
    form_name = "sds_germline"
    gene_mentions: Collection[SDSGermlineGeneMention]

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        yield min(len(self.gene_mentions), MAXIMUM_SDS_GERMLINES)
        for germline in padded(
            self.gene_mentions, n=MAXIMUM_SDS_GERMLINES, fillvalue=None
        ):
            yield from (
                SDSGermlineGeneMention.blank_row_fragment()
                if germline is None
                else germline.to_row_fragment()
            )

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SDS_GERMLINE_SCHEMA, orient="row")
