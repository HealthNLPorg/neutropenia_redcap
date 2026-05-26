from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from functools import partial
from itertools import chain

import polars as pl

from neutropenia_redcap.utils.iter import up_to_n
from neutropenia_redcap.variants.germline.sds import (
    MAXIMUM_SDS_GERMLINE_VARIANTS,
    MAXIMUM_SDS_GERMLINES,
    MINIMUM_SDS_GERMLINE_VARIANTS,
    MINIMUM_SDS_GERMLINES,
    SDSGermlineGeneMention,
)

from ..generic import (
    GermlineForm,
)


def germline_and_variant_index_to_columns(
    germline_index: int, variant_index: int
) -> Sequence[str]:
    return [
        f"sum_germ_var{variant_index}_cdna_{germline_index}",
        f"sum_germ_var{variant_index}_pro_{germline_index}",
        f"sum_germ_var{variant_index}_acmg_{germline_index}",
        f"sum_germ_var{variant_index}_comment_{germline_index}",
    ]


def germline_index_to_columns(germline_index: int) -> Iterable[str]:
    # Other way around from data labels format
    variant_index_to_columns = partial(
        germline_and_variant_index_to_columns, germline_index
    )
    return chain(
        (f"sum_germ_gene_{germline_index}", f"sum_germ_num_var_{germline_index}"),
        chain.from_iterable(
            map(
                variant_index_to_columns,
                range(MINIMUM_SDS_GERMLINE_VARIANTS, MAXIMUM_SDS_GERMLINE_VARIANTS + 1),
            )
        ),
    )


_SDS_GERMLINE_COLUMNS = chain(
    (
        "patient_id",
        "sum_germ",
        "sum_germ_num_gen",
    ),
    chain.from_iterable(
        germline_index_to_columns(germline_index)
        for germline_index in range(MINIMUM_SDS_GERMLINES, MAXIMUM_SDS_GERMLINES + 1)
    ),
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
