from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from itertools import chain
from operator import attrgetter

import polars as pl
from more_itertools import padded

from neutropenia_redcap.variants.somatic import (
    MAXIMUM_SOMATIC_VARIANTS,
    MINIMUM_SOMATIC_VARIANTS,
    SomaticGeneMention,
    SomaticVariant,
)

from .generic import REDCapForm


def variant_index_to_columns(variant_index: int) -> Sequence[str]:
    return [
        f"sum_som_cdna_{variant_index}",
        f"sum_som_pro_{variant_index}",
        f"sum_som_acmg_{variant_index}",
        f"sum_som_comment_{variant_index}",
    ]


COLUMNS = list(
    chain(
        (
            "patient_id",
            # "Were any variants identified?"
            "sum_som",
            # "How many?"
            "sum_som_num_var",
        ),
        chain.from_iterable(
            variant_index_to_columns(variant_index=variant_index)
            for variant_index in range(
                MINIMUM_SOMATIC_VARIANTS,
                MAXIMUM_SOMATIC_VARIANTS + 1,
            )
        ),
    )
)

SCHEMA = [(column_name, pl.String) for column_name in COLUMNS]


@dataclass
class SomaticForm(REDCapForm):
    gene_mentions: Collection[SomaticGeneMention]

    @staticmethod
    def collect_variants(
        gene_mentions: Collection[SomaticGeneMention],
    ) -> Sequence[SomaticVariant]:
        # In case we benefit from ordering later
        return list(chain.from_iterable(map(attrgetter("variants"), gene_mentions)))

    def to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # sum_germ, 1 == "Yes"
        yield 1
        # sum_germ_num_gen
        variants = SomaticForm.collect_variants(self.gene_mentions)
        yield min(len(variants), MAXIMUM_SOMATIC_VARIANTS)
        for variant in padded(variants, n=MAXIMUM_SOMATIC_VARIANTS, fillvalue=None):
            yield from (
                SomaticVariant.blank_row_fragment(SomaticVariant.total_variant_attrs)
                if variant is None
                else variant.to_row_fragment()
            )

    def to_data_frame(self) -> pl.DataFrame:
        data = [list(self.to_row())]
        return pl.DataFrame(data=data, schema=SCHEMA, orient="row")
