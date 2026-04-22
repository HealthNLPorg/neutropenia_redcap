from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from itertools import chain
from operator import attrgetter

import polars as pl

from neutropenia_redcap.utils.iter import up_to_n
from neutropenia_redcap.variants.somatic import (
    MAXIMUM_SOMATIC_VARIANTS,
    MINIMUM_SOMATIC_VARIANTS,
    SomaticGeneMention,
    SomaticVariant,
)

from .generic import REDCapForm


def variant_index_to_columns(variant_index: int) -> Sequence[str]:
    return [
        f"sum_som_gene_{variant_index}",
        f"sum_som_cdna_{variant_index}",
        f"sum_som_pro_{variant_index}",
        f"sum_som_path_{variant_index}",
        f"sum_som_vaf_{variant_index}",
        f"sum_som_comment_{variant_index}",
    ]


_SOMATIC_COLUMNS = chain(
    (
        "patient_id",
        # Form type (in this case "somatic_testing_form")
        "redcap_repeat_instrument",
        # On our end will be 1
        # unless maybe with using the REDCap API
        # we can check whether a given MRN somatic testing
        # form is already populated
        "redcap_repeat_instance",
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


@dataclass
class SomaticForm(REDCapForm):
    gene_mentions: Collection[SomaticGeneMention]
    schema = [(column_name, pl.String) for column_name in _SOMATIC_COLUMNS]

    @staticmethod
    def collect_variants(
        gene_mentions: Collection[SomaticGeneMention],
    ) -> Sequence[SomaticVariant]:
        # In case we benefit from ordering later
        return list(chain.from_iterable(map(attrgetter("variants"), gene_mentions)))

    def to_row(self) -> Sequence[Sequence[str | bool | None]]:
        return [list(self._to_row())]

    def _to_row(self) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # redcap_repeat_instrument
        yield "somatic_testing_form"
        # redcap_repeat_instance
        yield 1
        # sum_som, 1 == "Yes"
        yield 1
        # sum_som_num_var
        variants = SomaticForm.collect_variants(self.gene_mentions)
        yield min(len(variants), MAXIMUM_SOMATIC_VARIANTS)
        for variant in up_to_n(variants, n=MAXIMUM_SOMATIC_VARIANTS, fillvalue=None):
            yield from (
                SomaticVariant.blank_row_fragment(SomaticVariant.total_variant_attrs)
                if variant is None
                else variant.to_row_fragment()
            )
