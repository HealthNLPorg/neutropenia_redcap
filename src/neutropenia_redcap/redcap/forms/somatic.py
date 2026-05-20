from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass
from itertools import batched, chain
from operator import attrgetter, methodcaller

import polars as pl
from more_itertools import partition

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
    ) -> tuple[Iterable[SomaticVariant], bool]:
        # In case we benefit from ordering later
        all_variants = list(
            chain.from_iterable(map(attrgetter("variants"), gene_mentions))
        )
        empty_variants_iter, non_empty_variants_iter = partition(
            methodcaller("is_non_empty"), all_variants
        )
        non_empty_variants = list(non_empty_variants_iter)
        # Keep most informative variants if they are available
        if len(non_empty_variants) > 0:
            return non_empty_variants, True
        # If only empties retain them for the form logic
        return empty_variants_iter, False

    def to_rows(self) -> Sequence[Sequence[str | bool | None]]:
        variants, variants_are_non_empty = SomaticForm.collect_variants(
            self.gene_mentions
        )
        return [
            list(self._to_row(variants, variants_are_non_empty))
            for variants in batched(
                variants,
                n=MAXIMUM_SOMATIC_VARIANTS,
            )
        ]

    def _to_row(
        self, variants: Sequence[SomaticVariant], variants_are_non_empty: bool
    ) -> Iterable[str | bool | None]:
        # patient_id
        yield self.mrn
        # redcap_repeat_instrument
        yield "somatic_testing_form"
        # redcap_repeat_instance
        yield 1
        # sum_som, 1 == "Yes"
        yield 1 if variants_are_non_empty else None
        # sum_som_num_var
        yield len(variants) if variants_are_non_empty else None
        for variant in up_to_n(variants, n=MAXIMUM_SOMATIC_VARIANTS, fillvalue=None):
            yield from (
                SomaticVariant.blank_row_fragment(SomaticVariant.total_variant_attrs)
                if variant is None or variants_are_non_empty
                else variant.to_row_fragment()
            )
        # TODO figure out form level comments for summary
