from __future__ import annotations

from collections.abc import Collection, Iterable
from dataclasses import dataclass, field

from more_itertools import padded
from redcap.forms.germline.scnir import (
    MAXIMUM_SCNIR_GERMLINE_VARIANTS,
    MINIMUM_SCNIR_GERMLINE_VARIANTS,
)

from .generic import GermlineVariant


@dataclass(eq=True, frozen=True)
class SCNIRGermlineGeneMention:
    gene: str
    variants: Collection[GermlineVariant] = field(compare=False)

    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SCNIRGermlineGeneMention.blank_row_fragment()
        # sum_germ_gene_{germline_index}
        yield self.gene
        # sum_germ_num_var_{germline_index}
        yield min(len(self.variants), MAXIMUM_SCNIR_GERMLINE_VARIANTS)
        for variant in padded(
            self.variants, n=MAXIMUM_SCNIR_GERMLINE_VARIANTS, fillvalue=None
        ):
            yield from (
                variant.to_row_fragment()
                if variant is not None
                else GermlineVariant.blank_row_fragment(
                    total_variant_attrs=GermlineVariant.total_variant_attrs
                )
            )

    @staticmethod
    def blank_row_fragment() -> Iterable[None]:
        yield None
        yield None
        for _ in range(
            MINIMUM_SCNIR_GERMLINE_VARIANTS, MAXIMUM_SCNIR_GERMLINE_VARIANTS + 1
        ):
            yield from GermlineVariant.blank_row_fragment(
                total_variant_attrs=GermlineVariant.total_variant_attrs
            )
