from collections.abc import Iterable
from dataclasses import dataclass

from neutropenia_redcap.utils.iter import up_to_n

from .generic import GermlineGeneMention, GermlineVariant

MINIMUM_SDS_GERMLINES = 1
MAXIMUM_SDS_GERMLINES = 3

MINIMUM_SDS_GERMLINE_VARIANTS = 1
MAXIMUM_SDS_GERMLINE_VARIANTS = 4

MINIMUM_SDS_ALLELES = 1
MAXIMUM_SDS_ALLELES = 2


@dataclass(eq=True, frozen=True)
class SDSGermlineGeneMention(GermlineGeneMention):
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SDSGermlineGeneMention.blank_row_fragment()
        # sum_germ_gene_{germline_index}
        yield self.gene
        # sum_germ_num_var_{germline_index}
        yield min(len(self.variants), MAXIMUM_SDS_GERMLINE_VARIANTS)
        for variant in up_to_n(
            self.variants, n=MAXIMUM_SDS_GERMLINE_VARIANTS, fillvalue=None
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
            MINIMUM_SDS_GERMLINE_VARIANTS, MAXIMUM_SDS_GERMLINE_VARIANTS + 1
        ):
            yield from GermlineVariant.blank_row_fragment(
                total_variant_attrs=GermlineVariant.total_variant_attrs
            )
