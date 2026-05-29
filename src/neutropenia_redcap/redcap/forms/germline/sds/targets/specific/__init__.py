from collections.abc import Set
from enum import IntEnum

from attrs import Factory, define

from neutropenia_redcap.variants.germline.generic import GermlineVariant


class TargetNotFound(IntEnum):
    other = 2
    unknown = 98
    none_detected = 99


@define(frozen=True)
class Target:
    radio_button_value: int
    gene: str
    attributes: Set[str] = Factory(frozenset)

    def match_germline_variant(
        self, germline_variant: GermlineVariant, literal=False
    ) -> bool:
        if literal:

            def _normalize(s: str) -> str:
                return s
        else:

            def _normalize(s: str) -> str:
                return " ".join(s.strip().split()).lower()

        return _normalize(self.gene) == _normalize(germline_variant.gene) and all(
            germline_variant.search_attrs(keyword=attribute, literal=literal)
            for attribute in self.attributes
        )


@define(frozen=True)
class TargetCollection:
    targets: Set[Target] = Factory(frozenset)

    def match_target(self, germline_variant: GermlineVariant, literal=False) -> int:
        matches = {
            target
            for target in self.targets
            if target.match_germline_variant(
                germline_variant=germline_variant, literal=literal
            )
        }
        match len(matches):
            case 0:
                return TargetNotFound.other.value
            case 1:
                return next(iter(matches)).radio_button_value
            case _:
                raise ValueError(
                    f"Match for germline variant {germline_variant} should be unique, found {matches}"
                )
