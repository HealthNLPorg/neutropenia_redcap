from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable
from dataclasses import dataclass, field
from datetime import date
from enum import Enum, auto
from typing import ClassVar

from .sources import TextSource


class SupportedVariantTypes(Enum):
    GERMLINE = 0
    SOMATIC = 1
    NA = auto()


class SupportedMentionTypes(Enum):
    GERMLINE = 0
    SOMATIC = 1
    NA = auto()


@dataclass(eq=True, frozen=True)
class Variant(ABC):
    gene: str
    syntax_p: str | None
    syntax_n: str | None
    variant_type: str | None
    vaf: str | None
    heterozygous: (
        bool | None
    )  # True for is heterozygous, False for definitely isn't, None for unknown
    text_sources: Collection[TextSource] = field(compare=False)
    specimen_collection_dates: Collection[date] = field(compare=False)
    sample_sources: Collection[str] = field(compare=False)
    # protein syntax, nucleotide syntax, variant type, comment
    total_variant_attrs: ClassVar[int] = 4

    @abstractmethod
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        return []


@dataclass(eq=True, frozen=True)
class GeneMention(ABC):
    gene: str
    variants: Collection[Variant] = field(compare=False)

    @abstractmethod
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        return []
