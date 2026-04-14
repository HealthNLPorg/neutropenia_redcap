from abc import ABC
from dataclasses import dataclass
from enum import Enum, auto


class SupportedVariantTypes(Enum):
    SCNIR_GERMLINE = 0
    SDS_GERMLINE = 1
    SOMATIC = 2
    NA = auto()


@dataclass(eq=True, frozen=True)
class Variant(ABC):
    pass


@dataclass(eq=True, frozen=True)
class GeneMention(ABC):
    pass
