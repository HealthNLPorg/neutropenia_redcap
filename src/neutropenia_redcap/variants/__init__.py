from abc import ABC
from dataclasses import dataclass


@dataclass(eq=True, frozen=True)
class Variant(ABC):
    raise NotImplementedError


@dataclass(eq=True, frozen=True)
class GeneMention(ABC):
    raise NotImplementedError
