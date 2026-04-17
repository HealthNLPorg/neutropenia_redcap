from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto

import polars as pl


class SupportedFormTypes(Enum):
    SCNIR_GERMLINE = 0
    SDS_GERMLINE = 1
    SOMATIC = 2
    NA = auto()


@dataclass
class REDCapForm(ABC):
    mrn: int

    @abstractmethod
    def to_data_frame(self) -> pl.DataFrame:
        return pl.DataFrame([])
