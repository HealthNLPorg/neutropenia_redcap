from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum, auto
from typing import Any, ClassVar

import polars as pl

PolarsDataType = pl.datatypes.classes.DataTypeClass


class SupportedFormTypes(Enum):
    SCNIR_GERMLINE = 0
    SDS_GERMLINE = 1
    SOMATIC = 2
    NA = auto()


@dataclass
class REDCapForm(ABC):
    mrn: int
    schema: ClassVar[Sequence[tuple[str, PolarsDataType]]]

    @abstractmethod
    def to_rows(self) -> Sequence[Sequence[Any]]:
        return []

    def to_data_frame(self) -> pl.DataFrame:
        try:
            return pl.DataFrame(data=self.to_rows(), schema=self.schema, orient="row")
        except Exception:
            print(self.schema)
            print(self.to_rows()[0])
            raise ValueError(f"{len(self.schema)} {len(self.to_rows()[0])}")
