from abc import ABC, abstractmethod
from dataclasses import dataclass

import polars as pl


@dataclass
class REDCapForm(ABC):
    mrn: int
    form_name: str

    @abstractmethod
    def to_data_frame(self) -> pl.DataFrame:
        raise NotImplementedError("Need to implement for your particular form")
