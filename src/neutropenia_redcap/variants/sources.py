import datetime
from dataclasses import dataclass


@dataclass(eq=True, frozen=True)
class TextSource:
    filename: str
    section: str
    sentence: str
    file_date: datetime.date | None
