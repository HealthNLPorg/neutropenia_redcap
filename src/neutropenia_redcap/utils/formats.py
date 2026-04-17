from enum import StrEnum, auto


class Formats(StrEnum):
    RAW_TSV = "RAW_TSV"
    REDCAP = "REDCAP"
    NA = auto()
