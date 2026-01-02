from enum import IntEnum


class Formats(IntEnum):
    RAW_TSV = 0
    REDCAP = 1


valid_format_choices = [_format.name for _format in Formats]
