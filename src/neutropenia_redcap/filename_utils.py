from itertools import takewhile


def get_pdf_process(fn: str) -> str:
    return fn.split("_")[0]


def get_original_filename(fn: str) -> str:
    return "_".join(fn.split("_")[1:])


def get_mrn(fn: str) -> int:
    return int("".join(takewhile(str.isnumeric, fn.split("_")[1])))
