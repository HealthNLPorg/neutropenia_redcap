import datetime
import logging
from functools import cache
from itertools import takewhile

logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


@cache
def get_pdf_process(filename: str) -> str:
    return filename.split("_")[0]


@cache
def get_original_filename(filename: str) -> str:
    return "_".join(filename.split("_")[1:])


@cache
def get_mrn(filename: str) -> str:  # int |:
    potential_mrn = "".join(takewhile(str.isnumeric, filename.split("_")[1]))
    if potential_mrn.isnumeric():
        return potential_mrn
    else:
        logger.error("Bad filename: %s", filename)
        return filename.split("_")[1].split("-")[0]


@cache
def parse_type_from_filename(filename: str) -> str:
    try:
        return filename.split("-")[1]
    except Exception:
        raise ValueError(f"Bad filename: {filename}")


@cache
def parse_date_from_filename(filename: str) -> datetime.date | None:
    raw_date = (
        filename.split("-")[-2] if filename.endswith("lab") else filename.split("-")[-1]
    )
    try:
        month, day, year = (int(attr) for attr in raw_date.split("_"))
        return datetime.date(year=year, month=month, day=day)
    except Exception:
        logger.error(
            "Bad filename, either no date information or incorrectly formatted: %s",
            filename,
        )
        return None
