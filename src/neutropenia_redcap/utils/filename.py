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


def get_pdf_process(fn: str) -> str:
    return fn.split("_")[0]


def get_original_filename(fn: str) -> str:
    return "_".join(fn.split("_")[1:])


def get_mrn(fn: str) -> str:  # int |:
    potential_mrn = "".join(takewhile(str.isnumeric, fn.split("_")[1]))
    if potential_mrn.isnumeric():
        return potential_mrn
    else:
        logger.error("Bad fn: %s", fn)
        return fn.split("_")[1].split("-")[0]


@cache
def parse_date_from_filename(filename: str) -> datetime.date | None:
    raw_date = (
        filename.split("-")[-2] if filename.endswith("lab") else filename.split("-")[-1]
    )
    try:
        month, day, year = (int(attr) for attr in raw_date.split("_"))
        return datetime.date(year=year, month=month, day=day)
    except Exception:
        logger.error("Bad filename, no date information: %s", filename)
        return None
