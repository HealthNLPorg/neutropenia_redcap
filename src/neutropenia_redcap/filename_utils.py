import logging
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
    if "".join(takewhile(str.isnumeric, fn.split("_")[1])).isnumeric():
        return "".join(takewhile(str.isnumeric, fn.split("_")[1]))
    else:
        logger.error("Bad fn: %s", fn)
        return fn.split("_")[1].split("-")[0]
