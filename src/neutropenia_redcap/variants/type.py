import logging
from functools import cache

logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)
VARIANT_TYPES = [
    "Pathogenic",
    "Likely Pathogenic",
    "Benign",
    "Likely Benign",
    "Uncertain",
]


@cache
def map_variant_type(variant_type: str | None) -> int | None:
    if variant_type is None:
        return None
    normalized_variant_type = " ".join(variant_type.lower().split())
    is_likely = "likely" in normalized_variant_type
    is_benign = "benign" in normalized_variant_type
    is_pathogenic = "patho" in normalized_variant_type
    is_uncertain = (
        "uncertain" in normalized_variant_type
        or "vus" in normalized_variant_type
        or "unknown significance" in normalized_variant_type
    )
    match is_likely, is_benign, is_pathogenic, is_uncertain:
        # Unmodified pathogenic
        case False, False, True, False:
            return 1
        # Likely pathogenic
        case True, False, True, False:
            return 2
        # Unmodified benign
        case False, True, False, False:
            return 3
        # Likely benign
        case True, True, False, False:
            return 4
        # Unmodified uncertain
        case False, False, False, True:
            return 5
    logger.warning("Cannot currently map: %s", variant_type)
    return None
