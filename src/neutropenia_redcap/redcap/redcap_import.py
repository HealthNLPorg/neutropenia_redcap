from collections.abc import Sequence
from functools import partial
from itertools import chain

MINIMUM_GERMLINES = 1
MAXIMUM_GERMLINES = 3

MINIMUM_VARIANTS = 1
MAXIMUM_VARIANTS = 4


def germline_and_variant_index_to_columns(
    germline_index: int, variant_index: int
) -> Sequence[str]:
    return [
        f"sum_germ_var{variant_index}_cdna_{germline_index}",
        f"sum_germ_var{variant_index}_pro_{germline_index}",
        f"sum_germ_var{variant_index}_acmg_{germline_index}",
        f"sum_germ_var{variant_index}_comment_{germline_index}",
    ]


def germline_index_to_columns(germline_index: int) -> Sequence[str]:
    # Other way around from data labels format
    variant_index_to_columns = partial(
        germline_and_variant_index_to_columns, germline_index
    )
    return [
        f"sum_germ_gene_{germline_index}",
        f"sum_germ_num_var_{germline_index}",
        *chain.from_iterable(
            map(variant_index_to_columns, range(MINIMUM_VARIANTS, MAXIMUM_VARIANTS + 1))
        ),
    ]


SCNIR_COLUMNS = [
    "patient_id",
    "sum_germ",
    "sum_germ_num_gen",
    *chain.from_iterable(
        germline_index_to_columns(germline_index)
        for germline_index in range(MINIMUM_GERMLINES, MAXIMUM_GERMLINES + 1)
    ),
]
