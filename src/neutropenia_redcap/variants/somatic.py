from collections.abc import Collection
from dataclasses import dataclass, field
from datetime import date

from .sources import TextSource


@dataclass(eq=True, frozen=True)
class SomaticVariant:
    gene: str
    syntax_p: str | None
    syntax_n: str | None
    variant_type: str | None
    vaf: str | None
    heterozygous: (
        bool | None
    )  # True for is heterozygous, False for definitely isn't, None for unknown
    text_sources: Collection[TextSource] = field(compare=False)
    specimen_collection_dates: Collection[date] = field(compare=False)
    sample_sources: Collection[str] = field(compare=False)

    # TODO a lot of the potential methods might be coordinated with the
    # germline variant code
