from abc import abstractmethod
from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass

from neutropenia_redcap.variants.germline.generic import GermlineGeneMention

from ..generic import REDCapForm


@dataclass
class GermlineForm(REDCapForm):
    gene_mentions: Collection[GermlineGeneMention]

    def to_rows(self) -> Sequence[Sequence[str | bool | None]]:
        # Singleton row for now
        return [list(self.to_row())]

    @abstractmethod
    def to_row(self) -> Iterable[str | bool | None]:
        return []
