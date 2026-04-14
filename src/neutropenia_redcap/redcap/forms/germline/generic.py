from collections.abc import Collection
from dataclasses import dataclass

from variants.germline.generic import GermlineGeneMention

from ..generic import REDCapForm


@dataclass
class GermlineForm(REDCapForm):
    gene_mentions: Collection[GermlineGeneMention]
