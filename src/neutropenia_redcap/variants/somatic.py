import logging
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from operator import attrgetter, is_not_none

from more_itertools import partition

from . import GeneMention, Variant
from .sources import TextSource
from .type import VARIANT_TYPES, map_variant_type

MINIMUM_SOMATIC_VARIANTS = 1
MAXIMUM_SOMATIC_VARIANTS = 5  # Per REDCap - caveat is this is total, not per gene

logger = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s -   %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


@dataclass(eq=True, frozen=True)
class SomaticVariant(Variant):
    # Another weird thing is I can't find the field where the specimen collection date
    # would go
    total_variant_attrs = 6

    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        if blank:
            yield from SomaticVariant.blank_row_fragment(
                total_variant_attrs=SomaticVariant.total_variant_attrs
            )
        variant_type, source = self.select_variant_type()
        # sum_som_gene_{variant_index}
        yield self.gene
        # sum_som_cdna_{variant_index}
        yield self.syntax_n
        # sum_som_pro_{variant_index}
        yield self.syntax_p
        # sum_som_path_{variant_index}
        yield variant_type
        # sum_som_vaf_{variant_index}
        yield self.vaf
        # sum_som_comment_{variant_index}
        yield self.build_comment(variant_type=variant_type, source=source)

    # Heterozygosity in "Clinical and Research Sequencing Summary Form Comments"
    # or variant level comments
    def build_comment(self, variant_type: int | None, source: str) -> str:
        mention_summary = self.build_mention_summary(
            variant_type=variant_type, source=source
        )
        source_guide = self.build_source_guide()
        return mention_summary + source_guide

    def build_mention_summary(self, variant_type: int | None, source: str) -> str:
        normalized_syntax_n = (
            self.syntax_n.lower() if self.syntax_n is not None else "None found"
        )
        normalized_syntax_p = (
            self.syntax_p.lower() if self.syntax_p is not None else "None found"
        )
        normalized_vaf = self.vaf.lower() if self.vaf is not None else "None found"
        match self.heterozygous:
            case False:
                heterozygous = "No"
            case True:
                heterozygous = "Yes"
            case _:
                heterozygous = "Unknown"
        normalized_variant_type = (
            VARIANT_TYPES[variant_type - 1] if variant_type is not None else "Unknown"
        )
        return (
            "Mention Summary:\n"
            f"Gene: {self.gene.upper()} Nucleotide Syntax: {normalized_syntax_n} Protein Syntax: {normalized_syntax_p} VAF: {normalized_vaf} Variant Type (Parsed from {source.title()}): {normalized_variant_type} Heterozygous: {heterozygous}\n\n"
        )

    def select_variant_type(self) -> tuple[int | None, str]:
        model_parsed_variant_type = map_variant_type(self.variant_type)
        if model_parsed_variant_type is not None:
            return (
                model_parsed_variant_type
                if model_parsed_variant_type % 2
                else model_parsed_variant_type - 1,
                "sentence",
            )
        return self.select_variant_type_from_sources()

    def select_variant_type_from_sources(self) -> tuple[int | None, str]:
        sections = list(map(attrgetter("section"), self.text_sources))
        mapped_sections = list(
            filter(is_not_none, map(map_variant_type, sections)),
        )
        match len(mapped_sections):
            case 0:
                variant_type = None
            case 1:
                variant_type = mapped_sections[0]
            case _:
                if all(vtype == mapped_sections[0] for vtype in mapped_sections):
                    variant_type = mapped_sections[0]
                else:
                    logger.error(
                        "Variant type inconsistencies with sections, using None: %s",
                        ", ".join(sections),
                    )
                    variant_type = None
        return variant_type, "sections"

    def build_source_guide(self) -> str:
        match len(self.text_sources):
            case 0:
                raise ValueError("Mention has no text sources")
            case 1:
                return SomaticVariant.build_single_source_guide(
                    text_source=next(iter(self.text_sources))
                )
            case _:
                return SomaticVariant.build_multi_source_guide(
                    text_sources=self.text_sources
                )

    @staticmethod
    def build_single_source_guide(text_source: TextSource) -> str:
        return f"Found in sentence:\n{text_source.sentence}\nFrom file {text_source.filename} - section '{text_source.section}'"

    @staticmethod
    def build_multi_source_guide(text_sources: Collection[TextSource]) -> str:
        no_date_info, has_date_info = partition(
            lambda source: source.file_date is not None, text_sources
        )
        sorted_sources = sorted(has_date_info, key=attrgetter("file_date"))
        others = (
            f"From file {text_source.filename} - section '{text_source.section}'"
            for text_source in sorted_sources[1:]
        )
        proper = f"Multiple sources found, showing earliest:\n{SomaticVariant.build_single_source_guide(sorted_sources[0])}\n\nOthers least to most recent:\n{'\n'.join(others)}\n"
        no_date_info = list(no_date_info)
        if len(no_date_info) > 0:
            no_date_info_mentions = "\n".join(
                f"From file {text_source.filename} - section '{text_source.section}'"
                for text_source in no_date_info
            )
            return f"{proper}Others with no date information:\n{no_date_info_mentions}"
        return proper

    @staticmethod
    def blank_row_fragment(total_variant_attrs: int) -> Iterable[None]:
        for _ in range(total_variant_attrs):
            yield None


@dataclass(eq=True, frozen=True)
class SomaticGeneMention(GeneMention):
    def to_row_fragment(self, blank: bool = False) -> Iterable[str | bool | None]:
        raise ValueError(
            "Do not use this method, somatic forms use variants rather than full mentions"
        )
