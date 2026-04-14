from collections.abc import Collection

from variants import SupportedVariantTypes
from variants.sources import TextSource

from .filename import parse_date_from_filename, parse_type_from_filename


def attributes_to_text_source(sentence: str, section: str, filename: str) -> TextSource:
    return TextSource(
        filename=filename,
        section=section,
        sentence=sentence,
        file_date=parse_date_from_filename(filename),
    )


# Don't worry about whether a section header contains "somatic" since the only case
# identified is primarily germline, i.e.:
# RELEVANT Potential Germline with Possible Secondary Somatic Variants||POTENTIAL\s+GERMLINE\s+\(WITH\s+POSSIBLE\s+SECONDARY\s+SOMATIC\)\s+VARIANTS\s+\(SBDS,\s+TERC,\s+TERT,\s+DKC[L1],\s+AND\s+DDX41\s+ONLY\)\*\:\s*
def attributes_to_variant_type(
    corpus: str,
    filename: str,
    section: str,
    valid_corpora: Collection[str] = ("sds", "scnir"),
) -> SupportedVariantTypes:
    report_or_note_type = parse_type_from_filename(filename)
    if corpus not in valid_corpora:
        raise ValueError(f"Corpus {corpus} not in supported corpora {valid_corpora}.")
    is_germline = "germ" in section.strip().lower()
    is_rapid_heme = "RHP" == report_or_note_type.strip().upper()
    match corpus, is_germline, is_rapid_heme:
        case "sds", True, _:
            return SupportedVariantTypes.SDS_GERMLINE
        case "scnir", True, _:
            return SupportedVariantTypes.SDS_GERMLINE
        case _, False, True:
            return SupportedVariantTypes.SOMATIC
        case "sds", _, False:
            return SupportedVariantTypes.SDS_GERMLINE
        case "scnir", _, False:
            return SupportedVariantTypes.SCNIR_GERMLINE
        case _:
            return SupportedVariantTypes.NA
