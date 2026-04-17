from neutropenia_redcap.variants import SupportedVariantTypes
from neutropenia_redcap.variants.sources import TextSource

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
    filename: str,
    section: str,
) -> SupportedVariantTypes:
    report_or_note_type = parse_type_from_filename(filename)
    is_germline = "germ" in section.strip().lower()
    is_rapid_heme = "RHP" == report_or_note_type.strip().upper()
    match is_germline, is_rapid_heme:
        case True, _:
            return SupportedVariantTypes.GERMLINE
        case False, True:
            return SupportedVariantTypes.SOMATIC
        case _, False:
            return SupportedVariantTypes.GERMLINE
        case _:
            return SupportedVariantTypes.NA
