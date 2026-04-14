from variants.sources import TextSource

from .filename import parse_date_from_filename


def attributes_to_text_source(sentence: str, section: str, filename: str) -> TextSource:
    return TextSource(
        filename=filename,
        section=section,
        sentence=sentence,
        file_date=parse_date_from_filename(filename),
    )


# def attributes_to_variant_type(corpus: str, filename: str, section: str, valid_corpora: Collection[str] = ("sds", "scnir")) -> REDCapForm:
#     report_or_note_type = parse_type_from_filename(filename)
#     if corpus not in valid_corpora:
#         raise ValueError(f"Corpus {corpus} not in supported corpora {valid_corpora}.")
#     is_germline = "germ" in section.strip().lower()
#     is_rapid_heme = "RHP" == report_or_note_type.strip().upper()
#     match corpus, is_germline, is_rapid_heme:
#         case "sds", True, _:
#             return SDSGermlineForm()
#         case "scnir", True, _:
#             return SCNIRGermlineForm
#         case _, False, True:
#             return SomaticForm
