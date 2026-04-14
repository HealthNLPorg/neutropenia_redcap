from variants.sources import TextSource

from .filename import parse_date_from_filename


def attributes_to_text_source(sentence: str, section: str, filename: str) -> TextSource:
    return TextSource(
        filename=filename,
        section=section,
        sentence=sentence,
        file_date=parse_date_from_filename(filename),
    )
