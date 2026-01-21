from collections.abc import Iterable, Sequence
from functools import lru_cache
from itertools import chain

MINIMUM_GERMLINES = 1
MAXIMUM_GERMLINES = 3

MINIMUM_VARIANTS = 1
MAXIMUM_VARIANTS = 4


def nth_germline_testing_information(n: int) -> list[str]:
    return [
        f"Were the germline variants found by clinical testing? [{n}]",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=single gene sequencing)",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=gene panel sequencing)",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=whole exome sequencing (wes, exome))",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=whole genome sequencing (wgs, genome))",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=mitochondrial dna sequencing)",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=mitochondrial dna deletion/duplication (del/dup) analysis)",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=deletion/duplication (del/dup) analysis (targeted gene array))",
        f"type of clinical and/or research tests in which the germline variants were identified [{n}]: (choice=copy number analysis (chromosomal microarray))",
        f"sample types in which the germline variants were identified [{n}]: (choice=blood)",
        f"sample types in which the germline variants were identified [{n}]: (choice=bone marrow)",
        f"sample types in which the germline variants were identified [{n}]: (choice=skin fibroblasts)",
        f"sample types in which the germline variants were identified [{n}]: (choice=marrow fibroblasts)",
        f"sample types in which the germline variants were identified [{n}]: (choice=saliva)",
        f"sample types in which the germline variants were identified [{n}]: (choice=hair follicles)",
        f"sample types in which the germline variants were identified [{n}]: (choice=other)",
        f"sample types in which the germline variants were identified [{n}]: (choice=unknown)",
        f"Specify other sample type(s) in which the germline variants were identified [{n}]:",
        # Variant stuff
        *nth_germline_information(n),
        *chain.from_iterable(
            nth_germline_parental_variant_id_information(n, variant_id + 1)
            for variant_id in range(MINIMUM_VARIANTS - 1, MAXIMUM_VARIANTS)
        ),
        f"Significance of the germline variant genetic findings [{n}]:",
        f"Has a definitive or probable causative germline variant result been confirmed by clinical testing? [{n}]",
    ]


def nth_germline_information(n: int) -> Iterable[str]:
    return [
        f"Germline variant(s) gene [{n}]:",
        f"Germline variant(s) gene reference transcript [{n}]:",
        f"Number of germline variants in gene [{n}]:",
    ]


def nth_germline_parental_variant_id_information(
    n: int, variant_id: int
) -> Iterable[str]:
    variant = None
    match variant_id:
        case 1:
            variant = "1/maternal allele variant"
        case 2:
            variant = "2/paternal allele variant"
        case _:
            variant = variant_id
    return [
        f"Germline variant {variant} (genomic or mitochondrial/mtDNA) [{n}]:",
        f"Germline variant {variant} (cDNA) [{n}]:",
        f"Germline variant {variant} (protein) [{n}]:",
        f"Was germline variant {variant_id} confirmed by parental genetic testing? [{n}]",
        f"Type of germline variant {variant} [{n}]: (choice=Synonymous (silent))",
        f"Type of germline variant {variant} [{n}]: (choice=Missense (MS))",
        f"Type of germline variant {variant} [{n}]: (choice=Stop gained)",
        f"Type of germline variant {variant} [{n}]: (choice=Small INDEL - frameshift (FS))",
        f"Type of germline variant {variant} [{n}]: (choice=Small INDEL - in-frame)",
        f"Type of germline variant {variant} [{n}]: (choice=Large INDEL)",
        f"Type of germline variant {variant} [{n}]: (choice=Splicing (SPL))",
        f"Type of germline variant {variant} [{n}]: (choice=Chromosomal aneuploidy)",
        f"Type of germline variant {variant} [{n}]: (choice=Mitochondrial)",
        f"Type of germline variant {variant} [{n}]: (choice=Insertion)",
        f"Type of germline variant {variant} [{n}]: (choice=Deletion)",
        f"Type of germline variant {variant} [{n}]: (choice=Other)",
        f"Specify other type of germline variant {variant} [{n}]:",
        f"Germline variant {variant} PolyPhen2 score [{n}]:",
        f"Mitochondrial variant {variant} heteroplasmy (%) [{n}]:",
        f"Germline variant {variant} ACMG classification [{n}]:",
        f"Germline variant {variant} comment [{n}]:",  # Where to put text context etc.
    ]


TESTING_INFORMATION = [
    "Subject ID:",
    "Repeat Instrument",
    "Repeat Instance",
    "Date the Clinical and Research Genetic Sequencing Summary Form was last updated:",
    "Were any potentially disease-causing GERMLINE (including mitochondrial and copy number) variants identified on any clinical or research sequencing tests?",
    "How many GENES with germline variants and/or copy number variants were identified?",
    "Other, please specify",
]

FORM_LEVEL_INFORMATION = [
    "Clinical and Research Sequencing Summary Form Comments",
    "Complete?",
]


@lru_cache
def build_scnir_columns() -> Sequence[str]:
    return list(
        chain(
            TESTING_INFORMATION,
            # Germline mentions 1 - 3
            chain.from_iterable(
                nth_germline_testing_information(n + 1)
                for n in range(MINIMUM_GERMLINES - 1, MAXIMUM_GERMLINES)
            ),
            FORM_LEVEL_INFORMATION,
        )
    )
