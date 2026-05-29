from neutropenia_redcap.redcap.forms.germline.sds.targets.specific import (
    Target,
    TargetCollection,
)

FirstProteinTargets = TargetCollection(
    targets=frozenset(
        {
            Target(
                radio_button_value=0, gene="SDBS", attributes=frozenset({"c.258+2T>C"})
            ),
            Target(
                radio_button_value=3,
                gene="SDBS",
                attributes=frozenset({"c.258+2T>C", "Homozygous"}),
            ),
            Target(
                radio_button_value=1,
                gene="SDBS",
                attributes=frozenset({"c.183_184TA>CT"}),
            ),
        }
    )
)
