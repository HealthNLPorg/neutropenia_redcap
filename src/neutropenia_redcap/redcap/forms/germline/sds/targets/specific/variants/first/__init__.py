from neutropenia_redcap.redcap.forms.germline.sds.targets.specific import (
    Target,
    TargetCollection,
)

FirstAlleleTargets = TargetCollection(
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
            Target(
                radio_button_value=4,
                gene="DNAJC21",
                attributes=frozenset({"c.983+1G>A"}),
            ),
            Target(
                radio_button_value=5,
                gene="EFL1",
                attributes=frozenset({"c2305A>G"}),
            ),
            Target(
                radio_button_value=6,
                gene="EFL1",
                attributes=frozenset({"Deletion", "15q25.2", "(82422596_82456360)"}),
            ),
            Target(
                radio_button_value=7,
                gene="SRP54",
                attributes=frozenset({"c572C>G"}),
            ),
            Target(
                radio_button_value=8,
                gene="SRP54",
                attributes=frozenset({"c676AG>A"}),
            ),
            Target(
                radio_button_value=9,
                gene="SRP54",
                attributes=frozenset({"c349_351del"}),
            ),
            Target(
                radio_button_value=10,
                gene="SAMD9L",
                attributes=frozenset({"c2606T>C"}),
            ),
        }
    )
)

FirstProteinTargets = TargetCollection(
    targets=frozenset(
        {
            Target(
                radio_button_value=0,
                gene="SDBS",
                attributes=frozenset({"IVS+2", "T>C"}),
            ),
            Target(
                radio_button_value=3,
                gene="SDBS",
                attributes=frozenset({"Homozygous", "T>C", "IVS+2"}),
            ),
            Target(
                radio_button_value=1,
                gene="SDBS",
                attributes=frozenset({"Lys26Term"}),
            ),
            Target(
                radio_button_value=4,
                gene="DNAJC21",
                attributes=frozenset({"IVS7", "splice", "donor", "+1"}),
            ),
            Target(
                radio_button_value=5,
                gene="EFL1",
                attributes=frozenset({"Thr1068Ala"}),
            ),
            Target(
                radio_button_value=6,
                gene="EFL1",
                attributes=frozenset({"deletion", "15q25.2", "(82422596_82456360)"}),
            ),
            Target(
                radio_button_value=7,
                gene="SRP54",
                attributes=frozenset({"Thr191Arg"}),
            ),
            Target(
                radio_button_value=8,
                gene="SRP54",
                attributes=frozenset({"p.Gly226Arg"}),
            ),
            Target(
                radio_button_value=9,
                gene="SRP54",
                attributes=frozenset({"Thr117del"}),
            ),
            Target(
                radio_button_value=10,
                gene="SAMD9L",
                attributes=frozenset({"p.Phe869Ser"}),
            ),
        }
    )
)
