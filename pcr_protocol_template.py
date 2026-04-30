from opentrons import protocol_api
import math
from collections import defaultdict

metadata = {
    "protocolName": "User_Configurable_PCR_Setup",
    "author": "Generated from PCR Protocol GUI",
    "description": "Automated PCR setup with master mix, primer mixes, DNA samples, reaction map, seal pause, and thermocycler heated lid.",
}

requirements = {
    "robotType": "OT-2",
    "apiLevel": "2.16",
}

# ======================================================================================
# USER CONFIGURATION
# The GUI replaces this CONFIG block when exporting a protocol.
# You can also edit it manually.
# ======================================================================================

CONFIG = {
    "replicates": 1,
    "layout_mode": "sample_first",

    "volumes": {
        "master_mix_ul": 23.0,
        "primer_ul": 1.0,
        "dna_ul": 1.0,
        "master_mix_overage_percent": 10,
    },

    "thermocycler": {
        "lid_temp_c": 105,
        "initial_denaturation_temp_c": 94,
        "initial_denaturation_seconds": 90,
        "denature_temp_c": 94,
        "denature_seconds": 90,
        "anneal_temp_c": 56,
        "anneal_seconds": 45,
        "extend_temp_c": 59,
        "extend_seconds": 60,
        "cycles": 30,
        "final_extension_temp_c": 68,
        "final_extension_seconds": 300,
        "final_hold_temp_c": 4,
    },

    "sample_genes": {
        "M29": ["ARO8", "NIT1", "AMD"],
    },

    "primer_wells": {
        "ARO8": "A1",
        "NIT1": "B1",
        "AMD": "C1",
    },

    "labware": {
        "pcr_plate": "opentrons_96_aluminumblock_generic_pcr_strip_200ul",
        "primer_rack": "opentrons_96_aluminumblock_generic_pcr_strip_200ul",
        "dna_plate": "opentrons_96_aluminumblock_generic_pcr_strip_200ul",
        "master_mix_tuberack": "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap",
        "p300_tiprack": "opentrons_96_tiprack_300ul",
        "p20_tiprack": "opentrons_96_tiprack_20ul",
    },

    "deck": {
        "master_mix_tuberack": "2",
        "p20_tiprack": "3",
        "dna_plate": "4",
        "primer_rack": "5",
        "p300_tiprack": "6",
    },

    "source_wells": {
        "master_mix_tube": "A1",
    },

    "pipetting": {
        "primer_tip_strategy": "safer",
        "dna_tip_strategy": "always",
        "touch_tip": True,
        "mix_after_dna": False,
        "mix_repetitions": 3,
        "mix_volume_ul": 10,
    },
}


def run(protocol: protocol_api.ProtocolContext):
    """Run a validated PCR setup and thermocycler program."""

    tc_mod = protocol.load_module(module_name="thermocyclerModuleV1")
    pcr_plate = tc_mod.load_labware(CONFIG["labware"]["pcr_plate"])

    master_mix_tuberack = protocol.load_labware(
        CONFIG["labware"]["master_mix_tuberack"],
        CONFIG["deck"]["master_mix_tuberack"],
    )
    primer_rack = protocol.load_labware(
        CONFIG["labware"]["primer_rack"],
        CONFIG["deck"]["primer_rack"],
    )
    dna_plate = protocol.load_labware(
        CONFIG["labware"]["dna_plate"],
        CONFIG["deck"]["dna_plate"],
    )

    p300_tiprack = protocol.load_labware(
        CONFIG["labware"]["p300_tiprack"],
        CONFIG["deck"]["p300_tiprack"],
    )
    p20_tiprack = protocol.load_labware(
        CONFIG["labware"]["p20_tiprack"],
        CONFIG["deck"]["p20_tiprack"],
    )

    p300 = protocol.load_instrument("p300_single_gen2", "left", tip_racks=[p300_tiprack])
    p20 = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[p20_tiprack])

    master_mix_tube = master_mix_tuberack.wells_by_name()[CONFIG["source_wells"]["master_mix_tube"]]

    volumes = CONFIG["volumes"]
    thermo = CONFIG["thermocycler"]
    sample_genes = CONFIG["sample_genes"]
    primer_wells = CONFIG["primer_wells"]

    dna_sources = {}
    for i, sample in enumerate(sample_genes.keys()):
        dna_sources[sample] = dna_plate.wells()[i]

    primer_sources = {
        gene: primer_rack.wells_by_name()[well_name]
        for gene, well_name in primer_wells.items()
    }

    reaction_list = build_reaction_list(
        sample_genes=sample_genes,
        replicates=CONFIG["replicates"],
        layout_mode=CONFIG["layout_mode"],
    )

    reaction_assignments = []
    for i, rxn in enumerate(reaction_list):
        dest = pcr_plate.wells()[i]
        sample = rxn["sample"]
        gene = rxn["gene"]
        replicate = rxn["replicate"]
        reaction_assignments.append({
            "dest": dest,
            "well_name": dest.display_name.split(" ")[0],
            "sample": sample,
            "gene": gene,
            "replicate": replicate,
            "dna_source": dna_sources[sample],
            "primer_source": primer_sources[gene],
        })

    validate_run(protocol, reaction_assignments, sample_genes, primer_sources, volumes, thermo)

    target_wells = [r["dest"] for r in reaction_assignments]

    protocol.comment("=== PCR SETUP SUMMARY ===")
    protocol.comment(f"Total reactions: {len(reaction_assignments)}")
    protocol.comment(f"Layout mode: {CONFIG['layout_mode']}")
    reaction_volume = volumes["master_mix_ul"] + volumes["primer_ul"] + volumes["dna_ul"]
    protocol.comment(f"Reaction volume: {reaction_volume:.1f} uL")

    required_mm = (
        len(reaction_assignments)
        * volumes["master_mix_ul"]
        * (1 + volumes["master_mix_overage_percent"] / 100)
    )
    protocol.comment(
        f"Prepare at least {required_mm:.1f} uL master mix "
        f"including {volumes['master_mix_overage_percent']}% overage."
    )

    protocol.comment("=== SOURCE MAP ===")
    protocol.comment(f"Master mix: {master_mix_tube.display_name}")
    for sample, well in dna_sources.items():
        protocol.comment(f"DNA sample {sample}: {well.display_name}")
    for gene, well in primer_sources.items():
        protocol.comment(f"Primer mix {gene}: {well.display_name}")

    protocol.comment("=== REACTION MAP ===")
    for r in reaction_assignments:
        protocol.comment(
            f"{r['dest'].display_name}: sample={r['sample']} | gene={r['gene']} | "
            f"replicate={r['replicate']} | DNA={r['dna_source'].display_name} | "
            f"primer={r['primer_source'].display_name}"
        )

    tc_mod.open_lid()

    distribute_master_mix(
        protocol=protocol,
        pipette=p300,
        source=master_mix_tube,
        dest_wells=target_wells,
        volume_ul=volumes["master_mix_ul"],
    )

    add_primers(
        protocol=protocol,
        pipette=p20,
        reactions=reaction_assignments,
        volume_ul=volumes["primer_ul"],
        strategy=CONFIG["pipetting"]["primer_tip_strategy"],
        touch_tip=CONFIG["pipetting"]["touch_tip"],
    )

    add_dna(
        protocol=protocol,
        pipette=p20,
        reactions=reaction_assignments,
        volume_ul=volumes["dna_ul"],
        touch_tip=CONFIG["pipetting"]["touch_tip"],
        mix_after=CONFIG["pipetting"]["mix_after_dna"],
        mix_repetitions=CONFIG["pipetting"]["mix_repetitions"],
        mix_volume_ul=CONFIG["pipetting"]["mix_volume_ul"],
    )

    protocol.pause(
        "Seal or cap the PCR plate now. Confirm every well is fully sealed before continuing. "
        "This pause is included to reduce evaporation risk during thermocycling."
    )

    run_pcr(protocol, tc_mod, thermo)

    protocol.comment("Protocol complete.")


def build_reaction_list(sample_genes, replicates, layout_mode):
    """Return a list of dicts: sample, gene, replicate."""
    if layout_mode not in {"sample_first", "gene_first"}:
        raise RuntimeError("layout_mode must be 'sample_first' or 'gene_first'.")

    reactions = []

    if layout_mode == "sample_first":
        for sample, genes in sample_genes.items():
            for rep in range(1, replicates + 1):
                for gene in genes:
                    reactions.append({"sample": sample, "gene": gene, "replicate": rep})

    if layout_mode == "gene_first":
        all_genes_in_order = []
        for genes in sample_genes.values():
            for gene in genes:
                if gene not in all_genes_in_order:
                    all_genes_in_order.append(gene)

        for gene in all_genes_in_order:
            for sample, genes in sample_genes.items():
                if gene in genes:
                    for rep in range(1, replicates + 1):
                        reactions.append({"sample": sample, "gene": gene, "replicate": rep})

    return reactions


def validate_run(protocol, reaction_assignments, sample_genes, primer_sources, volumes, thermo):
    """Raise errors for problems that should be fixed before the robot moves."""

    if not reaction_assignments:
        raise RuntimeError("No reactions were generated. Check sample and gene inputs.")

    if len(reaction_assignments) > 96:
        raise RuntimeError(
            f"{len(reaction_assignments)} reactions were generated, but this template supports one 96-well PCR plate. "
            "Reduce reactions or extend the template to multi-plate mode."
        )

    missing_primer_genes = []
    for sample, genes in sample_genes.items():
        for gene in genes:
            if gene not in primer_sources:
                missing_primer_genes.append(gene)

    if missing_primer_genes:
        raise RuntimeError(f"Missing primer well assignments for: {sorted(set(missing_primer_genes))}")

    reaction_volume = volumes["master_mix_ul"] + volumes["primer_ul"] + volumes["dna_ul"]
    if reaction_volume <= 0:
        raise RuntimeError("Reaction volume must be greater than zero.")

    if reaction_volume > 50:
        protocol.comment(
            f"WARNING: reaction volume is {reaction_volume:.1f} uL. "
            "Confirm this is appropriate for the selected PCR labware and seal."
        )

    if thermo["lid_temp_c"] is None or thermo["lid_temp_c"] <= 0:
        raise RuntimeError("A heated lid temperature is required to reduce evaporation risk.")

    if thermo["cycles"] < 1:
        raise RuntimeError("PCR cycles must be at least 1.")


def distribute_master_mix(protocol, pipette, source, dest_wells, volume_ul):
    """Distribute master mix with refill chunks that do not exceed pipette capacity."""
    protocol.comment("Adding master mix.")
    pipette.pick_up_tip()

    max_aspirate = pipette.max_volume * 0.90
    wells_per_aspirate = max(1, int(max_aspirate // volume_ul))

    for start in range(0, len(dest_wells), wells_per_aspirate):
        chunk = dest_wells[start:start + wells_per_aspirate]
        total_volume = volume_ul * len(chunk)

        pipette.aspirate(total_volume, source)
        for dest in chunk:
            pipette.dispense(volume_ul, dest)
        pipette.blow_out(source.top())

    pipette.drop_tip()


def add_primers(protocol, pipette, reactions, volume_ul, strategy, touch_tip=True):
    """Add primer mixes with either safest or grouped tip strategy."""
    protocol.comment("Adding primer mixes.")

    if strategy == "safer":
        for r in reactions:
            pipette.transfer(
                volume_ul,
                r["primer_source"],
                r["dest"],
                new_tip="always",
                touch_tip=touch_tip,
            )
        return

    if strategy == "grouped":
        grouped = defaultdict(list)
        for r in reactions:
            grouped[r["gene"]].append(r)

        for gene, gene_reactions in grouped.items():
            source = gene_reactions[0]["primer_source"]
            protocol.comment(f"Adding primer mix for {gene}.")
            pipette.pick_up_tip()
            for r in gene_reactions:
                pipette.aspirate(volume_ul, source)
                pipette.dispense(volume_ul, r["dest"])
                if touch_tip:
                    pipette.touch_tip(r["dest"])
            pipette.drop_tip()
        return

    raise RuntimeError("primer_tip_strategy must be 'safer' or 'grouped'.")


def add_dna(protocol, pipette, reactions, volume_ul, touch_tip=True, mix_after=False, mix_repetitions=3, mix_volume_ul=10):
    """Add DNA with a fresh tip for every sample/reaction."""
    protocol.comment("Adding DNA samples.")

    for r in reactions:
        pipette.transfer(
            volume_ul,
            r["dna_source"],
            r["dest"],
            new_tip="always",
            touch_tip=touch_tip,
            mix_after=(mix_repetitions, mix_volume_ul) if mix_after else None,
        )


def run_pcr(protocol, tc_mod, thermo):
    """Run thermocycler program with heated lid enabled."""
    protocol.comment("Starting thermocycler program.")

    tc_mod.close_lid()
    tc_mod.set_lid_temperature(thermo["lid_temp_c"])

    if thermo.get("initial_denaturation_seconds", 0) > 0:
        protocol.comment("Initial denaturation.")
        tc_mod.set_block_temperature(
            temperature=thermo["initial_denaturation_temp_c"],
            hold_time_seconds=thermo["initial_denaturation_seconds"],
        )

    for cycle in range(1, thermo["cycles"] + 1):
        protocol.comment(f"Cycle {cycle} of {thermo['cycles']}")

        tc_mod.set_block_temperature(
            temperature=thermo["denature_temp_c"],
            hold_time_seconds=thermo["denature_seconds"],
        )

        tc_mod.set_block_temperature(
            temperature=thermo["anneal_temp_c"],
            hold_time_seconds=thermo["anneal_seconds"],
        )

        tc_mod.set_block_temperature(
            temperature=thermo["extend_temp_c"],
            hold_time_seconds=thermo["extend_seconds"],
        )

    protocol.comment("Final extension.")
    tc_mod.set_block_temperature(
        temperature=thermo["final_extension_temp_c"],
        hold_time_seconds=thermo["final_extension_seconds"],
    )

    protocol.comment(f"Final hold at {thermo['final_hold_temp_c']} °C.")
    tc_mod.set_block_temperature(thermo["final_hold_temp_c"])
