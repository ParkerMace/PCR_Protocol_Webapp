import csv
import io
import json
from pathlib import Path

import pandas as pd
import streamlit as st


st.set_page_config(
    page_title="PCR OT-2 Protocol Builder",
    page_icon="🧬",
    layout="wide",
)

st.title("PCR OT-2 Protocol Builder")
st.caption(
    "Edit PCR setup parameters, preview sample/primer/well maps, validate common issues, "
    "and export an OT-2 protocol plus reaction map."
)


DEFAULT_CONFIG = {
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


WELL_ORDER_96 = [f"{row}{col}" for col in range(1, 13) for row in "ABCDEFGH"]


def parse_sample_table(text: str):
    sample_genes = {}
    reader = csv.DictReader(io.StringIO(text.strip()))
    for row in reader:
        sample = str(row.get("sample", "")).strip()
        genes_text = str(row.get("genes", "")).strip()
        if not sample:
            continue
        genes = [g.strip() for g in genes_text.split(",") if g.strip()]
        sample_genes[sample] = genes
    return sample_genes


def parse_primer_table(text: str):
    primer_wells = {}
    reader = csv.DictReader(io.StringIO(text.strip()))
    for row in reader:
        gene = str(row.get("gene", "")).strip()
        primer_well = str(row.get("primer_well", "")).strip().upper()
        if gene and primer_well:
            primer_wells[gene] = primer_well
    return primer_wells


def build_reaction_rows(sample_genes, primer_wells, replicates, layout_mode):
    rows = []

    dna_source_map = {
        sample: WELL_ORDER_96[i]
        for i, sample in enumerate(sample_genes.keys())
    }

    if layout_mode == "sample_first":
        base = []
        for sample, genes in sample_genes.items():
            for rep in range(1, replicates + 1):
                for gene in genes:
                    base.append((sample, gene, rep))
    else:
        all_genes = []
        for genes in sample_genes.values():
            for gene in genes:
                if gene not in all_genes:
                    all_genes.append(gene)

        base = []
        for gene in all_genes:
            for sample, genes in sample_genes.items():
                if gene in genes:
                    for rep in range(1, replicates + 1):
                        base.append((sample, gene, rep))

    for i, (sample, gene, rep) in enumerate(base):
        pcr_well = WELL_ORDER_96[i] if i < len(WELL_ORDER_96) else "OVERFLOW"
        rows.append({
            "plate": 1,
            "pcr_well": pcr_well,
            "sample": sample,
            "gene": gene,
            "replicate": rep,
            "dna_source_well": dna_source_map.get(sample, ""),
            "primer_source_well": primer_wells.get(gene, "MISSING"),
        })

    return pd.DataFrame(rows)


def validate_config(config, reaction_df):
    warnings = []
    errors = []

    total_rxns = len(reaction_df)
    rxn_vol = (
        config["volumes"]["master_mix_ul"]
        + config["volumes"]["primer_ul"]
        + config["volumes"]["dna_ul"]
    )

    if total_rxns == 0:
        errors.append("No reactions were generated. Check your sample table.")

    if total_rxns > 96:
        errors.append(
            f"{total_rxns} reactions were generated. This starter template supports one 96-well PCR plate."
        )

    missing_primers = sorted(
        set(reaction_df.loc[reaction_df["primer_source_well"] == "MISSING", "gene"].tolist())
    ) if not reaction_df.empty else []
    if missing_primers:
        errors.append(f"Missing primer well assignments for: {', '.join(missing_primers)}")

    if rxn_vol <= 0:
        errors.append("Reaction volume must be greater than zero.")

    if rxn_vol > 50:
        warnings.append(
            f"Reaction volume is {rxn_vol:.1f} µL. Confirm this is appropriate for your PCR labware and seal."
        )

    if config["thermocycler"]["lid_temp_c"] <= 0:
        errors.append("Heated lid temperature must be greater than 0 °C to reduce evaporation risk.")

    if config["thermocycler"]["cycles"] < 1:
        errors.append("PCR cycles must be at least 1.")

    if config["pipetting"]["primer_tip_strategy"] == "grouped":
        warnings.append(
            "Grouped primer transfer saves tips but reuses one tip per primer source. "
            "Use the safer strategy if contamination risk is a concern."
        )

    required_mm = (
        total_rxns
        * config["volumes"]["master_mix_ul"]
        * (1 + config["volumes"]["master_mix_overage_percent"] / 100)
    )

    return errors, warnings, rxn_vol, required_mm


def config_to_protocol_text(config):
    template_path = Path(__file__).with_name("pcr_protocol_template.py")
    if template_path.exists():
        template = template_path.read_text()
    else:
        st.error("pcr_protocol_template.py was not found in the same folder as this app.")
        return ""

    start = template.index("CONFIG = ")
    end = template.index("\n\n\ndef run", start)
    config_block = "CONFIG = " + json.dumps(config, indent=4)
    return template[:start] + config_block + template[end:]


with st.sidebar:
    st.header("Protocol settings")

    replicates = st.number_input("Replicates", min_value=1, max_value=10, value=DEFAULT_CONFIG["replicates"], step=1)

    layout_mode_label = st.radio(
        "Plate layout",
        ["Sample-first", "Gene-first"],
        help=(
            "Sample-first groups all genes for each sample. "
            "Gene-first groups reactions by target gene, which can make primer additions easier to follow."
        ),
    )
    layout_mode = "sample_first" if layout_mode_label == "Sample-first" else "gene_first"

    primer_tip_strategy = st.radio(
        "Primer tip strategy",
        ["safer", "grouped"],
        help="Safer uses a fresh tip for each primer transfer. Grouped uses one tip per primer source.",
    )

    touch_tip = st.checkbox("Touch tip after primer/DNA dispense", value=True)
    mix_after_dna = st.checkbox("Mix after DNA addition", value=False)

    st.divider()
    st.subheader("Volumes")
    master_mix_ul = st.number_input("Master mix per reaction (µL)", min_value=0.0, max_value=200.0, value=23.0, step=0.5)
    primer_ul = st.number_input("Primer mix per reaction (µL)", min_value=0.0, max_value=20.0, value=1.0, step=0.1)
    dna_ul = st.number_input("DNA per reaction (µL)", min_value=0.0, max_value=20.0, value=1.0, step=0.1)
    overage = st.number_input("Master mix overage (%)", min_value=0, max_value=50, value=10, step=1)

    st.divider()
    st.subheader("Thermocycler")
    lid_temp_c = st.number_input("Heated lid temperature (°C)", min_value=0, max_value=120, value=105, step=1)
    cycles = st.number_input("Cycles", min_value=1, max_value=50, value=30, step=1)

tabs = st.tabs(["Samples & primers", "Thermocycler", "Preview & validation", "Export"])


with tabs[0]:
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Sample table")
        sample_text = st.text_area(
            "CSV: sample,genes",
            value='sample,genes\nM29,"ARO8,NIT1,AMD"',
            height=220,
            help='Use one row per sample. Put multiple genes in quotes separated by commas.',
        )

    with col2:
        st.subheader("Primer source table")
        primer_text = st.text_area(
            "CSV: gene,primer_well",
            value="gene,primer_well\nARO8,A1\nNIT1,B1\nAMD,C1",
            height=220,
            help="Each gene should point to the physical well containing its primer mix.",
        )

    sample_genes = parse_sample_table(sample_text)
    primer_wells = parse_primer_table(primer_text)

    st.write("Parsed sample setup")
    st.dataframe(
        pd.DataFrame(
            [{"sample": s, "genes": ", ".join(g)} for s, g in sample_genes.items()]
        ),
        use_container_width=True,
        hide_index=True,
    )

    st.write("Parsed primer setup")
    st.dataframe(
        pd.DataFrame(
            [{"gene": g, "primer_well": w} for g, w in primer_wells.items()]
        ),
        use_container_width=True,
        hide_index=True,
    )


with tabs[1]:
    col1, col2, col3 = st.columns(3)

    with col1:
        initial_denaturation_temp_c = st.number_input("Initial denaturation temp (°C)", value=94, min_value=0, max_value=120)
        initial_denaturation_seconds = st.number_input("Initial denaturation time (sec)", value=90, min_value=0, max_value=600)

    with col2:
        denature_temp_c = st.number_input("Denature temp (°C)", value=94, min_value=0, max_value=120)
        denature_seconds = st.number_input("Denature time (sec)", value=90, min_value=0, max_value=600)
        anneal_temp_c = st.number_input("Anneal temp (°C)", value=56, min_value=0, max_value=120)
        anneal_seconds = st.number_input("Anneal time (sec)", value=45, min_value=0, max_value=600)

    with col3:
        extend_temp_c = st.number_input("Extension temp (°C)", value=59, min_value=0, max_value=120)
        extend_seconds = st.number_input("Extension time (sec)", value=60, min_value=0, max_value=600)
        final_extension_temp_c = st.number_input("Final extension temp (°C)", value=68, min_value=0, max_value=120)
        final_extension_seconds = st.number_input("Final extension time (sec)", value=300, min_value=0, max_value=1800)
        final_hold_temp_c = st.number_input("Final hold temp (°C)", value=4, min_value=0, max_value=25)


config = json.loads(json.dumps(DEFAULT_CONFIG))
config["replicates"] = int(replicates)
config["layout_mode"] = layout_mode
config["volumes"]["master_mix_ul"] = float(master_mix_ul)
config["volumes"]["primer_ul"] = float(primer_ul)
config["volumes"]["dna_ul"] = float(dna_ul)
config["volumes"]["master_mix_overage_percent"] = int(overage)
config["sample_genes"] = sample_genes
config["primer_wells"] = primer_wells
config["thermocycler"].update({
    "lid_temp_c": int(lid_temp_c),
    "initial_denaturation_temp_c": int(initial_denaturation_temp_c),
    "initial_denaturation_seconds": int(initial_denaturation_seconds),
    "denature_temp_c": int(denature_temp_c),
    "denature_seconds": int(denature_seconds),
    "anneal_temp_c": int(anneal_temp_c),
    "anneal_seconds": int(anneal_seconds),
    "extend_temp_c": int(extend_temp_c),
    "extend_seconds": int(extend_seconds),
    "cycles": int(cycles),
    "final_extension_temp_c": int(final_extension_temp_c),
    "final_extension_seconds": int(final_extension_seconds),
    "final_hold_temp_c": int(final_hold_temp_c),
})
config["pipetting"]["primer_tip_strategy"] = primer_tip_strategy
config["pipetting"]["touch_tip"] = bool(touch_tip)
config["pipetting"]["mix_after_dna"] = bool(mix_after_dna)

reaction_df = build_reaction_rows(sample_genes, primer_wells, int(replicates), layout_mode)
errors, warnings, rxn_vol, required_mm = validate_config(config, reaction_df)


with tabs[2]:
    st.subheader("Validation")

    metric_cols = st.columns(4)
    metric_cols[0].metric("Total reactions", len(reaction_df))
    metric_cols[1].metric("Reaction volume", f"{rxn_vol:.1f} µL")
    metric_cols[2].metric("Master mix needed", f"{required_mm:.1f} µL")
    metric_cols[3].metric("Primer strategy", primer_tip_strategy)

    if errors:
        for error in errors:
            st.error(error)
    else:
        st.success("No blocking validation errors found.")

    for warning in warnings:
        st.warning(warning)

    st.info(
        "Evaporation safeguard is enabled in the exported protocol: the robot pauses for plate sealing/capping "
        "before closing the lid and setting the heated lid temperature."
    )

    st.subheader("Reaction map")
    st.dataframe(reaction_df, use_container_width=True, hide_index=True)

    st.subheader("PCR plate grid")
    if not reaction_df.empty:
        grid_data = []
        for row in "ABCDEFGH":
            grid_row = {"row": row}
            for col in range(1, 13):
                well = f"{row}{col}"
                hit = reaction_df[reaction_df["pcr_well"] == well]
                if hit.empty:
                    grid_row[str(col)] = ""
                else:
                    r = hit.iloc[0]
                    grid_row[str(col)] = f"{r['sample']} | {r['gene']} | R{r['replicate']}"
            grid_data.append(grid_row)
        st.dataframe(pd.DataFrame(grid_data), use_container_width=True, hide_index=True)


with tabs[3]:
    st.subheader("Export files")

    protocol_text = config_to_protocol_text(config)
    reaction_csv = reaction_df.to_csv(index=False)
    config_json = json.dumps(config, indent=2)

    st.download_button(
        "Download OT-2 protocol.py",
        data=protocol_text,
        file_name="generated_pcr_protocol.py",
        mime="text/x-python",
        disabled=bool(errors) or not bool(protocol_text),
    )

    st.download_button(
        "Download reaction_map.csv",
        data=reaction_csv,
        file_name="reaction_map.csv",
        mime="text/csv",
    )

    st.download_button(
        "Download config.json",
        data=config_json,
        file_name="pcr_config.json",
        mime="application/json",
    )

    with st.expander("Preview generated CONFIG"):
        st.code(config_json, language="json")
