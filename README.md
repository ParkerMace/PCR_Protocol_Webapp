# PCR OT-2 Protocol Builder Template

This starter package contains:

- `pcr_protocol_template.py` - OT-2 PCR setup and thermocycler template.
- `pcr_protocol_gui.py` - Streamlit web app for editing parameters, previewing reaction maps, validating common issues, and exporting a generated protocol.
- `requirements.txt` - Python packages for the GUI.

## How to run the GUI

From this folder:

```bash
pip install -r requirements.txt
streamlit run pcr_protocol_gui.py
```

The app opens in your browser.

## Recommended workflow

1. Enter your sample table.
2. Enter your primer source table.
3. Pick sample-first or gene-first layout.
4. Confirm volumes and thermocycler settings.
5. Check the validation tab.
6. Download:
   - `generated_pcr_protocol.py`
   - `reaction_map.csv`
   - `pcr_config.json`

## What changed from the original PCR script

The exported OT-2 protocol includes:

- A generated reaction map in protocol comments.
- Master mix overage calculation.
- Missing primer checks.
- One-plate 96-reaction validation.
- Manual pause to seal/cap the PCR plate before thermocycling.
- Heated-lid setting before cycling.
- Safer primer transfer option using a fresh tip per transfer.
- Optional grouped primer transfer if you want fewer tips.

The template intentionally supports one 96-well PCR plate to avoid silent multi-plate assignment errors. Extend deliberately if you need multi-plate mode.


## PCR presets added

The GUI now includes editable PCR starting presets:

- Original script / current lab starting point
- NEB Taq / routine PCR starting point
- NEB OneTaq / colony PCR starting point
- NEB Q5 high-fidelity starting point
- Two-step PCR starting point for high-Tm primers

These are starting points, not final assay guarantees. Always adjust annealing temperature based on primer Tm and extension time based on expected amplicon length/polymerase recommendations.
