# Membrane Transport Trace Classifier

This tool classifies particle transport traces across the nuclear membrane using trajectory data from **ImageJ** and a user-defined membrane reference. It segments, analyzes, and labels particle movement (e.g., import, export, docking, abortive transport) based on spatial landmarks along the nuclear pore complex (NPC).

---

## 📦 Features

- Supports both **import** and **export** directionality
- Classifies traces as:
  - `successful`
  - `docking_event`
  - `abortive_central_channel`
  - `abortive_cytoplasmic_fibril`
  - `abortive_nuclear_basket`
  - `incomplete`
- Calculates:
  - Regional travel distances
  - Dwell time in non-terminal regions
  - Mean step distance in homogeneous regions
- Outputs:
  - Per-trace classification
  - Trace visualizations (optional)
  - Regional statistics

---

## 📁 Inputs

1. **Membrane File (CSV)**  
   Exported from ImageJ. Contains line segments describing membrane boundaries across six regions:
   - `bsk` – basal scaffold
   - `nucChnl` – nuclear channel
   - `midline` – nuclear midline
   - `cytChnl` – cytoplasmic channel
   - `npcCyto` – cytoplasmic edge

   **Format:** A CSV with 10 columns, each pair representing X/Y for each segment:
   ```
   bsk_X, bsk_Y, nucChnl_X, nucChnl_Y, ..., npcCyto_Y
   ```

2. **Trace File (CSV)**  
   A list of particle positions over time generated by ImageJ's particle tracker.

   **Expected columns:**
   - `Frame` – Frame number
   - `X`, `Y` – Coordinates of the particle

---

## 🚀 Usage

```bash
python controller.py <membrane.csv> <traces.csv> <pixel_size> [OPTIONS]
```

### Required arguments:
- `<membrane.csv>`: Path to the membrane boundary CSV
- `<traces.csv>`: Path to the particle trajectory CSV
- `<pixel_size>`: Scalar for converting pixel to nanometer scale

### Optional arguments:
| Flag | Description |
|------|-------------|
| `--imp`, `-i` | Set direction as **import** |
| `--exp`, `-e` | Set direction as **export** |
| `--normalize`, `-n` | Normalize coordinates by pixel size |
| `--plot`, `-p` | Save trace plots for complete traces |
| `--plot_incomplete`, `-pi` | Save plots for incomplete traces |
| `--frame_time`, `-ft` | Specify frame duration in ms for dwell time reporting |

---

## 📊 Output

Outputs are saved under:

```bash
./output/<trace_filename>_import/  # or _export/
```

Directory contents:
- `out.csv` — Classification and metrics per trace
- `traces/` — Visualizations of categorized traces:
  - `successful/`
  - `docking_event/`
  - `abortive_central_channel/`
  - `abortive_cytoplasmic_fibril/`
  - `abortive_nuclear_basket/`
  - `incomplete/`
  - `midline/` — Traces that crossed the NPC midline
- `AllTraces.png` — Combined plot of all complete traces

---

## 📦 Dependencies

- `numpy`
- `pandas`
- `matplotlib`
- `argparse`
- Python ≥ 3.7

Install dependencies via pip:
```bash
pip install numpy pandas matplotlib
```

---

## 🧠 Classification Logic Summary

- Traces are segmented based on whether they span defined NPC subregions.
- Direction is user-specified (`import` or `export`).
- Traces are pruned to include only continuous and valid segments.
- Based on the deepest region reached (e.g., `nuclear_basket` vs `central_scaffold1`), traces are labeled as abortive, successful, or docking.

---

## 📝 Notes

- ImageJ must be used to:
  - Define membrane segments for boundary input
  - Export particle tracking coordinates

Ensure coordinate systems and scale (pixel size) are consistent across inputs.


