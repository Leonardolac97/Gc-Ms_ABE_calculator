# ██████████████████████████████████████████████████████████████████████
# █  GC-MS calculator for ABE mixtures reaction                        █
# █                                                                    █
# █  AI-assisted: OpenAI ChatGPT (GPT-5.3)                             █
# █  04.2026                                                           █
# ██████████████████████████████████████████████████████████████████████


# -----------------------------
# Libraries
# -----------------------------

import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd

from PyQt6.QtCore import QUrl
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib as mpl


# -----------------------------
# Matplotlib global styling
# -----------------------------
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 14


# -------------------------------------------------------
# Target products + aliases --> to be changed accordingly
# -------------------------------------------------------
TARGET_PRODUCTS: Dict[str, List[str]] = {
    "2-Pentanone": ["2-pentanone"],
    "2-Heptanone": ["2-heptanone"],
    "4-Nonanone": ["4-nonanone"],
    "6-Undecanone": ["6-undecanone"],
    "2-Heptanol": ["2-heptanol"],
    "6-Undecanol": ["6-undecanol"],
    "2-Ethyl-1-hexanol": ["2-ethyl-1-hexanol", "2-ethylhexanol"],
    "Isopropanol": ["isopropanol", "isopropyl alcohol"],
    "Butyl butyrate": ["butyl butyrate", "butanoic acid, butyl ester"],
    "n-Butyl acetate": ["n-butyl acetate", "butyl acetate", "acetic acid, butyl ester"],
    "3-Hepten-2-one": ["3-hepten-2-one"],
    "3,4-Dimethyl-3-penten-2-one": [
        "3,4-dimethyl-3-penten-2-one",
        "3-penten-2-one, 3,4-dimethyl-",
    ],
}


def parse_carbon_number_from_formula(formula: str) -> Optional[int]:
    if not isinstance(formula, str):
        return None
    m = re.search(r"C(\d+)", formula.strip())
    return int(m.group(1)) if m else None


def _norm(s: str) -> str:
    return str(s).strip().lower()


def extract_reactant_area_percents(df: pd.DataFrame) -> Dict[str, float]:
    """
    Returns dict with keys: ethanol, acetone, butanol
    Uses Family == Reactants first; else falls back to Compound Name matching.
    Sums Area % if multiple rows match.
    """
    out = {"ethanol": 0.0, "acetone": 0.0, "butanol": 0.0}

    df = df.copy()
    df["Area %"] = pd.to_numeric(df.get("Area %", 0.0), errors="coerce").fillna(0.0)
    df["Compound Name"] = df.get("Compound Name", "").astype(str)
    df["Family"] = df.get("Family", "").astype(str)

    react = df[df["Family"].astype(str).str.strip().str.lower() == "reactants"].copy()
    if not react.empty:
        for _, r in react.iterrows():
            name = _norm(r["Compound Name"])
            if "ethanol" in name:
                out["ethanol"] += float(r["Area %"])
            elif "acetone" in name:
                out["acetone"] += float(r["Area %"])
            elif "1-butanol" in name or name == "butanol" or ("butanol" in name and "2-" not in name):
                out["butanol"] += float(r["Area %"])
        return out

    for _, r in df.iterrows():
        name = _norm(r["Compound Name"])
        if "ethanol" in name:
            out["ethanol"] += float(r["Area %"])
        elif "acetone" in name:
            out["acetone"] += float(r["Area %"])
        elif "1-butanol" in name or name == "butanol" or ("butanol" in name and "2-" not in name):
            out["butanol"] += float(r["Area %"])

    return out


def compute_product_selectivity(
    df_keep_products: pd.DataFrame,
    target_products: Dict[str, List[str]],
) -> pd.Series:
    """
    Computes selectivity (%) for selected products based on the same basis
    used for family/carbon selectivity:
      - Classification == KEEP
      - Family != Reactants

    Matching is done through aliases in Compound Name.
    """
    if df_keep_products.empty:
        return pd.Series(dtype=float)

    df = df_keep_products.copy()
    df["Compound Name"] = df["Compound Name"].astype(str)
    df["Area %"] = pd.to_numeric(df["Area %"], errors="coerce").fillna(0.0)

    basis_total = float(df["Area %"].sum())
    if basis_total <= 0:
        return pd.Series(dtype=float)

    compound_names_norm = df["Compound Name"].str.strip().str.lower()
    results = {}

    for product_label, aliases in target_products.items():
        aliases_norm = [a.strip().lower() for a in aliases]
        mask = compound_names_norm.apply(lambda x: any(alias in x for alias in aliases_norm))
        area_sum = float(df.loc[mask, "Area %"].sum())
        results[product_label] = area_sum / basis_total * 100.0

    return pd.Series(results, dtype=float)


@dataclass
class ProcessedData:
    df_keep_products: pd.DataFrame
    fam_sel: pd.Series
    c_sel: pd.Series
    bin_sel: pd.Series
    product_sel: pd.Series
    basis_total_area: float
    reactants_at_area: Dict[str, float]


def load_and_process_csv(path: str) -> ProcessedData:
    df = pd.read_csv(path)

    required = {"RT", "Compound Name", "Formula", "Area %", "Classification", "Notes", "Family"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

    if "Carbon Number" in df.columns:
        df["Carbon Number"] = pd.to_numeric(df["Carbon Number"], errors="coerce")
    else:
        df["Carbon Number"] = df["Formula"].apply(parse_carbon_number_from_formula)

    df["Area %"] = pd.to_numeric(df["Area %"], errors="coerce").fillna(0.0)
    df["Classification"] = df["Classification"].astype(str).str.strip().str.upper()
    df["Family"] = df["Family"].astype(str).str.strip()

    reactants_at = extract_reactant_area_percents(df)

    df_keep = df[df["Classification"] == "KEEP"].copy()
    df_keep_products = df_keep[df_keep["Family"].str.lower() != "reactants"].copy()

    if df_keep_products.empty:
        return ProcessedData(
            df_keep_products=df_keep_products,
            fam_sel=pd.Series(dtype=float),
            c_sel=pd.Series(dtype=float),
            bin_sel=pd.Series(dtype=float),
            product_sel=pd.Series(dtype=float),
            basis_total_area=0.0,
            reactants_at_area=reactants_at,
        )

    basis_total = float(df_keep_products["Area %"].sum())

    fam_sum = df_keep_products.groupby("Family", dropna=False)["Area %"].sum().sort_values(ascending=False)
    fam_sel = fam_sum / basis_total * 100.0

    df_c = df_keep_products.dropna(subset=["Carbon Number"]).copy()
    df_c["Carbon Number"] = df_c["Carbon Number"].astype(int)

    c_sum = df_c.groupby("Carbon Number")["Area %"].sum().sort_index()
    c_sel = c_sum / basis_total * 100.0

    def bin_label(cn: int) -> Optional[str]:
        if 2 <= cn <= 6:
            return "C2–C6"
        if 7 <= cn <= 15:
            return "C7–C15"
        return None

    df_c["Cbin"] = df_c["Carbon Number"].apply(bin_label)
    df_bins = df_c.dropna(subset=["Cbin"])
    bin_sum = df_bins.groupby("Cbin")["Area %"].sum()
    bin_sum = bin_sum.reindex(["C2–C6", "C7–C15"]).fillna(0.0)
    bin_sel = bin_sum / basis_total * 100.0

    product_sel = compute_product_selectivity(df_keep_products, TARGET_PRODUCTS)

    return ProcessedData(
        df_keep_products=df_keep_products,
        fam_sel=fam_sel,
        c_sel=c_sel,
        bin_sel=bin_sel,
        product_sel=product_sel,
        basis_total_area=basis_total,
        reactants_at_area=reactants_at,
    )


def export_all_csvs(output_dir: Path, data: ProcessedData, conversion_df: Optional[pd.DataFrame]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    data.df_keep_products.to_csv(output_dir / "keep_products_table.csv", index=False)

    fam_df = data.fam_sel.rename("Selectivity_%").to_frame().reset_index().rename(columns={"index": "Family"})
    fam_df["Basis_Total_Area%"] = data.basis_total_area
    fam_df.to_csv(output_dir / "family_selectivity.csv", index=False)

    c_df = data.c_sel.rename("Selectivity_%").to_frame().reset_index().rename(columns={"index": "Carbon_Number"})
    c_df["Basis_Total_Area%"] = data.basis_total_area
    c_df.to_csv(output_dir / "carbon_selectivity.csv", index=False)

    b_df = data.bin_sel.rename("Selectivity_%").to_frame().reset_index().rename(columns={"index": "Carbon_Bin"})
    b_df["Basis_Total_Area%"] = data.basis_total_area
    b_df.to_csv(output_dir / "carbon_bins_selectivity.csv", index=False)

    p_df = data.product_sel.rename("Selectivity_%").to_frame().reset_index().rename(columns={"index": "Product"})
    p_df["Basis_Total_Area%"] = data.basis_total_area
    p_df.to_csv(output_dir / "product_selectivity.csv", index=False)

    if conversion_df is not None:
        conversion_df.to_csv(output_dir / "conversion.csv", index=False)


class PlotPanel(QWidget):
    def __init__(self, parent: QWidget = None):
        super().__init__(parent)
        self.fig = Figure(constrained_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.fig.add_subplot(111)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def clear(self):
        self.ax.clear()

    def draw(self):
        self.canvas.draw_idle()


# ----------------
# USER INTERFACE
# ----------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GC Selectivity + Conversion Viewer (KEEP products only)")

        self.data: Optional[ProcessedData] = None
        self.current_path: Optional[Path] = None
        self.output_dir: Optional[Path] = None

        self.baseline_a0 = {"ethanol": None, "acetone": None, "butanol": None}

        self.btn_load = QPushButton("Load sample CSV")
        self.btn_load.clicked.connect(self.load_sample_csv)

        self.btn_load_baseline = QPushButton("Load baseline (0 h)")
        self.btn_load_baseline.clicked.connect(self.load_baseline_csv)

        self.cmap_label = QLabel("Color scheme:")
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems([
            "plasma", "tab10", "Set2", "Set3", "Paired",
            "viridis", "cividis",
            "Blues", "Greens", "Oranges", "Purples"
        ])
        self.cmap_combo.currentTextChanged.connect(self.redraw_all)

        self.btn_open_folder = QPushButton("Open results folder")
        self.btn_open_folder.clicked.connect(self.open_results_folder)
        self.btn_open_folder.setEnabled(False)

        self.a0_eth = QLineEdit()
        self.a0_ace = QLineEdit()
        self.a0_but = QLineEdit()
        for w in (self.a0_eth, self.a0_ace, self.a0_but):
            w.setFixedWidth(90)
            w.setPlaceholderText("A0")

        self.a0_eth.editingFinished.connect(self._read_manual_a0)
        self.a0_ace.editingFinished.connect(self._read_manual_a0)
        self.a0_but.editingFinished.connect(self._read_manual_a0)

        top_bar = QHBoxLayout()
        top_bar.addWidget(self.btn_load)
        top_bar.addWidget(self.btn_load_baseline)
        top_bar.addSpacing(15)
        top_bar.addWidget(QLabel("Initial A% (optional): Ethanol"))
        top_bar.addWidget(self.a0_eth)
        top_bar.addWidget(QLabel("Acetone"))
        top_bar.addWidget(self.a0_ace)
        top_bar.addWidget(QLabel("1-Butanol"))
        top_bar.addWidget(self.a0_but)
        top_bar.addSpacing(15)
        top_bar.addWidget(self.cmap_label)
        top_bar.addWidget(self.cmap_combo)
        top_bar.addSpacing(15)
        top_bar.addWidget(self.btn_open_folder)
        top_bar.addStretch(1)

        self.tabs = QTabWidget()
        self.tab_family = PlotPanel()
        self.tab_carbon = PlotPanel()
        self.tab_bins = PlotPanel()
        self.tab_products = PlotPanel()
        self.tab_conv = PlotPanel()

        self.tabs.addTab(self.tab_family, "Family")
        self.tabs.addTab(self.tab_carbon, "Carbon #")
        self.tabs.addTab(self.tab_bins, "C-bins")
        self.tabs.addTab(self.tab_products, "Products")
        self.tabs.addTab(self.tab_conv, "Conversion")

        root = QWidget()
        root_layout = QVBoxLayout()
        root_layout.addLayout(top_bar)
        root_layout.addWidget(self.tabs)
        root.setLayout(root_layout)
        self.setCentralWidget(root)

        self.redraw_all()

    def _read_manual_a0(self):
        def parse_float(x: str) -> Optional[float]:
            x = str(x).strip()
            if not x:
                return None
            try:
                return float(x)
            except ValueError:
                return None

        self.baseline_a0["ethanol"] = parse_float(self.a0_eth.text())
        self.baseline_a0["acetone"] = parse_float(self.a0_ace.text())
        self.baseline_a0["butanol"] = parse_float(self.a0_but.text())
        self.redraw_all()
        self._export_if_ready()

    def load_baseline_csv(self):
        path_str, _ = QFileDialog.getOpenFileName(self, "Select baseline (0 h) CSV", "", "CSV files (*.csv);;All files (*.*)")
        if not path_str:
            return
        try:
            df0 = pd.read_csv(path_str)
            a0 = extract_reactant_area_percents(df0)
            self.baseline_a0 = {k: float(v) for k, v in a0.items()}

            self.a0_eth.setText(f"{self.baseline_a0['ethanol']:.6g}")
            self.a0_ace.setText(f"{self.baseline_a0['acetone']:.6g}")
            self.a0_but.setText(f"{self.baseline_a0['butanol']:.6g}")

            self.redraw_all()
            self._export_if_ready()

        except Exception as e:
            QMessageBox.critical(self, "Failed to load baseline CSV", str(e))

    def load_sample_csv(self):
        path_str, _ = QFileDialog.getOpenFileName(self, "Select sample CSV", "", "CSV files (*.csv);;All files (*.*)")
        if not path_str:
            return

        path = Path(path_str)

        try:
            self.data = load_and_process_csv(str(path))
            self.current_path = path

            sample_name = path.stem
            self.output_dir = path.parent / f"{sample_name}_results"
            self.output_dir.mkdir(parents=True, exist_ok=True)

            self.btn_open_folder.setEnabled(True)

        except Exception as e:
            QMessageBox.critical(self, "Failed to load sample CSV", str(e))
            self.data = None
            self.current_path = None
            self.output_dir = None
            self.btn_open_folder.setEnabled(False)
            return

        self.redraw_all()
        self._export_if_ready()
        self.open_results_folder()

    def open_results_folder(self):
        if not self.output_dir:
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.output_dir)))

    def _bar_colors(self, n: int):
        cmap_name = self.cmap_combo.currentText()
        cmap = mpl.colormaps.get(cmap_name)
        if n <= 1:
            return [cmap(0.6)]
        return [cmap(i / (n - 1)) for i in range(n)]

    def _style_axes(self, ax):
        ax.tick_params(axis="both", labelsize=14)
        ax.set_title("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def redraw_all(self):
        self.draw_family_plot()
        self.draw_carbon_plot()
        self.draw_bins_plot()
        self.draw_products_plot()
        self.draw_conversion_plot()

    def compute_conversion_df(self) -> Optional[pd.DataFrame]:
        if not self.data:
            return None

        at = self.data.reactants_at_area
        a0 = self.baseline_a0

        rows = []
        for key, label in [("ethanol", "Ethanol"), ("acetone", "Acetone"), ("butanol", "1-Butanol")]:
            A0 = a0.get(key)
            At = float(at.get(key, 0.0))
            if A0 is None or A0 <= 0:
                conv = None
            else:
                conv = 100.0 * (A0 - At) / A0
            rows.append({"Reactant": label, "A0_Area%": A0, "At_Area%": At, "Conversion_%": conv})

        return pd.DataFrame(rows)

    def draw_conversion_plot(self):
        panel = self.tab_conv
        panel.clear()

        conv_df = self.compute_conversion_df()
        if conv_df is None:
            panel.ax.text(0.5, 0.5, "Load a sample CSV first.", ha="center", va="center", transform=panel.ax.transAxes)
            self._style_axes(panel.ax)
            panel.draw()
            return

        if conv_df["Conversion_%"].isna().all():
            panel.ax.text(
                0.5, 0.5,
                "Conversion needs initial A% (A0).\nLoad baseline (0 h) or type A0 values above.",
                ha="center", va="center", transform=panel.ax.transAxes
            )
            self._style_axes(panel.ax)
            panel.draw()
            return

        plot_df = conv_df.dropna(subset=["Conversion_%"]).copy()
        labels = plot_df["Reactant"].tolist()
        values = plot_df["Conversion_%"].astype(float).tolist()

        colors = self._bar_colors(len(values))
        panel.ax.bar(labels, values, color=colors)
        panel.ax.set_ylabel("Conversion (%)")
        panel.ax.set_xlabel("")
        panel.ax.set_ylim(0, max(5, float(max(values)) * 1.15))
        self._style_axes(panel.ax)
        panel.draw()

    def _export_if_ready(self):
        if not (self.data and self.output_dir):
            return

        conv_df = self.compute_conversion_df()
        if conv_df is not None and conv_df["Conversion_%"].notna().any():
            export_all_csvs(self.output_dir, self.data, conv_df)
        else:
            export_all_csvs(self.output_dir, self.data, None)

    def draw_family_plot(self):
        panel = self.tab_family
        panel.clear()

        if not self.data or self.data.fam_sel.empty:
            panel.ax.text(
                0.5, 0.5,
                "Load a CSV to plot.\n(KEEP products only; Reactants/FLAG/DISCARD ignored)",
                ha="center", va="center", transform=panel.ax.transAxes
            )
            self._style_axes(panel.ax)
            panel.draw()
            return

        s = self.data.fam_sel
        labels = list(s.index)
        values = s.values
        colors = self._bar_colors(len(values))

        panel.ax.bar(labels, values, color=colors)
        panel.ax.set_ylabel("Selectivity (%)")
        panel.ax.set_xlabel("")
        panel.ax.set_ylim(0, max(5, float(values.max()) * 1.15))
        panel.ax.tick_params(axis="x", rotation=45, labelsize=14)
        self._style_axes(panel.ax)
        panel.draw()

    def draw_carbon_plot(self):
        panel = self.tab_carbon
        panel.clear()

        if not self.data or self.data.c_sel.empty:
            panel.ax.text(
                0.5, 0.5,
                "No carbon-number data to plot.\n(Check 'Formula' looks like C11H22O or include 'Carbon Number'.",
                ha="center", va="center", transform=panel.ax.transAxes
            )
            self._style_axes(panel.ax)
            panel.draw()
            return

        s = self.data.c_sel
        labels = [f"C{int(i)}" for i in s.index]
        values = s.values
        colors = self._bar_colors(len(values))

        panel.ax.bar(labels, values, color=colors)
        panel.ax.set_ylabel("Selectivity (%)")
        panel.ax.set_xlabel("")
        panel.ax.set_ylim(0, max(5, float(values.max()) * 1.15))
        self._style_axes(panel.ax)
        panel.draw()

    def draw_bins_plot(self):
        panel = self.tab_bins
        panel.clear()

        if not self.data or self.data.bin_sel.empty:
            panel.ax.text(0.5, 0.5, "No binned data to plot.", ha="center", va="center", transform=panel.ax.transAxes)
            self._style_axes(panel.ax)
            panel.draw()
            return

        s = self.data.bin_sel
        labels = list(s.index)
        values = s.values
        colors = self._bar_colors(len(values))

        panel.ax.bar(labels, values, color=colors)
        panel.ax.set_ylabel("Selectivity (%)")
        panel.ax.set_xlabel("")
        panel.ax.set_ylim(0, max(5, float(values.max()) * 1.15))
        self._style_axes(panel.ax)
        panel.draw()

    def draw_products_plot(self):
        panel = self.tab_products
        panel.clear()

        if not self.data or self.data.product_sel.empty:
            panel.ax.text(0.5, 0.5, "No selected-product data to plot.", ha="center", va="center", transform=panel.ax.transAxes)
            self._style_axes(panel.ax)
            panel.draw()
            return

        s = self.data.product_sel
        labels = list(s.index)
        values = s.values
        colors = self._bar_colors(len(values))

        panel.ax.bar(labels, values, color=colors)
        panel.ax.set_ylabel("Selectivity (%)")
        panel.ax.set_xlabel("")
        panel.ax.set_ylim(0, max(5, float(values.max()) * 1.15))
        panel.ax.tick_params(axis="x", rotation=75, labelsize=12)
        self._style_axes(panel.ax)
        panel.draw()


def main():
    app = QApplication(sys.argv)
    w = MainWindow()
    w.resize(1350, 820)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
