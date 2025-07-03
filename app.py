import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol                         # pip install py3Dmol
from streamlit.components.v1 import html
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# 1. Page configuration
st.set_page_config(page_title="PFAS Database", layout="wide")

# 2. Load data
df = pd.read_csv("pfas_fina.csv")
df.insert(0, "ID", range(1, len(df) + 1))

# 3. Title
st.title("🔬 PFAS Chemical Database")

# 4. Sidebar filters
st.sidebar.header("🔍 Filter Options")
is_pfas_options        = df["Is_PFAS"].dropna().unique()
pfas_structure_classes = df["PFAS_Structure_Class"].dropna().unique()
structure_classes      = df["Structure_Class"].dropna().unique()
use_categories         = df["Use_Category"].dropna().unique()

selected_is_pfas  = st.sidebar.multiselect("PFAS Status (Is_PFAS)", is_pfas_options)
selected_pfas     = st.sidebar.multiselect("PFAS Structure Class", pfas_structure_classes)
selected_struct   = st.sidebar.multiselect("Structure Class",      structure_classes)
selected_use      = st.sidebar.multiselect("Use Category",         use_categories)

# 5. Apply filters
filtered_df = df.copy()
if selected_is_pfas:
    filtered_df = filtered_df[filtered_df["Is_PFAS"].isin(selected_is_pfas)]
if selected_struct:
    filtered_df = filtered_df[filtered_df["Structure_Class"].isin(selected_struct)]
if selected_use:
    filtered_df = filtered_df[filtered_df["Use_Category"].isin(selected_use)]
if selected_pfas:
    filtered_df = filtered_df[filtered_df["PFAS_Structure_Class"].isin(selected_pfas)]

# 6. Configure AgGrid
gb = GridOptionsBuilder.from_dataframe(filtered_df)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
grid_options = gb.build()

# 7. Column descriptions
with st.expander("ℹ️ Column Descriptions and PFAS Classification"):
    st.markdown("""
- **Use Category**: Automatically assigned category based on name keywords (e.g. Pharmaceutical, Pesticide, etc.)  
- **Structure Class**: Broad structural type (e.g. Fluorinated aromatic, Fluoroalkyl chain, Heterocycle)  
- **PFAS Structure Class**: PFAS sub-classification (e.g. CF₃-containing organics, CF₂-containing organics)  
- **PFAS Status (Is_PFAS)**: “Yes” if contains ≥1 per-fluorinated methyl (–CF₃) or methylene (–CF₂–) not bound to H; else “No”

> **Reference:** OECD (2021), *Reconciling Terminology of the Universe of Per- and Polyfluoroalkyl Substances (PFASs)*, OECD Environment, Health and Safety Publications.
    """)

# 8. Render table
grid_response = AgGrid(
    filtered_df,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=500,
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# 3D helper
def show_3d(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=350, height=250)
    view.addModel(block, format='sdf')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view

# 9. Show details & 3D preview
selected = pd.DataFrame(grid_response["selected_rows"])
if not selected.empty:
    row = selected.iloc[0]
    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1, 2])
    with col1:
        # 2D image
        mol2d = Chem.MolFromSmiles(row["SMILES"])
        img = Draw.MolToImage(mol2d, size=(300, 300))
        st.image(img, caption=f"2D Structure of {row['Name']}")

        # 3D model
        view = show_3d(row["SMILES"])
        html(view._make_html(), height=300)
    with col2:
        st.markdown(f"""
**ID:** {row['ID']}  
**Name:** {row['Name']}  
**Identifier:** {row.get('CAS_or_Identifier','')}  
**Exact Mass:** {row.get('Exact_Mass','')}  
**PFAS Status:** {row.get('Is_PFAS','')}  
**PFAS Structure Class:** {row.get('PFAS_Structure_Class','')}  
**Structure Class:** {row.get('Structure_Class','')}  
**Use Category:** {row.get('Use_Category','')}
""")
else:
    st.info("Click a row in the table to view compound structure and details.")

# 10. Creator label
st.markdown(
    """
    <div style="
        position: fixed;
        bottom: 10px;
        right: 10px;
        font-size: 18px;
        color: #888;
        opacity: 0.7;
    ">
        Created by Josslyn
    </div>
    """,
    unsafe_allow_html=True
)
