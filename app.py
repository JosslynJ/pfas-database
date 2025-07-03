import streamlit as st
import pandas as pd
import requests
import pubchempy as pcp
import py3Dmol                            # pip install py3Dmol
from streamlit.components.v1 import html
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# 1. Page config
st.set_page_config(page_title="PFAS Database", layout="wide")

# 2. Load data
df = pd.read_csv("pfas_fina.csv")
df.insert(0, "ID", range(1, len(df) + 1))

# 3. Title
st.title("🔬 PFAS Chemical Database")

# 4. Sidebar filters
st.sidebar.header("🔍 Filter Options")
opts_pfas  = df["Is_PFAS"].dropna().unique()
opts_cls   = df["PFAS_Structure_Class"].dropna().unique()
opts_str   = df["Structure_Class"].dropna().unique()
opts_use   = df["Use_Category"].dropna().unique()

sel_pfas = st.sidebar.multiselect("PFAS Status", opts_pfas)
sel_cls  = st.sidebar.multiselect("PFAS Structure Class", opts_cls)
sel_str  = st.sidebar.multiselect("Structure Class", opts_str)
sel_use  = st.sidebar.multiselect("Use Category", opts_use)

# 5. Apply filters
f = df.copy()
if sel_pfas: f = f[f["Is_PFAS"].isin(sel_pfas)]
if sel_cls:  f = f[f["PFAS_Structure_Class"].isin(sel_cls)]
if sel_str:  f = f[f["Structure_Class"].isin(sel_str)]
if sel_use:  f = f[f["Use_Category"].isin(sel_use)]

# 6. Render table w/ AgGrid
gb = GridOptionsBuilder.from_dataframe(f)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
gridOptions = gb.build()
resp = AgGrid(f, gridOptions=gridOptions,
              update_mode=GridUpdateMode.SELECTION_CHANGED,
              height=500, allow_unsafe_jscode=True)

# 7. Column descriptions
with st.expander("ℹ️ Column Descriptions"):
    st.markdown("""
- **Use Category**: e.g. Pharmaceutical, Pesticide…  
- **Structure Class**: Aromatic, Heterocycle, Chain…  
- **PFAS Structure Class**: CF₃-containing, CF₂-containing…  
- **PFAS Status**: Yes/No per OECD PFAS definition  
    """)

# 8. Helper: get CID from PubChemPy
def get_cid(smiles, name):
    # 优先 SMILES 查询
    try:
        comps = pcp.get_compounds(smiles, namespace="smiles")
        if comps: return comps[0].cid
    except: pass
    # 再按 Name
    try:
        comps = pcp.get_compounds(name, namespace="name")
        if comps: return comps[0].cid
    except: pass
    return None

# 9. Show details + 2D/3D
sel = pd.DataFrame(resp["selected_rows"])
if sel.shape[0]:
    row = sel.iloc[0]
    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1,2])

    # 取 PubChem CID
    cid = None
    if str(row.get("CAS_or_Identifier","")).startswith("CID:"):
        cid = int(row["CAS_or_Identifier"].split("CID:")[1])
    else:
        cid = get_cid(row["SMILES"], row["Name"])

    with col1:
        # 2D image from PubChem
        if cid:
            png_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
            st.image(png_url, caption="2D structure", use_column_width=True)
        else:
            st.warning("No CID → cannot fetch 2D image")

        # 3D model via py3Dmol
        if cid:
            sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
            r = requests.get(sdf_url)
            if r.ok:
                view = py3Dmol.view(width=350, height=250)
                view.addModel(r.text, "sdf")
                view.setStyle({"stick": {}})
                view.zoomTo()
                html(view._make_html(), height=300)
            else:
                st.error("Failed to fetch 3D SDF")
        else:
            st.info("3D preview unavailable")

    with col2:
        st.markdown(f"""
**ID:** {row['ID']}  
**Name:** {row['Name']}  
**CAS / CID:** {row.get('CAS_or_Identifier','')}  
**Exact Mass:** {row.get('Exact_Mass','')}  
**PFAS Status:** {row.get('Is_PFAS','')}  
**PFAS Structure Class:** {row.get('PFAS_Structure_Class','')}  
**Structure Class:** {row.get('Structure_Class','')}  
**Use Category:** {row.get('Use_Category','')}
""")
else:
    st.info("Click a table row to view details")

# 10. Footer
st.markdown("""
<div style="
    position: fixed; bottom:10px; right:10px;
    font-size:14px; color:#888; opacity:0.7;
">
    Created by Josslyn
</div>
""", unsafe_allow_html=True)
