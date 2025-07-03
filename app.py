import os
import streamlit as st
import pandas as pd
import requests
import pubchempy as pcp
import py3Dmol                        # pip install py3Dmol
from streamlit.components.v1 import html
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# 1. 页面配置
st.set_page_config(page_title="PFAS Database", layout="wide")

# 2. 读取 CSV
CSV_FILE = "pfas_fina.csv"
if not os.path.exists(CSV_FILE):
    st.error(f"❌ 找不到文件: {CSV_FILE}")
    st.stop()

df = pd.read_csv(CSV_FILE)
df.insert(0, "ID", range(1, len(df) + 1))

# 3. 标题
st.title("🔬 PFAS Chemical Database")

# 4. Sidebar：筛选 + Debug Info
with st.sidebar:
    st.header("🔍 Filter Options")
    opts_pfas = df["Is_PFAS"].dropna().unique().tolist()
    opts_psc  = df["PFAS_Structure_Class"].dropna().unique().tolist()
    opts_sc   = df["Structure_Class"].dropna().unique().tolist()
    opts_use  = df["Use_Category"].dropna().unique().tolist()

    sel_pfas = st.multiselect("PFAS Status", opts_pfas)
    sel_psc  = st.multiselect("PFAS Structure Class", opts_psc)
    sel_sc   = st.multiselect("Structure Class", opts_sc)
    sel_use  = st.multiselect("Use Category", opts_use)

    with st.expander("🔧 Debug Info", expanded=False):
        st.markdown(f"- loaded rows: **{df.shape[0]}**\n- cols: **{df.shape[1]}**")
        st.markdown(f"**Columns:** {', '.join(df.columns)}")

# 5. 应用筛选
fdf = df.copy()
if sel_pfas: fdf = fdf[fdf["Is_PFAS"].isin(sel_pfas)]
if sel_psc:  fdf = fdf[fdf["PFAS_Structure_Class"].isin(sel_psc)]
if sel_sc:   fdf = fdf[fdf["Structure_Class"].isin(sel_sc)]
if sel_use:  fdf = fdf[fdf["Use_Category"].isin(sel_use)]

# 6. 渲染主表格
gb = GridOptionsBuilder.from_dataframe(fdf)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
grid = gb.build()

grid_response = AgGrid(
    fdf,
    gridOptions=grid,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=500,
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# 7. 选中行详情 & 2D/3D 预览
def get_cid(smiles, name):
    try:
        c = pcp.get_compounds(smiles, namespace="smiles")
        if c: return c[0].cid
    except: pass
    try:
        c = pcp.get_compounds(name, namespace="name")
        if c: return c[0].cid
    except: pass
    return None

selected = pd.DataFrame(grid_response["selected_rows"])
if not selected.empty:
    row = selected.iloc[0]
    st.markdown("### 🧬 Selected Compound Info")
    c1, c2 = st.columns([1, 2])

    # —— 2D 图像, 缩小到 200×200 ——  
    cid = get_cid(row["SMILES"], row["Name"])
    if cid:
        png = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
        c1.image(png, caption="2D Structure", width=200)
    else:
        c1.warning("2D 图不可用")

    # —— 3D 预览, 缩小到 250×200 ——  
    if cid:
        sdf = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        r = requests.get(sdf)
        if r.ok:
            view = py3Dmol.view(width=250, height=200)
            view.addModel(r.text, "sdf")
            view.setStyle({"stick": {}})
            view.zoomTo()
            c1.components(html(view._make_html(), height=200))
        else:
            c1.error("3D 下载失败")
    else:
        c1.info("3D 不可用")

    # —— 文本信息 ——  
    c2.markdown(f"""
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
    st.info("👈 Click a row to view details.")

# 8. 主视图底部：列说明
with st.expander("ℹ️ Column Descriptions"):
    st.markdown("""
- **Use Category**：按名称关键词自动分类  
- **Structure Class**：芳香 / 杂环 / 链烷  
- **PFAS Structure Class**：CF₃ / CF₂ / Sulfonate 等子类  
- **PFAS Status**：是否含 ≥1 个 –CF₃ 或 –CF₂–  
    """)

# 9. Creator label
st.markdown("""
<div style="
  position: fixed; bottom: 10px; right: 10px;
  font-size: 14px; color: #888; opacity: 0.7;
">
  Created by Josslyn
</div>
""", unsafe_allow_html=True)
