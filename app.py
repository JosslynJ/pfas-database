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

# 4. 侧栏筛选
st.sidebar.header("🔍 Filter Options")
opts_pfas = df["Is_PFAS"].dropna().unique().tolist()
opts_psc  = df["PFAS_Structure_Class"].dropna().unique().tolist()
opts_sc   = df["Structure_Class"].dropna().unique().tolist()
opts_use  = df["Use_Category"].dropna().unique().tolist()

sel_pfas = st.sidebar.multiselect("PFAS Status", opts_pfas)
sel_psc  = st.sidebar.multiselect("PFAS Structure Class", opts_psc)
sel_sc   = st.sidebar.multiselect("Structure Class", opts_sc)
sel_use  = st.sidebar.multiselect("Use Category", opts_use)

# 5. 应用筛选
fdf = df.copy()
if sel_pfas: fdf = fdf[fdf["Is_PFAS"].isin(sel_pfas)]
if sel_psc:  fdf = fdf[fdf["PFAS_Structure_Class"].isin(sel_psc)]
if sel_sc:   fdf = fdf[fdf["Structure_Class"].isin(sel_sc)]
if sel_use:  fdf = fdf[fdf["Use_Category"].isin(sel_use)]

# 6. 渲染表格
gb = GridOptionsBuilder.from_dataframe(fdf)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
grid_options = gb.build()

grid_response = AgGrid(
    fdf,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=500,               # 调整可视高度
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# 7. 选中行的详情 & 2D+3D 预览
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
    col1, col2 = st.columns([1, 2])

    # 2D image
    cid = get_cid(row["SMILES"], row["Name"])
    if cid:
        png_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
        st.image(png_url, caption="2D Structure", use_column_width=True)
    else:
        st.warning("无法获取 PubChem CID，2D 图像不可用。")

    # 3D 模型
    if cid:
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        resp = requests.get(sdf_url)
        if resp.ok:
            view = py3Dmol.view(width=350, height=250)
            view.addModel(resp.text, "sdf")
            view.setStyle({"stick": {}})
            view.zoomTo()
            html(view._make_html(), height=300)
        else:
            st.error("3D SDF 下载失败。")
    else:
        st.info("3D 预览不可用。")

    # 文字详情
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
    st.info("👈 Click a row to view details.")

# 8. 在表格下面放 Debug 和 列说明
with st.expander("🔧 Debug Info"):
    st.write("✅ loaded rows:", df.shape[0], "cols:", df.shape[1])
    st.write("cols:", df.columns.tolist())

with st.expander("ℹ️ Column Descriptions"):
    st.markdown("""
- **Use Category**：根据名称关键词自动分配  
- **Structure Class**：芳香、杂环还是链状  
- **PFAS Structure Class**：CF₃ / CF₂ / sulfonate 等子分类  
- **PFAS Status**：是否符合 ≥1 个 –CF₃ 或 –CF₂– 的 PFAS 定义  
    """)

# 9. 页脚
st.markdown("""
<div style="
  position: fixed; bottom: 10px; right: 10px;
  font-size: 14px; color: #888; opacity: 0.7;
">
  Created by Josslyn
</div>
""", unsafe_allow_html=True)
