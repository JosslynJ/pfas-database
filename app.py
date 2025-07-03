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

# 2. 检查文件是否存在
CSV_FILE = "pfas_fina.csv"
if not os.path.exists(CSV_FILE):
    st.error(f"❌ 找不到文件: {CSV_FILE}，请确认文件名和路径完全匹配。")
    st.stop()

# 3. 读取 CSV 并打印行列数
df = pd.read_csv(CSV_FILE)
st.write("✅ 已加载 CSV 行数:", df.shape[0], "，列数:", df.shape[1])
st.write("列名:", df.columns.tolist())

# 4. 插入 ID 列
df.insert(0, "ID", range(1, len(df) + 1))

# 5. 标题
st.title("🔬 PFAS Chemical Database")

# 6. 侧边栏筛选
st.sidebar.header("🔍 Filter Options")
opts_pfas = df["Is_PFAS"].dropna().unique().tolist()
opts_psc  = df["PFAS_Structure_Class"].dropna().unique().tolist()
opts_sc   = df["Structure_Class"].dropna().unique().tolist()
opts_use  = df["Use_Category"].dropna().unique().tolist()

sel_pfas = st.sidebar.multiselect("PFAS Status", opts_pfas)
sel_psc  = st.sidebar.multiselect("PFAS Structure Class", opts_psc)
sel_sc   = st.sidebar.multiselect("Structure Class", opts_sc)
sel_use  = st.sidebar.multiselect("Use Category", opts_use)

# 7. 应用筛选
fdf = df.copy()
if sel_pfas: fdf = fdf[fdf["Is_PFAS"].isin(sel_pfas)]
if sel_psc:  fdf = fdf[fdf["PFAS_Structure_Class"].isin(sel_psc)]
if sel_sc:   fdf = fdf[fdf["Structure_Class"].isin(sel_sc)]
if sel_use:  fdf = fdf[fdf["Use_Category"].isin(sel_use)]

# 8. 构建 AgGrid（固定高度 400px，内部滚动）
gb = GridOptionsBuilder.from_dataframe(fdf)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
grid_options = gb.build()

grid_response = AgGrid(
    fdf,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=400,
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# 9. 列说明
with st.expander("ℹ️ Column Descriptions"):
    st.markdown("""
- **Use Category**: 自动类别 (e.g. Pharmaceutical, Pesticide…)  
- **Structure Class**: 结构类型 (Fluorinated aromatic, Chain…)  
- **PFAS Structure Class**: CF₃/CF₂ 分类  
- **PFAS Status**: Yes/No  
    """)

# 10. PubChem CID 获取函数
def get_cid(smiles, name):
    try:
        comps = pcp.get_compounds(smiles, namespace="smiles")
        if comps: return comps[0].cid
    except: pass
    try:
        comps = pcp.get_compounds(name, namespace="name")
        if comps: return comps[0].cid
    except: pass
    return None

# 11. 详情 & 2D/3D 预览
selected = pd.DataFrame(grid_response["selected_rows"])
if not selected.empty:
    row = selected.iloc[0]
    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1, 2])

    # 确定 CID
    raw = str(row.get("CAS_or_Identifier", ""))
    if raw.startswith("CID:"):
        cid = int(raw.split("CID:")[1])
    else:
        cid = get_cid(row["SMILES"], row["Name"])

    with col1:
        # 2D
        if cid:
            png_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
            st.image(png_url, caption="2D Structure", use_column_width=True)
        else:
            st.warning("No CID → cannot fetch 2D image")

        # 3D
        if cid:
            sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
            r = requests.get(sdf_url)
            if r.ok and r.text:
                view = py3Dmol.view(width=350, height=250)
                view.addModel(r.text, "sdf")
                view.setStyle({"stick": {}})
                view.zoomTo()
                html(view._make_html(), height=300)
            else:
                st.error("3D SDF 下载失败")
        else:
            st.info("3D preview unavailable")

    with col2:
        st.markdown(f"""
**ID:** {row['ID']}  
**Name:** {row['Name']}  
**CAS/CID:** {row.get('CAS_or_Identifier','')}  
**Exact Mass:** {row.get('Exact_Mass','')}  
**PFAS Status:** {row.get('Is_PFAS','')}  
**PFAS Structure Class:** {row.get('PFAS_Structure_Class','')}  
**Structure Class:** {row.get('Structure_Class','')}  
**Use Category:** {row.get('Use_Category','')}
""")
else:
    st.info("点击一行查看分子信息。")

# 12. 页脚
st.markdown("""
<div style="
    position: fixed; bottom: 10px; right: 10px;
    font-size: 14px; color: #888; opacity: 0.7;
">
    Created by Josslyn
</div>
""", unsafe_allow_html=True)
