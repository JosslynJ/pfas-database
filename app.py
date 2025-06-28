import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import math

st.set_page_config(page_title="PFAS Database", layout="wide")
st.title("🔬 PFAS Chemical Database")

# 载入数据
df = pd.read_csv("pfas_data.csv")

# 过滤非法 SMILES 行，避免 RDKit 报错
def is_valid_smiles(smiles):
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except:
        return False

df = df[df["SMILES"].apply(is_valid_smiles)]

# 侧边栏筛选
st.sidebar.header("🔍 Filter Options")
compound_classes = df["Compound_class"].dropna().unique()
potential_uses = df["Potential_use"].dropna().unique()

selected_class = st.sidebar.multiselect("Compound Class", compound_classes)
selected_use = st.sidebar.multiselect("Potential Use", potential_uses)

# 应用筛选条件
filtered_df = df.copy()
if selected_class:
    filtered_df = filtered_df[filtered_df["Compound_class"].isin(selected_class)]
if selected_use:
    filtered_df = filtered_df[filtered_df["Potential_use"].isin(selected_use)]

# 分页功能
PAGE_SIZE = 50
total_pages = max(1, math.ceil(len(filtered_df) / PAGE_SIZE))
page_num = st.sidebar.number_input("Page", min_value=1, max_value=total_pages, value=1, step=1)

start_idx = (page_num - 1) * PAGE_SIZE
end_idx = start_idx + PAGE_SIZE
page_df = filtered_df.iloc[start_idx:end_idx]

# 配置 AgGrid
gb = GridOptionsBuilder.from_dataframe(page_df)
gb.configure_selection("single", use_checkbox=False)
gb.configure_column("SMILES", hide=True)
grid_options = gb.build()

# 显示分页表格
grid_response = AgGrid(
    page_df,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=500,
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# 显示详细信息
selected = pd.DataFrame(grid_response["selected_rows"])
if len(selected) > 0:
    row = selected.iloc[0]
    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1, 2])
    with col1:
        mol = Chem.MolFromSmiles(row["SMILES"])
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, caption=f"Structure of {row['Name']}")
    with col2:
        st.markdown(f"""
        **Name:** {row['Name']}  
        **PubChem CID:** {row['PubChem_CID']}  
        **Exact Mass:** {row['Exact_Mass']}  
        **m/z:** {row['mz']}  
        **Compound Class:** {row['Compound_class']}  
        **Potential Use:** {row['Potential_use']}
        """)
else:
    st.info("Click a row in the table to view compound structure and details.")

# Debug 用，开发时可显示所有数据
# st.dataframe(filtered_df)
