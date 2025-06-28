import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

df = pd.read_csv("pfas_data.csv")

st.set_page_config(page_title="PFAS Chemical Database", layout="wide")
st.title("🔬 PFAS Chemical Database")

# 侧边栏筛选
st.sidebar.header("🔍 Filter Options")
compound_classes = df["Compound_class"].dropna().unique()
potential_uses = df["Potential_use"].dropna().unique()

selected_class = st.sidebar.multiselect("Compound Class", compound_classes)
selected_use = st.sidebar.multiselect("Potential Use", potential_uses)

filtered_df = df.copy()
if selected_class:
    filtered_df = filtered_df[filtered_df["Compound_class"].isin(selected_class)]
if selected_use:
    filtered_df = filtered_df[filtered_df["Potential_use"].isin(selected_use)]

# ------- 这块用 AgGrid 展示表格 -------
if filtered_df.empty:
    st.warning("没有符合条件的数据。请调整筛选条件。")
else:
    gb = GridOptionsBuilder.from_dataframe(filtered_df)
    gb.configure_selection('single', use_checkbox=False)
    grid_options = gb.build()

    grid_response = AgGrid(
        filtered_df,
        gridOptions=grid_options,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        height=400,
        fit_columns_on_grid_load=True,
        allow_unsafe_jscode=True,
        theme="streamlit"
    )

    selected = grid_response["selected_rows"]

    if selected:
        row = pd.DataFrame(selected).iloc[0]
        st.markdown("### 🧬 Selected Compound Info")
        col1, col2 = st.columns([1, 2])
        with col1:
            smiles = str(row["SMILES"])
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    st.error("SMILES 无法识别：" + smiles)
                else:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, caption=f"Structure of {row['Name']}")
            except Exception as e:
                st.error(f"RDKit 错误: {e}")
        with col2:
            st.markdown(f"""
            **Name:** {row['Name']}  
            **PubChem CID:** {row['PubChem_CID']}  
            **Exact Mass:** {row['Exact_Mass']}  
            **m/z:** {row['mz'] if 'mz' in row else ''}  
            **Compound Class:** {row['Compound_class']}  
            **Potential Use:** {row['Potential_use']}
            """)
    else:
        st.info("点击上表任意一行，即可查看结构式和详细信息！")
