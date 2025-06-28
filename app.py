import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

st.set_page_config(page_title="PFAS Database", layout="wide")
df = pd.read_csv("pfas_data.csv")
st.title("🔬 PFAS Chemical Database")

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

# 关键：配置 AgGrid 可以整行选中
gb = GridOptionsBuilder.from_dataframe(filtered_df)
gb.configure_selection(selection_mode="single", use_checkbox=True)  # 用checkbox，用户体验好！
grid_options = gb.build()

st.write("Data shape:", filtered_df.shape)

grid_response = AgGrid(
    filtered_df,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=400,
    allow_unsafe_jscode=True,
    theme="streamlit"
)

selected = grid_response.get("selected_rows", [])

st.write("DEBUG-selected:", selected)

if selected and isinstance(selected, list) and len(selected) > 0:
    row = selected[0]
    st.write("DEBUG-row:", row)
    smiles = row.get("SMILES", None)
    st.write("DEBUG-SMILES:", smiles)

    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1, 2])
    with col1:
        try:
            if not smiles or pd.isna(smiles) or smiles == "":
                st.warning("SMILES 字段为空！")
            else:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is None:
                    st.error("SMILES 无法识别：" + str(smiles))
                else:
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, caption=f"Structure of {row.get('Name', '')}")
        except Exception as e:
            st.error(f"RDKit 错误: {e}")
    with col2:
        st.markdown(f"""
        **Name:** {row.get('Name', '')}  
        **PubChem CID:** {row.g
