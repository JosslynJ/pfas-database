import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# Set wide layout
st.set_page_config(page_title="PFAS Database", layout="wide")

# Load data
df = pd.read_csv("pfas_data.csv")

st.title("🔬 PFAS Chemical Database")

# Sidebar filters
st.sidebar.header("🔍 Filter Options")
compound_classes = df["Compound_class"].dropna().unique()
potential_uses = df["Potential_use"].dropna().unique()

selected_class = st.sidebar.multiselect("Compound Class", compound_classes)
selected_use = st.sidebar.multiselect("Potential Use", potential_uses)

# Apply filters
filtered_df = df.copy()
if selected_class:
    filtered_df = filtered_df[filtered_df["Compound_class"].isin(selected_class)]
if selected_use:
    filtered_df = filtered_df[filtered_df["Potential_use"].isin(selected_use)]

# Configure AgGrid
gb = GridOptionsBuilder.from_dataframe(filtered_df)
gb.configure_selection("single", use_checkbox=False)
# gb.configure_column("SMILES", hide=True)
grid_options = gb.build()

# Debug: Show table and shape
st.write("Data shape:", filtered_df.shape)
st.dataframe(filtered_df)

# Render table
grid_response = AgGrid(
    filtered_df,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=500,
    fit_columns_on_grid_load=True,
    allow_unsafe_jscode=True
)

# Handle selected row
selected = grid_response["selected_rows"]

if len(selected) > 0:
    # 兼容 dict 或 DataFrame
    row = selected[0]
    st.write("DEBUG-选中的 row:", row)  # 调试代码，部署后可删除

    smiles = row["SMILES"] if "SMILES" in row else None
    st.write("DEBUG-SMILES:", smiles)  # 调试代码，部署后可删除

    st.markdown("### 🧬 Selected Compound Info")
    col1, col2 = st.columns([1, 2])
    with col1:
        try:
            if smiles is None or pd.isna(smiles) or smiles == "":
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
        **PubChem CID:** {row.get('PubChem_CID', '')}  
        **Exact Mass:** {row.get('Exact_Mass', '')}  
        **m/z:** {row.get('mz', '')}  
        **Compound Class:** {row.get('Compound_class', '')}  
        **Potential Use:** {row.get('Potential_use', '')}
        """)
else:
    st.info("Click a row in the table to view compound structure and details.")
