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
# gb.configure_column("SMILES", hide=True)  # 可以选择隐藏SMILES列
grid_options = gb.build()

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
selected = pd.DataFrame(grid_response["selected_rows"])

if len(selected) > 0:
    row = selected.iloc[0]
    st.write("你选中的row:", row)  # 新加的调试代码
    st.markdown("### 🧬 Selected Compound Info")

    col1, col2 = st.columns([1, 2])

    with col1:
        try:
            mol = Chem.MolFromSmiles(row["SMILES"])
            if mol is None:
                st.error("SMILES 无法识别：" + str(row["SMILES"]))
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
        **m/z:** {row['mz']}  
        **Compound Class:** {row['Compound_class']}  
        **Potential Use:** {row['Potential_use']}
        """)
else:
    st.info("Click a row in the table to view compound structure and details.")
