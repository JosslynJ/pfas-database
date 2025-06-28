import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# 加载数据
df = pd.read_csv("pfas_data.csv")

st.set_page_config(page_title="PFAS Chemical Database", layout="wide")
st.title("🔬 PFAS Chemical Database")

# 显示原始数据
st.write("数据总览（调试用，可隐藏）:")
st.dataframe(df)

# 侧边栏筛选
st.sidebar.header("🔍 Filter Options")
compound_classes = df["Compound_class"].dropna().unique()
potential_uses = df["Potential_use"].dropna().unique()

selected_class = st.sidebar.multiselect("Compound Class", compound_classes)
selected_use = st.sidebar.multiselect("Potential Use", potential_uses)

# 筛选数据
filtered_df = df.copy()
if selected_class:
    filtered_df = filtered_df[filtered_df["Compound_class"].isin(selected_class)]
if selected_use:
    filtered_df = filtered_df[filtered_df["Potential_use"].isin(selected_use)]

option = st.selectbox(
    '👇 从下拉菜单中选择一个化合物（支持筛选后列表）',
    filtered_df['Name'].values if not filtered_df.empty else ['无可选项']
)

if filtered_df.empty or option == '无可选项':
    st.warning("没有符合条件的化合物。")
else:
    # 展示详情
    row = filtered_df[filtered_df['Name'] == option].iloc[0]
    # ... 下面不变 ...

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

