import streamlit as st
import requests
import pandas as pd
import io

BASE = "https://lincs-reverse-search-dashboard.dev.maayanlab.cloud/api/table"

def lincs_cp(gene: str, direction: str = "up", top_n: int = 10000):
    """Query Chemical Perturbagens that up/down regulate a gene."""
    url = f"{BASE}/cp/{direction}/{gene}"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    df = pd.DataFrame(data)
    # Sort by CD Coefficient (descending = stronger effect)
    df_sorted = df.sort_values(by="CD Coefficient", ascending=False)
    return df_sorted.head(top_n)

def convert_df_to_csv(df):
    """Convert DataFrame to CSV for download"""
    return df.to_csv(index=False).encode('utf-8')

# Streamlit UI
st.set_page_config(page_title="LINCS Gene Query Tool", page_icon="🧬", layout="wide")

st.title("🧬 LINCS Gene Query Tool")
st.markdown("Query chemical perturbagens that regulate your gene of interest")

# Sidebar for inputs
with st.sidebar:
    st.header("Query Parameters")
    gene_name = st.text_input("Enter Gene Name", value="BRAF", help="Enter a valid gene symbol (e.g., BRAF, TP53)")
    
    col1, col2 = st.columns(2)
    with col1:
        get_up = st.checkbox("Up-regulating", value=True)
    with col2:
        get_down = st.checkbox("Down-regulating", value=True)
    
    top_n = st.number_input("Max Results", min_value=10, max_value=10000, value=1000, step=50)
    
    query_button = st.button("🔍 Query Gene", type="primary")

# Main content area
if query_button and gene_name:
    gene_upper = gene_name.upper()
    
    with st.spinner(f'Querying data for {gene_upper}...'):
        try:
            # Initialize containers for results
            results = {}
            
            # Query up-regulating compounds
            if get_up:
                with st.status("Fetching up-regulating compounds..."):
                    up_df = lincs_cp(gene_name, "up", top_n)
                    results['up'] = up_df
                    st.write(f"✅ Found {len(up_df)} up-regulating compounds")
            
            # Query down-regulating compounds  
            if get_down:
                with st.status("Fetching down-regulating compounds..."):
                    down_df = lincs_cp(gene_name, "down", top_n)
                    results['down'] = down_df
                    st.write(f"✅ Found {len(down_df)} down-regulating compounds")
            
            # Display results
            st.success(f"Query completed for gene: **{gene_upper}**")
            
            # Create tabs for up and down regulation
            if len(results) > 1:
                tab_up, tab_down = st.tabs(["📈 Up-regulating", "📉 Down-regulating"])
            else:
                if 'up' in results:
                    tab_up = st.container()
                else:
                    tab_down = st.container()
            
            # Up-regulating results
            if get_up and 'up' in results:
                with (tab_up if len(results) > 1 else st.container()):
                    st.subheader(f"📈 Compounds that UP-regulate {gene_upper}")
                    st.dataframe(results['up'], use_container_width=True, height=400)
                    
                    # Download button
                    csv_up = convert_df_to_csv(results['up'])
                    st.download_button(
                        label="💾 Download Up-regulating Data",
                        data=csv_up,
                        file_name=f"{gene_upper}_Up-regulating_perturbations.csv",
                        mime="text/csv",
                        key="download_up"
                    )
                    
                    # Summary stats
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Compounds", len(results['up']))
                    with col2:
                        st.metric("Max CD Coefficient", f"{results['up']['CD Coefficient'].max():.4f}")
                    with col3:
                        st.metric("Avg Fold Change", f"{results['up']['Fold Change'].mean():.3f}")
            
            # Down-regulating results
            if get_down and 'down' in results:
                with (tab_down if len(results) > 1 else st.container()):
                    st.subheader(f"📉 Compounds that DOWN-regulate {gene_upper}")
                    st.dataframe(results['down'], use_container_width=True, height=400)
                    
                    # Download button
                    csv_down = convert_df_to_csv(results['down'])
                    st.download_button(
                        label="💾 Download Down-regulating Data",
                        data=csv_down,
                        file_name=f"{gene_upper}_Down-regulating_perturbations.csv",
                        mime="text/csv",
                        key="download_down"
                    )
                    
                    # Summary stats
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Compounds", len(results['down']))
                    with col2:
                        st.metric("Min CD Coefficient", f"{results['down']['CD Coefficient'].min():.4f}")
                    with col3:
                        st.metric("Avg Fold Change", f"{results['down']['Fold Change'].mean():.3f}")
        
        except requests.exceptions.RequestException as e:
            st.error(f"❌ API request failed: {str(e)}")
        except Exception as e:
            st.error(f"❌ An error occurred: {str(e)}")

elif query_button and not gene_name:
    st.warning("⚠️ Please enter a gene name to query")

# Footer
with st.sidebar:
    st.markdown("---")
    st.markdown("**About this tool:**")
    st.markdown("Query the LINCS database for chemical compounds that regulate your gene of interest.")
    st.markdown("Data source: [LINCS Reverse Search](https://lincs-reverse-search-dashboard.dev.maayanlab.cloud/)")

# Instructions
if not query_button:
    st.info("👆 Enter a gene name in the sidebar and click 'Query Gene' to get started!")
    
    with st.expander("ℹ️ How to use this tool"):
        st.markdown("""
        1. **Enter a gene name** (e.g., BRAF, TP53, MYC)
        2. **Select regulation type** (up-regulating, down-regulating, or both)
        3. **Set max results** (how many compounds to return)
        4. **Click 'Query Gene'** to fetch data
        5. **Download CSV files** with the gene name included in filename
        
        The tool will show you chemical compounds that affect your gene's expression,
        along with their effectiveness scores and experimental details.
        """)
