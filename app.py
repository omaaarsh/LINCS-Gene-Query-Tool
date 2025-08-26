import streamlit as st
import requests
import pandas as pd
import io
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

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

def create_volcano_plot(up_df, down_df, gene_name):
    """Create a volcano plot showing fold change vs significance"""
    # Combine data
    up_df_copy = up_df.copy()
    down_df_copy = down_df.copy()
    up_df_copy['Regulation'] = 'Up-regulated'
    down_df_copy['Regulation'] = 'Down-regulated'
    
    combined_df = pd.concat([up_df_copy, down_df_copy])
    
    # Create volcano plot
    fig = px.scatter(
        combined_df, 
        x='Log2(Fold Change)', 
        y='CD Coefficient',
        color='Regulation',
        hover_data=['Perturbagen', 'Dose', 'Cell Line'],
        title=f'Differential Expression of {gene_name} in RNA-seq-like Chemical Signatures',
        labels={'Log2(Fold Change)': 'Log2(Fold Change)', 'CD Coefficient': 'CD Coefficient (Significance)'},
        color_discrete_map={'Up-regulated': '#ff4444', 'Down-regulated': '#4444ff'}
    )
    
    fig.update_layout(
        width=800, 
        height=600,
        template='plotly_white'
    )
    
    return fig

def create_dose_response_plot(df, regulation_type, gene_name):
    """Create dose-response analysis plot"""
    # Extract numeric dose values where possible
    df_copy = df.copy()
    df_copy['Dose_Numeric'] = df_copy['Dose'].str.extract('(\d+\.?\d*)').astype(float)
    
    # Filter out rows where dose couldn't be parsed
    dose_df = df_copy.dropna(subset=['Dose_Numeric'])
    
    if len(dose_df) > 10:  # Only show if we have enough data points
        fig = px.scatter(
            dose_df.head(50),  # Top 50 for readability
            x='Dose_Numeric',
            y='CD Coefficient',
            size='Fold Change',
            hover_data=['Perturbagen', 'Cell Line'],
            title=f'Dose-Response Analysis: {regulation_type} of {gene_name}',
            labels={'Dose_Numeric': 'Dose (μM)', 'CD Coefficient': 'CD Coefficient'},
            log_x=True
        )
        fig.update_layout(template='plotly_white')
def convert_df_to_csv(df):
    """Convert DataFrame to CSV for download"""
    return df.to_csv(index=False).encode('utf-8')
    return None

def create_cell_line_analysis(df, gene_name, regulation_type):
    """Analyze effects across different cell lines"""
    cell_line_stats = df.groupby('Cell Line').agg({
        'CD Coefficient': ['mean', 'count'],
        'Fold Change': 'mean'
    }).round(4)
    
    cell_line_stats.columns = ['Mean_CD_Coeff', 'Count', 'Mean_Fold_Change']
    cell_line_stats = cell_line_stats.reset_index()
    cell_line_stats = cell_line_stats[cell_line_stats['Count'] >= 3].sort_values('Mean_CD_Coeff', ascending=False)
    
    if len(cell_line_stats) > 0:
        fig = px.bar(
            cell_line_stats.head(15),
            x='Cell Line',
            y='Mean_CD_Coeff',
            title=f'Cell Line Sensitivity: {regulation_type} of {gene_name}',
            labels={'Mean_CD_Coeff': 'Mean CD Coefficient', 'Cell Line': 'Cell Line'},
            hover_data=['Count', 'Mean_Fold_Change']
        )
        fig.update_layout(
            xaxis_tickangle=-45,
            template='plotly_white',
            height=500
        )
        return fig, cell_line_stats
    return None, None

def create_timepoint_analysis(df, gene_name, regulation_type):
    """Analyze temporal effects"""
    timepoint_stats = df.groupby('Timepoint').agg({
        'CD Coefficient': ['mean', 'count'],
        'Fold Change': 'mean'
    }).round(4)
    
    timepoint_stats.columns = ['Mean_CD_Coeff', 'Count', 'Mean_Fold_Change']
    timepoint_stats = timepoint_stats.reset_index()
    timepoint_stats = timepoint_stats[timepoint_stats['Count'] >= 5].sort_values('Mean_CD_Coeff', ascending=False)
    
    if len(timepoint_stats) > 0:
        fig = px.bar(
            timepoint_stats,
            x='Timepoint',
            y='Mean_CD_Coeff',
            title=f'Temporal Response: {regulation_type} of {gene_name}',
            labels={'Mean_CD_Coeff': 'Mean CD Coefficient', 'Timepoint': 'Treatment Duration'},
            hover_data=['Count', 'Mean_Fold_Change']
        )
        fig.update_layout(template='plotly_white')
        return fig, timepoint_stats
    return None, None

def create_top_compounds_plot(df, gene_name, regulation_type, top_n=20):
    """Show top compounds by effectiveness"""
    top_compounds = df.head(top_n)
    
    fig = px.bar(
        top_compounds,
        x='CD Coefficient',
        y='Perturbagen',
        orientation='h',
        title=f'Top {top_n} Most Effective Compounds: {regulation_type} of {gene_name}',
        labels={'CD Coefficient': 'CD Coefficient', 'Perturbagen': 'Compound'},
        hover_data=['Dose', 'Cell Line', 'Fold Change']
    )
    fig.update_layout(
        height=max(400, top_n * 25),
        template='plotly_white'
    )
    return fig

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
    print(gene_upper)
    
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
            
            # Analysis Section
            st.markdown("---")
            st.header(f"📊 Analysis: Differential Expression of {gene_upper} in RNA-seq-like Chemical Signatures")
            
            # Create analysis tabs
            if len(results) > 1:
                analysis_tabs = st.tabs(["🌋 Volcano Plot", "📈 Up-regulation Analysis", "📉 Down-regulation Analysis", "📋 Summary"])
                
                with analysis_tabs[0]:
                    st.subheader(f"Volcano Plot: {gene_upper} Expression Changes")
                    volcano_fig = create_volcano_plot(results['up'], results['down'], gene_upper)
                    st.plotly_chart(volcano_fig, use_container_width=True)
                    
                    st.markdown("""
                    **Interpretation:**
                    - **X-axis**: Log2(Fold Change) - magnitude of gene expression change
                    - **Y-axis**: CD Coefficient - significance/confidence of the effect
                    - **Red points**: Up-regulating compounds
                    - **Blue points**: Down-regulating compounds
                    """)
                
                with analysis_tabs[1]:
                    if 'up' in results:
                        st.subheader(f"📈 Up-regulation Analysis for {gene_upper}")
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            top_up_fig = create_top_compounds_plot(results['up'], gene_upper, "Up-regulation")
                            st.plotly_chart(top_up_fig, use_container_width=True)
                        
                        with col2:
                            dose_up_fig = create_dose_response_plot(results['up'], "Up-regulation", gene_upper)
                            if dose_up_fig:
                                st.plotly_chart(dose_up_fig, use_container_width=True)
                        
                        # Cell line analysis
                        cell_fig_up, cell_stats_up = create_cell_line_analysis(results['up'], gene_upper, "Up-regulation")
                        if cell_fig_up:
                            st.plotly_chart(cell_fig_up, use_container_width=True)
                
                with analysis_tabs[2]:
                    if 'down' in results:
                        st.subheader(f"📉 Down-regulation Analysis for {gene_upper}")
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            top_down_fig = create_top_compounds_plot(results['down'], gene_upper, "Down-regulation")
                            st.plotly_chart(top_down_fig, use_container_width=True)
                        
                        with col2:
                            dose_down_fig = create_dose_response_plot(results['down'], "Down-regulation", gene_upper)
                            if dose_down_fig:
                                st.plotly_chart(dose_down_fig, use_container_width=True)
                        
                        # Cell line analysis
                        cell_fig_down, cell_stats_down = create_cell_line_analysis(results['down'], gene_upper, "Down-regulation")
                        if cell_fig_down:
                            st.plotly_chart(cell_fig_down, use_container_width=True)
                
                with analysis_tabs[3]:
                    st.subheader(f"📋 Summary Statistics for {gene_upper}")
                    
                    col1, col2 = st.columns(2)
                    
                    if 'up' in results:
                        with col1:
                            st.markdown("### 📈 Up-regulation Summary")
                            st.metric("Total Compounds", len(results['up']))
                            st.metric("Max CD Coefficient", f"{results['up']['CD Coefficient'].max():.4f}")
                            st.metric("Mean Fold Change", f"{results['up']['Fold Change'].mean():.3f}")
                            st.metric("Unique Cell Lines", results['up']['Cell Line'].nunique())
                            st.metric("Unique Timepoints", results['up']['Timepoint'].nunique())
                    
                    if 'down' in results:
                        with col2:
                            st.markdown("### 📉 Down-regulation Summary")
                            st.metric("Total Compounds", len(results['down']))
                            st.metric("Min CD Coefficient", f"{results['down']['CD Coefficient'].min():.4f}")
                            st.metric("Mean Fold Change", f"{results['down']['Fold Change'].mean():.3f}")
                            st.metric("Unique Cell Lines", results['down']['Cell Line'].nunique())
                            st.metric("Unique Timepoints", results['down']['Timepoint'].nunique())
            
            # Data Tables Section
            st.markdown("---")
            st.header("📊 Raw Data Tables")
            
            # Create tabs for up and down regulation data tables
            if len(results) > 1:
                data_tab_up, data_tab_down = st.tabs(["📈 Up-regulating Data", "📉 Down-regulating Data"])
            else:
                if 'up' in results:
                    data_tab_up = st.container()
                else:
                    data_tab_down = st.container()
            
            # Up-regulating results
            if get_up and 'up' in results:
                with (data_tab_up if len(results) > 1 else st.container()):
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
            
            # Down-regulating results
            if get_down and 'down' in results:
                with (data_tab_down if len(results) > 1 else st.container()):
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
        5. **Explore the analysis** with interactive plots and statistics
        6. **Download CSV files** with the gene name included in filename
        
        **Analysis Features:**
        - **Volcano Plot**: Overview of all compounds and their effects
        - **Top Compounds**: Most effective perturbagens ranked by CD coefficient
        - **Dose-Response**: Relationship between compound dose and effect magnitude
        - **Cell Line Sensitivity**: Which cell lines show strongest responses
        - **Summary Statistics**: Key metrics and data overview
        """)

# Update requirements.txt note
st.sidebar.markdown("---")
st.sidebar.markdown("**💡 New in this version:**")
st.sidebar.markdown("• Interactive analysis plots")
st.sidebar.markdown("• Volcano plot visualization") 
st.sidebar.markdown("• Dose-response analysis")
st.sidebar.markdown("• Cell line comparisons")
