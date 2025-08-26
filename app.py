import streamlit as st
import requests
import pandas as pd
import io
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import time
from typing import Optional, Tuple, Dict, Any
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BASE = "https://lincs-reverse-search-dashboard.dev.maayanlab.cloud/api/table"

def validate_gene_name(gene: str) -> bool:
    """Validate gene name format"""
    if not gene or not isinstance(gene, str):
        return False
    # Basic validation: alphanumeric, allow hyphens and underscores
    gene = gene.strip()
    return len(gene) > 0 and gene.replace('-', '').replace('_', '').isalnum()

def lincs_cp(gene: str, direction: str = "up", top_n: int = 10000) -> Optional[pd.DataFrame]:
    """Query Chemical Perturbagens that up/down regulate a gene with error handling."""
    try:
        # Validate inputs
        if not validate_gene_name(gene):
            raise ValueError(f"Invalid gene name format: {gene}")
        
        if direction not in ["up", "down"]:
            raise ValueError(f"Invalid direction: {direction}. Must be 'up' or 'down'")
        
        if not isinstance(top_n, int) or top_n <= 0:
            raise ValueError(f"Invalid top_n value: {top_n}. Must be positive integer")
        
        url = f"{BASE}/cp/{direction}/{gene.strip()}"
        logger.info(f"Querying URL: {url}")
        
        # Add timeout and retry logic
        max_retries = 3
        retry_delay = 2
        
        for attempt in range(max_retries):
            try:
                r = requests.get(url, timeout=30)
                r.raise_for_status()
                break
            except requests.exceptions.Timeout:
                if attempt < max_retries - 1:
                    st.warning(f"Request timeout, retrying in {retry_delay} seconds... (Attempt {attempt + 1}/{max_retries})")
                    time.sleep(retry_delay)
                    retry_delay *= 2  # Exponential backoff
                else:
                    raise requests.exceptions.Timeout("Request timed out after multiple attempts")
            except requests.exceptions.ConnectionError:
                if attempt < max_retries - 1:
                    st.warning(f"Connection error, retrying in {retry_delay} seconds... (Attempt {attempt + 1}/{max_retries})")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                else:
                    raise requests.exceptions.ConnectionError("Failed to connect after multiple attempts")
        
        data = r.json()
        
        # Validate response data
        if not data:
            st.warning(f"No data found for gene {gene} in {direction} direction")
            return pd.DataFrame()
        
        if not isinstance(data, list):
            raise ValueError(f"Unexpected data format: expected list, got {type(data)}")
        
        df = pd.DataFrame(data)
        
        # Validate DataFrame structure
        required_columns = ['CD Coefficient']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Handle missing or invalid numeric data
        numeric_columns = ['CD Coefficient', 'Fold Change', 'Log2(Fold Change)']
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Sort by CD Coefficient (descending = stronger effect)
        df_sorted = df.sort_values(by="CD Coefficient", ascending=False)
        result = df_sorted.head(top_n)
        
        logger.info(f"Successfully retrieved {len(result)} records for {gene} ({direction})")
        return result
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            st.error(f"❌ Gene '{gene}' not found in database. Please check the gene symbol.")
        elif e.response.status_code == 500:
            st.error("❌ Server error. The LINCS database may be temporarily unavailable.")
        else:
            st.error(f"❌ HTTP Error {e.response.status_code}: {str(e)}")
        return None
    except requests.exceptions.RequestException as e:
        st.error(f"❌ Network error: {str(e)}")
        return None
    except ValueError as e:
        st.error(f"❌ Data validation error: {str(e)}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error in lincs_cp: {str(e)}", exc_info=True)
        st.error(f"❌ Unexpected error: {str(e)}")
        return None

def safe_plot_creation(plot_func, *args, **kwargs):
    """Wrapper for safe plot creation with error handling"""
    try:
        return plot_func(*args, **kwargs)
    except Exception as e:
        logger.error(f"Error creating plot {plot_func.__name__}: {str(e)}", exc_info=True)
        st.error(f"❌ Error creating {plot_func.__name__}: {str(e)}")
        return None

def create_volcano_plot(up_df: pd.DataFrame, down_df: pd.DataFrame, gene_name: str) -> Optional[go.Figure]:
    """Create a volcano plot showing fold change vs significance with error handling"""
    try:
        if up_df.empty and down_df.empty:
            st.warning("No data available for volcano plot")
            return None
        
        # Combine data with safety checks
        combined_data = []
        
        if not up_df.empty:
            up_df_copy = up_df.copy()
            up_df_copy['Regulation'] = 'Up-regulated'
            combined_data.append(up_df_copy)
        
        if not down_df.empty:
            down_df_copy = down_df.copy()
            down_df_copy['Regulation'] = 'Down-regulated'
            combined_data.append(down_df_copy)
        
        if not combined_data:
            return None
            
        combined_df = pd.concat(combined_data, ignore_index=True)
        
        # Ensure required columns exist
        required_cols = ['Log2(Fold Change)', 'CD Coefficient']
        for col in required_cols:
            if col not in combined_df.columns:
                st.error(f"Missing required column: {col}")
                return None
        
        # Remove rows with NaN values in critical columns
        combined_df = combined_df.dropna(subset=required_cols)
        
        if combined_df.empty:
            st.warning("No valid data points for volcano plot after filtering")
            return None
        
        # Create volcano plot
        hover_cols = ['Perturbagen']
        optional_cols = ['Dose', 'Cell Line']
        for col in optional_cols:
            if col in combined_df.columns:
                hover_cols.append(col)
        
        fig = px.scatter(
            combined_df, 
            x='Log2(Fold Change)', 
            y='CD Coefficient',
            color='Regulation',
            hover_data=hover_cols,
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
        
    except Exception as e:
        logger.error(f"Error creating volcano plot: {str(e)}", exc_info=True)
        st.error(f"❌ Error creating volcano plot: {str(e)}")
        return None

def create_dose_response_plot(df: pd.DataFrame, regulation_type: str, gene_name: str) -> Optional[go.Figure]:
    """Create dose-response analysis plot with error handling"""
    try:
        if df.empty:
            return None
            
        # Extract numeric dose values where possible
        df_copy = df.copy()
        
        if 'Dose' not in df_copy.columns:
            return None
            
        # Try to extract numeric values from dose column
        df_copy['Dose_Numeric'] = pd.to_numeric(
            df_copy['Dose'].astype(str).str.extract(r'(\d+\.?\d*)')[0], 
            errors='coerce'
        )
        
        # Filter out rows where dose couldn't be parsed
        dose_df = df_copy.dropna(subset=['Dose_Numeric', 'CD Coefficient'])
        
        if len(dose_df) < 5:  # Need minimum data points
            return None
        
        # Ensure we have required columns
        if 'Fold Change' not in dose_df.columns:
            dose_df['Fold Change'] = 1  # Default value
        
        hover_cols = ['Perturbagen']
        if 'Cell Line' in dose_df.columns:
            hover_cols.append('Cell Line')
        
        fig = px.scatter(
            dose_df.head(50),  # Top 50 for readability
            x='Dose_Numeric',
            y='CD Coefficient',
            size='Fold Change',
            hover_data=hover_cols,
            title=f'Dose-Response Analysis: {regulation_type} of {gene_name}',
            labels={'Dose_Numeric': 'Dose (μM)', 'CD Coefficient': 'CD Coefficient'},
            log_x=True
        )
        fig.update_layout(template='plotly_white')
        return fig
        
    except Exception as e:
        logger.error(f"Error creating dose-response plot: {str(e)}", exc_info=True)
        return None

def convert_df_to_csv(df: pd.DataFrame) -> Optional[bytes]:
    """Convert DataFrame to CSV for download with error handling"""
    try:
        if df.empty:
            st.warning("No data to download")
            return None
        return df.to_csv(index=False).encode('utf-8')
    except Exception as e:
        logger.error(f"Error converting DataFrame to CSV: {str(e)}", exc_info=True)
        st.error(f"❌ Error preparing download: {str(e)}")
        return None

def create_cell_line_analysis(df: pd.DataFrame, gene_name: str, regulation_type: str) -> Tuple[Optional[go.Figure], Optional[pd.DataFrame]]:
    """Analyze effects across different cell lines with error handling"""
    try:
        if df.empty or 'Cell Line' not in df.columns:
            return None, None
            
        # Ensure numeric columns
        df_copy = df.copy()
        numeric_cols = ['CD Coefficient', 'Fold Change']
        for col in numeric_cols:
            if col in df_copy.columns:
                df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce')
        
        cell_line_stats = df_copy.groupby('Cell Line').agg({
            'CD Coefficient': ['mean', 'count'],
            'Fold Change': 'mean'
        }).round(4)
        
        cell_line_stats.columns = ['Mean_CD_Coeff', 'Count', 'Mean_Fold_Change']
        cell_line_stats = cell_line_stats.reset_index()
        cell_line_stats = cell_line_stats[cell_line_stats['Count'] >= 3].sort_values('Mean_CD_Coeff', ascending=False)
        
        if len(cell_line_stats) == 0:
            return None, None
            
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
        
    except Exception as e:
        logger.error(f"Error creating cell line analysis: {str(e)}", exc_info=True)
        return None, None

def create_timepoint_analysis(df: pd.DataFrame, gene_name: str, regulation_type: str) -> Tuple[Optional[go.Figure], Optional[pd.DataFrame]]:
    """Analyze temporal effects with error handling"""
    try:
        if df.empty or 'Timepoint' not in df.columns:
            return None, None
            
        # Ensure numeric columns
        df_copy = df.copy()
        numeric_cols = ['CD Coefficient', 'Fold Change']
        for col in numeric_cols:
            if col in df_copy.columns:
                df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce')
        
        timepoint_stats = df_copy.groupby('Timepoint').agg({
            'CD Coefficient': ['mean', 'count'],
            'Fold Change': 'mean'
        }).round(4)
        
        timepoint_stats.columns = ['Mean_CD_Coeff', 'Count', 'Mean_Fold_Change']
        timepoint_stats = timepoint_stats.reset_index()
        timepoint_stats = timepoint_stats[timepoint_stats['Count'] >= 5].sort_values('Mean_CD_Coeff', ascending=False)
        
        if len(timepoint_stats) == 0:
            return None, None
            
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
        
    except Exception as e:
        logger.error(f"Error creating timepoint analysis: {str(e)}", exc_info=True)
        return None, None

def create_top_compounds_plot(df: pd.DataFrame, gene_name: str, regulation_type: str, top_n: int = 20) -> Optional[go.Figure]:
    """Show top compounds by effectiveness with error handling"""
    try:
        if df.empty:
            return None
            
        # Ensure we have required columns
        if 'CD Coefficient' not in df.columns or 'Perturbagen' not in df.columns:
            st.error("Missing required columns for top compounds plot")
            return None
        
        top_compounds = df.head(top_n)
        
        hover_cols = []
        optional_cols = ['Dose', 'Cell Line', 'Fold Change']
        for col in optional_cols:
            if col in top_compounds.columns:
                hover_cols.append(col)
        
        fig = px.bar(
            top_compounds,
            x='CD Coefficient',
            y='Perturbagen',
            orientation='h',
            title=f'Top {top_n} Most Effective Compounds: {regulation_type} of {gene_name}',
            labels={'CD Coefficient': 'CD Coefficient', 'Perturbagen': 'Compound'},
            hover_data=hover_cols if hover_cols else None
        )
        fig.update_layout(
            height=max(400, top_n * 25),
            template='plotly_white'
        )
        return fig
        
    except Exception as e:
        logger.error(f"Error creating top compounds plot: {str(e)}", exc_info=True)
        return None

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
if query_button:
    if not gene_name or not gene_name.strip():
        st.error("⚠️ Please enter a gene name to query")
    elif not validate_gene_name(gene_name):
        st.error("⚠️ Invalid gene name format. Please use alphanumeric characters only.")
    elif not get_up and not get_down:
        st.error("⚠️ Please select at least one regulation type (up or down)")
    else:
        gene_upper = gene_name.strip().upper()
        
        with st.spinner(f'Querying data for {gene_upper}...'):
            try:
                # Initialize containers for results
                results = {}
                
                # Query up-regulating compounds
                if get_up:
                    with st.status("Fetching up-regulating compounds..."):
                        up_df = lincs_cp(gene_name, "up", top_n)
                        if up_df is not None and not up_df.empty:
                            results['up'] = up_df
                            st.write(f"✅ Found {len(up_df)} up-regulating compounds")
                        else:
                            st.write("⚠️ No up-regulating compounds found")
                
                # Query down-regulating compounds  
                if get_down:
                    with st.status("Fetching down-regulating compounds..."):
                        down_df = lincs_cp(gene_name, "down", top_n)
                        if down_df is not None and not down_df.empty:
                            results['down'] = down_df
                            st.write(f"✅ Found {len(down_df)} down-regulating compounds")
                        else:
                            st.write("⚠️ No down-regulating compounds found")
                
                # Check if we have any results
                if not results:
                    st.warning(f"No data found for gene {gene_upper}. Please check the gene symbol and try again.")
                else:
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
                            volcano_fig = create_volcano_plot(
                                results.get('up', pd.DataFrame()), 
                                results.get('down', pd.DataFrame()), 
                                gene_upper
                            )
                            if volcano_fig:
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
                                    if top_up_fig:
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
                                    if top_down_fig:
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
                                    df_up = results['up']
                                    st.metric("Total Compounds", len(df_up))
                                    if 'CD Coefficient' in df_up.columns and not df_up['CD Coefficient'].isna().all():
                                        st.metric("Max CD Coefficient", f"{df_up['CD Coefficient'].max():.4f}")
                                    if 'Fold Change' in df_up.columns and not df_up['Fold Change'].isna().all():
                                        st.metric("Mean Fold Change", f"{df_up['Fold Change'].mean():.3f}")
                                    if 'Cell Line' in df_up.columns:
                                        st.metric("Unique Cell Lines", df_up['Cell Line'].nunique())
                                    if 'Timepoint' in df_up.columns:
                                        st.metric("Unique Timepoints", df_up['Timepoint'].nunique())
                            
                            if 'down' in results:
                                with col2:
                                    st.markdown("### 📉 Down-regulation Summary")
                                    df_down = results['down']
                                    st.metric("Total Compounds", len(df_down))
                                    if 'CD Coefficient' in df_down.columns and not df_down['CD Coefficient'].isna().all():
                                        st.metric("Min CD Coefficient", f"{df_down['CD Coefficient'].min():.4f}")
                                    if 'Fold Change' in df_down.columns and not df_down['Fold Change'].isna().all():
                                        st.metric("Mean Fold Change", f"{df_down['Fold Change'].mean():.3f}")
                                    if 'Cell Line' in df_down.columns:
                                        st.metric("Unique Cell Lines", df_down['Cell Line'].nunique())
                                    if 'Timepoint' in df_down.columns:
                                        st.metric("Unique Timepoints", df_down['Timepoint'].nunique())
                    
                    # Handle single result type
                    elif len(results) == 1:
                        result_type = list(results.keys())[0]
                        direction_name = "Up-regulation" if result_type == 'up' else "Down-regulation"
                        st.subheader(f"📊 {direction_name} Analysis for {gene_upper}")
                        
                        # Create plots for single result
                        col1, col2 = st.columns(2)
                        with col1:
                            top_fig = create_top_compounds_plot(results[result_type], gene_upper, direction_name)
                            if top_fig:
                                st.plotly_chart(top_fig, use_container_width=True)
                        
                        with col2:
                            dose_fig = create_dose_response_plot(results[result_type], direction_name, gene_upper)
                            if dose_fig:
                                st.plotly_chart(dose_fig, use_container_width=True)
                    
                    # Data Tables Section
                    st.markdown("---")
                    st.header("📊 Raw Data Tables")
                    
                    # Create tabs for up and down regulation data tables
                    if len(results) > 1:
                        data_tab_up, data_tab_down = st.tabs(["📈 Up-regulating Data", "📉 Down-regulating Data"])
                    
                    # Up-regulating results
                    if 'up' in results:
                        with (data_tab_up if len(results) > 1 else st.container()):
                            st.subheader(f"📈 Compounds that UP-regulate {gene_upper}")
                            st.dataframe(results['up'], use_container_width=True, height=400)
                            
                            # Download button
                            csv_up = convert_df_to_csv(results['up'])
                            if csv_up:
                                st.download_button(
                                    label="💾 Download Up-regulating Data",
                                    data=csv_up,
                                    file_name=f"{gene_upper}_Up-regulating_perturbations.csv",
                                    mime="text/csv",
                                    key="download_up"
                                )
                    
                    # Down-regulating results
                    if 'down' in results:
                        with (data_tab_down if len(results) > 1 else st.container()):
                            st.subheader(f"📉 Compounds that DOWN-regulate {gene_upper}")
                            st.dataframe(results['down'], use_container_width=True, height=400)
                            
                            # Download button
                            csv_down = convert_df_to_csv(results['down'])
                            if csv_down:
                                st.download_button(
                                    label="💾 Download Down-regulating Data",
                                    data=csv_down,
                                    file_name=f"{gene_upper}_Down-regulating_perturbations.csv",
                                    mime="text/csv",
                                    key="download_down"
                                )
                
            except Exception as e:
                logger.error(f"Unexpected error in main execution: {str(e)}", exc_info=True)
                st.error(f"❌ An unexpected error occurred: {str(e)}")
                st.error("Please try again or contact support if the problem persists.")

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
        
        **Troubleshooting:**
        - Ensure gene symbols are correct (e.g., use "TP53" not "p53")
        - Check spelling and capitalization
        - Some genes may not have data in the LINCS database
        """)

# Update requirements.txt note
st.sidebar.markdown("---")
st.sidebar.markdown("**💡 Features:**")
st.sidebar.markdown("• Comprehensive error handling")
st.sidebar.markdown("• Automatic retry on timeouts")
st.sidebar.markdown("• Input validation") 
st.sidebar.markdown("• Robust data processing")
