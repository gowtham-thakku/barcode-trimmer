#!/usr/bin/env python3
"""
Streamlit web application for barcode trimming.
Save this as streamlit_app.py and run with: streamlit run streamlit_app.py
"""

import streamlit as st
import tempfile
import json
from datetime import datetime
from pathlib import Path
import io
import zipfile
import time

# Import your filtering function (adjust import as needed)
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    import parasail
    DEPENDENCIES_AVAILABLE = True
except ImportError:
    DEPENDENCIES_AVAILABLE = False
    st.warning("⚠️ Warning: BioPython and/or parasail not available. Using simplified version.")

# Configure page
st.set_page_config(
    page_title="Barcode Trimmer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main > div {
        padding-top: 2rem;
    }
    
    .stAlert > div {
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        border-radius: 0.5rem;
    }
    
    .upload-section {
        border: 2px dashed #ccc;
        border-radius: 10px;
        padding: 20px;
        text-align: center;
        margin: 10px 0;
    }
    
    .results-section {
        background-color: #f8f9fa;
        padding: 20px;
        border-radius: 10px;
        border: 1px solid #e9ecef;
        margin: 20px 0;
    }
</style>
""", unsafe_allow_html=True)

# Your filtering functions (same as before)
def load_adapters(fasta_content):
    """Return list of adapter sequences + reverse complements (upper-case)."""
    adapters = []
    lines = fasta_content.strip().split('\n')
    
    current_seq = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                seq = current_seq.upper()
                adapters.append(seq)
                adapters.append(str(Seq(seq).reverse_complement()))
            current_seq = ''
        else:
            current_seq += line
    
    # Handle last sequence
    if current_seq:
        seq = current_seq.upper()
        adapters.append(seq)
        adapters.append(str(Seq(seq).reverse_complement()))
    
    return adapters

def file_format_from_name(filename):
    """Infer file format from extension: 'fastq' or 'fasta'."""
    ext = Path(filename).suffix.lower()
    if ext in {".fa", ".fasta"}:
        return "fasta"
    if ext in {".fq", ".fastq"}:
        return "fastq"
    return "fastq"  # default

def filter_reads_web(adapter_content, reads_content, reads_filename, params, progress_placeholder=None):
    """
    Web version of your filter_reads_parasail function.
    Returns (kept_content, discarded_content, log_content).
    """
    
    if not DEPENDENCIES_AVAILABLE:
        # Fallback implementation without parasail
        return simple_filter_fallback(adapter_content, reads_content, reads_filename, params, progress_placeholder)
    
    # Load adapters
    adapters = load_adapters(adapter_content)
    
    # Build scoring matrix
    matrix = parasail.matrix_create("ACGT", params['match'], params['mismatch'])
    
    # Parse reads
    fmt = file_format_from_name(reads_filename)
    reads_io = io.StringIO(reads_content)
    
    # First pass to count total reads for progress tracking
    total_reads = 0
    for _ in SeqIO.parse(io.StringIO(reads_content), fmt):
        total_reads += 1
    
    # Second pass for actual processing
    reads_io = io.StringIO(reads_content)
    kept_reads = []
    discarded_reads = []
    processed = 0
    
    for rec in SeqIO.parse(reads_io, fmt):
        processed += 1
        seq = str(rec.seq).upper()
        contaminated = False
        
        for adapter in adapters:
            if parasail.sw_scan_16(seq, adapter, 
                                 params['gap_open'], params['gap_extend'], matrix).score >= params['min_score']:
                contaminated = True
                break
        
        if contaminated:
            discarded_reads.append(rec)
        else:
            kept_reads.append(rec)
        
        # Update progress every 100 reads
        if progress_placeholder and processed % 100 == 0:
            progress_placeholder.text(f"Processing reads: {processed:,}/{total_reads:,}")
    
    # Final progress update
    if progress_placeholder:
        progress_placeholder.text(f"Processing reads: {total_reads:,}/{total_reads:,}")
    
    # Generate output content
    kept_io = io.StringIO()
    discarded_io = io.StringIO()
    
    SeqIO.write(kept_reads, kept_io, fmt)
    SeqIO.write(discarded_reads, discarded_io, fmt)
    
    kept_content = kept_io.getvalue()
    discarded_content = discarded_io.getvalue()
    
    # Generate log
    log_content = f"""Mid-Read Barcode Trimming Log
Started: {datetime.now().isoformat()}
Tool: Parasail SIMD Smith-Waterman Filter

Parameters:
{json.dumps(params, indent=2)}

Input File: {reads_filename}
Format: {fmt.upper()}
Total reads: {total_reads:,}
Adapters loaded: {len(adapters) // 2}

Results:
Reads kept: {len(kept_reads):,}
Reads discarded: {len(discarded_reads):,}
Completed: {datetime.now().isoformat()}
"""
    
    return kept_content, discarded_content, log_content

def simple_filter_fallback(adapter_content, reads_content, reads_filename, params, progress_placeholder=None):
    """Simplified filtering when parasail is not available."""
    
    # Parse adapters (simplified)
    adapters = []
    lines = adapter_content.strip().split('\n')
    current_seq = ''
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                adapters.append(current_seq.upper())
            current_seq = ''
        else:
            current_seq += line
    if current_seq:
        adapters.append(current_seq.upper())
    
    # Simple substring-based filtering
    fmt = file_format_from_name(reads_filename)
    reads_lines = reads_content.strip().split('\n')
    
    kept_lines = []
    discarded_lines = []
    total = 0
    kept_count = 0
    discarded_count = 0
    processed = 0
    
    if fmt == 'fastq':
        total_reads = len(reads_lines) // 4
        for i in range(0, len(reads_lines), 4):
            if i + 3 < len(reads_lines):
                processed += 1
                seq = reads_lines[i + 1].strip().upper()
                
                # Simple contamination check
                contaminated = any(adapter in seq for adapter in adapters)
                
                read_block = reads_lines[i:i+4]
                if contaminated:
                    discarded_lines.extend(read_block)
                    discarded_count += 1
                else:
                    kept_lines.extend(read_block)
                    kept_count += 1
                
                # Update progress every 100 reads
                if progress_placeholder and processed % 100 == 0:
                    progress_placeholder.text(f"Processing reads: {processed:,}/{total_reads:,}")
        
        total = processed
        # Final progress update
        if progress_placeholder:
            progress_placeholder.text(f"Processing reads: {total:,}/{total:,}")
            
    else:  # FASTA
        # First pass to count reads
        total_reads = sum(1 for line in reads_lines if line.strip().startswith('>'))
        
        current_header = ''
        current_seq = ''
        
        for line in reads_lines:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    processed += 1
                    contaminated = any(adapter in current_seq.upper() for adapter in adapters)
                    
                    if contaminated:
                        discarded_lines.extend([current_header, current_seq])
                        discarded_count += 1
                    else:
                        kept_lines.extend([current_header, current_seq])
                        kept_count += 1
                    
                    # Update progress every 100 reads
                    if progress_placeholder and processed % 100 == 0:
                        progress_placeholder.text(f"Processing reads: {processed:,}/{total_reads:,}")
                
                current_header = line
                current_seq = ''
            else:
                current_seq += line
        
        # Handle last sequence
        if current_header and current_seq:
            processed += 1
            contaminated = any(adapter in current_seq.upper() for adapter in adapters)
            
            if contaminated:
                discarded_lines.extend([current_header, current_seq])
                discarded_count += 1
            else:
                kept_lines.extend([current_header, current_seq])
                kept_count += 1
        
        total = processed
        # Final progress update
        if progress_placeholder:
            progress_placeholder.text(f"Processing reads: {total:,}/{total:,}")
    
    kept_content = '\n'.join(kept_lines) + '\n' if kept_lines else ''
    discarded_content = '\n'.join(discarded_lines) + '\n' if discarded_lines else ''
    
    log_content = f"""Mid-Read Barcode Trimming Log (Simplified Mode)
Started: {datetime.now().isoformat()}
Tool: Simple substring-based filter

Parameters:
{json.dumps(params, indent=2)}

Input File: {reads_filename}
Format: {fmt.upper()}
Total reads: {total:,}
Adapters loaded: {len(adapters)}

Results:
Reads kept: {kept_count:,}
Reads discarded: {discarded_count:,}
Completed: {datetime.now().isoformat()}

Note: This is a simplified version. Install BioPython and parasail for full Smith-Waterman alignment.
"""
    
    return kept_content, discarded_content, log_content

def create_download_zip(kept_content, discarded_content, log_content, reads_filename):
    """Create a ZIP file with all results."""
    zip_buffer = io.BytesIO()
    fmt = file_format_from_name(reads_filename)
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        zip_file.writestr(f'filtered_reads.{fmt}', kept_content)
        zip_file.writestr(f'discarded_reads.{fmt}', discarded_content)
        zip_file.writestr('filtering_log.txt', log_content)
    
    return zip_buffer.getvalue()

# Main Streamlit App
def main():
    # Header
    st.title("🧬 Mid-Read Barcode Trimmer")
    st.markdown("### Filter sequencing reads using Smith-Waterman alignment against adapter/barcode sequences")
    
    # Info box
    with st.container():
        st.info("""
        **How it works:** This tool uses Parasail SIMD Smith-Waterman alignment to identify and filter reads 
        containing adapter/barcode sequences. Reads with alignment scores ≥ threshold are discarded as contaminated. If using for Oxford Nanopore data, use the default parameters.
        """)
    
    # Check dependencies status
    if DEPENDENCIES_AVAILABLE:
        st.success("✅ Full functionality available with BioPython and parasail")
    else:
        st.warning("⚠️ Running in simplified mode (BioPython/parasail not available)")
    
    # Create two columns for file uploads
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### 📄 Upload Sequencing Reads")
        reads_file = st.file_uploader(
            "Choose your reads file",
            type=['fastq', 'fq', 'fasta', 'fa'],
            help="Supported formats: FASTQ, FASTA (.fastq, .fq, .fasta, .fa)"
        )
        
        if reads_file:
            st.success(f"✅ Loaded: {reads_file.name} ({reads_file.size:,} bytes)")
    
    with col2:
        st.markdown("#### 🔬 Upload Adapter/Barcode FASTA")
        adapter_file = st.file_uploader(
            "Choose your adapter FASTA file",
            type=['fasta', 'fa'],
            help="FASTA format with adapter sequences (.fasta, .fa)"
        )
        
        if adapter_file:
            st.success(f"✅ Loaded: {adapter_file.name} ({adapter_file.size:,} bytes)")
    
    # Parameters section
    st.markdown("#### ⚙️ Smith-Waterman Parameters")
    
    # Create columns for parameters
    param_cols = st.columns(5)
    
    with param_cols[0]:
        min_score = st.number_input("Min Score Threshold", value=30, min_value=1, help="Minimum alignment score to consider a read contaminated")
    
    with param_cols[1]:
        match = st.number_input("Match Score", value=2, help="Score for nucleotide matches")
    
    with param_cols[2]:
        mismatch = st.number_input("Mismatch Penalty", value=-1, help="Penalty for nucleotide mismatches")
    
    with param_cols[3]:
        gap_open = st.number_input("Gap Open Penalty", value=5, min_value=0, help="Penalty for opening a gap")
    
    with param_cols[4]:
        gap_extend = st.number_input("Gap Extend Penalty", value=1, min_value=0, help="Penalty for extending a gap")
    
    # Process button
    st.markdown("---")
    
    if st.button("🚀 Filter Reads", type="primary", use_container_width=True):
        if not reads_file or not adapter_file:
            st.error("❌ Please upload both files before processing")
        else:
            # Create parameters dictionary
            params = {
                'min_score': min_score,
                'match': match,
                'mismatch': mismatch,
                'gap_open': gap_open,
                'gap_extend': gap_extend
            }
            
            # Show processing message with progress placeholder
            progress_placeholder = st.empty()
            with st.spinner('🔄 Processing your files... This may take a few minutes.'):
                try:
                    # Read file contents
                    reads_content = reads_file.read().decode('utf-8')
                    adapter_content = adapter_file.read().decode('utf-8')
                    
                    # Process files with progress tracking
                    kept_content, discarded_content, log_content = filter_reads_web(
                        adapter_content, reads_content, reads_file.name, params, progress_placeholder
                    )
                    
                    # Clear progress display
                    progress_placeholder.empty()
                    
                    # Store results in session state
                    st.session_state.results = {
                        'kept_content': kept_content,
                        'discarded_content': discarded_content,
                        'log_content': log_content,
                        'reads_filename': reads_file.name,
                        'processed_at': datetime.now()
                    }
                    
                    st.success("✅ Processing complete!")
                    st.rerun()  # Refresh to show results
                    
                except Exception as e:
                    st.error(f"❌ Error processing files: {str(e)}")
    
    # Results section
    if 'results' in st.session_state:
        st.markdown("---")
        st.markdown("### 🎉 Results")
        
        results = st.session_state.results
        
        # Parse log to extract statistics
        log_lines = results['log_content'].split('\n')
        stats = {}
        for line in log_lines:
            if 'Total reads:' in line:
                stats['total'] = line.split(':')[1].strip()
            elif 'Reads kept:' in line:
                stats['kept'] = line.split(':')[1].strip()
            elif 'Reads discarded:' in line:
                stats['discarded'] = line.split(':')[1].strip()
        
        # Display statistics
        if stats:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Reads", stats.get('total', 'N/A'))
            with col2:
                st.metric("Reads Kept", stats.get('kept', 'N/A'), delta=None)
            with col3:
                st.metric("Reads Discarded", stats.get('discarded', 'N/A'), delta=None)
        
        # Download section
        st.markdown("#### 📥 Download Results")
        
        # Create download buttons
        col1, col2, col3, col4 = st.columns(4)
        
        fmt = file_format_from_name(results['reads_filename'])
        
        with col1:
            st.download_button(
                label="📥 Filtered Reads",
                data=results['kept_content'],
                file_name=f"filtered_reads.{fmt}",
                mime="text/plain",
                help="Download reads that passed the filter"
            )
        
        with col2:
            st.download_button(
                label="🗑️ Discarded Reads",
                data=results['discarded_content'],
                file_name=f"discarded_reads.{fmt}",
                mime="text/plain",
                help="Download reads that were filtered out"
            )
        
        with col3:
            st.download_button(
                label="📋 Log File",
                data=results['log_content'],
                file_name="filtering_log.txt",
                mime="text/plain",
                help="Download processing log with parameters and statistics"
            )
        
        with col4:
            # Create ZIP download
            zip_data = create_download_zip(
                results['kept_content'], 
                results['discarded_content'], 
                results['log_content'],
                results['reads_filename']
            )
            st.download_button(
                label="📦 Download All",
                data=zip_data,
                file_name="barcode_trimmer_results.zip",
                mime="application/zip",
                help="Download all results in a ZIP file"
            )
        
        # Show log preview
        with st.expander("📋 View Processing Log"):
            st.text(results['log_content'])
        
        # Clear results button
        if st.button("🔄 Process New Files", help="Clear current results to process new files"):
            if 'results' in st.session_state:
                del st.session_state.results
            st.rerun()

if __name__ == "__main__":
    main()
