# barcode-trimmer
Web app for filtering sequencing reads using Smith-Waterman alignment

# 🧬 Mid-Read Barcode Trimmer

A web application for filtering sequencing reads using Smith-Waterman alignment against adapter/barcode sequences, designed for CZID.org analysis preparation.

## 🚀 Quick Start - Use Online

**👉 [Launch the App](https://your-app-name.streamlit.app)** *(No installation required!)*

## 📖 What it does

This tool helps bioinformaticians and researchers clean up their sequencing data by:

- **Identifying contaminated reads** using advanced Smith-Waterman alignment
- **Filtering out adapter/barcode sequences** that can interfere with downstream analysis
- **Preparing clean data** for platforms like CZID.org
- **Providing detailed logs** of the filtering process

## ✨ Features

- 🌐 **Web-based interface** - No command line required
- 🚀 **Fast processing** with Parasail SIMD optimization
- 📊 **Real-time statistics** showing filtering results
- 📥 **Multiple download options** (individual files or ZIP bundle)
- 🔧 **Customizable parameters** for different use cases
- 📱 **Mobile-friendly** responsive design

## 🛠️ Supported Formats

### Input Files
- **Sequencing reads**: FASTQ (.fastq, .fq), FASTA (.fasta, .fa)
- **Adapters/Barcodes**: FASTA format (.fasta, .fa)

### Output Files
- **Filtered reads**: Same format as input
- **Discarded reads**: Same format as input  
- **Processing log**: Detailed text report

## 🎯 How to Use

1. **Upload your sequencing reads file** (FASTQ or FASTA)
2. **Upload your adapter/barcode FASTA file**
3. **Adjust parameters** if needed (defaults work well for most cases)
4. **Click "Filter Reads"** and wait for processing
5. **Download your results** - filtered reads, discarded reads, or everything as ZIP

### Parameters Explained

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Min Score Threshold** | 25 | Minimum alignment score to mark read as contaminated |
| **Match Score** | 2 | Points awarded for nucleotide matches |
| **Mismatch Penalty** | -1 | Points deducted for mismatches |
| **Gap Open Penalty** | 5 | Cost to start a gap in alignment |
| **Gap Extend Penalty** | 1 | Cost to extend an existing gap |

## 🖥️ Local Installation (Optional)

If you prefer to run locally or need to process very large files:

### Prerequisites
```bash
pip install streamlit biopython parasail
