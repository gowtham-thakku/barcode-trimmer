# barcode-trimmer*
barcode-trimmer filters sequencing reads that contain barcodes or adapters using Smith–Waterman alignment. It helps clean raw data before downstream analysis. It was specifically designed to remove mid-read barcodes from Oxford Nanopore data to address the technical issue of barcode bleeding. It is however generalizable for other types of inputs/filtering.

## Usage Options

### Option 1: Web App (Easiest)
**Access at: [barcode-trimmer.streamlit.app](https://barcode-trimmer.streamlit.app)**

**Pros:**
- No installation required
- Works in any browser

**Cons:**
- **File size limit: 200MB per file** (Streamlit Community Cloud restriction)
- Processing may be slower on shared resources

**Best for:** Small to medium-sized files, quick analyses, testing the tool

---

### Option 2: Local Installation (For Large Files)
**Run on your own computer with custom limits**

**Pros:**
- **File size limit: 2GB+ (configurable up to your system's memory)**
- Faster processing on dedicated hardware
- Full control over resources
- Private - your data never leaves your computer

**Cons:**
- Local installation needed (but only 3-4 lines of command line code! A local instance of the app will be launched on your browser)

**Best for:** Large files (>200MB), sensitive data, high-performance processing

### Installation (local)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/barcode-trimmer.git
   cd barcode-trimmer
2. **Install dependencies:**
   ```bash
   pip install streamlit biopython parasail
3. **Run with desired upload limit:**
   ```bash
   # For 5GB limit
   streamlit run barcode_trimmer_streamlit_app.py --server.maxUploadSize=5120

   # For 10GB limit  
   streamlit run barcode_trimmer_streamlit_app.py --server.maxUploadSize=10240


## How to Use

1. **Upload your sequencing reads file** (FASTQ or FASTA)
2. **Upload your adapter/barcode FASTA file** (list of Oxford Nanopore barcodes used in the Rapid PCR Barcoding kit is provided in the repo: ont_rlb_barcodes.fasta)
3. **Adjust parameters** if needed (defaults work well for most cases)
4. **Click "Filter Reads"** and wait for processing
5. **Download your results** - filtered reads, discarded reads, or everything as ZIP

### Parameters Explained

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Min Score Threshold** | 30 | Minimum alignment score to mark read as contaminated |
| **Match Score** | 2 | Points awarded for nucleotide matches |
| **Mismatch Penalty** | -1 | Points deducted for mismatches |
| **Gap Open Penalty** | 5 | Cost to start a gap in alignment |
| **Gap Extend Penalty** | 1 | Cost to extend an existing gap |

*Note: This code was generated using the assistance of Claude Code, but included specific prompts to which bioinformatics tools to use (e.g Parasail); the code has been validated using multiple datasets. Default thresholds were optimized using these datasets.
