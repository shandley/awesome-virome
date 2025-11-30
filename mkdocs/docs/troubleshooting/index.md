# Troubleshooting Guide

> **Last Updated:** November 29, 2025

This guide helps you diagnose and fix common problems in virome analysis. Issues are organized by workflow stage.

## Quick Problem Finder

**Jump to your issue:**

- Installation Problems
- Quality Control Issues
- Assembly Problems
- Viral Identification Issues
- CheckV Errors
- Performance Issues

## Installation Problems

### Tool Won't Install via Conda
```bash
# Update conda
conda update -n base conda

# Try mamba (faster)
conda install mamba
mamba install -c bioconda tool_name

# Create fresh environment
conda create -n fresh_env python=3.9
conda install -c bioconda tool_name
```

### Database Download Fails
```bash
# Check disk space
df -h

# Resume download
wget -c URL

# Use different download location
export TMPDIR=/scratch/tmp
```

## Quality Control Issues

### Low QC Pass Rate (<70%)

**Solutions:**
```bash
# Relax quality threshold
fastp ... --qualified_quality_phred 15

# Lower length requirement
fastp ... --length_required 30
```

### High Duplication (>60%)

**Normal for viromes** (20-70% expected)

If problematic:
```bash
clumpify.sh in=reads.fq out=dedup.fq dedupe
```

## Assembly Problems

### Out of Memory

**Solutions:**
```bash
# Increase memory
spades.py -m 250

# Subsample reads
seqtk sample reads.fq 5000000 > subset.fq

# Fewer k-mer sizes
spades.py -k 21,33,55
```

### Low N50 (<3 kb)

**Causes:** Low coverage, high diversity, poor quality

**Solutions:**
- Sequence deeper
- Try different k-mer sizes
- Improve QC

## Viral Identification Issues

### Very Few Viruses

**If <1%:**
```bash
# Lower threshold
virsorter run --min-score 0.3 ...

# Use multiple tools
# VirSorter2 + VIBRANT + geNomad

# Check if viral reads present
blastn -query contigs.fa -db viral_refseq
```

### Too Many Viruses (>50%)

**Validate with CheckV:**
```bash
checkv end_to_end predicted.fa checkv_out/

# Check contamination
awk -F'\t' '$10 > 10' quality_summary.tsv | wc -l
```

## CheckV Errors

### High Contamination

**Solutions:**
```bash
# Use cleaned sequences
cp checkv_out/viruses.fna cleaned.fa

# Or filter
awk -F'\t' '$10 < 10 {print $1}' quality_summary.tsv > clean_ids.txt
seqkit grep -f clean_ids.txt contigs.fa > filtered.fa
```

## Host Prediction Issues

### No Predictions

**Expected:** Only 30-60% get predictions

**Improve:**
```bash
# Lower threshold
iphop predict --min_score 70 ...

# Add more MAGs from environment
# Use multiple methods (CRISPR + WIsH + ML)
```

## Performance Issues

### Job Killed

**Solutions:**
```bash
# Request more resources
#SBATCH --mem=250G
#SBATCH --time=48:00:00

# Parallelize
for sample in S*; do
    process.sh $sample &
done
wait
```

### Analysis Too Slow

**Optimize:**
```bash
# More threads
tool -t 32

# Use faster tools
diamond  # instead of blastp
bowtie2  # instead of bwa

# Subsample for pilot
seqtk sample reads.fq 1000000 > pilot.fq
```

## Common Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| `std::bad_alloc` | Out of memory | Increase RAM |
| `Killed` | Out of memory | Request more in job |
| `Permission denied` | File permissions | Check with `ls -la` |
| `No such file` | Wrong path | Verify path exists |
| `Command not found` | Not installed | Install tool, activate env |

## Getting Help

**If problems persist:**

1. Check tool GitHub issues
2. Search error message on Google/Biostars
3. [Awesome-Virome Discussions](https://github.com/shandley/awesome-virome/discussions)
4. Ask on [Biostars](https://biostars.org)

**When asking for help, include:**
- Complete error message
- Tool version
- Input file characteristics
- Steps to reproduce

## Prevention Checklist

- [ ] Test on small dataset first
- [ ] Record all commands
- [ ] Check input formats
- [ ] Monitor resources (RAM, disk)
- [ ] Use version control
- [ ] Include controls
- [ ] Validate each step

## See Also

- [Tutorial Troubleshooting Sections](../tutorials/basic-metagenome-virome.md#troubleshooting)
- [Best Practices: Quality Control](../best-practices/quality-control.md)
- [Tool Selection Guide](../tools/selection-guide.md)
