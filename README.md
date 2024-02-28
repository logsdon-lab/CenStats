# Centromere Status Checker
Determine the status of centromeric contigs based on RepeatMasker annotations.

### Setup
```bash
pip install
```

### Usage
Takes RepeatMasker output from centromeric contigs to check and centromeric contigs from a reference.
```bash
cens-status -i input_cens.out -r ref_cens.out
```

Return tab-delimited output with the following fields:
1. Original contig name.
2. Remapped contig name.
3. Correct orientation of contig with `fwd` indicating an already correctly oriented contig.
4. Contig is a partial centromeric contig.

```
HG00171_chr21_haplotype2-0000154:33297526-38584922      HG00171_chr21_haplotype2-0000154:33297526-38584922      rev     false
HG00171_chr21_haplotype1-0000029:4196961-10792258       HG00171_chr13_haplotype1-0000029:4196961-10792258       fwd     false
HG00171_chr21_haplotype1-0000019:33327260-37313906      HG00171_chr21_haplotype1-0000019:33327260-37313906      rev     true
chm13_chr21:7700001-11850000    chm13_chr21:7700001-11850000    fwd     false
```

### Build
```bash
make build && pip install dist/*.whl
```
