# vcf_to_hapmap
Convert large gzipped VCF to double-bit HapMap format

This Rust program was primarily designed to be a more memory-efficient solution than TASSEL to convert a large VCF to a double-bit HapMap for GAPIT

```
git clone https://github.com/jlboat/vcf_to_hapmap
cd vcf_to_hapmap

cargo build --release
./vcf_to_hapmap/target/release/vcf_to_hapmap 

## prints the following help
Convert large gzipped VCF to hapmap
vcf_to_hapmap Input.vcf.gz Output.hapmap.txt
```
