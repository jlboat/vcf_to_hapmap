use std::env;
use std::process;
use std::fs::File;
use std::io::Write;
use vcf::VCFReader;
use bstr::ByteSlice;
use std::io::BufReader;
use std::io::BufWriter;
use flate2::read::MultiGzDecoder;

pub(crate) fn vcf_to_hapmap(input_filename: &str, output_filename: &str) {
    let mut reader = VCFReader::new(
        BufReader::new(
            MultiGzDecoder::new(
                File::open(input_filename,
                ).expect("Failed to read file.")))).expect("Failed to read VCF");

    let hapmap_vector: [String; 11] = ["rs".to_string(), "alleles".to_string(), "chrom".to_string(),
        "pos".to_string(), "strand".to_string(), "assembly".to_string(),
        "center".to_string(), "protLSID".to_string(),
        "assayLSID".to_string(), "panelLSID".to_string(), "QCcode".to_string()];

    let output = File::create(output_filename)
        .expect("Failed to create output file");
    let mut stream = BufWriter::new(output);
    let buffer_out = hapmap_vector.join("\t");
    let bytes_amount = stream.write(buffer_out.as_ref()).expect("Failed to write header");
    if buffer_out.len() != bytes_amount {
        panic!("Failed to write header to BufWriter");
    }

    for sample in reader.header().samples() {
        let buffer_out = format!("\t{}", sample.to_str()
            .expect("Failed to convert sample to string")
        );
        let bytes_amount = stream.write(buffer_out.as_ref())
            .expect("Failed to write to output");
        if buffer_out.len() != bytes_amount {
            panic!("Failed to write header to BufWriter");
        }
    }
    let bytes_amount = stream.write(b"\n")
        .expect("Failed to write newline after header");
    if b"\n".len() != bytes_amount {
        panic!("Failed to write header to BufWriter");
    }

    let mut vcf_record = reader.empty_record();
    while reader.next_record(&mut vcf_record).expect("Failed to read record") {
        let reference = String::from_utf8(vcf_record.reference.to_vec()
        ).expect("Failed ref");
        let alternative = String::from_utf8(vcf_record.alternative
            .to_vec().get(0).expect("Failed alt").to_vec()
        ).expect("Failed to convert alt");
        let chromosome = String::from_utf8(vcf_record.chromosome.to_vec())
            .expect("Failed to convert chromosome to string [snp name]");
        let buffer_out = format!("{}_{}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                             &chromosome,
                             &vcf_record.position.to_string(),
                             &reference,
                             &alternative,
                             &chromosome.replace("Chr0", "").replace("Chr", ""),
                             &vcf_record.position.to_string(),
                             "+", "NA", "NA", "NA", "NA", "NA", "NA"
        );
        let bytes_amount = stream.write(buffer_out.as_ref())
            .expect("Failed to write to output");
        if buffer_out.len() != bytes_amount {
            panic!("Failed to write marker to BufWriter");
        }

        for sample in reader.header().samples() {
            let genotype = String::from_utf8(
                vcf_record.genotype(sample, b"GT")
                    .expect("Failed to get genotype")
                    .to_vec().get(0).expect("").to_vec())
                .expect("Failed to convert genotype to string");
            let bytes_amount = stream.write(b"\t").expect("Genotype tab failed");
            if b"\t".len() != bytes_amount{
                panic!("Failed to write genotype tab");
            }

            let split_genotype;
            if genotype.contains('/'){
                split_genotype = genotype.split('/');
            } else if genotype.contains('|'){
                split_genotype = genotype.split('|');
            } else {
                panic!("Uncertain genotype encoding");
            }
            for item in split_genotype {
                if item.eq("0") {
                    let bytes_amount = stream.write((&reference).as_ref())
                        .expect("Failed to write reference allele to output");
                    if reference.len() != bytes_amount{
                        panic!("Ref allele write failed");
                    }
                } else if item.eq("1") {
                    let bytes_amount = stream.write((&alternative).as_ref())
                        .expect("Failed to write alternative allele to output");
                    if alternative.len() != bytes_amount{
                        panic!("Alt allele write failed");
                    }
                }
            }
        }
        let bytes_amount = stream.write(b"\n").expect("");
        if b"\n".len() != bytes_amount{
            panic!("Failed to write newline");
        }

    }
}

fn main(){
    let args: Vec<String> = env::args().collect();
    if args.len() < 3{
        println!("\nConvert large gzipped VCF to hapmap");
        println!("vcf_to_hapmap Input.vcf.gz Output.hapmap.txt\n");
        process::exit(1);
    }
    let input = &args[1];
    let output = &args[2];
    vcf_to_hapmap(input, output);
}
