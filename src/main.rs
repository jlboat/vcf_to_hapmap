use std::env;
use std::process;
use std::fs::File;
use std::io::Write;
use vcf::VCFReader;
use bstr::ByteSlice;
use std::io::BufReader;
use std::io::BufWriter;
use flate2::read::MultiGzDecoder;

fn open_gzvcf(filename: &str) -> VCFReader<BufReader<MultiGzDecoder<File>>> {
    let vcf = VCFReader::new(
        BufReader::new(
            MultiGzDecoder::new(
                File::open(filename,
                ).expect("Failed to read file.")))).expect("Failed to read VCF");
    return vcf
}

fn write_check(mut stream: BufWriter<File>, output: &[u8]) -> BufWriter<File>{
    let bytes_amount = stream.write(output)
        .expect("Failed to write to hapmap");
    if output.len() != bytes_amount {
        panic!("Failed to write to BufWriter");
    }
    return stream
}

fn get_header_prefix() -> String {
    let hapmap_vector: [String; 11] = ["rs".to_string(), "alleles".to_string(), "chrom".to_string(),
        "pos".to_string(), "strand".to_string(), "assembly".to_string(),
        "center".to_string(), "protLSID".to_string(),
        "assayLSID".to_string(), "panelLSID".to_string(), "QCcode".to_string()];
    return hapmap_vector.join("\t")
}

pub(crate) fn vcf_to_hapmap(input_filename: &str, output_filename: &str) {
    let mut reader = open_gzvcf(input_filename);
    let output = File::create(output_filename)
        .expect("Failed to create output file");
    let mut stream = BufWriter::new(output);
    stream = write_check(stream, get_header_prefix().as_ref());

    for sample in reader.header().samples() {
        stream = write_check(stream, b"\t");
        stream = write_check(stream, sample.to_str().expect("").as_ref());
    }
    stream = write_check(stream, b"\n");

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
        stream = write_check(stream, buffer_out.as_ref());

        for sample in reader.header().samples() {
            let genotype = String::from_utf8(
                vcf_record.genotype(sample, b"GT")
                    .expect("Failed to get genotype")
                    .to_vec().get(0).expect("").to_vec())
                .expect("Failed to convert genotype to string");
            stream = write_check(stream, b"\t");

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
                    stream = write_check(stream, (reference).as_ref());
                } else if item.eq("1") {
                    stream = write_check(stream, (alternative).as_ref());
                }
            }
        }
        stream = write_check(stream, b"\n");
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
