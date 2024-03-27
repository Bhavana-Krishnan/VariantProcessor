import csv
from collections import defaultdict


class Variant:
    def __init__(self, chrom, pos, vid, ref, alt, qual, filt, info):
        # Initializing variant attributes
        self.chrom = chrom
        self.pos = pos
        self.vid = vid
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filt = filt
        self.info = info

    @classmethod
    def from_line(cls, line):
        # Creates Variant instance from a line of text
        fields = line.strip().split("\t")
        return cls(*fields[:8])

    def process_info(self):
        # Extracts relevant information from the INFO field
        info_dict = dict(item.split("=") for item in self.info.split(";"))
        allele_id = info_dict.get("ALLELEID", "-")
        clnhgvs = info_dict.get("CLNHGVS", "-")
        clnsig = info_dict.get("CLNSIG", "-")
        clnvc = info_dict.get("CLNVC", "-")
        origin = info_dict.get("ORIGIN", "-")
        rs = info_dict.get("RS", "-")
        gene_info = info_dict.get("GENEINFO", "-").split(":")
        gene_id = gene_info[1] if len(gene_info) > 1 else "-"
        gene_symbol = gene_info[0] if len(gene_info) > 0 else "-"
        consequence = "-"
        if "MC" in info_dict:
            consequence = info_dict["MC"].split("|")[1]
        return allele_id, clnhgvs, clnsig, clnvc, origin, rs, gene_id, gene_symbol, consequence


class VariantProcessor:
    def __init__(self, input_file):
        # Initialize VariantProcessor with input file
        self.input_file = input_file
        self.variants = []

    def process_variants(self):
        # Process variants from input file
        with open(self.input_file, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            for row in reader:
                if row and not row[0].startswith("#"):
                    variant = Variant(*row[:8])
                    self.variants.append(variant)

    def write_processed_data_to_csv(self, output_file):
        # Write processed variant data to CSV file
        with open(output_file, "w") as file:
            writer = csv.writer(file)
            writer.writerow(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNHGVS',
                            'CLNSIG', 'CLNVC', 'ORIGIN', 'RS', 'Gene_ID', 'Gene_symbol', 'Consequence'])
            for variant in self.variants:
                allele_id, clnhgvs, clnsig, clnvc, origin, rs, gene_id, gene_symbol, consequence = variant.process_info()
                writer.writerow([variant.chrom, variant.pos, variant.vid, variant.ref, variant.alt,
                                allele_id, clnhgvs, clnsig, clnvc, origin, rs, gene_id, gene_symbol, consequence])

    def count_origins(self):
        # Count occurrences of different origin types
        origins = {0: 'unknown',
                   1: 'germline',
                   2: 'somatic',
                   4: 'inherited',
                   8: 'paternal',
                   16: 'maternal',
                   32: 'de-novo',
                   64: 'biparental',
                   128: 'uniparental',
                   256: 'not-tested',
                   512: 'tested-inconclusive'
                   }
        origin_counts = defaultdict(int)
        for variant in self.variants:
            origin = variant.process_info()[4]
            origin_counts[origin] += 1
        mut_counts = {}
        for k1, v1 in origins.items():
            for k2, v2 in origin_counts.items():
                # print(k1, k2)
                if str(k1) == k2:
                    mut_counts[v1] = v2
        return mut_counts

    def identify_mutation_types(self):
        # Identify different mutation types
        mutation_types = {'Substitution': 0, 'Insertion': 0, 'Deletion': 0}
        for variant in self.variants:
            ref_length = len(variant.ref)
            alt_length = len(variant.alt)
            if ref_length == 1 and alt_length == 1:
                mutation_types['Substitution'] += 1
            elif alt_length > ref_length:
                mutation_types['Insertion'] += 1
            elif alt_length < ref_length:
                mutation_types['Deletion'] += 1
        return mutation_types

    def write_mutation_files(self):
        # Write mutation files for different mutation types
        substitution_variants = []
        insertion_variants = []
        deletion_variants = []
        for variant in self.variants:
            ref_length = len(variant.ref)
            alt_length = len(variant.alt)
            if ref_length == 1 and alt_length == 1:
                substitution_variants.append(variant)
            elif alt_length > ref_length:
                insertion_variants.append(variant)
            elif alt_length < ref_length:
                deletion_variants.append(variant)
        self._write_variant_file(
            substitution_variants, working_dir + 'Substitution.tsv')
        self._write_variant_file(
            insertion_variants, working_dir + 'Insertion.tsv')
        self._write_variant_file(
            deletion_variants, working_dir + 'Deletion.tsv')

    def _write_variant_file(self, variants, filename):
        # Helper method to write variant file
        with open(filename, "w") as file:
            writer = csv.writer(file, delimiter="\t")
            for variant in variants:
                writer.writerow([variant.chrom, variant.pos,
                                variant.vid, variant.ref, variant.alt])

    def count_ref_alt(self, variants):
        # Count occurrences of ref-alt pairs for substitution variants
        ref_alt_counts = defaultdict(int)
        for variant in variants:
            ref = variant.ref
            alt = variant.alt
            ref_alt_counts[(ref, alt)] += 1
        return ref_alt_counts


if __name__ == "__main__":
    input_file = ".\\VariantProcessor\\input\\exampleVCF.vcf"
    output_file = ".\\VariantProcessor\\output\\processed_data.csv"
    working_dir = ".\\VariantProcessor\\output\\"
    processor = VariantProcessor(input_file)
    processor.process_variants()
    processor.write_processed_data_to_csv(output_file)

    # Count origins
    origin_counts = processor.count_origins()
    print(origin_counts)

    # Identify mutation types
    mutation_types = processor.identify_mutation_types()
    print(mutation_types)

    # Write mutation files
    processor.write_mutation_files()

    # Count ref-alt pairs for substitution file
    substitution_variants = [variant for variant in processor.variants if len(
        variant.ref) == len(variant.alt) == 1]
    ref_alt_counts = processor.count_ref_alt(substitution_variants)
    print(ref_alt_counts)
