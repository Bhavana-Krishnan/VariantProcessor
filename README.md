# Variant Processor

This repository is focused on developing a tool called Variant Processor. This repo contains a Python script "VariantProcessor.py" designed in Object Oriented Programmimg style to derive meaningful insights from Variant Calling Files in vcf format by processing them. It provides functionalities to process variants and gain meaningful insights by writing processed data to CSV files, counting different types of origins, identifying mutation types, writing mutation files based on Substitution, Insertion or Deletion, and counting ref-alt pairs for substitution variants.

To learn more about VCF files, please refer to: [VCF-Variant-Call-Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format)

## Features

- **Variant Processing**: Parses variant calling files in VCF format and extracts relevant information.
- **Data Export**: Writes processed variant data to a CSV file with specified columns.
- **Origin Counting**: Counts occurrences of different origin types in the processed data.
- **Mutation Identification**: Identifies and counts different types of mutations (Substitution, Insertion, Deletion) based on reference and alternate alleles.
- **Mutation File Generation**: Writes mutation files (Substitution.tsv, Insertion.tsv, Deletion.tsv) containing records classified under different mutation types.
- **Ref-Alt Pair Counting**: Counts occurrences of reference and alternate allele pairs for substitution variants.

## Usage

1. **Input File**: Provide the path to the input variant calling file in VCF format. An example input vcf file is given in this repository in the input folder.
2. **Output File**: Specify the path for the output CSV file to store processed data. Example output files obtained can be viewed in the output folder.
3. **Working Directory**: Provide the path to the Working Directory. It could be the same path where the input and output folders can be placed.
3. **Processing**: Run the script to process variants and generate the output CSV file.
4. **Analysis**:
   - Use the generated CSV file to analyze variant data.
   - Count different types of origins using the `count_origins()` method.
   - Identify mutation types and their counts using the `identify_mutation_types()` method.
   - Write mutation files for further analysis using the `write_mutation_files()` method.
   - Count ref-alt pairs for substitution variants using the `count_ref_alt()` method.

## Dependencies

- Python 3.x

## Installation

1. Clone this repository to your local machine:

    ```bash
    git clone https://github.com/yourusername/variant-processor.git
    ```

2. Navigate to the project directory:

    ```bash
    cd variant-processor
    ```

3. Make the necessary modifications:
    
    -Modify the paths for input file, output file and Working Directory.
    - Run the script using the command "python VariantProcessor.py"

    ```bash
    python VariantProcessor.py
    ```

## Example

```python
# Example usage of VariantProcessor
if __name__ == "__main__":
    input_file = ".\\VariantProcessor\\input\\exampleVCF.vcf"
    output_file = ".\\VariantProcessor\\output\\Exampleprocessed_data.csv"
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
    substitution_variants = [variant for variant in processor.variants if len(variant.ref) == len(variant.alt) == 1]
    ref_alt_counts = processor.count_ref_alt(substitution_variants)
    print(ref_alt_counts)
```