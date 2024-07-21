
from argparse import ArgumentParser
from meta.sample_data.qiime import split_and_dump_main_metadata_table


def parse_args():
    parser = ArgumentParser(
        description="Split QIIME2 main metadata file",
        epilog=""
    )
    parser.add_argument("-i", "--input", required=True, help="Main metadata file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    _namespace = parser.parse_args()
    return (
        _namespace.input,
        _namespace.output,
    )


if __name__ == '__main__':
    input_file, output_dir = parse_args()
    split_and_dump_main_metadata_table(input_file, output_dir)
