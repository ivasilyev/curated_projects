#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description="Generate sampledata based on directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Input directory (directories)")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-r", "--regex", default=DEFAULT_REGEX,
                        help="Regular expression to extract sample name(s) from file name(s) without extension")
    parser.add_argument("-t", "--taxa", default="",
                        help="(Optional) Taxonomy to add to all samples")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.extension, _namespace.regex,  _namespace.taxa, _namespace.output


if __name__ == '__main__':
    inputDirs, inputExtension, inputRegex, inputTaxa, outputFile = parse_args()

    pair2dArray = []
    for input_dir in inputDirs:
        pair2dArray.extend(Utilities.get_most_similar_word_pairs(Utilities.find_file_by_tail(
            input_dir, inputExtension, multiple=True)))

    sampleDataArray = SampleDataArray.generate(pair2dArray, inputRegex)

    for sampleDataLine in sampleDataArray.values:
        sampleDataLine.taxa = inputTaxa

    sampleDataArray.dump(outputFile)
"""


from meta.utils.io import dump_dict
from meta.utils.language import regex_based_tokenization
from meta.utils.file_system import find_file_by_tail


def tokenize_reads_file_name(s: str):
    d = regex_based_tokenization({
        "extension": ["\.(.{2,8})$", "(\..{2,8})$"],  # E.g. '.fastq.gz'
        "last_segment": ["[^A-Za-z0-9]([A-Za-z0-9]+)$", "([^A-Za-z0-9][A-Za-z0-9]+)$"],  # The last segment is always 001,
        "read_index": ["[^A-Za-z0-9](R[0-9]+)$", "([^A-Za-z0-9]R[0-9]+)$"],
        "lane_number": ["[^A-Za-z0-9](L[0-9]+)$", "([^A-Za-z0-9]L[0-9]+)$"],
        "sample_sheet_number": ["[^A-Za-z0-9](S[0-9]+)$", "([^A-Za-z0-9]S[0-9]+)$"],
        "sample_name": ["(.+)", "(.+)"],
    }, os.path.basename(s))
    d["reads_file"] = s
    return d


def create_sampledata_dict_from_list(reads_files: list):
    tokenized_reads_files = [
        tokenize_reads_file_name(i) for i in sorted(sorted(reads_files), key=len, reverse=True)
    ]
    out = dict()
    for token_dict in tokenized_reads_files:
        sample_name = token_dict["sample_name"]
        if sample_name in out.keys():
            out[sample_name]["reads"].append(token_dict["reads_file"])
        else:
            out[sample_name] = {
                "name": sample_name,
                "reads": [token_dict["reads_file"], ],
                "taxa": ""
            }
        out[sample_name]["reads"] = sorted(out[sample_name]["reads"])
    return out


def create_sampledata_dict_from_dir(directory: str, reads_extension: str = DEFAULT_READS_EXTENSION):
    reads_files = find_file_by_tail(
        directory, ".{}".format(reads_extension.strip(".")), multiple=True
    )
    return create_sampledata_dict_from_list(reads_files)


def parse_args():
    parser = ArgumentParser(
        description="Generate sample data based on scanned raw reads files in the directory".strip(),
        epilog=""
    )
    parser.add_argument("-i", "--input", nargs="+",
                        help="Input directory (directories)")
    parser.add_argument("-e", "--extension", default=DEFAULT_READS_EXTENSION,
                        help="Extension of reads files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file")
    _namespace = parser.parse_args()
    return _namespace.input, _namespace.extension, _namespace.output


if __name__ == '__main__':
    inputDirs, inputExtension, outputFile = parse_args()

    sampledata_dict = dict()
    for input_dir in inputDirs:
        sampledata_dict.update(create_sampledata_dict_from_dir(input_dir, inputExtension))

    dump_dict(sampledata_dict, outputFile)
    print(f"Sampledata successfully created in: '{outputFile}'")
