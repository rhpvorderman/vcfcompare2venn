# Copyright (c) 2022 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import argparse
import os
from typing import TextIO

from matplotlib import pyplot as plt
import matplotlib_venn


def parse_label_filename(label: str):
    filename, percentage = label.split()
    return filename


def subsets_from_vcf_compare_file(file: TextIO):
    lines = [line for line in file if line.startswith("VN\t")]
    if len(lines) != 7:
        raise NotImplementedError("Only venn3 diagrams are implemnted")
    singles = []
    duos = []
    all = []
    for line in lines:
        spl = line.split("\t")
        count = int(spl[1])
        filenames = [parse_label_filename(x) for x in spl[2:]]
        if len(filenames) == 1:
            singles.append((filenames[0], count))
        elif len(filenames) == 2:
            duos.append((filenames, count))
        elif len(filenames) == 3:
            all.append((filenames, count))
        else:
            raise RuntimeError("Corrupt vcf-compare output.")
    if len(singles) != 3 or len(duos) != 3 or len(all) != 1:
        raise RuntimeError("Corrupt vcf-compare output")
    duos_dict = {frozenset(filenames): count for filenames, count in duos}
    _, all_count = all[0]
    singles.sort()
    a_name, a_count = singles[0]
    b_name, b_count = singles[1]
    c_name, c_count = singles[2]
    a_and_b = frozenset((a_name, b_name))
    a_and_c = frozenset((a_name, c_name))
    b_and_c = frozenset((b_name, c_name))
    subsets = (
        a_count,
        b_count,
        duos_dict[a_and_b],
        c_count,
        duos_dict[a_and_c],
        duos_dict[b_and_c],
        all_count
    )
    labels = (
        os.path.basename(a_name),
        os.path.basename(b_name),
        os.path.basename(c_name)
    )
    return subsets, labels


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", metavar="INPUT",
                        help="A vcf-compare output file")
    parser.add_argument("output", metavar="PLOT", nargs="?",
                        help="Save location of the plot. If not given shows "
                             "the interactive matplotlib interface.")
    parser.add_argument("--title", "-t",
                        help="Title for the plot.")
    return parser


def main():
    args = argument_parser().parse_args()
    with open(args.input) as f:
        subsets, lables = subsets_from_vcf_compare_file(f)
    matplotlib_venn.venn3(subsets=subsets, set_labels=lables,)
    plt.title(args.title)
    if not args.output:
        plt.show()
    else:
        plt.savefig(args.output)


if __name__ == "__main__":
    main()
