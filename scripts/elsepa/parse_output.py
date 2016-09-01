import pandas as pd
import re


def arg_first(pred, s):
    return next(i for i, v in enumerate(s) if pred(v))


def extract_header(l1_, l2_):
    m = max(len(l1_), len(l2_))
    l1 = l1_.ljust(m+1)
    l2 = l2_.ljust(m+1)
    c = list(zip(l1, l2))

    x1 = x2 = 0
    while True:
        x1 = x2 + arg_first(lambda v: v[0] != ' ' or v[1] != ' ', c[x2:])
        x2 = x1 + arg_first(lambda v: v[0] == ' ' and v[1] == ' ', c[x1:])
        yield ' '.join([l1[x1:x2].strip(), l2[x1:x2].strip()])


def parse_most_elscata_output(lines):
    comments = []
    i = 0
    for line in lines:
        if line.startswith(" #----"):
            if comments[-2].strip() != '':
                h = list(extract_header(comments[-2], comments[-1]))
            else:
                h = [a.strip() for a in comments[-1].split('  ') if a]
            data = pd.DataFrame(columns=h)
        elif line.startswith(" #"):
            comments.append(line[2:])
        else:
            values = [float(v) for v in line.split()]
            if len(values) < len(h):
                continue
            i += 1
            data.loc[i] = values
    return data, comments


class RegexDict(dict):
    """A dictionary that uses regular expressions to match a
    requested key to a value."""
    def __init__(self, *args):
        super(RegexDict, self).__init__(*args)

    def __getitem__(self, query):
        for key, value in self.items():
            if re.fullmatch(key, query):
                return value
        raise KeyError("None of the patterns matched the query.")


elsepa_output_parsers = RegexDict({
    "dpwa":     None,
    "input":    None,
    "dcs_.*":   parse_most_elscata_output,
    "scatamp":  parse_most_elscata_output,
    "scfield":  parse_most_elscata_output,
    "tcstable": parse_most_elscata_output
})
