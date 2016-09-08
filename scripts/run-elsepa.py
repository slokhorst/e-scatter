from elsepa import units
from elsepa.generate_input import (generate_elscata_input, Settings)
from elsepa.parse_output import (elsepa_output_parsers)
from elsepa.executable import (DockerContainer, Archive)

import numpy as np
import re

energies = np.array([
    # 10, 20, 30, 40, 50, 60, 70, 80, 90,
    # 100, 200, 300, 400, 500, 600, 700, 800, 900,
    # 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
    # 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
    10, 100, 100000]) * units.eV


def elscata(settings: Settings):
    with DockerContainer('elsepa', working_dir='/opt/elsepa') as elsepa:
        elsepa.put_archive(
            Archive('w')
            .add_text_file('input.dat', generate_elscata_input(s))
            .close())

        elsepa.sh('./elscata < input.dat',
                  'mkdir result && mv *.dat result')

        result_tar = elsepa.get_archive('result')
        result = {}

        for info in result_tar:
            if not info.isfile():
                continue

            name = re.match("result/(.*)\\.dat", info.name).group(1)
            parser = elsepa_output_parsers[name]

            if parser:
                lines = result_tar.get_text_file(info).split('\n')
                result[name] = parser(lines)

    return result


if __name__ == "__main__":
    s = Settings(IZ=80, MNUCL=3, MELEC=4, MUFFIN=0, IELEC=-1,
                 MEXCH=1, MCPOL=2, IHEF=0, MABS=0, EV=energies)

    result = elscata(s)
    print(result['tcstable'])
