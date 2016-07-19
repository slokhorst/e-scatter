#!/usr/bin/env python3

import noodles
import subprocess
import json

from itertools import (chain)
from noodles.display import SimpleDisplay

@noodles.has_scheduled_methods
class Command(object):
    def __init__(self, name, executable, *fixed_args, **kwargs):
        self.name = name
        self.executable = executable
        self.fixed_args = fixed_args
        self.parameters = kwargs

    @noodles.schedule_hint(
        display='{self.name} {kwargs}',
        confirm=True)
        # annotated=True)
    def __call__(self, **kwargs):
        def generate_argument(k, v):
            if isinstance(v, bool):
                return [self.parameters[k]]
            else:
                return [self.parameters[k], str(v)]

        command = list(chain(
            [self.executable], 
            self.fixed_args,
            *(generate_argument(k, v) for k, v in kwargs.items())))

        p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
            universal_newlines=True)
        p.check_returncode()

        result = json.loads(p.stdout)
        return result # , p.stderr

def error_filter(xcptn):
    if isinstance(xcptn, subprocess.CalledProcessError):
        return xcptn.stderr
    else:
        return None

if __name__ == '__main__':
    command = Command('e⁻Scatter', './escat', 'sp', '-no-png', '-machine',
        depth='-n', 
        multiplicity='-m',
        scaling='-x', 
        ball_size='-s',
        hilbert='-H',
        random='-R')

    sizes = [i * 0.05 for i in range(20)]
    a = command(depth=10, multiplicity=3, scaling=0.5, ball_size=0.01, hilbert=True)
    b = command(depth=10, multiplicity=3, scaling=0.5, ball_size=0.01, random=True)

    work = noodles.gather(a, b)
    with SimpleDisplay(error_filter) as display:
        result = noodles.run_logging(work, 1, display)

    for line in result:
        print(json.dumps(line, indent=2))


