#!/usr/bin/env python3

import subprocess
import argparse

BIN_DIR = '../../bin/'
ROOT_DIR = '../../'

class Benchmark:
    def __init__(self, program, n_runs, dataset_size, test, show_diff) -> None:
        if 'deriche' in program:
            self.type = 'deriche'
        elif 'seidel-2d' in program:
            self.type = 'seidel-2d'
        elif 'heat-3d' in program:
            self.type = 'heat-3d'
        self.program = BIN_DIR + program
        try:
            self.n_runs = int(n_runs)
        except TypeError:
            self.n_runs = 1
        self.dataset_size = dataset_size
        self.test = test
        self.show_diff = show_diff

    def make(self):
        p = subprocess.run(['make', f'SIZE={self.dataset_size}', 'DUMP=' + '-DPOLYBENCH_DUMP_ARRAYS' if self.test else ''],
            cwd=ROOT_DIR, capture_output=True)

    def check_implementation(self):
        self.make()
        p = subprocess.run([self.program], capture_output=True, text=True)

        with open('diff_output', 'w') as text_file:
            text_file.write(p.stderr)

        diff = subprocess.run(['diff', 'diff_output', ROOT_DIR + 'test/' + self.type + '/' + self.type + '_' + self.dataset_size + '_dataset.out'],
            capture_output=True, text=True)

        if self.show_diff:
            if (diff.stdout != ''):
                print(diff.stdout)
        else:
            print('Correct output' if diff.stdout == '' else 'Wrong output')
        
    
    def create_csv(self):
        self.make()
        # TODO
        raise NotImplementedError()
        runtimes = []
        for _ in range(self.n_runs):
            p = subprocess.run([*self.program.split(' ')], capture_output=True, text=True)
            time = float(p.stdout)
            output = p.stderr
            runtimes.append(time)

        average = sum(runtimes) / len(runtimes)
        print(average)
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('-s', '--dataset-size', required=True)
    parser.add_argument('-n', '--n-runs', required=False)
    parser.add_argument('-d', '--show-diff', required=False, action='store_true')
    parser.add_argument('-t', '--test', required=False, action='store_true')
   
    args = parser.parse_args()
    program = args.executable
    n_runs = args.n_runs
    dataset_size = args.dataset_size
    test = args.test
    show_diff = args.show_diff

    bench = Benchmark(program, n_runs, dataset_size, test, show_diff)
    if test:
        bench.check_implementation()
    else:
        bench.create_csv()


if __name__ == '__main__':
    main()
