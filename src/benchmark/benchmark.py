#!/usr/bin/env python3

import subprocess
import argparse

BIN_DIR = '../../bin/'
ROOT_DIR = '../../'

class Benchmark:
    def __init__(self, program, n_runs, dataset_size, test, show_diff, n_procs) -> None:
        if 'deriche' in program:
            self.type = 'deriche'
        elif 'seidel-2d' in program:
            self.type = 'seidel-2d'
        elif 'heat-3d' in program:
            self.type = 'heat-3d'

        try:
            self.n_runs = int(n_runs)
        except TypeError:
            self.n_runs = 1

        self.program = program
        self.dataset_size = dataset_size
        self.test = test
        self.show_diff = show_diff
        self.n_procs = n_procs
        

    def make(self):
        subprocess.run(['make', 'clean'], cwd=ROOT_DIR, capture_output=True)
        command = ['make', 'DERICHE_DIM=64', f'SIZE={self.dataset_size}']
        if self.test:
            command += 'DUMP=-DPOLYBENCH_DUMP_ARRAYS'
        p = subprocess.run(command, cwd=ROOT_DIR, capture_output=True)


    def check_implementation(self):
        self.make()
        p = subprocess.run(['./'+self.program], capture_output=True, text=True, cwd=BIN_DIR)

        with open('diff_output', 'w') as text_file:
            text_file.write(p.stderr)

        diff = subprocess.run(['diff', 'diff_output', ROOT_DIR + 'test/' + self.type + '/' + self.type + '_' + self.dataset_size + '_dataset.out'],
            capture_output=True, text=True)

        if self.show_diff:
            if (diff.stdout != ''):
                print(diff.stdout)
        else:
            print('Correct output' if diff.stdout == '' else 'Wrong output')
        
    
    def bench(self):
        # ref impl
        P = ROOT_DIR + 'polybench/medley/deriche/'
        subprocess.run(['make', 'clean'], cwd=P, capture_output=True)
        subprocess.run(['make', f'EXTRA_FLAGS=-D{self.dataset_size}_DATASET'], cwd=P, capture_output=True)
        runtimes = []
        for _ in range(self.n_runs):
            ref_impl = subprocess.run(['./deriche'], cwd=P, capture_output=True)
            runtimes.append(float(ref_impl.stdout))
        average = sum(runtimes) / len(runtimes)
        print('Ref impl runtime:', average)

        # bench impl
        self.make()
        averages = []
        for np in self.n_procs:
            runtimes = []
            for _ in range(self.n_runs):
                p = subprocess.run(['mpiexec', '-np', np, f'./{self.program}'], capture_output=True, text=True, cwd=BIN_DIR)
                try:
                    time = float(p.stdout)
                except ValueError:
                    print('Program quit unexpectedly:')
                    print(p.stderr)
                    print(p.stdout)

                runtimes.append(time)
            average = sum(runtimes) / len(runtimes)
            averages.append(average)
        print('Runtimes:', averages)
        print('With #cores:', self.n_procs)

        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('-s', '--dataset-size', required=True)
    parser.add_argument('-n', '--runs', required=False)
    parser.add_argument('-np', '--procs', required=False, action='append')
    parser.add_argument('-d', '--show-diff', required=False, action='store_true')
    parser.add_argument('-t', '--test', required=False, action='store_true')
   
    args = parser.parse_args()
    program = args.executable
    n_runs = args.runs
    dataset_size = args.dataset_size
    test = args.test
    show_diff = args.show_diff
    n_procs = args.procs

    bench = Benchmark(program, n_runs, dataset_size, test, show_diff, n_procs)
    if test:
        bench.check_implementation()
    else:
        bench.bench()


if __name__ == '__main__':
    main()
