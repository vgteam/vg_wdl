import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-w', help='WDL file', required=True)

args = parser.parse_args()

wdl_f = open(args.w, 'rt')

print_param = False

for line in wdl_f:
    if print_param:
        if '}' in line:
            break
        else:
            line = line.rstrip().split(':')
            par_n = line[0].lstrip()
            par_d = ':'.join(line[1:])
            par_d = par_d.rstrip('"').lstrip().lstrip('"')
            print("- *{}*: {}".format(par_n, par_d))
    elif 'parameter_meta' in line:
        print_param = True

wdl_f.close()
