#!/usr/bin/env python3



import sys

def concatenator(hydro, average, aacomp):
    with open(hydro, 'r') as inhydro, open(average, 'r') as inaverage, open(aacomp, 'r') as inaacomp:
        #print(inhydro, inaverage, inaacomp)
        for line in inhydro:
            identifier = line.split()[0] + ' ' + line.split()[1]
            avg_line = inaverage.readline()
            aacomp_line = inaacomp.readline()
            if identifier in avg_line and identifier in aacomp_line:
                #print(line.rstrip(), avg_line.rstrip(), aacomp_line.rstrip())
                print(identifier, line.split()[2], aacomp_line.split()[3], line.split()[4], line.split()[3], end=' [ ')
                for elem in avg_line.split()[3:]:
                    if '[' in elem:
                        break
                    else:
                        print(elem, end=' ')
                print(']', end=' [ ')
                for freq in aacomp_line.split()[4:-1]:
                    print(freq, end=' ')
                print(']')
            else:
                print('The protein id ', identifier, ' found in <', hydro, '> is not present in either <', average, ' or ', aacomp, '>', file=sys.stderr)
                raise SystemExit

if __name__ == "__main__":
    try:
        hydro_input = sys.argv[1]
        average_input = sys.argv[2]
        aacomp_input = sys.argv[3]
    except:
        print('Program usage: text.py', file=sys.stderr)
        raise SystemExit
    else:
        concatenator(hydro_input, average_input, aacomp_input)
