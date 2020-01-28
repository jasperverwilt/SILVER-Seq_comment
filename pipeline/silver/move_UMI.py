import sys
import itertools



line_no = 0
with open(sys.argv[1],"r") as myfile, open(sys.argv[2], "r") as umifile:
    for line, umi in itertools.zip_longest(myfile, umifile):
        line = line.rstrip()
        umi = umi.rstrip()
        if line_no == 1 :
            print(umi + line)
        elif line_no == 2 :
            print("+")
        elif line_no == 3 :
            print(umi + line)
        else :
            print(line)
            line_no = 0
        line_no += 1

