import sys
import os
print(sys.argv)

for path in sys.argv[1:]:
    print(path,':')
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not os.path.isfile(line):
                print(line)


