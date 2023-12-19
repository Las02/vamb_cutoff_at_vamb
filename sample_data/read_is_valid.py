import sys
import os
print(sys.argv)

for path in sys.argv[1:]:
    print(path,':')
    path_append = ''
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line.split()) == 1:
                path_append = line
            else:
                line_split = line.split()   
                for sub_path in line_split[1:]:
                    full_path =path_append + sub_path 
                    if not os.path.isfile(full_path):
                        print(full_path)


