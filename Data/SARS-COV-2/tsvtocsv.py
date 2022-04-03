import re
import os
import glob


def main():
    files = glob.glob('./*/*.txt')

    for f in files:
        with open(f, 'r') as tsv_file:
            out = f.replace('.txt', '.csv')
            if os.path.exists(out):
                os.remove(out)
            with open(out, 'a') as csv_file:
                fileContent = ""
                for line in tsv_file:
                    if line.startswith("TYPE"):
                        continue
                    fileContent += re.sub('\t', ',', line)
                csv_file.write(fileContent)

if __name__ == '__main__':
    main()