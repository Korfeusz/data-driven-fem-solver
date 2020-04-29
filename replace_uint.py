import os


def inplace_change(in_folder, out_folder, filename, old_string, new_string):
    with open('{}/{}'.format(in_folder, filename), "rt") as fin:
        with open('{}/{}'.format(out_folder, filename), "wt") as fout:
            for line in fin:
                fout.write(line.replace(old_string, new_string))


in_f = 'test'
out_f = 'test_corrected'
for f in os.listdir(in_f):
    print(f)
    inplace_change(in_f, out_f, f, 'UInt32', 'Int32')
