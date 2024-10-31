import sys

cur_name = None
cur_len = None
for line in sys.stdin:
    if line.startswith(">"):
        if cur_name:
            print(cur_name, cur_len)
        cur_name = line.strip()
        cur_len = 0
    else:
        cur_len += len(line.strip())
print(cur_name, cur_len)
