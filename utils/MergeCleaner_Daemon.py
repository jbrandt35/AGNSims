import os
from time import sleep
from Setup import print_to_stdout

for i in range(0, 30):
    print_to_stdout(f"Initializing Merge Cleaner for {i}th time...")
    os.system("python MergeCleaner.py")
    sleep(60 * 60)
