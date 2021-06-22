import sys


with open(sys.argv[1],'r') as file:
  print(file.readline().strip().split('\t'))