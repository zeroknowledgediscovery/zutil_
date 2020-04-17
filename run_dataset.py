
import subprocess as sp
import sys

run_command = './dtw ntb/seq.dat 10 dtw.dst'
sp.call('./dtw ntb/seq.dat 10 dtw.dst',shell=True)