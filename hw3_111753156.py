#<陳羽暉, 111753156>

import argparse
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="inputpath, pls enter a path")
parser.add_argument("--score", type=str, help="inputscorepath, pls enter a path")
parser.add_argument("--aln", type=str, help="inputaln, pls enter global or local")
parser.add_argument("--gap", type=int, help="input gap score, pls enter a score")
parser.add_argument("--output", type=str, help="outputpath, pls enter a path")
args = parser.parse_args()
print(vars(args))

# Read file
file = open(args.score, "r")
data = file.readlines()
file.close()

pamx = np.zeros((24, 24))
for i in range(10, 34):
    s = data[i].rstrip()
    rowData = s.split()
    del rowData[0]
    for j in range(len(rowData)):
        pamx[i - 10][j] = int(rowData[j])
aminoAcid = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']

file = open(args.input, 'r')
ls = []
name = []
for line in file:
    if line.startswith('>'):
        name.append(line)
    if not line.startswith('>'):
        ls.append(line.replace('\n', ''))
file.close()

#Caculate the dynamic programming
gapScore = args.gap
dp_score = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
dp_direction = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))

direction = ['Left', 'Upper Left', 'Up']

dp_score[0][0] = 1
for i in range(1, len(ls[0]) + 1):
    dp_score[i][0] = dp_score[i-1][0] + gapScore
for j in range(1, len(ls[1]) + 1):
    dp_score[0][j] = dp_score[0][j-1] + gapScore

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        aminoIndex1 = aminoAcid.index(ls[0][i - 1])
        aminoIndex2 = aminoAcid.index(ls[1][j - 1])
        left = dp_score[i][j-1] + gapScore
        upperLeft = dp_score[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        up = dp_score[i][j-1] + gapScore
        if (max(upperLeft, left, up) == upperLeft):
            dp_score[i][j] = upperLeft
            dp_direction[i][j] = 2
        elif (max(upperLeft, left, up) == left):
            dp_score[i][j] = left
            dp_direction[i][j] = 1
        elif (max(upperLeft, left, up) == up):
            dp_score[i][j] = up
            dp_direction[i][j] = 3
            
dp_score_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))
dp_direction_local = np.zeros((len(ls[0]) + 1, len(ls[1]) + 1))

direction = ['Left', 'Upper Left', 'Up']

dp_score_local[0][0] = 1
for i in range(1, len(ls[0]) + 1):
    dp_score_local[i][0] = dp_score_local[i-1][0] + gapScore
    if dp_score_local[i][0] < 0:
        dp_score_local[i][0] = 0
for j in range(1, len(ls[1]) + 1):
    dp_score_local[0][j] = dp_score_local[0][j-1] + gapScore
    if dp_score_local[0][j] < 0:
        dp_score_local[0][j] = 0

for i in range(1, len(ls[0]) + 1):
    for j in range(1, len(ls[1]) + 1):
        aminoIndex1 = aminoAcid.index(ls[0][i - 1])
        aminoIndex2 = aminoAcid.index(ls[1][j - 1])
        left = dp_score_local[i][j-1] + gapScore
        upperLeft = dp_score_local[i-1][j-1] + pamx[aminoIndex1][aminoIndex2]
        up = dp_score_local[i][j-1] + gapScore
        if (max(upperLeft, left, up) < 0):
            dp_score_local[i][j] = 0
        elif (max(upperLeft, left, up) == left):
            dp_score_local[i][j] = left
            dp_direction_local[i][j] = 1
        elif (max(upperLeft, left, up) == upperLeft):
            dp_score_local[i][j] = upperLeft
            dp_direction_local[i][j] = 2
        elif (max(upperLeft, left, up) == up):
            dp_score_local[i][j] = up
            dp_direction_local[i][j] = 3

#for i in range(len(ls[0]) + 1):
#    for j in range(len(ls[1]) + 1):
#        print(i, j, ":", dp_score_local[i][j], end=" ")
#    print('\n')
            
#Choose global or local alignment and trace back, then save the result to fasta

if args.aln == 'global':
    aminoAcid1 = []
    aminoAcid2 = []
    row = len(ls[0])
    column = len(ls[1])

    traceback = dp_direction[row][column]

    while (traceback != 0):
        if traceback == 3:
            aminoAcid1.append(ls[0][row-1])
            aminoAcid2.append('-')
            row -= 1
        elif traceback == 2:
            aminoAcid1.append(ls[0][row-1])
            aminoAcid2.append(ls[1][column-1])
            row -= 1
            column -= 1
        elif traceback == 1:
            aminoAcid1.append('-')
            aminoAcid2.append(ls[1][column-1])
            column -= 1
        traceback = dp_direction[row][column]

    aminoAcid1.reverse()
    aminoAcid2.reverse()
    result = [aminoAcid1, aminoAcid2]

    file = open(args.output, 'w')
    lines = []
    for i in range(len(result)):
        lines.append(name[i])
        a = ''.join(result[i]) + "\n"
        lines.append(a)
    file.writelines(lines)
    file.close()


elif (args.aln == 'local'):
    maxScore = 0
    for i in range(len(ls[0]) + 1):
        for j in range(len(ls[1]) + 1):
            if dp_score_local[i][j] > maxScore:
                maxScore = dp_score_local[i][j]
    maxIndex = []
    for i in range(len(ls[0]) + 1):
        for j in range(len(ls[1]) + 1):
            if dp_score_local[i][j] == maxScore:
                maxIndex.append((i, j))
    #print(maxIndex)
    result = []
    for i in range(len(maxIndex)):
        aminoAcid1 = []
        aminoAcid2 = []
        (row, column) = maxIndex[i]

        traceback = dp_direction_local[row][column]

        while (traceback != 0):
            if traceback == 3:
                aminoAcid1.append(ls[0][row-1])
                aminoAcid2.append('-')
                row -= 1
            elif traceback == 2:
                aminoAcid1.append(ls[0][row-1])
                aminoAcid2.append(ls[1][column-1])
                row -= 1
                column -= 1
            elif traceback == 1:
                aminoAcid1.append('-')
                aminoAcid2.append(ls[1][column-1])
                column -= 1
            traceback = dp_direction_local[row][column]
        aminoAcid1.reverse()
        aminoAcid2.reverse()
        result.append(aminoAcid1)
        result.append(aminoAcid2)
    
    finalResult = []
    maxLen = 0
    for i in range(len(result)):
        if len(result[i]) > maxLen:
            finalResult = []
            finalResult.append(result[i])
            maxLen = len(result[i])
        elif len(result[i]) == maxLen:
            finalResult.append(result[i])
    
    file = open(args.output, 'w')
    lines = []
    for i in range(0, len(finalResult), 2):
        lines.append(name[0])
        a = ''.join(finalResult[i]) + "\n"
        lines.append(a)
        lines.append(name[1])
        b = ''.join(finalResult[i+1]) + "\n"
        lines.append(b)
    file.writelines(lines)
    file.close()