import numpy
import time
# import textfile for protein sequences of dengue 1 - 4 sequences.
# import the scoring matrix.
# define the gap penalty

# compute a distance matrix

# create a guide tree

# align sequences

class InvalidPairException(Exception):
  pass

class Matrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0)
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix

  def lookup_score(self, a, b):
    a = a.upper()
    b = b.upper()
    if a == ' ':
        a = '*'
    if b == ' ':
        b = '*'

    if a not in self._matrix or b not in self._matrix[a]:
      raise InvalidPairException('[%s, %s]' % (a, b))
    return self._matrix[a][b]

def getSequence(id):
    sequencePath = 'DengueV' + str(id) + 'AminoAcidSequences'
    with open(sequencePath, 'r') as sequenceFile:
        data = sequenceFile.read().replace('\n', '')
    return data

def pairwiseAlignment(s1, s2):
    matrix_filename = 'blosum62'
    gapPenalty  = -2
    matrix = Matrix(matrix_filename)
    if len(s1) > len(s2):
        diff = len(s1) - len(s2)
        for i in range(diff):
            s2 = ' ' + s2

    elif len(s1) < len(s2):
        diff = len(s2) - len(s1)
        for i in range(diff):
            s1 = ' ' + s1

    s1 = ' ' + s1
    s2 = ' ' + s2
    valueMatrix = numpy.zeros((len(s1), len(s1)))
    directionMatrix = {}

    # forward pass to build the matrix of values and matrix of directions of arrows.
    for row in range(len(s1)):
        for col in range(len(s1)):
            if row == 0 and col == 0:
                # directionMatrix[row][col] = ['nw']
                directionMatrix[str(row) + str(col)] = ['nw']
                valueMatrix[row][col] = int(matrix.lookup_score(s1[col], s2[row]))
            elif row == 0:
                # directionMatrix[row][col] = ['w']
                directionMatrix[str(row) + str(col)] = ['w']
                # int(matrix.lookup_score(s1[col], s2[row]))
                valueMatrix[row][col] = valueMatrix[row][col - 1] + gapPenalty
            elif col == 0:
                # directionMatrix[row][col] = ['n']
                directionMatrix[str(row) + str(col)] = ['n']
                # int(matrix.lookup_score(s1[col], s2[row]))
                valueMatrix[row][col] = valueMatrix[row - 1][col] + gapPenalty
            else:
                map = {}
                map['nw'] = valueMatrix[row - 1][col - 1]
                map['w'] = valueMatrix[row][col - 1]
                map['n'] = valueMatrix[row - 1][col]
                maxValue, maxList = findMax(map)
                if 'nw' in maxList:
                    valueMatrix[row][col] = valueMatrix[row - 1][col - 1] + int(matrix.lookup_score(s1[col], s2[row]))
                    directionMatrix[str(row) + str(col)] = maxList
                elif 'w' in maxList:
                    # int(matrix.lookup_score(s1[col], s2[row]))
                    valueMatrix[row][col] = valueMatrix[row][col - 1] + gapPenalty
                    directionMatrix[str(row) + str(col)] = maxList
                elif 'n' in maxList:
                    # int(matrix.lookup_score(s1[col], s2[row]))
                    valueMatrix[row][col] = valueMatrix[row - 1][col] + gapPenalty
                    directionMatrix[str(row) + str(col)] = maxList

    # print 'value matrix'
    # printMatrix(valueMatrix)
    #
    # print 'direction matrix'
    # printDirection(directionMatrix, len(s1))

    #backward pass to find the best alignment.
    pass


def printMatrix(m):
    for i in range(len(m)):
        for j in range(len(m[0])):
            print m[i][j],
        print '\n'

def printDirection(m, size):
    for row in range(size):
        for col in range(size):
            print str(row) + str(col)+" " + str(m[str(row) + str(col)]),
        print "\n"

def findMax(map):
    maxList = ['nw']
    maxValue = map['nw']
    if map['n'] == maxValue:
        maxList.append('n')
    elif map['n'] > maxValue:
        maxValue = map['n']
        maxList = ['n']
    if map['w'] == maxValue:
        maxList.append('w')
    elif map['w'] > maxValue:
        maxValue = map['w']
        maxList = ['w']
    return maxValue, maxList

def findMin(map):
    min = ['nw']
    minValue = map['nw']
    if map['n'] == minValue:
        min.append('n')
    elif map['n'] < minValue:
        minValue = map['n']
        min = ['n']
    if map['w'] == minValue:
        min.append('w')
    elif map['w'] < minValue:
        minValue = map['w']
        min = ['w']
    return minValue, min

def constructDistanceMatrix(listOfSequences):
    distanceMatrix = numpy.zeros((4,4))
    for i in range(len(listOfSequences)):
        for j in range(i,len(listOfSequences)):
            if j > i:
                s1 = listOfSequences[i]
                s2 = listOfSequences[j]
                distanceMatrix[j][i] = pairwiseAlignment(s1, s2)

listOfSequences = []
s1 = getSequence(1)
s2 = getSequence(2)
s3 = getSequence(3)
s4 = getSequence(4)
listOfSequences.append(s1)
listOfSequences.append(s2)
listOfSequences.append(s3)
listOfSequences.append(s4)

print s1 + "\n"
print s2 + "\n"
print s3 + "\n"
print s4 + "\n"

start = time.time()
pairwiseAlignment(s1[:10], s2[:10])
end = time.time()
print "first 10 " + str(end-start)

start = time.time()
pairwiseAlignment(s1[:100], s2[:100])
end = time.time()
print "first 100 " + str(end-start)

start = time.time()
pairwiseAlignment(s1[:1000], s2[:1000])
end = time.time()
print "first 1000 " + str(end-start)

start = time.time()
pairwiseAlignment(s1, s2)
end = time.time()
print "full sequence " + str(end-start)

# pairwiseAlignment("www", "www")