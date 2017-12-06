import numpy

class InvalidPairException(Exception):
  pass

class NonExistantDistanceMatrix(Exception):
    pass

'''
Distance Matrix Helper Class
https://gist.github.com/jwintersinger/1870047
'''
class DistanceMatrix:
  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)

  def _load_matrix(self, matrix_filename):
    matrix_filename = 'DistanceMatrices/' + matrix_filename
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


'''
Needleman-Wunsch Alignment
'''
class NeedlemanWunsch:
    def __init__(self, distanceMatrixName, gapPenalty):
        availableDistanceMatrices = ['blosum62', 'blosum80']
        #we should be initializing this with some parameters...

        # gap penalty
        self.gapPenalty = gapPenalty

        # type of distance matrix used for alignment
        if distanceMatrixName in availableDistanceMatrices:
            self.distanceMatrix = DistanceMatrix(distanceMatrixName)
        else:
            raise NonExistantDistanceMatrix('%s' % (distanceMatrixName))
        # extension penalty (?)
        pass

    def align(self, s1, s2):
        # we could be producing some better metrics here and instead of just printing it.
        # it would be cleaner to return the results as an object to the caller
        metrics = {}

        # if len(s1) > len(s2):
        #     diff = len(s1) - len(s2)
        #     for i in range(diff):
        #         s2 = ' ' + s2
        #
        # elif len(s1) < len(s2):
        #     diff = len(s2) - len(s1)
        #     for i in range(diff):
        #         s1 = ' ' + s1

        s1 = ' ' + s1
        s2 = ' ' + s2
        valueMatrix = numpy.zeros((len(s2), len(s1)))
        directionMatrix = {}

        highestValue = -999999
        highestValueIndex = (0, 0)
        # forward pass to build the matrix of values and matrix of directions of arrows.
        for row in range(len(s2)):
            for col in range(len(s1)):
                if row == 0 and col == 0:
                    # directionMatrix[row][col] = ['nw']
                    directionMatrix[str(row) + str(col)] = ['nw']
                    # int(matrix.lookup_score(s1[col], s2[row]))
                    valueMatrix[row][col] = 1
                elif row == 0:
                    # directionMatrix[row][col] = ['w']
                    directionMatrix[str(row) + str(col)] = ['w']
                    # int(matrix.lookup_score(s1[col], s2[row]))

                    # self.gapPenalty would turn into an self.extensionPenalty if used
                    # for testing
                    newValue = valueMatrix[row][col - 1] + (-0.5)

                    # newValue = valueMatrix[row][col - 1] + self.gapPenalty
                    if (newValue > highestValue):
                        highestValue = newValue
                        highestValueIndex = (row, col)
                    valueMatrix[row][col] = newValue
                elif col == 0:
                    # directionMatrix[row][col] = ['n']
                    directionMatrix[str(row) + str(col)] = ['n']
                    # int(matrix.lookup_score(s1[col], s2[row]))

                    # self.gapPenalty would turn into an self.extensionPenalty if used
                    # for testing
                    newValue = valueMatrix[row - 1][col] + (-0.5)

                    # newValue = valueMatrix[row - 1][col] + self.gapPenalty
                    if (newValue > highestValue):
                        highestValue = newValue
                        highestValueIndex = (row, col)
                    valueMatrix[row][col] = newValue
                else:
                    map = {}
                    map['nw'] = valueMatrix[row - 1][col - 1]
                    map['w'] = valueMatrix[row][col - 1]
                    map['n'] = valueMatrix[row - 1][col]
                    matchScore = int(self.distanceMatrix.lookup_score(s1[col], s2[row]))
                    maxValue, maxList = self.findMax(map, matchScore)
                    if (maxValue > highestValue):
                        highestValue = maxValue
                        highestValueIndex = (row, col)
                    if 'nw' in maxList:
                        # valueMatrix[row - 1][col - 1] + int(matrix.lookup_score(s1[col], s2[row]))
                        valueMatrix[row][col] = maxValue
                        directionMatrix[str(row) + str(col)] = maxList
                    elif 'w' in maxList:
                        # int(matrix.lookup_score(s1[col], s2[row]))
                        # valueMatrix[row][col - 1] + gapPenalty
                        valueMatrix[row][col] = maxValue
                        directionMatrix[str(row) + str(col)] = maxList
                    elif 'n' in maxList:
                        # int(matrix.lookup_score(s1[col], s2[row]))
                        # valueMatrix[row - 1][col] + gapPenalty
                        valueMatrix[row][col] = maxValue
                        directionMatrix[str(row) + str(col)] = maxList

        # print 'value matrix'
        # self.printMatrix(valueMatrix)
        #
        # print 'direction matrix'
        # self.printDirection(directionMatrix, len(s2), len(s1))

        # backward pass to find the best global alignment.
        newS1 = ""
        newS2 = ""
        pathDirection = []
        foundAlignment = False
        row = len(s2) - 1
        col = len(s1) - 1
        matchScore = valueMatrix[row][col]
        lastPick = 'w'
        while not foundAlignment:
            position = str(row) + str(col)
            pathDirection.append(position)
            if position == "00":
                foundAlignment = True
            if 'nw' in directionMatrix[position]:
                newDir = 'nw'
            elif len(directionMatrix[position]) == 2:
                if lastPick == 'w':
                    newDir = 'n'
                    lastPick = 'n'
                else:
                    newDir = 'w'
                    lastPick = 'w'
            else:
                newDir = directionMatrix[position][0]
            if newDir == 'w':
                newS1 = s1[col] + newS1
                newS2 = "-" + newS2
                col = col - 1
            elif newDir == 'n':
                newS1 = "-" + newS1
                newS2 = s2[row] + newS2
                row = row - 1
            elif newDir == 'nw':
                newS1 = s1[col] + newS1
                newS2 = s2[row] + newS2
                col = col - 1
                row = row - 1
        # print pathDirection
        metrics['pathDirection'] = pathDirection
        metrics['newS1'] = newS1
        metrics['newS2'] = newS2
        metrics['score'] = matchScore
        self.printValueMatrix(valueMatrix)
        self.printDirectionMatrix(directionMatrix, len(s2), len(s1))
        print newS1
        print newS2
        print "score: " + str(matchScore)
        print "\n"

    def printValueMatrix(self, m):
        for i in range(len(m)):
            for j in range(len(m[0])):
                print m[i][j],
            print '\n'

    def printDirectionMatrix(self, m, rows, cols):
        for row in range(rows):
            for col in range(cols):
                print str(row) + str(col) + " " + str(m[str(row) + str(col)]),
            print "\n"

    def findMax(self, map, matchScore):
        map['nw'] = map['nw'] + matchScore
        map['n'] = map['n'] + self.gapPenalty
        map['w'] = map['w'] + self.gapPenalty
        # you are adding because gapPenalty is a negative value

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

def getSequence(id):
    sequencePath = 'DengueVirusFullGenome/AminoAcidSequences/DengueV' + str(id) + 'AminoAcidSequence'
    with open(sequencePath, 'r') as sequenceFile:
        data = sequenceFile.read().replace('\n', '')
    return data


dengue1_fullGenome = getSequence(1)
dengue2_fullGenome = getSequence(2)
dengue3_fullGenome = getSequence(3)
dengue4_fullGenome = getSequence(4)

#
NeedlemanWunsch = NeedlemanWunsch('blosum62', -10)
# dengue1_NS3 = "SGVLWDTPSPPEVERAVLDNGIYRILQRGLLGRSQVGVGVFQEGVFHTMWHVTRGAVLMYQGKRLEPSWASVKKDLISYGGGWRFQGSWNTGEEVQVIAVEPGKNPKNVQTTPGTFKTPEGEIGAIALDFKPGTSGSPIVNREGKIVGLYGNGVVTTSGTYVSAIAQAKVSQEGPLPEIEDEVFRK"

# NeedlemanWunsch.align(s1, dengue1_NS3)

'''
http://www.uniprot.org/uniprot/A0A0K0TQ28
>tr|A0A0K0TQ28|A0A0K0TQ28_9FLAV NS3 protease (Fragment) OS=Dengue virus 1 PE=4 SV=1
SGVLWDTPSPPEVERAVLDNGIYRILQRGLLGRSQVGVGVFQEGVFHTMWHVTRGAVLMY
QGKRLEPSWASVKKDLISYGGGWRFQGSWNTGEEVQVIAVEPGKNPKNVQTTPGTFKTPE
GEIGAIALDFKPGTSGSPIVNREGKIVGLYGNGVVTTSGTYVSAIAQAKVSQEGPLPEIE
DEVFRK
'''

dengue1_NS3 = "SGVLWDTPSPPEVERAVLDNGIYRILQRGLLGRSQVGVGVFQEGVFHTMWHVTRGAVLMYQGKRLEPSWASVKKDLISYGGGWRFQGSWNTGEEVQVIAVEPGKNPKNVQTTPGTFKTPEGEIGAIALDFKPGTSGSPIVNREGKIVGLYGNGVVTTSGTYVSAIAQAKVSQEGPLPEIEDEVFRK"

'''
https://www.ncbi.nlm.nih.gov/protein/AAL07258.1
DKKGKVVGLYGNGVVTRSGTYVSAIAQTEKSIEDNPEIEDDIFRKKRLTIMD
LHPGAGKTKRYLPAIVREAIKRGLRTLILAPTRVVAAEMEEALRGLPIRYQ
'''
dengue2_NS3 = "DKKGKVVGLYGNGVVTRSGTYVSAIAQTEKSIEDNPEIEDDIFRKKRLTIMDLHPGAGKTKRYLPAIVREAIKRGLRTLILAPTRVVAAEMEEALRGLPIRYQ"


'''
http://www.uniprot.org/uniprot/Q90213
>tr|Q90213|Q90213_9FLAV Nonstructural protein NS3 (Fragment) OS=Dengue virus 3 PE=4 SV=1
KVVGLYGNGVVTKNGGYVSGIAQTNAEPDGPTPELEEEMFKKRNLTIMDLHPGSGKTRKY
LPAIVREAIKRRLRTLILAPTRVVAAEMEEALKGLPIRYQTTATKSEHTGREIVDLMCHA
TFTMRLLSPVRV
'''

dengue3_NS3 = 'KVVGLYGNGVVTKNGGYVSGIAQTNAEPDGPTPELEEEMFKKRNLTIMDLHPGSGKTRKYLPAIVREAIKRRLRTLILAPTRVVAAEMEEALKGLPIRYQTTATKSEHTGREIVDLMCHATFTMRLLSPVRV'


'''
http://www.uniprot.org/uniprot/Q66422
>tr|Q66422|Q66422_9FLAV Nonstructural protein NS3 (Fragment) OS=Dengue virus 4 PE=4 SV=1
RKGKVIGLYGNGVVTKSGDYVSAITQAERTGEPDYEVDEDIFRKKRLTIMDLHPGAGKTK
RILPSIVREALKRRLRTLILAPTRVVAAEMEEALRGLPIRYQTPAVKSEHTGREIVDLMC
HATFTTRLLSSTRVPNYNLI
'''

dengue4_NS3 = 'RKGKVIGLYGNGVVTKSGDYVSAITQAERTGEPDYEVDEDIFRKKRLTIMDLHPGAGKTKRILPSIVREALKRRLRTLILAPTRVVAAEMEEALRGLPIRYQTPAVKSEHTGREIVDLMCHATFTTRLLSSTRVPNYNLI'

#
# NeedlemanWunsch = NeedlemanWunsch()
# NeedlemanWunsch.align(dengue1_NS3, dengue2_NS3)

# NeedlemanWunsch = NeedlemanWunsch()
# NeedlemanWunsch.align("DKKKKKKKKKKKGW", "DDKKWW")
#
# NeedlemanWunsch.align(s1, dengue1_NS3)
# NeedlemanWunsch.align(s2, dengue2_NS3)


# NeedlemanWunsch.align(dengue1_NS3, dengue2_NS3)

# print "Dengue 1 & 2 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue1_fullGenome, dengue2_fullGenome)
# print "Dengue 1 & 3 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue1_fullGenome, dengue3_fullGenome)
# print "Dengue 1 & 4 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue1_fullGenome, dengue4_fullGenome)
# print "Dengue 2 & 3 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue2_fullGenome, dengue3_fullGenome)
# print "Dengue 2 & 4 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue2_fullGenome, dengue4_fullGenome)
# print "Dengue 3 & 4 Virus Full Genome Pairwise Sequence Alignment"
# NeedlemanWunsch.align(dengue3_fullGenome, dengue4_fullGenome)

print "Dengue 1 & 2 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue1_NS3, dengue2_NS3)
print "Dengue 1 & 3 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue1_NS3, dengue3_NS3)
print "Dengue 1 & 4 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue1_NS3, dengue4_NS3)
print "Dengue 2 & 3 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue2_NS3, dengue3_NS3)
print "Dengue 2 & 4 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue2_NS3, dengue4_NS3)
print "Dengue 3 & 4 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue3_NS3, dengue4_NS3)