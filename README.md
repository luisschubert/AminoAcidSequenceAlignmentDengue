# AminoAcidSequenceAlignmentDengue
Program for aligning Amino Acid Sequences using Needleman-Wunsch Algorithm

https://docs.google.com/presentation/d/1yMQ2yGutQMkOrLSY-CqGH5ZAxHhQDcnZM-GQ-c66BL4/edit?usp=sharing

### alignment of 2 amino acid sequences
```
# instantiate Needleman-Wunsch with blosum 62 and a gap penalty of -10
NeedlemanWunsch = NeedlemanWunsch(‘blosum62’, -10)

# declare 2 sequences
sequence1 = “KVVGLYGNGVVTNG”
sequence2 = “DKKGKVVGLYGNGVVTR”

# align sequence1 with sequence2
NeedlemanWunsch.align(sequence1, sequence2)
```

### Ouput from alignment
```

Value Matrix

 1.0  0.5   0.0  -0.5  -1.0  -1.5  -2.0  -2.5  -3.0  -3.5  -4.0  -4.5  -5.0  -5.5  -6.0 
 0.5  0.0  -2.5  -3.0  -1.5  -5.0  -4.5  -3.0  -1.5  -4.0  -6.5  -7.0  -5.5  -4.0  -6.5 
 0.0  5.5  -2.0  -4.5  -5.0  -3.5  -7.0  -6.5  -3.0  -3.5  -6.0  -8.5  -8.0  -5.5  -6.0 
-0.5  5.0   3.5  -4.0  -6.5  -7.0  -5.5  -9.0  -6.5  -5.0  -5.5  -8.0  -9.5  -8.0  -7.5 
-1.0 -2.5   2.0   0.5   2.0  -8.0 -10.0   0.5  -9.0  -0.5  -8.0  -8.5 -10.0  -9.5  -2.0 
-1.5  4.0  -4.5   0.0  -1.5   0.0 -10.0  -9.5   0.5  -9.5  -2.5 -10.0  -9.5 -10.0 -11.5 
-2.0 -3.5   8.0  -0.5  -3.0  -0.5  -1.0 -11.0  -9.5  -2.5  -5.5   1.5  -8.5 -12.5 -13.0 
-2.5 -4.0   0.5  12.0   2.0  -2.0  -1.5  -4.0 -14.0 -12.5   1.5  -1.5   1.5  -8.5 -15.5 
-3.0 -4.5  -7.0   2.0  18.0   8.0  -2.0   4.5  -4.0  -8.0  -8.5  -1.5  -3.5   1.5  -2.5 
-3.5 -5.0  -3.5  -6.0   8.0  22.0  12.0   2.0   1.5  -8.0  -7.0  -7.5  -2.5  -6.5  -2.5 
-4.0 -5.5  -6.0  -4.5  -2.0  12.0  29.0  19.0   9.0  -1.0  -9.0  -8.0  -9.5  -4.5  -9.5 
-4.5 -6.0  -8.5  -9.0   1.5   2.0  19.0  35.0  25.0  15.0   5.0  -5.0 -10.0  -9.5   1.5 
-5.0 -4.5  -9.0 -11.5  -8.5  -1.5   9.0  25.0  41.0  31.0  21.0  11.0   1.0  -4.0  -8.5 
-5.5 -7.0  -7.5 -12.0  -5.5 -11.5  -1.0  15.0  31.0  47.0  37.0  27.0  17.0   7.0   2.0 
-6.0 -7.5  -3.0  -3.5 -13.5  -4.5 -11.0   5.0  21.0  37.0  51.0  41.0  31.0  21.0  11.0 
-6.5 -8.0  -3.5   1.0  -6.5 -12.5  -5.5  -5.0  11.0  27.0  41.0  55.0  45.0  35.0  25.0 
-7.0 -7.5  -8.0  -3.5  -1.0  -7.5 -14.5  -7.5   1.0  17.0  31.0  45.0  60.0  50.0  40.0 
-7.5 -5.0 -10.5 -11.0  -5.5  -3.0  -9.5 -16.5  -7.5   7.0  21.0  35.0  50.0  60.0  50.0


Direction Matrix:

00  ['nw'] 01 ['w']   02 ['w']   03 ['w']   04 ['w']   05 ['w']   06 ['w']       07 ['w']        08 ['w']       09 ['w']        010 ['w']        011 ['w']        012 ['w']   013 ['w']   014 ['w'] 
10  ['n']  11 ['nw']  12 ['nw']  13 ['nw']  14 ['nw']  15 ['nw']  16 ['nw']      17 ['nw']       18 ['nw']      19 ['nw']       110 ['n']        111 ['nw']       112 ['nw']  113 ['nw']  114 ['nw'] 
20  ['n']  21 ['nw']  22 ['nw']  23 ['nw']  24 ['nw']  25 ['nw']  26 ['nw']      27 ['nw']       28 ['nw']      29 ['nw']       210 ['nw']       211 ['nw']       212 ['nw']  213 ['nw']  214 ['nw'] 
30  ['n']  31 ['nw']  32 ['nw']  33 ['nw']  34 ['nw']  35 ['nw']  36 ['nw']      37 ['nw']       38 ['nw']      39 ['nw']       310 ['nw']       311 ['nw']       312 ['nw']  313 ['nw']  314 ['nw'] 
40  ['n']  41 ['nw']  42 ['nw']  43 ['nw']  44 ['nw']  45 ['w']   46 ['nw']      47 ['nw']       48 ['nw']      49 ['nw']       410 ['nw']       411 ['nw']       412 ['nw']  413 ['nw']  414 ['nw'] 
50  ['n']  51 ['nw']  52 ['nw']  53 ['nw']  54 ['nw']  55 ['nw']  56 ['nw', 'w'] 57 ['n']        58 ['nw']      59 ['w']        510 ['nw']       511 ['nw']       512 ['nw']  513 ['nw']  514 ['nw'] 
60  ['n']  61 ['nw']  62 ['nw']  63 ['nw']  64 ['nw']  65 ['nw']  66 ['nw']      67 ['w']        68 ['n']       69 ['nw']       610 ['nw']       611 ['nw']       612 ['w']   613 ['nw']  614 ['nw'] 
70  ['n']  71 ['nw']  72 ['nw']  73 ['nw']  74 ['w']   75 ['nw']  76 ['nw']      77 ['nw']       78 ['nw', 'w'] 79 ['nw', 'n']  710 ['nw']       711 ['nw']       712 ['nw']  713 ['w']   714 ['nw'] 
80  ['n']  81 ['nw']  82 ['nw']  83 ['n']   84 ['nw']  85 ['w']   86 ['w']       87 ['nw']       88 ['nw']      89 ['nw']       810 ['n']        811 ['nw']       812 ['nw']  813 ['nw']  814 ['nw'] 
90  ['n']  91 ['nw']  92 ['nw']  93 ['nw']  94 ['n']   95 ['nw']  96 ['w']       97 ['w']        98 ['nw']      99 ['nw']       910 ['nw']       911 ['nw']       912 ['nw']  913 ['nw']  914 ['nw'] 
100 ['n'] 101 ['nw'] 102 ['nw'] 103 ['nw'] 104 ['n']  105 ['n']  106 ['nw']     107 ['w']       108 ['w']      109 ['w']       1010 ['nw']      1011 ['nw']      1012 ['nw'] 1013 ['nw'] 1014 ['nw'] 
110 ['n'] 111 ['nw'] 112 ['nw'] 113 ['nw'] 114 ['nw'] 115 ['n']  116 ['n']      117 ['nw']      118 ['w']      119 ['nw', 'w'] 1110 ['w']       1111 ['w']       1112 ['nw'] 1113 ['nw'] 1114 ['nw'] 
120 ['n'] 121 ['nw'] 122 ['nw'] 123 ['nw'] 124 ['n']  125 ['nw'] 126 ['n']      127 ['n']       128 ['nw']     129 ['w']       1210 ['w']       1211 ['w']       1212 ['w']  1213 ['nw'] 1214 ['n'] 
130 ['n'] 131 ['nw'] 132 ['nw'] 133 ['nw'] 134 ['nw'] 135 ['n']  136 ['n']      137 ['nw', 'n'] 138 ['n']      139 ['nw']      1310 ['w']       1311 ['w']       1312 ['w']  1313 ['w']  1314 ['nw'] 
140 ['n'] 141 ['nw'] 142 ['nw'] 143 ['nw'] 144 ['w']  145 ['nw'] 146 ['n']      147 ['n']       148 ['n']      149 ['n']       1410 ['nw']      1411 ['nw', 'w'] 1412 ['w']  1413 ['w']  1414 ['w'] 
150 ['n'] 151 ['nw'] 152 ['nw'] 153 ['nw'] 154 ['nw'] 155 ['nw'] 156 ['nw']     157 ['n']       158 ['n']      159 ['n']       1510 ['nw', 'n'] 1511 ['nw']      1512 ['w']  1513 ['w']  1514 ['w'] 
160 ['n'] 161 ['nw'] 162 ['nw'] 163 ['nw'] 164 ['nw'] 165 ['nw'] 166 ['nw']     167 ['nw']      168 ['n']      169 ['n']       1610 ['n']       1611 ['n']       1612 ['nw'] 1613 ['w']  1614 ['w'] 
170 ['n'] 171 ['nw'] 172 ['nw'] 173 ['nw'] 174 ['nw'] 175 ['nw'] 176 ['nw']     177 ['nw']      178 ['nw']     179 ['n']       1710 ['n']       1711 ['n']       1712 ['n']  1713 ['nw'] 1714 ['w'] 


Alignment:

 ----KVVGLYGNGVVTNG
 DKKGKVVGLYGNGVVTR-
score: 50.0
```








### alignment of Dengue 3 & 4 NS3 protein fragments
```
# instantiate Needleman-Wunsch with blosum 62 and a gap penalty of -10
NeedlemanWunsch = NeedlemanWunsch(‘blosum62’, -10)

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

print "Dengue 3 & 4 Virus NS3 Pairwise Sequence Alignment"
NeedlemanWunsch.align(dengue3_NS3, dengue4_NS3)
```
