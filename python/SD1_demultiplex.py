# Left side barcodes
lbc = {"AACCA":"V6L-0",
"CCAAC":"V6L-1",
"GGTTG":"V6L-2",
"TTGGT":"V6L-3",
"AGGTG":"V6L-4",
"CTTAC":"V6L-5",
"GAACT":"V6L-6",
"TCCGA":"V6L-7",
"ACAGT":"V6L-8",
"CACTG":"V6L-9",
"GTGAC":"V6L-10",
"TGTCA":"V6L-11",
"AAACC":"V6L-12",
"CCCAA":"V6L-13",
"GGGTT":"V6L-14",
"TTTAA":"V6L-15",
"GAGGT":"V6L-16",
"CCTTA":"V6L-17",
"TGAAC":"V6L-18",}

# Right side barcodes
rbc = {"TGGTT":"V6R-0",
"GTTGG":"V6R-1",
"CAACC":"V6R-2",
"ACCAA":"V6R-3",
"CACCT":"V6R-4",
"GTAAG":"V6R-5",
"AGTTC":"V6R-6",
"TCGGA":"V6R-7",
"ACTGT":"V6R-8",
"CAGTG":"V6R-9",
"GTCAC":"V6R-10",
"TGACA":"V6R-11",
"GGTTT":"V6R-12",
"TTGGG":"V6R-13",
"AACCC":"V6R-14",
"TTAAA":"V6R-15",
"ACCTC":"V6R-16",
"TAAGG":"V6R-17",
"GTTCA":"V6R-18",}

samples = {
"V6L-0V6R-0":"SD1.CK4D1.CK4.B",
"V6L-0V6R-1":"SD1.CK4D1.CK4.G",
"V6L-0V6R-10":"SD1.CK4D2.CK4.B",
"V6L-0V6R-11":"SD1.CK4D2.CK4.G",
"V6L-0V6R-12":"SD1.CK4D2.CK4.O",
"V6L-0V6R-13":"SD1.CK4D2.CK4.P",
"V6L-0V6R-14":"SD1.CK4D2.CK4.W",
"V6L-0V6R-17":"SD1.CK4D2.CK14.O",
"V6L-0V6R-18":"SD1.CK4D2.CK14.P",
"V6L-0V6R-2":"SD1.CK4D1.CK4.O",
"V6L-0V6R-3":"SD1.CK4D1.CK4.P",
"V6L-0V6R-4":"SD1.CK4D1.CK4.W",
"V6L-0V6R-5":"H6",
"V6L-0V6R-6":"SD1.CK4D1.CK14.G",
"V6L-0V6R-7":"SD1.CK4D1.CK14.O",
"V6L-0V6R-8":"SD1.CK4D1.CK14.P",
"V6L-10V6R-0":"SD1.CK4D2.CK4.G.T10",
"V6L-10V6R-1":"SD1.CK4D2.CK4.O.T10",
"V6L-10V6R-10":"SD1.CK4.D3.CK4.B.T10",
"V6L-10V6R-11":"SD1.CK4.D3.CK4.G.T10",
"V6L-10V6R-12":"SD1.CK4.D3.CK4.O.T10",
"V6L-10V6R-13":"SD1.CK4.D3.CK4.P.T10",
"V6L-10V6R-14":"SD1.CK4.D3.CK4.W.T10",
"V6L-10V6R-15":"SD1.CK4.D3.CK14.B.T10",
"V6L-10V6R-16":"A1.+.D2.GA",
"V6L-10V6R-17":"SD1.CK4.D3.CK14.O.T10",
"V6L-10V6R-18":"SD1.CK4.D3.CK14.W.T10",
"V6L-10V6R-2":"SD1.CK4D2.CK4.P.T10",
"V6L-10V6R-3":"SD1.CK4D2.CK4.W.T10",
"V6L-10V6R-4":"SD1.CK4D2.CK4.W.T10",
"V6L-10V6R-5":"SD1.CK4.D2.CK14.B.T10",
"V6L-10V6R-6":"SD1.CK4.D2.CK14.G.T10",
"V6L-10V6R-7":"SD1.CK4.D2.CK14.O.T10",
"V6L-10V6R-8":"SD1.CK4.D2.CK14.P.T10",
"V6L-10V6R-9":"SD1.CK4.D2.CK14.W.T10",
"V6L-11V6R-0":"SD1.CK4.D3.CK14.P.T10",
"V6L-11V6R-1":"SD1.CK14.D1.CK4.B.T10",
"V6L-11V6R-10":"SD1.CK14.D1.CK14.W.T10",
"V6L-11V6R-11":"SD1.CK14.D2.CK4.B.T10",
"V6L-11V6R-12":"SD1.CK14.D2.CK4.G.T10",
"V6L-11V6R-13":"SD1.CK14.D2.CK4.O.T10",
"V6L-11V6R-14":"SD1.CK14.D2.CK4.P.T10",
"V6L-11V6R-15":"SD1.CK4.D2.CK4.W.T10",
"V6L-11V6R-16":"SD1.CK14.D2.CK14.B.T10",
"V6L-11V6R-17":"SD1.CK14.D2.CK14.G.T10",
"V6L-11V6R-18":"A1.+.H3.GA",
"V6L-11V6R-2":"SD1.CK14.D1.CK4.G.T10",
"V6L-11V6R-3":"SD1.CK14.D1.CK4.O.T10",
"V6L-11V6R-4":"SD1.CK14.D1.CK4.P.T10",
"V6L-11V6R-5":"SD1.CK14.D1.CK4.W.T10",
"V6L-11V6R-6":"SD1.CK14.D1.CK14.B.T10",
"V6L-11V6R-7":"SD1.CK14.D1.CK14.G.T10",
"V6L-11V6R-8":"SD1.CK14.D1.CK14.O.T10",
"V6L-11V6R-9":"SD1.CK14.D1.CK14.P.T10",
"V6L-12V6R-0":"SD1.CK14.D2.CK14.P.T10",
"V6L-12V6R-1":"SD1.CK4.D2.CK14.W.T10",
"V6L-12V6R-10":"SD1.CK14.D3.CK14.P.T10",
"V6L-12V6R-11":"SD1.CK4.D3.CK14.W.T10",
"V6L-12V6R-12":"SD1.CK4.H3.CK4.B.T22.H",
"V6L-12V6R-13":"SD1.CK4.H3.CK4.G.T22.H",
"V6L-12V6R-14":"SD1.CK4.H3.CK4.O.T22.H",
"V6L-12V6R-15":"SD1.CK4.H3.CK4.W.T22.H",
"V6L-12V6R-16":"SD1.CK4.H3.CK14.B.T22.H",
"V6L-12V6R-17":"SD1.CK4.H3.CK14.W.T22.H",
"V6L-12V6R-18":"SD1.CK4.H2.CK4.P.T34.H",
"V6L-12V6R-2":"SD1.CK14.D3.CK4.B.T10",
"V6L-12V6R-3":"SD1.CK14.D3.CK4.G.T10",
"V6L-12V6R-4":"SD1.CK14.D3.CK4.O.T10",
"V6L-12V6R-5":"SD1.CK14.D3.CK4.P.T10",
"V6L-12V6R-6":"SD1.CK14.D3.CK4.W.T10",
"V6L-12V6R-7":"SD1.CK14.D3.CK14.B.T10",
"V6L-12V6R-8":"A1.+.H1.BB",
"V6L-12V6R-9":"SD1.CK14.D3.CK14.O.T10",
"V6L-13V6R-0":"SD1.CK4.H2.CK4.W.T34.H",
"V6L-13V6R-1":"SD1.CK4.H2.CK14.G.T34.H",
"V6L-13V6R-10":"SD1.CK14.H1.CK14.W.T37.H",
"V6L-13V6R-11":"SD1.CK4.H2.CK14.B.T37.H",
"V6L-13V6R-12":"SD1.CK14.H3.CK14.B.T37.H",
"V6L-13V6R-13":"SD1.CK4.H1.CK4.O.T37.H",
"V6L-13V6R-14":"SD1.CK4.D3.CK4.B.T22.D",
"V6L-13V6R-15":"SD1.CK4.D3.CK4.G.T22.D",
"V6L-13V6R-16":"SD1.CK4.D3.CK4.O.T22.D",
"V6L-13V6R-17":"A1.+.D2.BA",
"V6L-13V6R-18":"SD1.CK4.D3.CK14.B.T22.D",
"V6L-13V6R-2":"SD1.CK4.H2.CK14.O.T34.H",
"V6L-13V6R-3":"SD1.CK4.H3.CK4.P.T34.H",
"V6L-13V6R-4":"SD1.CK4.H3.CK14.G.T34.H",
"V6L-13V6R-5":"SD1.CK4.H3.CK14.O.T34.H",
"V6L-13V6R-6":"SD1.CK4.H3.CK14.P.T34.H",
"V6L-13V6R-7":"SD1.CK4.H1.CK14.W.T34.H",
"V6L-13V6R-8":"SD1.CK4.H2.CK4.O.T37.H",
"V6L-13V6R-9":"SD1.CK14.H1.CK14.B.T37.H",
"V6L-14V6R-0":"SD1.CK4.D3.CK14.W.T22.D",
"V6L-14V6R-1":"SD1.CK4.D2.CK4.P.T34.D",
"V6L-14V6R-10":"SD1.CK4.D2.CK4.O.T37.D",
"V6L-14V6R-11":"SD1.CK14.D1.CK14.B.T37.D",
"V6L-14V6R-12":"SD1.CK14.D1.CK14.W.T37.D",
"V6L-14V6R-13":"SD1.CK14.D2.CK14.B.T37.D",
"V6L-14V6R-14":"SD1.CK14.D3.CK14.B.T37.D",
"V6L-14V6R-15":"SD1.CK4.D1.CK4.O.T37.D",
"V6L-14V6R-16":"SD1.CK4.H1.CK4.B.T57.H",
"V6L-14V6R-17":"SD1.CK4.H1.CK4.G.T57.H",
"V6L-14V6R-18":"SD1.CK4.H1.CK4.P.T57.H",
"V6L-14V6R-2":"SD1.CK4.D2.CK4.W.T34.D",
"V6L-14V6R-3":"SD1.CK4.D2.CK14.G.T34.D",
"V6L-14V6R-4":"SD1.CK4.D2.CK14.O.T34.D",
"V6L-14V6R-5":"SD1.CK4.D3.CK4.P.T34.D",
"V6L-14V6R-6":"SD1.CK4.D3.CK14.G.T34.D",
"V6L-14V6R-7":"SD1.CK4.D3.CK14.O.T34.D",
"V6L-14V6R-8":"SD1.CK4.D3.CK14.P.T34.D",
"V6L-14V6R-9":"SD1.CK4.D1.CK14.W.T34.D",
"V6L-15V6R-0":"SD1.CK4.H1.CK4.W.T57.H",
"V6L-15V6R-1":"A1.+.D3.BA",
"V6L-15V6R-10":"SD1.CK14.H1.CK4.G.T57.H",
"V6L-15V6R-11":"SD1.CK14.H1.CK4.O.T57.H",
"V6L-15V6R-12":"SD1.CK14.H1.CK4.P.T57.H",
"V6L-15V6R-13":"SD1.CK14.H1.CK4.W.T57.H",
"V6L-15V6R-14":"SD1.CK14.H1.CK14.G.T57.H",
"V6L-15V6R-15":"SD1.CK14.H1.CK14.O.T57.H",
"V6L-15V6R-16":"SD1.CK14.H1.CK14.P.T57.H",
"V6L-15V6R-17":"SD1.CK14.H3.CK4.B.T57.H",
"V6L-15V6R-18":"SD1.CK14.H3.CK4.O.T57.H",
"V6L-15V6R-2":"SD1.CK4.H1.CK14.G.T57.H",
"V6L-15V6R-3":"SD1.CK4.H1.CK14.P.T57.H",
"V6L-15V6R-4":"SD1.CK4.H1.CK4.B.T57.H",
"V6L-15V6R-5":"SD1.CK4.H2.CK4.G.T57.H",
"V6L-15V6R-6":"SD1.CK4.H2.CK14.B.T57.H",
"V6L-15V6R-7":"SD1.CK4.H2.CK14.P.T57.H",
"V6L-15V6R-8":"SD1.CK4.H2.CK14.W.T57.H",
"V6L-15V6R-9":"SD1.CK14.H1.CK4.B.T57.H",
"V6L-16V6R-0":"SD1.CK14.H3.CK4.P.T57.H",
"V6L-16V6R-1":"SD1.CK14.H3.CK4.W.T57.H",
"V6L-16V6R-10":"SD1.CK14.H3.CK4.G.T57.H",
"V6L-16V6R-11":"SD1.CK14.H3.CK14.G.T57.H",
"V6L-16V6R-12":"SD1.CK14.H3.CK14.O.T57.H",
"V6L-16V6R-13":"SD1.CK14.H3.CK14.P.T57.H",
"V6L-16V6R-14":"SD1.CK14.H3.CK14.W.T57.H",
"V6L-16V6R-15":"SD1.CK4.D1.CK4.G.T57.D",
"V6L-16V6R-16":"STR8UP.NEG.CONTROL",
"V6L-16V6R-17":"SD1.CK4.D1.CK14.G.T57.H",
"V6L-16V6R-18":"STR8UP.NEG.CONTROL",
"V6L-16V6R-2":"SD1.CK14.H2.CK4.W.T57.H",
"V6L-16V6R-3":"SD1.CK14.H2.CK14.G.T57.H",
"V6L-16V6R-4":"SD1.CK14.H2.CK14.O.T57.H",
"V6L-16V6R-5":"SD1.CK14.H2.CK14.P.T57.H",
"V6L-16V6R-6":"SD1.CK14.H2.CK4.B.T57.H",
"V6L-16V6R-7":"SD1.CK14.H2.CK4.G.T57.H",
"V6L-16V6R-8":"SD1.CK14.H2.CK4.O.T57.H",
"V6L-16V6R-9":"SD1.CK14.H2.CK4.P.T57.H",
"V6L-17V6R-0":"SD1.CK4.D2.CK4.B.T57.H",
"V6L-17V6R-1":"SD1.CK14.D1.CK4.B.T57.H",
"V6L-17V6R-10":"SD1.CK14.D2.CK4.G.T57.D",
"V6L-17V6R-11":"SD1.CK14.D2.CK4.P.T57.D",
"V6L-17V6R-12":"SD1.CK14.D2.CK4.W.T57.H",
"V6L-17V6R-13":"SD1.CK14.D2.CK14.G.T57.D",
"V6L-17V6R-14":"SD1.CK14.D2.CK14.O.T57.D",
"V6L-17V6R-15":"SD1.CK14.D2.CK14.P.T57.H",
"V6L-17V6R-16":"STR8UP.NEG.CONTROL",
"V6L-17V6R-17":"SD1.CK14.D3.CK14.B.T57.H",
"V6L-17V6R-18":"STR8UP NEG CONTROL",
"V6L-17V6R-2":"SD1.CK14.D1.CK4.G.T57.D",
"V6L-17V6R-3":"SD1.CK14.D1.CK4.O.T57.H",
"V6L-17V6R-4":"SD1.CK14.D1.CK4.P.T57.D",
"V6L-17V6R-5":"SD1.CK14.D1.CK4.W.T57.D",
"V6L-17V6R-6":"SD1.CK14.D1.CK14.G.T57.D",
"V6L-17V6R-7":"SD1.CK14.D1.CK14.O.T57.H",
"V6L-17V6R-8":"SD1.CK14.D1.CK14.P.T57.H",
"V6L-17V6R-9":"SD1.CK14.D2.CK4.B.T57.D",
"V6L-18V6R-0":"SD1.CK14.D3.CK4.O.T57.H",
"V6L-18V6R-1":"SD1.CK14.D3.CK4.P.T57.H",
"V6L-18V6R-10":"SD1.CK4.H1",
"V6L-18V6R-11":"SD1.CK4.H2",
"V6L-18V6R-12":"SD1.CK4.H3",
"V6L-18V6R-13":"SD1.CK14.D1",
"V6L-18V6R-14":"SD1.CK14.D2",
"V6L-18V6R-15":"SD1.CK14.D3",
"V6L-18V6R-16":"SD1.CK4.D1",
"V6L-18V6R-17":"SD1.CK4.D2",
"V6L-18V6R-18":"SD1.CK4.D3",
"V6L-18V6R-2":"SD1.CK14.D3.CK4.W.T57.H",
"V6L-18V6R-3":"SD1.CK14.D3.CK14.G.T57.D",
"V6L-18V6R-4":"SD1.CK14.D3.CK14.O.T57.H",
"V6L-18V6R-5":"SD1.CK14.D3.CK14.P.T57.H",
"V6L-18V6R-6":"SD1.CK14.D2.CK4.W.T57.D",
"V6L-18V6R-7":"SD1.CK14.H1",
"V6L-18V6R-8":"SD1.CK14.H2",
"V6L-18V6R-9":"SD1.CK14.H3",
"V6L-1V6R-0":"H7",
"V6L-1V6R-1":"SD1.CK4D3.CK4.B",
"V6L-1V6R-10":"H9",
"V6L-1V6R-11":"SD1.CK14D1.CK4.B",
"V6L-1V6R-12":"SD1.CK14D1.CK4.G",
"V6L-1V6R-13":"SD1.CK14D1.CK4.O",
"V6L-1V6R-14":"SD1.CK14D1.CK4.P",
"V6L-1V6R-15":"SD1.CK14D1.CK4.W",
"V6L-1V6R-16":"SD1.CK14D1.CK14.B",
"V6L-1V6R-17":"SD1.CK14D1.CK14.G",
"V6L-1V6R-18":"D6",
"V6L-1V6R-2":"SD1.CK4D3.CK4.G",
"V6L-1V6R-3":"SD1.CK4D3.CK4.O",
"V6L-1V6R-4":"SD1.CK4D3.CK4.P",
"V6L-1V6R-5":"SD1.CK4D3.CK14.O",
"V6L-1V6R-6":"SD1.CK4D3.CK14.P",
"V6L-1V6R-8":"SD1.CK4D3.CK4.W",
"V6L-1V6R-9":"H8",
"V6L-2V6R-0":"SD1.CK14D1.CK14.P",
"V6L-2V6R-1":"SD1.CK14D1.CK14.W",
"V6L-2V6R-10":"SD1.CK14D2.CK14.P",
"V6L-2V6R-11":"SD1.CK14D2.CK14.W",
"V6L-2V6R-12":"SD1.CK14D3.CK4.B",
"V6L-2V6R-13":"SD1.CK14D3.CK4.G",
"V6L-2V6R-14":"SD1.CK14D3.CK4.O",
"V6L-2V6R-15":"SD1.CK14D3.CK4.P",
"V6L-2V6R-16":"SD1.CK14D3.CK4.W",
"V6L-2V6R-17":"SD1.CK14D3.CK14.B",
"V6L-2V6R-18":"D7",
"V6L-2V6R-2":"SD1.CK14D2.CK4.B",
"V6L-2V6R-3":"SD1.CK14D2.CK4.G",
"V6L-2V6R-4":"SD1.CK14D2.CK4.O",
"V6L-2V6R-5":"SD1.CK14D2.CK4.P",
"V6L-2V6R-6":"SD1.CK14D2.CK4.W",
"V6L-2V6R-7":"SD1.CK14D2.CK14.B",
"V6L-2V6R-8":"SD1.CK14D2.CK14.G",
"V6L-2V6R-9":"SD1.CK14D2.CK14.O",
"V6L-3V6R-0":"SD1.CK14D3.CK14.O",
"V6L-3V6R-1":"SD1.CK14D3.CK14.P",
"V6L-3V6R-10":"SD1.CK4H1.CK14.O",
"V6L-3V6R-11":"SD1.CK4H1.CK14.P",
"V6L-3V6R-12":"D8",
"V6L-3V6R-13":"SD1.CK4H2.CK4.B",
"V6L-3V6R-14":"SD1.CK4H2.CK4.O",
"V6L-3V6R-15":"SD1.CK4H2.CK4.G",
"V6L-3V6R-16":"SD1.CK4H2.CK4.P",
"V6L-3V6R-17":"SD1.CK4H2.CK4.W",
"V6L-3V6R-18":"D9",
"V6L-3V6R-2":"SD1.CK14D3.CK14.W",
"V6L-3V6R-3":"SD1.CK4H1.CK4.B",
"V6L-3V6R-4":"SD1.CK4H1.CK4.G",
"V6L-3V6R-5":"SD1.CK4H1.CK4.O",
"V6L-3V6R-6":"SD1.CK4H1.CK4.P",
"V6L-3V6R-7":"SD1.CK4H1.CK4.W",
"V6L-3V6R-8":"SD1.CK4H1.CK14.B",
"V6L-3V6R-9":"SD1.CK4H1.CK14.G",
"V6L-4V6R-0":"SD1.CK4H2.CK14.G",
"V6L-4V6R-1":"SD1.CK4H2.CK14.O",
"V6L-4V6R-10":"SD1.CK4H3.CK4.O",
"V6L-4V6R-11":"SD1.CK4H3.CK4.P",
"V6L-4V6R-12":"A1.DOSE.DISEASE",
"V6L-4V6R-13":"SD1.CK4H3.CK14.W",
"V6L-4V6R-14":"SD1.CK14H1.CK4.B",
"V6L-4V6R-15":"SD1.CK14H1.CK4.G",
"V6L-4V6R-16":"SD1.CK14H1.CK4.O",
"V6L-4V6R-17":"SD1.CK14H1.CK4.P",
"V6L-4V6R-18":"SD1.CK14H1.CK4.W",
"V6L-4V6R-2":"SD1.CK4H2.CK14.P",
"V6L-4V6R-3":"D10",
"V6L-4V6R-4":"SD1.CK4H3.CK14.B",
"V6L-4V6R-5":"SD1.CK4H3.CK14.G",
"V6L-4V6R-6":"A1.+B.PD ",
"V6L-4V6R-7":"SD1.CK4H3.CK14.P",
"V6L-4V6R-8":"SD1.CK4H3.CK4.B",
"V6L-4V6R-9":"SD1.CK4H3.CK4.G",
"V6L-5V6R-0":"SD1.CK14H1.CK14.B",
"V6L-5V6R-1":"SD1.CK14H1.CK14.G",
"V6L-5V6R-10":"SD1.CK14H2.CK14.G",
"V6L-5V6R-11":"A1.PREDOSE.+AG",
"V6L-5V6R-12":"SD1.CK14H2.CK14.P",
"V6L-5V6R-13":"SD1.CK14H2.CK14.W",
"V6L-5V6R-14":"SD1.CK14H3.CK4.B",
"V6L-5V6R-15":"SD1.CK14H3.CK4.G",
"V6L-5V6R-16":"SD1.CK14H3.CK4.P",
"V6L-5V6R-17":"SD1.CK14H3.CK4.O",
"V6L-5V6R-18":"SD1.CK14H3.CK4.W",
"V6L-5V6R-2":"SD1.CK14H1.CK14.P",
"V6L-5V6R-3":"SD1.CK14H1.CK14.W",
"V6L-5V6R-4":"A1.PREDOSE.+P",
"V6L-5V6R-5":"SD1.CK14H2.CK4.G",
"V6L-5V6R-6":"SD1.CK14H2.CK4.O",
"V6L-5V6R-7":"SD1.CK14H2.CK4.P",
"V6L-5V6R-8":"SD1.CK14H2.CK4.W",
"V6L-5V6R-9":"A1.PREDOSE.+AB",
"V6L-6V6R-0":"SD1.CK14H1.CK14.B",
"V6L-6V6R-1":"SD1.CK14H1.CK14.G",
"V6L-6V6R-10":"SD1.CK4H1.CK14.B.T10",
"V6L-6V6R-11":"SD1.CK4H1.CK14.G.T10",
"V6L-6V6R-12":"A1.+.D2.BB",
"V6L-6V6R-13":"SD1.CK4H1.CK14.P.T10",
"V6L-6V6R-14":"SD1.CK4H1.CK14.W.T10",
"V6L-6V6R-15":"SD1.CK4H2.CK4.B.T10",
"V6L-6V6R-16":"SD1.CK4H2.CK4.G.T10",
"V6L-6V6R-17":"SD1.CK4H2.CK4.O.T10",
"V6L-6V6R-18":"SD1.CK4H2.CK4.P.T10",
"V6L-6V6R-2":"SD1.CK14H3.CK14.O",
"V6L-6V6R-3":"SD1.CK14H3.CK14.P",
"V6L-6V6R-4":"A1.+.D3.BB",
"V6L-6V6R-5":"A1.+.D3.PA",
"V6L-6V6R-6":"SD1.CK4H1.CK4.G.T10",
"V6L-6V6R-7":"SD1.CK4H1.CK4.O.T10",
"V6L-6V6R-8":"SD1.CK4H1.CK4.P.T10",
"V6L-6V6R-9":"SD1.CK4H1.CK4.W.T10",
"V6L-7V6R-0":"SD1.CK4H2.CK4.W.T10",
"V6L-7V6R-1":"SD1.CK4H2.CK14.B.T10",
"V6L-7V6R-10":"SD1.CK4H3.CK4.W.T10",
"V6L-7V6R-11":"SD1.CK4H3.CK14.B.T10",
"V6L-7V6R-12":"SD1.CK4H3.CK14.G.T10",
"V6L-7V6R-13":"SD1.CK4H3.CK14.O.T10",
"V6L-7V6R-14":"SD1.CK4H3.CK14.P.T10",
"V6L-7V6R-15":"SD1.CK4H3.CK14.W.T10",
"V6L-7V6R-16":"SD1.CK14H1.CK4.B.T10",
"V6L-7V6R-17":"SD1.CK14H1.CK4.G.T10",
"V6L-7V6R-18":"SD1.CK14H1.CK4.O.T10",
"V6L-7V6R-2":"A1.+.H1.BA",
"V6L-7V6R-3":"SD1.CK4H2.CK14.O.T10",
"V6L-7V6R-4":"SD1.CK4H2.CK14.P.T10",
"V6L-7V6R-5":"SD1.CK4H2.CK14.W.T10",
"V6L-7V6R-6":"SD1.CK4H3.CK4.B.T10",
"V6L-7V6R-7":"SD1.CK4H3.CK4.G.T10",
"V6L-7V6R-8":"SD1.CK4H3.CK4.O.T10",
"V6L-7V6R-9":"SD1.CK4H3.CK4.P.T10",
"V6L-8V6R-0":"SD1.CK14H1.CK4.P.T10",
"V6L-8V6R-1":"SD1.CK14H1.CK4.W.T10",
"V6L-8V6R-10":"SD1.CK14H2.CK4.P.T10",
"V6L-8V6R-11":"SD1.CK14H2.CK4.W.T10",
"V6L-8V6R-12":"SD1.CK14H2.CK14.B.T10",
"V6L-8V6R-13":"SD1.CK14H2.CK14.G.T10",
"V6L-8V6R-14":"SD1.CK14H2.CK14.O.T10",
"V6L-8V6R-15":"SD1.CK14H2.CK14.P.T10",
"V6L-8V6R-16":"SD1.CK14H2.CK14.W.T10",
"V6L-8V6R-17":"SD1.CK14H3.CK4.B.T10",
"V6L-8V6R-18":"SD1.CK14H3.CK4.G.T10",
"V6L-8V6R-2":"SD1.CK14H1.CK14.B.T10",
"V6L-8V6R-3":"SD1.CK14H1.CK14.G.T10",
"V6L-8V6R-4":"SD1.CK14H1.CK14.O.T10",
"V6L-8V6R-5":"SD1.CK14H1.CK14.P.T10",
"V6L-8V6R-6":"A1.+.D2.BB",
"V6L-8V6R-7":"SD1.CK14H2.CK4.B.T10",
"V6L-8V6R-8":"SD1.CK14H2.CK4.G.T10",
"V6L-8V6R-9":"SD1.CK14H2.CK4.O.T10",
"V6L-9V6R-0":"SD1.CK14H3.CK4.O.T10",
"V6L-9V6R-1":"SD1.CK14H3.CK4.P.T10",
"V6L-9V6R-10":"SD1.CK4D1.CK4.O.T10",
"V6L-9V6R-11":"SD1.CK4D1.CK4.P.T10",
"V6L-9V6R-12":"SD1.CK4D1.CK4.W.T10",
"V6L-9V6R-13":"SD1.CK4D1.CK14.B.T10",
"V6L-9V6R-14":"SD1.CK4D1.CK14.G.T10",
"V6L-9V6R-15":"SD1.CK4D1.CK14.O.T10",
"V6L-9V6R-16":"SD1.CK4D1.CK14.P.T10",
"V6L-9V6R-17":"SD1.CK4D1.CK14.W.T10",
"V6L-9V6R-18":"SD1.CK4D2.CK4.B.T10",
"V6L-9V6R-2":"SD1.CK14H3.CK4.W.T10",
"V6L-9V6R-3":"SD1.CK14H3.CK14.B.T10",
"V6L-9V6R-4":"SD1.CK14H3.CK14.G.T10",
"V6L-9V6R-5":"SD1.CK14H3.CK14.O.T10",
"V6L-9V6R-6":"SD1.CK14H3.CK14.P.T10",
"V6L-9V6R-7":"SD1.CK14H3.CK14.W.T10",
"V6L-9V6R-8":"SD1.CK4D1.CK4.B.T10",
"V6L-9V6R-9":"SD1.CK4D1.CK4.G.T10"
}

# Can also be passed as arguments
fadaptor = "ACGCTCTTCCGATCT"
radaptor = "AGATCGGAAGAGC"

files = {}

import gc
from argparse import ArgumentParser
from Bio import SeqIO
from time import time

            
def prepare_files(o):
    """
    Creates a set of files for Demultiplexing.
    """
    for key in samples:
        name = samples[key]
        
        ######
        # RENAME THESE
        ######
        t = o + name + "_SD1.fasta"
        
        # Directory path, lines list, lines written
        files[name] = [t,[], 1]
        
    files["NO_MATCH"] = [o + "NO_MATCH.fasta", [], 1]
    files["LOW_QUALITY_READS"] = [o + "LOW_QUALITY_READS.fasta", [], 1]
    
    for key in files:
        f = open(files[key][0], 'w')
        f.close()


def dump_files():
    """
    Dumps queued content into the appropriate file and garbage collects.
    """
    for key in files:
        t = files[key]
        f = open(t[0],'a')
        
        for line in t[1]:
            f.write(line)

        files[key][1] = []
        f.close()
    
    gc.collect()


def get_sample_name(read, filter_limit):
    """
    Returns the filename in which to dump the record
    """
    d = str(read.seq)
    leftBC = d[0:5]
    rightBC = d[-5:]
    
    if not ((leftBC in lbc) and (rightBC in rbc)):
        return "NO_MATCH"
    
    if compute_phred_score(read) <= filter_limit:
        return "LOW_QUALITY_READS"
        
    ids = lbc[leftBC] + rbc[rightBC]
    
    return samples[ids]
 
 
def compute_phred_score(rec):
    """
    Computes the phred score of the record by determning the lowest score in the barcodes
    """
    phreds = rec.letter_annotations["phred_quality"]
    
    # Maximum possible score
    low_score = 60
    
    for i in xrange(len(phreds)):
        low_score = min(phreds[i], low_score)

    return low_score
        
        
def trim_adaptors(record):
    """
    Removes the adaptors from the front and end of the read, if they exist.
    """
    indexf = record.seq.find(fadaptor)
    indexr = record.seq.find(radaptor)
    
    # Trim the rear adaptor first, obvi
    if indexr >= 0:
        record = record[:indexr]
    if indexf >= 0: 
        len_fadaptor = len(fadaptor)
        record = record[indexf+len_fadaptor:]

    return record
 

def format_line(read, tupl, name):
    """
    Converts a fastq line into a fasta line, removes the adapters, and appends the header name to
    the line.
    """
    read = read[5:-5]
    line = format(read, "fasta")
    line = "{0}{1}_{2}:{3}".format(line[0], name, tupl[2], line[1:])
    tupl[1].append(line)
    tupl[2] += 1
    return line


def process(infile, outdir, filter_limit):
    """
    Process the fastq file. For each record (read), we do the following:
    
    1) Trim the adaptors
    2) Remove the barcodes
    3) Lookup the appropriate file, filter if necessary
    4) Queue the file for writing
    """
    
    print "reading from: " + infile
    print "writing to: " + outdir
    
    prepare_files(outdir)
    
    count = 0
    filtered = 0
    data = SeqIO.parse(infile,"fastq")
    start = time()
    
    for read in data:
       
        count += 1
        if count % 25000 == 0 and count > 1:
            print "processed {0} reads in {1} seconds ({2}% have been filtered)".format(\
                str(count), str(time() - start), str(((filtered+0.0)/count)*100)) 
        
        # Trim the adaptor
        read = trim_adaptors(read)
        
        # Filter and find file
        name = get_sample_name(read, filter_limit)
        
        # Inefficient
        if (name == "LOW_QUALITY_READS"):
            filtered += 1

        # format and append
        format_line(read, files[name], name)
        
        # Arbitrary... but seems to prevent memory overloads.
        if count % 500000 == 0:
            print "DUMPING FILES... "
            dump_files()
            
    print "analysis complete, successfully processed {0} reads... writing to files".format(str(count))
    dump_files()
    
    print "files written, {0} reads processed, filtered {1}".format(str(count),str(filtered))
 
        
if __name__ == "__main__":
    
    parser = ArgumentParser(description="Process a fastq file according to gws.")
    parser.add_argument("infile", help="Input file. Should be fastq", type=str)
    
    #####
    #Putting in constants for now
    #parser.add_argument("front_adaptor", help="front adaptor for reads", type=str)
    #parser.add_argument("rear_adaptor", help="rear adaptor for reads", type=str)
    #####
    parser.add_argument("-o", "--outdir", help="Output filename directory", type=str, default="")
    parser.add_argument("-f", "--filter_limit", help="filter limit", type=int, default=20)
    args = parser.parse_args()
    
    infile = args.infile
    outdir = args.outdir if args.outdir[-1] == "/" else args.outdir + "/"
    filter_limit = args.filter_limit
    #fadaptor = args.front_adaptor
    #radaptor = args.rear_adaptor
    
    process(infile, outdir, filter_limit)