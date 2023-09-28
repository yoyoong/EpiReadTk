from enum import Enum

class CODE(Enum):
    IGNORE= 0
    UNMETHYLATED = 10
    METHYLATED = 11
    DELETION = -1
    A_REPLACE = 21
    T_REPLACE = 22
    C_REPLACE = 23
    G_REPLACE = 24
    AorG_REPLACE = 25
    CorT_REPLACE = 26
    A_INSERT = 31
    T_INSERT = 32
    C_INSERT = 33
    G_INSERT = 34
    UNKNOWN_INSERT = 35
CODE.IGNORE.label = 'Fx'
CODE.UNMETHYLATED.label = 'U'
CODE.METHYLATED.label = 'M'
CODE.DELETION.label1 = 'D'
CODE.DELETION.label2 = 'd'
CODE.A_REPLACE.label = 'A'
CODE.T_REPLACE.label = 'T'
CODE.C_REPLACE.label = 'C'
CODE.G_REPLACE.label = 'G'
CODE.AorG_REPLACE.label = 'R'
CODE.CorT_REPLACE.label = 'Y'
CODE.A_INSERT.label = 'a'
CODE.T_INSERT.label = 't'
CODE.C_INSERT.label = 'c'
CODE.G_INSERT.label = 'g'
CODE.UNKNOWN_INSERT.label = 'i'

CPG_DICT = {CODE.UNMETHYLATED.label: CODE.UNMETHYLATED.value, CODE.METHYLATED.label: CODE.METHYLATED.value}
SNP_DICT = {CODE.DELETION.label1: CODE.DELETION.value, CODE.DELETION.label2: CODE.DELETION.value,
            CODE.A_REPLACE.label: CODE.A_REPLACE.value, CODE.T_REPLACE.label: CODE.T_REPLACE.value,
            CODE.C_REPLACE.label: CODE.C_REPLACE.value, CODE.G_REPLACE.label: CODE.G_REPLACE.value,
            CODE.AorG_REPLACE.label: CODE.AorG_REPLACE.value, CODE.CorT_REPLACE.label: CODE.CorT_REPLACE.value,
            CODE.A_INSERT.label: CODE.A_INSERT.value, CODE.T_INSERT.label: CODE.T_INSERT.value,
            CODE.C_INSERT.label: CODE.C_INSERT.value, CODE.G_INSERT.label: CODE.G_INSERT.value,
            CODE.UNKNOWN_INSERT.label: CODE.UNKNOWN_INSERT.value}