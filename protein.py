"""
Класс объектов белков
Предоставляет структуру данных аминокислотной последовательности
Supported by Bespalov Andrei
"""

import torch
import torch.nn as nn
from PIL import Image
import math
import numpy as np

def sigmoid(x):
    return 1 / (1 + math.exp(-x))

class Resolve_Dicts:
    simple_numbers = {
        'A': 3,
        'R': 5,
        'N': 7,
        'D': 11,
        'B': 13,
        'C': 17,
        'E': 23,
        'Q': 29,
        'Z': 31,
        'G': 37,
        'H': 41,
        'I': 43,
        'L': 47,
        'K': 53,
        'M': 59,
        'F': 61,
        'P': 67,
        'S': 71,
        'T': 73,
        'W': 79,
        'Y': 83,
        'V': 89,
        'U': 97,
        'O': 101
    }
    backward_simple_numbers = {
        3: 'A',
        5: 'R',
        7: 'N',
        11: 'D',
        13: 'B',
        17: 'C',
        23: 'E',
        29: 'Q',
        31: 'Z',
        37: 'G',
        41: 'H',
        43: 'I',
        47: 'L',
        53: 'K',
        59: 'M',
        61: 'F',
        67: 'P',
        71: 'S',
        73: 'T',
        79: 'W',
        83: 'Y',
        89: 'V',
        97: 'U',
        101: 'O'
    }
    three_to_simple = {
        'GLY': 'G',
        'ALA': 'A',
        'VAL': 'V',
        'ILE': 'I',
        'LEU': 'L',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'CYS': 'C',
        'MET': 'M',
        'ASP': 'D',
        'ASN': 'N',
        'GLU': 'E',
        'GLN': 'Q',
        'LYS': 'K',
        'ARG': 'R',
        'HIS': 'H',
        'PHE': 'F',
        'TYR': 'Y',
        'TRP': 'W',
        'SEC': 'U',
        'PYR': 'O'
    }
    main_dict = {
        'A': 3,
        'R': 5,
        'N': 7,
        'D': 11,
        'B': 13,
        'C': 17,
        'E': 23,
        'Q': 29,
        'Z': 31,
        'G': 37,
        'H': 41,
        'I': 43,
        'L': 47,
        'K': 53,
        'M': 59,
        'F': 61,
        'P': 67,
        'S': 71,
        'T': 73,
        'W': 79,
        'Y': 83,
        'V': 89,
        'U': 97,
        'O': 101
    }
    coordinates = None
    atype = None
    numerical = None

class Residue:
    atoms = None

class Protein:

    content = None

    # protein entities
    sequence = None
    CA_trace = None
    AASequence = None
    CALen = 0
    distanceMatrix = None
    distanceTensor = None

    # Параметры преобразования данных
    cutoff = 5  # дистанция отсечки

    # Константы описания
    _on_ = 107
    _to_ = 109
    _is_ = 113
    _drop_ = 127
    _breakline_ = 131
    _next_number_ = 137

    def __init__(self, pdb_link=None, pdb_content=None):
        if pdb_link is not None:
            f = open(pdb_link, 'r')
            self.content = f.readlines()
            self.CA_trace = []
            self.__parse_content()
        elif pdb_content is not None:
            if type(pdb_content) == type([]):
                self.content = pdb_content
            elif type(pdb_content) == type(""):
                self.content = pdb_content.split("\n")
            self.CA_trace = []
            self.__parse_content()
        else:
            del self

    def __parse_content(self):
        self.AASequence = []
        for line in self.content:
            record = line[0:4]
            if record == "ATOM":
                chain = line[21]
                if chain == "A":
                    res_name = line[17:20]
                    name = line[12:16]
                    name = name.replace(' ', '')
                    if res_name != "HOH" and name == "CA":
                        self.AASequence.append(Resolve_Dicts.three_to_simple[res_name])
                        simple_number = Resolve_Dicts.simple_numbers[Resolve_Dicts.three_to_simple[res_name]]
                        serial = int(line[6:11])
                        seq_num = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[47:54])
                        self.CA_trace.append({
                            'simple_number': simple_number,
                            'res_name': Resolve_Dicts.three_to_simple[res_name],
                            'serial': serial,
                            'seq_num': seq_num,
                            'x': x,
                            'y': y,
                            'z': z
                        })
        self.CALen = len(self.CA_trace)

    def calculate_direction_value(self, seq):
        power_dict = {}
        direction_coef = 1
        for i, x in enumerate(seq):
            if x['simple_number'] in power_dict.keys():
                power_dict[x['simple_number']] += i
            else:
                power_dict[x['simple_number']] = i
        for base, power in power_dict.items():
            direction_coef *= sigmoid(base**sigmoid(power))
        return direction_coef

    def getCATrace(self):
        return self.CA_trace

    def getCATraceTensor(self):
        result = []
        for a in self.CA_trace:
            result.append(Resolve_Dicts.simple_numbers[a['res_name']])
        return torch.tensor(result)

    def printCATrace(self):
        for ca in self.CA_trace:
            print(ca)

    def getCATraceLen(self):
        return self.CALen

    # Генерация матрицы расстояний
    # с заданным размером начиная со
    # стартовой позиции в аминокислотной
    # последовательности
    def generateDistanceMatrix(self, coreLen, startPos):
        toRightPos = startPos + coreLen
        if toRightPos > self.CALen:
            raise IOError('Wrong core length')
        coreArray = []
        self.distanceMatrix = np.zeros([coreLen, coreLen])
        for i1, ca1 in enumerate(self.CA_trace):
            for i2, ca2 in enumerate(self.CA_trace):
                if i1 >= startPos and i1 < toRightPos and i2 >= startPos and i2 < toRightPos:
                    d = math.sqrt(((ca2['x'] - ca1['x'])**2) + ((ca2['y'] - ca1['y'])**2) + ((ca2['z'] - ca1['z'])**2))
                    self.distanceMatrix[i1-startPos, i2-startPos] = d
        self.distanceTensor = torch.from_numpy(self.distanceMatrix)
        return self.distanceTensor

    # Генерация текстового описания расстояний и контактов
    # Для последовательностей любой длинны
    def generateCompleteProteinFormLanguage(self):
        self.complete_sentence = ""
        for i1, ca1 in enumerate(self.CA_trace):
            first_line_sentence = f"{ca1['res_name']} "
            second_line_sentence = ""
            for i2, ca2 in enumerate(self.CA_trace):
                d = math.sqrt(
                    ((ca2['x'] - ca1['x']) ** 2) + ((ca2['y'] - ca1['y']) ** 2) + ((ca2['z'] - ca1['z']) ** 2))
                if d <= self.cutoff and abs(i1 - i2) > 4:
                    second_line_sentence += f"to {i2} is {ca2['res_name']} on {round(d)} "
            self.complete_sentence += first_line_sentence + second_line_sentence + "drop \n "
        return self.complete_sentence

    # Генерация тензорного описания расстояний и контактов
    # Для последовательностей любой длинны
    def generateCompleteProteinTensorFormLanguage(self):
        self.complete_sentence = []
        for i1, ca1 in enumerate(self.CA_trace):
            first_line_sentence = [Resolve_Dicts.simple_numbers[ca1['res_name']]]  # f"{ca1['res_name']} "
            second_line_sentence = []
            for i2, ca2 in enumerate(self.CA_trace):
                d = math.sqrt(
                    ((ca2['x'] - ca1['x']) ** 2) + ((ca2['y'] - ca1['y']) ** 2) + ((ca2['z'] - ca1['z']) ** 2))
                if d <= self.cutoff and abs(i1 - i2) > 4:
                    second_line_sentence += [self._to_, self._next_number_, i2, self._is_,
                                             Resolve_Dicts.simple_numbers[ca1['res_name']],
                                             self._on_,
                                             self._next_number_,
                                             round(d)]  # f"to {i2} is {ca2['res_name']} on {round(d)} "
            self.complete_sentence += first_line_sentence + second_line_sentence + [self._drop_, self._breakline_]
        return torch.tensor(self.complete_sentence)

    def restoreSequenceFromTensor(self, tensor):
        result_text = ""
        next_number = False
        for number in tensor:
            if next_number:
                result_text += " "+str(number.item())
                next_number = False
            else:
                if number.item() == self._on_:
                    result_text += " on"
                elif number.item() == self._next_number_:
                    next_number = True
                elif number.item() == self._drop_:
                    result_text += " drop"
                elif number.item() == self._is_:
                    result_text += " is"
                elif number.item() == self._to_:
                    result_text += " to"
                elif number.item() == self._breakline_:
                    result_text += " \n"
                else:
                    result_text += " "+Resolve_Dicts.backward_simple_numbers[number.item()]
        return result_text

    def generateCompleteDistanceMatrix(self):
        coreArray = []
        coreLen = len(self.CA_trace)
        self.distanceMatrix = np.zeros([coreLen, coreLen])
        for i1, ca1 in enumerate(self.CA_trace):
            for i2, ca2 in enumerate(self.CA_trace):
                d = math.sqrt(
                    ((ca2['x'] - ca1['x']) ** 2) + ((ca2['y'] - ca1['y']) ** 2) + ((ca2['z'] - ca1['z']) ** 2))
                if d <= self.cutoff:
                    self.distanceMatrix[i1, i2] = d
                else:
                    self.distanceMatrix[i1, i2] = 0
        self.distanceTensor = torch.from_numpy(self.distanceMatrix)
        return self.distanceTensor

    def generateInputData(self, coreLen, startPos):
        toRightPos = startPos + coreLen
        if toRightPos > self.CALen:
            raise IOError('Wrong core length')
        leftArray = []
        coreArray = []
        rightArray = []
        for i, a in enumerate(self.CA_trace):
            if i < startPos:
                leftArray.append(a)
            elif i >= startPos and i < toRightPos:
                coreArray.append(a)
            else:
                rightArray.append(a)
        left_value = self.calculate_direction_value(leftArray)
        right_value = self.calculate_direction_value(rightArray)
        content_simple_numbers_array = []
        for x in coreArray:
            content_simple_numbers_array.append(x['simple_number'])
        return torch.tensor(np.concatenate(([left_value], content_simple_numbers_array, [right_value])))

    def generateInput2DData(self, coreLen, startPos):
        toRightPos = startPos + coreLen
        if toRightPos > self.CALen:
            raise IOError('Wrong core length')
        leftArray = []
        coreArray = []
        rightArray = []
        for i, a in enumerate(self.CA_trace):
            if i < startPos:
                leftArray.append(a)
            elif i >= startPos and i < toRightPos:
                coreArray.append(a)
            else:
                rightArray.append(a)
        left_value = self.calculate_direction_value(leftArray)
        right_value = self.calculate_direction_value(rightArray)
        content_simple_numbers_array = []
        for x in coreArray:
            content_simple_numbers_array.append(x['simple_number'])
        t2 = torch.tensor(np.concatenate(([left_value], content_simple_numbers_array, [right_value])))
        t1 = torch.zeros(52, 52)
        for i in range(52):
            t1[i] = t2
        return t1


p = Protein(pdb_link="./training-set/1A73.pdb")
print(p.getCATraceTensor())
a = p.generateCompleteProteinTensorFormLanguage()
print(a)
# print(p.restoreSequenceFromTensor(a))