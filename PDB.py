import requests
import threading
from multiprocessing import *
import os
import torch
import time
from protein import Protein
"""
Класс базы данных PDB
1) Предоставляет интерфейс взаимодействия с удаленным хранилищем белковых структур
2) Предоставляет интерфейс подготовки сырых данных
"""
path_delimiter = '/'
ptStoragePath = "./easy_ptstorage"
pdb_data_path = "./PDB_DATABASE/pdb/"

class PDB:

    link = None
    path = None

    # Initial method of PDB class
    # if init method is local_storage class will use local replica of PDB database instead of original

    def __init__(self, link=None, init_method=None, local_path=None):
        global path_delimiter
        if init_method is None:
            self.link = link
        elif init_method == "local_storage":
            self.path = local_path
            if self.path[-1] != path_delimiter:
                self.path += path_delimiter

    def pdb_to_nlp_training_data(self, target_file_path):
        objects = open(target_file_path, "r").readlines()
        pdb_path = "/home/main/V-Fold/"
        for object in objects:
            pass

    # Метод преобразованиея данных из файлы PDB к карте расстояний с отсечкой
    # Отсечка задается параметром one_dimension_input_size
    def prepare_data(self, one_dimension_input_size=50):
        global ptStoragePath, path_delimiter
        if not os.path.exists(ptStoragePath):
            os.makedirs(ptStoragePath)
        if ptStoragePath[:-1] != path_delimiter:
            ptStoragePath += path_delimiter
        if not os.path.exists(ptStoragePath+"pt"):
            os.makedirs(ptStoragePath+"pt")
        if not os.path.exists(ptStoragePath+"des"):
            os.makedirs(ptStoragePath+"des")
        fileList = os.listdir(self.path)
        for fileName in fileList:
            p = Protein(pdb_link=self.path + fileName)
            if one_dimension_input_size < p.getCATraceLen():
                different = p.getCATraceLen() - one_dimension_input_size
                for i in range(different):
                    input_tensor = p.generateInputData(one_dimension_input_size, i)
                    ptfname = ptStoragePath+"pt"+path_delimiter+fileName+".{}.pt".format(i)
                    desfname = ptStoragePath + "des" + path_delimiter + fileName + ".{}.pt.des".format(i)
                    torch.save(input_tensor, ptfname)
                    desired = p.generateDistanceMatrix(one_dimension_input_size, i)
                    torch.save(desired, desfname)
                    print(f"Save: {ptfname}")
            elif one_dimension_input_size == p.getCATraceLen():
                input_tensor = p.generateInputData(one_dimension_input_size, 0)
                ptfname = ptStoragePath+"pt"+ path_delimiter + fileName + ".pt"
                desfname = ptStoragePath + "des"+ path_delimiter + fileName + ".pt.des"
                torch.save(input_tensor, ptfname)
                desired = p.generateDistanceMatrix(one_dimension_input_size, 0)
                torch.save(desired, desfname)
                print(f"Save: {ptfname}")
            else:
                print(f"Flush input: {fileName}")

    # Подготовка данных без отсечки. Сохраняется вся карта расстояний любой длины
    def prepare_data(self):
        global ptStoragePath, path_delimiter
        if not os.path.exists(ptStoragePath):
            os.makedirs(ptStoragePath)
        if ptStoragePath[:-1] != path_delimiter:
            ptStoragePath += path_delimiter
        if not os.path.exists(ptStoragePath+"pt"):
            os.makedirs(ptStoragePath+"pt")
        if not os.path.exists(ptStoragePath+"des"):
            os.makedirs(ptStoragePath+"des")
        fileList = os.listdir(self.path)
        for fileName in fileList:
            p = Protein(pdb_link=self.path + fileName)
            for i in range(p.getCATraceLen()):
                input_tensor = p.generateCompleteProteinTensorFormLanguage()
                ptfname = ptStoragePath+"pt"+path_delimiter+fileName+".{}.pt".format(i)
                desfname = ptStoragePath + "des" + path_delimiter + fileName + ".{}.pt.des".format(i)
                torch.save(input_tensor, ptfname)
                desired = p.generateDistanceMatrix(one_dimension_input_size, i)
                torch.save(desired, desfname)
                print(f"Save: {ptfname}")

    def printAllCATraceLen(self):
        fileList = os.listdir(self.path)
        for fileName in fileList:
            p = Protein(pdb_link=self.path + fileName)
            print(f"CA Trace len of {fileName}: {p.getCATraceLen()}")

    # download files which noted in files list from PDB
    # abs_path is absolute_path of files list
    # save_path is folder which will be contain saved files
    def download_from_file_list(self, abs_path, save_path):
        global path_delimiter
        pdb_file_link = "https://files.rcsb.org/download/"
        f = open(abs_path, "r")
        content = f.readlines()
        all_count = len(content)
        for i, pdb_id in enumerate(content):
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            if save_path[-1] != path_delimiter:
                save_path += path_delimiter
            pdb_id = pdb_id.replace('\n', '')
            pdb_id = pdb_id[0:4]
            complete = float(i*100)/all_count
            tmp_link = pdb_file_link+pdb_id+".pdb"
            print(f"Load: {pdb_id} complete: {complete} link: {tmp_link}")
            pdb_content = requests.get(tmp_link)
            pdb_file = open(save_path+pdb_id+".pdb", 'wb+')
            pdb_file.write(pdb_content.content)
            pdb_file.close()

    def prepare_data(self, one_dimension_input_size=50):
        global ptStoragePath, path_delimiter
        if not os.path.exists(ptStoragePath):
            os.makedirs(ptStoragePath)
        if ptStoragePath[:-1] != path_delimiter:
            ptStoragePath += path_delimiter
        if not os.path.exists(ptStoragePath+"pt"):
            os.makedirs(ptStoragePath+"pt")
        if not os.path.exists(ptStoragePath+"des"):
            os.makedirs(ptStoragePath+"des")
        fileList = os.listdir(self.path)
        for fileName in fileList:
            p = Protein(pdb_link=self.path + fileName)
            if one_dimension_input_size < p.getCATraceLen():
                different = p.getCATraceLen() - one_dimension_input_size
                for i in range(different):
                    input_tensor = p.generateInputData(one_dimension_input_size, i)
                    ptfname = ptStoragePath+"pt"+path_delimiter+fileName+".{}.pt".format(i)
                    desfname = ptStoragePath + "des" + path_delimiter + fileName + ".{}.pt.des".format(i)
                    torch.save(input_tensor, ptfname)
                    desired = p.generateDistanceMatrix(one_dimension_input_size, i)
                    torch.save(desired, desfname)
                    print(f"Save: {ptfname}")
            elif one_dimension_input_size == p.getCATraceLen():
                input_tensor = p.generateInputData(one_dimension_input_size, 0)
                ptfname = ptStoragePath+"pt"+ path_delimiter + fileName + ".pt"
                desfname = ptStoragePath + "des"+ path_delimiter + fileName + ".pt.des"
                torch.save(input_tensor, ptfname)
                desired = p.generateDistanceMatrix(one_dimension_input_size, 0)
                torch.save(desired, desfname)
                print(f"Save: {ptfname}")
            else:
                print(f"Flush input: {fileName}")

    def fileIteration(self, fileName, one_dimension_input_size):
        global ptStoragePath, path_delimiter
        p = Protein(pdb_link=self.path + fileName)
        if one_dimension_input_size < p.getCATraceLen():
            different = p.getCATraceLen() - one_dimension_input_size
            for i in range(different):
                input_tensor = p.generateInput2DData(one_dimension_input_size, i)
                ptfname = ptStoragePath + "pt" + path_delimiter + fileName + ".{}.pt".format(i)
                desfname = ptStoragePath + "des" + path_delimiter + fileName + ".{}.pt.des".format(i)
                torch.save(input_tensor, ptfname)
                desired = p.generateDistanceMatrix(one_dimension_input_size, i)
                torch.save(desired, desfname)
                print(f"Save: {ptfname}")
        elif one_dimension_input_size == p.getCATraceLen():
            input_tensor = p.generateInput2DData(one_dimension_input_size, 0)
            ptfname = ptStoragePath + "pt" + path_delimiter + fileName + ".pt"
            desfname = ptStoragePath + "des" + path_delimiter + fileName + ".pt.des"
            torch.save(input_tensor, ptfname)
            desired = p.generateDistanceMatrix(one_dimension_input_size, 0)
            torch.save(desired, desfname)
            print(f"Save: {ptfname}")
        else:
            print(f"Flush input: {fileName}")

    def prepare_2d_data(self, one_dimension_input_size=50):
        global ptStoragePath, path_delimiter
        if not os.path.exists(ptStoragePath):
            os.makedirs(ptStoragePath)
        if ptStoragePath[:-1] != path_delimiter:
            ptStoragePath += path_delimiter
        if not os.path.exists(ptStoragePath+"pt"):
            os.makedirs(ptStoragePath+"pt")
        if not os.path.exists(ptStoragePath+"des"):
            os.makedirs(ptStoragePath+"des")
        fileList = os.listdir(self.path)
        for fileName in fileList:
            while threading.active_count() >= cpu_count()*2:
                time.sleep(1)
            t = threading.Thread(target=self.fileIteration, args=(fileName,one_dimension_input_size, ))
            t.start()
        while threading.active_count() > 1:
            time.sleep(1)
            print(f"Active threads{threading.active_count()}", end='\r')



pdb = PDB()
pdb.download_from_file_list('./targets.txt', './targets')
# pdb = PDB(init_method="local_storage", local_path=".\\easy_target")
# pdb.prepare_2d_data()
# pdb.printAllCATraceLen()

