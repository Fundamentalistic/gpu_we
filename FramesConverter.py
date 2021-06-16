import os
import pymysql
import time
import re
from threading import Thread

mutex = False
water_id = 0

class FramesConverter:

    _CA = 1
    _N = 2
    _C = 3
    _CB = 4
    _OTHER = 5
    _OH2 = 6
    _H1 = 7
    _H2 = 8

    def __init__(self, path, frames_path):
        self.frames_path = frames_path
        if not os.path.exists(path):
            os.makedirs(path)
        self.path = path
        self.db = pymysql.connect(
            host='localhost',
            user='models',
            password='zxcasdqwe123'
        )
        model_index_store_filepath = 'model.index'
        if os.path.isfile(model_index_store_filepath):
            with open(model_index_store_filepath, 'r') as f:
                self.model_index = int(f.read())
                f.close()
            with open(model_index_store_filepath, 'w') as f:
                f.truncate()
                f.write(str(self.model_index + 1))
                f.close()
        else:
            self.model_index = 1
            with open(model_index_store_filepath, 'a+') as f:
                f.write(str(self.model_index + 1))
                f.close()



    def create_database(self):
        create_database_query_path = "create_database.sql"
        with open(create_database_query_path, 'r') as f:
            query = f.read()
        complete_query = query.format(self.model_index)
        dbname = f'modelling{self.model_index}'
        commands = complete_query.split(';')
        cursor = self.db.cursor()
        for command in commands:
            if command.strip() != '':
                cursor.execute(command)
        return dbname

    def fill_protein(self):
        self.protein_cursor = 0
        with open(self.frames_path) as f:
            for line in f:
                if self.is_water(line):
                    break
                if not 'ATOM' in line:
                    continue
                x = float(line[30:38].replace(' ', ''))
                y = float(line[38:46].replace(' ', ''))
                z = float(line[46:55].replace(' ', ''))
                res_id = [self.protein_cursor] = int(line[22:27].replace(' ', ''))
                atom_id = [self.protein_cursor] = self.protein_cursor
                name = line[13:17].replace(' ', '')
                if name == 'CA':
                    atom_name_id = self._CA
                elif name == 'C':
                    atom_name_id = self._C
                elif name == 'CB':
                    atom_name_id = self._CB
                elif name == 'HA2':
                    atom_name_id = self._CB
                elif name == 'N':
                    atom_name_id = self._N
                else:
                    atom_name_id = self._OTHER
                self.protein_cursor += 1

    def fill_water(self, file_path):
        self.water_cursor = 0
        with open(file_path) as f:
            for line in f:
                if "END" in line:
                    timestamps()
                    self.db.commit()
                if self.is_water(line):
                    if not 'ATOM' in line:
                        continue
                    x = float(line[30:38].replace(' ', ''))
                    y = float(line[38:46].replace(' ', ''))
                    z = float(line[46:55].replace(' ', ''))
                    res_id = int(line[22:27].replace(' ', ''))
                    atom_id = self.water_cursor
                    res = line[13:17].replace(' ', '')
                    if res == 'H1':
                        atom_name_id = self._H1
                    elif res == 'H2':
                        atom_name_id = self._H2
                    elif res == 'OH2':
                        atom_name_id = self._OH2
                    query = "INSERT INTO `water` (`res_id`,`atom_name_id`,`atom_id`,`x`,`y`,`z` )"
                    query += "VALUES ({0}, {1}, {2}, {3}, {4}, {5});"
                    query = query.format(res_id, atom_name_id, atom_id, x, y, z)
                    cursor = self.db.cursor()
                    cursor.execute(query)
                    self.water_cursor += 1


    def is_water(self, line):
        if line[17:22] == 'TIP3W':
            return True
        else:
            return False

    def set_database(self, host='localhost', user='models', password='zxcasdqwe123', database=''):
        self.db = pymysql.connect(
            host=host,
            user=user,
            password=password,
            database=database
        )

    def split_file(self):
        timestamps()
        file_index = 0
        line_counter = 0
        with open(self.frames_path, 'r') as main_file:
            fn = open(self.frames_path+f".{file_index}.pdb", 'a+')
            for line in main_file:
                fn.write(line)
                line_counter += 1
                if line_counter == 10000000:
                    timestamps()
                    fn.close()
                    file_index += 1
                    fn = open(self.frames_path + f".{file_index}.pdb", "a+")
                    line_counter = 0
            fn.close()
            main_file.close()

    def fill_water_from_decomposition(self, decomposition_path, db_name):
        timestamps()
        file_index = 0
        line_counter = 0
        if decomposition_path[-1] != '\\':
            decomposition_path += '\\'
        preg = re.compile('frames.pdb...')
        file_list = os.listdir(decomposition_path)
        target_list = []
        threads = []
        for file in file_list:
            if preg.match(file) != None:
                target_list.append(decomposition_path+file)
        for i, target in enumerate(target_list):
            index_of_water = i*10**7
            t = Thread(target=self.water_iteration, args=(target, index_of_water, i, db_name, ))
            t.start()
            threads.append(t)
        for t in threads:
            t.join()


    def water_iteration(self, target, water_id, i, db_name):
        global mutex
        print(f"Thread {water_id} is started")
        db = pymysql.connect(
            host='localhost',
            user='models',
            password='zxcasdqwe123',
            database=db_name,
        )
        with open(target) as f:
            for line in f:
                if "END" in line:
                    timestamps()
                    db.commit()
                if self.is_water(line):
                    if not 'ATOM' in line:
                        continue
                    x = float(line[30:38].replace(' ', ''))
                    y = float(line[38:46].replace(' ', ''))
                    z = float(line[46:55].replace(' ', ''))
                    res_id = int(line[22:27].replace(' ', ''))
                    atom_id = water_id
                    res = line[13:17].replace(' ', '')
                    if res == 'H1':
                        atom_name_id = self._H1
                    elif res == 'H2':
                        atom_name_id = self._H2
                    elif res == 'OH2':
                        atom_name_id = self._OH2
                    query = "INSERT INTO `water` (`res_id`,`atom_name_id`,`atom_id`,`x`,`y`,`z` )"
                    query += "VALUES ({0}, {1}, {2}, {3}, {4}, {5});"
                    query = query.format(res_id, atom_name_id, atom_id, x, y, z)
                    cursor = db.cursor()
                    cursor.execute(query)
                    water_id += 1
                    if water_id % 100000 == 0:
                        print(f"Thread {i} have proceedeng. Water identificator is: {water_id}")
            db.commit()
            db.close()
            f.close()


    


def timestamps():
    global current_time
    if not 'current_time' in globals():
        current_time = time.time()
    delta = time.time() - current_time
    current_time = time.time()
    print(f"TIMESTAMP: {delta}")
    return delta

timestamps()

path = 'C:\\frames\\frames_decomposition\\'
frames_path = 'C:\\frames\\frames.pdb'
f = FramesConverter(path=path, frames_path=frames_path)
dbname = f.create_database()
f.fill_water_from_decomposition('C:\\frames\\', dbname)