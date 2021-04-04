import torch
import math as m
import os

naigbour_distance = 0.3

class EntropyAnalyzer:

    frame_number = 1

    # 4 dimension of tensor
    _CA = 1
    _N = 2
    _C = 3
    _CB = 4
    _OTHER = 5
    _OH2 = 6
    _H1 = 7
    _H2 = 8

    # 5 dimension of tensor
    _first_line = 1
    _second_line = 2
    _third_line = 3
    _excessive = 4

    _x = 0
    _y = 1
    _z = 2
    _name = 3

    deltaB = 3.5
    hydrogen_bond_len = 4

    def __init__(self):
        self.log_path = "C:\\logs\\"
        self.folder_delimiter = '\\'
        if self.log_path[-1] != self.folder_delimiter:
            self.log_path += self.folder_delimiter
        if not os.path.exists(self.log_path):
            os.makedirs(self.log_path)

    def analyze(self, pdb_frames_path):
        self.protein_points_count = 0
        self.protein_cursor = 0
        self.water_points_count = 0
        self.water_cursor = 0
        with open(pdb_frames_path, 'r') as frames_file:
            for line in frames_file:
                if self.is_water(line):
                    self.water_points_count += 1
                else:
                    self.protein_points_count += 1
                if "END" in line:
                    break
            self.protein_data = torch.zeros(self.protein_points_count, 5, device='cuda:0')
            self.water_data = torch.zeros(self.water_points_count, 5, device='cuda:0')
            for line in frames_file:
                if "END" in line:
                    self.proceed_data()
                    self.truncate_data()
                    continue
                if self.is_water(line):
                    self.update_water_data(line)
                else:
                    self.update_protein_data(line)
            frames_file.close()

    def update_water_data(self, line):
        self.water_data[self.water_cursor, self._x] = float(line[31:38].replace(' ', ''))
        self.water_data[self.water_cursor, self._y] = float(line[38:46].replace(' ', ''))
        self.water_data[self.water_cursor, self._z] = float(line[46:55].replace(' ', ''))
        name = line[13:17].replace(' ', '')
        if name == 'OH2':
            self.water_data[self.water_cursor, self._name] = self._OH2
        elif name == 'H1':
            self.water_data[self.water_cursor, self._name] = self._H1
        elif name == 'H2':
            self.water_data[self.water_cursor, self._name] = self._H2

    def update_protein_data(self, line):
        self.water_data[self.water_cursor, self._x] = float(line[30:38].replace(' ', ''))
        self.water_data[self.water_cursor, self._y] = float(line[38:46].replace(' ', ''))
        self.water_data[self.water_cursor, self._z] = float(line[46:55].replace(' ', ''))
        name = line[13:17].replace(' ', '')
        if name == 'N':
            self.water_data[self.water_cursor, self._name] = self._N
        elif name == 'CA':
            self.water_data[self.water_cursor, self._name] = self._CA
        elif name == 'C':
            self.water_data[self.water_cursor, self._name] = self._C
        elif name == 'CB' or name == 'CH2':
            self.water_data[self.water_cursor, self._name] = self._CB

    def is_water(self, line):
        if line[17:22] == 'TIP3W':
            return True
        else:
            return False

    def is_end_of_frame(self, line):
        if line.find("END") != -1:
            print(line)
            return True
        return False

    def proceed_data(self):
        pass

    def distance(self, a, b):
        return m.sqrt( (a[self._x] - b[self._x])**2 + (a[self._y] - b[self._y])**2 + (a[self._z] - b[self._z])**2 )

    def update_water_array_if_first_line(self, line):
        x = float(line[31:38].replace(' ', ''))
        y = float(line[38:46].replace(' ', ''))
        z = float(line[46:55].replace(' ', ''))
        for protein_point in self.protein_data:
            if self.distance(protein_point, (x, y, z)) < self.hydrogen_bond_len:
                self.water_data[self.water_cursor, self._x] = x
                self.water_data[self.water_cursor, self._y] = y
                self.water_data[self.water_cursor, self._z] = z
                name = line[13:17].replace(' ', '')
                if name == 'OH2':
                    self.water_data[self.water_cursor, self._name] = self._OH2
                elif name == 'H1':
                    self.water_data[self.water_cursor, self._name] = self._H1
                elif name == 'H2':
                    self.water_data[self.water_cursor, self._name] = self._H2
                break

    def is_first_line_molecula(self, line):
        x = float(line[31:38].replace(' ', ''))
        y = float(line[38:46].replace(' ', ''))
        z = float(line[46:55].replace(' ', ''))
        for protein_point in self.protein_data:
            if self.distance(protein_point, (x, y, z)) < self.hydrogen_bond_len:
                return True

    def update_protein_coordinate_indexes(self, line):
        if not 'ATOM' in line:
            return
        self.x[self.protein_cursor] = float(line[30:38].replace(' ', ''))
        self.y[self.protein_cursor] = float(line[38:46].replace(' ', ''))
        self.z[self.protein_cursor] = float(line[46:55].replace(' ', ''))
        self.protein_amino_acid_id[self.protein_cursor] = int(line[23:27].replace(' ', ''))
        name = line[13:17].replace(' ', '')
        if name == 'CA':
            self.protein_atom_type[self.protein_cursor] = self._CA
        elif name == 'C':
            self.protein_atom_type[self.protein_cursor] = self._C
        elif name == 'CB':
            self.protein_atom_type[self.protein_cursor] = self._CB
        elif name == 'HA2':
            self.protein_atom_type[self.protein_cursor] = self._CB
        elif name == 'N':
            self.protein_atom_type[self.protein_cursor] = self._N
        else:
            self.protein_atom_type[self.protein_cursor] = self._OTHER


    def update_water_coordinate_indexes(self, line):
        if self.water_cursor == 623938:
            print(line)
        if not 'ATOM' in line:
            return
        self.hoh_x[self.water_cursor] = float(line[30:38].replace(' ', ''))
        self.hoh_y[self.water_cursor] = float(line[38:46].replace(' ', ''))
        self.hoh_z[self.water_cursor] = float(line[46:55].replace(' ', ''))
        self.hoh_indexes[self.water_cursor] = int(line[23:27].replace(' ', ''))
        res = line[17:21].replace(' ', '')
        if res == 'H1':
            self.hoh_atom_names[self.water_cursor] = self._H1
        elif res == 'H2':
            self.hoh_atom_names[self.water_cursor] = self._H2
        elif res == 'OH2':
            self.hoh_atom_names[self.water_cursor] = self._OH2

    def binary_search(self, val, target):
        target_len = int(len(target) / 2)
        if target_len == 0:
            return -1
        sp = target[:target_len]
        if val < sp[-1]:
            res = self.binary_search(val, sp)
            if res < 0:
                return -1
            else:
                return res
        elif val > sp[-1]:
            sp = target[target_len:]
            res = self.binary_search(val, sp)
            if res < 0:
                return -1
            else:
                return res
        return sp[-1]

    def truncate_data(self):
        self.protein_data = torch.zeros(self.protein_points_count, 5, device='cuda:0')
        self.water_data = torch.zeros(self.water_points_count, 5, device='cuda:0')

    def scalar_multiplication(self, ax, ay, az, bx, by, bz):
        return ax*bx + ay*by + az*bz

    def module(self, ax, ay, az):
        return torch.sqrt(ax.pow(2) + ay.pow(2) + az.pow(2))

    def calculate_cosa(self, ax, ay, az, bx, by, bz):
        return self.scalar_multiplication(ax, ay, az, bx, by, bz) / self.module(ax, ay, az) * self.module(bx, by, bz)

    def x_rotation(self, ax, ay, az, cosa, sina):
        x = ax
        y = ay*cosa - az*sina
        z = ay*sina + az*cosa
        return x, y, z

    def y_rotation(self, ax, ay, az, cosa, sina):
        x = ax*cosa - az*sina
        y = ay
        z = ax*sina + az*cosa
        return x, y, z

    def z_rotation(self, ax, ay, az, cosa, sina):
        x = ax * cosa - ay * sina
        y = ax * sina + ay * cosa
        z = az
        return x, y, z


    def cuda_fill_protein(self, pdb_frames_path):

        self.protein_points_count = 0
        self.water_points_count = 0

        self.protein_cursor = 0
        self.water_cursor = 0

        with open(pdb_frames_path, 'r') as frames_file:
            for line in frames_file:
                if self.is_water(line):
                    break
                else:
                    self.protein_points_count += 1

            frames_file.seek(0)
            self.x = torch.zeros(self.protein_points_count).cuda()
            self.y = torch.zeros(self.protein_points_count).cuda()
            self.z = torch.zeros(self.protein_points_count).cuda()

            self.protein_atom_type = torch.zeros(self.protein_points_count).cuda()
            self.protein_amino_acid_id = torch.zeros(self.protein_points_count).cuda()

            self.protein_length = self.protein_points_count

            print(f"protein length: {self.protein_length}")

            for line in frames_file:
                if self.is_water(line):
                    break
                else:
                    self.update_protein_coordinate_indexes(line)
                    self.protein_cursor += 1

            print(self.x, self.y, self.z)
            frames_file.seek(0)
            for line in frames_file:
                if self.is_water(line):
                    self.water_points_count += 1
                if self.is_end_of_frame(line):
                    break

            frames_file.seek(0)

            self.water_length = self.water_points_count

            self.hoh_x = torch.zeros(self.water_length).cuda()
            self.hoh_y = torch.zeros(self.water_length).cuda()
            self.hoh_z = torch.zeros(self.water_length).cuda()
            self.hoh_indexes = torch.zeros(self.water_length).cuda()
            self.hoh_atom_names = torch.zeros(self.water_length).cuda()

            print(f"water length: {self.water_length}")
            self.water_cursor = 0
            self.water_initiated = False
            for line in frames_file:
                if self.is_water(line):
                    self.update_water_coordinate_indexes(line)
                    self.water_cursor += 1
                    self.water_initiated = True
                elif self.is_end_of_frame(line) and self.water_initiated:
                    select_condition = self.check_neighbour()
                    self.main_calculations()
                    self.water_statistic_calculation()
                    self.water_cursor = 0

    def check_neighbour(self):
        global naigbour_distance
        print_counter = 0
        for atom_index in range(self.protein_length):
            d = ((self.hoh_x.pow(2) - self.x[atom_index].pow(2)) +
                (self.hoh_y.pow(2) - self.y[atom_index].pow(2)) +
                (self.hoh_z.pow(2) - self.z[atom_index].pow(2))).sqrt()
            d = d[d < naigbour_distance]
            print(d.size()[0], end=' ')
            print_counter += 1
            if print_counter % 50 == 0:
                print('\n', end='')
        print('\n', end='')
        print("Iteration complete")

    def main_calculations(self):
        CA_x = self.x[self.protein_atom_type == self._CA]
        CB_x = self.x[self.protein_atom_type == self._CB]
        C_x = self.x[self.protein_atom_type == self._C]
        N_x = self.x[self.protein_atom_type == self._N]

        CA_y = self.y[self.protein_atom_type == self._CA]
        CB_y = self.y[self.protein_atom_type == self._CB]
        C_y = self.y[self.protein_atom_type == self._C]
        N_y = self.y[self.protein_atom_type == self._N]

        CA_z = self.z[self.protein_atom_type == self._CA]
        CB_z = self.z[self.protein_atom_type == self._CB]
        C_z = self.z[self.protein_atom_type == self._C]
        N_z = self.z[self.protein_atom_type == self._N]

        # Запись контрольных расстояний
        start_distance = torch.sqrt(
            (CB_x - C_x).pow(2) + (CB_y - C_y).pow(2) + (CB_z - C_z).pow(2)
        )

        # Сдвиг начала координат в точку CA

        CB_x = CB_x - CA_x
        CB_y = CB_y - CA_y
        CB_z = CB_z - CA_z

        C_x = C_x - CA_x
        C_y = C_y - CA_y
        C_z = C_z - CA_z

        N_x = N_x - CA_x
        N_y = N_y - CA_y
        N_z = N_z - CA_z

        # Находим проекцию вектора C на плоскость XOY

        C_proj_x = C_x
        C_proj_y = C_y
        C_proj_z = torch.zeros(C_z.size()[0]).cuda()

        # Получаем значения синуса и косинуса угла между проекцией вектора С на плоскость ХОУ и осью Х
        self.cos_1 = self.calculate_cosa(
            C_x, C_y, torch.zeros(C_z.size()[0]).cuda(), # XOY_proj
            torch.zeros(C_z.size()[0]).cuda() + 1, # X_Direction
            torch.zeros(C_z.size()[0]).cuda(),
            torch.zeros(C_z.size()[0]).cuda()
        )

        self.sin_1 = torch.sqrt(1 - self.cos_1.pow(2))

        # Если Y > 0, необходимо вращение противоположное ориентации плоскости
        self.sin_1[C_y > 0] = self.sin_1[C_y > 0] * (-1)

        C_x, C_y, C_z = self.z_rotation(C_x, C_y, C_z, self.cos_1, self.sin_1)
        CB_x, CB_y, CB_z = self.z_rotation(CB_x, CB_y, CB_z, self.cos_1, self.sin_1)

        # Находим синус и косинус угла между проекцией осью Х и прямой вектором C
        self.cos_2 = self.calculate_cosa(
            C_x, C_y, C_z,
            torch.zeros(C_z.size()[0]).cuda() + 1,  # X_Direction
            torch.zeros(C_z.size()[0]).cuda(),
            torch.zeros(C_z.size()[0]).cuda()
        )
        self.sin_2 = torch.sqrt(1 - self.cos_2.pow(2))
        self.sin_2[C_z > 0] = self.sin_2[C_z > 0] * (-1)

        # Выполянем вращение вокруг Y
        C_x, C_y, C_z = self.y_rotation(C_x, C_y, C_z, self.cos_2, self.sin_2)
        CB_x, CB_y, CB_z = self.y_rotation(CB_x, CB_y, CB_z, self.cos_2, self.sin_2)

        # Ось Х совмещена с прямой вектором С. Требуется совмещение плоскости ХОУ с точкой CB вращением вокруг X
        # Находим проекцию вектора СВ на плоскость YOZ
        # Вычисляем синус и косинус угла между проекцией вектора CB на плоскость YOZ и осью Y
        self.cos_3 = self.calculate_cosa(
            torch.zeros(CB_x.size()[0]).cuda(), CB_y, CB_z,  # YOZ_Proj
            torch.zeros(CB_x.size()[0]).cuda(),
            torch.zeros(CB_x.size()[0]).cuda() + 1,  # Y_Direction
            torch.zeros(CB_x.size()[0]).cuda()
        )
        self.sin_3 = torch.sqrt(1 - self.cos_3.pow(2))
        self.sin_3[C_z > 0] = self.sin_3[C_z > 0] * (-1)

        # Выполянем вращение вокруг X
        CB_x, CB_y, CB_z = self.x_rotation(CB_x, CB_y, CB_z, self.cos_3, self.sin_3)

        control_distance = torch.sqrt(
            (CB_x - C_x).pow(2) + (CB_y - C_y).pow(2) + (CB_z - C_z).pow(2)
        )
        check_vector = torch.abs(start_distance - control_distance)
        check = check_vector[check_vector > 1e-5].size()[0]
        if check > 0:
            print(f"Обнаружены ошибки в рассчетах: Размерность проверочного вектора не нулевая {check}")
            print(check_vector[check_vector > 1e-5])

    def water_statistic_calculation(self):  # Написать метод вычисления статистики водной ориентации
        global naigbour_distance
        print_counter = 0
        print_string = ""
        for atom_index in range(self.protein_length):
            # Вычисляем расстояние до каждой молекулы воды
            d = ((self.hoh_x.pow(2) - self.x[atom_index].pow(2)) +
                 (self.hoh_y.pow(2) - self.y[atom_index].pow(2)) +
                 (self.hoh_z.pow(2) - self.z[atom_index].pow(2))).sqrt()

            # Получаем молекулярные индексы каждого атома, попавшего в сферу
            indexes = self.hoh_indexes[d < naigbour_distance]
            indexes = indexes.unique()

            # По полученным индексам получаем целиком молекулы воды
            w_x = torch.tensor([]).cuda()
            w_y = torch.tensor([]).cuda()
            w_z = torch.tensor([]).cuda()
            w_n = torch.tensor([]).cuda()
            for i in indexes:
                w_x = torch.cat((w_x, self.hoh_x[self.hoh_indexes == i]))
                w_y = torch.cat((w_y, self.hoh_y[self.hoh_indexes == i]))
                w_z = torch.cat((w_z, self.hoh_z[self.hoh_indexes == i]))
                w_n = torch.cat((w_n, self.hoh_atom_names[self.hoh_indexes == i]))

            # Разделяем полученные массивы по именам атомов для начала вычисления ориентации
            h1_x = w_x[w_n == self._H1]
            h1_y = w_y[w_n == self._H1]
            h1_z = w_z[w_n == self._H1]

            h2_x = w_x[w_n == self._H2]
            h2_y = w_y[w_n == self._H2]
            h2_z = w_z[w_n == self._H2]

            oh2_x = w_x[w_n == self._OH2]
            oh2_y = w_y[w_n == self._OH2]
            oh2_z = w_z[w_n == self._OH2]

            # Начинаем вычисление водной ориентации
            # Выполняем сдвиг начала координат в атом кислорода для каждого водорода
            h1_x = h1_x - oh2_x
            h1_y = h1_y - oh2_y
            h1_z = h1_z - oh2_z

            h2_x = h2_x - oh2_x
            h2_y = h2_y - oh2_y
            h2_z = h2_z - oh2_z

            # Вычисляем расстояние для каждого водорода
            r1 = torch.sqrt(h1_x.pow(2) + h1_y.pow(2) + h1_z.pow(2))
            r2 = torch.sqrt(h2_x.pow(2) + h2_y.pow(2) + h2_z.pow(2))

            # Вычисляем зенит для каждого водорода
            teta1 = torch.arccos(h1_z / r1)
            teta2 = torch.arccos(h2_z / r2)

            # Вычисляем азимут для каждого водорода
            phi1 = torch.arctan(h1_y / h1_x)
            phi2 = torch.arctan(h2_y / h2_x)

            # Пишем лог вычисления
            self.print_log(teta1, phi1, teta2, phi2)

    def print_log(self, teta1, phi1, teta2, phi2):
        i = 0
        if not teta1.size()[0] == phi1.size()[0] == teta2.size()[0] == phi2.size()[0]:
            print("Wrong dimensions: log can`t be created. Calculation error was occured")
            quit()

        while i < teta1.size()[0]:
            log_line = f"{teta1[i]}\t{phi1[i]}\n{teta2[i]}\t{phi2[i]}\n"
            log_file = self.log_path + f"atom_{i}_log.txt"
            with open(log_file, 'a+') as f:
                f.write(log_line)
                f.close()
            i += 1



e = EntropyAnalyzer()
# e.fill_first_line_molecules("C:\\Users\\softc\\Desktop\\Data\\Models\\3rjp\\frames\\frames.pdb")
e.cuda_fill_protein("D:\\Data\\Models\\3rjp\\frames\\frames.pdb")
print("Complete")