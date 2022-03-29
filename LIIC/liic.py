import os
import re
import shutil
import sys
from itertools import islice
import copy 

import numpy as np
import pandas as pd

"""
This script is a linear interpolation between the two types of coordinates.
Meanwhile, Cartesian coordinates or internal coordinates can be given.
author: wuzb_1996@outlook.com 
"""

def coord_rot(coord1, coord2):
    """
    url=https://en.wikipedia.org/wiki/Kabsch_algorithm
    """
    # natom = len(coord1)
    coord_1 = copy.deepcopy(coord1)
    coord_2 = copy.deepcopy(coord2)
    coord_1 -= np.mean(coord_1, axis=0)
    coord_2 -= np.mean(coord_2, axis=0)
    cov_matrix = coord_1.T @ coord_2
    U, sigma, VT = np.linalg.svd(cov_matrix)
    d = np.sign(np.linalg.det(U @ VT))
    diag = np.identity(3)
    diag[2][2] = d
    rot_mat = VT.T @ diag @ U.T
    # rot_coord = coord_2 @ rot_matrix
    # RMSD = np.sqrt(np.sum((coord_1 - rot_coord) ** 2) / natom)
    return rot_mat

def translation(geom1, geom2):
    rot_mat = coord_rot(geom1, geom2)
    trans_mat =  geom1 - geom2 @ rot_mat
    return rot_mat, trans_mat



def rodrigues_rotation(axis, angle, vector):
    angle = np.radians(angle)
    axis = np.array(axis) / np.linalg.norm(axis)
    vector = np.array(vector) / np.linalg.norm(vector)
    V_rot = np.cos(angle) * vector + (1-np.cos(angle)) * np.dot(axis, vector) * axis \
        + np.sin(angle) * np.cross(axis, vector)
    return V_rot


def zmatrix_to_cartesian(atom, value, zmatrix):
    # the format of zmatrix : List(str)
    # C              
    # C                  1            B1
    # C                  2            B2    1            A1
    # C                  3            B3    2            A2      1        D1
    # ....
    # value :List(float)
    # B1 = .....
    # B2 = .....
    # .......
    # A1 = .....
    # A2 = .....
    # ........
    # D1 = .....
    # D2 = .....
    coord = []
    distance = 0.0
    angle = 0.0
    dihedral = 0.0
    for i in range(atom):
        if i >= 1:
            distance = value[i-1]
            if i >= 2:
                angle = value[i+atom - 3]
            if i >= 3:
                dihedral = value[i+atom*2 - 6]
        if i == 0:
            xyz = np.array([0.000000, 0.000000, 0.000000])
        elif i == 1:
            xyz = np.array([0.000000, 0.000000, distance])
        elif i == 2:
            atom1 = coord[int(zmatrix[i][1]) - 1]
            atom2 = coord[int(zmatrix[i][3]) - 1]
            axis = [0, 1, 0]
            u1 = rodrigues_rotation(axis, -angle, (atom1 - atom2))*distance
            u2 = atom1 - u1
            xyz = u2
        else:
            atom1 = coord[int(zmatrix[i][1]) - 1]
            atom2 = coord[int(zmatrix[i][3]) - 1]
            atom3 = coord[int(zmatrix[i][5]) - 1]
            q1 = atom2 - atom1
            q2 = atom3 - atom2
            u1 = rodrigues_rotation(np.cross(q1, q2), angle, q1)
            u2 = rodrigues_rotation(q1, -dihedral, u1) * distance + atom1
            xyz = u2
        coord.append(xyz)
    return coord

def get_geom(filename):
    with open(filename, 'r') as f:
        regex = re.compile('[A-Za-z]{1,2}\s*(\s*(-?[0-9]+\.[0-9]*)){3}')
        geom = []
        for line in f:
            if regex.search(line):
                xyz = list(map(float, line.split()[1:4]))
                geom.append(xyz)
        return np.array(geom)

def read_zmatrix(file, file_1):
    # file ->end.gjf file_1 = Gauss
    with open(file, 'r') as f, open(file_1, 'w+') as Gauss:
        flag = False
        zmatrix = []
        # regex=re.compile("(0\s*1)|(1\s*1)|(-1\s*1)|(0\s*2)(1\s*2)(-1\s*2)")
        regex = re.compile("-?[0-9]\s*[1-9]")
        l = 0
        for line in f:
            if flag:
                if line.split():
                    if l > 3:
                        zmatrix.append(line.split()[:7])
                    else:
                        n = 2 * l + 1
                        zmatrix.append(line.split()[:n])
                    l += 1
                else:
                    break
            else:
                Gauss.write(line)
            if regex.match(line):
                flag = True
    return l, zmatrix


def main():
    if not (os.path.exists(begin) and os.path.exists(end)):
        print(f'The file {begin} or {end} do not exist')
        sys.exit(1)
    # 读取gjf 格式内坐标内容
    l = -1
    with open('Gauss', 'w+') as Gauss, open(end, 'r') as A:
        for line in A:
            l += 1
            Lines = line.lstrip()
            if Lines[0:2] != 'B1':
                Gauss.write(line)
            else:
                break

    # 从第l+1行读取内坐标数据,保存为end.txt 和begin.txt
    tmp_file = [os.path.splitext(i)[0] for i in [begin, end]]
    file = ['begin', 'end']
    for i,j in zip(tmp_file,file):
        suffix_txt = str(j) + '.txt'
        suffix_gjf = str(i) + '.gjf'
        with open(suffix_txt, 'w+') as txt, open(suffix_gjf, 'r') as gjf:
            for index, line in enumerate(islice(gjf, l, None)):
                txt.write(line)

    name_begin = ['key', 'begin']
    name_end = ['key', 'end']
    data_begin = pd.read_csv('begin.txt', names=name_begin, sep='\s+')
    data_end = pd.read_csv('end.txt', names=name_end, sep='\s+')
    data = pd.merge(data_begin, data_end).set_index(
        ['key'])
    data['Delta'] = (data['end'] - data['begin'])/18  # 计算差分18次内坐标的Delta

    def adj(i):
        data_adjust = data[np.abs(data.loc[:, 'Delta']) > i]
        data_adjust_index = data_adjust.index
        if np.abs(i) > 11.5:
            for index in data_adjust_index:
                if data['Delta'][index] > 0:
                    data['begin'][index] = data['begin'][index] + 360
                else:
                    data['end'][index] = data['end'][index] + 360
        else:
            for index in data_adjust.index:
                if data['Delta'][index] < 0:
                    data['begin'][index] = data['begin'][index] - 360
        data['Delta'] = (data['end'] - data['begin'])/18
        print("Adjusted value")
        print(data_adjust)
    adj(12)  # 调整差值　
    adj(9.0)  # 使得关键二面角差值正负号一致。
    #          begin         end      Delta
    # key
    # D3  -179.315640  178.305304  19.867830
    # D4  -179.687982  179.129581  19.934309
    # D21 -160.096475  163.151186  17.958203
    # D5   179.818898 -179.302698 -19.951200

    def rotation(i):  # 调整 关键二面角差值的正负号，正负号对于两种旋转不同方向.
        data_adjust = data[np.abs(data.loc[:, 'Delta']) > 9]
        data_adjust_index = data_adjust.index
        if i == 1:
            for index in data_adjust_index:
                data['begin'][index] = data['begin'][index] + 360
        data['Delta'] = (data['end'] - data['begin'])/18
    if direction:
        rotation(1)
    data['Delta'] = (data['end'] - data['begin'])/Num_inter
    for i in range(0, Num_inter+1):
        num = i * Num_coef + Num_begin
        data[num] = data['begin'] + data['Delta'] * i

    Columns = data.columns.drop(
        ['begin', 'end', 'Delta'])  # 丢弃begin end Delta 三列
    
    if flag_remove:
        if os.path.exists("simulation.xyz"):
            print("The file 'simulation.xyz' exist and will be remove")
            os.remove("simulation.xyz")
    
    refer_geom = None
    if refer:
        if not os.path.exists(refer):
            raise Exception("The refer cartesian file %s does not exist" %(refer))
        else:
            refer_geom = get_geom(refer)
    

    for count, C in enumerate(Columns):  # Columns = [1,2...,N]
        i_str = str(C)
        head = i_str + '.gjf'
        tail = i_str + '.txt'
        if not flag_Z_C:
            shutil.copy('Gauss', head)
            data[C].to_csv(tail, header=False, sep=" ",
                           float_format='%.8f')  # 读取每列数据
        else:
            natom, zmatrix = read_zmatrix(end, 'Gauss')
            shutil.copy('Gauss', head)
            cartesian = zmatrix_to_cartesian(natom, data[C], zmatrix)
            
            if count == 0 and refer_geom:
                basis_geom = np.array(cartesian)
                rot_mat, trans_mat = translation(refer_geom, basis_geom)
            if refer_geom:
                # 调整分子趋向，与参考分子趋向一致。
                cartesian = np.array(cartesian) @ rot_mat + trans_mat
            
            with open(tail, 'w+') as file, open('simulation.xyz', 'a+') as file_1:
                file_1.write(str(natom)+'\n' + f'Time {str(count):>5s}'+'\n')
                for i in range(natom):
                    line = format(
                        zmatrix[i][0], '5s')+''.join(format(x, '>18.8f') for x in cartesian[i])
                    file.write(line+'\n')
                    file_1.write(line+'\n')
                file.write('\n')

        with open(head, 'a+') as f_1, open(tail, 'r') as f_2:
            f_1.write(f_2.read())
        os.remove(tail)

    # 删除临时文件
    for i in ['begin.txt', 'Gauss', 'end.txt']: 
        os.remove(i)

    # print(data_re)
    # print(data_Delta)
    #             begin       end     Delta
    # key
    # D9  -174.469371   -0.320910  9.674915
    # D10    3.409999  179.933096  9.806839
    #print(data.loc[np.abs(data.loc[:, 'Delta']) > 11 ][['begin', 'end', 'Delta']])
    print(data.loc[np.abs(data.loc[:, 'Delta']) > 1]
          [['begin', 'end', 'Delta']])
    # print(data['Delta'])


if __name__ == '__main__':
    direction = True  # 两种旋转不同方向 "True or False"
    Num_inter = 18  # 差分次数
    Num_begin = 0  # 起始编号
    Num_coef = 1  # 编号系数( 文件编号：a+bx  (a-> 起始编号，b->差分次数 ,x 编号系数)
    flag_Z_C = True  # zmatrix->cartesian
    begin = 'Np-Fu.gjf' # 初始构型
    end = 'Np-Fu-S1S0.gjf' #结束构型
    flag_remove = True # 是否移除文件"simulation.xyz"
    refer = None   # 初始构型未转化为Z坐标后的笛卡尔坐标
    main()
