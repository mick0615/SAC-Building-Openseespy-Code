import openseespy.opensees as ops
import opsvis as opsv
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import time

class OpenseesModel:
    def __init__(self, Type, Analysis_Type, LCol, LCol_1, LGird, Story_num, Span_num, nps, Floor_Mass, damping_ratio, Earthquake):
        self.Type = Type
        self.Analysis_Type = Analysis_Type
        self.LCol = LCol
        self.LCol_1 = LCol_1
        self.LGird = LGird
        self.Story_num = Story_num
        self.Span_num = Span_num
        self.Node = []
        self.nps = nps
        self.Floor_Mass = Floor_Mass
        self.IDColTransf = 1
        self.IDGirdTransf = 2
        self.story_disp = []
        self.damping_ratio = damping_ratio
        self.Earthquake = Earthquake
        self.plot_title = ('Displacement' if self.Analysis_Type == 'disp' else 'Acceleration' if self.Analysis_Type == 'accel'
                            else 'Velocity' if self.Analysis_Type == 'vel' else 'Unknown Data Type')
        self.strain_stress_folder_name = os.path.join(os.getcwd(), 'Opensees Model', 'Inelastic Hysteresis Loop')
        os.makedirs(self.strain_stress_folder_name, exist_ok=True)
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        ops.geomTransf('Linear', self.IDColTransf)
        ops.geomTransf('Linear', self.IDGirdTransf)
        
    def create_node(self, i, j, x, y):
        nodeID = 10 * i + j
        ops.node(nodeID, x, y)
        return nodeID

    def build_nodes(self):
        self.Node = [self.create_node(1, j, (j - 1) * self.LGird, 0) for j in range(1, self.Span_num + 2)]
        [ops.fix(nodeID, 1, 1, 1) for nodeID in self.Node]
        
        for i in range(2, self.Story_num + 1):
            for j in range(1, self.Span_num + 2):
                self.create_node(i, j, (j - 1) * self.LGird, (i - 2) * self.LCol + self.LCol_1)
        
        [self.create_node(10, i, (i - 1) * self.LGird, self.LCol_1 + self.LCol * 8) for i in range(1, self.Span_num + 2)]
        
    
    def setup_perpendicular_direction(self):
        perpDirn = 1
        for i in range(2, self.Story_num + 2):
            nodeID = 1000 + i
            ops.node(nodeID, self.Span_num * self.LGird / 2, (i - 2) * self.LCol + self.LCol_1)
            ops.fix(nodeID, 0, 1, 1)
            
            node_floor = []
            for j in range(1, self.Span_num + 2):
                node_ID = i * 10 + j
                node_floor.append(node_ID)
            ops.rigidDiaphragm(perpDirn, nodeID, *node_floor)
            
        #print("Nodes have been added.")
            
    def define_materials_and_sections(self, Steel_Type):
        self.Steel_Type = Steel_Type
        Fy_A36 = 36.0 * 4448.22 * 39.37 * 39.37    # N/m^2
        Fy_A572 = 50.0 * 4448.22 * 39.37 * 39.37   # N/m^2
        self.Es = 2040 * 9800 * 10000              # N/m^2
        nu = 0.3
        Gs = self.Es / (2 * (1 + nu))
        b = 0.01
        R0 = 18.0
        cR1 = 0.925
        cR2 = 0.15
        a1 = 0.0
        a2 = 1.0
        a3 = 0.0
        a4 = 1.0
        sigInit = 0.0

        matIDSteel_Gird = 3
        matIDSteel_Col = 4

        if self.Steel_Type == 'Steel01':
            ops.uniaxialMaterial('Steel01', matIDSteel_Gird, Fy_A36, self.Es, b)
            ops.uniaxialMaterial('Steel01', matIDSteel_Col, Fy_A572, self.Es, b)
        elif self.Steel_Type == 'Steel02':
            ops.uniaxialMaterial('Steel02', matIDSteel_Gird, Fy_A36, self.Es, b, R0, cR1, cR2, a1, a2, a3, a4, sigInit)
            ops.uniaxialMaterial('Steel02', matIDSteel_Col, Fy_A572, self.Es, b, R0, cR1, cR2, a1, a2, a3, a4, sigInit)

        def RHsection(secID, matID, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ):
            ops.section('Fiber', secID, '-GJ', GJ)
            ops.patch('rect', matID, nftf, nfbf, -0.5 * d, -0.5 * bf, -0.5 * d + tf, 0.5 * bf)
            ops.patch('rect', matID, nfdw, nftw, -0.5 * d + tf, -0.5 * tw, 0.5 * d - tf, 0.5 * tw)
            ops.patch('rect', matID, nftf, nfbf, 0.5 * d - tf, -0.5 * bf, 0.5 * d, 0.5 * bf)
        
        SoF_dw = 0.01     # Fiber Spacing dw (m)
        SoF_tw = 0.0025   # Fiber Spacing tw (m)
        SoF_bf = 0.0175   # Fiber Spacing bf (m)
        SoF_tf = 0.00275  # Fiber Spacing tf (m)

        self.Col_sections = [
            {"tag": 370, "d": 17.9 * 0.0254, "bf": 16.5 * 0.0254, "tw": 1.66 * 0.0254, "tf": 2.66 * 0.0254, "J": 222.0 * 4.16e-7, "A": 109 * 0.0254 ** 2, "I": 5440 * 0.0254 ** 4},
            {"tag": 283, "d": 16.7 * 0.0254, "bf": 16.1 * 0.0254, "tw": 1.29 * 0.0254, "tf": 2.07 * 0.0254, "J": 104.0 * 4.16e-7, "A": 83.3 * 0.0254 ** 2, "I": 3840 * 0.0254 ** 4},
            {"tag": 257, "d": 16.4 * 0.0254, "bf": 16.0 * 0.0254, "tw": 1.18 * 0.0254, "tf": 1.89 * 0.0254, "J": 79.1 * 4.16e-7, "A": 75.6 * 0.0254 ** 2, "I": 3400 * 0.0254 ** 4},
            {"tag": 233, "d": 16.0 * 0.0254, "bf": 15.9 * 0.0254, "tw": 1.07 * 0.0254, "tf": 1.72 * 0.0254, "J": 59.5 * 4.16e-7, "A": 68.5 * 0.0254 ** 2, "I": 3010 * 0.0254 ** 4},
            {"tag": 500, "d": 19.6 * 0.0254, "bf": 17.0 * 0.0254, "tw": 2.19 * 0.0254, "tf": 3.5 * 0.0254, "J": 514.0 * 4.16e-7, "A": 147 * 0.0254 ** 2, "I": 8210 * 0.0254 ** 4},
            {"tag": 455, "d": 19.0 * 0.0254, "bf": 16.8 * 0.0254, "tw": 2.02 * 0.0254, "tf": 3.21 * 0.0254, "J": 395.0 * 4.16e-7, "A": 134 * 0.0254 ** 2, "I": 7190 * 0.0254 ** 4},
        ]
        
        for col_sec in self.Col_sections:
            ColSecTag = col_sec["tag"]
            d = col_sec["d"]
            bf = col_sec["bf"]
            tw = col_sec["tw"]
            tf = col_sec["tf"]
            dw = d - 2 * tf

            nfdw = int(dw // SoF_dw)
            nftw = int(tw // SoF_tw)
            nfbf = int(bf // SoF_bf)
            nftf = int(tf // SoF_tf)
            
            #nfdw, nftw, nfbf, nftf = 2, 2, 2, 2

            JColumn = col_sec["J"]
            GJ = Gs * JColumn
            A = col_sec["A"]
            I = col_sec["I"]

            if self.Type == "Elastic":
                ops.section('Elastic', ColSecTag, self.Es, A, I)
            elif self.Type == "Inelastic":
                RHsection(ColSecTag, matIDSteel_Col, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ)

        self.Gird_sections = [
            {"tag": 160, "d": 36.0 * 0.0254, "bf": 12.0 * 0.0254, "tw": 0.65 * 0.0254, "tf": 1.02 * 0.0254, "J": 12.4 * 4.16e-7, "A": 47 * 0.0254 ** 2, "I": 9760 * 0.0254 ** 4},
            {"tag": 135, "d": 35.6 * 0.0254, "bf": 12.0 * 0.0254, "tw": 0.6 * 0.0254, "tf": 0.79 * 0.0254, "J": 7.0 * 4.16e-7, "A": 39.9 * 0.0254 ** 2, "I": 7800 * 0.0254 ** 4},
            {"tag": 99, "d": 29.7 * 0.0254, "bf": 10.5 * 0.0254, "tw": 0.52 * 0.0254, "tf": 0.67 * 0.0254, "J": 3.77 * 4.16e-7, "A": 29 * 0.0254 ** 2, "I": 3990 * 0.0254 ** 4},
            {"tag": 84, "d": 26.7 * 0.0254, "bf": 10.0 * 0.0254, "tw": 0.46 * 0.0254, "tf": 0.64 * 0.0254, "J": 2.81 * 4.16e-7, "A": 24.7 * 0.0254 ** 2, "I": 2850 * 0.0254 ** 4},
            {"tag": 68, "d": 23.7 * 0.0254, "bf": 8.97 * 0.0254, "tw": 0.415 * 0.0254, "tf": 0.585 * 0.0254, "J": 1.87 * 4.16e-7, "A": 20.1 * 0.0254 ** 2, "I": 1830 * 0.0254 ** 4},
        ]

        for gird_sec in self.Gird_sections:
            GirdSecTag = gird_sec["tag"]
            d = gird_sec["d"]
            bf = gird_sec["bf"]
            tw = gird_sec["tw"]
            tf = gird_sec["tf"]
            dw = d - 2 * tf

            nfdw = int(dw // SoF_dw)
            nftw = int(tw // SoF_tw)
            nfbf = int(bf // SoF_bf)
            nftf = int(tf // SoF_tf)

            JGird = gird_sec["J"]
            GJ = Gs * JGird
            A = col_sec["A"]
            I = col_sec["I"]
            
            if self.Type == "Elastic":
                ops.section('Elastic', GirdSecTag, self.Es, A, I)
            elif self.Type == "Inelastic":
                RHsection(GirdSecTag, matIDSteel_Gird, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ)

    def add_elements(self):
        elements = {
            "Exterior Column": {
                "tags": [1121, 1626, 2131, 2636, 3141, 3646, 4151, 4656, 5161, 5666, 6171, 6676, 7181, 7686, 8191, 8696, 91101, 96106],
                "nodes": [(11, 21), (16, 26), (21, 31), (26, 36), (31, 41), (36, 46), (41, 51), (46, 56), (51, 61), (56, 66), (61, 71),
                          (66, 76), (71, 81), (76, 86), (81, 91), (86, 96), (91, 101), (96, 106)],
                "sections": [370] * 8 + [283] * 4 + [257] * 4 + [233] * 2,
                "transf_id": self.IDColTransf
            },

            "Interior Column": {
                "tags": [1222, 1323, 1424, 1525, 2232, 2333, 2434, 2535, 3242, 3343, 3444, 3545, 4252, 4353, 4454, 4555, 5262, 5363,
                          5464, 5565, 6272, 6373, 6474, 6575, 7282, 7383, 7484, 7585, 8292, 8393, 8494, 8595, 92102, 93103, 94104, 95105],
                "nodes": [(12, 22), (13, 23), (14, 24), (15, 25), (22, 32), (23, 33), (24, 34), (25, 35), (32, 42), (33, 43), (34, 44), 
                          (35, 45), (42, 52), (43, 53), (44, 54), (45, 55), (52, 62), (53, 63), (54, 64), (55, 65), (62, 72), (63, 73),
                          (64, 74), (65, 75), (72, 82), (73, 83), (74, 84), (75, 85), (82, 92), (83, 93), (84, 94), (85, 95), (92, 102), (93, 103), (94, 104), (95, 105)],
                "sections": [500] * 8 + [455] * 8 + [370] * 8 + [283] * 8 + [257] * 4,
                "transf_id": self.IDColTransf
            },

            "Girder": {
                "tags": [2122, 2223, 2324, 2425, 2526, 3132, 3233, 3334, 3435, 3536, 4142, 4243, 4344, 4445, 4546, 5152, 5253, 5354, 5455,
                          5556, 6162, 6263, 6364, 6465, 6566, 7172, 7273, 7374, 7475, 7576, 8182, 8283, 8384, 8485, 8586, 9192, 9293, 9394,
                            9495, 9596, 101102, 102103, 103104, 104105, 105106],
                "nodes": [(21, 22), (22, 23), (23, 24), (24, 25), (25, 26), (31, 32), (32, 33), (33, 34), (34, 35), (35, 36), (41, 42), (42, 43), 
                          (43, 44), (44, 45), (45, 46), (51, 52), (52, 53), (53, 54), (54, 55), (55, 56), (61, 62), (62, 63), (63, 64), (64, 65),
                            (65, 66), (71, 72), (72, 73), (73, 74), (74, 75), (75, 76), (81, 82), (82, 83), (83, 84), (84, 85), (85, 86), (91, 92),
                              (92, 93), (93, 94), (94, 95), (95, 96), (101, 102), (102, 103), (103, 104), (104, 105), (105, 106)],
                "sections": [160] * 10 + [135] * 20 + [99] * 5 + [84] * 5 + [68] * 5,
                "transf_id": self.IDGirdTransf
            }
        }

        if self.Type == "Elastic":
            for element_data in elements.values():
                for tag, node_pair, section_tag in zip(element_data["tags"], element_data["nodes"], element_data["sections"]):
                    if section_tag in [sec["tag"] for sec in self.Col_sections]:
                        A = next(sec["A"] for sec in self.Col_sections if sec["tag"] == section_tag)
                        I = next(sec["I"] for sec in self.Col_sections if sec["tag"] == section_tag)
                    elif section_tag in [sec["tag"] for sec in self.Gird_sections]:
                        A = next(sec["A"] for sec in self.Gird_sections if sec["tag"] == section_tag)
                        I = next(sec["I"] for sec in self.Gird_sections if sec["tag"] == section_tag)
                    ops.element('elasticBeamColumn', tag, *node_pair, A, self.Es, I, element_data["transf_id"])
        elif self.Type == "Inelastic":
            for element_data in elements.values():
                for tag, node_pair, section in zip(element_data["tags"], element_data["nodes"], element_data["sections"]):
                    ops.element('nonlinearBeamColumn', tag, *node_pair, self.nps, section, element_data["transf_id"])

        #print("Elements have been added.")

    def add_leaning_column(self):

        ops.node(17, 6 * self.LGird, 0)
        ops.fix(17, 1, 1, 0)

        for i in range(2, 10):
          ops.node(i * 100 + 71, 6 * self.LGird, self.LCol_1 + self.LCol * (i - 2))
          ops.node(i * 10 + 7, 6 * self.LGird, self.LCol_1 + self.LCol * (i - 2))
        ops.node(107, 6 * self.LGird, self.LCol_1 + self.LCol * 8)

        rigidMaterial = 5
        ops.uniaxialMaterial('Elastic', rigidMaterial, 1e12)
        rotarySpring = 6
        ops.uniaxialMaterial('ElasticPP', rotarySpring, 1, 1)
        rigidA = 1

        for i in range(2, 10):
            lower_node = (i - 1) * 10 + 7
            upper_node = lower_node * 10 + 100 + 1
            ops.element('elasticBeamColumn', upper_node + lower_node * 1000, lower_node, upper_node, rigidA, 1e12, 1e-4, self.IDColTransf)
            ops.element('zeroLength',  upper_node * 100 + ((upper_node - 1)//10), upper_node, (upper_node - 1)//10, '-mat', rotarySpring, '-dir', 3)
            ops.equalDOF(upper_node, (upper_node - 1)//10, 1, 2)
        ops.element('elasticBeamColumn', 97107, 97, 107, rigidA, 1e12, 1e-4, self.IDColTransf)

        for i in range(2, 11):
            lower_node = i * 10 + 6
            upper_node = i * 10 + 7
            ops.element('Truss', lower_node * 100 + upper_node, lower_node, upper_node, rigidA, rigidMaterial)

        leaning_column_Xfactor = 0.0
        leaning_column_Yfactor = 1.0
        leaning_column_Zfactor = 0.0
        ops.mass(27, self.Floor_Mass[0] * leaning_column_Xfactor, self.Floor_Mass[0] * leaning_column_Yfactor, self.Floor_Mass[0] * leaning_column_Zfactor)
        for i in range(3, 10):
            ops.mass(i * 10 + 7, self.Floor_Mass[1] * leaning_column_Xfactor, self.Floor_Mass[1] * leaning_column_Yfactor, self.Floor_Mass[1] * leaning_column_Zfactor)
        ops.mass(107, self.Floor_Mass[2] * leaning_column_Xfactor, self.Floor_Mass[2] * leaning_column_Yfactor, self.Floor_Mass[2] * leaning_column_Zfactor)

    def assign_element_mass(self):
        Ext_Col_Weight, Int_Col_Weight, Girder_Weight  = {}, {}, {}

        for col_sec in self.Col_sections:
            tag = col_sec["tag"] 
            d = col_sec["d"]                                      # m
            bf = col_sec["bf"]                                    # m
            tw = col_sec["tw"]                                    # m
            tf = col_sec["tf"]                                    # m
            Area = ((d - 2 * tf) * tw) + 2 * (bf * tf)            # m^2
            unit_length_weight = Area * 76872.25                  # N/m  (unit volumn of steel : 76872.25 N/m3)         
            Col_Weight = unit_length_weight * 13 * 0.3048 / 9.81  # N    (1 ft = 0.3048 m)    
            Ext_Col_Weight[str(tag)] = Col_Weight
            Int_Col_Weight[str(tag)] = Col_Weight

        for gird_sec in self.Gird_sections:
            tag = gird_sec["tag"]
            d = gird_sec["d"]                                      # m
            bf = gird_sec["bf"]                                    # m
            tw = gird_sec["tw"]                                    # m
            tf = gird_sec["tf"]                                    # m^2
            Area = ((d - 2 * tf) * tw) + 2 * (bf * tf)             # N/m  (unit volumn of steel : 76872.25 N/m3)  
            unit_length_weight = Area * 76872.25                   # N    (1 ft = 0.3048 m)
            Gird_Weight = unit_length_weight * 30 * 0.3048 / 9.81  # N    (1 ft = 0.3048 m)    
            Girder_Weight[str(tag)] = Gird_Weight 

        Ext_Col_Weight_1F = {key: value * (18/13) for key, value in Ext_Col_Weight.items()}
        Int_Col_Weight_1F = {key: value * (18/13) for key, value in Int_Col_Weight.items()}

        floor_data = [
            (2, "370", "370", "500","500","160"),
            (3, "370", "370", "500","455","160"),
            (4, "370", "370", "455","455","135"),
            (5, "370", "283", "455","370","135"),
            (6, "283", "283", "370","370","135"),
            (7, "283", "257", "370","283","135"),
            (8, "257", "257", "283","283","99"),
            (9, "257", "233", "283","257","84"),
            (10, "233", "233", "257","257", "68"),   
        ]

        node_mass_1F, node_mass_middle_floor_ext, node_mass_middle_floor_int, node_mass_9F = [], [], [], []

        for floor, ext_col_tag_below, ext_col_tag_above, int_col_tag_below, int_col_tag_above, gird_tag in floor_data:
            for i in range(1, self.Span_num + 2):
                node_tag = i + floor * 10

                if floor == 2:  # 1F

                    if i == 1 or i == self.Span_num + 1:  # 外柱
                        node_mass = Ext_Col_Weight_1F[ext_col_tag_below] / 2 + Ext_Col_Weight[ext_col_tag_above] / 2 + Girder_Weight[gird_tag] / 2
                        ops.mass(node_tag, node_mass, node_mass, node_mass)
                        node_mass_1F.append(node_mass)

                    else:  # 内柱
                        node_mass = Int_Col_Weight_1F[int_col_tag_below] / 2 + Int_Col_Weight[int_col_tag_above] / 2 + Girder_Weight[gird_tag]
                        ops.mass(node_tag, node_mass, node_mass, node_mass)
                        node_mass_1F.append(node_mass)

                elif floor == 10:  # 9F
                    if i == 1 or i == self.Span_num + 1:  # 外柱
                        node_mass = Ext_Col_Weight[ext_col_tag_below] / 2 + Girder_Weight[gird_tag] / 2
                        ops.mass(node_tag, node_mass, node_mass, node_mass)
                        node_mass_9F.append(node_mass)
                    else:  # 内柱
                        node_mass = Int_Col_Weight[int_col_tag_below] / 2 + Girder_Weight[gird_tag]
                        ops.mass(node_tag, node_mass, node_mass, node_mass)    
                        node_mass_9F.append(node_mass)

                else:  # 3F~8F
                    if i == 1 or i == self.Span_num + 1:  # 外柱
                        node_mass = Ext_Col_Weight[ext_col_tag_below] / 2 + Ext_Col_Weight[ext_col_tag_above] / 2 + Girder_Weight[gird_tag] / 2
                        ops.mass(node_tag, node_mass, node_mass, node_mass)
                        node_mass_middle_floor_ext.append(node_mass)
                    else:  # 内柱
                        node_mass = Int_Col_Weight[int_col_tag_below] / 2 + Int_Col_Weight[int_col_tag_above] / 2 + Girder_Weight[gird_tag]
                        ops.mass(node_tag, node_mass, node_mass, node_mass)
                        node_mass_middle_floor_int.append(node_mass)
        
    def assign_node_mass(self):
        Xfactor = 1.25
        Yfactor = 1.25
        Zfactor = 1.25
        for i in range(2, self.Story_num + 2):
            node_id = 1000 + i
            if i == 2:
                node_mass_x = self.Floor_Mass[0] * Xfactor
                node_mass_y = self.Floor_Mass[0] * Yfactor
                node_mass_z = self.Floor_Mass[0] * Zfactor

            elif i > 2 and i < 10:
                node_mass_x = self.Floor_Mass[1] * Xfactor
                node_mass_y = self.Floor_Mass[1] * Yfactor
                node_mass_z = self.Floor_Mass[1] * Zfactor
            elif i == 10:
                node_mass_x = self.Floor_Mass[2] * Xfactor
                node_mass_y = self.Floor_Mass[2] * Yfactor
                node_mass_z = self.Floor_Mass[2] * Zfactor
            ops.mass(node_id, node_mass_x, node_mass_y, node_mass_z)
            
        #print("Mass have been added.")
        
    def calculate_natural_frequencies(self):
        self.num_modes = 9
        eigenValues = ops.eigen('-fullGenLapack', self.num_modes)
        self.natural_frequencies = [np.sqrt(fre) for fre in eigenValues]
        self.frequencies = [omega / (2 * np.pi) for omega in self.natural_frequencies]
        print()
        print("Natural Circular frequencies (rad/s):"," ".join(f"{round(na_freq, 3)} (rad/s)" for na_freq in self.natural_frequencies))
        print("Natural Cyclic frequencies (Hz):"," ".join(f"{round(freq, 3)} Hz" for freq in self.frequencies))
    
    def calculate_mode_shape(self):
        modeshape_value = []
        mode_shapes = {}
        each_node = [i * 10 + 1 for i in range(1, self.num_modes + 2)]
        for mode in range(1, self.num_modes + 1):
            mode_shapes[mode] = {}
            for node in each_node:
                mode_shapes[mode][node] = ops.nodeEigenvector(node, mode)
                
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))

        for mode in range(1, self.num_modes + 1):
            ax = axes[(mode-1)//3, (mode-1)%3]

            x_base = [ops.nodeCoord(i)[0] for i in range(11, 111, 10)]
            y_base = [i for i in range(self.num_modes + 1)]
            
            mode_x_disp = [mode_shapes[mode][node][0] for node in each_node]
            max_disp = max(np.abs(mode_x_disp))
            mode_x_disp = [disp / max_disp for disp in mode_x_disp] 

            ax.plot(x_base, y_base, 'b--', label="Original Structure", marker='o')
            deformed_x = [x + dx for x, dx in zip(x_base, mode_x_disp)]
            ax.plot(deformed_x, y_base, 'r-', label=f"Mode {mode} Shape", marker='o')

            ax.set_title(f'Mode {mode} Shape (Freq: {self.frequencies[mode - 1]:.3f} Hz)')
            ax.set_ylabel('Floor Number')
            ax.grid(True)
            ax.set_xlim(-1.05, 1.05)
            ax.set_yticks(np.arange(0, self.num_modes + 1, 1))
            ax.set_xticks(np.arange(-1.0, 1.1, 0.2))
            modeshape_value.append(deformed_x)
        
        modeshape_data_1_to_5 = np.column_stack([modeshape_value[i] for i in range(5)])
        modeshape_data_6_to_9 = np.column_stack([modeshape_value[i] for i in range(5,9)])
        file_name = 'mode shape'
        file_path = os.path.join(os.getcwd(), 'Opensees Model', file_name) 
        os.makedirs(file_path, exist_ok=True)
        full_path_1_to_5 = os.path.join(file_path, 'modeshape_data_1_to_5.txt')
        full_path_6_to_9 = os.path.join(file_path, 'modeshape_data_6_to_9.txt')
        np.savetxt(full_path_1_to_5, modeshape_data_1_to_5, delimiter='\t')
        np.savetxt(full_path_6_to_9, modeshape_data_6_to_9, delimiter='\t')
        
        plt.tight_layout()
        mode_shape_img_path = os.path.join(file_path, 'Mode Shape.png')
        plt.savefig(mode_shape_img_path)
        #fig.subplots_adjust(hspace=1.2)
              
    def visualize_model(self):
        opsv.plot_model(node_labels=1, element_labels=0, node_supports=True, gauss_points=False)
        plt.title('9 Story Steel Structure', fontsize = 18)
        plt.ylim(ops.nodeCoord(11)[1] - 1.8, ops.nodeCoord(101)[1] + 2.3)
        #plt.show()
    
    def add_rayleigh_damping(self, mode1, mode2):
        omega1 = self.natural_frequencies[mode1 - 1]
        omega2 = self.natural_frequencies[mode2 - 1]

        betaK = 2 * self.damping_ratio / (omega1 + omega2)
        alphaM = betaK * omega1 * omega2
        betaKinit = 0
        betaKcomm = 4 * self.damping_ratio / (omega1 + omega2)
        ops.rayleigh(alphaM, betaK, betaKinit, betaKcomm)
    
    def plot_stress_strain_curve(self, Target_elements, section_num):
        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        fig.suptitle(f'Hysteresis Curves for Elements (Section {section_num})', fontsize=16)

        for idx, element_tag in enumerate(Target_elements):

            file_path = os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt')
            stress_strain_data = np.loadtxt(file_path)
        
            strain = stress_strain_data[:, 0]
            stress = stress_strain_data[:, 1]
        
            row = idx // 3
            col = idx % 3
            axs[row, col].plot(strain, stress, label=f'Element {element_tag}')
            axs[row, col].set_xlabel('Strain')
            axs[row, col].set_ylabel('Stress')
            axs[row, col].set_title(f'Element {element_tag}')
            axs[row, col].grid(True)
    
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        img_path = os.path.join(self.hysteresis_loop_folder, 'Stress_Strain_Hysteresis_Curves.png')
        plt.savefig(img_path)
        #plt.show()

    def dynamic_analysis(self, data, DataType, add_leaning_column, Target_elements, section_num, scale_factor, record_stress_strain):
        time = data[:, 0]
        acc = data[:, 1]
        
        self.acc = acc if self.Earthquake == 'Groundmotion' else (acc * 9.81 / np.max(np.abs(acc * 9.81)))
        
        dt = time[1] - time[0]
        self.DataType = DataType    
        ops.timeSeries('Path', 1, '-dt', dt, '-values', *self.acc, '-factor', scale_factor)
        ops.pattern('MultipleSupport',1)
        ops.groundMotion(1,'Plain','-accel',1)

        value = 2 if add_leaning_column else 1
            
        for Base_node in range(11, 11 + self.Span_num + value):
            ops.imposedMotion(Base_node,1,1)

        self.main_folder_name = os.path.join(os.getcwd(),"Opensees Model", f"{self.Type} Opensees Structure Response", self.Earthquake)
        if not os.path.exists(self.main_folder_name):
            os.makedirs(self.main_folder_name)

        self.abs_folder_name = os.path.join(self.main_folder_name, f'abs_{self.DataType}')

        if os.path.exists(self.abs_folder_name):
             shutil.rmtree(self.abs_folder_name)
        os.makedirs(self.abs_folder_name)

        for i in range(1, self.Story_num + 2):
            ops.recorder('Node', '-file', os.path.join(self.abs_folder_name, f'node{i * 10 + 1}{self.DataType}.txt'), '-time', '-node', i * 10 + 1, '-dof', 1, self.DataType)
        
        if record_stress_strain:
            if self.Type == 'Inelastic':
                self.hysteresis_loop_folder = os.path.join(self.strain_stress_folder_name, f'{self.Earthquake}') 
                if os.path.exists(self.hysteresis_loop_folder):
                    shutil.rmtree(self.hysteresis_loop_folder)
                os.makedirs(self.hysteresis_loop_folder)

                for element_tag in Target_elements:
                    ops.recorder('Element', '-file', os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt'), 
                                 '-ele', element_tag, 'section', section_num, 'fiber', 1, 1, 'stressStrain')
            
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        tolerance = 1 if scale_factor > 9.81 else 1e-1
        ops.test('EnergyIncr', tolerance, 100)
        ops.algorithm('Newton')
        #ops.algorithm('NewtonLineSearch')
        ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')
        result = ops.analyze(len(self.acc), dt)
        
        if record_stress_strain:
            if self.Type == 'Inelastic':
                self.plot_stress_strain_curve(Target_elements, section_num)

    def plot_dynamic_result(self, data, add_leaning_column):
        plt.figure()
        plt.plot(data[:, 0], self.acc)
        plt.xlabel('Time (sec)', fontsize = 14)
        plt.ylabel('Acceleration (m/s2)', fontsize = 14)
        plt.xlim(data[:,0][0], data[:,0][-1])
        plt.title(self.Earthquake, fontsize = 18)
        plt.grid(True)

        self.rela_folder_name = os.path.join(self.main_folder_name, f'rela_{self.DataType}')
        if not os.path.exists(self.rela_folder_name):
            os.makedirs(self.rela_folder_name)

        Ground_data = np.loadtxt(os.path.join(self.abs_folder_name, f'node{11}{self.DataType}.txt'))
        ground_data = Ground_data[:,1]

        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        fig.suptitle(f'{self.Earthquake} Earthquake {self.plot_title} Response', fontsize=12)

        for i in range(2, Story_num + 2):
            data = np.loadtxt(os.path.join(self.abs_folder_name, f'node{i * 10 + 1}{self.DataType}.txt'))
            relative_data = data[:, 1] - ground_data
            time_data = np.insert(data[:, 0], 0, 0)
            relative_data = np.insert(relative_data, 0, 0) 
            np.savetxt(os.path.join(self.rela_folder_name, f'node{i * 10 + 1}_relative_{self.DataType}.txt'), np.column_stack((time_data, relative_data)), delimiter='\t')
            row = (i - 2) // 3  
            col = (i - 2) % 3   
    
            axs[row, col].plot(time_data, relative_data)
            axs[row, col].set_xlabel('Time', fontsize=10)
            axs[row, col].set_ylabel(f'Displacement (m)' if self.DataType == 'disp' else f'Acceleration (m/s²)' if self.DataType == 'accel' else f'Velocity (m/s)', fontsize = 10)
            axs[row, col].set_xlim(time_data[0], time_data[-1]) 
            axs[row, col].set_title(f'{self.DataType} of Node {i * 10 + 1} Over Time', fontsize=12)
            axs[row, col].grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(self.rela_folder_name, f'rela_{self.DataType}.png'))
        #plt.show()
        
        add = 'with' if add_leaning_column else 'no'
        print(f"Complete the generation of {self.plot_title} data for {self.Earthquake} Earthquake ({self.Type} {add} leaning column)")

    def calculate_floor_stiffness_method(self,push_floor, push_force):

        ops.wipe() 
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        ops.geomTransf('PDelta', self.IDColTransf)
        ops.geomTransf('Linear', self.IDGirdTransf)
        self.build_nodes()  
        self.setup_perpendicular_direction() 
        self.define_materials_and_sections(self.Steel_Type)  
        self.add_elements()  
        self.assign_node_mass() 

        push_node = 1000 + push_floor

        ops.wipeAnalysis()  
        ops.remove('loadPattern', 1)  
        ops.remove('timeSeries', 1)   

        for i in range(2, self.Story_num + 2):
            rigid_nodeID = 1000 + i
            ops.remove('sp', rigid_nodeID, 1)
            ops.remove('sp', rigid_nodeID, 2)
            ops.remove('sp', rigid_nodeID, 3)

            if i < push_floor:
                   ops.fix(rigid_nodeID, 1, 1, 1)
            else:
                   ops.fix(rigid_nodeID, 0, 1, 1)

        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
    
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1,1)
        ops.load(push_node, push_force, 0.0, 0.0)
        
        analyze_result = ops.analyze(1)
        
        top_disp = ops.nodeDisp(push_node, 1)
        
        if top_disp != 0:
            stiffness = push_force / top_disp
        else:
            stiffness = float('inf')
        
        print(f"Floor {push_floor - 1} (Node {push_node}) Displacement: {round(top_disp * 100, 6)} cm, Lateral Stiffness: {round(stiffness / 100000, 3)} kN/cm")
    
    def perform_cyclic_pushover_analysis(self, Target_elements, max_displacement, cycles):
        hysteresis_data = []
        element_hysteresis_data = []

        for floor, element in enumerate(Target_elements, start=2):
            ops.wipe() 
            ops.model('basic', '-ndm', 2, '-ndf', 3)
            ops.geomTransf('PDelta', self.IDColTransf)
            ops.geomTransf('Linear', self.IDGirdTransf)
            self.build_nodes()
            self.setup_perpendicular_direction() 
            self.define_materials_and_sections(Steel_Type)
            self.add_elements()
            self.assign_node_mass()

            for i in range(2, self.Story_num + 2):
                nodeID = 1000 + floor
                ops.remove('sp', nodeID, 1)
                ops.remove('sp', nodeID, 2)
                ops.remove('sp', nodeID, 3)
                if i < floor:
                    ops.fix(nodeID, 1, 1, 1)
                else:
                    ops.fix(nodeID, 0, 1, 1)

            ops.wipeAnalysis()
            ops.constraints('Plain')
            ops.numberer('Plain')
            ops.system('BandGeneral')
            ops.test('EnergyIncr', 1e-1, 10)
            ops.algorithm('Newton')

            push_node = 1000 + floor

            amplitudes = np.linspace(0, max_displacement, cycles + 1)[1:]  
            displacement_history = []
            for amp in amplitudes:
                displacement_history.extend([amp, -amp])

            ops.timeSeries('Linear', 1)
            ops.pattern('Plain', 1, 1)
            ops.load(push_node, 1.0, 0.0, 0.0)  

            recorder_folder = os.path.join('Opensees Model', 'hysteresis_data', f'floor_{floor}')
            if not os.path.exists(recorder_folder):
                os.makedirs(recorder_folder)

            element_recorder = os.path.join(recorder_folder, f'element_{element}_stressStrain.txt')
            ops.recorder('Element', '-file', element_recorder, '-time', '-ele', element, 'section', 1, 'fiber', 1, 1, 'stressStrain')

            ops.integrator('DisplacementControl', push_node, 1, 0.0001) 
            ops.analysis('Static')

            ok = 0
            for target_disp in displacement_history:
                current_disp = ops.nodeDisp(push_node, 1)
                disp_increment = target_disp - current_disp
                num_steps = int(abs(disp_increment) / 0.0001)
                if num_steps == 0:
                     num_steps = 1
                disp_incr_per_step = disp_increment / num_steps
                ops.integrator('DisplacementControl', push_node, 1, disp_incr_per_step)
                for step in range(num_steps):
                    ok = ops.analyze(1)
                    if ok != 0:
                        print(f"Analysis failed on floor {floor} at target displacement {target_disp}")
                        break
                if ok != 0:
                    break

            if ok == 0:
                print(f"Cyclic pushover analysis successfully completed on floor {floor}")
            else:
                print(f"Cyclic pushover analysis failed on floor {floor}")

            element_data = np.loadtxt(element_recorder)
            element_hysteresis_data.append((floor, element_data))

            ops.remove('recorders')
            ops.wipe()

        self.plot_element_hysteresis_curves(element_hysteresis_data)  
    
    def plot_element_hysteresis_curves(self, element_hysteresis_data):
        num_floors = len(element_hysteresis_data)
        num_cols = 3  
        num_rows = (num_floors + num_cols - 1) // num_cols  

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols * 5, num_rows * 4))
        axs = axs.flatten()  

        for idx, (floor, element_data) in enumerate(element_hysteresis_data):
            strain = element_data[:, 1]  
            stress = element_data[:, 0] 

            axs[idx].plot(strain, stress)
            axs[idx].set_xlabel('Strain')
            axs[idx].set_ylabel('Stress (Pa)')
            axs[idx].set_title(f'Stress-Strain Curve for Floor {floor}')
            axs[idx].grid(True)

        plt.tight_layout()
        plt.show()
    
if __name__ == "__main__":

    start_time = time.time()
    
    # Elastic or Inelastic
    Type = "Inelastic"
    
    # Steel Type 
    Steel_Type = 'Steel01' # Steel01 or Steel02
    
    # Earthquake Type
    Earthquake_dict = {
    0: "Groundmotion.txt",
    1: "El Centro.txt",
    2: "Kobe.txt",
    3: "Chi-Chi.txt",
    4: "Northridge.txt",
    }
    scale_factor = 9.81 * 0.1
    Earthquake = Earthquake_dict[4]
    Earthquake_Data = np.loadtxt(os.path.join(os.getcwd(),"Earthquake", Earthquake))
    
    # Response Type
    Analysis_Type = 'accel'  # 'disp' 'vel' 'accel' 'rayleighForces'
    
    # Leaning Column ?
    add_leaning_column = True  # True or False
    
    # Parameter
    lCol = 3.9624    # m
    lCol1 = 5.4864   # m
    lGird = 9.144    # m
    Story_num = 9       
    Span_num = 5
    nps = 5
    #Floor_Dead_Load = 96 * 150 * 150 * 0.453 / 2    # kg
    #Floor_Mass = 86 * 150 * 150 * 0.453 / 2         # kg
    Mass_2F = (1007304 / 2)                          # kg
    Mass_3_to_9F = (990088 / 2)                      # kg
    Mass_Roof = (1066540 / 2)                        # kg
    Floor_Mass = [Mass_2F, Mass_3_to_9F, Mass_Roof]  # kg
    Damping_ratio = 0.02

    section_num = 1  # nps section
    Target_elements = [1121, 2131, 3141, 4151, 5161, 6171, 7181, 8191, 91101]
    #[1424, 2434, 3444, 4454, 5464, 6474, 7484, 8494, 94104]

    model = OpenseesModel(Type, Analysis_Type, lCol, lCol1, lGird, Story_num, Span_num, nps, Floor_Mass, Damping_ratio, Earthquake.replace(".txt", ""))
    model.build_nodes()
    model.setup_perpendicular_direction()
    model.define_materials_and_sections(Steel_Type)
    model.add_elements()
    if add_leaning_column:
        model.add_leaning_column()
    model.assign_element_mass()
    model.assign_node_mass()
    model.calculate_natural_frequencies()
    model.calculate_mode_shape()
    model.visualize_model()
    model.add_rayleigh_damping(1, 3)
    model.dynamic_analysis(Earthquake_Data, Analysis_Type, add_leaning_column, Target_elements,
                           section_num, scale_factor, record_stress_strain = False)
    model.plot_dynamic_result(Earthquake_Data, add_leaning_column)
    
    Calculate_Lateral_Stiffness = False  # True or False
    
    if Calculate_Lateral_Stiffness:
        push_force = 10000 # N
        for floor in range(2, Story_num + 2):
            model.calculate_floor_stiffness_method(floor, push_force = 10000)
            
    #model.perform_cyclic_pushover_analysis(Target_elements, max_displacement=0.05, cycles=15)
    plt.show()
    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)  
    minutes, seconds = divmod(remainder, 60)
    print(f"Total runtime: {int(hours)} hours {int(minutes)} minutes {int(seconds)} seconds")