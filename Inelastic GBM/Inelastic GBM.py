import openseespy.opensees as ops
import opsvis as opsv
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import shutil

class Inelastic_GBM:
    def __init__(self, Analysis_Type, Steel_Type, LCol, LCol_1, LGird, Story_num, nps, Floor_Mass, damping_ratio, 
                 shear_Area_factor, flexural_Area_factor, Earthquake,  yielding_factor):
        self.Analysis_Type = Analysis_Type
        self.Steel_Type = Steel_Type
        self.LCol = LCol
        self.LCol_1 = LCol_1
        self.LGird = LGird
        self.Story_num = Story_num
        self.Node = []
        self.nps = nps
        self.Floor_Mass = Floor_Mass
        self.IDColTransf = 1
        self.IDGirdTransf = 2
        self.IDLColTransf = 3
        self.story_disp = []
        self.damping_ratio = damping_ratio
        self.Earthquake = Earthquake
        self.shear_Area_factor = shear_Area_factor
        self.flexural_Area_factor = flexural_Area_factor
        self.yielding_factor =  yielding_factor
        self.plot_title = ('Displacement' if self.Analysis_Type == 'disp' else 'Acceleration' if self.Analysis_Type == 'accel'
                           else 'Velocity' if self.Analysis_Type == 'vel' else 'Unknown Data Type')
        self.strain_stress_folder_name = os.path.join(os.getcwd(), 'Inelastic GBM', 'Inelastic Hysteresis Loop')
        os.makedirs(self.strain_stress_folder_name, exist_ok=True)
        self.initialize_model()

    def initialize_model(self):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        ops.geomTransf('Linear', self.IDColTransf)
        ops.geomTransf('Linear', self.IDGirdTransf)
        ops.geomTransf('PDelta', self.IDLColTransf)

    def create_node(self, i, j, x, y):
        nodeID = 10 * i + j
        ops.node(nodeID, x, y)
        return nodeID
    
    def build_nodes(self):
        # Base node
        self.Node = [self.create_node(1, j, (j - 1) * self.LGird, 0) for j in range(1, 3)]
        [ops.fix(nodeID, 1, 1, 1) for nodeID in self.Node]
        
        # Stick node
        for i in range(2, self.Story_num + 2):
            for j in range(1, 3):
                nodeID = self.create_node(i, j, (j - 1) * self.LGird, (i - 2) * self.LCol + self.LCol_1)
                if j == 1:
                    ops.fix(nodeID, 0, 1, 1)   # Shear node
                elif j == 2:
                    ops.fix(nodeID, 0, 1, 0)   # Flexural node

    def define_materials_and_sections(self):
        Fy = 50.0 * 4448.22 * 39.37 * 39.37       # N/m^2
        self.Es = 2040 * 9800 * 10000 * 0.815     # N/m^2
        b = 0.02
        nu = 0.3  
        Gs = self.Es / (2 * (1 + nu))
        R0, cR1, cR2 = 18.0, 0.8, 0.15
        a1, a2, a3, a4, sigInit = 0.0, 1.0, 0.0, 1.0, 0.0

        matIDSteel_list = [11, 22, 33, 44, 55, 66, 77, 88, 99]

        if self.Steel_Type == 'Steel01':
            for i, matID in enumerate(matIDSteel_list):
                ops.uniaxialMaterial('Steel01', matID, Fy * self.yielding_factor[i], self.Es, b)

        elif self.Steel_Type == 'Steel02':
            for i, matID in enumerate(matIDSteel_list):
                ops.uniaxialMaterial('Steel02', matID, Fy * self.yielding_factor[i], self.Es, b, R0, cR1, cR2, a1, a2, a2, a4, sigInit)


        def RHsection(secID, matID, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ):
            ops.section('Fiber', secID, '-GJ', GJ)
            ops.patch('rect', matID, nftf, nfbf, -0.5 * d, -0.5 * bf, -0.5 * d + tf, 0.5 * bf)
            ops.patch('rect', matID, nfdw, nftw, -0.5 * d + tf, -0.5 * tw, 0.5 * d - tf, 0.5 * tw)
            ops.patch('rect', matID, nftf, nfbf, 0.5 * d - tf, -0.5 * bf, 0.5 * d, 0.5 * bf)
            
        SoF_dw = 0.01     # Fiber Spacing dw (m)
        SoF_tw = 0.0025   # Fiber Spacing tw (m)
        SoF_bf = 0.0175   # Fiber Spacing bf (m)
        SoF_tf = 0.00275  # Fiber Spacing tf (m)
        
        self.shear_Col_sections = [
            {"tag": 3701, "matID": matIDSteel_list[0], "d": 17.9 * 0.0254 * self.shear_Area_factor[0], "bf": 16.5 * 0.0254 * self.shear_Area_factor[0], "tw": 1.66 * 0.0254 * self.shear_Area_factor[0], "tf": 2.66 * 0.0254 * self.shear_Area_factor[0], "J": 222.0 * 4.16e-7 * (self.shear_Area_factor[0] ** 4), "A": 109 * (0.0254 ** 2) * (self.shear_Area_factor[0] ** 2), "I": 5440 * 0.0254 ** 4 * (self.shear_Area_factor[0] ** 4)},
            {"tag": 3702, "matID": matIDSteel_list[1], "d": 17.9 * 0.0254 * self.shear_Area_factor[1], "bf": 16.5 * 0.0254 * self.shear_Area_factor[1], "tw": 1.66 * 0.0254 * self.shear_Area_factor[1], "tf": 2.66 * 0.0254 * self.shear_Area_factor[1], "J": 222.0 * 4.16e-7 * (self.shear_Area_factor[1] ** 4), "A": 109 * (0.0254 ** 2) * (self.shear_Area_factor[1] ** 2), "I": 5440 * 0.0254 ** 4 * (self.shear_Area_factor[1] ** 4)},
            {"tag": 3703, "matID": matIDSteel_list[2], "d": 17.9 * 0.0254 * self.shear_Area_factor[2], "bf": 16.5 * 0.0254 * self.shear_Area_factor[2], "tw": 1.66 * 0.0254 * self.shear_Area_factor[2], "tf": 2.66 * 0.0254 * self.shear_Area_factor[2], "J": 222.0 * 4.16e-7 * (self.shear_Area_factor[2] ** 4), "A": 109 * (0.0254 ** 2) * (self.shear_Area_factor[2] ** 2), "I": 5440 * 0.0254 ** 4 * (self.shear_Area_factor[2] ** 4)},
            {"tag": 3704, "matID": matIDSteel_list[3], "d": 17.9 * 0.0254 * self.shear_Area_factor[3], "bf": 16.5 * 0.0254 * self.shear_Area_factor[3], "tw": 1.66 * 0.0254 * self.shear_Area_factor[3], "tf": 2.66 * 0.0254 * self.shear_Area_factor[3], "J": 222.0 * 4.16e-7 * (self.shear_Area_factor[3] ** 4), "A": 109 * (0.0254 ** 2) * (self.shear_Area_factor[3] ** 2), "I": 5440 * 0.0254 ** 4 * (self.shear_Area_factor[3] ** 4)},
            {"tag": 2835, "matID": matIDSteel_list[4], "d": 16.7 * 0.0254 * self.shear_Area_factor[4], "bf": 16.1 * 0.0254 * self.shear_Area_factor[4], "tw": 1.29 * 0.0254 * self.shear_Area_factor[4], "tf": 2.07 * 0.0254 * self.shear_Area_factor[4], "J": 104.0 * 4.16e-7 * (self.shear_Area_factor[4] ** 4), "A": 83.3 * (0.0254 ** 2) * (self.shear_Area_factor[4] ** 2), "I": 3840 * 0.0254 ** 4 * (self.shear_Area_factor[4] ** 4)},
            {"tag": 2836, "matID": matIDSteel_list[5], "d": 16.7 * 0.0254 * self.shear_Area_factor[5], "bf": 16.1 * 0.0254 * self.shear_Area_factor[5], "tw": 1.29 * 0.0254 * self.shear_Area_factor[5], "tf": 2.07 * 0.0254 * self.shear_Area_factor[5], "J": 104.0 * 4.16e-7 * (self.shear_Area_factor[5] ** 4), "A": 83.3 * (0.0254 ** 2) * (self.shear_Area_factor[5] ** 2), "I": 3840 * 0.0254 ** 4 * (self.shear_Area_factor[5] ** 4)},       
            {"tag": 2577, "matID": matIDSteel_list[6], "d": 16.4 * 0.0254 * self.shear_Area_factor[6], "bf": 16.0 * 0.0254 * self.shear_Area_factor[6], "tw": 1.18 * 0.0254 * self.shear_Area_factor[6], "tf": 1.89 * 0.0254 * self.shear_Area_factor[6], "J": 79.1 * 4.16e-7 * (self.shear_Area_factor[6] ** 4), "A": 75.6 * (0.0254 ** 2) * (self.shear_Area_factor[6] ** 2), "I": 3400 * 0.0254 ** 4 * (self.shear_Area_factor[6] ** 4)},
            {"tag": 2578, "matID": matIDSteel_list[7], "d": 16.4 * 0.0254 * self.shear_Area_factor[7], "bf": 16.0 * 0.0254 * self.shear_Area_factor[7], "tw": 1.18 * 0.0254 * self.shear_Area_factor[7], "tf": 1.89 * 0.0254 * self.shear_Area_factor[7], "J": 79.1 * 4.16e-7 * (self.shear_Area_factor[7] ** 4), "A": 75.6 * (0.0254 ** 2) * (self.shear_Area_factor[7] ** 2), "I": 3400 * 0.0254 ** 4 * (self.shear_Area_factor[7] ** 4)},
            {"tag": 2339, "matID": matIDSteel_list[8], "d": 16.0 * 0.0254 * self.shear_Area_factor[8], "bf": 15.9 * 0.0254 * self.shear_Area_factor[8], "tw": 1.07 * 0.0254 * self.shear_Area_factor[8], "tf": 1.72 * 0.0254 * self.shear_Area_factor[8], "J": 59.5 * 4.16e-7 * (self.shear_Area_factor[8] ** 4), "A": 68.5 * (0.0254 ** 2) * (self.shear_Area_factor[8] ** 2), "I": 3010 * 0.0254 ** 4 * (self.shear_Area_factor[8] ** 4)},
        ]
        
        for col_sec in self.shear_Col_sections:
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

            nfdw, nftw, nfbf, nftf = 2, 2, 2, 2
            
            JColumn = col_sec["J"]
            GJ = Gs * JColumn
            A = col_sec["A"]
            I = col_sec["I"]
            matIDSteel = col_sec["matID"]
            RHsection(ColSecTag, matIDSteel, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ)

        self.flexural_Col_sections = [
            {"tag": 37012, "matID": matIDSteel_list[0], "d": 17.9 * 0.0254 * self.flexural_Area_factor[0], "bf": 16.5 * 0.0254 * self.flexural_Area_factor[0], "tw": 1.66 * 0.0254 * self.flexural_Area_factor[0], "tf": 2.66 * 0.0254 * self.flexural_Area_factor[0], "J": 222.0 * 4.16e-7 * (self.flexural_Area_factor[0] ** 4), "A": 109 * (0.0254 ** 2) * (self.flexural_Area_factor[0] ** 2), "I": 5440 * 0.0254 ** 4 * (self.flexural_Area_factor[0] ** 4)},
            {"tag": 37022, "matID": matIDSteel_list[1], "d": 17.9 * 0.0254 * self.flexural_Area_factor[1], "bf": 16.5 * 0.0254 * self.flexural_Area_factor[1], "tw": 1.66 * 0.0254 * self.flexural_Area_factor[1], "tf": 2.66 * 0.0254 * self.flexural_Area_factor[1], "J": 222.0 * 4.16e-7 * (self.flexural_Area_factor[1] ** 4), "A": 109 * (0.0254 ** 2) * (self.flexural_Area_factor[1] ** 2), "I": 5440 * 0.0254 ** 4 * (self.flexural_Area_factor[1] ** 4)},
            {"tag": 37032, "matID": matIDSteel_list[2], "d": 17.9 * 0.0254 * self.flexural_Area_factor[2], "bf": 16.5 * 0.0254 * self.flexural_Area_factor[2], "tw": 1.66 * 0.0254 * self.flexural_Area_factor[2], "tf": 2.66 * 0.0254 * self.flexural_Area_factor[2], "J": 222.0 * 4.16e-7 * (self.flexural_Area_factor[2] ** 4), "A": 109 * (0.0254 ** 2) * (self.flexural_Area_factor[2] ** 2), "I": 5440 * 0.0254 ** 4 * (self.flexural_Area_factor[2] ** 4)},
            {"tag": 37042, "matID": matIDSteel_list[3], "d": 17.9 * 0.0254 * self.flexural_Area_factor[3], "bf": 16.5 * 0.0254 * self.flexural_Area_factor[3], "tw": 1.66 * 0.0254 * self.flexural_Area_factor[3], "tf": 2.66 * 0.0254 * self.flexural_Area_factor[3], "J": 222.0 * 4.16e-7 * (self.flexural_Area_factor[3] ** 4), "A": 109 * (0.0254 ** 2) * (self.flexural_Area_factor[3] ** 2), "I": 5440 * 0.0254 ** 4 * (self.flexural_Area_factor[3] ** 4)},
            {"tag": 28352, "matID": matIDSteel_list[4], "d": 16.7 * 0.0254 * self.flexural_Area_factor[4], "bf": 16.1 * 0.0254 * self.flexural_Area_factor[4], "tw": 1.29 * 0.0254 * self.flexural_Area_factor[4], "tf": 2.07 * 0.0254 * self.flexural_Area_factor[4], "J": 104.0 * 4.16e-7 * (self.flexural_Area_factor[4] ** 4), "A": 83.3 * (0.0254 ** 2) * (self.flexural_Area_factor[4] ** 2), "I": 3840 * 0.0254 ** 4 * (self.flexural_Area_factor[4] ** 4)},
            {"tag": 28362, "matID": matIDSteel_list[5], "d": 16.7 * 0.0254 * self.flexural_Area_factor[5], "bf": 16.1 * 0.0254 * self.flexural_Area_factor[5], "tw": 1.29 * 0.0254 * self.flexural_Area_factor[5], "tf": 2.07 * 0.0254 * self.flexural_Area_factor[5], "J": 104.0 * 4.16e-7 * (self.flexural_Area_factor[5] ** 4), "A": 83.3 * (0.0254 ** 2) * (self.flexural_Area_factor[5] ** 2), "I": 3840 * 0.0254 ** 4 * (self.flexural_Area_factor[5] ** 4)},       
            {"tag": 25772, "matID": matIDSteel_list[6], "d": 16.4 * 0.0254 * self.flexural_Area_factor[6], "bf": 16.0 * 0.0254 * self.flexural_Area_factor[6], "tw": 1.18 * 0.0254 * self.flexural_Area_factor[6], "tf": 1.89 * 0.0254 * self.flexural_Area_factor[6], "J": 79.1 * 4.16e-7 * (self.flexural_Area_factor[6] ** 4), "A": 75.6 * (0.0254 ** 2) * (self.flexural_Area_factor[6] ** 2), "I": 3400 * 0.0254 ** 4 * (self.flexural_Area_factor[6] ** 4)},
            {"tag": 25782, "matID": matIDSteel_list[7], "d": 16.4 * 0.0254 * self.flexural_Area_factor[7], "bf": 16.0 * 0.0254 * self.flexural_Area_factor[7], "tw": 1.18 * 0.0254 * self.flexural_Area_factor[7], "tf": 1.89 * 0.0254 * self.flexural_Area_factor[7], "J": 79.1 * 4.16e-7 * (self.flexural_Area_factor[7] ** 4), "A": 75.6 * (0.0254 ** 2) * (self.flexural_Area_factor[7] ** 2), "I": 3400 * 0.0254 ** 4 * (self.flexural_Area_factor[7] ** 4)},
            {"tag": 23392, "matID": matIDSteel_list[8], "d": 16.0 * 0.0254 * self.flexural_Area_factor[8], "bf": 15.9 * 0.0254 * self.flexural_Area_factor[8], "tw": 1.07 * 0.0254 * self.flexural_Area_factor[8], "tf": 1.72 * 0.0254 * self.flexural_Area_factor[8], "J": 59.5 * 4.16e-7 * (self.flexural_Area_factor[8] ** 4), "A": 68.5 * (0.0254 ** 2) * (self.flexural_Area_factor[8] ** 2), "I": 3010 * 0.0254 ** 4 * (self.flexural_Area_factor[8] ** 4)},
        ]
        for col_sec in self.flexural_Col_sections:
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

            nfdw, nftw, nfbf, nftf = 2, 2, 2, 2
            
            JColumn = col_sec["J"]
            GJ = Gs * JColumn
            A = col_sec["A"]
            I = col_sec["I"]
            matIDSteel = col_sec["matID"]
            RHsection(ColSecTag, matIDSteel, d, bf, tf, tw, nfdw, nftw, nfbf, nftf, GJ)

    def add_elements(self):
        shear_elements = {
            "Exterior Column": {
                "tags": [1121, 2131, 3141, 4151, 5161, 6171, 7181, 8191, 91101],
                "nodes": [(11, 21), (21, 31), (31, 41), (41, 51), (51, 61), (61, 71), (71, 81), (81, 91), (91, 101)],
                "sections": [3701, 3702, 3703, 3704, 2835, 2836, 2577, 2578, 2339],
                "transf_id": self.IDColTransf
            }
        }
        for element_data in shear_elements.values():
                for tag, node_pair, section in zip(element_data["tags"], element_data["nodes"], element_data["sections"]):
                    ops.element('nonlinearBeamColumn', tag, *node_pair, self.nps, section, element_data["transf_id"])
        
        flexural_elements = {
            "Exterior Column": {
                "tags": [1222, 2232, 3242, 4252, 5262, 6272, 7282, 8292, 92102],
                "nodes": [(12, 22), (22, 32), (32, 42), (42, 52), (52, 62), (62, 72), (72, 82), (82, 92), (92, 102)],
                "sections": [37012, 37022, 37032, 37042, 28352, 28362, 25772, 25782, 23392],
                "transf_id": self.IDColTransf
            }
        }
        for element_data in flexural_elements.values():
                for tag, node_pair, section in zip(element_data["tags"], element_data["nodes"], element_data["sections"]):
                    ops.element('nonlinearBeamColumn', tag, *node_pair, self.nps, section, element_data["transf_id"])


    def add_rigid_link(self):
        self.rigidMaterial = 5
        ops.uniaxialMaterial('Elastic', self.rigidMaterial, 1e12) 
        rigidA = 1     
        
        for i in range(2, self.Story_num + 2):
            left_node = i * 10 + 1
            right_node = i * 10 + 2
            ops.element('Truss', left_node * 100 + right_node, left_node, right_node, rigidA, self.rigidMaterial) 

    def assign_node_mass(self):
        for i in range(2, self.Story_num + 2):
            node_id_shear = 10 * i + 1
            node_id_flexural = 10 * i + 2
            if i == 2:
                node_mass = self.Floor_Mass[0]
            elif i > 2 and i < 10:
                node_mass = self.Floor_Mass[1]
            elif i == 10:
                node_mass = self.Floor_Mass[2]
            ops.mass(node_id_shear, node_mass / 2, 0.0, 0.0)
            ops.mass(node_id_flexural, node_mass/ 2, 0.0, 0.0)

    
    def add_leaning_column(self):
        ops.node(13, 2 * self.LGird, 0)
        ops.fix(13, 1, 1, 0)
        
        for i in range(2, 10):
            ops.node(i * 100 + 31, 2 * self.LGird, self.LCol_1 + self.LCol * (i - 2))
            ops.node(i * 10 + 3, 2 * self.LGird, self.LCol_1 + self.LCol * (i - 2))
        ops.node(103, 2 * self.LGird, self.LCol_1 + self.LCol * 8)

        rigidMaterial = 6
        ops.uniaxialMaterial('Elastic', rigidMaterial, 1e12)
        rotarySpring = 7
        ops.uniaxialMaterial('ElasticPP', rotarySpring, 1, 1)
        rigidA = 1

        for i in range(2, 10):
            lower_node = (i - 1) * 10 + 3
            upper_node = lower_node * 10 + 100 + 1
            ops.element('elasticBeamColumn', upper_node + lower_node * 1000, lower_node, upper_node, rigidA, 1e12, 1e-4, self.IDColTransf)
            ops.element('zeroLength',  upper_node * 100 + ((upper_node - 1)//10), upper_node, (upper_node - 1)//10, '-mat', rotarySpring, '-dir', 3)
            ops.equalDOF(upper_node, (upper_node - 1)//10, 1, 2)
        ops.element('elasticBeamColumn', 93103, 93, 103, rigidA, 1e12, 1e-4, self.IDColTransf)

        for i in range(2, 11):
            lower_node = i * 10 + 2
            upper_node = i * 10 + 3
            ops.element('Truss', lower_node * 100 + upper_node, lower_node, upper_node, rigidA, rigidMaterial)

        leaning_column_Xfactor = 0.0
        leaning_column_Yfactor = 1.0
        leaning_column_Zfactor = 0.0
        ops.mass(23, self.Floor_Mass[0] * leaning_column_Xfactor, self.Floor_Mass[0] * leaning_column_Yfactor, self.Floor_Mass[0] * leaning_column_Zfactor)
        for i in range(3, 10):
            ops.mass(i * 10 + 3, self.Floor_Mass[1] * leaning_column_Xfactor, self.Floor_Mass[1] * leaning_column_Yfactor, self.Floor_Mass[1] * leaning_column_Zfactor)
        ops.mass(103, self.Floor_Mass[2] * leaning_column_Xfactor, self.Floor_Mass[2] * leaning_column_Yfactor, self.Floor_Mass[2] * leaning_column_Zfactor)
    
    def calculate_natural_frequencies(self):
        self.num_modes = 9
        eigenValues = ops.eigen('-fullGenLapack', self.num_modes)
        self.natural_frequencies = [np.sqrt(fre) for fre in eigenValues]
        self.frequencies = [omega / (2 * np.pi) for omega in self.natural_frequencies]
        print()
        print("Natural Circular frequencies (rad/s):"," ".join(f"{round(na_freq, 3)} (rad/s)" for na_freq in self.natural_frequencies))
        print("Natural Cyclic frequencies (Hz):"," ".join(f"{round(freq, 3)} Hz" for freq in self.frequencies))

    def calculate_mode_shape(self, Store_data):
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

        if Store_data:
            modeshape_data_1_to_5 = np.column_stack([modeshape_value[i] for i in range(5)])
            modeshape_data_6_to_9 = np.column_stack([modeshape_value[i] for i in range(5,9)])
            file_name = 'modeshape'
            file_path = os.path.join(os.getcwd(), 'Opensees Model', file_name) 
            os.makedirs(file_path, exist_ok=True)
            full_path_1_to_5 = os.path.join(file_path, 'modeshape_data_1_to_5.txt')
            full_path_6_to_9 = os.path.join(file_path, 'modeshape_data_6_to_9.txt')
            np.savetxt(full_path_1_to_5, modeshape_data_1_to_5, delimiter='\t')
            np.savetxt(full_path_6_to_9, modeshape_data_6_to_9, delimiter='\t')
            
        #fig.subplots_adjust(hspace=1.2)
        plt.tight_layout()

    def visualize_model(self, add_leaning_column):
        opsv.plot_model(node_labels=0, element_labels=0, node_supports=False, gauss_points=False)
        fig = plt.gcf() 
        ax = plt.gca()  
        
        for i in range(1, self.Story_num + 2):
            for j in range(1, 3):
                nodeID = 10 * i + j
                coords = ops.nodeCoord(nodeID)
                if j == 1:
                    ax.text(coords[0] - 0.5, coords[1] + 0.5, str(nodeID), fontsize=10, color='blue', ha='right')
                else:
                    ax.text(coords[0] - 1.5 , coords[1] + 0.5, str(nodeID), fontsize=10, color='blue', ha='left')
        
        if add_leaning_column:
            leaning_column_nodes = [13] + [i * 100 + 31 for i in range(2, 10)] + [103]
            for nodeID in leaning_column_nodes:
                coords = ops.nodeCoord(nodeID)
                ax.text(coords[0] + 3.7, coords[1] + 0.5, str(nodeID), fontsize=10, color='red', ha='right')      
        plt.title('9 Story Enhanced GBM', fontsize=18)
        plt.ylabel('Height (m)', fontsize = 12)
        plt.xlabel('Spacing (m)', fontsize = 12)
        plt.grid(True)
        #plt.show()

    def add_rayleigh_damping(self, mode1, mode2):
        omega1 = self.natural_frequencies[mode1 - 1]
        omega2 = self.natural_frequencies[mode2 - 1]

        betaK = 2 * self.damping_ratio / (omega1 + omega2)
        alphaM = betaK * omega1 * omega2
        betaKinit = 0
        betaKcomm = 4 * self.damping_ratio / (omega1 + omega2)
        ops.rayleigh(alphaM, betaK, betaKinit, betaKcomm)

    def plot_stress_strain_curve(self, Target_elements_shear, Target_elements_flexural, section_num):
        # Plot Shear Column Hysteresis Loop
        fig_shear, axs_shear = plt.subplots(3, 3, figsize=(15, 15))
        fig_shear.suptitle(f'Hysteresis Curves for Shear Column Elements (Section {section_num})', fontsize=16)

        for idx, element_tag in enumerate(Target_elements_shear):
            file_path = os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt')
            stress_strain_data = np.loadtxt(file_path)

            strain = stress_strain_data[:, 0]
            stress = stress_strain_data[:, 1]

            row = idx // 3
            col = idx % 3
            axs_shear[row, col].plot(strain, stress, label=f'Element {element_tag}')
            axs_shear[row, col].set_xlabel('Strain')
            axs_shear[row, col].set_ylabel('Stress')
            axs_shear[row, col].set_title(f'Element {element_tag}')
            axs_shear[row, col].grid(True)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        img_path_shear = os.path.join(self.hysteresis_loop_folder, 'Stress_Strain_Hysteresis_Curves_Shear_Column.png')
        plt.savefig(img_path_shear)

        # Plot Flexural Column Hysteresis Loop
        fig_flexural, axs_flexural = plt.subplots(3, 3, figsize=(15, 15))
        fig_flexural.suptitle(f'Hysteresis Curves for Flexural Column Elements (Section {section_num})', fontsize=16)

        for idx, element_tag in enumerate(Target_elements_flexural):
            file_path = os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt')
            stress_strain_data = np.loadtxt(file_path)

            strain = stress_strain_data[:, 0]
            stress = stress_strain_data[:, 1]

            row = idx // 3
            col = idx % 3
            axs_flexural[row, col].plot(strain, stress, label=f'Element {element_tag}')
            axs_flexural[row, col].set_xlabel('Strain')
            axs_flexural[row, col].set_ylabel('Stress')
            axs_flexural[row, col].set_title(f'Element {element_tag}')
            axs_flexural[row, col].grid(True)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        img_path_flexural = os.path.join(self.hysteresis_loop_folder, 'Stress_Strain_Hysteresis_Curves_Flexural_Column.png')
        plt.savefig(img_path_flexural)

    def dynamic_analysis(self, data, DataType, add_leaning_column, Target_elements_shear, Target_elements_flexural,
                          section_num,  scale_factor, record_stress_strain):
        time = data[:, 0]
        acc = data[:, 1]

        self.acc = acc if self.Earthquake == 'Groundmotion' else acc * 9.81 / np.max(np.abs(acc * 9.81))

        dt = time[1] - time[0]
        self.DataType = DataType    
        ops.timeSeries('Path', 1, '-dt', dt, '-values', *self.acc, '-factor', scale_factor)
        ops.pattern('MultipleSupport',1)
        ops.groundMotion(1,'Plain','-accel',1)

        value = 3 if add_leaning_column else 2
            
        for Base_node in range(11, 11 + value):
            ops.imposedMotion(Base_node,1,1)

        self.main_folder_name = os.path.join(os.getcwd(),"Inelastic GBM", "Inelastic GBM Structure Response", self.Earthquake)
        if not os.path.exists(self.main_folder_name):
            os.makedirs(self.main_folder_name)

        self.abs_folder_name = os.path.join(self.main_folder_name, f'abs_{self.DataType}')

        if os.path.exists(self.abs_folder_name):
             shutil.rmtree(self.abs_folder_name)
        os.makedirs(self.abs_folder_name)

        for i in range(1, self.Story_num + 2):
            ops.recorder('Node', '-file', os.path.join(self.abs_folder_name, f'node{i * 10 + 1}{self.DataType}.txt'), '-time', '-node', i * 10 + 1, '-dof', 1, self.DataType)

        if record_stress_strain:    
            self.hysteresis_loop_folder = os.path.join(self.strain_stress_folder_name, f'{self.Earthquake}') 
            if os.path.exists(self.hysteresis_loop_folder):
                shutil.rmtree(self.hysteresis_loop_folder)
            os.makedirs(self.hysteresis_loop_folder)

            for element_tag in Target_elements_shear:
                ops.recorder('Element', '-file', os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt'), 
                            '-ele', element_tag, 'section', section_num, 'fiber', 1, 1, 'stressStrain')
    
            for element_tag in Target_elements_flexural:
                ops.recorder('Element', '-file', os.path.join(self.hysteresis_loop_folder, f'element_{element_tag}_stressStrain.txt'), 
                             '-ele', element_tag, 'section', section_num, 'fiber', 1, 1, 'stressStrain')
        
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        tolerance = 1 if scale_factor > 9.81 else 1e-1
        ops.test('EnergyIncr', tolerance, 100)
        ops.algorithm('Newton')
        ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')

        result = ops.analyze(len(self.acc), dt)
        if record_stress_strain:
            self.plot_stress_strain_curve(Target_elements_shear, Target_elements_flexural, section_num)

    def plot_dynamic_result(self, data):
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

        print(f"Complete the generation of {self.plot_title} data for {self.Earthquake} Earthquake (with leaning column)")

    def calculate_shear_stiffness(self, push_force, stick_type):
 
        print(f"\nShear Stiffness of each Floor ({stick_type.capitalize()} Stick)：")
    
        for floor in range(2, self.Story_num + 2):
            self.initialize_model()
            self.build_nodes()
            self.define_materials_and_sections()
            self.add_elements()
            self.assign_node_mass()

            node_offset = 1 if stick_type == 'shear' else 2
        
            for i in range(2, self.Story_num + 2):
                nodeID = 10 * i + node_offset
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
            ops.test('NormUnbalance', 1e-8, 10)
            ops.algorithm('Newton')
            ops.integrator('LoadControl', 1.0)
            ops.analysis('Static')

            push_node = 10 * floor + node_offset
            ops.timeSeries('Linear', 1)
            ops.pattern('Plain', 1, 1)
            ops.load(push_node, push_force, 0.0, 0.0)

            analyze_result = ops.analyze(1)
            if analyze_result != 0:
                print(f"Floor {floor - 1} analysis didn't converge")
                continue

            disp = ops.nodeDisp(push_node, 1)

            if disp != 0:
                stiffness = push_force / disp
                print(f"Floor {floor - 1}，Displacement：{disp:.6f} m，Shear Stiffness：{stiffness:.2f} N/m")
            else:
                print(f"Floor {floor - 1} total displacement is equal to zero，cannot calculate shear stiffness")

    def perform_cyclic_pushover_analysis(self, max_displacement, cycles, stick_type='shear'):
        elements_shear = [1121, 2131, 3141, 4151, 5161, 6171, 7181, 8191, 91101]
        hysteresis_data = []
        element_hysteresis_data = []

        for floor, element in enumerate(elements_shear, start=2):
            self.initialize_model()
            self.build_nodes()
            self.define_materials_and_sections()
            self.add_elements()
            self.assign_node_mass()

            node_offset = 1 if stick_type == 'shear' else 2

            for i in range(2, self.Story_num + 2):
                nodeID = 10 * i + node_offset
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
            ops.test('EnergyIncr', 1e-1, 100)
            ops.algorithm('Newton')

            push_node = 10 * floor + node_offset

            amplitudes = np.linspace(0, max_displacement, cycles + 1)[1:]  
            displacement_history = []
            for amp in amplitudes:
                displacement_history.extend([amp, -amp])

            ops.timeSeries('Linear', 1)
            ops.pattern('Plain', 1, 1)
            ops.load(push_node, 1.0, 0.0, 0.0)  

            recorder_folder = os.path.join('Inelastic GBM', 'hysteresis_data', f'floor_{floor}')
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
        #plt.show()

    
if __name__ == "__main__":

    start_time = time.time()
    Steel_Type = 'Steel01' # Steel01 or Steel02
    LCol = 3.9624     # m
    LCol_1 = 5.4864   # m
    LGird = 1.5       # m
    Story_num = 9       
    nps = 5
    Mass_2F = (1007304 / 2)                          # kg
    Mass_3_to_9F = (990088 / 2)                      # kg
    Mass_Roof = (1066540 / 2)                        # kg
    Floor_Mass = [Mass_2F, Mass_3_to_9F, Mass_Roof]  # kg
    Damping_ratio = 0.02
    section_num = 1  # nps section
    Target_elements_shear = [1121, 2131, 3141, 4151, 5161, 6171, 7181, 8191, 91101]
    Target_elements_flexural = [1222, 2232, 3242, 4252, 5262, 6272, 7282, 8292, 92102]

    # Flexural Stiffness (N/m) of Each Floor Obtain From Elastic GBM (From 9F to 1F)
    shear_stiffness_per_floor = [75171408, 86327220, 117553069, 168256700, 158692299, 193112556, 206133463, 223983896, 139013238]
    flexural_stiffness_per_floor = [31216417, 134002357, 151888442, 129929474, 164399055, 195467799, 165231962, 75308144, 12717634]

    Analysis_Type = 'disp'  # 'disp' 'vel' 'accel' 'rayleighForces'

    # Earthquake Type
    Earthquake_dict = {
    0: "Groundmotion.txt",
    1: "El Centro.txt",
    2: "Kobe.txt",
    3: "Chi-Chi.txt",
    4: "Northridge.txt"
    }
    scale_factor = 9.81 * 0.4
    Earthquake = Earthquake_dict[4]
    Earthquake_Data = np.loadtxt(os.path.join(os.getcwd(),"Earthquake", Earthquake))
    shear_Area_factor = [1.441, 1.285, 1.246, 1.225, 1.275, 1.294, 1.216, 1.126, 1.124]
    flexural_Area_factor = [0.792, 0.968, 1.178, 1.229, 1.286, 1.212, 1.297, 1.257, 1.605]
    
    yielding_factor = [1.7, 1.4, 1.4, 1.5, 1.5, 2, 3, 2.5, 3.5]

    add_leaning_column = True  # True or False

    model = Inelastic_GBM(Analysis_Type, Steel_Type, LCol, LCol_1, LGird, Story_num, nps, Floor_Mass, Damping_ratio, 
                          shear_Area_factor, flexural_Area_factor, Earthquake.replace(".txt", ""),  yielding_factor)
    model.build_nodes()
    model.define_materials_and_sections()
    model.add_elements()
    model.add_rigid_link()
    model.assign_node_mass()
    if add_leaning_column:
        model.add_leaning_column()
    model.calculate_natural_frequencies()
    model.calculate_mode_shape(Store_data = False)
    model.visualize_model(add_leaning_column)
    model.add_rayleigh_damping(1, 3)
    model.dynamic_analysis(Earthquake_Data, Analysis_Type, add_leaning_column,Target_elements_shear, Target_elements_flexural, 
                           section_num, scale_factor, record_stress_strain = False)
    model.plot_dynamic_result(Earthquake_Data)

    #model.calculate_shear_stiffness(push_force = 10000, stick_type='shear')    
    #model.calculate_shear_stiffness(push_force = 10000,stick_type='flexural')
    #model.perform_cyclic_pushover_analysis(max_displacement=0.05, cycles=15, stick_type='shear')
    plt.show()

    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)  
    minutes, seconds = divmod(remainder, 60)
    print(f"Total runtime: {int(hours)} hours {int(minutes)} minutes {int(seconds)} seconds")
